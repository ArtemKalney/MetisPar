/*
 * Copyright 1997-2015, Regents of the University of Minnesota
 *
 * ometis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id: ometis.c 18208 2015-01-17 18:02:38Z dominique $
 *
 *
 * Dominique LaSalle - 2015-02-28
 * Modified to support parallel task scheduling.
 *
 */

#include <amd.h>
#include "metislib.h"
// #include "dlthread_pool.h"


//#define OMP_TASK 1

static idx_t const MIN_TASK_SIZE = 2048;


static ctrl_t * __duplicate_ctrl(
    ctrl_t const * const ctrl,
    graph_t * const graph)
{
  ctrl_t * myctrl;

  myctrl = gk_malloc(sizeof(ctrl_t),"X"); 
  memcpy(myctrl,ctrl,sizeof(ctrl_t));

  myctrl->mcore = NULL;
  myctrl->cnbrpool = NULL;
  myctrl->vnbrpool = NULL;
  myctrl->maxnads = NULL;
  myctrl->nads = NULL;
  myctrl->adids = NULL;
  myctrl->adwgts = NULL;
  myctrl->pvec1 = NULL;
  myctrl->pvec2 = NULL;

  myctrl->maxvwgt = ismalloc(ctrl->ncon, 0, "SetupCtrl: maxvwgt");

  myctrl->tpwgts = rsmalloc(2, .5, "ctrl->tpwgts");

  myctrl->ubfactors = rsmalloc(ctrl->ncon, \
      I2RUBFACTOR(ctrl->ufactor), "SetupCtrl: ubfactors");

  myctrl->pijbm = rmalloc(ctrl->nparts, "SetupCtrl: ctrl->pijbm");

  AllocateWorkSpace(myctrl, graph);

  return myctrl;
}


/*************************************************************************/
/*! This function is the entry point for the multilevel nested dissection 
    ordering code. At each bisection, a node-separator is computed using
    a node-based refinement approach.

    \param nvtxs is the number of vertices in the graph.
    \param xadj is of length nvtxs+1 marking the start of the adjancy 
           list of each vertex in adjncy.
    \param adjncy stores the adjacency lists of the vertices. The adjnacy
           list of a vertex should not contain the vertex itself.
    \param vwgt is an array of size nvtxs storing the weight of each 
           vertex. If vwgt is NULL, then the vertices are considered 
           to have unit weight.
    \param numflag is either 0 or 1 indicating that the numbering of 
           the vertices starts from 0 or 1, respectively.
    \param options is an array of size METIS_NOPTIONS used to pass 
           various options impacting the of the algorithm. A NULL
           value indicates use of default options.
    \param perm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A'[i] = A[perm[i]].
    \param iperm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A[i] = A'[iperm[i]].
*/
/*************************************************************************/
int METIS_NodeND(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
          idx_t *options, idx_t *perm, idx_t *iperm) 
{
  int sigrval=0, renumber=0;
  idx_t i, ii, j, l, nnvtxs=0;
  graph_t *graph=NULL;
  ctrl_t *ctrl;
  idx_t *cptr, *cind, *piperm;
  int numflag = 0;

  /* set up the run time parameters */
  ctrl = SetupCtrl(METIS_OP_OMETIS, options, 1, 3, NULL, NULL);
  if (!ctrl) {
    return METIS_ERROR_INPUT;
  }

  /* if required, change the numbering to 0 */
  if (ctrl->numflag == 1) {
    Change2CNumbering(*nvtxs, xadj, adjncy);
    renumber = 1;
  }

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->TotalTmr));

  /* prune the dense columns */
  if (ctrl->pfactor > 0.0) { 
    piperm = imalloc(*nvtxs, "OMETIS: piperm");

    graph = PruneGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, piperm, ctrl->pfactor);
    if (graph == NULL) {
      /* if there was no prunning, cleanup the pfactor */
      gk_free((void **)&piperm, LTERM);
      ctrl->pfactor = 0.0;
    }
    else {
      nnvtxs = graph->nvtxs;
      ctrl->compress = 0;  /* disable compression if prunning took place */
    }
  }

  /* compress the graph; note that compression only happens if not prunning 
     has taken place. */
  if (ctrl->compress) { 
    cptr = imalloc(*nvtxs+1, "OMETIS: cptr");
    cind = imalloc(*nvtxs, "OMETIS: cind");

    graph = CompressGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, cptr, cind);
    if (graph == NULL) {
      /* if there was no compression, cleanup the compress flag */
      gk_free((void **)&cptr, &cind, LTERM);
      ctrl->compress = 0; 
    }
    else {
      nnvtxs = graph->nvtxs;
      ctrl->cfactor = 1.0*(*nvtxs)/nnvtxs;
      if (ctrl->cfactor > 1.5 && ctrl->nseps == 1)
        ctrl->nseps = 2;
      //ctrl->nseps = (idx_t)(ctrl->cfactor*ctrl->nseps);
    }
  }

  /* if no prunning and no compression, setup the graph in the normal way. */
  if (ctrl->pfactor == 0.0 && ctrl->compress == 0) 
    graph = SetupGraph(ctrl, *nvtxs, 1, xadj, adjncy, vwgt, NULL, NULL);

  ASSERT(CheckGraph(graph, ctrl->numflag, 1));

  /* allocate workspace memory */
  AllocateWorkSpace(ctrl, graph);

  /* do the nested dissection ordering  */
  if (ctrl->ccorder) 
  {
    MlevelNestedDissectionCC(ctrl, graph, iperm, graph->nvtxs);
  }
  else 
  {
    #pragma omp parallel 
    {
        #pragma omp single 
        {
            MlevelNestedDissection(ctrl, graph, iperm, graph->nvtxs);
        }
    }
  }

  if (ctrl->pfactor > 0.0) { /* Order any prunned vertices */
    icopy(nnvtxs, iperm, perm);  /* Use perm as an auxiliary array */
    for (i=0; i<nnvtxs; i++)
      iperm[piperm[i]] = perm[i];
    for (i=nnvtxs; i<*nvtxs; i++)
      iperm[piperm[i]] = i;

    gk_free((void **)&piperm, LTERM);
  }
  else if (ctrl->compress) { /* Uncompress the ordering */
    /* construct perm from iperm */
    for (i=0; i<nnvtxs; i++)
      perm[iperm[i]] = i; 
    for (l=ii=0; ii<nnvtxs; ii++) {
      i = perm[ii];
      for (j=cptr[i]; j<cptr[i+1]; j++)
        iperm[cind[j]] = l++;
    }

    gk_free((void **)&cptr, &cind, LTERM);
  }

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopwctimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

  /* clean up */
  FreeCtrl(&ctrl);

  /* if required, change the numbering back to 1 */
  if (renumber)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);

  return metis_rcode(sigrval);
}

void LowLevelReordering(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t lastvtx)
{
    if (ctrl->ltype == METIS_MMD)
    {
        MMDOrder(ctrl, graph, order, lastvtx);
    }

    if (ctrl->ltype == METIS_AMD)
    {
        AMDOrder(ctrl, graph, order, lastvtx);
    }

    if (ctrl->ltype == METIS_RCM)
    {
        RCMOrder(ctrl, graph, order, lastvtx);
    }

    if (ctrl->ltype == METIS_EMP)
    {
        EmptyOrder(ctrl, graph, order, lastvtx);
    }

    FreeGraph(&graph);
}

/*************************************************************************/
/*! This is the driver for the recursive tri-section of a graph into the
    left, separator, and right partitions. The graphs correspond to the 
    left and right parts are further tri-sected in a recursive fashion.
    The nodes in the separator are ordered at the end of the left & right
    nodes.
 */
/*************************************************************************/
void MlevelNestedDissection(ctrl_t* ctrl, graph_t* graph, idx_t* order,
    idx_t lastvtx)
{
    idx_t i, j, nvtxs, nbnd, rvtx, lvtx;
    idx_t* label, * bndind;
    graph_t* lgraph, * rgraph;

    nvtxs = graph->nvtxs;

    MlevelNodeBisectionMultiple(ctrl, graph);

    IFSET(ctrl->dbglvl, METIS_DBG_SEPINFO,
        printf("Nvtxs: %6"PRIDX", [%6"PRIDX" %6"PRIDX" %6"PRIDX"]\n",
            graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));


    /* Order the nodes in the separator */
    nbnd = graph->nbnd;
    bndind = graph->bndind;
    label = graph->label;
    for (i = 0; i < nbnd; i++)
        order[label[bndind[i]]] = --lastvtx;

    SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

    /* Free the memory of the top level graph */
    FreeGraph(&graph);

    rvtx = lastvtx;
    lvtx = lastvtx - rgraph->nvtxs;

    /* Recurse on lgraph first, as its lastvtx depends on rgraph->nvtxs, which
       will not be defined upon return from MlevelNestedDissection. */
    if (lgraph->nedges > ctrl->minnedges && lgraph->nvtxs > ctrl->minnvtxs)
    {
        #pragma omp task default(none) firstprivate(lgraph, ctrl, order, lvtx)
        {
            ctrl_t *myctrl;

            myctrl = __duplicate_ctrl(ctrl, lgraph);

            if (lgraph->nvtxs > ctrl->llswitch && lgraph->nedges > 0)
            {
                MlevelNestedDissection(myctrl, lgraph, order, lvtx);
            }
            else
            {
                LowLevelReordering(myctrl, lgraph, order, lvtx);
            }

            // FreeCtrl(&myctrl);
        }
    }
    else
    {
        if (lgraph->nvtxs > ctrl->llswitch && lgraph->nedges > 0)
        {
            MlevelNestedDissection(ctrl, lgraph, order, lvtx);
        }
        else
        {
            LowLevelReordering(ctrl, lgraph, order, lvtx);
        }
    }

    // if (rgraph->nedges > ctrl->minnedges && rgraph->nvtxs > ctrl->minnvtxs)
    // {
    //     #pragma omp task default(none) shared(ctrl, rgraph, order, rvtx)
    //     {
    //         ctrl_t *myctrl;

    //         myctrl = __duplicate_ctrl(ctrl, rgraph);

    //         if (rgraph->nvtxs > ctrl->llswitch && rgraph->nedges > 0)
    //         {
    //             MlevelNestedDissection(myctrl, rgraph, order, rvtx);
    //         }
    //         else
    //         {
    //             LowLevelReordering(myctrl, rgraph, order, rvtx);
    //         }

    //         FreeCtrl(&myctrl);
    //     }
    // }
    // else
    {
        if (rgraph->nvtxs > ctrl->llswitch && rgraph->nedges > 0)
        {
            MlevelNestedDissection(ctrl, rgraph, order, rvtx);
        }
        else
        {
            LowLevelReordering(ctrl, rgraph, order, rvtx);
        }
    }

    // #pragma omp taskwait
}


/*************************************************************************/
/*! This routine is similar to its non 'CC' counterpart. The difference is
    that after each tri-section, the connected components of the original
    graph that result after removing the separator vertises are ordered
    independently (i.e., this may lead to more than just the left and 
    the right subgraphs).
*/
/*************************************************************************/
void MlevelNestedDissectionCC(ctrl_t* ctrl, graph_t* graph, idx_t* order,
    idx_t lastvtx)
{
    idx_t i, j, nvtxs, nbnd, ncmps, rnvtxs, snvtxs;
    idx_t* label, * bndind, * offset;
    idx_t* cptr, * cind;
    graph_t** sgraphs;

    nvtxs = graph->nvtxs;

    MlevelNodeBisectionMultiple(ctrl, graph);

    IFSET(ctrl->dbglvl, METIS_DBG_SEPINFO,
        printf("Nvtxs: %6"PRIDX", [%6"PRIDX" %6"PRIDX" %6"PRIDX"]\n",
            graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

    /* Order the nodes in the separator */
    nbnd = graph->nbnd;
    bndind = graph->bndind;
    label = graph->label;
    for (i = 0; i < nbnd; i++)
        order[label[bndind[i]]] = --lastvtx;

    WCOREPUSH;
    cptr = iwspacemalloc(ctrl, nvtxs + 1);
    cind = iwspacemalloc(ctrl, nvtxs);
    ncmps = FindSepInducedComponents(ctrl, graph, cptr, cind);

    if (ctrl->dbglvl & METIS_DBG_INFO) {
        if (ncmps > 2)
            printf("  Bisection resulted in %"PRIDX" connected components\n", ncmps);
    }

    sgraphs = SplitGraphOrderCC(ctrl, graph, ncmps, cptr, cind);

    WCOREPOP;

    /* Free the memory of the top level graph */
    FreeGraph(&graph);

    omp_set_nested(1);

    offset = imalloc(ncmps, "X");
    offset[0] = 0;
    for (i = 1; i < ncmps; ++i)
    {
        offset[i] = sgraphs[i - 1]->nvtxs + offset[i - 1];
    }

    for (i = 0; i < ncmps; ++i)
    {
/* Go and process the subgraphs */
#pragma omp task shared(ctrl, sgraphs, order, offset, lastvtx, i) default(none)
        {
            idx_t rnvtxs;

            /* Save the number of vertices in sgraphs[i] because sgraphs[i] is freed
               inside MlevelNestedDissectionCC, and as such it will be undefined. */
            ctrl_t *myctrl;

            /* create my control */
            myctrl = __duplicate_ctrl(ctrl, sgraphs[i]);

            rnvtxs = offset[i];

            if (sgraphs[i]->nvtxs > ctrl->llswitch && sgraphs[i]->nedges > 0)
            {
                MlevelNestedDissectionCC(myctrl, sgraphs[i], order, lastvtx - rnvtxs);
            }
            else
            {
                LowLevelReordering(myctrl, sgraphs[i], order, lastvtx - rnvtxs);
            }

            FreeCtrl(&myctrl);
        }
    }
#pragma omp taskwait

    gk_free((void **)&sgraphs, &offset, LTERM);
}


/*************************************************************************/
/*! This function performs multilevel node bisection (i.e., tri-section).
    It performs multiple bisections and selects the best. */
/*************************************************************************/
void MlevelNodeBisectionMultiple(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, mincut;
  idx_t *bestwhere;

  /* if the graph is small, just find a single vertex separator */
  if (ctrl->nseps == 1 || graph->nvtxs < (ctrl->compress ? 1000 : 2000)) {
    MlevelNodeBisectionL2(ctrl, graph, ctrl->nIparts);
    return;
  }

  WCOREPUSH;

  bestwhere = iwspacemalloc(ctrl, graph->nvtxs);

  mincut = graph->tvwgt[0];
  for (i=0; i<ctrl->nseps; i++) {
    MlevelNodeBisectionL2(ctrl, graph, ctrl->nIparts);

    if (i == 0 || graph->mincut < mincut) {
      mincut = graph->mincut;
      if (i < ctrl->nseps-1)
        icopy(graph->nvtxs, graph->where, bestwhere);
    }

    if (mincut == 0)
      break;

    if (i < ctrl->nseps-1) 
      FreeRData(graph);
  }

  if (mincut != graph->mincut) {
    icopy(graph->nvtxs, bestwhere, graph->where);
    Compute2WayNodePartitionParams(ctrl, graph);
  }

  WCOREPOP;
}


/*************************************************************************/
/*! This function performs multilevel node bisection (i.e., tri-section).
    It performs multiple bisections and selects the best. */
/*************************************************************************/
void MlevelNodeBisectionL2(ctrl_t *ctrl, graph_t *graph, idx_t niparts)
{
  idx_t i, mincut, nruns=ctrl->nrunsmlndl1;
  graph_t *cgraph; 
  idx_t *bestwhere;

  /* if the graph is small, just find a single vertex separator */
  if (graph->nvtxs < ctrl->smallgr) {
    MlevelNodeBisectionL1(ctrl, graph, niparts);
    return;
  }

  WCOREPUSH;

  ctrl->CoarsenTo = gk_max(100, graph->nvtxs/30);

  cgraph = CoarsenGraphNlevels(ctrl, graph, ctrl->ncoarsengr);

  bestwhere = iwspacemalloc(ctrl, cgraph->nvtxs);

  mincut = graph->tvwgt[0];
  for (i=0; i<nruns; i++) {
    MlevelNodeBisectionL1(ctrl, cgraph, 0.7*niparts);

    if (i == 0 || cgraph->mincut < mincut) {
      mincut = cgraph->mincut;
      if (i < nruns-1)
        icopy(cgraph->nvtxs, cgraph->where, bestwhere);
    }

    if (mincut == 0)
      break;

    if (i < nruns-1) 
      FreeRData(cgraph);
  }

  if (mincut != cgraph->mincut) 
    icopy(cgraph->nvtxs, bestwhere, cgraph->where);

  WCOREPOP;

  Refine2WayNode(ctrl, graph, cgraph);

}


/*************************************************************************/
/*! The top-level routine of the actual multilevel node bisection */
/*************************************************************************/
void MlevelNodeBisectionL1(ctrl_t *ctrl, graph_t *graph, idx_t niparts)
{
  graph_t *cgraph;

  ctrl->CoarsenTo = graph->nvtxs/8;
  if (ctrl->CoarsenTo > 100)
    ctrl->CoarsenTo = 100;
  else if (ctrl->CoarsenTo < 40)
    ctrl->CoarsenTo = 40;

  cgraph = CoarsenGraph(ctrl, graph);

  niparts = gk_max(1, (cgraph->nvtxs <= ctrl->CoarsenTo ? niparts/2: niparts));
  /*niparts = (cgraph->nvtxs <= ctrl->CoarsenTo ? SMALLNIPARTS : LARGENIPARTS);*/
  InitSeparator(ctrl, cgraph, niparts);

  Refine2WayNode(ctrl, graph, cgraph);
}


/*************************************************************************/
/*! This function takes a graph and a tri-section (left, right, separator)
    and splits it into two graphs. 
    
    This function relies on the fact that adjwgt is all equal to 1.
*/
/*************************************************************************/
void SplitGraphOrder(ctrl_t *ctrl, graph_t *graph, graph_t **r_lgraph, 
         graph_t **r_rgraph)
{
  idx_t i, ii, j, k, l, istart, iend, mypart, nvtxs, snvtxs[3], snedges[3];
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where, *bndptr, *bndind;
  idx_t *sxadj[2], *svwgt[2], *sadjncy[2], *sadjwgt[2], *slabel[2];
  idx_t *rename;
  idx_t *auxadjncy;
  graph_t *lgraph, *rgraph;

  WCOREPUSH;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->SplitTmr));

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  label   = graph->label;
  where   = graph->where;
  bndptr  = graph->bndptr;
  bndind  = graph->bndind;
  ASSERT(bndptr != NULL);

  rename = iwspacemalloc(ctrl, nvtxs);
  
  snvtxs[0] = snvtxs[1] = snvtxs[2] = snedges[0] = snedges[1] = snedges[2] = 0;
  for (i=0; i<nvtxs; i++) {
    k = where[i];
    rename[i] = snvtxs[k]++;
    snedges[k] += xadj[i+1]-xadj[i];
  }

  lgraph      = SetupSplitGraph(graph, snvtxs[0], snedges[0]);
  sxadj[0]    = lgraph->xadj;
  svwgt[0]    = lgraph->vwgt;
  sadjncy[0]  = lgraph->adjncy; 
  sadjwgt[0]  = lgraph->adjwgt; 
  slabel[0]   = lgraph->label;

  rgraph      = SetupSplitGraph(graph, snvtxs[1], snedges[1]);
  sxadj[1]    = rgraph->xadj;
  svwgt[1]    = rgraph->vwgt;
  sadjncy[1]  = rgraph->adjncy; 
  sadjwgt[1]  = rgraph->adjwgt; 
  slabel[1]   = rgraph->label;

  /* Go and use bndptr to also mark the boundary nodes in the two partitions */
  for (ii=0; ii<graph->nbnd; ii++) {
    i = bndind[ii];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      bndptr[adjncy[j]] = 1;
  }

  snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
  sxadj[0][0] = sxadj[1][0] = 0;
  for (i=0; i<nvtxs; i++) {
    if ((mypart = where[i]) == 2)
      continue;

    istart = xadj[i];
    iend   = xadj[i+1];
    if (bndptr[i] == -1) { /* This is an interior vertex */
      auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
      for(j=istart; j<iend; j++) 
        auxadjncy[j] = adjncy[j];
      snedges[mypart] += iend-istart;
    }
    else {
      auxadjncy = sadjncy[mypart];
      l = snedges[mypart];
      for (j=istart; j<iend; j++) {
        k = adjncy[j];
        if (where[k] == mypart) 
          auxadjncy[l++] = k;
      }
      snedges[mypart] = l;
    }

    svwgt[mypart][snvtxs[mypart]]    = vwgt[i];
    slabel[mypart][snvtxs[mypart]]   = label[i];
    sxadj[mypart][++snvtxs[mypart]]  = snedges[mypart];
  }

  for (mypart=0; mypart<2; mypart++) {
    iend = snedges[mypart];
    iset(iend, 1, sadjwgt[mypart]);

    auxadjncy = sadjncy[mypart];
    for (i=0; i<iend; i++) 
      auxadjncy[i] = rename[auxadjncy[i]];
  }

  lgraph->nvtxs  = snvtxs[0];
  lgraph->nedges = snedges[0];
  rgraph->nvtxs  = snvtxs[1];
  rgraph->nedges = snedges[1];

  SetupGraph_tvwgt(lgraph);
  SetupGraph_tvwgt(rgraph);

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopwctimer(ctrl->SplitTmr));

  *r_lgraph = lgraph;
  *r_rgraph = rgraph;

  WCOREPOP;
}


/*************************************************************************/
/*! This function takes a graph and generates a set of graphs, each of 
    which is a connected component in the original graph.

    This function relies on the fact that adjwgt is all equal to 1.

    \param ctrl stores run state info.
    \param graph is the graph to be split.
    \param ncmps is the number of connected components.
    \param cptr is an array of size ncmps+1 that marks the start and end
           locations of the vertices in cind that make up the respective
           components (i.e., cptr, cind is in CSR format).
    \param cind is an array of size equal to the number of vertices in 
           the original graph and stores the vertices that belong to each
           connected component.

    \returns an array of subgraphs corresponding to the extracted subgraphs.
*/
/*************************************************************************/
graph_t **SplitGraphOrderCC(ctrl_t *ctrl, graph_t *graph, idx_t ncmps, 
              idx_t *cptr, idx_t *cind)
{
  idx_t i, ii, iii, j, k, l, istart, iend, mypart, nvtxs, snvtxs, snedges;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *label, *where, *bndptr, *bndind;
  idx_t *sxadj, *svwgt, *sadjncy, *sadjwgt, *slabel;
  idx_t *rename;
  idx_t *auxadjncy;
  graph_t **sgraphs;

  WCOREPUSH;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->SplitTmr));

  nvtxs   = graph->nvtxs;
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  label   = graph->label;
  where   = graph->where;
  bndptr  = graph->bndptr;
  bndind  = graph->bndind;
  ASSERT(bndptr != NULL);

  /* Go and use bndptr to also mark the boundary nodes in the two partitions */
  for (ii=0; ii<graph->nbnd; ii++) {
    i = bndind[ii];
    for (j=xadj[i]; j<xadj[i+1]; j++)
      bndptr[adjncy[j]] = 1;
  }

  rename = iwspacemalloc(ctrl, nvtxs);
  
  sgraphs = (graph_t **)gk_malloc(sizeof(graph_t *)*ncmps, "SplitGraphOrderCC: sgraphs");

  /* Go and split the graph a component at a time */
  for (iii=0; iii<ncmps; iii++) {
    my_irandArrayPermute_r(cptr[iii+1]-cptr[iii], cind+cptr[iii], cptr[iii+1]-cptr[iii], 0, &ctrl->curseed);
    snvtxs = snedges = 0;
    for (j=cptr[iii]; j<cptr[iii+1]; j++) {
      i = cind[j];
      rename[i] = snvtxs++;
      snedges += xadj[i+1]-xadj[i];
    }

    sgraphs[iii] = SetupSplitGraph(graph, snvtxs, snedges);

    sxadj    = sgraphs[iii]->xadj;
    svwgt    = sgraphs[iii]->vwgt;
    sadjncy  = sgraphs[iii]->adjncy;
    sadjwgt  = sgraphs[iii]->adjwgt;
    slabel   = sgraphs[iii]->label;

    snvtxs = snedges = sxadj[0] = 0;
    for (ii=cptr[iii]; ii<cptr[iii+1]; ii++) {
      i = cind[ii];

      istart = xadj[i];
      iend   = xadj[i+1];
      if (bndptr[i] == -1) { /* This is an interior vertex */
        auxadjncy = sadjncy + snedges - istart;
        for(j=istart; j<iend; j++) 
          auxadjncy[j] = adjncy[j];
        snedges += iend-istart;
      }
      else {
        l = snedges;
        for (j=istart; j<iend; j++) {
          k = adjncy[j];
          if (where[k] != 2) 
            sadjncy[l++] = k;
        }
        snedges = l;
      }

      svwgt[snvtxs]    = vwgt[i];
      slabel[snvtxs]   = label[i];
      sxadj[++snvtxs]  = snedges;
    }

    iset(snedges, 1, sadjwgt);
    for (i=0; i<snedges; i++) 
      sadjncy[i] = rename[sadjncy[i]];

    sgraphs[iii]->nvtxs  = snvtxs;
    sgraphs[iii]->nedges = snedges;

    SetupGraph_tvwgt(sgraphs[iii]);
  }

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopwctimer(ctrl->SplitTmr));

  WCOREPOP;

  return sgraphs;
}


/*************************************************************************/
/*! This function uses MMD to order the graph. The vertices are numbered
    from lastvtx downwards. */
/*************************************************************************/
void MMDOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t lastvtx)
{
  idx_t i, j, k, nvtxs, nofsub, firstvtx;
  idx_t *xadj, *adjncy, *label;
  idx_t *perm, *iperm, *head, *qsize, *list, *marker;

  WCOREPUSH;

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->LowLevelTmr))
	IFSET(ctrl->dbglvl, METIS_DBG_TIME, ctrl->LowLevelCount++)

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;

  /* Relabel the vertices so that it starts from 1 */
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]++;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]++;

  perm   = iwspacemalloc(ctrl, nvtxs+5);
  iperm  = iwspacemalloc(ctrl, nvtxs+5);
  head   = iwspacemalloc(ctrl, nvtxs+5);
  qsize  = iwspacemalloc(ctrl, nvtxs+5);
  list   = iwspacemalloc(ctrl, nvtxs+5);
  marker = iwspacemalloc(ctrl, nvtxs+5);

  genmmd(nvtxs, xadj, adjncy, iperm, perm, 1, head, qsize, list, marker, IDX_MAX, &nofsub);

  label = graph->label;
  firstvtx = lastvtx-nvtxs;
  for (i=0; i<nvtxs; i++)
    order[label[i]] = firstvtx+iperm[i]-1;

  /* Relabel the vertices so that it starts from 0 */
  for (i=0; i<nvtxs+1; i++)
    xadj[i]--;
  k = xadj[nvtxs];
  for (i=0; i<k; i++)
    adjncy[i]--;

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->LowLevelTmr))

  WCOREPOP;
}

int compare_idx_t(const void* a, const void* b)
{
    idx_t arg1 = *(const idx_t*)a;
    idx_t arg2 = *(const idx_t*)b;

    if (arg1 < arg2)
    {
	    return -1;
    }

    if (arg1 > arg2)
    {
	    return 1;
    }

    return 0;
}

void AMDOrder(ctrl_t* ctrl, graph_t* graph, idx_t* order, idx_t lastvtx)
{
    idx_t i, adjncy_index, adjncy_count, nvtxs, firstvtx, nedges;
    idx_t* xadj, * adjncy, * label;
    idx_t* perm, * iperm;

    idx_t status;
    idx_t* xadjProcessed, * adjncyProcessed;
    double control[AMD_CONTROL], info[AMD_INFO];

    idx_t * degrees, * mask_lp, * mask_ep, * mask, *w;
    my_graph_t* my_graph, * my_graph_e, * my_graph_l, * my_graph_h, * my_graph_s;

    WCOREPUSH;

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startwctimer(ctrl->LowLevelTmr));
    IFSET(ctrl->dbglvl, METIS_DBG_TIME, ctrl->LowLevelCount++);

    nvtxs = graph->nvtxs;
    xadj = graph->xadj;
    adjncy = graph->adjncy;
    nedges = graph->nedges;

    perm = iwspacemalloc(ctrl, nvtxs + 5);
    iperm = iwspacemalloc(ctrl, nvtxs + 5);

    xadjProcessed = iwspacemalloc(ctrl, nvtxs + 5);
    adjncyProcessed = iwspacemalloc(ctrl, nedges + 5);

    amd_defaults(control);
    amd_preprocess(ctrl, nvtxs, xadj, adjncy, xadjProcessed, adjncyProcessed);
    // using qsort for preprocess csr matrix
    /*adjncy_index = 0;
    for (i = 0; i < nvtxs; i++)
    {
        adjncy_count = xadj[i + 1] - xadj[i];
        qsort(adjncy + adjncy_index, adjncy_count, sizeof(idx_t), compare_idx_t);
        adjncy_index += adjncy_count;
    }*/
    status = amd_order(ctrl, nvtxs, xadjProcessed, adjncyProcessed, perm, control, info);

    switch (status) {
    case AMD_OK: IFSET(ctrl->dbglvl, METIS_DBG_LOW_LEVEL, printf("amd: Success\n")) break;
    case AMD_INVALID: IFSET(ctrl->dbglvl, METIS_DBG_LOW_LEVEL, printf("amd: Input error!\n")) break;
    case AMD_OUT_OF_MEMORY: IFSET(ctrl->dbglvl, METIS_DBG_LOW_LEVEL, printf("amd: Memory allocation error!\n")) break;
    default: IFSET(ctrl->dbglvl, METIS_DBG_LOW_LEVEL, printf("amd: Unknown error!\n")) break;
    }

    IFSET(ctrl->dbglvl, METIS_DBG_LOW_LEVEL, amd_control(control));
    IFSET(ctrl->dbglvl, METIS_DBG_LOW_LEVEL, amd_info(info));

	for (i = 0; i < nvtxs; i++)
    {
        iperm[perm[i]] = i;
    }

    // custom amd algorithm
    /*my_graph = CreateMyGraphWspace(ctrl, nvtxs, nedges + 1);
    my_graph_e = CreateMyGraphWspace(ctrl, nvtxs, nedges / 2 + 1);
    my_graph_l = CreateMyGraphWspace(ctrl, nvtxs, nedges / 2 + 1);
    my_graph_h = CreateMyGraphWspace(ctrl, nvtxs, nvtxs + 1);
    my_graph_s = CreateMyGraphWspace(ctrl, nvtxs, nvtxs + 1);

    degrees = iwspacemalloc(ctrl, nvtxs);
    mask_lp = iwspacemalloc(ctrl, nvtxs);
    mask_ep = iwspacemalloc(ctrl, nvtxs);
    mask = iwspacemalloc(ctrl, nvtxs);
    w = iwspacemalloc(ctrl, nvtxs);

    amd(xadj, adjncy, perm, iperm, my_graph, my_graph_e, my_graph_l, my_graph_h, my_graph_s, degrees, mask_lp, mask_ep, mask, w);*/

    label = graph->label;
    firstvtx = lastvtx - nvtxs;
    for (i = 0; i < nvtxs; i++)
    {
	    order[label[i]] = firstvtx + iperm[i];
    }

    IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->LowLevelTmr));

    WCOREPOP;
}