/*!
\file 
\brief Functions dealing with memory allocation and workspace management

\date Started 2/24/96
\author George
\author Copyright 1997-2009, Regents of the University of Minnesota 
\version $Id: wspace.c 17623 2014-09-09 03:51:44Z dominique $
*/

#include "metislib.h"


/*************************************************************************/
/*! This function allocates memory for the workspace */
/*************************************************************************/
void AllocateWorkSpace(ctrl_t *ctrl, graph_t *graph)
{
  size_t coresize;

  /*
  switch (ctrl->optype) {
    case METIS_OP_PMETIS:
      coresize = 3*(graph->nvtxs+1)*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(real_t);
      break;
    default:
      coresize = 4*(graph->nvtxs+1)*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(idx_t) + 
                 5*(ctrl->nparts+1)*graph->ncon*sizeof(real_t);
  }*/
  coresize = 0;
  ctrl->mcore = gk_mcoreCreate(coresize);

  ctrl->nbrpoolsize = 0;
  ctrl->nbrpoolcpos = 0;
}


/*************************************************************************/
/*! This function frees the workspace */
/*************************************************************************/
void FreeWorkSpace(ctrl_t *ctrl)
{
  gk_mcoreDestroy(&ctrl->mcore, ctrl->dbglvl&METIS_DBG_INFO);

  IFSET(ctrl->dbglvl, METIS_DBG_INFO,
      printf(" nbrpool statistics\n" 
             "        nbrpoolsize: %12zu   nbrpoolcpos: %12zu\n"
             "    nbrpoolreallocs: %12zu\n\n",
             ctrl->nbrpoolsize,  ctrl->nbrpoolcpos, 
             ctrl->nbrpoolreallocs));

  gk_free((void **)&ctrl->cnbrpool, &ctrl->vnbrpool, LTERM);
  ctrl->nbrpoolsize = 0;
  ctrl->nbrpoolcpos = 0;

  if (ctrl->minconn) {
    iFreeMatrix(&(ctrl->adids),  ctrl->nparts, INIT_MAXNAD);
    iFreeMatrix(&(ctrl->adwgts), ctrl->nparts, INIT_MAXNAD);

    gk_free((void **)&ctrl->pvec1, &ctrl->pvec2, 
        &ctrl->maxnads, &ctrl->nads, LTERM);
  }
}


/*************************************************************************/
/*! This function allocate space from the workspace/heap */
/*************************************************************************/
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes)
{
  return gk_mcoreMalloc(ctrl->mcore, nbytes);
}


/*************************************************************************/
/*! This function sets a marker in the stack of malloc ops to be used
    subsequently for freeing purposes */
/*************************************************************************/
void wspacepush(ctrl_t *ctrl)
{
  gk_mcorePush(ctrl->mcore);
}


/*************************************************************************/
/*! This function frees all mops since the last push */
/*************************************************************************/
void wspacepop(ctrl_t *ctrl)
{
  gk_mcorePop(ctrl->mcore);
}


/*************************************************************************/
/*! This function allocate space from the core  */
/*************************************************************************/
idx_t *iwspacemalloc(ctrl_t *ctrl, idx_t n)
{
  return (idx_t *)wspacemalloc(ctrl, n*sizeof(idx_t));
}


/*************************************************************************/
/*! This function allocate space from the core */
/*************************************************************************/
real_t *rwspacemalloc(ctrl_t *ctrl, idx_t n)
{
  return (real_t *)wspacemalloc(ctrl, n*sizeof(real_t));
}


/*************************************************************************/
/*! This function allocate space from the core  */
/*************************************************************************/
ikv_t *ikvwspacemalloc(ctrl_t *ctrl, idx_t n)
{
  return (ikv_t *)wspacemalloc(ctrl, n*sizeof(ikv_t));
}