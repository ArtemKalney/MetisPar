/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h 17622 2014-09-09 03:27:49Z dominique $
 *
 */

#ifndef _LIBMETIS_PROTO_H_
#define _LIBMETIS_PROTO_H_

/* auxapi.c */

/* balance.c */
void Balance2Way(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts);
void Bnd2WayBalance(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts);
void General2WayBalance(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts);
void McGeneral2WayBalance(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts);


/* bucketsort.c */
void BucketSortKeysInc(ctrl_t *ctrl, idx_t n, idx_t max, idx_t *keys,
         idx_t *tperm, idx_t *perm);

/* coarsen.c */
graph_t *CoarsenGraph(ctrl_t *ctrl, graph_t *graph);
graph_t *CoarsenGraphNlevels(ctrl_t *ctrl, graph_t *graph, idx_t nlevels);
idx_t Match_RM(ctrl_t *ctrl, graph_t *graph);
idx_t Match_SHEM(ctrl_t *ctrl, graph_t *graph);
idx_t Match_2Hop(ctrl_t *ctrl, graph_t *graph, idx_t *perm, idx_t *match,
          idx_t cnvtxs, size_t nunmatched);
idx_t Match_2HopAny(ctrl_t *ctrl, graph_t *graph, idx_t *perm, idx_t *match,
          idx_t cnvtxs, size_t *r_nunmatched, size_t maxdegree);
idx_t Match_2HopAll(ctrl_t *ctrl, graph_t *graph, idx_t *perm, idx_t *match,
          idx_t cnvtxs, size_t *r_nunmatched, size_t maxdegree);
void PrintCGraphStats(ctrl_t *ctrl, graph_t *graph);
void CreateCoarseGraph(ctrl_t *ctrl, graph_t *graph, idx_t cnvtxs, 
         idx_t *match);
void CreateCoarseGraphNoMask(ctrl_t *ctrl, graph_t *graph, idx_t cnvtxs, 
         idx_t *match);
graph_t *SetupCoarseGraph(graph_t *graph, idx_t cnvtxs, idx_t dovsize);
void ReAdjustMemory(ctrl_t *ctrl, graph_t *graph, graph_t *cgraph);



/* compress.c */
graph_t *CompressGraph(ctrl_t *ctrl, idx_t nvtxs, idx_t *xadj, idx_t *adjncy, 
             idx_t *vwgt, idx_t *cptr, idx_t *cind);
graph_t *PruneGraph(ctrl_t *ctrl, idx_t nvtxs, idx_t *xadj, idx_t *adjncy, 
             idx_t *vwgt, idx_t *iperm, real_t factor);


/* contig.c */
idx_t FindSepInducedComponents(ctrl_t *, graph_t *, idx_t *, idx_t *);


/* fm.c */
void FM_2WayRefine(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, idx_t niter);
void FM_2WayCutRefine(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, idx_t niter);
void FM_Mc2WayCutRefine(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, idx_t niter);
void SelectQueue(graph_t *graph, real_t *pijbm, real_t *ubfactors, rpq_t **queues, 
         idx_t *from, idx_t *cnum);
void Print2WayRefineStats(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, 
         real_t deltabal, idx_t mincutorder);


/* fortran.c */
void Change2CNumbering(idx_t, idx_t *, idx_t *);
void Change2FNumberingOrder(idx_t, idx_t *, idx_t *, idx_t *, idx_t *);


/* graph.c */
graph_t *SetupGraph(ctrl_t *ctrl, idx_t nvtxs, idx_t ncon, idx_t *xadj, 
             idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt);
void SetupGraph_tvwgt(graph_t *graph);
void SetupGraph_label(graph_t *graph);
graph_t *SetupSplitGraph(graph_t *graph, idx_t snvtxs, idx_t snedges);
graph_t *CreateGraph(void);
void InitGraph(graph_t *graph);
void FreeRData(graph_t *graph);
void FreeGraph(graph_t **graph);
my_graph_t* CreateMyGraph(idx_t n, idx_t m_max);
my_graph_t* CreateMyGraphWspace(ctrl_t* ctrl, idx_t n, idx_t m_max);
void FreeMyGraph(my_graph_t* graph);

/* initpart.c */
void InitSeparator(ctrl_t *ctrl, graph_t *graph, idx_t niparts);
void RandomBisection(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, idx_t niparts);
void GrowBisection(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, idx_t niparts);
void GrowBisectionNode(ctrl_t *ctrl, graph_t *graph, real_t *ntpwgts, idx_t niparts);


/* mcutil.c */
int ivecle(idx_t n, idx_t *x, idx_t *z);
int ivecaxpylez(idx_t n, idx_t a, idx_t *x, idx_t *y, idx_t *z);
int BetterVBalance(idx_t ncon, real_t *itvwgt, idx_t *v_vwgt, idx_t *u1_vwgt,
            idx_t *u2_vwgt);
int BetterBalance2Way(idx_t n, real_t *x, real_t *y);
real_t ComputeLoadImbalance(graph_t *graph, idx_t nparts, real_t *pijbm);
real_t ComputeLoadImbalanceDiff(graph_t *graph, idx_t nparts, real_t *pijbm, 
           real_t *ubvec);
real_t ComputeLoadImbalanceDiffVec(graph_t *graph, idx_t nparts, real_t *pijbm, 
         real_t *ubfactors, real_t *diffvec);


/* mmd.c */
void genmmd(idx_t, idx_t *, idx_t *, idx_t *, idx_t *, idx_t , idx_t *, idx_t *, idx_t *, idx_t *, idx_t, idx_t *);
void mmdelm(idx_t, idx_t *xadj, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t, idx_t);
idx_t mmdint(idx_t, idx_t *xadj, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *);
void mmdnum(idx_t, idx_t *, idx_t *, idx_t *);
void mmdupd(idx_t, idx_t, idx_t *, idx_t *, idx_t, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, idx_t, idx_t *tag);


/* ometis.c */
void MlevelNestedDissection(ctrl_t *ctrl, graph_t *graph, idx_t *order,
         idx_t lastvtx);
void MlevelNestedDissectionCC(ctrl_t *ctrl, graph_t *graph, idx_t *order,
         idx_t lastvtx);
void MlevelNodeBisectionMultiple(ctrl_t *ctrl, graph_t *graph);
void MlevelNodeBisectionL2(ctrl_t *ctrl, graph_t *graph, idx_t niparts);
void MlevelNodeBisectionL1(ctrl_t *ctrl, graph_t *graph, idx_t niparts);
void SplitGraphOrder(ctrl_t *ctrl, graph_t *graph, graph_t **r_lgraph, 
         graph_t **r_rgraph);
graph_t **SplitGraphOrderCC(ctrl_t *ctrl, graph_t *graph, idx_t ncmps,
              idx_t *cptr, idx_t *cind);
void MMDOrder(ctrl_t *ctrl, graph_t *graph, idx_t *order, idx_t lastvtx);
void AMDOrder(ctrl_t* ctrl, graph_t* graph, idx_t* order, idx_t lastvtx);
void RCMOrder(ctrl_t* ctrl, graph_t* graph, idx_t* order, idx_t lastvtx);
void EmptyOrder(ctrl_t* ctrl, graph_t* graph, idx_t* order, idx_t lastvtx);


/* options.c */
ctrl_t *SetupCtrl(moptype_et optype, idx_t *options, idx_t ncon, idx_t nparts, 
            real_t *tpwgts, real_t *ubvec);
void Setup2WayBalMultipliers(ctrl_t *ctrl, graph_t *graph, real_t *tpwgts);
void PrintCtrl(ctrl_t *ctrl);
int CheckParams(ctrl_t *ctrl);
void FreeCtrl(ctrl_t **r_ctrl);


/* refine.c */
void Allocate2WayPartitionMemory(ctrl_t *ctrl, graph_t *graph);
void Compute2WayPartitionParams(ctrl_t *ctrl, graph_t *graph);


/* separator.c */
void ConstructSeparator(ctrl_t *ctrl, graph_t *graph);


/* sfm.c */
void FM_2WayNodeRefine2Sided(ctrl_t *ctrl, graph_t *graph, idx_t niter);
void FM_2WayNodeRefine1Sided(ctrl_t *ctrl, graph_t *graph, idx_t niter);
void FM_2WayNodeBalance(ctrl_t *ctrl, graph_t *graph);


/* srefine.c */
void Refine2WayNode(ctrl_t *ctrl, graph_t *orggraph, graph_t *graph);
void Allocate2WayNodePartitionMemory(ctrl_t *ctrl, graph_t *graph);
void Compute2WayNodePartitionParams(ctrl_t *ctrl, graph_t *graph);
void Project2WayNodePartition(ctrl_t *ctrl, graph_t *graph);


/* timing.c */
void InitTimers(ctrl_t *);
void PrintTimers(ctrl_t *);

/* util.c */
idx_t my_irandInRange_r(size_t range, idx_t * seed);
idx_t my_irandArrayPermute_r(size_t n, idx_t * p, size_t nshuffles,
    idx_t flag, idx_t * seed);
idx_t iargmax_nrm(size_t n, idx_t *x, real_t *y);
idx_t iargmax2_nrm(size_t n, idx_t *x, real_t *y);
void InitRandom(idx_t);
int metis_rcode(int sigrval);



/* wspace.c */
void AllocateWorkSpace(ctrl_t *ctrl, graph_t *graph);
void FreeWorkSpace(ctrl_t *ctrl);
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes);
void wspacepush(ctrl_t *ctrl);
void wspacepop(ctrl_t *ctrl);
idx_t *iwspacemalloc(ctrl_t *, idx_t);
real_t *rwspacemalloc(ctrl_t *, idx_t);
ikv_t *ikvwspacemalloc(ctrl_t *, idx_t);

/* amd.c */
void mmd(idx_t* xadj, idx_t* adjncy, idx_t* perm, idx_t* iperm, my_graph_t* graph, idx_t* degrees, idx_t* mask);
void amd(idx_t* xadj, idx_t* adjncy, idx_t* perm, idx_t* iperm, my_graph_t* graph, my_graph_t* graph_e, my_graph_t* graph_l, my_graph_t* graph_h, my_graph_t* graph_s, idx_t* degrees, idx_t* mask_lp, idx_t* mask_ep, idx_t* mask, idx_t* w);

#endif
