/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * fortran.c
 *
 * This file contains code for the fortran to C interface
 *
 * Started 8/19/97
 * George
 *
 */

#include "metislib.h"


/*************************************************************************/
/*! This function changes the numbering to start from 0 instead of 1 */
/*************************************************************************/
void Change2CNumbering(idx_t nvtxs, idx_t *xadj, idx_t *adjncy)
{
  idx_t i;

  for (i=0; i<=nvtxs; i++)
    xadj[i]--;

  for (i=0; i<xadj[nvtxs]; i++)
    adjncy[i]--;
}



/*************************************************************************/
/*! This function changes the numbering to start from 1 instead of 0 */
/*************************************************************************/
void Change2FNumberingOrder(idx_t nvtxs, idx_t *xadj, idx_t *adjncy, 
         idx_t *v1, idx_t *v2)
{
  idx_t i, nedges;

  for (i=0; i<nvtxs; i++) {
    v1[i]++;
    v2[i]++;
  }

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}