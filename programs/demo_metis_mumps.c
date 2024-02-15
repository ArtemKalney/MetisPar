#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include<dmumps_c.h>
#include<metis.h>
// #include<mpi.h>
#include<omp.h>

// #include "affinity.h"
#include "data_formats.h"

// Some macro magic for time measurement
#define macro_concat(v1,v2) v1 ## v2
#define macro_concat2(v1,v2) macro_concat(v1,v2)
#define macro_var(name) macro_concat2(name,__LINE__)
#define measure_time(dt) \
  struct timespec macro_var(_t0_),macro_var(_t1_); for ( \
  int i##__LINE__ = (clock_gettime(CLOCK_MONOTONIC,&macro_var(_t0_)), 0); \
  !i##__LINE__; \
  i##__LINE__++, \
    clock_gettime(CLOCK_MONOTONIC,&macro_var(_t1_)), \
    ((dt)=(macro_var(_t1_).tv_sec-macro_var(_t0_).tv_sec)+1e-9*(macro_var(_t1_).tv_nsec-macro_var(_t0_).tv_nsec)) \
  )

//---------------------------------------------------------------------------------
void change_affinity() {
  int rank;
//   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank == 0) printf("Affinity changed: ");
//   if (rank == 0) print_current_cpu_set();
//   set_affinity_all_cores();
  if (rank == 0) printf(" -> ");
//   if (rank == 0) print_current_cpu_set();
  if (rank == 0) printf("\n");
}
//---------------------------------------------------------------------------------
double solution(int i,int n) {
  static const double PI = 3.1415926535897932;
  return sin(2*PI*i/n);
}
//---------------------------------------------------------------------------------
// Distance between vectors
double norm2_delta(int n,const double *x,const double *y) {
  double s = 0.0;
  for (int i=0;i<n;i++) {
    const double d = x[i] - y[i];
    s += d * d;
  }
  return sqrt(s);
}

// y[] += dia[][] * x[]
void coo_matrix_vector_multiply_add(const coo_matrix *coo,const double *x,double *y) {
  const int base = coo->base;
  const int nnz = coo->nnz;
  if (coo->type == MT_SYMMETRIC) {
    for (int i=0;i<nnz;i++) {
      const int row = coo->row[i] - base;
      const int col = coo->col[i] - base;
      const double value = coo->val[i];
      y[row] += value * x[col];
      if (row != col) y[col] += value * x[row];
    }
  }
  else {
    for (int i=0;i<nnz;i++) {
      const int row = coo->row[i] - base;
      const int col = coo->col[i] - base;
      const double value = coo->val[i];
      y[row] += value * x[col];
    }
  }
}
//---------------------------------------------------------------------------------
// METIS usage

char *metis_error_message(int code) {
  switch (code) {
    case METIS_OK: return "success";
    case METIS_ERROR_INPUT: return "input error";
    case METIS_ERROR_MEMORY: return "memory allocation error";
    case METIS_ERROR: return "other error";
  }
  return "unknown error";
}

void run_metis_ND(const csr_matrix *csr,int *perm,int *iperm) {
  assert(sizeof(idx_t) == sizeof(int));

  int n = csr->n;
  int *xadj = csr->row;
  int *adjncy = csr->col;

  int options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NUMBERING] = csr->base;
    // 0 - C-style numbering is assumed that starts from 0.
    // 1 - Fortran-style numbering is assumed that starts from 1.

  options[METIS_OPTION_DBGLVL] = 0
    //| METIS_DBG_INFO       // (1) Prints various diagnostic messages.
    //| METIS_DBG_TIME       // (2) Performs timing analysis.
    //| METIS_DBG_COARSEN    // (4) Displays various statistics during coarsening.
    //| METIS_DBG_REFINE     // (8) Displays various statistics during refinement.
    //| METIS_DBG_IPART      // (16) Displays various statistics during initial partitioning.
    //| METIS_DBG_MOVEINFO   // (32) Displays detailed information about vertex moves during refinement.
    //| METIS_DBG_SEPINFO    // (64) Displays information about vertex separators.
    //| METIS_DBG_CONNINFO   // (128) Displays information related to the minimization of subdomain connectivity.
    //| METIS_DBG_CONTIGINFO // (256) Displays information related to the elimination of connected components.
    ;

  int ret = METIS_NodeND(&n, xadj, adjncy, NULL, options, perm, iperm);
  printf("METIS_NodeND: %s\n",metis_error_message(ret));
}
//---------------------------------------------------------------------------------
// METIS usage

char *mumps_reorder_method(int method) {
  // see icntl(7) description in MUMPS userguide
  switch(method) {
    case 0: return "AMD"; // Approximate Minimum Degree (AMD)
    case 1: return "user"; // User defined
    case 2: return "AMF"; // Approximate Minimum Fill (AMF)
    case 3: return "SCOTCH";
    case 4: return "PORD";
    case 5: return "METIS";
    case 6: return "QAMD"; // Approximate Minimum Degree with automatic quasi-dense row detection (QAMD)
    case 7: return "auto"; // Automatic choice by the software during analysis phase
  }
  return "unknown";
}

// coo - (in) matrix in coordinate format
// b - (in/out) in: right-hand side vector, out: solution vector
// iperm - (in) reordering vector / NULL
void coo_matrix_solve_mumps(const coo_matrix *coo,double *b,int *iperm) {

  int rank;
//   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  const int JOB_INIT = -1;
  const int JOB_END = -2;
  const int USE_COMM_WORLD = -987654;

  DMUMPS_STRUC_C id = {
    .job = JOB_INIT,
    .par = 1,
    .sym = coo->type == MT_SYMMETRIC ? 2 : 0, // matrix symmetry
    .comm_fortran = USE_COMM_WORLD,
  };

  dmumps_c(&id);

  if (rank == 0) {
    assert(coo->base == 1); // MUMPS requirement
    printf("Job %2d. MUMPS initialization: ",id.job);
    if (id.infog[0] == 0) {
      if (iperm == NULL) printf("Done. Reordering: auto\n");
      else printf("Done. Reordering: user defined\n");
    }
    else {
      printf("Failed. infog(1:2): %d %d\n",id.infog[0],id.infog[1]);
      return;
    }
    id.n = coo->n;
    id.nnz = coo->nnz;
    id.irn = coo->row;
    id.jcn = coo->col;
    id.a = coo->val;
    id.rhs = b;
  }
  id.icntl[ 0] = 6; // output stream for errors (0 - none, 6 - standatd output)
  id.icntl[ 1] = 6; // output stream for diagnostics (0 - none, 6 - standatd output)
  id.icntl[ 2] = 6; // output stream for global info (0 - none, 6 - standatd output)
  id.icntl[ 3] = 0; // message level (<=0 - none, 1 - errors only, 2 - errors,warnings,main statstics, )
  id.icntl[15] = 1; // set number of OpenMP threads to 1
  id.icntl[ 4] = 0; // Matrix in assembled format
  id.icntl[17] = 0; // Matrix is centralized on the host
  id.icntl[19] = 0; // right-hand side is dense and centralized on the host
  id.icntl[20] = 0; // solution vector is gathered on the host
  id.icntl[27] = 1; // use sequential analysis (for user-defined reordering)
  id.icntl[10] =-2; // iterative refinement: 0 - no, -n - perform exactly n steps, n - perform maximum n steps
  id.icntl[13] = 100; // controls the percentage increase in the estimated working space

  if (iperm != NULL) {
    id.icntl[ 6] = 1; // reordering: 7 - automatic choice, 1 - user-defined
    id.perm_in = iperm; // vector with user-defined reordering
  }
  else
    id.icntl[ 6] = 7; // reordering: 7 - automatic choice, 1 - user-defined

  id.job = 1;
  dmumps_c(&id);
  if (rank == 0) {
    printf("Job %2d. Analysis: ",id.job);
    if (id.infog[0] == 0) {
      unsigned long long nnz = id.infog[2] > 0 ? id.infog[2] : -id.infog[2] * 1000000;
      printf("Done. Reordering: %s  Factor NNZ (est.): %llu  Tree front: %d  Tree nodes: %d\n",
             mumps_reorder_method(id.infog[6]),nnz,id.infog[4],id.infog[5]);
    }
    else printf("Failed. infog(1:2): %d %d",id.infog[0],id.infog[1]);  
  }

  id.job = 2;
  dmumps_c(&id);
  if (rank == 0) {
    printf("Job %2d. Factorization: ",id.job);
    if (id.infog[0] == 0) {
      unsigned long long nnz = id.infog[28] > 0 ? id.infog[28] : -id.infog[28] * 1000000;
      printf("Done. Factor NNZ: %llu\n",nnz);
    }
    else printf("Failed. infog(1:2): %d %d\n",id.infog[0],id.infog[1]);
  }

  id.job = 3;
  dmumps_c(&id);
  if (rank == 0) {
    printf("Job %2d. Computing solution: ",id.job);
    if (id.infog[0] == 0) printf("Done.\n");
    else printf("Failed. infog(1:2): %d %d\n",id.infog[0],id.infog[1]);
  }

  id.job = JOB_END;
  dmumps_c(&id);
  if (rank == 0) {
    printf("Job %2d. MUMPS finalization: ",id.job);
    if (id.infog[0] == 0) printf("Done.\n");
    else printf("Failed. infog(1:2): %d %d\n",id.infog[0],id.infog[1]);
  }
}
//---------------------------------------------------------------------------------

int main(int argc,char *argv[]) {
  int rank,np;
//   MPI_Init(&argc,&argv);
//   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
//   MPI_Comm_size(MPI_COMM_WORLD,&np);

  int stop = 0;
  if (rank == 0 && argc < 2) { 
    printf("Usage: %s <filename.mtx>\n",argv[0]);
    stop = 1;
  }
//   MPI_Bcast(&stop,1,MPI_INT,0,MPI_COMM_WORLD);
//   if (stop) { MPI_Finalize(); return 0; }

  //change_affinity();
//   set_affinity_all_cores();

  coo_matrix coo = rank == 0 ? coo_matrix_from_mtx(argv[1]) : coo_matrix_empty();
  csr_matrix csr_metis = rank == 0 ? convert_coo_matrix_to_csr_metis_matrix(&coo) : csr_matrix_empty();
  const int n = coo.n;

  double *x = NULL;
  double *b = NULL;
  int *perm = NULL;
  int *iperm = NULL;

  double delta_time_metis,delta_time_mumps;
  if (rank == 0) {
    //printf("N: %d   NNZ: %d\n",n,coo.nnz);
    x = malloc(sizeof(double[n])); check_null(x,sizeof(double[n]));
    for (int i=0;i<n;i++) x[i] = solution(i,n);
    b = calloc(n,sizeof(double)); check_null(b,sizeof(double[n]));
    coo_matrix_vector_multiply_add(&coo,x,b);
    perm = malloc(sizeof(int[n])); check_null(perm,sizeof(int[n]));
    iperm = malloc(sizeof(int[n])); check_null(iperm,sizeof(int[n]));

    measure_time(delta_time_metis) {
      run_metis_ND(&csr_metis,perm,iperm);
    }
  }
//   MPI_Barrier(MPI_COMM_WORLD);
  measure_time(delta_time_mumps) {
    coo_matrix_solve_mumps(&coo,b,iperm);
  }

  if (rank == 0) {
    double delta = norm2_delta(n,x,b);
    printf("Cores: %d   %s   N: %d   NNZ: %d   Time METIS: %lf   Time MUMPS: %lf   Residual: %le\n",np,
           coo.type == MT_GENERAL ? "Gereral" : coo.type == MT_SYMMETRIC ? "Symmetric" : "???",
           n,coo.nnz,delta_time_metis,delta_time_mumps,delta);
    //vector_save_bin("x",n,b);
    coo_matrix_free(&coo);
    csr_matrix_free(&csr_metis);
    free(x);
    free(b);
    free(perm);
    free(iperm);
  }

//   MPI_Finalize();
  return 0;
}
