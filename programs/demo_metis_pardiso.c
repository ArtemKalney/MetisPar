#include<assert.h>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<metis.h>
#include<mkl_pardiso.h>
#include <time.h>

#include "data_formats.h"

//---------------------------------------------------------------------------------
// METIS usage

void run_metis_ND(const csr_matrix *csr,int *perm,int *iperm) {
  assert(sizeof(idx_t) == sizeof(int));

  int base = csr->base;
  int n = csr->n;
  int *xadj = csr->row;
  int *adjncy = csr->col;

  int options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NUMBERING] = base;
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

  options[METIS_OPTION_LLEVEL] = METIS_MMD;

  int ret = METIS_NodeND(&n, xadj, adjncy, NULL, options, perm, iperm);

  switch (ret) {
    case METIS_OK: printf("METIS_NodeND: Success\n"); break;
    case METIS_ERROR_INPUT: printf("METIS_NodeND: Input error!\n"); break;
    case METIS_ERROR_MEMORY: printf("METIS_NodeND: Memory allocation error!\n"); break;
    case METIS_ERROR: printf("METIS_NodeND: Other error!\n"); break;
    default: printf("METIS_NodeND: Unknown error!\n"); break;
  }
}

//---------------------------------------------------------------------------------
char *pardiso_error_message(int error) {
  switch (error) {
    case   0: return "no error";
    case  -1: return "input inconsistent";
    case  -2: return "no enough memory";
    case  -3: return "ordering problem";
    case  -4: return "zero pivot, numerical factorization or iterative refinement problem";
    case  -5: return "unclassified (internal) error";
    case  -6: return "reordering failed (matrix types 11 and 13 only)";
    case  -7: return "diagonal matrix is singular";
    case  -8: return "32-bit integer overflow problem";
    case  -9: return "not enough memory for OOC";
    case -10: return "error opening OOc files";
    case -11: return "read/write error with OOC files";
    case -12: return "pardiso_64 called from 32-bit library";
  }
  return "unknown error";
}

void csr_matrix_solve_pardiso(const csr_matrix* csr, double* b, int* iperm, int* nnz) {
    const int n = csr->n;
    const int base = csr->base;

    // Setup Pardiso control parameters
    void* pt[64] = { 0 };
    int iparm[64] = {
          [0] = 1,    // No solver default parameters
          [1] = 2,    // Fill-in reordering: 0 - minimum degree algorithm, 2 - nested dissection algorithm from METIS
          [5] = 1,    // Write solution to: 0 - vector x[], 1 - vector b[]
          [7] = 2,    // Max numbers of iterative refinement steps: 0 - default(2), >0 - user defined, <0 - user defined with extended precision
          [9] = csr->type == MT_SYMMETRIC ? 8 : 13,  // Perturb the pivot elements with eps: 8 - 1e-8 (symmetric), 13 - 1e-13 (non-symmetric)
          [10] = csr->type == MT_SYMMETRIC ? 0 : 1,   // Use nonsymmetric permutation and scaling MPS (symmetric indefinite - 0-off, non-symmetric - 1-on)
          [11] = 0,   // Non-transposed matrix: 0 - A x = b, 1 - A^H x  = b, 2 - A^T x = b
          [12] = csr->type == MT_SYMMETRIC ? 0 : 1,   // Maximum weighted matching algorithm (symmetric - 0-off, non-symmetric - 1-on)
          [17] = -1,  // Output: Number of nonzeros in the factor LU: <0 - enable, >=0 - disable
          [18] = 0,   // Output: Mflops for LU factorization: <0 - enable, >=0 - disable
          [26] = 0,   // Check matrix representation: 0 - no check, 1 - check
          [27] = 0,   // Precision: 0 - double, 1 - single
          [34] = 1 - base, // Array indexing: 0 - one-based (Fortran-style), 1 - zero-based (C-style)
          [36] = 0,   // 0 - CSR format, -2 - convert to VBSR format
    };

    int* perm;

    if (iperm != NULL) {
        iparm[4] = 1; // Permutation: 0 - auto, 1 - user (requires [30]=0 and [35]=0), 2 - auto & return it in perm[]
        perm = iperm;
    }
    else {
        perm = malloc(sizeof(int[n])); check_null(perm, sizeof(int[n]));
        iparm[4] = 2; // Permutation: 0 - auto, 1 - user (requires [30]=0 and [35]=0), 2 - auto & return it in perm[]
    }

    int maxfct = 1;  // Maximum number of numerical factorizations
    int mnum = 1;    // Factorization number to use
    int mtype = csr->type == MT_SYMMETRIC ? -2 : 11;  // Matrix type: 11 - real unsymmetric / -2 - real symmetric indefinite
    int nrhs = 1;    // Number of right hand sides
    int msglvl = 0;  // Print statistical information: 0 - no, 1 - yes
    int error = 0;   // Initialize error flag

    int phase = 11;  // Perform phase 1:Permute
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, csr->val, csr->row, csr->col, perm, &nrhs, iparm, &msglvl, NULL, NULL, &error);
    if (error != 0)
    {
        printf("PARDISO: Error during symbolic factorization: %s\n", pardiso_error_message(error));
    }
    else
    {
        *nnz = iparm[17];
    }

    return;

    printf("1. Reordering done. Memory (KB) peak: %d permanent: %d internals: %d   Nonzeros: %d   MFLOPS for factorization: %d\n",
        iparm[14], iparm[15], iparm[16], iparm[17], iparm[18]);

    vector_int_save_txt("perm.txt", n, perm);

    phase = 22;  // Perform phase 2:Factorize
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, csr->val, csr->row, csr->col, perm, &nrhs, iparm, &msglvl, NULL, NULL, &error);
    if (error != 0) { printf("PARDISO: Error during numerical factorization: %s\n", pardiso_error_message(error)); }
    printf("2. Factorization done. Perturbed pivots: %d\n", iparm[13]);

    double* x = malloc(sizeof(double[n]));
    phase = 33;  // Perform phase 3:Solve,Refine solution
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, csr->val, csr->row, csr->col, perm, &nrhs, iparm, &msglvl, b, x, &error);
    if (error != 0) { printf("PARDISO: Error during solution: %s\n", pardiso_error_message(error)); }
    printf("3. Solution done. Iterative refinement steps: %d\n", iparm[6]);

    phase = -1;  // Release internal memory
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, NULL, csr->row, csr->col, NULL, &nrhs, iparm, &msglvl, NULL, NULL, &error);
    if (error != 0) { printf("PARDISO: Error during memory release: %s\n", pardiso_error_message(error)); }

    if (iperm == NULL) free(perm);
    free(x);
}

//---------------------------------------------------------------------------------

void coo_matrix_vector_multiply_add(const coo_matrix* coo, const double* x, double* y) {
    const int base = coo->base;
    const int nnz = coo->nnz;
    if (coo->type == MT_SYMMETRIC) {
        for (int i = 0;i < nnz;i++) {
            const int row = coo->row[i] - base;
            const int col = coo->col[i] - base;
            const double value = coo->val[i];
            y[row] += value * x[col];
            if (row != col) y[col] += value * x[row];
        }
    }
    else {
        for (int i = 0;i < nnz;i++) {
            const int row = coo->row[i] - base;
            const int col = coo->col[i] - base;
            const double value = coo->val[i];
            y[row] += value * x[col];
        }
    }
}

double solution(int i, int n) {
    static const double PI = 3.1415926535897932;
    return sin(2 * PI * i / n);
}

int main(int argc, char* argv[]) {

    coo_matrix coo = coo_matrix_from_mtx(argv[1]);
    const int n = coo.n;

    int* perm = malloc(sizeof(int[n])); check_null(perm, sizeof(int[n]));
    int* iperm = malloc(sizeof(int[n])); check_null(iperm, sizeof(int[n]));

    csr_matrix csr_metis = convert_coo_matrix_to_csr_metis_matrix(&coo);

    clock_t start = clock();
    run_metis_ND(&csr_metis, perm, iperm);
    clock_t end = clock();
    double delta_time_metis = end - start;

    double* x = NULL;
    double* b = NULL;

    /*x = malloc(sizeof(double[n])); check_null(x, sizeof(double[n]));
    for (int i = 0;i < n;i++) x[i] = solution(i, n);
    b = calloc(n, sizeof(double)); check_null(b, sizeof(double[n]));
    coo_matrix_vector_multiply_add(&coo, x, b);*/

    csr_matrix csr;
    if (coo.type == MT_SYMMETRIC) {
        csr = csr_matrix_upper_from_coo_matrix(&coo);
    }
    else {
        csr = csr_matrix_from_coo_matrix(&coo);
    }
    int nnz = 0;

    start = clock();
    csr_matrix_solve_pardiso(&csr, b, perm, &nnz);
    end = clock();
    double delta_time_pardiso = end - start;

    FILE* f = fopen("out.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");

        return EXIT_FAILURE;
    }

    fprintf(f, "%s %f %u\n", argv[1], (float)delta_time_metis / CLOCKS_PER_SEC, nnz);
    //vector_save_bin("x", n, b);

    free(perm);
    free(iperm);
    free(b);
    csr_matrix_free(&csr);
    coo_matrix_free(&coo);
    csr_matrix_free(&csr_metis);
    fclose(f);

    return 0;
}