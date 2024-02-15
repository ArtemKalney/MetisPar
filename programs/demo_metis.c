#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<metis.h>
#include <time.h>
#include "data_formats.h"
#include <omp.h>

//---------------------------------------------------------------------------------
// METIS usage

void run_metis_ND(const csr_matrix* csr, int* perm, int* iperm) {
    assert(sizeof(idx_t) == sizeof(int));

    int base = csr->base;
    int n = csr->n;
    int* xadj = csr->row;
    int* adjncy = csr->col;

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
		//| METIS_DBG_LOW_LEVEL    // (4096) Displays information related to the low level reordering.
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

int main(int argc, char* argv[]) {
    coo_matrix coo = coo_matrix_from_mtx(argv[1]);
    csr_matrix csr_metis = convert_coo_matrix_to_csr_metis_matrix(&coo);
    const int n = coo.n;

    int* perm = malloc(sizeof(int[n])); check_null(perm, sizeof(int[n]));
    int* iperm = malloc(sizeof(int[n])); check_null(iperm, sizeof(int[n]));

    double start = omp_get_wtime();
	run_metis_ND(&csr_metis, perm, iperm);
	double end = omp_get_wtime();

    FILE* f = fopen("out.txt", "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");

        return EXIT_FAILURE;
    }

    fprintf(f, "%s %f\n", argv[1], end - start);

    coo_matrix_free(&coo);
    csr_matrix_free(&csr_metis);
    free(perm);
    free(iperm);
    fclose(f);

    return EXIT_SUCCESS;
}