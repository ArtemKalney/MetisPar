#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "data_formats.h"
#include "mmio.h"

//---------------------------------------------------------------------------------
// Helper functions

static int merge_sorted_uniq(int n1, int* x1, int n2, int* x2, int maxval, int* y);
static int cmp_int(const void* pa, const void* pb);

//---------------------------------------------------------------------------------
// Vector subroutines

void vector_save_bin(const char* filename, int n, const double* x) {
    FILE* f = fopen(filename, "wb");
    fwrite(x, sizeof(double), n, f);
    fclose(f);
}

void vector_save_txt(const char* filename, int n, const double* x) {
    FILE* f = fopen(filename, "w");
    for (int i = 0;i < n;i++) fprintf(f, "%d %lf\n", i, x[i]);
    fclose(f);
}

void vector_int_save_txt(const char* filename, int n, const int* x) {
    FILE* f = fopen(filename, "w");
    for (int i = 0;i < n;i++) fprintf(f, "%d %d\n", i, x[i]);
    fclose(f);
}

//---------------------------------------------------------------------------------
// Diagoinal matrix routines

dia_matrix dia_matrix_empty() {
    dia_matrix dia = {
      .n = 0,
      .ndiag = 0,
      .val = NULL,
      .dis = NULL,
      .type = MT_GENERAL,
    };
    return dia;
}

dia_matrix dia_matrix_allocate(int n, int ndiag, int dis[ndiag]) {
    double* dia_val = calloc(ndiag * n, sizeof(double)); check_null(dia_val, ndiag * n * sizeof(double));
    int* dia_dis = malloc(sizeof(int[ndiag])); check_null(dia_dis, sizeof(int[ndiag]));
    dia_matrix dia = {
      .n = n,
      .ndiag = ndiag,
      .val = dia_val,
      .dis = dia_dis,
      .type = MT_GENERAL,
    };
    memcpy(dia.dis, dis, ndiag * sizeof(int));
    return dia;
}

void dia_matrix_free(dia_matrix* dia) {
    dia->n = 0;
    dia->ndiag = 0;
    if (dia->val) { free(dia->val); dia->val = NULL; }
    if (dia->dis) { free(dia->dis); dia->dis = NULL; }
}

void dia_matrix_save_txt(const char* filename, const dia_matrix* dia) {
    const int n = dia->n;
    const int ndiag = dia->ndiag;
    const double(*dia_val)[ndiag][n] = (const double(*)[ndiag][n]) dia->val;
    const int(*dia_dis)[ndiag] = (const int(*)[ndiag]) dia->dis;
    FILE* f = fopen(filename, "w");
    for (int i = 0;i < n;i++)
        for (int j = 0;j < ndiag;j++) {
            const int dis = (*dia_dis)[j];
            if (i + dis >= 0 && (*dia_val)[j][i] != 0.0)
                fprintf(f, "%d %d %lf\n", i, i + dis, (*dia_val)[j][i]);
        }
    fclose(f);
}

void dia_matrix_save_as_dia_txt(const char* filename, const dia_matrix* dia) {
    const int n = dia->n;
    const int ndiag = dia->ndiag;
    const double(*dia_val)[ndiag][n] = (const double(*)[ndiag][n]) dia->val;
    const int(*dia_dis)[ndiag] = (const int(*)[ndiag]) dia->dis;
    FILE* f = fopen(filename, "w");
    fprintf(f, "%d %d\n", ndiag, n);
    for (int i = 0;i < n;i++)
        for (int j = 0;j < ndiag;j++) {
            const int dis = (*dia_dis)[j];
            if (i + dis >= 0 && (*dia_val)[j][i] != 0.0)
                fprintf(f, "%d %d %lf\n", i, i + dis, (*dia_val)[j][i]);
        }
    fclose(f);
}

//---------------------------------------------------------------------------------
// CSR matrix routines

csr_matrix csr_matrix_empty() {
    csr_matrix csr = {
      .n = 0,
      .base = 0,
      .row = NULL,
      .col = NULL,
      .val = NULL,
      .type = MT_GENERAL,
    };
    return csr;
}

csr_matrix csr_matrix_from_dia_matrix(const dia_matrix* dia, int base) {
    const int n = dia->n;
    const int ndiag = dia->ndiag;
    const double(*dia_val)[ndiag][n] = (const double(*)[ndiag][n]) dia->val;

    int k = 0;
    for (int j = 0;j < ndiag;j++) {
        const int dis = dia->dis[j];
        const int i0 = dis < 0 ? -dis : 0;
        const int i1 = dis < 0 ? n : n - dis;
        for (int i = i0;i < i1;i++) {
            const double value = (*dia_val)[j][i];
            if (value != 0.0) k++;
        }
    }
    const int nnz = k;

    int* csr_row = malloc(sizeof(int[n + 1])); check_null(csr_row, sizeof(int[n + 1]));
    int* csr_col = malloc(sizeof(int[nnz])); check_null(csr_col, sizeof(int[nnz]));
    double* csr_val = malloc(sizeof(double[nnz])); check_null(csr_val, sizeof(double[nnz]));

    csr_matrix csr = {
      .n = n,
      .base = base,
      .row = csr_row,
      .col = csr_col,
      .val = csr_val,
      .type = dia->type,
    };

    csr.row[0] = base;
    k = 0;
    for (int i = 0;i < n;i++) {
        for (int j = 0;j < ndiag;j++) {
            const int col = dia->dis[j] + i;
            if (col >= 0 && col < n) {
                const double value = (*dia_val)[j][i];
                if (value != 0.0) {
                    csr.col[k] = col + base;
                    csr.val[k] = value;
                    k++;
                }
            }
        }
        csr.row[i + 1] = base + k;
    }
    return csr;
}

csr_matrix csr_matrix_from_coo_matrix(const coo_matrix* coo) {
    const int base = coo->base;
    const int n = coo->n;
    const int nnz = coo->nnz;

    int* rs = calloc(n + 1, sizeof(int)); check_null(rs, sizeof(int[n]));
    int* row = malloc(sizeof(int[n + 1])); check_null(row, sizeof(int[n + 1]));
    int* col = malloc(sizeof(int[nnz])); check_null(col, sizeof(int[nnz]));
    double* val = malloc(sizeof(double[nnz])); check_null(val, sizeof(double[nnz]));
    char* need_sort = calloc(n, sizeof(char)); check_null(need_sort, sizeof(char[n]));

    { int prev_r = -1, prev_c = -1;
    for (int i = 0;i < nnz;i++) {
        const int r = coo->row[i] - base;
        const int c = coo->col[i] - base;
        rs[r]++;
        if (prev_r == r && prev_c > c) need_sort[r] = 1;
        prev_r = r;
        prev_c = c;
    }
    }

    { int k = base;
    for (int i = 0;i < n;i++) {
        const int current = k;
        k += rs[i];
        rs[i] = current;
        row[i] = current;
    }
    row[n] = k;
    }

    assert(row[n] == nnz + base);

    for (int i = 0;i < nnz;i++) {
        const int r = coo->row[i] - base;
        const int idxr = rs[r];
        col[idxr - base] = coo->col[i];
        val[idxr - base] = coo->val[i];
        rs[r] = idxr + 1;
    }

    for (int i = 0;i < n;i++) {
        const int ri = row[i] - base;
        if (need_sort[i]) {
            const int count = row[i + 1] - ri - base;
            int(*ca)[count] = (int(*)[count]) & col[ri];
            double(*va)[count] = (double(*)[count]) & val[ri];
            for (int i = 0;i < count - 1;i++) { // sort
                const int ci = (*ca)[i];
                int min_c = ci;
                int min_idx = i;
                for (int j = i + 1;j < count;j++) { // find min
                    const int cj = (*ca)[j];
                    if (min_c > cj) { min_c = cj; min_idx = j; }
                }
                if (min_idx != i) {
                    (*ca)[i] = min_c;
                    (*ca)[min_idx] = ci;
                    const double vi = (*va)[i];
                    (*va)[i] = (*va)[min_idx];
                    (*va)[min_idx] = vi;
                }
            }
        }
    }
    /*
    for (int r=0;r<n;r++) {
      int prev = -1;
      for (int j=row[r];j<row[r+1];++j) {
        int c = col[j-base];
        if (c <= prev) printf("row[%d] %d : col[%d] %d - [%d] %d\n",r,r+base,j-1,prev,j,c);
        prev = c;
      }
    }
    */
    free(rs);
    free(need_sort);

    csr_matrix csr = {
      .base = base,
      .n = n,
      .row = row,
      .col = col,
      .val = val,
      .type = coo->type,
    };
    return csr;
}

csr_matrix csr_matrix_upper_from_coo_matrix(const coo_matrix* coo) {
    const int base = coo->base;
    const int n = coo->n;
    const int nnz = coo->nnz;

    int* rs = calloc(n + 1, sizeof(int)); check_null(rs, sizeof(int[n]));
    int* row = malloc(sizeof(int[n + 1])); check_null(row, sizeof(int[n + 1]));
    int* col = malloc(sizeof(int[nnz])); check_null(col, sizeof(int[nnz]));
    double* val = malloc(sizeof(double[nnz])); check_null(val, sizeof(double[nnz]));
    char* need_sort = calloc(n, sizeof(char)); check_null(need_sort, sizeof(char[n]));

    { int prev_r = -1, prev_c = -1;
    for (int i = 0;i < nnz;i++) {
        int r = coo->row[i] - base;
        int c = coo->col[i] - base;
        if (r > c) { int tmp = r; r = c; c = tmp; }
        rs[r]++;
        if (prev_r == r && prev_c > c) need_sort[r] = 1;
        prev_r = r;
        prev_c = c;
    }
    }

    { int k = base;
    for (int i = 0;i < n;i++) {
        const int current = k;
        k += rs[i];
        rs[i] = current;
        row[i] = current;
    }
    row[n] = k;
    }

    assert(row[n] == nnz + base);

    for (int i = 0;i < nnz;i++) {
        int r = coo->row[i] - base;
        int c = coo->col[i] - base;
        if (r > c) { int tmp = r; r = c; c = tmp; }
        const int idxr = rs[r];
        col[idxr - base] = c + base;
        val[idxr - base] = coo->val[i];
        rs[r] = idxr + 1;
    }

    for (int i = 0;i < n;i++) {
        const int ri = row[i] - base;
        if (need_sort[i]) {
            const int count = row[i + 1] - ri - base;
            int(*ca)[count] = (int(*)[count]) & col[ri];
            double(*va)[count] = (double(*)[count]) & val[ri];
            for (int i = 0;i < count - 1;i++) { // sort
                const int ci = (*ca)[i];
                int min_c = ci;
                int min_idx = i;
                for (int j = i + 1;j < count;j++) { // find min
                    const int cj = (*ca)[j];
                    if (min_c > cj) { min_c = cj; min_idx = j; }
                }
                if (min_idx != i) {
                    (*ca)[i] = min_c;
                    (*ca)[min_idx] = ci;
                    const double vi = (*va)[i];
                    (*va)[i] = (*va)[min_idx];
                    (*va)[min_idx] = vi;
                }
            }
        }
    }
    /*
    for (int r=0;r<n;r++) {
      int prev = -1;
      for (int j=row[r];j<row[r+1];++j) {
        int c = col[j-base];
        if (c <= prev) printf("row[%d] %d : col[%d] %d - [%d] %d\n",r,r+base,j-1,prev,j,c);
        prev = c;
      }
    }
    */
    free(rs);
    free(need_sort);

    csr_matrix csr = {
      .base = base,
      .n = n,
      .row = row,
      .col = col,
      .val = val,
      .type = coo->type,
    };
    return csr;
}

void csr_matrix_free(csr_matrix* csr) {
    csr->n = 0;
    if (csr->row) { free(csr->row); csr->row = NULL; }
    if (csr->col) { free(csr->col); csr->col = NULL; }
    if (csr->val) { free(csr->val); csr->val = NULL; }
}

void csr_matrix_save_txt(const char* filename, const csr_matrix* csr) {
    const int base = csr->base;
    const int n = csr->n;
    const int nnz = csr->row[n] - base;
    FILE* f = fopen(filename, "w");
    fprintf(f, "%d %d %d\n", n, nnz, base);
    if (csr->val != NULL)
        for (int i = 0;i < n;i++)
            for (int j = csr->row[i] - base;j < csr->row[i + 1] - base;j++)
                fprintf(f, "%d %d %lf\n", i + base, csr->col[j], csr->val[j]);
    else
        for (int i = 0;i < n;i++)
            for (int j = csr->row[i] - base;j < csr->row[i + 1] - base;j++)
                fprintf(f, "%d %d\n", i + base, csr->col[j]);
    fclose(f);
}

void csr_matrix_save_as_csr_txt(const char* filename, const csr_matrix* csr) {
    const int base = csr->base;
    const int n = csr->n;
    const int nnz = csr->row[n] - base;
    FILE* f = fopen(filename, "w");
    fprintf(f, "%d %d %d\n", n, nnz, base);
    for (int i = 0;i < n + 1;i++) fprintf(f, "%d ", csr->row[i]);
    fprintf(f, "\n");
    for (int i = 0;i < nnz;i++) fprintf(f, "%d ", csr->col[i]);
    fprintf(f, "\n");
    if (csr->val != NULL) {
        for (int i = 0;i < nnz;i++) fprintf(f, "%le ", csr->val[i]);
        fprintf(f, "\n");
    }
    fclose(f);
}

//---------------------------------------------------------------------------------
// Coordinate matrix routines

coo_matrix coo_matrix_empty() {
    coo_matrix coo = {
      .n = 0,
      .nnz = 0,
      .base = 0,
      .type = MT_GENERAL,
      .row = NULL,
      .col = NULL,
      .val = NULL,
    };
    return coo;
}

coo_matrix coo_matrix_from_dia_matrix(const dia_matrix* dia, int base) {
    const int n = dia->n;
    const int ndiag = dia->ndiag;
    const double(*dia_val)[ndiag][n] = (const double(*)[ndiag][n]) dia->val;

    int k = 0;
    for (int j = 0;j < ndiag;j++) {
        const int dis = dia->dis[j];
        const int i0 = dis < 0 ? -dis : 0;
        const int i1 = dis < 0 ? n : n - dis;
        for (int i = i0;i < i1;i++) {
            const double value = (*dia_val)[j][i];
            if (value != 0.0) k++;
        }
    }
    const int nnz = k;

    int* coo_row = malloc(sizeof(int[nnz])); check_null(coo_row, sizeof(int[nnz]));
    int* coo_col = malloc(sizeof(int[nnz])); check_null(coo_col, sizeof(int[nnz]));
    double* coo_val = malloc(sizeof(double[nnz])); check_null(coo_val, sizeof(double[nnz]));

    coo_matrix coo = {
      .n = n,
      .nnz = nnz,
      .base = base,
      .type = dia->type,
      .row = coo_row,
      .col = coo_col,
      .val = coo_val,
    };

    k = 0;
    for (int i = 0;i < n;i++)
        for (int j = 0;j < ndiag;j++) {
            const int dis = dia->dis[j];
            if (i + dis >= 0 && i + dis < n) {
                const double value = (*dia_val)[j][i];
                if (value != 0.0) {
                    coo.row[k] = base + i;
                    coo.col[k] = base + i + dis;
                    coo.val[k] = value;
                    k++;
                }
            }
        }

    return coo;
}

char* mm_error(int code) {
    return
        code == MM_COULD_NOT_READ_FILE ? "could not read file" :
        code == MM_PREMATURE_EOF ? "premature end of file" :
        code == MM_NOT_MTX ? "not mtx file" :
        code == MM_NO_HEADER ? "the file does not begin with " MatrixMarketBanner :
        code == MM_UNSUPPORTED_TYPE ? "unsupported type" :
        code == MM_LINE_TOO_LONG ? "line in file is too long" :
        code == MM_COULD_NOT_WRITE_FILE ? "could not write file" :
        "unknown error";
}

coo_matrix coo_matrix_from_mtx(char* filename) {
    FILE* f = fopen(filename, "r");
    if (f == NULL) {
        printf("%s : Can't open file '%s' : ", __func__, filename); perror("");
        exit(1);
    }
    MM_typecode matcode;
    int ret = mm_read_banner(f, &matcode);
    if (ret != 0) {
        printf("%s : Can't read Matrix Market header in '%s' : %s\n", __func__, filename, mm_error(ret));
        exit(1);
    }
    if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) || !mm_is_coordinate(matcode)) {
        printf("%s : Not a sparse matrix in coordinate format in file '%s'\n", __func__, filename);
        exit(1);
    }
    if (!mm_is_real(matcode)) {
        printf("%s : Not a matrix with real values in '%s'\n", __func__, filename);
        exit(1);
    }
    int m, n, nnz;
    ret = mm_read_mtx_crd_size(f, &m, &n, &nnz);
    if (ret != 0) {
        printf("%s : Can't read matrix sizes in '%s' : %s\n", __func__, filename, mm_error(ret));
        exit(1);
    }
    if (m != n) {
        printf("%s : Not a squere matrix in '%s'\n", __func__, filename);
        exit(1);
    }
    int* row = malloc(sizeof(int[nnz])); check_null(row, sizeof(int[nnz]));
    int* col = malloc(sizeof(int[nnz])); check_null(col, sizeof(int[nnz]));
    double* val = malloc(sizeof(double[nnz])); check_null(val, sizeof(double[nnz]));

    for (int i = 0;i < nnz;i++) {
        if (fscanf(f, "%d %d %lg\n", &row[i], &col[i], &val[i]) != 3) {
            printf("%s : Unexpected EOF in file '%s' on number %d\n", __func__, filename, i);
            exit(1);
        }
    }
    fclose(f);

    coo_matrix coo = {
      .base = 1,
      .n = n,
      .nnz = nnz,
      .type = mm_is_symmetric(matcode) ? MT_SYMMETRIC : MT_GENERAL,
      .row = row,
      .col = col,
      .val = val,
    };
    return coo;
}

void coo_matrix_free(coo_matrix* coo) {
    coo->n = 0;
    coo->nnz = 0;
    if (coo->row) { free(coo->row); coo->row = NULL; }
    if (coo->col) { free(coo->col); coo->col = NULL; }
    if (coo->val) { free(coo->val); coo->val = NULL; }
}

void coo_matrix_save_txt(const char* filename, const coo_matrix* coo) {
    const int nnz = coo->nnz;
    FILE* f = fopen(filename, "w");
    fprintf(f, "%d %d %d\n", coo->n, nnz, coo->base);
    for (int i = 0;i < nnz;i++)
        fprintf(f, "%d %d %lf\n", coo->row[i], coo->col[i], coo->val[i]);
    fclose(f);
}

//---------------------------------------------------------------------------------
// Cenvertion for METIS

// helper function: merge sorted lists
static int merge_sorted_uniq(int n1, int* x1, int n2, int* x2, int maxval, int* y) {
    int i = 0;
    int j = 0;
    int v1 = x1[i];
    int v2 = x2[j];
    int go1 = i < n1&& v1 < maxval;
    int go2 = j < n2&& v2 < maxval;
    int k = 0;
    while (go1 || go2) {
        while (i < n1 && v1 < maxval && (!go2 || v1 <= v2)) {
            y[k++] = v1;
            if (go2 && v1 == v2) {
                j++;
                if (j < n2) v2 = x2[j];
            }
            i++;
            if (i < n1) v1 = x1[i];
        }
        go1 = i < n1&& v1 < maxval;
        while (j < n2 && v2 < maxval && (!go1 || v2 <= v1)) {
            y[k++] = v2;
            if (go1 && v2 == v1) {
                i++;
                if (i < n1) v1 = x1[i];
            }
            j++;
            if (j < n2) v2 = x2[j];
        }
        go2 = j < n2&& v2 < maxval;
    }
    return k;
}

// helper function: compare ints
static int cmp_int(const void* pa, const void* pb) {
    const int a = *(const int*)pa;
    const int b = *(const int*)pb;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

csr_matrix convert_coo_matrix_to_csr_metis_matrix(const coo_matrix* coo) {
    const int base = coo->base;
    const int n = coo->n;
    const int nnz = coo->nnz;
    const int nnz_max = 2 * nnz;

    int* rs = calloc(n + 1, sizeof(int));
    int* ri = calloc(n + 1, sizeof(int));
    int* rc = malloc(sizeof(int[nnz_max]));

    int count = 0;
    for (int i = 0;i < nnz;i++) {
        const int row = coo->row[i] - base;
        const int col = coo->col[i] - base;
        if (row != col) { rs[row]++; rs[col]++; }
    }

    { int k = 0;
    for (int i = 0;i < n;i++) {
        //printf("rs %3d : %2d -> %2d\n",i,rs[i],k);
        const int current = k;
        k += rs[i];
        rs[i] = current;
    }
    //printf("rs %3d : %2d -> %2d\n",n,rs[n],k);
    rs[n] = k;
    }
    if (rs[n] > nnz_max) printf("%d %d \n", rs[n], nnz_max);
    assert(rs[n] <= nnz_max);

    for (int i = 0;i < nnz;i++) {
        const int row = coo->row[i] - base;
        const int col = coo->col[i] - base;
        if (row != col) {
            const int idxr = ri[row];
            rc[rs[row] + idxr] = col;
            ri[row] = idxr + 1;
            const int idxc = ri[col];
            rc[rs[col] + idxc] = row;
            ri[col] = idxc + 1;
        }
    }
    int nnz_new = 0;
    int* row_new = malloc(sizeof(int[n + 1]));
    row_new[0] = base;
    for (int i = 0;i < n;i++) {
        const int count = ri[i];
        //if (rs[i]+ri[i] > rs[i+1]) printf("%3d ! res: %d cnt: %d\n",i,rs[i+1]-rs[i],count);
        if (count > 1) {
            int* cols = &rc[rs[i]];
            qsort(cols, count, sizeof(int), cmp_int);
            //{ printf("%3d :",i); for (int j=0;j<count;j++) printf(" %d",cols[j]); printf("\n"); }
            int k = 1;
            int prev = cols[0];
            for (int j = 1;j < count;j++) {
                const int next = cols[j];
                if (next != prev) {
                    cols[k] = next;
                    k++;
                    prev = next;
                }
            }
            ri[i] = k;
        }
        //else printf("%3d : %d\n",i,rc[rs[i]]);
        nnz_new += ri[i];
        row_new[i + 1] = row_new[i] + ri[i];
    }
    int* col_new = malloc(sizeof(int[nnz_new]));
    for (int i = 0;i < n;i++) {
        const int* from = &rc[rs[i]];
        int* to = &col_new[row_new[i] - base];
        for (int j = 0;j < ri[i];j++) to[j] = from[j] + base;
    }

    free(rs);
    free(ri);
    free(rc);

    csr_matrix csr_metis = {
      .base = base,
      .n = n,
      .row = row_new,
      .col = col_new,
      .val = NULL,
    };
    return csr_metis;
}

//---------------------------------------------------------------------------------
