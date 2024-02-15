#pragma once

#define check_null(ptr,size) \
  if (ptr == NULL) { printf("Error in %s(): error allocating %ld bytes for %s\n",__func__,(size_t)(size),(#ptr)); exit(1); }

//---------------------------------------------------------------------------------
// Data types

enum matrix_type { MT_GENERAL = 0, MT_SYMMETRIC = 1 };

// Diagoinal matrix
typedef struct {
	int n;         // matrix size
	int ndiag;     // number of stored diagonals
	double* val;   // diagonal values: val[ndiag][n]
	int* dis;      // diagonals displacements: [ndiags]
	int type;      // matrix_type
} dia_matrix;

// CSR matrix
typedef struct {
	int n;         // matrix size
	int base;      // 0 - C indexing / 1 - Fortran indexing
	int* row;      // row starting indices: [n+1]
	int* col;      // col indices: [row[n]]
	double* val;   // elements values: val[row[n]]
	int type;      // matrix_type
} csr_matrix;

// Coordinate matrix
typedef struct {
	int n;         // matrix size
	int nnz;       // number of stored elements
	int base;      // 0 - C indexing / 1 - Fortran indexing
	int type;      // matrix_type
	int* row;      // row numbers
	int* col;      // column numbers
	double* val;   // values
} coo_matrix;

//---------------------------------------------------------------------------------
// Vectors

void vector_save_bin(const char* filename, int n, const double* x);
void vector_save_txt(const char* filename, int n, const double* x);
void vector_int_save_txt(const char* filename, int n, const int* x);

//---------------------------------------------------------------------------------
// Diagoinal matrix subroutines

dia_matrix dia_matrix_empty();
dia_matrix dia_matrix_allocate(int n, int ndiag, int dis[ndiag]);
void dia_matrix_free(dia_matrix* dia);
void dia_matrix_save_txt(const char* filename, const dia_matrix* dia);
void dia_matrix_save_as_dia_txt(const char* filename, const dia_matrix* dia);

//---------------------------------------------------------------------------------
// CSR matrix subroutines

csr_matrix csr_matrix_empty();
csr_matrix csr_matrix_from_dia_matrix(const dia_matrix* dia, int base);
csr_matrix csr_matrix_from_coo_matrix(const coo_matrix* coo);
csr_matrix csr_matrix_upper_from_coo_matrix(const coo_matrix* coo);
void csr_matrix_free(csr_matrix* csr);
void csr_matrix_save_txt(const char* filename, const csr_matrix* csr);
void csr_matrix_save_as_csr_txt(const char* filename, const csr_matrix* csr);

//---------------------------------------------------------------------------------
// Coordinate matrix subroutines
coo_matrix coo_matrix_empty();
coo_matrix coo_matrix_from_dia_matrix(const dia_matrix* dia, int base);
coo_matrix coo_matrix_from_mtx(char* filename);
void coo_matrix_free(coo_matrix* coo);
void coo_matrix_save_txt(const char* filename, const coo_matrix* coo);

//---------------------------------------------------------------------------------
// Cenvertion for METIS
csr_matrix convert_coo_matrix_to_csr_metis_matrix(const coo_matrix* coo);

//---------------------------------------------------------------------------------
