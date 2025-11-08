#ifndef matrixfunc_h
#define matrixfunc_h

double** transpose(double** matrix, int rows, int cols);

double** multiply(double** matrixA, int rowsA, int colsA, double** matrixB, int rowsB, int colsB);

void freeMatrix(double** matrix, int rows);

double** createMatrix(int rows, int cols);

void copyMatrix(double** source, double** destination, int rows, int cols);

void swapcols(double** matrix, int rows, int col1, int col2);

double calc_error(double **A, double **B, int m, int n);

double dot_product(double *a, double *b, int n);

double* multiply_matrix_vector(double **C, int n, double *v);

double vector_norm(double *v, int n);

double square(double x);

#endif /* matrixfunc_h */