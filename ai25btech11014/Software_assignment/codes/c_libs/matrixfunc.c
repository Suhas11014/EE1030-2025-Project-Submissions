#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixfunc.h"

double** transpose(double** matrix, int rows, int cols){
double** result = createMatrix(cols, rows);

if (result == NULL) {
        return NULL; // allocation failed
    }

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            result[j][i] = matrix[i][j];
        }
    }
    return result;

}

double** multiply(double** matrixA, int rowsA, int colsA, double** matrixB, int rowsB, int colsB){

    if(colsA!=rowsB){
        return NULL ;
    }// A,B cannot be multiplied


    double** result = createMatrix(rowsA, colsB);
    for(int i = 0; i < rowsA; i++){
        for(int j = 0; j < colsB; j++){
            result[i][j] = 0;
            for(int k = 0; k < colsA; k++){
                result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
    return result;
}


void freeMatrix(double** matrix, int rows){
    for(int i = 0; i < rows; i++){
        free(matrix[i]);
    }
    free(matrix);
}// frees a dynamically allocated matrix


double** createMatrix(int rows, int cols){
    double** matrix = (double**)calloc(rows, sizeof(double*));
    if (matrix == NULL) {
        return NULL; // allocation failed
    }

    for(int i = 0; i < rows; i++){
        matrix[i] = (double*)calloc(cols, sizeof(double));
        if (matrix[i] == NULL) {
            // allocation failed
            for(int j = 0; j < i; j++){
                free(matrix[j]);
            }
            free(matrix);
            return NULL; // allocation failed
        }
    }
    return matrix;
}


void copyMatrix(double** source, double** destination, int rows, int cols){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            destination[i][j] = source[i][j];
        }
    }
}// pretty self-explanatory


void swapcols(double** matrix, int rows, int col1, int col2){
    for(int i = 0; i < rows; i++){
        double temp = matrix[i][col1];
        matrix[i][col1] = matrix[i][col2];
        matrix[i][col2] = temp;
    }
}// again, pretty self-explanatory


double square(double x) {
    return x * x;
}//seems obvious


double calc_error(double **A, double **B, int m, int n){
    double result ;
    double sum = 0.0;
    double** diff = createMatrix(m, n);
    if (diff == NULL) {
        return -1; // allocation failed
    }

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            diff[i][j] = A[i][j] - B[i][j];
        }
    }

    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            sum += square(diff[i][j]);
        }
    }
    result=sqrt(sum);
    freeMatrix(diff, m);
    return result;
   
}// calculates the Frobenius norm of the difference between the original and approximated img(matrix)

double dot_product(double *a, double *b, int n){
    double result = 0.0;
    for(int i = 0; i < n; i++){
        result += a[i] * b[i];
    }
    return result;
}//again, pretty obvious

double* multiply_matrix_vector(double **C, int n, double *v){
    double* result = (double*)calloc(n, sizeof(double));
    if (result == NULL) {
        return NULL; 
    }

    for(int i = 0; i < n; i++){
        result[i] = 0.0;
        for(int j = 0; j < n; j++){
            result[i] += C[i][j] * v[j];
        }
    }
    return result;
}// so the reason i defined this again is that whaen a vector is involved, the answer is represented as a 1D array, not a 2D array with one column

double vector_norm(double *v, int n){
    double result=dot_product(v, v, n);
    return sqrt(result);
}// calculates the "length"(or Euclidean norm) of a vector
