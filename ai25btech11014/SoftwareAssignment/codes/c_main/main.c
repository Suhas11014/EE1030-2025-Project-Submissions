#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrixfunc.h"
#include "svd.h"
#include "img.h"

int main(int argc, char *argv[]) {//allows and keeps track of arguments
    int m = 0;
    int n = 0;
    int k = 0;
    
    char *input_file;
    char *output_file;
    char *format;
    
    double **A = NULL;
    double **A_T = NULL;
    double **C = NULL;
    double **V = NULL;
    double **U = NULL;
    double *lambda = NULL;
    double **A_k = NULL;

    if (argc != 5) { 
        fprintf(stderr, "Usage: %s <input_file> <output_file> <rank_k> <format (png/jpg)>\n", argv[0]);
        return 1;
    }//makes sure 5 arguments are given
    
    input_file = argv[1];
    output_file = argv[2];
    k = atoi(argv[3]);
    format = argv[4];

    if (k <= 0) {
        fprintf(stderr, "Error: Rank k must be a positive integer.\n");
        return 1;
    }//rejects k<0

    A = read_img(input_file, &m, &n);
    if (A == NULL) {
        return 1; 
    }

    if (k > m || k > n) {
        int min_dim = (m < n) ? m : n;
        fprintf(stderr, "Warning: Rank k reduced from %d to min dimension %d.\n", k, min_dim);
        k = min_dim;
    }//fallback if k is out of bounds
    

    A_T = transpose(A, m, n);
    if (A_T == NULL) {
        goto cleanup_A;
    }//frees A,m

    C = multiply(A_T, n, m, A, m, n);
    freeMatrix(A_T, n); 
    if (C == NULL){
         goto cleanup_A;
    }//frees A,m

    V = createMatrix(n, n);
    lambda = (double*)calloc(n, sizeof(double)); 
    if (V == NULL || lambda == NULL) {
         goto cleanup_C; 
        }//frees C,n
    
    svd_power_deflation(C, n, k, V, lambda); 
    
    freeMatrix(C, n); 

    U = compute_U(A, V, lambda, m, n);
    if (U == NULL) { 
        goto cleanup_VL; 
    }//frees V,n,lambda
    
    A_k = final_image(U, V, lambda, m, n, k);
    if (A_k == NULL) {
         goto cleanup_UVL;
         }//frees U,m
    
    double error = calc_error(A, A_k, m, n); //calculates frobenius norm
    
    freeMatrix(A, m);

    if (error == -1.0) { 
        fprintf(stderr, "Error: Failed to calculate error due to allocation failure.\n");
        goto cleanup_UVL;
    }//frees U,m
    
    printf("Frobenius Norm Error (||A - A_k||): %.4f\n", error);
    write_img(output_file, A_k, m, n, format);
    
    // below is the code for freeing allocated memory using goto

cleanup_Ak: 
    if (A_k) freeMatrix(A_k, m);
cleanup_UVL: 
    if (U) freeMatrix(U, m);
cleanup_VL: 
    if (V) freeMatrix(V, n);
    if (lambda) free(lambda);
    return 0;

cleanup_C: 
    if (C) freeMatrix(C, n);
cleanup_A_T: 
    if (A_T) freeMatrix(A_T, n);
cleanup_A: 
    if (A) freeMatrix(A, m); 
    return 1;
}