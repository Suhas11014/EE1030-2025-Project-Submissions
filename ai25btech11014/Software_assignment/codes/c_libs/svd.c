#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixfunc.h"
#include "svd.h"
#define MAX_ITERATIONS 50 //prevents the no of loop from going to infinity
#define LIMIT 1.0e-9 //prevents x/0 type of cases

void power_iteration(double **C_temp, int n, double *v_out, double *lambda_out) {
    double initial_val = 1.0 / sqrt((double)n);
    for (int i = 0; i < n; i++) {
        v_out[i] = initial_val; 
    }

    double *w = NULL;

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        w = multiply_matrix_vector(C_temp, n, v_out); 
        if (w == NULL) {
            fprintf(stderr, "Error: Memory allocation failed in multiply_matrix_vector.\n");
            return; 
        }//check for failed memory allocation

        double norm = vector_norm(w, n); 
        
        if (norm < LIMIT) { 
            free(w);
            break; 
        }//prevents x/0 type of cases

        for (int i = 0; i < n; i++) {
            v_out[i] = w[i] / norm; 
        }

        free(w); 
        w = NULL;
    }

    double *Cv = multiply_matrix_vector(C_temp, n, v_out);
    if (Cv == NULL) return;

    double lambda_up = dot_product(v_out, Cv, n); //finds numerator
    double lambda_down = dot_product(v_out, v_out, n);  //finds denominator
    *lambda_out =  lambda_up/ lambda_down;//divides

    free(Cv);
}

void svd_power_deflation(double **C, int n, int k, double **V, double *lambda) {
    double **C_temp = createMatrix(n, n); 
    if (C_temp == NULL) return;
    copyMatrix(C, C_temp, n, n);
    
    double *v_temp = (double*)malloc(n * sizeof(double)); 
    if (v_temp == NULL) {
        freeMatrix(C_temp, n);
        return; 
    }//check for memory allocation fail

    for (int i = 0; i < k; i++) {
        double current_lambda;
        power_iteration(C_temp, n, v_temp, &current_lambda); 

        if (current_lambda <= LIMIT) { 
            k = i;
            fprintf(stderr, "Warning: Stopping SVD early at rank %d due to negligible eigenvalue.\n", i);
            break; 
        }//prevents x/0 type of cases

        lambda[i] = current_lambda;
        for (int j = 0; j < n; j++) {
            V[j][i] = v_temp[j]; 
        }

        for (int r = 0; r < n; r++) {
            for (int c = 0; c < n; c++) {
                C_temp[r][c] -= current_lambda * v_temp[r] * v_temp[c]; 
            }
        }
    }
    
    freeMatrix(C_temp, n);
    free(v_temp);
}



double** compute_U(double **A, double **V, double *lambda, int m, int n) {
    double **U = createMatrix(m, n); 
    if (U == NULL) return NULL;

    for (int i = 0; i < n; i++) {
        double sigma = sqrt(lambda[i]);
        if (sigma < 1.0e-9) { 
            for (int r = 0; r < m; r++) {
                U[r][i] = 0.0;
            }
            continue; 
        }

        double *v_i_vec = (double*)malloc(n * sizeof(double));
        if (v_i_vec == NULL) { 
            freeMatrix(U, m); 
            return NULL; 
        }
        for (int j = 0; j < n; j++) {
            v_i_vec[j] = V[j][i];
        }

        double *Av = (double*)malloc(m * sizeof(double));
        if (Av == NULL) { 
            free(v_i_vec); 
            freeMatrix(U, m); 
            return NULL; 
        }

        for (int r = 0; r < m; r++) {
            double sum = 0.0;
            for (int c = 0; c < n; c++) {
                sum += A[r][c] * v_i_vec[c];
            }
            Av[r] = sum;
        }

        for (int r = 0; r < m; r++) {
            U[r][i] = Av[r] / sigma; 
        }

        free(Av); 
        free(v_i_vec); 
    }//finds the corrosponding U columns

    return U;
}


double** final_image(double **U, double **V, double *lambda, int m, int n, int k) {
    double **A_k = createMatrix(m, n);
    if (A_k == NULL) return NULL;

    for (int i = 0; i < k; i++) {
        double sigma_i = sqrt(lambda[i]);
        if (sigma_i < 1.0e-9) {
            continue; 
        }

        for (int r = 0; r < m; r++) {
            for (int c = 0; c < n; c++) {
                A_k[r][c] += sigma_i * U[r][i] * V[c][i];
            }
        }
    }

    return A_k;
}//creates what we want, just in a matrix form
