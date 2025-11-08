#ifndef svd_h
#define svd_h

void power_iteration(double **C_temp, int n, double *v_out, double *lambda_out) ;

void svd_power_deflation(double **C, int n, int k, double **V, double *lambda);

double** compute_U(double **A, double **V, double *lambda, int m, int n);

double** final_image(double **U, double **V, double *lambda, int m, int n, int k);

#endif // svd_h
