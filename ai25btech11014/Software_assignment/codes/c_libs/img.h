#ifndef img_h
#define img_h

double** read_img(const char *file, int *m, int *n);

void write_img(const char *file, double **A, int m, int n, const char *format);


#endif // img_h