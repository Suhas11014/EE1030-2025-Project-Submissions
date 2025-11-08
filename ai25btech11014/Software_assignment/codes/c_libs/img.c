#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "matrixfunc.h"
#include "img.h"

double** read_img(const char *file, int *m, int *n) {
    int width, height, channels;
    double** A = NULL;

    unsigned char *data = stbi_load(file, &width, &height, &channels, 1); 
    
    if (data == NULL) {
        fprintf(stderr, "Error: Could not load image file %s. Check path and format.\n", file);
        return NULL;
    }//coudnt load image

    *n = width;
    *m = height;

    A = createMatrix(*m, *n);
    if (A == NULL) {
        fprintf(stderr, "Error: Matrix allocation failed during image reading.\n");
        stbi_image_free(data);
        return NULL;
    }//couldnt create matrix

    for (int i = 0; i < *m; i++) {
        for (int j = 0; j < *n; j++) {
            int index = (i * width + j);
            A[i][j] = (double)data[index];
        }
    }//give the matrix the img values

    stbi_image_free(data);
    return A;
}

void write_img(const char *file, double **A, int m, int n, const char *format) {
    const int CHANNELS = 3; 
    int data_size = m * n * CHANNELS;
    unsigned char *data = (unsigned char *)malloc(data_size);
    if (data == NULL) {
        fprintf(stderr, "Error: Failed to allocate output image buffer.\n");
        return;
    }//data buffer allocation failed

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int index = (i * n + j) * CHANNELS;
            double val = A[i][j];
            val = fmax(0.0, fmin(255.0, val));
            unsigned char pixel_val = (unsigned char)round(val);
            data[index + 0] = pixel_val; 
            data[index + 1] = pixel_val; 
            data[index + 2] = pixel_val; 
        }
    }//create the data buffer from the matrix

    int success = 0;
    
    if (strcmp(format, "png") == 0) {
        success = stbi_write_png(file, n, m, CHANNELS, data, n * CHANNELS);
    } else if (strcmp(format, "jpg") == 0 || strcmp(format, "jpeg") == 0) {
        success = stbi_write_jpg(file, n, m, CHANNELS, data, 90);
    } else {
        fprintf(stderr, "Error: Unsupported output format: %s. Use 'png' or 'jpg'.\n", format);
    }//check format and write accordingly

    if (!success) {
        fprintf(stderr, "Error: Failed to write image to %s.\n", file);
    } else {
    }

    free(data);
}



