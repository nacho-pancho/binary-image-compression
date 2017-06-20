#ifndef PNM_H
#define PNM_H
#include <cstdio>
//#include "img.h"

typedef unsigned int pixel_t;

int write_pgm(const pixel_t* pixels, int tipo, int ancho, int alto, int maxval, const char* ruta_archivo);
int read_pnm_header(FILE* f, int& tipo, int& ancho, int& alto, int &maxval);
int read_pgm_data(FILE* f, int tipo, int ancho, int alto, int maxval, pixel_t* buf);

int write_ppm(const pixel_t* pixels, int tipo, int ancho, int alto, int maxval, const char* ruta_archivo);
int write_ppm_header(int tipo, int ancho, int alto, int maxval, FILE* ruta_archivo);
int write_p2_data(const pixel_t* pixels, const int npixels, const int maxval, FILE* fw);
int write_p5_data(const pixel_t* pixels, const int npixels, const int maxval, FILE* fw);
int read_ppm_data(FILE* f, int tipo, int ancho, int alto, int maxval, pixel_t* buf);


#endif
