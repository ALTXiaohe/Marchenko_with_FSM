/* xiaohe.h - include file for general purpose */

#ifndef XIAOHE_H
#define XIAOHE_H
/* INCLUDES */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>

#ifndef size_t
#define size_t int
#endif

#ifndef PI
#define PI (3.141592653589793)
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

float **space2d(int nr,int nc);
void  free_space2d(float **a,int nr);
int   **space2dint(int nr,int nc);
void  free_space2dint(int **a,int nr);
float *FDmodulus(int Mm);
float source(int k,float f0,float t0,float dt,float AP);
float d(float x,int deta,float c);
float dd(float x,int deta,float c);
float **model(int Nxx,int Nzz, int nn, int M);
float **acoustic_modeling_2d(float *para, int num, int **receiver, int Nreceiver);
void  fourier(float _Complex XXf[],float XX[],int Nt);
void  ifourier(float _Complex XXt[], float _Complex XXf[], int Nt);
void  DFFT(float complex *X, float *x_in, int Nx, int Nfft);
void  IDFFT(float complex *X_in, float *x, int Nx, int Nfft);
void convolution(float *X, float *refl, int Nr, float *wavelet, int Nw, int flag);//--flag==0 ->convolution; flag==else ->correlation
void xcorr(float *r, float *x, float *y, int N);
#endif
