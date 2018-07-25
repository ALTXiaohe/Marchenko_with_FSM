#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
void  DFFT(float complex *X, float *x_in, int Nx, int Nfft)
{
	int i;
	fftwf_complex *out_refl,*in_refl;
	fftwf_plan p;
	out_refl = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * Nfft);
	in_refl  = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * Nfft);
	for(i=0;i<Nfft;i++) 
	{
		if (i < Nx)
			in_refl[i]=x_in[i] + 0*I;
		else 
			in_refl[i]=0.0 + 0.0 * I;
	}
	p = fftwf_plan_dft_1d(Nt, in_refl, out_refl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p);

	for(i=0;i<Nfft;i++) X[i] = crealf(out_refl[i]) + cimagf(out_refl[i]) * I ;

	free(out_refl);
	free(in_refl);
	fftwf_destroy_plan(p);
}
