/* Some small codes for daily life */
#include "xiaohe.h"
/*Data:10/10/2016*/
//***The application of 2D array in float
float **space2d(int nr,int nc)
{
	float **a;
	int i;
	a=(float **)calloc(nr,sizeof(float *));
	for (i=0; i<nr; i++)
		a[i]=(float *)calloc(nc,sizeof(float));
	return a;
}

//***The application of 2D array in int
int **space2dint(int nr,int nc)
{
	int **a;
	int i;
	a=(int **)calloc(nr,sizeof(int *));
	for (i=0; i<nr; i++)
		a[i]=(int *)calloc(nc,sizeof(int));
	return a;
}
//---To release a 2D array in float
void free_space2d(float **a,int nr)
{
	int i;
	for (i=0; i<nr; i++)
		free(a[i]);
	free(a);
}
//---To release a 2D array in int
void free_space2dint(int **a,int nr)
{
	int i;
	for (i=0; i<nr; i++)
		free(a[i]);
	free(a);
}

//////////////////////////////////////////////////////////////////////

//---DFFT---//
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
			in_refl[i]=x_in[i] + 0.0 * I;
		else 
			in_refl[i]=0.0 + 0.0 * I;
	}

	p = fftwf_plan_dft_1d(Nfft, in_refl, out_refl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p);

	for(i=0;i<Nfft;i++) X[i] = crealf(out_refl[i]) + cimagf(out_refl[i]) * I ; 
	
	free(out_refl);
	free(in_refl);
	fftwf_destroy_plan(p);
}

//---IDFFT---//
void  IDFFT(float complex *X_in, float *x, int Nx, int Nfft)
{
	int i;
	fftwf_complex *out_refl,*in_refl;
	fftwf_plan p;

	out_refl = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * Nfft);
	in_refl  = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * Nfft);
	
	for(i=0;i<Nfft;i++)   in_refl[i] = crealf(X_in[i]) + cimagf(X_in[i])*I;

	p = fftwf_plan_dft_1d(Nfft, in_refl, out_refl, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(p);

	for(i=0;i<Nx;i++) x[i] = crealf(out_refl[i])/Nfft;
	
	free(out_refl);
	free(in_refl);
	fftwf_destroy_plan(p);
}
//////////////////////////////////////////////////////////////////////

void xcorr(float *r, float *x, float *y, int N)
{
	float sxy;
	int delay, i, j;
	for(delay=-N+1;delay<N;delay++)
	{
		sxy=0;
		for(i=0; i<N; i++)
		{
			j=i+delay;
			if((j < 0) || (j >= N))
				continue;
			else
				sxy+=(x[j]*y[i]);
		}
		r[delay + N - 1] = sxy;
	}
}

void convolution(float *X, float *refl, int Nr, float *wavelet, int Nw, int flag)
{
	// ***covolution &cross-correlationin frequency domain*** //
	// ***If flag==0:convolution; else : cross-correlation*** //
	int i, Nlen, Nmax, Nmin;
	Nlen   = Nr + Nw - 1;
	Nmax = ((Nr) > (Nw) ? (Nr) : (Nw));
	Nmin = ((Nr) < (Nw) ? (Nr) : (Nw));

	fftwf_complex *con_out_w;
	fftwf_complex *in_w, *out_w;
	fftwf_complex *in_r, *out_r;
	fftwf_complex *in_x, *out_x;

	fftwf_plan p1,p2,p3;

	con_out_w = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);
	
	in_r  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);
	out_r = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);

	in_w  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);
	out_w = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);

	in_x  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);
	out_x = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * Nlen);
	
	for(i=0;i<Nlen;i++)
	{
		if(i<Nr) in_r[i] = refl[i] + 0.0*I;
		else 	 in_r[i] = 0.0 + 0.0*I;

		if(i<Nw) in_w[i] = wavelet[i] + 0*I;
		else	 in_w[i] = 0 + 0*I;
	}

	p1 = fftwf_plan_dft_1d(Nlen,in_r,out_r,FFTW_FORWARD,FFTW_ESTIMATE);
	fftwf_execute(p1);

	p2 = fftwf_plan_dft_1d(Nlen,in_w,out_w,FFTW_FORWARD,FFTW_ESTIMATE);
	fftwf_execute(p2);

	if (flag==0)
	{
		for(i=0;i<Nlen;i++)
			in_x[i] = out_r[i]*out_w[i];
		p3 = fftwf_plan_dft_1d(Nlen,in_x,out_x,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(p3);

		for(i=0; i<Nmax; i++)
			X[i] = crealf(out_x[i+Nmin-1])/Nlen;	
	}	
	else 
	{
		for(i=0;i<Nlen;i++){
			// ***con_out_w is the conjugation of out_w*** //
			con_out_w[i] = crealf(out_w[i]) - cimagf(out_w[i])*I;
			in_x[i]      = con_out_w[i] * out_r[i];}

		p3 = fftwf_plan_dft_1d(Nlen,in_x,out_x,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(p3);

		for(i=0;i<Nmax;i++)
			X[i] = crealf(out_x[i])/Nlen;
	}

	free(con_out_w);

	free(in_r);
	free(out_r);

	free(in_w);
	free(out_w);

	free(in_x);
	free(out_x);

	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);
	fftwf_destroy_plan(p3);
}
	
/*------Calculation of higher order differential coef    ficients:c(m)-----*/
float *FDmodulus(int Mm)
{
	int i,j,m;
	int imax;
	float s,max,t1;
	float *Cc;
	float **A;
	A=space2d(Mm+1,Mm+2);
	Cc=(float*)calloc(Mm+1,sizeof(float));

	for(i=1;i<=Mm;i++)
	{ 
		for(j=1;j<=Mm;j++)
			A[i][j]=pow(j,2*i);
		A[i][Mm+1]=0;
		A[1][Mm+1]=1; 
	}
for(j=1;j<=Mm-1;j++)         
	{
		max=0.0;
	   for(i=j;i<=Mm;i++)
		   if(fabs(A[i][j])>max)
		   {	
			   max=A[i][j];
		       imax=i;
		}
	    if(max==0)
		{
			printf("stop\n");
		    break;
	    }
	    if(imax!=j)
	    for(i=j;i<=Mm+1;i++)
       	{
			t1=A[j][i];
		    A[j][i]=A[imax][i];
		    A[imax][i]=t1;
        }
	    for(i=j+1;i<=Mm;i++)
	    {
			A[i][j]=A[i][j]/A[j][j];
	        for(m=j+1;m<=Mm+1;m++)
	            A[i][m]=A[i][m]-A[i][j]*A[j][m];
	    }
	}
	if(A[Mm][Mm]==0)
	    printf("stop\n");
	else
	for(i=Mm;i>=1;i--)
	{
		s=0.0;
		for(j=i+1;j<=Mm;j++)
			s=s+A[i][j]*Cc[j];
		Cc[i]=(A[i][Mm+1]-s)/A[i][i];
	}
	free_space2d(A,Mm);
	return Cc;
}

float source(int k,float f0,float t0,float dt,float AP)	//Richer wavelet
{
	float y=AP*(1-2*pow((PI*f0*(k*dt-t0)),2))*exp(-pow((PI*f0*(k*dt-t0)),2));
	return y;
}

float d(float x,int deta,float c)   	//For PML - d(x)
{
	float y=log(1000.0)*1.5*c*pow(x/deta,2)/deta;
	return y;
}

float dd(float x,int deta,float c)	//For PML - d'(x)
{
	float y=log(1000.0)*3*c*x/pow(deta,3);
	return y;
}

float **model(int Nxx,int Nzz ,int nn, int M)		//***Velocity model
{
	FILE *fp_model, *fp_out;
	int i,j;
	float **vel;
	int Nz, Nx;
	int upbrim=nn+M;
	Nx=Nxx+nn+nn+M+M;
	Nz=Nzz+nn+M+upbrim;
	fp_model=fopen("velocity.dat","rb");
	fp_out=fopen("velocity1.dat","wb");

	vel=space2d(Nx,Nz);

	for(i=0;i<Nxx;i++)
		for(j=0;j<Nzz;j++)                     
			fread(&vel[i+nn+M][j+upbrim],4L,1,fp_model);

	for(i=0;i<Nx;i++)
		for(j=0;j<Nz;j++)
		{
			if(j<upbrim)	//Up
				vel[i][j]=vel[nn+M][upbrim];
			if(i<nn+M && j>=upbrim && j<Nz-nn-M)	//Left
				vel[i][j]=vel[nn+M][j];
			if(i>=Nx-nn-M && j>=upbrim && j<Nz-nn-M)	//Right
				vel[i][j]=vel[Nx-nn-M-1][j];
			if(j>=Nz-nn-M)	//Bottom
				vel[i][j]=vel[i][Nz-nn-M-1];
			fwrite(&vel[i][j],4L,1,fp_out);
		}
	fclose(fp_model);
	fclose(fp_out);
	return vel;
}

float **acoustic_modeling_2d(float *para, int num, int **receiver, int Nreceiver)
{
	FILE *fpssp,*fpvsp,*fpsnap;
	float AP=1.0, dx, dz, f0, dt, t0=0.035, t1, max, s;
	int i, j, k, m, x, z, deta;	//deta is PML's width
	int Nx, Nz, Nx0, Nz0, T, M, nn, t_snap, flag, x0, z0;
	float px, pz;
	int Flag_M, upbrim, leftbrim;
	int i_receiver, receiver_x, receiver_z;
	Nx0 = (int) para[0];
	Nz0 = (int) para[1];
	T = (int) para[2];
	M = (int) para[3];
	nn = (int) para[4];
	upbrim=nn+M;
	leftbrim=nn+M;
	Nx = Nx0+nn+M+nn+M;
	Nz = Nz0+nn+M+upbrim;
	t_snap = (int) para[5];
	flag = (int) para[6];
	f0 = para[7];
	x0 = (int) para[8]+leftbrim;
	z0 = (int) para[9]+upbrim;
	dx = para[10];
	dz = para[11];
	dt = para[12];
	
	//Application for 2d arraies
	float *c, *cc;
	float **v, **v_model;	//Velocity
	float **V;	//Vertical Seismic Profile of SSP.      
	float **p;
	float **temp1,**temp2;
	float **q;
	
	float **u1,**u2,**u3;
	float **temp1u1,**temp1u2,**temp1u3;
	float **temp2u1,**temp2u2,**temp2u3;
	float **w1,**temp1w1,**temp2w1;

	float **u21,**u22,**u23;
	float **temp1u21,**temp1u22,**temp1u23;
	float **temp2u21,**temp2u22,**temp2u23;
	float **w2,**temp1w2,**temp2w2;

	float **u31,**u32,**u33;
	float **temp1u31,**temp1u32,**temp1u33;
	float **temp2u31,**temp2u32,**temp2u33;
	float **w3,**temp1w3,**temp2w3;

	float **u41,**u42,**u43;
	float **temp1u41,**temp1u42,**temp1u43;
	float **temp2u41,**temp2u42,**temp2u43;
	float **w4,**temp1w4,**temp2w4;
	if(flag==3)
	{	V=space2d(Nreceiver,T); q=space2d(Nreceiver,T);}
	if(flag==2)
	{	V=space2d(Nz,T); q=space2d(para[1],T);}
	if(flag==1)
	{	V=space2d(Nx,T); q=space2d(para[0],T);}
	if(flag==0)
	{	V=space2d(Nx,Nz); q=space2d(para[0],para[1]);}
	p=space2d(Nx,Nz);
	temp1=space2d(Nx,Nz);//For memorize the data from k & k-1 
	temp2=space2d(Nx,Nz);
	v=space2d(Nx,Nz);
	
	//------u1+u2+u3=p;w1,w2 is the medium of data u2 & u22----------*
	u1=space2d(Nx,Nz);
	u2=space2d(Nx,Nz);
	u3=space2d(Nx,Nz);
	temp1u1=space2d(Nx,Nz);
	temp1u2=space2d(Nx,Nz);
	temp1u3=space2d(Nx,Nz);
	temp2u1=space2d(Nx,Nz);
	temp2u2=space2d(Nx,Nz);
	temp2u3=space2d(Nx,Nz);
	w1=space2d(Nx,Nz);
	temp1w1=space2d(Nx,Nz);
	temp2w1=space2d(Nx,Nz);

	//------Left side-------*
	u21=space2d(Nx,Nz);
	u22=space2d(Nx,Nz);
	u23=space2d(Nx,Nz);
	temp1u21=space2d(Nx,Nz);
	temp1u22=space2d(Nx,Nz);
	temp1u23=space2d(Nx,Nz);
	temp2u21=space2d(Nx,Nz);
	temp2u22=space2d(Nx,Nz);
	temp2u23=space2d(Nx,Nz);
	w2=space2d(Nx,Nz);
	temp1w2=space2d(Nx,Nz);
	temp2w2=space2d(Nx,Nz);

	//-----Right side--------*
	u31=space2d(Nx,Nz);
	u32=space2d(Nx,Nz);
	u33=space2d(Nx,Nz);
	temp1u31=space2d(Nx,Nz);
	temp1u32=space2d(Nx,Nz);
	temp1u33=space2d(Nx,Nz);
	temp2u31=space2d(Nx,Nz);
	temp2u32=space2d(Nx,Nz);
	temp2u33=space2d(Nx,Nz);
	w3=space2d(Nx,Nz);
	temp1w3=space2d(Nx,Nz);
	temp2w3=space2d(Nx,Nz);

	//-----bottom side-----*
	u41=space2d(Nx,Nz);
	u42=space2d(Nx,Nz);
	u43=space2d(Nx,Nz);
	temp1u41=space2d(Nx,Nz);
	temp1u42=space2d(Nx,Nz);
	temp1u43=space2d(Nx,Nz);
	temp2u41=space2d(Nx,Nz);
	temp2u42=space2d(Nx,Nz);
	temp2u43=space2d(Nx,Nz);
	w4=space2d(Nx,Nz);
	temp1w4=space2d(Nx,Nz);
	temp2w4=space2d(Nx,Nz);	

	deta=nn*dx;
	c=FDmodulus(M);//Calculation of differential coefficients	
	v=model(Nx0,Nz0,nn,M);	////////////////////////////////////////////////////

	for(k=0;k<T;k++)
	{
	   	for(i=0;i<Nx;i++)
			for(j=0;j<Nz;j++)
			{
				temp1[i][j]=temp2[i][j];             //For exchange
				temp2[i][j]=p[i][j];

			   	//----------Up side-------*
				temp1u1[i][j]=temp2u1[i][j];
			   	temp2u1[i][j]=u1[i][j];
				temp1u2[i][j]=temp2u2[i][j];
			   	temp2u2[i][j]=u2[i][j];
			   	temp1u3[i][j]=temp2u3[i][j];
			   	temp2u3[i][j]=u3[i][j];
			   	temp1w1[i][j]=temp2w1[i][j];
			   	temp2w1[i][j]=w1[i][j];

			   	//----------Left side---------*
				temp1u21[i][j]=temp2u21[i][j];
			   	temp2u21[i][j]=u21[i][j];
			   	temp1u22[i][j]=temp2u22[i][j];
			   	temp2u22[i][j]=u22[i][j];
			   	temp1u23[i][j]=temp2u23[i][j];
			   	temp2u23[i][j]=u23[i][j];
			   	temp1w2[i][j]=temp2w2[i][j];
			   	temp2w2[i][j]=w2[i][j];

			  	//----------Right side---------*
				temp1u31[i][j]=temp2u31[i][j];
			   	temp2u31[i][j]=u31[i][j];
			   	temp1u32[i][j]=temp2u32[i][j];
			   	temp2u32[i][j]=u32[i][j];
			   	temp1u33[i][j]=temp2u33[i][j];
			   	temp2u33[i][j]=u33[i][j];
			   	temp1w3[i][j]=temp2w3[i][j];
			   	temp2w3[i][j]=w3[i][j];

               	//----------Bottom side---------*
			   	temp1u41[i][j]=temp2u41[i][j];
			   	temp2u41[i][j]=u41[i][j];
			   	temp1u42[i][j]=temp2u42[i][j];
			   	temp2u42[i][j]=u42[i][j];
			   	temp1u43[i][j]=temp2u43[i][j];
			   	temp2u43[i][j]=u43[i][j];
			   	temp1w4[i][j]=temp2w4[i][j];
			   	temp2w4[i][j]=w4[i][j];

			}
		for(j=M;j<Nz-M;j++)	
		{
			for(i=M;i<Nx-M;i++)
			{
				px=0;
				pz=0;
				for(m=1;m<=M;m++)
				{
					px=px+c[m]*(temp2[i+m][j]-2*temp2[i][j]+temp2[i-m][j]); 
					pz=pz+c[m]*(temp2[i][j+m]-2*temp2[i][j]+temp2[i][j-m]);                               
				}
				//--------------Middle space---------*
				if(j>=nn+M && i>=nn+M && i<=Nx-nn-M && j<=Nz-nn-M)
				{
					if(i==x0&&j==z0)  
					{	//------Wavelet: Richer-------*
						p[i][j]=pow((v[i][j]*dt),2)*(px/(dx*dx)+pz/(dz*dz))
								+source(k,f0,t0,dt,AP)+2*temp2[i][j]-temp1[i][j];   
//						pow((v[i][j]*dt),2)*source(k,f0,t0,dt,AP)+2*temp2[i][j]-temp1[i][j];   
					}
					else
						p[i][j]=pow((v[i][j]*dt),2)*(px/(dx*dx)+pz/(dz*dz))+2*temp2[i][j]-temp1[i][j];
				}

				//----------------Up---------------*
				else if(j<=nn+M && i>=j && i+j<=Nx)
				{
					z=(nn+M-j)*dz; 
					
					u1[i][j]=pow((v[i][j]*dt/dz),2)*pz-(pow((d(z,deta,v[i][j])*dt+1),2)-3)*temp2u1[i][j]
							+(2*d(z,deta,v[i][j])*dt-1)*temp1u1[i][j];
					w1[i][j]=(3-pow((1+dt*d(z,deta,v[i][j])),2))*temp2w1[i][j]+(2*d(z,deta,v[i][j])*dt-1)
						*temp1w1[i][j]+pow(v[i][j]*dt,2)*dd(z,deta,v[i][j])*(temp2[i][j]-temp2[i][j-1])/dz;
					u2[i][j]=dt*w1[i][j]+temp2u2[i][j]*(1-dt*d(z,deta,v[i][j]));
   					u3[i][j]=pow((v[i][j]*dt/dx),2)*px+2*temp2u3[i][j]-temp1u3[i][j];  
					p[i][j]=u1[i][j]+u2[i][j]+u3[i][j];
				}
				
				//-------------Left-----------------*
				else if(i<=nn+M && i+j<=Nz && j>=i)
				{
					x=(nn+M-i)*dx;

					u21[i][j]=pow((v[i][j]*dt/dx),2)*px-(pow((d(x,deta,v[i][j])*dt+1),2)-3)*temp2u21[i][j]
							+(2*d(x,deta,v[i][j])*dt-1)*temp1u21[i][j];
					w2[i][j]=(3-pow((1+dt*d(x,deta,v[i][j])),2))*temp2w2[i][j]+(2*d(x,deta,v[i][j])*dt-1)
						*temp1w2[i][j]+pow(v[i][j]*dt,2)*dd(x,deta,v[i][j])*(temp2[i][j]-temp2[i-1][j])/dx;
					u22[i][j]=dt*w2[i][j]+temp2u22[i][j]*(1-dt*d(x,deta,v[i][j]));
   					u23[i][j]=pow((v[i][j]*dt/dz),2)*pz+2*temp2u23[i][j]-temp1u23[i][j];  
					p[i][j]=u21[i][j]+u22[i][j]+u23[i][j];
				}
				//---------------Right------i+j>=Nx---------*
				else if(i>=Nx-nn-M && i>=j+Nx-Nz && i+j>=Nx)
				{
					x=(i-Nx+nn+M)*dx;
					
					u31[i][j]=pow((v[i][j]*dt/dx),2)*px-(pow((d(x,deta,v[i][j])*dt+1),2)-3)*temp2u31[i][j]
							+(2*d(x,deta,v[i][j])*dt-1)*temp1u31[i][j];
					w3[i][j]=(3-pow((1+dt*d(x,deta,v[i][j])),2))*temp2w3[i][j]+(2*d(x,deta,v[i][j])*dt-1)
						*temp1w3[i][j]-pow(v[i][j]*dt,2)*dd(x,deta,v[i][j])*(temp2[i][j]-temp2[i-1][j])/dx;
					u32[i][j]=dt*w3[i][j]+temp2u32[i][j]*(1-dt*d(x,deta,v[i][j]));
   					u33[i][j]=pow((v[i][j]*dt/dz),2)*pz+2*temp2u33[i][j]-temp1u33[i][j];  
					p[i][j]=u31[i][j]+u32[i][j]+u33[i][j];
				}
				//---------------Bottom--------j>N-nn-M&&j>=i&&i+j>=N-------*
				else if(j>=Nz-nn-M && j>=i+Nz-Nx && i+j>=Nz)
				{
					z=(j+nn+M-Nz)*dz;

					u1[i][j]=pow((v[i][j]*dt/dz),2)*pz-(pow((d(z,deta,v[i][j])*dt+1),2)-3)*temp2u1[i][j]
							+(2*d(z,deta,v[i][j])*dt-1)*temp1u1[i][j];
					w1[i][j]=(3-pow((1+dt*d(z,deta,v[i][j])),2))*temp2w1[i][j]+(2*d(z,deta,v[i][j])*dt-1)
						*temp1w1[i][j]-pow(v[i][j]*dt,2)*dd(z,deta,v[i][j])*(temp2[i][j]-temp2[i][j-1])/dz;
					u2[i][j]=dt*w1[i][j]+temp2u2[i][j]*(1-dt*d(z,deta,v[i][j]));
   					u3[i][j]=pow((v[i][j]*dt/dx),2)*px+2*temp2u3[i][j]-temp1u3[i][j];  
					p[i][j]=u1[i][j]+u2[i][j]+u3[i][j];
				}
				if(flag==0)
					if(k==t_snap)
						V[i][j]=p[i][j];		//Snap shot
				if(flag==1)
					if(j==upbrim+220)
						V[i][k]=p[i][j];  		//SSP
				if(flag==2)	
					if(i==x0)
						V[j][k]=p[i][1+upbrim];	//VSP
			}
		}
		if(flag==3)
			for(i_receiver=0; i_receiver<Nreceiver; i_receiver++ )
			{
				receiver_x = receiver[i_receiver][0] + leftbrim;
				receiver_z = receiver[i_receiver][1] + upbrim;
				V[i_receiver][k]=p[receiver_x][receiver_z];
			}
	}
	if(flag==0)
	{
		for(i=0;i<para[0];i++)
			for(j=0;j<para[1];j++)
				q[i][j]=V[i+leftbrim][j+upbrim];		//Snap shot
		free_space2d(V,para[0]);
	}
	if(flag==1)
	{
		for(i=0;i<para[0];i++)
			for(j=0;j<para[2];j++)
				q[i][j]=V[i+leftbrim][j];		//SSP
		free_space2d(V,para[0]);
	}
	if(flag==2)
	{
		for(i=0;i<para[1];i++)
			for(j=0;j<para[2];j++)
				q[i][j]=V[i+upbrim][j];		//VSP
		free_space2d(V,para[1]);
	}
	if(flag==3)
	{
		for(i=0;i<Nreceiver;i++)
			for(j=0;j<para[2];j++)
				q[i][j]=V[i][j];		//SSP
		free_space2d(V,Nreceiver);
	}
	
	free(c);
	free_space2d(v,Nx);
	free_space2d(p,Nx);
	free_space2d(temp1,Nx);
	free_space2d(temp2,Nx);
	
	free_space2d(u1,Nx);
	free_space2d(u2,Nx);
	free_space2d(u3,Nx);
	free_space2d(temp1u1,Nx);
	free_space2d(temp1u2,Nx);
	free_space2d(temp1u3,Nx);
	free_space2d(temp2u1,Nx);
	free_space2d(temp2u2,Nx);
	free_space2d(temp2u3,Nx);
	free_space2d(w1,Nx);
	free_space2d(temp1w1,Nx);
	free_space2d(temp2w1,Nx);
	
	//--------Left side-------*
	free_space2d(u21,Nx);
	free_space2d(u22,Nx);
	free_space2d(u23,Nx);
	free_space2d(temp1u21,Nx);
	free_space2d(temp1u22,Nx);
	free_space2d(temp1u23,Nx);
	free_space2d(temp2u21,Nx);
	free_space2d(temp2u22,Nx);
	free_space2d(temp2u23,Nx);
	free_space2d(w2,Nx);
	free_space2d(temp1w2,Nx);
	free_space2d(temp2w2,Nx);

	//-------Right side-------*
	free_space2d(u31,Nx);
	free_space2d(u32,Nx);
	free_space2d(u33,Nx);
	free_space2d(temp1u31,Nx);
	free_space2d(temp1u32,Nx);
	free_space2d(temp1u33,Nx);
	free_space2d(temp2u31,Nx);
	free_space2d(temp2u32,Nx);
	free_space2d(temp2u33,Nx);
	free_space2d(w3,Nx);
	free_space2d(temp1w3,Nx);
	free_space2d(temp2w3,Nx);
	
	//------Bottom side------*
	free_space2d(u41,Nx);
	free_space2d(u42,Nx);
	free_space2d(u43,Nx);
	free_space2d(temp1u41,Nx);
	free_space2d(temp1u42,Nx);
	free_space2d(temp1u43,Nx);
	free_space2d(temp2u41,Nx);
	free_space2d(temp2u42,Nx);
	free_space2d(temp2u43,Nx);
	free_space2d(w4,Nx);
	free_space2d(temp1w4,Nx);
	free_space2d(temp2w4,Nx);
	return q;

}

