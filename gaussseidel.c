#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"xiaohe.h"
void gaussseidel(float aa[],float r[],float b[],int n)
{ 
	/****Gauss-seidel iteration****/
	int *js,l,k,i,j,is;
	float d,t;
	float **R;

	js = (int *)calloc(n, sizeof(int));
	R  = space2d(n, n);
	l=1;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			R[i][j]=r[abs(i-j)];
	for (k=0; k<=n-2; k++)
	{ 
		d=0.0;
		for (i=k;i<=n-1;i++)
			for (j=k;j<=n-1;j++)
			{
				t=fabs(R[i][j]);
				if (t>d) { d=t; js[k]=j; is=i;}
			}
		if (d+1.0==1.0) l=0;
		else
		{ 
			if (js[k]!=k)
				for (i=0;i<=n-1;i++)
				{
					t=R[i][k]; 
					R[i][k]=R[i][js[k]]; 
					R[i][js[k]]=t;
				}
			if (is!=k)
			{
				for (j=k;j<=n-1;j++)
				{ 
					t=R[k][j]; 
					R[k][j]=R[is][j]; 
					R[is][j]=t;
				}
				t=b[k]; b[k]=b[is]; b[is]=t;
			}
		}
		if (l==0)
		{ 
			free(js);
			printf("\nÏµÊýŸØÕóÆæÒì£¡ÎÞœâ.");
			return;
		}
		d=R[k][k];
		for (j=k+1;j<=n-1;j++)
			R[k][j]=R[k][j]/d;
		b[k]=b[k]/d;
		for (i=k+1;i<=n-1;i++)
		{ 
			for (j=k+1;j<=n-1;j++)
				R[i][j]=R[i][j]-R[i][k]*R[k][j];
			b[i]=b[i]-R[i][k]*b[k];
		}
	}
	d=R[n-1][n-1];
	if (fabs(d)+1.0==1.0)
	{ 
		free(js);
		printf("\nÏµÊýŸØÕóÆæÒì£¡ÎÞœâ.");
		return;
	}
	b[n-1]=b[n-1]/d;
	for (i=n-2;i>=0;i--)
	{ 
		t=0.0;
		for (j=i+1;j<=n-1;j++)
			t=t+R[i][j]*b[j];
		b[i]=b[i]-t;
	}
	js[n-1]=n-1;
	for (k=n-1;k>=0;k--)
		if (js[k]!=k)
		{ 
			t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;
		}
	for(i=0;i<n;i++)
		aa[i]=b[i];
	free(js);
	free_space2d(R,n);
}
