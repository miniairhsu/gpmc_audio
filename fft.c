/* 
 * Free FFT and convolution (C)
 * 
 * Copyright (c) 2014 Project Nayuki
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 * 
 * (MIT License)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fft.h"


// Private function prototypes
static size_t reverse_bits(size_t x, unsigned int n);
static void *memdup(const void *src, size_t n);

#define SIZE_MAX ((size_t)-1)

/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
short FFT(short int dir,long m,double *x,double *y)
{
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return(1);
}

int transform(double real[], double imag[], size_t n) {
	if (n == 0)
		return 1;
	else if ((n & (n - 1)) == 0)  // Is power of 2
		return transform_radix2(real, imag, n);
	else  // More complicated algorithm for arbitrary sizes
		return transform_bluestein(real, imag, n);
}


int inverse_transform(double real[], double imag[], size_t n) {
	return transform(imag, real, n);
}


int transform_radix2(double real_in[], double imag_in[], size_t n) {
	// Variables
	int status = 0;
	unsigned int levels;
	double *cos_table, *sin_table;
	size_t size;
	size_t i;
	double *real;
	double *imag;
	real = malloc(n * sizeof(double));
	imag = malloc(n * sizeof(double));
	memmove(real, real_in, n*sizeof(double));
	memmove(imag, imag_in, n*sizeof(double));
	memset(real_in, 0, n*sizeof(double));
	memset(imag_in, 0, n*sizeof(double));
	// Compute levels = floor(log2(n))
	//{
		size_t temp = n;
		levels = 0;
		while (temp > 1) {
			levels++;
			temp >>= 1;
		}
		if (1u << levels != n)
			return 0;  // n is not a power of 2
	//}
	
	// Trignometric tables
	if (SIZE_MAX / sizeof(double) < n / 2)
		return 0;
	size = (n / 2) * sizeof(double);
	cos_table = malloc(size);
	sin_table = malloc(size);
	if (cos_table == NULL || sin_table == NULL)
		goto cleanup;
	for (i = 0; i < n / 2; i++) {
		cos_table[i] = cos(2 * M_PI * i / n);
		sin_table[i] = sin(2 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (i = 0; i < n; i++) {
		size_t j = reverse_bits(i, levels);
		if (j > i) {
			double temp = real[i];
			real[i] = real[j];
			real[j] = temp;
			temp = imag[i];
			imag[i] = imag[j];
			imag[j] = temp;
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (i = 0; i < n; i += size) {
			size_t j;
			size_t k;
			for (j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				double tpre =  real[j+halfsize] * cos_table[k] + imag[j+halfsize] * sin_table[k];
				double tpim = -real[j+halfsize] * sin_table[k] + imag[j+halfsize] * cos_table[k];
				real[j + halfsize] = real[j] - tpre;
				imag[j + halfsize] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
	status = 1;
	memmove(real_in, real, n*sizeof(double));
	memmove(imag_in, imag, n*sizeof(double));
	
cleanup:
	free(cos_table);
	free(sin_table);
    free(real);
	free(imag);
	return status;
}


int transform_bluestein(double real_in[], double imag_in[], size_t n) {
	// Variables
	int status = 0;
	double *cos_table, *sin_table;
	double *areal, *aimag;
	double *breal, *bimag;
	double *creal, *cimag;
	size_t m;
	size_t size_n, size_m;
	size_t i;
	double *real;
	double *imag;
	real = malloc(n * sizeof(double));
	imag = malloc(n * sizeof(double));
	memmove(real, real_in, n*sizeof(double));
	memmove(imag, imag_in, n*sizeof(double));
	memset(real_in, 0, n*sizeof(double));
	memset(imag_in, 0, n*sizeof(double));
	// Find a power-of-2 convolution length m such that m >= n * 2 + 1
	{
		size_t target;
		if (n > (SIZE_MAX - 1) / 2)
			return 0;
		target = n * 2 + 1;
		for (m = 1; m < target; m *= 2) {
			if (SIZE_MAX / 2 < m)
				return 0;
		}
	}
	
	// Allocate memory
	if (SIZE_MAX / sizeof(double) < n || SIZE_MAX / sizeof(double) < m)
		return 0;
	size_n = n * sizeof(double);
	size_m = m * sizeof(double);
	cos_table = malloc(size_n);
	sin_table = malloc(size_n);
	areal = calloc(m, sizeof(double));
	aimag = calloc(m, sizeof(double));
	breal = calloc(m, sizeof(double));
	bimag = calloc(m, sizeof(double));
	creal = malloc(size_m);
	cimag = malloc(size_m);
	if (cos_table == NULL || sin_table == NULL
			|| areal == NULL || aimag == NULL
			|| breal == NULL || bimag == NULL
			|| creal == NULL || cimag == NULL)
		goto cleanup;
	
	// Trignometric tables
	for (i = 0; i < n; i++) {
		double temp = M_PI * (size_t)((unsigned long long)i * i % ((unsigned long long)n * 2)) / n;
		// Less accurate version if long long is unavailable: double temp = M_PI * i * i / n;
		cos_table[i] = cos(temp);
		sin_table[i] = sin(temp);
	}
	
	// Temporary vectors and preprocessing
	for (i = 0; i < n; i++) {
		areal[i] =  real[i] * cos_table[i] + imag[i] * sin_table[i];
		aimag[i] = -real[i] * sin_table[i] + imag[i] * cos_table[i];
	}
	breal[0] = cos_table[0];
	bimag[0] = sin_table[0];
	for (i = 1; i < n; i++) {
		breal[i] = breal[m - i] = cos_table[i];
		bimag[i] = bimag[m - i] = sin_table[i];
	}
	
	// Convolution
	if (!convolve_complex(areal, aimag, breal, bimag, creal, cimag, m))
		goto cleanup;
	
	// Postprocessing
	for (i = 0; i < n; i++) {
		real[i] =  creal[i] * cos_table[i] + cimag[i] * sin_table[i];
		imag[i] = -creal[i] * sin_table[i] + cimag[i] * cos_table[i];
	}
	status = 1;
	memmove(real_in, real, n*sizeof(double));
	memmove(imag_in, imag, n*sizeof(double));
	// Deallocation
cleanup:
	free(cimag);
	free(creal);
	free(bimag);
	free(breal);
	free(aimag);
	free(areal);
	free(sin_table);
	free(cos_table);
	free(real);
	free(imag);
	return status;
}


int convolve_real(const double x[], const double y[], double out[], size_t n) {
	double *ximag, *yimag, *zimag;
	int status = 0;
	ximag = calloc(n, sizeof(double));
	yimag = calloc(n, sizeof(double));
	zimag = calloc(n, sizeof(double));
	if (ximag == NULL || yimag == NULL || zimag == NULL)
		goto cleanup;
	
	status = convolve_complex(x, ximag, y, yimag, out, zimag, n);
cleanup:
	free(zimag);
	free(yimag);
	free(ximag);
	return status;
}


int convolve_complex(const double xreal[], const double ximag[], const double yreal[], const double yimag[], double outreal[], double outimag[], size_t n) {
	int status = 0;
	size_t size;
	size_t i;
	double *xr, *xi, *yr, *yi;
	if (SIZE_MAX / sizeof(double) < n)
		return 0;
	size = n * sizeof(double);
	xr = memdup(xreal, size);
	xi = memdup(ximag, size);
	yr = memdup(yreal, size);
	yi = memdup(yimag, size);
	if (xr == NULL || xi == NULL || yr == NULL || yi == NULL)
		goto cleanup;
	
	if (!transform(xr, xi, n))
		goto cleanup;
	if (!transform(yr, yi, n))
		goto cleanup;
	for (i = 0; i < n; i++) {
		double temp = xr[i] * yr[i] - xi[i] * yi[i];
		xi[i] = xi[i] * yr[i] + xr[i] * yi[i];
		xr[i] = temp;
	}
	if (!inverse_transform(xr, xi, n))
		goto cleanup;
	for (i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
		outreal[i] = xr[i] / n;
		outimag[i] = xi[i] / n;
	}
	status = 1;
	
cleanup:
	free(yi);
	free(yr);
	free(xi);
	free(xr);
	return status;
}


static size_t reverse_bits(size_t x, unsigned int n) {
	size_t result = 0;
	unsigned int i;
	for (i = 0; i < n; i++, x >>= 1){
		result = (result << 1) | (x & 1);
	}
	return result;
}


static void *memdup(const void *src, size_t n) {
	void *dest = malloc(n);
	if (dest != NULL)
		memcpy(dest, src, n);
	return dest;
}

// Generate window function (Hanning)

void wHanning(float *w, int nSize) // size of the window
{
	int i;
	const double M = nSize-1;
	for (i = 0; i < nSize; i++) {
		w[i] = 0.5 * (1.0 - cos(2.0*PI*(float)i/M));
	}
	return;
}

// Generate window function (Hamming)
void wHamming(float *w, int nSize)
{
	int i;
	const float M = nSize-1;
	for (i = 0; i < nSize; i++) {
		w[i] = 0.54 - (0.46*cos(2.0*PI*(float)i/M));
	}
	return;
}

//////////////////////////////////////////////////////////
// convert frequency to mel 
//////////////////////////////////////////////////////////
double mel(float f){
	return 2595 * log10(1 + f/700); 
}

/////////////////////////////////////////////////////////
// Sorensen in-place radix-2 FFT for real values
// data: array of doubles:
// re(0),re(1),re(2),...,re(size-1)
// 
// output:
// re(0),re(1),re(2),...,re(size/2),im(size/2-1),...,im(1)
// normalized by array length
//
// Source: 
// Sorensen et al: Real-Valued Fast Fourier Transform Algorithms,
// IEEE Trans. ASSP, ASSP-35, No. 6, June 1987

void realfft_radix2(float *data, float *predata /*slide window*/,long n, float *window, float *dTotalPow, int nMode /*nMode = 0 ->FFT nMode = 1 ->Mel*/){
    float  xt,a,e, t1, t2, cc, ss;
    long  i, j, k, n1, n2, n3, n4, i1, i2, i3, i4;
	*dTotalPow = 0;
	for(i=0 ; i < n ; i++){
		//window overlapping
		if(i < n/2)
			data[i] = predata[i] * window[i];
		else 
			data[i] = data[i] * window[i];
	}
	//moving window
	memmove(&predata[0], &data[n/2], (n/2) * sizeof(float));
	n4=n-1;
  //data shuffling
      for (i=0,j=0,n2=n/2; i<n4 ; i++){
	  if (i<j){
				xt=data[j];
				data[j]=data[i];
				data[i]=xt;
				}
	  k=n2;
	  while (k<=j){
				j-=k;
				k>>=1;	
				}
	  j+=k;
      }
	
/* -------------------- */
    for (i=0; i<n; i += 2)  
      {
	 xt = data[i];
	 data[i] = xt + data[i+1];
	 data[i+1] = xt - data[i+1];
      }
/* ------------------------ */
    n2 = 1;
    for (k=n;k>2;k>>=1){ 
		n4 = n2;
		n2 = n4 << 1;
		n1 = n2 << 1;
		e = 2*PI/(n1);
		for (i=0; i<n; i+=n1){  
			xt = data[i];
			data[i] = xt + data[i+n2];
			data[i+n2] = xt-data[i+n2];
			data[i+n4+n2] = -data[i+n4+n2];
			a = e;
			n3=n4-1;
			for (j = 1; j <=n3; j++){
				i1 = i+j;
				i2 = i - j + n2;
				i3 = i1 + n2;
				i4 = i - j + n1;
				cc = cos(a);
				ss = sin(a);
				a += e;
				t1 = data[i3] * cc + data[i4] * ss;
				t2 = data[i3] * ss - data[i4] * cc;
				data[i4] = data[i2] - t2;
				data[i3] = -data[i2] - t2;
				data[i2] = data[i1] - t1;
				data[i1] += t1;
		  }
	  }
  }
	
	//division with array length
   for(i=0;i<n;i++){
	   if(nMode == 0)
		data[i] = abs(data[i])/n;
	   else {
		   data[i] = abs(data[i])/n;
		   //data[i] = 2595 * log10(1+data[i]/700);
	   }
	   *dTotalPow += data[i]; 
   }
	   //data[i]=(abs(data[i]))/n;
}




/////////////////////////////////////////////////////////
// Sorensen in-place inverse split-radix FFT for real values
// data: array of doubles:
// re(0),re(1),re(2),...,re(size/2),im(size/2-1),...,im(1)
// 
// output:
// re(0),re(1),re(2),...,re(size-1)
// NOT normalized by array length
//
// Source: 
// Sorensen et al: Real-Valued Fast Fourier Transform Algorithms,
// IEEE Trans. ASSP, ASSP-35, No. 6, June 1987

void irealfft_split(double *data,long n){

  long i,j,k,i5,i6,i7,i8,i0,id,i1,i2,i3,i4,n2,n4,n8,n1;
  double t1,t2,t3,t4,t5,a3,ss1,ss3,cc1,cc3,a,e,sqrt2; 
  sqrt2=sqrt(2.0);
  
  n1=n-1;
  n2=n<<1;
  for(k=n;k>2;k>>=1){  
	id=n2;
	n2>>=1;
	n4=n2>>2;
	n8=n2>>3;
	e = 2*PI/(n2);
	i1=0;
	do{ 
		for (; i1<n; i1+=id){
			i2=i1+n4;
			i3=i2+n4;
			i4=i3+n4;
			t1=data[i1]-data[i3];
			data[i1]+=data[i3];
			data[i2]*=2;
			data[i3]=t1-2*data[i4];
			data[i4]=t1+2*data[i4];
			if (n4!=1){
				i0=i1+n8;
				i2+=n8;
				i3+=n8;
				i4+=n8;
				t1=(data[i2]-data[i0])/sqrt2;
				t2=(data[i4]+data[i3])/sqrt2;
				data[i0]+=data[i2];
				data[i2]=data[i4]-data[i3];
				data[i3]=2*(-t2-t1);
				data[i4]=2*(-t2+t1);
			}
	     }
		 id<<=1;
	     i1=id-n2;
	     id<<=1;
	  } while ( i1<n1 );
	a=e;
	for (j=2; j<=n8; j++){  
	      a3=3*a;
	      cc1=cos(a);
	      ss1=sin(a);
	      cc3=cos(a3);
	      ss3=sin(a3);
	      a=j*e;
	      i=0;
	      id=n2<<1;
	      do{
		   for (; i<n; i+=id){  
			  i1=i+j-1;
			  i2=i1+n4;
			  i3=i2+n4;
			  i4=i3+n4;
			  i5=i+n4-j+1;
			  i6=i5+n4;
			  i7=i6+n4;
			  i8=i7+n4;
			  t1=data[i1]-data[i6];
			  data[i1]+=data[i6];
			  t2=data[i5]-data[i2];
			  data[i5]+=data[i2];
			  t3=data[i8]+data[i3];
			  data[i6]=data[i8]-data[i3];
			  t4=data[i4]+data[i7];
			  data[i2]=data[i4]-data[i7];
			  t5=t1-t4;
			  t1+=t4;
			  t4=t2-t3;
			  t2+=t3;
			  data[i3]=t5*cc1+t4*ss1;
			  data[i7]=-t4*cc1+t5*ss1;
			  data[i4]=t1*cc3-t2*ss3;
			  data[i8]=t2*cc3+t1*ss3;
			  }
		   id<<=1;
		   i=id-n2;
		   id<<=1;
		 } while(i<n1);
	   }
	}	

   /*----------------------*/
	i0=0;
	id=4;
   do{
       for (; i0<n1; i0+=id){ 
			i1=i0+1;
			t1=data[i0];
			data[i0]=t1+data[i1];
			data[i1]=t1-data[i1];
		}
	   id<<=1;
       i0=id-2;
       id<<=1;
    } while ( i0<n1 );

/*----------------------*/

//data shuffling
      for (i=0,j=0,n2=n/2; i<n1 ; i++){
	  if (i<j){
				t1=data[j];
				data[j]=data[i];
				data[i]=t1;
				}
	  k=n2;
	  while (k<=j){
				j-=k;
				k>>=1;	
				}
	  j+=k;
      }	
}

void preemphasize(float *fSample, int nSampleLen)  
{  
    /* Setting emphFac=0 turns off preemphasis. */   
    int i;   
    float emphFac = (float)0.9;   
   
    for (i = nSampleLen-1; i > 0; i--) {   
        fSample[i] = fSample[i] - emphFac * fSample[i-1];   
    }   
    fSample[0] = (float)(1.0 - emphFac) * fSample[0];     
}   

/////////////////////////////////////////////////////////
// mel ceptrum index 
/////////////////////////////////////////////////////////
#define MfccRound(x) ((int)((x)+0.5)) 
void filter_bank_int(float fStartFreq, int nSampleLen, float fSampleRate, int nNumOfCh, WfMelFB *melFb, float *fMelWeight){
	 int i, k;
	 int *nBinIdx;
	 float start_mel, fs_per_2_mel;
	 float *freq;
	 freq = malloc((nNumOfCh+2) * sizeof(float));
	 nBinIdx = malloc((nNumOfCh+2) * sizeof(int));
	 freq[0] = fStartFreq;
	 start_mel = (2595.0 * log10 (1.0 + fStartFreq / 700.0)); 
	 nBinIdx[0] = MfccRound(nSampleLen * freq[0] / fSampleRate); 
	 freq[nNumOfCh+1] = fSampleRate/8;
	 fs_per_2_mel = (2595.0 * log10 (1.0 + (fSampleRate /8.0) / 700.0));
	 nBinIdx[nNumOfCh+1] = MfccRound(nSampleLen * freq[nNumOfCh+1] / fSampleRate);  
	/* Calculating mel-scaled frequency and the corresponding FFT-bin */   
    /* number for the lower edge of the band                          */   
    for (k = 1; k <= nNumOfCh; k++) {   
        freq[k] = (float)(700 * (pow (10, (start_mel +  (float)k/(nNumOfCh + 1) * (fs_per_2_mel - start_mel)) / 2595.0) - 1.0));   
        nBinIdx[k] = MfccRound(nSampleLen * freq[k] / fSampleRate);
		//printf("Filter bank index %d is %f \r\n", k, freq[k]);		
    }
	for(i = 0; i < nBinIdx[0]; i++){   
        fMelWeight[i] = 0;   
    }   
    fMelWeight[nSampleLen/8]=1;   
	/* Initialize low, center, high indices to FFT-bin */   
    for (k = 0; k <= nNumOfCh; k++) {   
		if(k < nNumOfCh){   
			melFb[k].nLowIndx = nBinIdx[k];   
            melFb[k].nCentreIndx = nBinIdx[k+1];   
            melFb[k].nHighIndx = nBinIdx[k+2];  
			//printf("Filter bank center index %d is %d \r\n", k, melFb[k].nCentreIndx);			
        }   
        for(i = nBinIdx[k]; i < nBinIdx[k+1]; i++){   
            fMelWeight[i]=(i-nBinIdx[k]+1)/(float)(nBinIdx[k+1]-nBinIdx[k]+1); 
			//printf("Filter bank weidght %d is %f \r\n", i, fMelWeight[i]); 
        }   
    }
	free(freq);
	free(nBinIdx);
}

void filter_bank_out(float *fFFT, int nNumOfCh,float* fFilterOut, int normalize, WfMelFB *melFb, float *fMelWeight)   
{   
    float sum, wsum;   
    int i, k;     
       
    for (k=0;k<nNumOfCh;k++){      
        sum = fFFT[(&melFb[k])->nCentreIndx];   
        wsum=1;   
        for (i = (&melFb[k])->nLowIndx; i < (&melFb[k])->nCentreIndx; i++){   
            sum += fMelWeight[i] * pow(fFFT[i],2);   
            wsum += fMelWeight[i];   
        }   
        for (i = (&melFb[k])->nCentreIndx+1; i <= (&melFb[k])->nHighIndx; i++){   
            sum += (1 - fMelWeight[i-1]) * pow(fFFT[i],2);   
            wsum += (1 - fMelWeight[i-1]);   
        }   
        fFilterOut[k] = 0.5 * log(sum+1);
		//fFilterOut[k] = sum; 
				
        if(normalize) {    
            fFilterOut[k] /= wsum;  			
        }   
    }   
    return;   
}  

int dctInit(float *dctMatrix, int nNumOfCep, int nNumOfCh)   
{   
    int i, j;   
    for (i = 0; i <= nNumOfCep; i++){//12+1   
        for (j = 0; j < nNumOfCh; j++){//23   
            dctMatrix[i * nNumOfCh + j] = (float) cos (PI * (float) i / (float) nNumOfCh * ((float) j + 0.5));   
            if(i==0) dctMatrix[i * nNumOfCh + j]*=(float)sqrt(1/(float)nNumOfCh);   
            else     dctMatrix[i * nNumOfCh + j]*=(float)sqrt(2/(float)nNumOfCh); 
			//printf("DCT matrix is %f\r\n", dctMatrix[i * nNumOfCh + j]);
        }   
    }   
    return 1;   
}   

//direct cosine transform 
void dct(float *fFilterOut, float *dctMatrix, int nCepIndx, int nNumOfCh, float *mel_cep, float *delta_mel_cep, float *delta_delta_mel_cep)   
{   
    int i, j;   
	//memmove(&delta_delta_mel_cep[0], &delta_mel_cep[0], nCepIndx * sizeof(float)); 
	//memmove(&delta_mel_cep[0], &mel_cep[0], nCepIndx * sizeof(float)); 
    for (i = 0; i <= nCepIndx; i++) {   
        mel_cep[i] = 0.0;   
        for (j = 0; j < nNumOfCh; j++){   
            mel_cep[i] += fFilterOut[j] * dctMatrix[i * nNumOfCh + j];   
        }  
    }      
	return;   
}   

void mel_shift(float *mel_cep, float **mel_cep_frame, int nFrameNum, int nNumOfCep){
	int i, j;
	for(i = (nFrameNum-1); i >= 1; i--){
		for (j = 0; j < nNumOfCep; j++)
		{
			*(*(mel_cep_frame+i)+j) = *(*(mel_cep_frame+i-1)+j);
		}
	}
	
	for (j = 0; j < nNumOfCep; j++)
	{
		*(*(mel_cep_frame)+j) = mel_cep[j];
	}
	//*(mel_cep_frame) = mel_cep;
}

void vel_cep(float **mel_cep_frame, float *del_mel_cep, int nCepIndx){
	int i;
	for(i = 0; i < nCepIndx; i++){
		del_mel_cep[i] =  (2*(*(*(mel_cep_frame+4)+i) - (*(*(mel_cep_frame)+i))) + (*(*(mel_cep_frame+3)+i) - (*(*(mel_cep_frame+1)+i))));
		//printf("delta cepstrum is %f\r\n", del_mel_cep[i]);
		//printf("delta cepstrum is %f\r\n", (*(*(mel_cep_frame+1)+i)));
		//printf("delta cepstrum is %f\r\n", (*(*(mel_cep_frame+2)+i)));
		//printf("delta cepstrum is %f\r\n", (*(*(mel_cep_frame+3)+i)));
		//printf("delta cepstrum is %f\r\n", (*(*(mel_cep_frame+4)+i)));
	}		
	return;
}

//frame energy summation
float framEnergy(unsigned short *frameBUF, int nSize){
	int i;
	float fEnergy = 0;
	for(i = 0; i < nSize; i++)
		fEnergy += (float)pow(frameBUF[i],2) ;
	
	return fEnergy;
}






