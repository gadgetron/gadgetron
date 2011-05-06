/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/

#include "common_magma.h"

#if defined( _WIN32 )
  /* This must be included before INPUT is defined below, otherwise we
     have a name clash/problem  */
  #include <windows.h>
  #include <limits.h>
#else 
  #include <inttypes.h>
#endif

#if defined(ADD_)
#    define magma_gettime_f        magma_gettime_f_
#    define magma_gettimervalue_f  magma_gettimervalue_f_
#elif defined(NOCHANGE)
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Get current time
*/ 
extern "C"
TimeStruct get_current_time(void)
{
  static struct timeval  time_val;
  static struct timezone time_zone;

  TimeStruct time;

  cudaThreadSynchronize();
  gettimeofday(&time_val, &time_zone);

  time.sec  = time_val.tv_sec;
  time.usec = time_val.tv_usec;
  return (time);
}

extern "C"
void magma_gettime_f(unsigned int *time) {
  TimeStruct tmp = get_current_time();
  time[0] = tmp.sec;
  time[1] = tmp.usec;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- End elapsed time
*/ 
extern "C"
double GetTimerValue(TimeStruct time_1, TimeStruct time_2)
{
  int sec, usec;

  sec  = time_2.sec  - time_1.sec;
  usec = time_2.usec - time_1.usec;

  return (1000.*(double)(sec) + (double)(usec) * 0.001);
}

extern "C"
void magma_gettimervalue_f(unsigned int *start, unsigned int *end, double *result) {
  TimeStruct time1, time2;
  time1.sec  = start[0];
  time1.usec = start[1];
  time2.sec  = end[0];
  time2.usec = end[1];
  *result = GetTimerValue(time1, time2);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Print the available GPU devices
*/
extern "C"
void printout_devices( )
{
  int ndevices;
  cuDeviceGetCount( &ndevices );
  for( int idevice = 0; idevice < ndevices; idevice++ )
    {
      char name[200];
#if CUDA_VERSION > 3010 
      size_t totalMem;
#else
      unsigned int totalMem;
#endif

      int clock;
      CUdevice dev;

      cuDeviceGet( &dev, idevice );
      cuDeviceGetName( name, sizeof(name), dev );
      cuDeviceTotalMem( &totalMem, dev );
      cuDeviceGetAttribute( &clock,
                            CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev );
      printf( "device %d: %s, %.1f MHz clock, %.1f MB memory\n",
              idevice, name, clock/1000.f, totalMem/1024.f/1024.f );
    }
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Put 0s in the upper triangular part of a panel (and 1s on the diagonal)
      if uplo is 'U'/'u', or 0s in the lower triangular part of a panel (and 
      1s on the diagonal) if uplo is 'L'/'l'.
      This is auxiliary function used in geqrf and geqlf.  
*/
extern "C"
void spanel_to_q(char uplo, int ib, float *a, int lda, float *work){
  int i, j, k = 0;
  float *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<i; j++){
	work[k++] = col[j];
	col[j] = 0.;
      }
      work[k++] = col[i];
      col[j] = 1.;
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      work[k++] = col[i];
      col[i] = 1.;
      for(j=i+1; j<ib; j++){
        work[k++] = col[j];
        col[j] = 0.;
      }
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Restores a panel (after call to "panel_to_q").
      This isauxiliary function usedin geqrf and geqlf.
*/
extern "C"
void sq_to_panel(char uplo, int ib, float *a, int lda, float *work){
  int i, j, k = 0;
  float *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<=i; j++)
	col[j] = work[k++];
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=i; j<ib; j++)
        col[j] = work[k++];
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Put 0s in the upper triangular part of a panel (and 1s on the diagonal)
*/
extern "C"
void cpanel_to_q(char uplo, int ib, float2 *a, int lda, float2 *work){
  int i, j, k = 0;
  float2 *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<i; j++){
	work[k  ].x = col[j].x;
	work[k++].y = col[j].y;
	col[j].x = col[j].y = 0.;
      }
      work[k  ].x = col[i].x;
      work[k++].y = col[i].y;
      col[j].x = 1.;
      col[j].y = 0.;
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      work[k  ].x = col[i].x;
      work[k++].y = col[i].y;
      col[i].x = 1.;
      col[i].y = 0.;
      for(j=i+1; j<ib; j++){
        work[k  ].x = col[j].x;
	work[k++].y = col[j].y;
        col[j].x = col[j].y = 0.;
      }
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Restores a panel (after call to "panel_to_q")
*/
extern "C"
void cq_to_panel(char uplo, int ib, float2 *a, int lda, float2 *work){
  int i, j, k = 0;
  float2 *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<=i; j++){
	col[j].x = work[k  ].x;
	col[j].y = work[k++].y;
      }
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=i; j<ib; j++){
        col[j].x = work[k  ].x;
	col[j].y = work[k++].y;
      }
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Put 0s in the upper triangular part of a panel (and 1s on the diagonal)
*/
extern "C"
void dpanel_to_q(char uplo, int ib, double *a, int lda, double *work){
  int i, j, k = 0;
  double *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<i; j++){
	work[k++] = col[j];
	col[j] = 0.;
      }
      work[k++] = col[i];
      col[j] = 1.;
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      work[k++] = col[i];
      col[i] = 1.;
      for(j=i+1; j<ib; j++){
        work[k++] = col[j];
	col[j] = 0.;
      }
    }
  } 
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Restores a panel (after call to "panel_to_q")
*/
extern "C"
void dq_to_panel(char uplo, int ib, double *a, int lda, double *work){
  int i, j, k = 0;
  double *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<=i; j++)
	col[j] = work[k++];
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=i; j<ib; j++)
        col[j] = work[k++];
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Put 0s in the upper triangular part of a panel (and 1s on the diagonal)
*/
extern "C"
void zpanel_to_q(char uplo, int ib, cuDoubleComplex *a, int lda, cuDoubleComplex *work){
  int i, j, k = 0;
  cuDoubleComplex *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<i; j++){
	work[k  ].x = col[j].x;
	work[k++].y = col[j].y;
	col[j].x = col[j].y = 0.;
      }
      work[k  ].x = col[i].x;
      work[k++].y = col[i].y;
      col[j].x = 1.;
      col[j].y = 0.;
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      work[k  ].x = col[i].x;
      work[k++].y = col[i].y;
      col[i].x = 1.;
      col[i].y = 0.;
      for(j=i+1; j<ib; j++){
	work[k  ].x = col[j].x;
	work[k++].y = col[j].y;
	col[j].x = col[j].y = 0.;
      }
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Restores a panel (after call to "panel_to_q")
*/
extern "C"
void zq_to_panel(char uplo, int ib, cuDoubleComplex *a, int lda, cuDoubleComplex *work){
  int i, j, k = 0;
  cuDoubleComplex *col;

  if (uplo == 'U' || uplo == 'u'){
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=0; j<=i; j++){
	col[j].x = work[k  ].x;
	col[j].y = work[k++].y;
      }
    }
  }
  else {
    for(i=0; i<ib; i++){
      col = a + i*lda;
      for(j=i; j<ib; j++){
        col[j].x = work[k  ].x;
        col[j].y = work[k++].y;
      }
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Auxiliary function: ipiv(i) indicates that row i has been swapped with 
      ipiv(i) from top to bottom. This function rearranges ipiv into newipiv
      where row i has to be moved to newipiv(i). The new pivoting allows for
      parallel processing vs the original one assumes a specific ordering and
      has to be done sequentially.
*/
extern "C"
void swp2pswp(char trans, int n, int *ipiv, int *newipiv){
  int i, newind, ind;
  char            trans_[2] = {trans, 0};
  long int    notran = lapackf77_lsame(trans_, "N");

  for(i=0; i<n; i++)
    newipiv[i] = -1;
  
  if (notran){
    for(i=0; i<n; i++){
      newind = ipiv[i] - 1;
      if (newipiv[newind] == -1) {
	if (newipiv[i]==-1){
	  newipiv[i] = newind;
	  if (newind>i)
	    newipiv[newind]= i;
	}
	else
	  {
	    ind = newipiv[i];
	    newipiv[i] = newind;
	    if (newind>i)
	      newipiv[newind]= ind;
	  }
      }
      else {
	if (newipiv[i]==-1){
	  if (newind>i){
	    ind = newipiv[newind];
	    newipiv[newind] = i;
	    newipiv[i] = ind;
	  }
	  else
	    newipiv[i] = newipiv[newind];
	}
	else{
	  ind = newipiv[i];
	  newipiv[i] = newipiv[newind];
	  if (newind > i)
	    newipiv[newind] = ind;
	}
      }
    }
  } else {
    for(i=n-1; i>=0; i--){
      newind = ipiv[i] - 1;
      if (newipiv[newind] == -1) {
        if (newipiv[i]==-1){
          newipiv[i] = newind;
          if (newind>i)
            newipiv[newind]= i;
        }
        else
          {
            ind = newipiv[i];
            newipiv[i] = newind;
            if (newind>i)
              newipiv[newind]= ind;
          }
      }
      else {
        if (newipiv[i]==-1){
          if (newind>i){
            ind = newipiv[newind];
            newipiv[newind] = i;
            newipiv[i] = ind;
          }
          else
            newipiv[i] = newipiv[newind];
        }
        else{
          ind = newipiv[i];
          newipiv[i] = newipiv[newind];
          if (newind > i)
            newipiv[newind] = ind;
        }
      }
    }
  }
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Auxiliary function: used for debugging. Given a pointer to floating
      point number on the GPU memory, the function returns the value
      at that location.
*/
extern "C"
float getv(float *da){
  float res[1];
  cublasGetVector(1, sizeof(float), da, 1, res, 1);
  return res[0];
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Auxiliary function sp_cat
*/
extern "C"
int sp_cat(char *lp, char *rpp[], magma_int_t *rnp, magma_int_t*np, magma_int_t ll)
{
  magma_int_t i, n, nc;
  char *f__rp;

  n = (int)*np;
  for(i = 0 ; i < n ; ++i)
    {
      nc = ll;
      if(rnp[i] < nc)
        nc = rnp[i];
      ll -= nc;
      f__rp = rpp[i];
      while(--nc >= 0)
        *lp++ = *f__rp++;
    }
  while(--ll >= 0)
    *lp++ = ' ';

  return 0;
}
