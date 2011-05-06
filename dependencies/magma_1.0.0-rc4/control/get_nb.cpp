/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/

// ==== Definition of blocking sizes for Tesla ===============================
#if (GPUSHMEM < 200)

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for spotrf based on n
*/ 
extern "C"
int magma_get_spotrf_nb(int n){
  if (n <= 3328)
    return 128;
  else if (n<=4256)
    return 224;
  else 
    return 288;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgeqrf based on m
*/
extern "C"
int magma_get_sgeqrf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgetrf based on m;
      the return value should be multiple of 64
*/
extern "C"
int magma_get_sgetrf_nb(int m){
  if (m <= 2048)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for spotrf based on n
*/ 
extern "C"
int magma_get_dpotrf_nb(int n){
  if (n <= 3328)
    return 128;
  else if (n<=4256)
    return 128;
  else 
    return 256;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgeqrf based on m
*/
extern "C"
int magma_get_dgeqrf_nb(int m){
  if (m <= 2048)
    return 64;
  else 
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgeqlf based on m
*/
extern "C"
int magma_get_sgeqlf_nb(int m){
  if (m <= 1024)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgelqf based on m
*/
extern "C"
int magma_get_sgelqf_nb(int m){
  return magma_get_sgeqrf_nb(m);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgetrf based on m; 
      the return value should be multiple of 64
*/
extern "C"
int magma_get_dgetrf_nb(int m){
  if (m <= 2048)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgeqlf based on m
*/
extern "C"
int magma_get_dgeqlf_nb(int m){
  if (m <= 1024)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgelqf based on m
*/
extern "C"
int magma_get_dgelqf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_dgehrd_nb(int m){
  if (m <= 2048)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_sgehrd_nb(int m){
  if (m <= 1024)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_ssytrd_nb(int m){
  return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_sgebrd_nb(int m){
  return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_dgebrd_nb(int m){
  return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_dsytrd_nb(int m){
  return 32;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgeqrf based on m
*/
extern "C"
int magma_get_cgeqrf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgeqlf based on m
*/
extern "C"
int magma_get_cgeqlf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgeqrf based on m
*/
extern "C"
int magma_get_cgelqf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgetrf based on m;
      the return value should be multiple of 64
*/
extern "C"
int magma_get_cgetrf_nb(int m){
  if (m <= 2048)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cpotrf based on n
*/
extern "C"
int magma_get_cpotrf_nb(int n){
  return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zpotrf based on n
*/
extern "C"
int magma_get_zpotrf_nb(int n){
  return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for spotrf based on n
*/
extern "C"
int magma_get_zgetrf_nb(int n){
  return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgeqrf based on m
*/
extern "C"
int magma_get_zgeqrf_nb(int m){
  if (m <= 1024)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgeqlf based on m
*/
extern "C"
int magma_get_zgeqlf_nb(int m){
  if (m <= 1024)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for czeqlf based on m
*/
extern "C"
int magma_get_zgelqf_nb(int m){
  if (m <= 1024)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_zgehrd_nb(int m){
  if (m <= 2048)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_cgehrd_nb(int m){
  if (m <= 1024)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_chetrd_nb(int m){
  return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_zhetrd_nb(int m){
  return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_cgebrd_nb(int m){
  return 32;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_zgebrd_nb(int m){
  return 32;
}

// ==== End Definition of blocking sizes for Tesla ===========================
#else
// ====     Definition of blocking sizes for Fermi ===========================

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for spotrf based on n
*/ 
extern "C"
int magma_get_spotrf_nb(int n){
  if (n<=1024)
    return 160;
  else
    return 320;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgeqrf based on m
*/
extern "C"
int magma_get_sgeqrf_nb(int m){
  if (m<2000)
    return 96;
  else
    return 192;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgetrf based on m;
      the return value should be multiple of 64
*/
extern "C"
int magma_get_sgetrf_nb(int m){
  if (m<=3200)
    return 64;
  else
    return 192;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for spotrf based on n
*/ 
extern "C"
int magma_get_dpotrf_nb(int n){
  if (n<=4256)
    return 128;
  else 
    return 256;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgeqrf based on m
*/
extern "C"
int magma_get_dgeqrf_nb(int m){
  if (m <= 2048)
    return 64;
  else 
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgeqlf based on m
*/
extern "C"
int magma_get_sgeqlf_nb(int m){
  return magma_get_sgeqrf_nb(m);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgelqf based on m
*/
extern "C"
int magma_get_sgelqf_nb(int m){
  return magma_get_sgeqrf_nb(m);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgetrf based on m; 
      the return value should be multiple of 64
*/
extern "C"
int magma_get_dgetrf_nb(int m){
  if (m <= 2048)
    return 64;
  else if (m<7200)
    return 192;
  else
    return 256;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgeqlf based on m
*/
extern "C"
int magma_get_dgeqlf_nb(int m){
  return magma_get_dgeqrf_nb(m);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgelqf based on m
*/
extern "C"
int magma_get_dgelqf_nb(int m){
  return magma_get_dgeqrf_nb(m);
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_dgehrd_nb(int m){
  if (m <= 2048)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_sgehrd_nb(int m){
  if (m <= 1024)
    return 32;
  else
    return 96;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_ssytrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_sgebrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_dgebrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_dsytrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgeqrf based on m
*/
extern "C"
int magma_get_cgeqrf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgeqlf based on m
*/
extern "C"
int magma_get_cgeqlf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgeqrf based on m
*/
extern "C"
int magma_get_cgelqf_nb(int m){
  if (m <= 2048)
    return 32;
  else if (m<=4032)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgetrf based on m;
      the return value should be multiple of 64
*/
extern "C"
int magma_get_cgetrf_nb(int m){
  if (m <= 2048)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cpotrf based on n
*/
extern "C"
int magma_get_cpotrf_nb(int n){
  return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zpotrf based on n
*/
extern "C"
int magma_get_zpotrf_nb(int n){
  return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgetrf based on n
*/
extern "C"
int magma_get_zgetrf_nb(int n){
  if (n<=3072)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgeqrf based on m
*/
extern "C"
int magma_get_zgeqrf_nb(int m){
  if (m <= 1024)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgeqlf based on m
*/
extern "C"
int magma_get_zgeqlf_nb(int m){
  if (m <= 1024)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for czeqlf based on m
*/
extern "C"
int magma_get_zgelqf_nb(int m){
  if (m <= 1024)
    return 64;
  else
    return 128;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for dgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_zgehrd_nb(int m){
  if (m <= 2048)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for sgehrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_cgehrd_nb(int m){
  if (m <= 1024)
    return 32;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_chetrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for ssytrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_zhetrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for cgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_cgebrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgebrd based on m;
      the return value should be a multiple of 32
*/
extern "C"
int magma_get_zgebrd_nb(int m){
  //return 24;
  return 32;
  if (m <= 1024)
    return 64;
  else
    return 64;
}

// ==== End Definition of blocking sizes for Fermi ===========================
#endif
