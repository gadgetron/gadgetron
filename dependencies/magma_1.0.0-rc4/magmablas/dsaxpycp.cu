/*
  -- MAGMA (version 1.0) --
  Univ. of Tennessee, Knoxville
  Univ. of California, Berkeley
  Univ. of Colorado, Denver
  November 2010

  @generated ds

*/
#include "common_magma.h"

extern "C" __global__ void
dsaxpycp_special(float *R, double *X, magma_int_t M, double *B,double *W )
{
    const magma_int_t ibx = blockIdx.x * 64;
    const magma_int_t idt = threadIdx.x;
    X += ibx+idt;
    R += ibx+idt;
    B += ibx+idt;
    W += ibx+idt;
    X[0] = MAGMA_D_ADD( X[0], (double)(R[0]) );
    W[0] = B[0];
}

extern "C" __global__ void
daxpycp_special(double *R, double *X, magma_int_t M, double *B)
{
    const magma_int_t ibx = blockIdx.x * 64;
    const magma_int_t idt = threadIdx.x;
    X += ibx+idt;
    R += ibx+idt;
    B += ibx+idt;
    X[0] = MAGMA_D_ADD( X[0], R[0] );
    R[0] = B[0];
}

extern "C" __global__ void
dsaxpycp_generic(float *R, double *X, magma_int_t M, double *B,double *W )
{
    const magma_int_t ibx = blockIdx.x * 64;
    const magma_int_t idt = threadIdx.x;
    if( ( ibx + idt ) < M ) {
        X += ibx+idt;
        R += ibx+idt;
        B += ibx+idt;
        W += ibx+idt;
    }
    else{
        X +=(M-1);
        R +=(M-1);
        B +=(M-1);
        W +=(M-1);
    }
    X[0] = MAGMA_D_ADD( X[0], (double)( R[0] ) );
    W[0] = B[0];
}

extern "C" __global__ void
daxpycp_generic(double *R, double *X, magma_int_t M, double *B)
{
    const magma_int_t ibx = blockIdx.x * 64;
    const magma_int_t idt = threadIdx.x;
    if( ( ibx + idt ) < M ) {
        X += ibx+idt;
        R += ibx+idt;
        B += ibx+idt;
    }
    else{
        X +=(M-1);
        R +=(M-1);
        B +=(M-1);
    }
    X[0] = MAGMA_D_ADD( X[0], R[0] );
    R[0] = B[0];
}


extern "C" void
magmablas_dsaxpycp(float *R, double *X, magma_int_t M, double *B, double *W)
{
    dim3 threads( 64, 1 );
    dim3 grid(M/64+(M%64!=0),1);
    if( M %64 == 0 ) {
        dsaxpycp_special <<< grid, threads >>> ( R, X, M, B, W) ;
    }
    else{
        dsaxpycp_generic <<< grid, threads >>> ( R, X, M, B, W) ;
    }
}

extern "C" void
magmablas_daxpycp(double *R, double *X, magma_int_t M, double *B)
{
    dim3 threads( 64, 1 );
    dim3 grid(M/64+(M%64!=0),1);
    if( M %64 == 0 ) {
        daxpycp_special <<< grid, threads >>> ( R, X, M, B) ;
    }
    else{
        daxpycp_generic <<< grid, threads >>> ( R, X, M, B) ;
    }
}
