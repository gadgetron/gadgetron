#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2010
#//////////////////////////////////////////////////////////////////////////////

#
# GPU_TARGET specifies for which GPU you want to compile MAGMA
#      0: Tesla family
#      1: Fermi Family
#
GPU_TARGET = 1

CC        = gcc
NVCC      = /usr/local/cuda/bin/nvcc
FORT      = gfortran

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

OPTS      = -DADD_ -O3 -fPIC
NVOPTS    = --compiler-options -fno-strict-aliasing -DUNIX -DADD_ -O3
LDOPTS    = -fPIC -Xlinker -zmuldefs

LIB       = -lf77blas -latlas -lcblas -lf2c -lcublas -lm -llapack

CUDADIR   = /usr/local/cuda

LIBDIR    = -L$(CUDADIR)/lib64
INC       = -I$(CUDADIR)/include

LIBMAGMA     = ../lib/libmagma.a
LIBMAGMABLAS = ../lib/libmagmablas.a
