/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/

#ifndef _MAGMABLAS_
#define _MAGMABLAS_

typedef int magma_int_t;

#include <cublas.h>
#include <cuda.h>

#include "magmablas_z.h"
#include "magmablas_c.h"
#include "magmablas_d.h"
#include "magmablas_s.h"
#include "magmablas_zc.h"
#include "magmablas_ds.h"

#define magmablas_zgemm cublasZgemm
#define magmablas_cgemm cublasCgemm

#endif
