/**
 *
 * @file common_magma.h
 *
 *  MAGMA (version 1.0) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  November 2010
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11
 *
 * Based on PLASMA common.h
 *
 **/

/***************************************************************************//**
 *  MAGMA facilities of interest to both src and magmablas directories
 **/
#ifndef _MAGMA_COMMON_H_
#define _MAGMA_COMMON_H_

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include "magmawinthread.h"
#else
#include <pthread.h>
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

#include "magma.h"
#include "magma_lapack.h"
#include "operators.h"
#include "transpose.h"

/** ****************************************************************************
 *  Determine if weak symbol are allowed 
 */
#if defined(linux) || defined(__linux) || defined(__linux__)
#if defined(__GNUC_EXCL__) || defined(__GNUC__) 
#define MAGMA_HAVE_WEAK    1
#endif
#endif

/***************************************************************************//**
 *  Global utilities
 **/
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef roundup
#define roundup(a, b) (b <= 0) ? (a) : (((a) + (b)-1) & ~((b)-1))
#endif

/** ****************************************************************************
 *  Define magma_[sd]sqrt functions 
 *    - sqrt alone cannot be catched by the generation script because of tsqrt
 */

#define magma_dsqrt sqrt
#define magma_ssqrt sqrt

#endif /* _MAGMA_COMMON_H_ */
