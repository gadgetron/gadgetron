/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/

#include <stdio.h>
#include <string.h>

#include "magma.h"

extern "C" void 
auto_tune (char algorithm, char precision, magma_int_t ncores, magma_int_t ncorespsocket, 
           magma_int_t m, magma_int_t n, magma_int_t *nb, magma_int_t *ob, magma_int_t *ib, 
           magma_int_t *nthreads, magma_int_t *nquarkthreads)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======
    This function initializes tunable parameters to be used for 
    subsequent calls to hybrid routines in the MAGMA library.

    The idea is to use the matrix size together with the number of cores 
    and the number of cores per socket to do a table lookup for tunable 
    parameter values based on existing research results.

    Arguments
    =========
    algorithm     (input) CHAR
                  'q' QR
                  'l' LU
                  'c' Choleskey

    precision     (input) CHAR
                  's' Single
                  'd' Double
                  'c' Complex single
                  'z' Complex double

    ncores        (input) INTEGER
                  Number of cores
 
    ncorespsocket (intput) INTEGER
                  Number of cores per socket

    m             (input) INTEGER
                  Number of rows

    n             (intput) INTEGER
                  Number of columns

    nb            (output) INTEGER
                  Block size
 
    ob            (output) INTEGER
                  Outer block size

    ib            (output) INTEGER
                  Inner block size

    nthreads      (output) INTEGER
                  Number of MAMGMA threads
 
    nquarkthreads (output) INTEGER
                  Number of QUARK threads

    ===================================================================== */

  /* if QR */
  if (algorithm == 'q') {

    /* The best inner block size is always 12 */
    *ib = 12;

    /* The best number of QUARK threads is the number of cores per socket, in general */
    *nquarkthreads = ncorespsocket;
   
    /* 0 <= m <= 2080 */
    if ((m > 0) && (m <= 2080)) {

      *nb = 64;
      *ob = 64;

      *nthreads = 2;

    }

    /* 2080 < m <= 3360 */
    if ((m > 2080) && (m <= 3360)) {

      *nb = 128;
      *ob = 128;

      *nthreads = 6;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

    }

    /* 3360 < m <= 4640 */
    if ((m > 3360) && (m <= 4640)) {

      *nb = 128;
      *ob = 128;

      *nthreads = 14;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

    }
    
    /* 4640 < m <= 5920 */
    if ((m > 4640) && (m <= 5920)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 18;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 128;

        *nquarkthreads = 4;
        *nthreads = 8;

      }

    }

    /* 5920 < m <= 7200 */
    if ((m > 5920) && (m <= 7200)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 22;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 128;

        *nquarkthreads = 4;
        *nthreads = 8;

      }

    }

    /* 7200 < m <= 8480 */
    if ((m > 7200) && (m <= 8480)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 26;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 128;

        *nquarkthreads = 3;
        *nthreads = 9;

      }

    }

    /* 8480 < m <= 9760 */
    if ((m > 8480) && (m <= 9760)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 30;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 128;

        *nquarkthreads = 3;
        *nthreads = 9;

        if (precision == 's'){

          *nb = 192;
          *ob = 192;

        }

      }

    }

    /* 9760 < m <= 11040 */
    if ((m > 9760) && (m <= 11040)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 34;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 128;

        *nquarkthreads = 2;
        *nthreads = 10;

        if (precision == 's'){

          *nb = 192;
          *ob = 192;

        }

      }

    }

    /* 11040 < m <= 12320 */
    if ((m > 11040) && (m <= 12320)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 36;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 224;

        if (precision == 'c'){

          *ob = 128;

        }

        *nquarkthreads = 2;
        *nthreads = 10;

        if (precision == 's'){

          *nb = 192;
          *ob = 192;

        }

      }

    }

    /* 12320 < m <= 13600 */
    if ((m > 12320) && (m <= 13600)) {

      *nb = 128;
      *ob = 160;

      *nthreads = 42;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      if (precision == 'z'){

        *ob = 192;

      }

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 224;

        if (precision == 'c'){

          *ob = 128;

        }

        *nquarkthreads = 2;
        *nthreads = 10;

        if (precision == 's'){

          *nb = 192;
          *ob = 224;

        }

      }

    }

    /* 13600 < m <= 15220 */
    if ((m > 13600) && (m <= 15220)) {

      *nb = 128;
      *ob = 192;

      *nthreads = 42;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      if (precision == 'd'){

        *ob = 160;

      }

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 224;

        if (precision == 'c'){

          *ob = 160;

        }

        *nquarkthreads = 2;
        *nthreads = 10;

        if (precision == 's'){

          *nb = 192;
          *ob = 224;

        }

      }

    }

    /* 15220 < m <= 16800 */
    if ((m > 15220) && (m <= 16800)) {

      *nb = 128;
      *ob = 192;

      *nthreads = 42;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      if (precision == 'd'){

        *nb = 160;
        *ob = 200;

      }

      if (precision == 'c'){

        *ob = 160;

      }

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 224;

        if (precision == 'c'){

          *ob = 192;

        }

        *nquarkthreads = 2;
        *nthreads = 10;

        if (precision == 's'){

          *nb = 192;
          *ob = 224;

        }

      }

    }

    /* 16800 < m */
    if (m > 16800) {

      *nb = 128;
      *ob = 224;

      *nthreads = 42;
      if ((*nthreads + *nquarkthreads) > ncores)
        *nthreads = ncores - *nquarkthreads;

      if (precision == 'd'){

        *nb = 160;
        *ob = 200;

      }

      if (precision == 'c'){

        *ob = 192;

      }

      /* ncores = 12; ncorespsocket = 6 */
      if ((ncores == 12) && (ncorespsocket == 6)) {

        *ob = 256;

        if (precision == 'c'){

          *ob = 192;

        }

        *nquarkthreads = 2;
        *nthreads = 10;

        if (precision == 's'){

          *nb = 192;
          *ob = 224;

        }

      }

    }

  }
    
}

