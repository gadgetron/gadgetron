#include <stdlib.h>
#include <stdio.h>
#include "magma.h"

extern "C"
void magma_xerbla(const char *srname , magma_int_t *info)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

    Purpose
    =======

    XERBLA  is an error handler for the LAPACK routines.
    It is called by an LAPACK routine if an input parameter has an
    invalid value.  A message is printed and execution stops.

    Installers may consider modifying the STOP statement in order to
    call system-specific exception-handling facilities.

    Arguments
    =========

    SRNAME  (input) CHARACTER*(*)
            The name of the routine which called XERBLA.

    INFO    (input) INTEGER
            The position of the invalid parameter in the parameter list
            of the calling routine.

    =====================================================================   */

    printf("MAGMA Error: On Routine %s argument number %d had an illegal value\n",
	   srname , (-1)*(*info));
    exit(1);
}
