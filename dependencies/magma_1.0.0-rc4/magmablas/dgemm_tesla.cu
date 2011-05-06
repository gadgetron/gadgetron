/*
    -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010
*/
#include "common_magma.h"
#include "commonblas_d.h"

#define magmablas_dgemm_tesla magmablas_dgemm

extern "C" void 
magmablas_dgemm_tesla( char TRANSA, char TRANSB, int m , int n , int k , 
                       double alpha, const double *A, int lda, 
                                     const double *B, int ldb,
                       double beta,        double *C, int ldc)
{
/*  -- MAGMA (version 1.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2010

   Purpose
   =======
 
   DGEMM  performs one of the matrix-matrix operations
 
      C := alpha*op( A )*op( B ) + beta*C,
 
   where  op( X ) is one of
 
      op( X ) = X   or   op( X ) = X',
 
   alpha and beta are scalars, and A, B and C are matrices, with op( A )
   an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 
   Parameters
   ==========
 
   TRANSA - CHARACTER*1.
            On entry, TRANSA specifies the form of op( A ) to be used in
            the matrix multiplication as follows:
 
               TRANSA = 'N' or 'n',  op( A ) = A.
 
               TRANSA = 'T' or 't',  op( A ) = A'.
 
               TRANSA = 'C' or 'c',  op( A ) = A'.
 
            Unchanged on exit.
 
   TRANSB - CHARACTER*1.
            On entry, TRANSB specifies the form of op( B ) to be used in
            the matrix multiplication as follows:
 
               TRANSB = 'N' or 'n',  op( B ) = B.
 
               TRANSB = 'T' or 't',  op( B ) = B'.
 
               TRANSB = 'C' or 'c',  op( B ) = B'.
 
            Unchanged on exit.
 
   M      - INTEGER.
            On entry,  M  specifies  the number  of rows  of the  matrix
            op( A )  and of the  matrix  C.  M  must  be at least  zero.
            Unchanged on exit.
 
   N      - INTEGER.
            On entry,  N  specifies the number  of columns of the matrix
            op( B ) and the number of columns of the matrix C. N must be
            at least zero.
            Unchanged on exit.
 
   K      - INTEGER.
            On entry,  K  specifies  the number of columns of the matrix
            op( A ) and the number of rows of the matrix op( B ). K must
            be at least  zero.
            Unchanged on exit.
 
   ALPHA  - DOUBLE PRECISION.
            On entry, ALPHA specifies the scalar alpha.
            Unchanged on exit.
 
   A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
            k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
            Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
            part of the array  A  must contain the matrix  A,  otherwise
            the leading  k by m  part of the array  A  must contain  the
            matrix A.
            Unchanged on exit.
 
   LDA    - INTEGER.
            On entry, LDA specifies the first dimension of A as declared
            in the calling (sub) program. When  TRANSA = 'N' or 'n' then
            LDA must be at least  max( 1, m ), otherwise  LDA must be at
            least  max( 1, k ).
            Unchanged on exit.
 
   B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
            n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
            Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
            part of the array  B  must contain the matrix  B,  otherwise
            the leading  n by k  part of the array  B  must contain  the
            matrix B.
            Unchanged on exit.
 
   LDB    - INTEGER.
            On entry, LDB specifies the first dimension of B as declared
            in the calling (sub) program. When  TRANSB = 'N' or 'n' then
            LDB must be at least  max( 1, k ), otherwise  LDB must be at
            least  max( 1, n ).
            Unchanged on exit.
 
   BETA   - DOUBLE PRECISION.
            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
            supplied as zero then C need not be set on input.
            Unchanged on exit.
 
   C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
            Before entry, the leading  m by n  part of the array  C must
            contain the matrix  C,  except when  beta  is zero, in which
            case C need not be set on entry.
            On exit, the array  C  is overwritten by the  m by n  matrix
            ( alpha*op( A )*op( B ) + beta*C ).
 
   LDC    - INTEGER.
            On entry, LDC specifies the first dimension of C as declared
            in  the  calling  (sub)  program.   LDC  must  be  at  least
            max( 1, m ).
            Unchanged on exit.
   =====================================================================    */

   if(m==0 || n==0  || ( ( alpha==0 || k==0 ) && beta ==1 ) ){
	return ;
   }
   TRANSA = toupper( TRANSA ) ; 	
   TRANSB = toupper( TRANSB ) ; 	
   if( alpha == 0.0){
	if( beta == 0.0){
		magmablas_dgemm_kernel_ab_0( C,A,B, m, n,k,lda,ldb, ldc, alpha, beta);
		return ;
	}	
	else{
		magmablas_dgemm_kernel_a_0( C,A,B, m, n,k,lda,ldb, ldc, alpha, beta);
		return ;
	}		
   }
 
   if(ldc < m ) return ;
   if(TRANSA=='N'){
	if(TRANSB=='N')
	   { 
		if(lda < m ) return ;
		if(ldb < k ) return ;
		/*=======================================================================
		  ===================C = alpha * A * B + beta * C =======================
		  =======================================================================*/
		if( m > 512 && n > 512 ){
			if( m % 64 == 0 && n%16 == 0 && k%16 == 0 ) 
				magmablas_dgemm_kernel_N_N_64_16_16_16_4_special( C,A,B, 
                                                        m, n,k,lda,ldb, ldc, alpha, beta);
			else
				magmablas_dgemm_kernel_N_N_64_16_16_16_4(C,A,B, m, n, k,
                                                               lda,ldb, ldc, alpha, beta);
		}
		else{
			cublasDgemm(TRANSA, TRANSB, m, n, k, 
                                    alpha, A, lda, B, ldb, beta, C, ldc );
		}
	   }
	   else{ 
		if(lda < m ) return ;
		if(ldb < n ) return ;
		/*=======================================================================
		  ===================C = alpha * A * B^T + beta * C======================
		  =======================================================================*/
		if( m > 512 && n > 512 ){
     		        if( m%64 == 0 && n %16 ==0 && k%4==0) 
				magmablas_dgemm_kernel_N_T_64_16_4_16_4(C,A,B, m, n, k,
                                                                lda,ldb, ldc, alpha, beta);
			else 
				magmablas_dgemm_kernel_N_T_64_16_4_16_4(C,A,B, m, n, k,
                                                                lda,ldb, ldc, alpha, beta);
		}
		else{
			cublasDgemm(TRANSA, TRANSB, m, n, k, 
                                    alpha, A, lda, B, ldb, beta, C, ldc );
		}
           }
	}
        else{
	    if(TRANSB=='N'){
		if(lda < k ) return ;
		if(ldb < k ) return ;
		/*=======================================================================
		  ===================C = alpha * A^T * B + beta * C======================
		  =======================================================================*/
		if(m>512 && n > 512){
     			if( m%32 == 0 && n %32 ==0 && k%8==0) 
				magmablas_dgemm_kernel_T_N_32_32_8_8_8(C,A,B, m, n, k,
                                                                lda,ldb, ldc, alpha, beta);
	
			else
				magmablas_dgemm_kernel_T_N_32_32_8_8_8(C,A,B, m, n, k,
                                                                lda,ldb, ldc, alpha, beta);
		}
		else{
	       	        cublasDgemm(TRANSA, TRANSB, m, n, k,
                                    alpha, A, lda, B, ldb, beta, C, ldc );
		}	
            }
            else{
		if(lda < k) return ;
		if(ldb < n ) return ;
		/*=======================================================================
		  ===================C = alpha * A^T* B^T + beta * C=====================
		  =======================================================================*/
		if( m > 512 && n > 512 ){
	     		if( m%64 == 0 && n %16 ==0 && k%16==0) 
			      magmablas_dgemm_kernel_T_T_64_16_16_16_4_v2(C,B,A, n, m, k,
                                                                ldb,lda, ldc, alpha, beta);
			else 
			      magmablas_dgemm_kernel_T_T_64_16_16_16_4(C,B,A, n, m, k,
                                                                ldb,lda, ldc, alpha, beta);
		}
		else{
		        cublasDgemm(TRANSA, TRANSB, m, n, k, 
                                    alpha, A, lda, B, ldb, beta, C, ldc );
		}

	   }
        }
        
}
