!
!   -- MAGMA (version 1.0) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      November 2010
!
!   @generated d
!

module magma_dfortran
  
  implicit none
  
  !---- Fortran interfaces to MAGMA subroutines ----
  interface
     
     subroutine magma_dgetrf_gpu(m, n, A, lda, ipiv, info)
       integer,          intent(in)     :: m, n, lda
       integer,          intent(out)    :: ipiv(*), info
       real, dimension(4), intent(inout):: A
     end subroutine magma_dgetrf_gpu
     
     subroutine magma_dgetrs_gpu(trans, n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       character,        intent(in)     :: trans
       integer,          intent(in)     :: n, nrhs, ldda, ipiv(*), lddb
       integer,          intent(out)    :: info
       real, dimension(4), intent(in)   :: dA
       real, dimension(4), intent(inout):: dB
     end subroutine magma_dgetrs_gpu
     
  end interface
end module magma_dfortran
