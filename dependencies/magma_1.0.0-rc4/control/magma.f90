!
!   -- MAGMA (version 1.0) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      November 2010
!

module magma

  use magma_zfortran
  use magma_dfortran
  use magma_cfortran
  use magma_sfortran
  
  implicit none

  integer, parameter :: sizeof_complex_16 = 16
  integer, parameter :: sizeof_complex    = 8
  integer, parameter :: sizeof_double     = 8
  integer, parameter :: sizeof_real       = 4
  
end module magma
