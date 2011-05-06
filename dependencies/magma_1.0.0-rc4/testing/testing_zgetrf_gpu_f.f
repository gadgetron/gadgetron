!
!   -- MAGMA (version 1.0) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      November 2010
!
!  @precisions normal z -> c d s
!

      program testing_zgetrf_gpu_f

      use magma

      external cublas_init, cublas_set_matrix, cublas_get_matrix
      external cublas_shutdown, cublas_alloc
      external zlange, zgemm, zgesv

      double precision zlange
      integer cublas_alloc

      double precision              :: rnumber(2), Anorm, Bnorm, Rnorm 
      double precision, allocatable :: work(:)
      complex*16, allocatable       :: h_A(:), h_B(:), h_X(:)
      complex*16, allocatable       :: h_A2(:)
      real, dimension(4)            :: devptrA, devptrB
      integer,    allocatable       :: ipiv(:)

      complex*16                    :: zone, mzone
      integer                       :: i, n, info, stat, lda
      integer                       :: size_of_elt, nrhs
      real(kind=8)                  :: flops, t
      integer                       :: tstart(2), tend(2)

      PARAMETER          ( nrhs = 1, zone = 1., mzone = -1. )
      
      call cublas_init()

      n   = 2048
      lda  = n
      ldda = ((n+31)/32)*32
      size_of_elt = sizeof_complex_16
 
!------ Allocate CPU memory
      allocate(h_A(lda*n))
      allocate(h_A2(n*n))
      allocate(h_B(lda*nrhs))
      allocate(h_X(lda*nrhs))
      allocate(work(n))
      allocate(ipiv(n))

!------ Allocate GPU memory
      stat = cublas_alloc(ldda*n, size_of_elt, devPtrA)
      if (stat .ne. 0) then
         write(*,*) "device memory allocation failed"
         stop
      endif

      stat = cublas_alloc(ldda*nrhs, size_of_elt, devPtrB)
      if (stat .ne. 0) then
         write(*,*) "device memory allocation failed"
         stop
      endif

!---- Initializa the matrix
      do i=1,n*n
        call random_number(rnumber)
        h_A(i) = rnumber(1)
        h_A2(i) = h_A(i)
      end do

      do i=1,n*nrhs
        call random_number(rnumber)
        h_B(i) = rnumber(1)
        h_X(i) = h_B(i)
      end do

!---- devPtrA = h_A
      call cublas_set_matrix(n, n, size_of_elt, h_A, lda, devptrA, ldda)

!---- devPtrB = h_B
      call cublas_set_matrix(n, nrhs, size_of_elt, h_B, lda, devptrB, 
     $                       ldda)

!---- Call magma LU ----------------
      call magma_gettime_f(tstart)
      call magma_zgetrf_gpu(n, n, devptrA, ldda, ipiv, info)
      call magma_gettime_f(tend)

      if ( info .ne. 0 )  then
         write(*,*) "Info : ", info
      end if

!---- Call magma solve -------------
      call magma_zgetrs_gpu('n', n, nrhs, devptrA, ldda, ipiv, devptrB, 
     $                      ldda, info)

      if ( info .ne. 0 )  then
         write(*,*) "Info : ", info
      end if

!---- h_X = devptrB
      call cublas_get_matrix (n, nrhs, size_of_elt, devptrB, ldda, h_X, 
     $                        lda)

!---- Solve using LAPACK ----------------------------------------------
c      call zgesv(n, nrhs, h_A2, lda, ipiv, h_X, lda, info)

!---- Compare the two results ------
      Anorm = zlange('I', n, n,    h_A, lda, work)
      Bnorm = zlange('I', n, nrhs, h_B, lda, work)

      call zgemm('n', 'n', n,  nrhs, n, zone, h_A, lda, h_X, lda,
     $           mzone, h_B, lda)
      Rnorm = zlange('I', n, nrhs, h_B, lda, work)

      write(*,*)
      write(*,*  ) 'Solving A x = b using LU factorization:'
      write(*,105) '  || A || = ', Anorm
      write(*,105) '  || b || = ', Bnorm
      write(*,105) '  || b - A x || / (||A|| ||b||) = ', 
     $                Rnorm/(Anorm*Bnorm)

      flops = 2. * n * n * n / 3.  
      call magma_gettimervalue_f(tstart, tend, t)

      write(*,*)   '  Gflops  = ',  flops / t / 1e6
      write(*,*)

!---- Free CPU memory
      deallocate(h_A, h_A2, h_X, h_B, work, ipiv)

!---- Free GPU memory
      call cublas_free(devPtrA)
      call cublas_free(devPtrB)
      call cublas_shutdown()

 105  format((a35,es10.3))

      end
