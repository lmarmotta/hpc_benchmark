!
! This subroutine solves an eigenvalue problem which can be of complex type. 
! Valid for a square matrix.
!
subroutine eig(matrix,size_matrix,eigs)
!
    implicit none
!
!   inputs
!
    integer(kind=4), intent(in) :: size_matrix
    complex(kind=8), dimension(size_matrix) :: eigs
    complex(kind=8), dimension(size_matrix,size_matrix) :: matrix
!
!   parameters
!
    integer(kind=4) :: n
    integer(kind=4) :: lda
    integer(kind=4) :: ldvl
    integer(kind=4) :: ldvr
    integer(kind=4) :: lwmax
!
    integer(kind=4) :: i
!
!   local scalars
!
    integer(kind=4) :: info, lwork
!
!   local arrays
!
    real(kind=8), allocatable, dimension(:) :: rwork
!
!   Declare arrays.
!
    complex(kind=8), allocatable, dimension(:,:) :: vl
    complex(kind=8), allocatable, dimension(:,:) :: vr
    complex(kind=8), allocatable, dimension(:)   :: w
    complex(kind=8), allocatable, dimension(:)   :: work
!
!   Define the needed parameters.
!
    n     = size_matrix
    lda   = n
    ldvl  = n
    ldvr  = n
    lwmax = 1000
!
!   Allocate the working matrices.
!
    allocate(vl(ldvl, n))
    allocate(vr(ldvr, n))
    allocate(w(n))
    allocate(work(lwmax))
    allocate(rwork(2*N))
!
!   query the optimal workspace.
!
    lwork = -1
!
    call zgeev( 'vectors', 'vectors', n, matrix, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
!
!   Compute the proper size for lwork.
!
    lwork = min( lwmax, int( work( 1 ) ) )
!
!   solve eigenproblem.
!
    call zgeev( 'vectors', 'vectors', n, matrix, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
!
!   check for convergence.
!
    if( info.gt.0 ) then
        write(*,*)'the algorithm failed to compute eigenvalues.'
        stop
    end if
!
!   Fill up the answer.
!
    do i = 1, size_matrix
        eigs(i) = w(i)
    end do
!
!   No memory leaks.
!
    deallocate(vl)
    deallocate(vr)
    deallocate(w)
    deallocate(work)
    deallocate(rwork)
!
end subroutine eig

program eigenvalue
!
    implicit none
!
    complex(kind=8) :: eigs(2)
    complex(kind=8) :: a(2,2)
!
    a(1,1) = ( 3.0, 0.0) 
    a(1,2) = (-2.0, 0.0)
    a(2,1) = ( 4.0, 0.0)
    a(2,2) = (-1.0, 0.0)
!
!   Solve the eigenvalue problem.
!
    call eig(a,2,eigs)
!
!   Dump the solution.
!
    write(*,*) eigs(1)
    write(*,*) eigs(2)
!
end program eigenvalue

