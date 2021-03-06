! This is a POC to solve the overdefined system generated by least squares
! routine.
! It solves:
!   |  -0.57  -1.28  -0.39  |        | -2.67|
!   |  -1.93   1.08  -0.31  |        | -0.55|
! A |   2.30   0.24   0.40  |  x = B |  3.34|
!   |  -1.93   0.64  -0.66  |        | -0.77|
!   |   0.15   0.30   0.15  |        |  0.48|
!   |  -0.02   1.03  -1.43  |        |  4.10|

! x = | 1.5462|
!     | 1.8372| 
!     |-1.5642|
!
! Base code is from:
! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11
!        /mkl_lapack_examples/dgels_ex.f.htm
! 
! Example taken from:
! https://www.nag.com/lapack-ex/node45.html
!
! Result checkded with:
! http://calculator.vhex.net/calculator/linear-algebra/least-squares-solution
!
program wleast

    implicit none
!
!   Get the lhs matrix of weights.
!
    real(kind=8),dimension(6,3) :: w_rhs
!
!   Get the lhs vector of weights.
!
    real(kind=8),dimension(6) :: w_lhs
!
!   Put your shit together.
!
    w_rhs(:,:) = 0.0d0
    w_lhs(:)   = 0.0d0
!
!   Fill the first column of the matrix.
!
    w_rhs(1,1) = -0.57d0
    w_rhs(2,1) = -1.93d0
    w_rhs(3,1) =  2.30d0
    w_rhs(4,1) = -1.93d0
    w_rhs(5,1) =  0.15d0
    w_rhs(6,1) = -0.02d0
!
!   Fill the second column of the A matrix.
!
    w_rhs(1,2) = -1.28d0 
    w_rhs(2,2) =  1.08d0
    w_rhs(3,2) =  0.24d0
    w_rhs(4,2) =  0.64d0
    w_rhs(5,2) =  0.30d0
    w_rhs(6,2) =  1.03d0
!
!   Fill the third column of the A matrix.
!
    w_rhs(1,3) = -0.39d0
    w_rhs(2,3) = -0.31d0
    w_rhs(3,3) =  0.40d0
    w_rhs(4,3) = -0.66d0
    w_rhs(5,3) =  0.15d0
    w_rhs(6,3) = -1.43d0
!
!   Now, fill the first column of B matrix.
!
    w_lhs(1) = -2.67d0
    w_lhs(2) = -0.55d0
    w_lhs(3) =  3.34d0
    w_lhs(4) = -0.77d0
    w_lhs(5) =  0.48d0
    w_lhs(6) =  4.10d0
!
!
    call l_leasts(w_rhs(:,:),w_lhs(:),6,3)
!
    write(*,'(A,F10.4)') " * The correct result of first column is:  1.5462 | Yours is: ", w_lhs(1)
    write(*,'(A,F10.4)') " * The correct result of first column is:  1.8372 | Yours is: ", w_lhs(2)
    write(*,'(A,F10.4)') " * The correct result of first column is: -1.5642 | Yours is: ", w_lhs(3)

end program wleast

subroutine l_leasts(lhs,rhs,m_a,m_b)
!
!    Here are some tips about the variables.
!    rhs: The B matrix in Ax=B.
!    lhs: The A matrix in Ax=B.
!    ans: The x vector in Ax=B.
!    m_a: Number of "appartments" in your matrix building.
!    m_b: Number of "buildings"in your matrix condominium.
!
    implicit none
!
    real(kind=8), dimension(m_a,m_b) :: lhs
    real(kind=8), dimension(m_a) :: rhs
    real(kind=8), dimension(m_a) :: ans
    integer(kind=4) :: m_a
    integer(kind=4) :: m_b
!
!   LAPACK Specifics.
!
    integer(kind=4) :: NRHS = 1      ! Number of columns in RHS (B).
    integer(kind=4) :: LWMAX = 100   ! Number of columns in RHS (B).
    integer(kind=4) :: INFO,LWORK    ! Error handlers.
!
    real(kind=8), dimension(1) :: AUX_WORK
    real(kind=8), allocatable, dimension(:) :: WORK
!
!   Call lapack.
!
    LWORK = -1
    call dgels('N',m_a,m_b,NRHS,lhs(:,:),m_a,rhs(:),m_a,AUX_WORK, &
                LWORK,INFO)
!
    LWORK = MIN( LWMAX, INT( AUX_WORK( 1 ) ) )
!
    allocate(WORK(LWORK))
!
    call dgels('N',m_a,m_b,NRHS,lhs(:,:),m_a,rhs(:),m_a,WORK, &
                 LWORK,INFO)
!
!   Check for errors.
!
    if( info.gt.0 ) then
        write(*,*)'The diagonal element ',INFO,' of the triangular '
        write(*,*)'factor of A is zero, so that A does not have full '
        write(*,*)'rank; the least squares solution could not be'
        write(*,*)'computed.'
        stop
    end if
!
    deallocate(WORK)
! 
!   The answeres are in the lhs vector man !
!
end subroutine l_leasts
