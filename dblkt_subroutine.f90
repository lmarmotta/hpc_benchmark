program testando

    ! This program main goal is to test the block tridiagonal solver developed
    ! for CFD applications. There are a lot of good optimizations that can be
    ! done here but this subroutine works perfectly.

    ! Compile the code with gfortran [this file] -lblas -llapack
    ! Run ./a.out
    ! The screen output should be the correct answer hard coded below.

    implicit none

    integer(4) :: i,j
    real(8), dimension(2,2,3) :: A
    real(8), dimension(2,2,3) :: B,C
    real(8), dimension(2,3) :: xb
    real(8), dimension(2,3) :: x


    !
    ! Define the matrix to be tested. It's a block matrix.
    !

    do i = 1,3
        A(1,1,i) = 4.0d0
        A(1,2,i) = 3.0d0
        A(2,1,i) = 3.0d0
        A(2,2,i) = 2.0d0
    end do

    do i = 1,3

        B(1,1,i) = 8.0d0
        B(1,2,i) = 6.0d0
        B(2,1,i) = 6.0d0
        B(2,2,i) = 4.0d0

        C(1,1,i) = 16.0d0
        C(1,2,i) = 9.0d0
        C(2,1,i) = 9.0d0
        C(2,2,i) = 6.0d0
    end do


    !
    ! Define the B vector in AX = B system.
    !

    do i = 1, 2
        do j = 1, 3
            xb(i,j) = 3.0d0
        end do
    end do


    !
    ! Call solver.
    !

    call blktriad(A,B,C,2,3,xb,x)


    !
    ! Hard coded answer.
    !

    write(*,*) "x(1) = 5.8571429"
    write(*,*) "x(2) = -8.9220779"
    write(*,*) "x(3) = 0.5714286"
    write(*,*) "x(4) = -0.3116883"
    write(*,*) "x(5) = 1.8571429"
    write(*,*) "x(6) = -2.3766234"

    write(*,*) "Calculated..." 

    do j = 1, 3
        do i = 1, 2
            write(*,*) "x(",i,") = ",x(i,j)
        end do
    end do

end program testando

subroutine  blktriad(maind,lower,upper,id,md,xb,x)

     !|     B(1)    C(1)             | | x(1)  |
     !| A(2)  B(2)    C(2)           | |       |
     !|   A(3)  B(3)                 | |       |
     !|     .     .                  |*|       | = xb[1:mb*3]
     !|       .     .                | |       |
     !|               .       C(mb-1)| |       |
     !|         A(mb)  B(mb)         | |x(n*id)|

    ! id = inner matrices dimension.
    ! md = number matrices.
    ! maind = main diagonal of matrices  format: maind(id,id,md)
    ! lower = lower diagonal of matrices format: lower(id,id,2:md)
    ! upper = upper diagonal of matrices format: maind(id,id,md-1)
    ! xb    = B vector in Ax=B           format: xb(md*id)
    ! x     = x vector in Ax=B           format: x(md*id)
    !

    implicit none

    ! +++ Inputs +++
    !
    ! Scalar input variables.
    integer(kind=4) :: id,md

    ! Main diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: maind

    ! Lower diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: lower

    ! Upper diagonal of matrices.
    real(kind=8), dimension(id,id,md) :: upper

    ! Vector of equalties B in Ax=B.
    real(kind=8), dimension(id,md)    :: xb

    ! Vector of answers x in Ax=B.
    real(kind=8), dimension(id,md)    :: x

    ! ++ Inside variables ++
    !
    ! Scalar variables
    integer(kind=4) :: i,ii,jj

    ! Array of gamma coefficients.
    real(kind=8), dimension(id,id,md) :: gamm

    ! Array of beta coefficients.
    real(kind=8), dimension(id,md) :: beta

    ! Auxiliary arrays.
    real(kind=8), allocatable, dimension(:,:) :: aux_copy
    real(kind=8), allocatable, dimension(:,:) :: aux_mult
    real(kind=8), allocatable, dimension(:,:) :: aux_summ
    real(kind=8), allocatable, dimension(:)   :: aux_dumm

    allocate(aux_mult(id,id))
    allocate(aux_summ(id,id))
    allocate(aux_copy(id,id))
    allocate(aux_dumm(id))


    !--------------------------------------------------------------------------!
    !                  Step 1: BLOCK TRIANGULARIZATION                         !
    !--------------------------------------------------------------------------!


    ! Zero out the auxiliar vectors.

    beta = 0.0d0
    gamm = 0.0d0


    ! Lets first get our first gamma.

    aux_copy = maind(:,:,1)


    call inv(aux_copy,aux_copy,id)


    gamm(:,:,1) = matmul(aux_copy,upper(:,:,1))


    ! Now that we have our first gamma, lets get the rest of then.

    do i = 2, md-1

        aux_mult = 0.0d0
        aux_summ = 0.0d0

        aux_mult = matmul(lower(:,:,i),gamm(:,:,i-1))

        do jj = 1, id
            do ii = 1, id
                aux_summ(ii,jj) = maind(ii,jj,i) - aux_mult(ii,jj)
            end do
        end do

        call inv(aux_summ,aux_summ,id)

        gamm(:,:,i) = matmul(aux_summ,upper(:,:,i))
            
    end do

    
    ! Now that we have our gammas, lets get the betas, starting from the first
    ! ones. Note that now the calls done by the Lapack library will get a bit
    ! more complicated so lets use matmul...

    aux_copy = 0.0d0

    aux_copy = maind(:,:,1)

    call inv(aux_copy,aux_copy,id)

    beta(:,1) = matmul(aux_copy,xb(:,1))


    ! We now have our first beta, lets get the rest.

    do i = 2, md

        aux_mult = 0.0d0
        aux_summ = 0.0d0

        aux_mult(:,:) = matmul(lower(:,:,i),gamm(:,:,i-1))

        do jj = 1, id
            do ii = 1, id
                aux_summ(ii,jj) = maind(ii,jj,i) - aux_mult(ii,jj)
            end do
        end do

        call inv(aux_summ,aux_summ,id)

        aux_dumm(:) = xb(:,i) - matmul(lower(:,:,i),beta(:,i-1))

        beta(:,i) = matmul(aux_summ(:,:),aux_dumm(:))

    end do


    !--------------------------------------------------------------------------!
    !                  Step 2: BACKWARD SWEEP                                  !
    !--------------------------------------------------------------------------!


    ! How cool is that, lets start build our solution vector... iupiiii!

    x = 0.0d0

    x(:,md) = beta(:,md)

    do i = md-1,1,-1

        aux_dumm(:) = matmul(gamm(:,:,i),x(:,i+1))

        do ii = 1, id
            x(ii,i) = beta(ii,i) - aux_dumm(ii)
        end do

    end do


    deallocate(aux_mult)
    deallocate(aux_summ)
    deallocate(aux_copy)
    deallocate(aux_dumm)

end subroutine blktriad

subroutine inv(A,A_inv,m)

  Implicit none
  integer :: m
  real(8), dimension(m,m)::A, A_inv
  real(8),dimension(m)::WORK
  integer,dimension(m)::IPIV
  integer info

  A_inv = A

  call DGETRF(M,M,A_inv,M,IPIV,info)

  if (info /=  0) then
    write(*,*)"DGETRF: Failed during matrix factorization"
    stop
  end if

  call DGETRI(M,A_inv,M,IPIV,WORK,M,info)

  if (info /=  0) then
   write(*,*)"DGETRI: Failed during matrix inversion."
   stop
  end if

end subroutine inv

