program dgemm_bench


    ! This simple program aims to compute the matrix multiplication operation
    ! using lapack's dgemm in a big square matrix and in a series of 5x5 squa
    ! re matrices.
    ! Test with: -mkl and -llapack -lblas implementations.
    !
    ! Author: Leonardo Motta 20/09/2016.

    implicit none


    ! General benchmark variables.

    real(kind=8) :: t1
    real(kind=8) :: t2
    real(kind=8) :: tf

    integer(kind=4) :: matrix_size = 2500
    integer(kind=4) :: i


    ! Variables for the big matrix multiplication.

    real(kind=8), allocatable, dimension(:,:) :: matrix_A
    real(kind=8), allocatable, dimension(:,:) :: matrix_B
    real(kind=8), allocatable, dimension(:,:) :: matrix_R


    ! Output performance file.

    open(1,file="perf_dgemm.dat")

    !--------------------------------------------------------------------------!
    !          TESTING OUT THE MULTIPLICATION OF A SINGLE BIG MATRIX           !
    !--------------------------------------------------------------------------!

    ! Allocating matrices.

    t1 = 0.0d0
    t2 = 0.0d0
    tf = 0.0d0

    allocate(matrix_A(matrix_size,matrix_size))
    allocate(matrix_B(matrix_size,matrix_size))
    allocate(matrix_R(matrix_size,matrix_size))


    ! Putting randon numbers in the matrices.

    call RANDOM_NUMBER(matrix_A)
    call RANDOM_NUMBER(matrix_B)


    ! Zero out answer.

    matrix_R = 0.0d0


    ! Calling clock timer.

    call cpu_time(t1)


    ! Calling Lapack's DGEMM subroutine.

    call dgemm('N','N',matrix_size,matrix_size,matrix_size,1.0d0, &
        matrix_A,matrix_size,matrix_B,matrix_size,1.0d0,matrix_R,matrix_size)


    call cpu_time(t2)

    tf = t2-t1

    write(1,'(A,I7,A,F10.7)') " Elapse time running Lapack on ",matrix_size," size matrix = ", tf


    ! Testing matmul.

    t1 = 0.0d0
    t2 = 0.0d0
    tf = 0.0d0

    matrix_R = 0.0d0

    call cpu_time(t1)


    ! Calling matmul internal routine.

    matrix_R = matmul(matrix_A, matrix_B)

    call cpu_time(t2)

    tf = t2 - t1

    write(1,'(A,I7,A,F10.7)') " Elapse time running matmul on ",matrix_size," size matrix = ", tf


    ! Finish with the single big matrix multiplication.

    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(matrix_R)


    !--------------------------------------------------------------------------!
    !          TESTING OUT THE MULTIPLICATION MANY SMALL MATRICES              !
    !--------------------------------------------------------------------------!

    ! Allocating matrices.

    allocate(matrix_A(5,5))
    allocate(matrix_B(5,5))
    allocate(matrix_R(5,5))


    ! Putting random numbers in the matrice.

    call RANDOM_NUMBER(matrix_A)
    call RANDOM_NUMBER(matrix_B)


    t1 = 0.0d0
    t2 = 0.0d0
    tf = 0.0d0


    call cpu_time(t1)


    do i = 1,matrix_size*matrix_size
        
        call dgemm('N','N',5,5,5,1.0d0, &
            matrix_A,5,matrix_B,5,1.0d0,matrix_R,5)

    end do


    call cpu_time(t2)

    tf = t2-t1

    write(1,'(A,I7,A,F10.7)') " Elapse time running lapack on ",matrix_size," small 5x5 matrices = ", tf


    ! Do the same with matmul.

    ! Putting random numbers in the matrice.

    call RANDOM_NUMBER(matrix_A)
    call RANDOM_NUMBER(matrix_B)


    t1 = 0.0d0
    t2 = 0.0d0
    tf = 0.0d0


    call cpu_time(t1)


    do i = 1,matrix_size*matrix_size

        matrix_R = matmul(matrix_A, matrix_B)

    end do


    call cpu_time(t2)

    tf = t2-t1

    write(1,'(A,I7,A,F10.7)') " Elapse time running matmul on ",matrix_size," small 5x5 matrices = ", tf


    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(matrix_R)


    close(1)

end program dgemm_bench
