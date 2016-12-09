program thomas

    ! This benchmark program compares LAPACK's tridiagonal solver with an in-hou
    ! se one. Since our usual usage involves solving system of equation with
    ! order greater than 1000x1000 matrices, the tests were done with this size
    ! of matrix.

    ! Author: Leonardo Motta

    implicit none

    integer(kind=4) :: i,INFO
    integer(kind=4) :: matrix_size  = 10000
    integer(kind=4) :: number_calls = 10000

    real(kind=8) :: t1,t2,tf

    real(kind=8), allocatable, dimension(:) :: main
    real(kind=8), allocatable, dimension(:) :: below
    real(kind=8), allocatable, dimension(:) :: upper
    real(kind=8), allocatable, dimension(:) :: bvect
    real(kind=8), allocatable, dimension(:) :: answr


    open(1,file="perf_dgtsv.dat")

    !--------------------------------------------------------------------------!
    !           TESTING OUT THE IN-HOUSE SCALAR THOMAS SUBROUTINE              !
    !--------------------------------------------------------------------------!


    ! Allocating vectors for the test.

    allocate(main(matrix_size))
    allocate(below(matrix_size))
    allocate(upper(matrix_size))
    allocate(bvect(matrix_size))
    allocate(answr(matrix_size))
    

    ! Fill up vectors for the test.

    call RANDOM_NUMBER(main)
    call RANDOM_NUMBER(upper)
    call RANDOM_NUMBER(below)
    call RANDOM_NUMBER(bvect)


    ! Zero out the answer vector.

    answr = 0.0d0



    call cpu_time(t1)

    do i = 1,number_calls
        call solve_tridiag(below,main,upper,bvect,answr,matrix_size)
    end do

    call cpu_time(t2)

    tf = t2 - t1

    write(1,'(A,I7,A,I7,A,F10.7,A)') " The time spent running the in-house & 
        thomas solver in ",number_calls," matrices with dimension: ", matrix_size," was:  ",tf," seconds."

    deallocate(main)
    deallocate(below)
    deallocate(upper)
    deallocate(bvect)
    deallocate(answr)


    !--------------------------------------------------------------------------!
    !            TESTING OUT LAPACK'S SCALAR THOMAS SUBROUTINE                 !
    !--------------------------------------------------------------------------!

    allocate(main(matrix_size))
    allocate(below(matrix_size-1))
    allocate(upper(matrix_size-1))
    allocate(bvect(matrix_size))


    call RANDOM_NUMBER(main)
    call RANDOM_NUMBER(upper)
    call RANDOM_NUMBER(below)
    call RANDOM_NUMBER(bvect)


    t1 = 0.0d0
    t2 = 0.0d0
    tf = 0.0d0


    call cpu_time(t1)

    do i = 1,number_calls

        call dgtsv(matrix_size,1,below,main,upper,bvect,matrix_size,INFO)

        if (INFO /= 0) then
            write(*,*) "Lapack call done improperlly come back to manual"
            stop
        end if

    end do

    call cpu_time(t2)

    tf = t2 - tf

    write(1,'(A,I7,A,I7,A,F10.7,A)') " The time spent running the LAPACK dgtsv & 
        thomas solver in ",number_calls," matrices with dimension: ", matrix_size," was:  ",tf," seconds."


    deallocate(main)
    deallocate(below)
    deallocate(upper)
    deallocate(bvect)

    close(1)

end program thomas

subroutine solve_tridiag(below,main,upper,bvect,answr,n)

    implicit none

    ! below - sub-diagonal (means it is the diagonal below the main diagonal)
    ! main - the main diagonal
    ! upper - sup-diagonal (means it is the diagonal above the main diagonal)
    ! bvectd - right part
    ! answr - the answer
    ! n - order of the matrix a

    integer(kind=4) :: n
    real(kind=8),dimension(n) :: below,main,upper,bvect
    real(kind=8),dimension(n) :: answr
    real(kind=8),dimension(n) :: gamm, beta
    real(kind=8) :: m
    integer(kind=4) :: i

    gamm(1) = upper(1) / main(1)

    do i = 2,n-1
        gamm(i) = upper(i) / (main(i)-(below(i)*gamm(i-1)))
    end do

    beta(1) = bvect(1)/main(1)

    do i = 2,n
        beta(i) = (bvect(i)-below(i)*beta(i-1))/(main(i) - below(i)*gamm(i-1))
    end do

    answr(n) = beta(n)

    do i = n-1,1,-1
        answr(i) = beta(i) - gamm(i)*answr(i+1)
    end do

end subroutine solve_tridiag

