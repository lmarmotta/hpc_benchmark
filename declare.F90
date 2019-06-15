module datatypes

    ! Default FORTRAN intrinsic variables.
    use, intrinsic :: iso_fortran_env

    ! Quarter precision integer variable.
    integer, parameter :: i_08 = int8 

    ! Half precision integer variable.
    integer, parameter :: i_16 = int16

    ! Single precision integer variable.
    ! Equivalent to integer(4)
    integer, parameter :: i_32 = int32

    ! Double precision integer variable.
    ! Equivalent to integer(8)
    integer, parameter :: i_64 = int64

    ! Single precision floating point.
    ! Equivalent to real(4)
    integer, parameter :: f_32 = real32

    ! Double precision floating point.
    ! Equivalent to real(8)
    integer, parameter :: f_64 = real64

    ! Extended precision floating point.
    ! Equivalent to real(16)
    integer, parameter :: f_128 = real128

end module datatypes

program main

    use datatypes
    implicit none

    integer(kind=i_08) :: a = 10
    integer(kind=i_16) :: b = 10
    integer(kind=i_32) :: c = 10
    integer(kind=i_64) :: d = 10

    real(kind=f_32) :: e  = 10.0d0
    real(kind=f_64) :: f  = 10.0d0
    real(kind=f_128) :: g = 10.0d0

    real(kind=8) :: h = 10.0d0

    write(*,*) a 
    write(*,*) b 
    write(*,*) c 
    write(*,*) d 
    write(*,*) e
    write(*,*) f
    write(*,*) h
    write(*,*) g

end program main
