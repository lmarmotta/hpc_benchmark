! ******************************************************************************
!
! This code implements a comprehensive and organized structure for CFD codes.
!
! ******************************************************************************

MODULE datatypes_module

    ! Default FORTRAN intrinsic variables.

    USE, INTRINSIC :: iso_fortran_env

    ! Quarter precision integer variable.

    INTEGER, PARAMETER :: i_08 = int8 

    ! Half precision integer variable.

    INTEGER, PARAMETER :: i_16 = int16

    ! Single precision integer variable.
    ! Equivalent to integer(4)

    INTEGER, PARAMETER :: i_32 = int32

    ! Double precision integer variable.
    ! Equivalent to integer(8)

    INTEGER, PARAMETER :: i_64 = int64

    ! Single precision floating point.
    ! Equivalent to real(4)

    INTEGER, PARAMETER :: f_32 = real32

    ! Double precision floating point.
    ! Equivalent to real(8)

    INTEGER, PARAMETER :: f_64 = real64

    ! Extended precision floating point.
    ! Equivalent to real(16)

    INTEGER, PARAMETER :: f_128 = real128

    ! Cell datatype

    TYPE cell
        REAL(KIND=f_64) :: vol
        REAL(KIND=f_64), DIMENSION(:), ALLOCATABLE :: nodes
    END TYPE cell

END MODULE datatypes_module

MODULE functions_module

    USE datatypes_module

    IMPLICIT NONE

    CONTAINS

        REAL(f_64) FUNCTION sum_numbers(a,b)

            IMPLICIT NONE

            REAL(KIND=f_64) :: a
            REAL(KIND=f_64) :: b

            sum_numbers = a + b

        END FUNCTION sum_numbers 

        REAL(f_64) FUNCTION mul_numbers(a,b)

            IMPLICIT NONE

            REAL(KIND=f_64) :: a
            REAL(KIND=f_64) :: b

            mul_numbers = a + b

        END FUNCTION  mul_numbers

END MODULE functions_module

PROGRAM main

    USE datatypes_module
    USE functions_module, ONLY : mul_numbers 

    IMPLICIT NONE

    ! Declare each integer type to test the lenth.

    INTEGER(KIND=i_08) :: a = 10
    INTEGER(KIND=i_16) :: b = 10
    INTEGER(KIND=i_32) :: c = 10
    INTEGER(KIND=i_64) :: d = 10

    INTEGER(KIND=i_08) :: i

    ! Declare each real type to test the lenth.

    REAL(KIND=f_32) :: e  = 10.0d0
    REAL(KIND=f_64) :: f  = 10.0d0
    REAL(KIND=f_128) :: g = 10.0d0

    REAL(KIND=f_64) :: h = 10.0d0

    TYPE(cell), DIMENSION(:), ALLOCATABLE :: volume_x

    ALLOCATE(volume_x(10))

    DO i =1, 10
        ALLOCATE(volume_x(i)%nodes(2))
    END DO

    WRITE(*,*) a 
    WRITE(*,*) b 
    WRITE(*,*) c 
    WRITE(*,*) d 
    WRITE(*,*) e
    WRITE(*,*) f
    WRITE(*,*) h
    WRITE(*,*) g

    WRITE(*,*) mul_numbers(f,f)

    volume_x(1)%vol = 3.0d0
    volume_x(2)%vol = 4.0d0
    volume_x(2)%nodes(1) = 4.0d0

    WRITE(*,*) volume_x(1)%vol
    WRITE(*,*) volume_x(2)%vol
    WRITE(*,*) volume_x(2)%nodes(2)

    DEALLOCATE(volume_x)

end program main
