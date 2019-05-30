!*********************************************************!
!*********************  ANT.G-2.4.0  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   David Jacob                                           !
!                                                         !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!                                                         !
!*********************************************************!
  MODULE AntCommon
!*********************************************************!
! module for sharing a *few* common variables among       !
! several modules                                         !
! i.e a kind of common block for ANT                      !
!*********************************************************!
!  use, intrinsic :: ISO_C_BINDING
  implicit none

PUBLIC :: C_F_STRING
  ! Job name from title section of Gaussian input file
  ! Used for output file names
  CHARACTER(len=50) :: jobname
!  character(jobname_len,kind=C_CHAR), intent(IN) :: jobname
  ! basic file name for ant1d input files
  CHARACTER(len=50) :: ant1dname
  
  ! *** Unique file unit identifiers ***
  ! To avoid chaos every input/output file has a unique identifier

  integer, parameter :: ifu_log = 6

  integer, parameter :: ifu_nam = 10
  integer, parameter :: ifu_ini = 11

  integer, parameter :: ifu_xyz = 101
  integer, parameter :: ifu_tra = 102
  integer, parameter :: ifu_red = 103
  integer, parameter :: ifu_dos = 104
  integer, parameter :: ifu_ham = 105
  integer, parameter :: ifu_diagham = 110
  integer, parameter :: ifu_mul = 106
  integer, parameter :: ifu_msh = 107
  integer, parameter :: ifu_hyb = 108
  integer, parameter :: ifu_ac  = 109

  integer, parameter :: ifu_bl  = 20

  integer, parameter :: ifu_dm  = 30
  integer, parameter :: ifu_dmx = 31
  integer, parameter :: ifu_fm  = 32

  integer, parameter :: ifu_ant = 40
  integer, parameter :: ifu_lead = 41


CONTAINS
  FUNCTION C_F_STRING(c_str) RESULT(f_str)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_CHAR
    TYPE(C_PTR), INTENT(IN) :: c_str
    CHARACTER(:,KIND=C_CHAR), POINTER :: f_str
    CHARACTER(KIND=C_CHAR), POINTER :: arr(:)
    INTERFACE
      ! Steal std C library function rather than writing our own.
      FUNCTION strlen(s) BIND(C, NAME='strlen')
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_SIZE_T
        IMPLICIT NONE
        !----
        TYPE(C_PTR), INTENT(IN), VALUE :: s
        INTEGER(C_SIZE_T) :: strlen
      END FUNCTION strlen
    END INTERFACE
    !****
    CALL C_F_POINTER(c_str, arr, [strlen(c_str)])
    CALL get_scalar_pointer(SIZE(arr), arr, f_str)
  END FUNCTION C_F_STRING
  SUBROUTINE get_scalar_pointer(scalar_len, scalar, ptr)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
    INTEGER, INTENT(IN) :: scalar_len
    CHARACTER(KIND=C_CHAR,LEN=scalar_len), INTENT(IN), TARGET :: scalar(1)
    CHARACTER(:,KIND=C_CHAR), INTENT(OUT), POINTER :: ptr
    !***
    ptr => scalar(1)
  END SUBROUTINE get_scalar_pointer

!function f_to_c_string(str,strlen) result(s)
!  use iso_c_binding
!  character(kind=c_char,len=:), allocatable :: s(*)
!  character(len=strlen), intent(in) :: str
!  integer i, nchars
!  allocate(character(len=strlen) :: s
!  s = transfer(str(1:nchars), s)
!end function f_to_c_string

function c_to_f_string(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  str = transfer(s(1:nchars), str)
end function c_to_f_string
 
!subroutine pstr(s) bind(c,name='pstr')
!  use iso_c_binding
!!  use cstr
!  character(kind=c_char,len=1), intent(in) :: s(*)
!  character(len=:), allocatable :: str
!  integer i, nchars
!  write(*,'(a)') c_to_f_string(s)
!end subroutine pstr


end module AntCommon
