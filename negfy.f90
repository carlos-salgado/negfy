

module linked_list_mod
   implicit none
   type node
!      integer, pointer :: array(:) => NULL()
      real, pointer :: array(:) => NULL()
      type(node), pointer :: next => NULL()
   end type node
end module linked_list_mod

module get_record_mod
   use linked_list_mod
   implicit none
   type char_node
      character c
      type(char_node), pointer :: next => NULL()
   end type char_node
   contains
      subroutine get_record(iunit,temp)
         integer, intent(in) :: iunit
         type(node) :: temp
         integer numc
         character c
         integer iostat
         type(char_node), pointer :: head
         type(char_node), pointer :: cursor

         numc = 0
         nullify(head)
         do
            read(iunit,'(a)',advance='no',iostat=iostat) c
            if(iostat /= 0) exit
            if(numc == 0) then
               allocate(head)
               cursor => head
            else
               allocate(cursor%next)
               cursor => cursor%next
            end if
            numc = numc+1
            cursor%c = c
         end do
         call gr_1(head,numc,temp%array)
      end subroutine get_record

      subroutine gr_1(head,numc,array)
         type(char_node), pointer :: head
         integer, intent(in) :: numc
!         integer, pointer :: array(:)
         real, pointer :: array(:)
         character(numc) line
         integer i
         type(char_node), pointer :: cursor
         character(3), parameter :: set = ', '//achar(9)
         integer numi
         integer deltai

         do i = 1, numc
            cursor => head
            head => head%next
            line(i:i) = cursor%c
            deallocate(cursor)
         end do
         numi = 0
         i = 1
         do
            deltai = verify(line(i:),set)
            if(deltai == 0) exit
            i = i+deltai-1
            numi = numi+1
            deltai = scan(line(i:),set)
            if(deltai == 0) exit
            i = i+deltai-1
         end do
         allocate(array(numi))
         read(line,*) array
      end subroutine gr_1
!end module get_record_mod

!program test
SUBROUTINE test
!   use get_record_mod
   implicit none
   integer iunit
   type(node), pointer :: head
   type(node), pointer :: cursor
!   integer, allocatable :: array(:,:)
   real, allocatable :: array(:,:)
   integer rows
   integer, allocatable :: temp(:)
   integer i
   character(80) fmt
   integer iostat

   allocate(head)
   iunit = 10
   open(iunit,file='test.dat',status='old')
   call get_record(iunit,head)
   rows = 1
   allocate(temp(size(head%array)))
   cursor => head
   do
      read(iunit,*,iostat=iostat) temp
      if(iostat /= 0) exit
      rows = rows+1
      allocate(cursor%next)
      cursor => cursor%next
      allocate(cursor%array(size(temp)))
      cursor%array = temp
   end do
   allocate(array(rows,size(temp)))
   deallocate(temp)
   do i = 1, rows
      cursor => head
      head => head%next
      array(i,:) = cursor%array
      deallocate(cursor)
   end do
! At this point array contains the values from the input file.
   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,i0))'
   write(*,fmt) transpose(array)
!end program test
END SUBROUTINE test

! character(kind=c_char) function f_to_c_string(str,strlen)! result(s)
function f_to_c_string(str,strlen) result(s)
  use iso_c_binding
  integer, intent(in) :: strlen
  character(kind=c_char,len=:), allocatable :: s!(*)
  character(len=strlen), intent(in) :: str
  integer i, nchars
  allocate(character(len=strlen) :: s)
  s = transfer(str(1:nchars), s)
  !f_to_c_string = s
end function f_to_c_string

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
 
subroutine pcstr(s) bind(c,name='pcstr')
  use iso_c_binding
!  use cstr
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars
  write(*,'(a)') c_to_f_string(s)
end subroutine pcstr

subroutine print_cstring_array(n, cstring) bind(C)

  use iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_loc, c_null_char
  implicit none

  integer(kind=c_int),               intent(in) :: n
  type(c_ptr), dimension(n), target, intent(in) :: cstring
  character, pointer                            :: fstring(:)
  integer                                       :: slen, i

  do i = 1, n
    call c_f_pointer(cstring(i), fstring, [4])
    write(*,*) fstring
  end do

end subroutine print_cstring_array

!subroutine print_cstring_array_beta(n, cstring) bind(C)
!
!!  use iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_loc, c_null_char, c_char
!  use iso_c_binding
!  implicit none
!
!  integer(kind=c_int),                 intent(in) :: n
!  type(c_ptr), target,                 intent(in) :: cstring
!  character(kind=c_char), dimension(:,:), pointer :: fptr
!  character(len=3), dimension(n)                  :: fstring
!  integer(kind=c_int) :: i
!
!  call c_f_pointer(c_loc(cstring), fptr, [3, n])
!  do i = 1, n
!     slen = 0
!     do while(fptr(slen+1,i) /= c_null_char)
!        slen = slen + 1
!     end do
!     fstring(i) = transfer(fptr(1:slen,i), fstring(i))
!     write(*,*) trim(fstring(i))
!  end do                                                
!
!end subroutine print_cstring_array_beta

end module get_record_mod

!--
!write(*,*) transfer((/17.392111325966148d0,3.6351694777236872d228, &
!6.0134700169991705d-154/),(/'x'/)); end



PROGRAM AntG
    !use fortranc
    use linked_list_mod ! ADDED BY C.SALGADO ON 2018-04-18
    use get_record_mod
    use antmod
    use, intrinsic :: ISO_C_BINDING
    implicit none
    integer, parameter :: ifu_log = 109
    integer, parameter :: ifu_com = 110
    integer :: ios
    CHARACTER(len=32) :: logfile, basfile
    real, dimension(:,:), allocatable :: S, H
    real, dimension(:,:), allocatable :: PSO
!---------------------------------------------------------------------------
!------------- GET RECORD MOD PART -----------------------------------------
!---------------------------------------------------------------------------
!IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig
   integer :: IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig
   real :: denerrj, Crit
  integer, PARAMETER :: zero=0, two=2, one=1
  logical    :: ANTOn
   integer :: JCycle
  logical :: UHF

   integer iunit
   type(node), pointer :: head
   type(node), pointer :: cursor
!   integer, allocatable :: array(:,:)
   integer :: inputNAtoms,inputNShell,NBasis,NAtoms,inputNSpin,inputNAE,inputNBE,inputNE
   integer, allocatable :: inputJAN(:), inputAOS(:), inputShellT(:), inputShellC(:)
   integer, allocatable :: inputIAN(:)
   real, allocatable :: inputAtmChg(:)
   real, allocatable :: inputAtmCo(:,:)
   Character (len=2), allocatable :: strElem(:)
   real, allocatable :: array(:,:)
   real, allocatable :: dens(:,:,:)
   real, allocatable :: olddens(:,:,:)
   real, allocatable :: singlespindens(:,:)
   real, allocatable :: ham(:,:)
   real, allocatable :: pivham(:,:)
   real :: EE, fockerror, denserror
   real, allocatable :: fock(:,:)
   real, allocatable :: ufock(:,:,:)
   real, allocatable :: oldufock(:,:,:)
   real, allocatable :: overlap(:,:)
   real, allocatable :: Coulomb(:,:)
   real, allocatable :: Exchange(:,:,:)
   real, allocatable :: eig(:)
   real, allocatable :: ERI(:,:,:,:)
   real, allocatable :: ERIsparse(:), bigERIsparse(:)
   real :: pivERIsparse
   integer, allocatable :: ERIsparseindexes(:,:), bigERIsparseindexes(:,:)
   integer :: pivERIsparseindexes(4)
   integer :: isizeERIsparse, ERIsize, isparse
   CHARACTER(LEN=256) :: newline
   integer rows
!   integer, allocatable :: temp(:)
   real, allocatable :: temp(:)
   integer :: i, j, k, l
   integer :: iSpin
   character(80) fmt
   integer iostat
   real, dimension(2) :: Q
   integer :: inputNCycles
   real :: ALPHA, BETA ! TO USE WITH DGEMM 
   real :: Hart = 27.2113834d0
!---------------------------------------------------------------------------
!------------- GETARG PART -------------------------------------------------
!---------------------------------------------------------------------------
!    integer :: i
    integer :: narg!, natoms
!    CHARACTER(len=32) :: arg
    CHARACTER(len=32), dimension(:), allocatable :: arg
    CHARACTER(len=50) :: antname, inputjobname
!    integer(C_INT), VALUE :: inputjobname_len
!    character(inputjobname_len,kind=C_CHAR) :: inputjobname
!    character(kind=C_CHAR,len=:) :: inputcjobname
!    character(kind=C_CHAR,len=:) :: inputcjobname(:)
!    character(kind=C_CHAR) :: inputcjobname(:)
!    character(kind=C_CHAR,len=:), allocatable :: inputcjobname(*)
!    character(kind=C_CHAR,len=:), allocatable :: inputcjobname(:)
!    character(kind=C_CHAR,len=:), allocatable :: inputcjobname(:)

! ORIGINAL WORKED I 2017 IN PORTATIL.
!     character(len=52,kind=C_CHAR) :: inputcjobname ! SIMILAR TO ANT.f90 BUT WITHOUT SPECIFIC LENGTH.
! NEW CODE 2018-12-31 QUITAMOS C_CHAR COMO EN mainOnlyEri.f90
     character(len=52) :: inputcjobname


    logical :: SOCtrue
    print *, "PROGRAM STARTED!!!"
    narg = iargc()
    allocate (arg(narg))
    DO i = 1,narg
      CALL getarg(i, arg(i))
      WRITE (*,*) arg(i)
    END DO
    inputjobname = arg(1)
    print *, "inputjobname = ",inputjobname

!    if(narg>1)then
!        antname = arg(2)
!        open(1,file='name',STATUS='unknown')
!        Write(1,'(A)')antname
!        close(1)
!    end if

    if(narg>1)then
        !NAtoms = arg(3)
        read(arg(2),*)NAtoms
        Write(1,'(A,I4)')"Read NAtoms = ",NAtoms
    end if
!    ifu_com = 110
!---------------------------------------------------------------------------

    !inputjobname="Au9-Au4-Au-Au4-Au9"
    print *, "inputjobname = ", inputjobname

!    open(1,file=inputjobname)
    IRwH=2001
    IRwPA=2002
    IRwPB=2003
    IRwFA=2004
    IRwFB=2005
    IRwS1=2006
    IRwEig=2007
    denerrj = 1.0d-5
    Crit = 1.0d-5
    AntOn = .TRUE.
    UHF = .FALSE.
    JCycle = 0
    ! inputNAtoms,inputNShell,NBasis,inputNAE,inputNBE,inputNE
    ! FOR 10 ATOMS
!    inputNAtoms = 10
!    inputNShell = 20
!    allocate(inputAOS(inputNShell),&
!             inputJAN(inputNShell),&
!             inputShellT(inputNShell),&
!             inputShellC(inputNShell))

    ! FOR 2 ATOMS.
    inputNAtoms = NAtoms

    inputNShell = 2*NAtoms
    allocate(inputAOS(inputNShell),&
             inputJAN(inputNShell),&
             inputShellT(inputNShell),&
             inputShellC(inputNShell))
    if(NAtoms==10)then
      inputAOS = (/ 1, 2, 6, 7,&
                    11, 12, 16, 17,&
                    21, 22, 26, 27,&
                    31, 32, 36, 37,&
                    41, 42, 46, 47 /)
      inputJAN = (/ 1, 1, 2, 2,&
                    3, 3, 4, 4,&
                    5, 5, 6, 6,&
                    7, 7, 8, 8,&
                    9, 9, 10, 10 /)
      inputShellT = (/ 0, 1, 0, 1,&
                       0, 1, 0, 1,&
                       0, 1, 0, 1,&
                       0, 1, 0, 1,&
                       0, 1, 0, 0 /)
      inputShellT = (/ 0, 0, 0, 0,&
                       0, 0, 0, 0,&
                       0, 0, 0, 0,&
                       0, 0, 0, 0,&
                       0, 0, 0, 0 /)
    else
      inputAOS = (/ 1, 2, 6, 7 /)
      inputJAN = (/ 1, 1, 2, 2 /)
      inputShellT = (/ 0, 1, 0, 1 /)
      inputShellT = (/ 0, 0, 0, 0 /)
    end if
    
!---------------------------------------------------------------------------
!------------- GET RECORD MOD PART -----------------------------------------
!---------------------------------------------------------------------------
!------------- READ DENS ---------------------------------------------------
   print *, "-------------------------------------------------------------\n"
   print *, "------------------------- READ DENS -------------------------\n"
   print *, "-------------------------------------------------------------\n"
   allocate(head)
   iunit = 10
!   open(iunit,file='test.dat',status='old')
   open(iunit,file='DENS.' // trim(inputjobname) // '.dat',status='old')
   call get_record(iunit,head)
   rows = 1
   allocate(temp(size(head%array)))
   cursor => head
   do
      read(iunit,*,iostat=iostat) temp
      if(iostat /= 0) exit
      rows = rows+1
      print *, 'rows',rows
      allocate(cursor%next)
      cursor => cursor%next
      allocate(cursor%array(size(temp)))
      cursor%array = temp
   end do
   allocate(array(rows,size(temp)))
   deallocate(temp)
   do i = 1, rows
      cursor => head
      head => head%next
      array(i,:) = cursor%array
      deallocate(cursor)
      print *, array(i,1)
   end do
! At this point array contains the values from the input file.
!   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,i0))'
   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,f12.6))'
   print *, fmt
   write(*,fmt) transpose(array)

   inputNSpin = 1
   IF(UHF) inputNSpin = 2
   allocate(dens(inputNSpin,SIZE(array,1),SIZE(array,2)))
   dens(1,:,:)=transpose(array)
   IF(UHF)dens(2,:,:)=dens(2,:,:)

   deallocate(array)
   deallocate(head)
   print *, "-------------------------------------------------------------\n"
   print *, "----------------------- END READ DENS -----------------------\n"
   print *, "-------------------------------------------------------------\n"
!------------- READ HAM ----------------------------------------------------
   print *, "-------------------------------------------------------------\n"
   print *, "------------------------- READ HAM --------------------------\n"
   print *, "-------------------------------------------------------------\n"
   allocate(head)
   iunit = 10
!   open(iunit,file='test.dat',status='old')
   open(iunit,file='HAM.' // trim(inputjobname) // '.dat',status='old')
   call get_record(iunit,head)
   rows = 1
   allocate(temp(size(head%array)))
   cursor => head
   do
      read(iunit,*,iostat=iostat) temp
      if(iostat /= 0) exit
      rows = rows+1
      allocate(cursor%next)
      cursor => cursor%next
      allocate(cursor%array(size(temp)))
      cursor%array = temp
   end do
   allocate(array(rows,size(temp)))
   deallocate(temp)
   do i = 1, rows
      cursor => head
      head => head%next
      array(i,:) = cursor%array
      deallocate(cursor)
   end do
! At this point array contains the values from the input file.
!   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,i0))'
   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,f12.6))'
   write(*,fmt) transpose(array)

   allocate(ham(SIZE(array,1),SIZE(array,2)))
   ham=transpose(array)

   deallocate(array)
   deallocate(head)
   print *, "-------------------------------------------------------------\n"
   print *, "----------------------- END READ HAM ------------------------\n"
   print *, "-------------------------------------------------------------\n"
!------------- READ FOCK ---------------------------------------------------
   print *, "-------------------------------------------------------------\n"
   print *, "------------------------- READ FOCK -------------------------\n"
   print *, "-------------------------------------------------------------\n"
   allocate(head)
   iunit = 10
!   open(iunit,file='test.dat',status='old')
   open(iunit,file='FOCK.' // trim(inputjobname) // '.dat',status='old')
   call get_record(iunit,head)
   rows = 1
   allocate(temp(size(head%array)))
   cursor => head
   do
      read(iunit,*,iostat=iostat) temp
      if(iostat /= 0) exit
      rows = rows+1
      allocate(cursor%next)
      cursor => cursor%next
      allocate(cursor%array(size(temp)))
      cursor%array = temp
   end do
   allocate(array(rows,size(temp)))
   deallocate(temp)
   do i = 1, rows
      cursor => head
      head => head%next
      array(i,:) = cursor%array
      deallocate(cursor)
   end do
! At this point array contains the values from the input file.
!   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,i0))'
   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,f12.6))'
   write(*,fmt) transpose(array)

   allocate(fock(SIZE(array,1),SIZE(array,2)))
   fock=transpose(array)

   deallocate(array)
   deallocate(head)
   print *, "-------------------------------------------------------------\n"
   print *, "----------------------- END READ FOCK -----------------------\n"
   print *, "-------------------------------------------------------------\n"
!------------- READ OVERLAP ------------------------------------------------
   print *, "-------------------------------------------------------------\n"
   print *, "------------------------ READ OVERLAP -----------------------\n"
   print *, "-------------------------------------------------------------\n"
   allocate(head)
   iunit = 10
!   open(iunit,file='test.dat',status='old')
   open(iunit,file='OV.' // trim(inputjobname) // '.dat',status='old')
   call get_record(iunit,head)
   rows = 1
   allocate(temp(size(head%array)))
   cursor => head
   do
      read(iunit,*,iostat=iostat) temp
      if(iostat /= 0) exit
      rows = rows+1
      allocate(cursor%next)
      cursor => cursor%next
      allocate(cursor%array(size(temp)))
      cursor%array = temp
   end do
   allocate(array(rows,size(temp)))
   deallocate(temp)
   do i = 1, rows
      cursor => head
      head => head%next
      array(i,:) = cursor%array
      deallocate(cursor)
   end do
! At this point array contains the values from the input file.
!   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,i0))'
   write(fmt,'(a,i0,a)') '(',size(array,2),'(1x,f12.6))'
   write(*,fmt) transpose(array)

   allocate(overlap(SIZE(array,1),SIZE(array,2)))
   overlap=transpose(array)

   deallocate(array)
   deallocate(head)
   print *, "-------------------------------------------------------------\n"
   print *, "---------------------- END READ OVERLAP ---------------------\n"
   print *, "-------------------------------------------------------------\n"
!------------------------------------------------------------------------------------------
!--------------------- Z-MATRIX -----------------------------------------------------------
!------------------------------------------------------------------------------------------
   print *, "-------------------------------------------------------------\n"
   print *, "-------------------------- READ ERI -------------------------\n"
   print *, "-------------------------------------------------------------\n"
   open(iunit,file='ERI.' // trim(inputjobname) // '.dat',status='old')
    DO
      READ (iunit,fmt="(A)",iostat=ios) newline
!      print *, "newline = ", newline
!      IF((newline(1:1) == '!') .OR. (newline(1:1) == '*') .OR. &
!        (newline(1:1) == '%') .OR. (newline(1:1) == '#') ) CYCLE
!      IF( (len_trim(newline) == 0)  .OR. (newline(1:1) == '\n') .OR. &
!          (newline(1:1) == '') .OR. (newline(1:1) == ' ') )THEN
!        write (*, *) "blank line AFTER Z-MATRIX"
!        EXIT
!      END IF

!      IF((newline(1:1) == ' ') .OR. (newline(1:1) == '\n') .OR. &
!         (newline(1:1) == '')) EXIT
       IF (ios < 0) THEN
         Write(*,'(A)')"END OF ERI FILE REACHED"
         EXIT
       ENDIF
      READ (newline,fmt=*) pivERIsparseindexes(1), pivERIsparseindexes(2), pivERIsparseindexes(3), pivERIsparseindexes(4), pivERIsparse
      Write(*,'(I4,I4,I4,I4,F18.6,A)') pivERIsparseindexes(1), pivERIsparseindexes(2), pivERIsparseindexes(3), pivERIsparseindexes(4),pivERIsparse

      if(allocated(ERIsparse)) then
        isizeERIsparse = size(ERIsparse)
        allocate(bigERIsparse(isizeERIsparse+1))
        allocate(bigERIsparseindexes(isizeERIsparse+1,4))
        do i=1,isizeERIsparse
          bigERIsparseindexes(i,:) = ERIsparseindexes(i,:)     
          bigERIsparse(i) = ERIsparse(i)
        end do
        bigERIsparseindexes(isizeERIsparse+1,:) = pivERIsparseindexes(:)
        bigERIsparse(isizeERIsparse+1) = pivERIsparse ! HERE IS WHERE I HAVE TO ADD THE NEW COORD TO THE FULL ERI(:).
        deallocate(ERIsparse)
        call move_alloc(bigERIsparse, ERIsparse)
        call move_alloc(bigERIsparseindexes, ERIsparseindexes)
      else
        allocate(ERIsparse(1))
        allocate(ERIsparseindexes(1,4))
        ERIsparseindexes(1,:) = pivERIsparseindexes(:)
        ERIsparse(1) = pivERIsparse ! HERE IS WHERE I HAVE TO ADD THE NEW COORD TO THE FULL ERIsparse(:).
      end if
      !PAUSE
    END DO
   print *, "-------------------------------------------------------------\n"
   print *, "------------------------ END READ ERI -----------------------\n"
   print *, "-------------------------------------------------------------\n"
   print *, "-------------------------------------------------------------\n"
   print *, "-------------------- CONVERT ERI TO TENSOR ------------------\n"
   print *, "-------------------------------------------------------------\n"
!   ERIsize=MAXVAL(MAXVAL(ERIsparseindexes))
   ERIsize=MAXVAL(ERIsparseindexes)
   Write(*,'(A,I4)')"ERIsize =",ERIsize
   NBasis=SIZE(ham,1)
   IF(ERIsize<NBasis)ERISize=NBasis
   IF(ERIsize>NBasis)THEN
     Write(*,'(A)')"ERI size is bigger than NBasis. Truncating it."
   ENDIF
   Allocate(ERI(ERIsize,ERIsize,ERIsize,ERIsize))
   ERI=0.0d0
   do isparse=1,size(ERIsparse)
     IF((ERIsparseindexes(isparse,1)>NBasis).or.(ERIsparseindexes(isparse,2)>NBasis).or.(ERIsparseindexes(isparse,3)>NBasis).or.(ERIsparseindexes(isparse,4)>NBasis))CYCLE ! TRUNCATING ERI.
     ERI(ERIsparseindexes(isparse,1),ERIsparseindexes(isparse,2),ERIsparseindexes(isparse,3),ERIsparseindexes(isparse,4))=ERIsparse(isparse)
   end do

   Allocate(Coulomb(ERIsize,ERIsize),&
            Exchange(inputNSpin,ERIsize,ERIsize))
   Coulomb=0.0d0
   Exchange=0.0d0
   do iSpin=1,inputNSpin
     do i=1,ERIsize
       do j=1,ERIsize
         do k=1,ERIsize
           do l=1,ERIsize
             Coulomb(i,j)=Coulomb(i,j)+eri(i,j,k,l)*dens(iSpin,k,l)
             Exchange(iSpin,i,k)=Exchange(iSpin,i,k)-0.5*eri(i,k,j,l)*dens(iSpin,j,l)
           end do
         end do
       end do
     end do
   end do
   do i=1,NBasis
     write(*, '(A,I2,A,10F8.4)')"Coulomb(",i,",:) = ",( Coulomb(i,j) ,j=1,NBasis)
   end do
   do iSpin=1,inputNSpin
     do i=1,NBasis
       write(*, '(A,I1,A,I2,A,10F8.4)')"Exchange(",iSpin,",",i,",:) = ",( Exchange(iSpin,i,j) ,j=1,NBasis)
     end do
   end do
   print *, "-------------------------------------------------------------\n"
   print *, "------------------ END CONVERT ERI TO TENSOR ------------------\n"
   print *, "-------------------------------------------------------------\n"
!------------------------------------------------------------------------------------------
!------------ FINISHED Z-MATRIX -----------------------------------------------------------
!------------------------------------------------------------------------------------------
   print *, "-------------------------------------------------------------\n"
   print *, "-------------------------- READ GEOM ------------------------\n"
   print *, "-------------------------------------------------------------\n"
!c
!c  This program reads n points from a data file and stores them in 
!c  3 arrays x, y, z.
!c
!      integer nmax, u
!      parameter (nmax=1000, u=20)
!      real x(nmax), y(nmax), z(nmax)

!c  Open the data file
      open (iunit, FILE=trim(inputjobname) // '.xyz', STATUS='OLD')

!c  Read the number of points
      read(iunit,*) inputNAtoms
      Write(*,'(A,I3)')"InputNAtoms read from xyz header: ", inputNAtoms
      read(iunit,*)
!      if (n.GT.nmax) then
!         write(*,*) 'Error: n = ', n, 'is larger than nmax =', nmax
!         goto 9999
!      endif

      allocate(inputIAN(inputNAtoms))
      allocate (strElem(inputNAtoms))
      allocate(inputAtmCo(3,inputNAtoms))
      allocate(inputAtmChg(inputNAtoms))

      write(*,*) 'inputNAtoms',inputNAtoms
!c  Loop over the data points
      do i=1,inputNAtoms
         !read(iunit,100) x(i), y(i), z(i)
         read(iunit,100) strElem(i), inputAtmCo(1,i), inputAtmCo(2,i), inputAtmCo(3,i)
         write(*,100) strElem(i), inputAtmCo(1,i), inputAtmCo(2,i), inputAtmCo(3,i)
         inputIAN(i) = 3
         inputAtmChg(i) = 3.0d0
   10 enddo
!  100 format (3(F10.4))
!  100 format ((A2),3(F14.6))
  100 format ((A2),3(F12.6))

!c  Close the file
      close (iunit)

      inputNE = SUM(inputIAN(:))
      IF (MOD(inputNE,2)==0) THEN
            inputNAE = inputNE/2
            inputNBE = inputNE/2
            PRINT*, 'THE NUMBER IS AN EVEN NUMBER'
            PRINT*, 'input(NAE,NBE)=(',inputNAE,',',inputNBE,')'
      ELSE
            inputNAE = FLOOR(REAL(inputNE/2))+1
            inputNBE = FLOOR(REAL(inputNE/2))
            PRINT*, 'THE NUMBER IS AN ODD NUMBER'
            PRINT*, 'input(NAE,NBE)=(',inputNAE,',',inputNBE,')'
      END IF

!c  Now we should process the data somehow...
!c  (missing part)

! 9999 stop
   print *, "-------------------------------------------------------------\n"
   print *, "------------------------ END READ GEOM ----------------------\n"
   print *, "-------------------------------------------------------------\n"
!---------------------------------------------------------------------------
!------------- END GET RECORD MOD PART -------------------------------------
!---------------------------------------------------------------------------






open(unit=2009, file='CONVER.' // trim(inputjobname) // '.dat', ACTION="write", STATUS="replace")
!    CALL ANT (UHF,JCycle,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn,NBasis)
    ! THE LINE BELOW WORKS.
    !CALL ANT (UHF,JCycle,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn,SIZE(ham,1),dens,ham,fock,overlap)
    !CALL ANT (UHF,JCycle,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn,NBasis,pivDens,pivHam,pivFock,pivOver)
    NBasis=SIZE(ham,1)
    allocate(ufock(inputNSpin,SIZE(fock,1),SIZE(fock,2)))
    allocate(oldufock(inputNSpin,SIZE(fock,1),SIZE(fock,2)))
    allocate(pivham(SIZE(fock,1),SIZE(fock,2)))

!    IF(UHF)THEN
!      !ufock = MATMUL(dens(1,:,:)+dens(2,:,:),(MATMUL(OVERLAP,fock)))
!    ELSE
!      !ufock = MATMUL(2*dens(1,:,:),(MATMUL(OVERLAP,fock)))
!    END IF
    do iSpin=1,inputNSpin
      ufock(iSpin,:,:)=Coulomb(:,:)+Exchange(iSpin,:,:)
    end do

    allocate(olddens(inputNSpin,NBasis,NBasis))
    allocate(singlespindens(NBasis,NBasis))
    olddens=dens
    oldufock=ufock


!    allocate(inputcjobname(LEN(inputjobname)))
!    allocate(character(len=len_trim(inputjobname)) :: inputcjobname)

    !inputcjobname=f_to_c_string(inputjobname,LEN(inputjobname))
    !inputcjobname=transfer(inputcjobname,inputjobname, size = len_trim(inputjobname))
    inputcjobname=transfer(inputjobname,inputcjobname)
    !string = transfer(text, ' ', size = len_trim(text))
    !CALL pcstr(inputcjobname)
    !CALL print_chararray_c(inputcjobname,LEN(inputcjobname))
    !call print_cstring_array(len(inputjobname), inputcjobname)
    Write(*,*)"STRING PASSED AS INPUTCJOBNAME TO ANT IS : ",inputcjobname
    Write(*,'(A)')"PRESS ANY KEY TO CONTINUE"
    PAUSE

    JCycle=0
    inputNCycles=999999

    DO JCycle=0,inputNCycles,1
!     DO WHILE (JCycle .LE. inputNCycles)
      olddens=dens
      !oldufock = 0.9d0*oldufock+0.1d0*ufock
      oldufock = 0.9d0*oldufock+0.1d0*ufock  
      !PRINT *, "UHF,JCycle,IRwH,IRwPA,IRwPB,IRwFA,IRwFB",UHF,JCycle,IRwH,IRwPA,IRwPB,IRwFA,IRwFB
      !PRINT *, "IRwS1,IRwEig,denerrj,Crit,ANTOn",IRwS1,IRwEig,denerrj,Crit,ANTOn
      !PRINT *, "inputNAtoms,inputNShell,NBasis",inputNAtoms,inputNShell,NBasis
      !PRINT *, "inputNAE,inputNBE,inputNE",inputNAE,inputNBE,inputNE
!      CALL ANT(UHF,JCycle,inputNCycles,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn,&
!                 inputNAtoms,inputNShell,NBasis,inputNSpin,inputJAN,inputAOS,inputShellT,inputShellC,inputNAE,&
!                 inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,dens,ham,oldufock,overlap)

 


! COMMENTED ON 2018-04-18 TO REMOVE IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig.
!      CALL ANT(UHF,JCycle,inputNCycles,inputcjobname,LEN(inputcjobname),dens,ham,oldufock,Coulomb,Exchange,overlap,&
!                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
!                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
!                 IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn)
      CALL ANT(UHF,JCycle,inputNCycles,inputcjobname,LEN(inputcjobname),dens,ham,oldufock,Coulomb,Exchange,overlap,&
                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
                 denerrj,Crit,ANTOn) 




     IF(JCycle > 0)THEN
        do i=1,NBasis
          write(*, '(A,I2,A,10F8.4)')"HAM(",i,",:) = ",( ham(i,j) ,j=1,NBasis)
        end do
        Coulomb=0.0d0
        Exchange=0.0d0
        do iSpin=1,inputNSpin
          do i=1,ERIsize
            do j=1,ERIsize
              do k=1,ERIsize
                do l=1,ERIsize
                  Coulomb(i,j)=Coulomb(i,j)+eri(i,j,k,l)*dens(iSpin,k,l)
                  Exchange(iSpin,i,k)=Exchange(iSpin,i,k)-0.5*eri(i,k,j,l)*dens(iSpin,j,l)
                end do
              end do
            end do
          end do
        end do
        do iSpin=1,inputNSpin
          ufock(iSpin,:,:)=Coulomb(:,:)+Exchange(iSpin,:,:)
        end do
        do i=1,NBasis
          write(*, '(A,I2,A,10F8.4)')"Coulomb(",i,",:) = ",( Coulomb(i,j) ,j=1,NBasis)
        end do
        do iSpin=1,inputNSpin
          do i=1,NBasis
            write(*, '(A,I1,A,I2,A,10F8.4)')"Exchange(",iSpin,",",i,",:) = ",( Exchange(iSpin,i,j) ,j=1,NBasis)
          end do
        end do
        !print *, dens(1,1:10,1:10)
        !print *, "= dens(1,1:3,1:3)"
        do iSpin=1,inputNSpin
          do i=1,NBasis
            write(*, '(A,I1,A,I2,A,10F8.4)')"D(",iSpin,",",i,",:) = ",( dens(iSpin,i,j) ,j=1,NBasis)
          end do
        end do
        print *, "= dens(1,1:10,1:10)"
        print *, pivham(1:3,1:3)
        do i=1,NBasis
          write(*, '(A,I1,A,I2,A,10F8.4)')"UFOCK(",iSpin,",",i,",:) = ",( ufock(iSpin,i,j) ,j=1,NBasis)
        end do
        print *, "= ufock(1:10,1:10)"
        EE=0.0d0
        do iSpin=1,inputNSpin
          DO i=1,NBasis
            EE = EE + pivham(i,i) + ufock(iSpin,i,i)
          END DO
        end do
        fockerror = MAXVAL(MAXVAL(MAXVAL(ufock-oldufock,1),2))
        denserror = MAXVAL(MAXVAL(MAXVAL(dens-olddens,1),2))
        if((.NOT.(ISNAN(fockerror)) .AND. (fockerror < 5.0d-3)).AND.(.NOT.(ISNAN(denserror)) .AND. (denserror < 5.0d-4)))EXIT
        write(*,*)"*******************************************************************"
        write(*,*)"********** PRINTING ENERGY AT CYCLE: ",JCycle," *****"
        write(*,'(3(A,F16.8))')" ENERGY =",EE,"; FERR =",fockerror,"; DERR =",denserror
        write(*,*)"*******************************************************************"
        write(2009, '(A,I10,A,F20.10,A,F20.10,A,F20.10)')"JCycle:",JCycle,"; Energy:",EE,"; FERR:",fockerror,"; DERR:",denserror
      END IF

    END DO
    write(*,*)"EXITED SCF LOOP IN AntG.f90"


        write(*,*)"*******************************************************************"
        write(*,*)"********************** POSTPROCESSING CYCLE ***********************"
        write(*,*)"*******************************************************************"
      JCycle=inputNCycles

!      CALL ANT(UHF,JCycle,inputNCycles,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn,&
!                 inputNAtoms,inputNShell,NBasis,inputNSpin,inputJAN,inputAOS,inputShellT,inputShellC,inputNAE,&
!                 inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,dens,ham,oldufock,overlap)
!       CALL ANT(UHF,JCycle,inputNCycles,inputcjobname,LEN(inputcjobname),dens,ham,oldufock,Coulomb,Exchange,overlap,&
!                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
!                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
!                 IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn)

 


! COMMENTED ON 2018-04-18 TO REMOVE IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig.
!      CALL ANT(UHF,JCycle,inputNCycles,inputcjobname,LEN(inputcjobname),dens,ham,oldufock,Coulomb,Exchange,overlap,&
!                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
!                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
!                 IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn)
      CALL ANT(UHF,JCycle,inputNCycles,inputcjobname,LEN(inputcjobname),dens,ham,oldufock,Coulomb,Exchange,overlap,&
                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
                 denerrj,Crit,ANTOn) 




        do i=1,NBasis
          write(*, '(A,I2,A,10F8.4)')"HAM(",i,",:) = ",( ham(i,j) ,j=1,NBasis)
        end do
        do i=1,NBasis
          write(*, '(A,I2,A,10F8.4)')"HAM(",i,",:) = ",( ham(i,j) ,j=1,NBasis)
        end do
        Coulomb=0.0d0
        Exchange=0.0d0
        do iSpin=1,inputNSpin
          do i=1,ERIsize
            do j=1,ERIsize
              do k=1,ERIsize
                do l=1,ERIsize
                  Coulomb(i,j)=Coulomb(i,j)+eri(i,j,k,l)*dens(iSpin,k,l)
                  Exchange(iSpin,i,k)=Exchange(iSpin,i,k)-0.5*eri(i,k,j,l)*dens(iSpin,j,l)
                end do
              end do
            end do
          end do
        end do
        do i=1,NBasis
          write(*, '(A,I2,A,10F8.4)')"Coulomb(",i,",:) = ",( Coulomb(i,j) ,j=1,NBasis)
        end do
        do iSpin=1,inputNSpin
          do i=1,NBasis
            write(*, '(A,I1,A,I2,A,10F8.4)')"Exchange(",iSpin,",",i,",:) = ",( Exchange(iSpin,i,j) ,j=1,NBasis)
          end do
        end do
        !print *, dens(1,1:10,1:10)
        !print *, "= dens(1,1:3,1:3)"
        do iSpin=1,inputNSpin
          do i=1,NBasis
            write(*, '(A,I1,A,I2,A,10F8.4)')"D(",iSpin,",",i,",:) = ",( dens(iSpin,i,j) ,j=1,NBasis)
          end do
        end do
        print *, "= dens(1,1:10,1:10)"
        print *, pivham(1:3,1:3)
        do i=1,NBasis
          write(*, '(A,I2,A,10F8.4)')"UFOCK(",i,",:) = ",( ufock(iSpin,i,j) ,j=1,NBasis)
        end do
        print *, "= ufock(1:10,1:10)"
        EE=0.0d0
        do iSpin=1,inputNSpin
          DO i=1,NBasis
            EE = EE + pivham(i,i) + ufock(iSpin,i,i)
          END DO
        end do
        fockerror = MAXVAL(MAXVAL(MAXVAL(ufock-oldufock,1),2))
        denserror = MAXVAL(MAXVAL(MAXVAL(dens-olddens,1),2))
        !if((.NOT.(ISNAN(fockerror)) .AND. (fockerror < 5.0d-3)).AND.(.NOT.(ISNAN(denserror)) .AND. (denserror < 5.0d-4)))EXIT
        write(*,*)"*******************************************************************"
        write(*,*)"********** PRINTING ENERGY AT CYCLE: ",JCycle," *****"
        write(*,'(3(A,F16.8))')" ENERGY =",EE,"; FERR =",fockerror,"; DERR =",denserror
        write(*,*)"*******************************************************************"
        write(2009, '(A,I10,A,F20.10,A,F20.10,A,F20.10)')"JCycle:",JCycle,"; Energy:",EE,"; FERR:",fockerror,"; DERR:",denserror

        write(*,*)"*******************************************************************"
        write(*,*)"******************** END POSTPROCESSING CYCLE *********************"
        write(*,*)"*******************************************************************"



    Q=0.0d0
    do iSpin=1,inputNSpin
      do i=1,NBasis
        do j=1,NBasis
          Q(iSpin)=Q(iSpin)+dens(iSpin,i,j)*overlap(j,i)
        end do
      end do
      write(*, '(A,I1,A,F8.4)')"Q(",iSpin,") = ",Q(iSpin)
    end do
    write(*, '(A,F8.4)')"QTot = ",SUM(Q(:))
  open(unit=2008, file='DENS.' // trim(inputjobname) // '_final.dat', ACTION="write", STATUS="replace")
  !open(iunit,file='DENS.' // trim(inputjobname) // '.dat',status='old')
  do iSpin=1,inputNSpin
    write(*,'(A,I1)')"SPIN",iSpin
    write(2008,'(A,I1)')"SPIN",iSpin
    do i=1,NBasis
      !write(*,'(A,I1,A,I2,A,I2,A,F14.6)')"D(",iSpin,",",i,",",j,") = ",D(iSpin,i,j)
      write(2008, '(10F8.4)')( dens(iSpin,i,j) ,j=1,NBasis)
      write(*, '(10F8.4)')( dens(iSpin,i,j) ,j=1,NBasis)
      !write(*, 10(1600F14.7))( dens(iSpin,i,j) ,j=1,NBasis)
    end do
  end do
  print *, "*******************************************************************"
  print *, "*********************** PROGRAM FINISHED!!! ***********************"
  print *, "*******************************************************************"
  close(2001)
  close(2002)
  close(2003)
  close(2004)
  close(2005)
  close(2006)
  close(2007)
  close(2008)
  close(2009)
END PROGRAM
