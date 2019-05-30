!#include "stdout.h"
!**********************************************************
!*********************  ANT.G-2.4.1  **********************
!**********************************************************
!*                                                        *
!*  Copyright (c) by                                      *
!*                                                        *
!*  Juan Jose Palacios (1)                                *
!*  David Jacob (2)                                       *
!*  Maria Soriano (1)                                     *
!*  Angel J. Perez-Jimenez (3)                            *
!*  Emilio SanFabian (3)                                  *
!*  Jose Antonio Antonio Verges (4)                       *
!*  Enrique Louis (5)                                     *
!*                                                        *
!* (1) Departamento de Fisica de la Materia Condensada    *
!*     Universidad Autonoma de Madrid                     *      
!*     28049 Madrid (SPAIN)                               *
!* (2) Theory Department                                  *
!*     Max-Planck-Institute for Microstructure Physics    *
!*     Halle, 06120 (GERMANY)                             *
!* (3) Departamento de Quimica Fisica                     *
!*     Universidad de Alicante                            *
!*     03690 Alicante (SPAIN)                             *
!* (4) Insto. de Ciencias de Materiales de Madrid (ICMM)  *
!*     Consejo Superior de Investigacion y Ciencia (CSIC) *
!*     28049 Madrid (SPAIN)                               *
!* (5) Departamento de Fisica Aplicada                    *
!*     Universidad de Alicante                            *      
!*     03690 Alicante (SPAIN)                             *
!*                                                        *
!**********************************************************
  MODULE AntMod
!*********************************************************!
!  Module for analysis of system provided in gaussian     !
!  input file .com                                        !
!**********************************************************
!  use cluster ! ADDED TO MATCH ARECONNECTED, WHICH IS A FUNCTION IN CLUSTER
  USE preproc !, ONLY: MaxAtm, MaxSh
  USE parameters, ONLY: NEmbed,NAtomEl
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: antc
  PUBLIC :: ANT, readmat

  integer, PARAMETER :: nvec = 12

!  real  :: vpb1(3,nvec), vpb2(3,nvec) 
!  integer :: ndir(MaxAtm), nvbet(MaxAtm,nvec) 
!  integer :: ifrpl(MaxAtm), norbmt
!  integer :: ANLead1, ANLead2, nlead1, nmol, nlead2, nneig1, nneig2, NAtomData
!  integer :: nn1, nn2
!
!  integer, DIMENSION(MaxAtm) :: ANMol, NAO, LAO, HAO, NSh, AN
!  CHARACTER, Dimension(100,MaxAtm) :: ShAtm  
!  real,  DIMENSION(MaxAtm) :: XL1, YL1, ZL1, XM, YM, ZM, XL2, YL2, ZL2

CONTAINS
! THE DEFINITION BELOW WORKS WITH FORTRAN TYPES BUT NOT ISO_C_BINDING ONES.
!  SUBROUTINE ANT(UHF,JCycle,inputNCycles,inputjobname,inputjobname_len,D,pivHam,pivFock,pivCoulomb,pivExchange,pivOver,&
!                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
!                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
!                 IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn) bind(c)
!  SUBROUTINE antc(UHF,JCycle,inputNCycles,Cinputjobname,inputjobname_len,D,pivHam,pivFock,pivCoulomb,pivExchange,pivOver,&
!                 NBasis,inputNSpin,inputNAtoms,inputNShell,CinputJAN,CinputAOS,CinputShellT,CinputShellC,&
!                 inputNAE,inputNBE,inputNE,CinputIAN,inputAtmCo,inputAtmChg,&
!                 denerrj,Crit,ANTOn) bind(c)
                 SUBROUTINE antc(UHF,JCycle,inputNCycles,Cinputjobname,Cinputjobname_len,&
                 NBasis,inputNSpin,inputNAtoms,inputNShell,&
                 DA,DB,pivHam,pivFockA,pivFockB,pivCoulomb,pivExchangeA,pivExchangeB,pivOver,&
                 outHWFockA,outHWFockB,outGibbsYA,outGibbsYB,CinputIAN,CinputJAN,CinputAOS,CinputShellT,CinputShellC,&
                 inputNAE,inputNBE,inputNE,inputAtmCo,inputAtmChg,&
                 denerrj,Crit,ANTOn) bind(c)
!**********************************************************************************************************************
!* Interface subroutine with Gaussian                                                                                 *
!**********************************************************************************************************************
  USE Parameters, ONLY: SL, SwOffSPL, alpha, Read_Parameters, Write_Parameters, NSpinLock, npulay
  USE Parameters, ONLY: ChargeAcc,ChargeA,FermiAcc,FermiA,PAcc,PA,FullAcc,RedTransmB,RedTransmE,ElType,LDOS_Beg,LDOS_End
  USE Parameters, ONLY: Mulliken, Hamilton, PFix, DFTU, FMixing, IntEnergy, DiagFock, SpinMu, CompFock
  USE constants, ONLY: Hart
  USE preproc
  USE device, ONLY: InitDevice, DevFockMat, DevDensMat, ReadDensMat, LeadsOn, DevShift, SwitchOnLeads, &
       SwitchOnEvaluation, SwitchOnSecant, SwitchOffSecant, SwitchOnSpinLock, SwitchOffSpinLock, &
       SwitchOnChargeCntr, SwitchOffChargeCntr, transport, CleanUpDevice, SetDevDensMat, SetDevFockMat, ReadFockMat, &
       DevHWFockMat, SetDevHWFockMat, DevDGibbsYMat, SetDevDGibbsYMat
!  USE g09Common, ONLY: GetNShell, GetAtm4Sh, Get1stAO4Sh, GetNBasis, GetAN, GetAtmChg, GetAtmCo, GetNAtoms
  USE g09Common
  use ANTCommon
  use util
!  use, intrinsic :: ISO_C_BINDING
  use iso_c_binding
  IMPLICIT NONE
  ! dummy arguments
!  logical, INTENT(in)      :: UHF
  logical(c_bool), INTENT(in)      :: UHF
!  integer, INTENT(inout) :: JCycle, inputNCycles
  integer(c_int), INTENT(inout) :: JCycle, inputNCycles
!  integer(C_INT), intent(IN), VALUE :: Cinputjobname_len
  integer(C_INT), intent(IN) :: Cinputjobname_len
!  character(inputjobname_len,kind=C_CHAR), intent(IN) :: inputjobname
!  character(inputjobname_len,kind=C_CHAR) :: Cinputjobname
  character, dimension(*), intent(in) :: Cinputjobname
  character(Cinputjobname_len) :: inputjobname
!  integer, INTENT(in)    :: inputNSpin,inputNAtoms,inputNShell,NBasis,inputNAE,inputNBE,inputNE
  integer(c_int), INTENT(in)    :: NBasis,inputNSpin,inputNAtoms,inputNShell
!  integer(c_int) :: NBasis, inputNSpin, inputNAtoms, inputNShell
!  real,DIMENSION(inputNSpin,NBasis,NBasis),INTENT(inout) :: D !, pivDens
!  real(c_double),DIMENSION(inputNSpin,NBasis,NBasis),INTENT(inout) :: D !, pivDens
  real(c_double),DIMENSION(NBasis,NBasis),INTENT(inout) :: DA, DB !, pivDens
!  real(c_double),DIMENSION(inputNSpin,NBasis,NBasis) :: D !, pivDens
  real,DIMENSION(inputNSpin,NBasis,NBasis) :: D
!  real,DIMENSION(NBasis,NBasis),INTENT(in) :: pivHam, pivOver, pivCoulomb
  real(c_double),DIMENSION(NBasis,NBasis),INTENT(in) :: pivHam, pivOver, pivCoulomb
!  real(c_double),DIMENSION(NBasis,NBasis) :: pivHam, pivOver, pivCoulomb
!  real,DIMENSION(inputNSpin,NBasis,NBasis),INTENT(in) :: pivFock, pivExchange
!  real(c_double),DIMENSION(inputNSpin,NBasis,NBasis),INTENT(in) :: pivFock, pivExchange
  real(c_double),DIMENSION(NBasis,NBasis),INTENT(in) :: pivFockA, pivFockB, pivExchangeA, pivExchangeB
  real(c_double),DIMENSION(inputNSpin,NBasis,NBasis) :: pivFock, pivExchange
  real(c_double),DIMENSION(NBasis,NBasis),INTENT(out) :: outHWFockA, outHWFockB
  real(c_double),DIMENSION(NBasis,NBasis),INTENT(out) :: outGibbsYA, outGibbsYB
  real(c_double),DIMENSION(inputNSpin,NBasis,NBasis) :: outHWFock
!  real(c_double),DIMENSION(inputNSpin,NBasis,NBasis) :: pivFock, pivExchange
!  integer, DIMENSION(inputNAtoms), INTENT(in):: inputIAN
  integer(c_long) :: CinputIAN(inputNAtoms)
  integer, DIMENSION(inputNAtoms):: inputIAN
  integer(c_long) :: CinputJAN(inputNShell), CinputAOS(inputNShell), CinputShellT(inputNShell), CinputShellC(inputNShell)
  integer, DIMENSION(inputNShell):: inputJAN, inputAOS, inputShellT, inputShellC
  integer(c_long), INTENT(in)    :: inputNAE,inputNBE,inputNE
!  integer(c_int)   :: inputNAE,inputNBE,inputNE
!  real, DIMENSION(inputNAtoms), INTENT(in) :: inputAtmChg
  real(c_double), DIMENSION(inputNAtoms), INTENT(in) :: inputAtmChg
!  real(c_double), DIMENSION(inputNAtoms) :: inputAtmChg
!  real, DIMENSION(3,inputNAtoms), INTENT(in) :: inputAtmCo
  real(c_double), DIMENSION(3,inputNAtoms), INTENT(in) :: inputAtmCo
!  real(c_double), DIMENSION(3,inputNAtoms) :: inputAtmCo
!  real,INTENT(in)        :: denerrj, Crit
!  real(c_double),INTENT(in)        :: denerrj, Crit
  real(c_double)        :: denerrj, Crit
!  logical,INTENT(inout)    :: ANTOn
!  logical(c_bool),INTENT(inout)    :: ANTOn
  logical(c_bool)    :: ANTOn
!  integer, INTENT(in)    :: IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig!,NBasis
!  integer(c_int), INTENT(in)    :: IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig!,NBasis
!  integer :: IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig!,NBasis
  integer :: inputNBasis ! LATER OVERWRITTEN BY NBasis, ADDED BY C.SALGADO ON 2017-03-23.
!  integer, DIMENSION(inputNShell), INTENT(in):: inputJAN, inputAOS, inputShellT, inputShellC
!  integer(c_int), DIMENSION(inputNShell), INTENT(in):: inputJAN, inputAOS, inputShellT, inputShellC
!  integer, DIMENSION(NBasis), INTENT(in):: 

  logical :: ADDP
  logical, SAVE :: FInit
  ! local variables
  real    :: val, density, fock, exspin, norma
  integer :: ic, ndim, i, j, info, ii, jj, ispin, acount, AllocErr, ios ,k, iatom, jatom, ntot, ia
  integer, DIMENSION(MaxAtm) :: NAO
  CHARACTER(len=80), SAVE :: densitymatrix, fockmatrix
  integer, SAVE :: isw, NCycLeadsOn, NSpin
  integer, SAVE :: ntt, len
  integer, PARAMETER :: zero=0, two=2, one=1
  real, DIMENSION(:),ALLOCATABLE   :: TS
  real, DIMENSION(:,:),ALLOCATABLE :: S
  real, DIMENSION(:,:,:),ALLOCATABLE :: PtimesS

  ! Matrices in lower trinagular form (alpha and beta) 
  ! for communication with gaussian
!  real, DIMENSION(:,:),ALLOCATABLE,SAVE ::  MT 
  !
  ! Fock and density matrices in regular form
  !
  real, DIMENSION(:,:,:),ALLOCATABLE,SAVE :: F!, D ! COMMENTED MATRIX D ON 2017-03-22.
  real, DIMENSION(:,:,:,:),ALLOCATABLE,SAVE ::  DD
  real, DIMENSION(:,:,:),ALLOCATABLE,SAVE ::  DDD
  real, DIMENSION(npulay+1,npulay+1) ::  b
  real, DIMENSION(npulay+1) ::  c
  integer, DIMENSION(npulay+1) ::  ipiv

  real :: fmix
  real, DIMENSION(2) :: TraceCharge, TraceChargeNoOv
  real, DIMENSION(:,:),ALLOCATABLE ::  SingleSpinDens
  character(80) strfmt

  Write(*,*)"**********************************************"
  Write(*,*)"**************** ENTERED ANTC ****************"
  Write(*,*)"**********************************************"
  ! VERY IMPORTANT ADDED BY C.SALGADO ON 2017-03-23 TO OVERWRITE inputNBasis by NBasis.
  Write(*,'(A)')"UPDATE INPUTNBASIS AND NSPIN TO AVOID ERRORS WITH THE SIZE OF ARRAYS"
  inputNBasis = NBasis
  NSpin = inputNSpin
  IF(JCYCLE==0)THEN
  Write(*,'(A)')"**************************************************************************"
  Write(*,'(A)')"****** FILL DENSITYD(:,:,:) WITH DA(:,:) AND DB(:,:) FOR JCYCLE = 0 ******"
  Write(*,'(A)')"**************************************************************************"
  do i=1,NBasis
    do j=1,NBasis
      D(1,i,j)=DA(i,j)
      pivFock(1,i,j)=pivFockA(i,j)
      !  pivCoulomb(1,:,:)=pivCoulombA(:,:)
      pivExchange(1,i,j)=pivExchangeA(i,j)
    end do
  end do
  if(NSpin.eq.2)then
    do i=1,NBasis
      do j=1,NBasis
        D(2,i,j)=DB(i,j)
        pivFock(2,i,j)=pivFockB(i,j)
!        pivCoulomb(2,i,j)=pivCoulombB(i,j)
        pivExchange(2,i,j)=pivExchangeB(i,j)
      end do
    end do
  end if
  Write(*,*)"JCycle = ",JCycle
  Write(*,*)"inputNCycles = ",inputNCycles
  Write(*,'(A)')"**************************************************************************"
  Write(*,'(A)')"**** END FILL DENSITYD(:,:,:) WITH DA(:,:) AND DB(:,:) FOR JCYCLE = 0 ****"
  Write(*,'(A)')"**************************************************************************"
  END IF

  IF(JCycle>0) THEN
  Write(*,'(A)')"*****************************************************************************"
  Write(*,'(A)')"************************ UPDATE EVERY CYCLE EXCEPT 0 ************************"
  Write(*,'(A)')"*****************************************************************************"
  Write(*,'(A)')"**************************************************************************"
  Write(*,'(A)')"****** FILL DENSITYD(:,:,:) WITH DA(:,:) AND DB(:,:) FOR JCYCLE > 0 ******"
  Write(*,'(A)')"**************************************************************************"
  Write(*,'(A)')"FILL ALPHA DENSITY MATRIX D(1,:,:)"
  do i=1,NBasis
    do j=1,NBasis
      D(1,i,j)=DA(i,j)
      pivFock(1,i,j)=pivFockA(i,j)
      !  pivCoulomb(1,:,:)=pivCoulombA(:,:)
      pivExchange(1,i,j)=pivExchangeA(i,j)
    end do
    Write(*, '(10F8.4)')( D(1,i,j) ,j=1,NBasis)
  end do
  if(NSpin.eq.2)then
    Write(*,'(A)')"FILL BETA DENSITY MATRIX D(2,:,:)"
    do i=1,NBasis
      do j=1,NBasis
        D(2,i,j)=DB(i,j)
        pivFock(2,i,j)=pivFockB(i,j)
!        pivCoulomb(2,i,j)=pivCoulombB(i,j)
        pivExchange(2,i,j)=pivExchangeB(i,j)
      end do
      Write(*, '(10F8.4)')( D(2,i,j) ,j=1,NBasis)
    end do
  end if
  Write(*,*)"JCycle = ",JCycle
  Write(*,*)"inputNCycles = ",inputNCycles
  Write(*,'(A)')"**************************************************************************"
  Write(*,'(A)')"**** END FILL DENSITYD(:,:,:) WITH DA(:,:) AND DB(:,:) FOR JCYCLE > 0 ****"
  Write(*,'(A)')"**************************************************************************"

  Write(*,'(A)')"*****************************************************************************"
  Write(*,'(A)')"********************** END UPDATE EVERY CYCLE EXCEPT 0 **********************"
  Write(*,'(A)')"*****************************************************************************"
  end if
!  do iSpin=1,NSpin
!    do i=1,NBasis
!      do j=1,NBasis
!        if(isnan(D(iSpin,i,j)))D(iSpin,i,j)=0.0d0
!        !write(*,'(A,I1,A,I2,A,I2,A,F14.6)')"D(",iSpin,",",i,",",j,") = ",D(iSpin,i,j)
!      end do
!    end do
!  end do

  Write(*,'(A)')"*******************************************************************"
  Write(*,'(A)')"******************** FOR USE WITH HARTREE-FOCK ********************"
  Write(*,'(A)')"*******************************************************************"
  do ia=1,inputNAtoms
    Write(*,'(A,I2,A,I2)')"CinputIAN(",ia,") = ",CinputIAN(ia)
  end do
  inputIAN = CinputIAN
  inputNBasis = NBasis
  NSpin = inputNSpin
  inputAOS=CinputAOS
  inputJAN=CinputJAN
  inputShellT=CinputShellT
  inputShellC=CinputShellC
  Write(*, '(A,10I2)')"inputIAN = ",( inputIAN(j) ,j=1,inputNAtoms)
  Write(*, '(A,10I2)')"inputAOS = ",( inputAOS(j) ,j=1,inputNShell)
  Write(*, '(A,10I2)')"inputJAN = ",( inputJAN(j) ,j=1,inputNShell)
  Write(*, '(A,10I2)')"inputShellT = ",( inputShellT(j) ,j=1,inputNShell)
  Write(*, '(A,10I2)')"inputShellC = ",( inputShellC(j) ,j=1,inputNShell)

  Write(*,*)"Cinputjobname = ",Cinputjobname(1:Cinputjobname_len)
  Write(*,*)"Cinputjobname_len = ",Cinputjobname_len
  Write(*,*)"NBasis = ",NBasis
  Write(*,*)"inputNSpin = ",inputNSpin
  Write(*,*)"inputNAToms = ",inputNAtoms
  Write(*,*)"inputNShell = ",inputNShell

  Write(*,'(A)')"************************************************"
  Write(*,'(A)')"****** COORDINATES RECEIVED IN inputAtmCo ******"
  Write(*,'(A)')"************************************************"
  do ia=1,inputNAtoms
    write(*, '(A,I2,A,10F8.4)')"inputAtmCo(1:3,",ia,") = ",( inputAtmCo(j,ia) ,j=1,3)
  end do
  Write(*,'(A)')"************************************************"
  Write(*,'(A)')"**** END COORDINATES RECEIVED IN inputAtmCo ****"
  Write(*,'(A)')"************************************************"
  Write(*,'(A,I4)')"JCycle = ",JCycle
  Write(*,'(A)')"ALPHA DENSITY MATRIX DA(:,:)"
  do i=1,NBasis
    Write(*, '(10F16.10)')( DA(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"BETA DENSITY MATRIX DB(:,:)"
  do i=1,NBasis
    Write(*, '(10F16.10)')( DB(i,j) ,j=1,NBasis)
  end do
  do iSpin=1,NSpin
    if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY MATRIX D(1,:,:)"
    if(iSpin==2)Write(*,'(A)')"BETA DENSITY MATRIX D(2,:,:)"
    do i=1,NBasis
      Write(*, '(10F16.10)')( D(iSpin,i,j) ,j=1,NBasis)
    end do
  end do
  Write(*,'(A)')"HAMILTONIAN pivHam(:,:)"
  do i=1,NBasis
    Write(*, '(18F16.10)')( pivHam(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"FOCK MATRIX ALPHA pivFockA(:,:)"
  do i=1,NBasis
    Write(*, '(10F16.10)')( pivFockA(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"COULOMB MATRIX pivCoulomb(:,:)"
  do i=1,NBasis
    Write(*, '(10F16.10)')( pivCoulomb(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"EXCHANGE MATRIX ALPHA pivExchangeA(:,:)"
  do i=1,NBasis
    Write(*, '(10F16.10)')( pivExchangeA(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"OVERLAP MATRIX pivOver(:,:)"
  do i=1,NBasis
    Write(*, '(18F16.10)')( pivOver(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A,I2,A,I2,A,I2)')"inputNAE=",inputNAE,"; inputNBE=",inputNBE,"; inputNE=",inputNE
  Write(*,'(A)')"*******************************************************************"
  Write(*,'(A)')"****************** END FOR USE WITH HARTREE-FOCK ******************"
  Write(*,'(A)')"*******************************************************************"
  !Write(*,'(A)')"PRESS ANY KEY..."
  !Pause
  !*******************************************************************
  ! Before first cycle: Initialize module device, allocate memory etc
  !*******************************************************************
  IF(JCycle.EQ.0) THEN
     NSpin = 1
     IF(UHF) NSpin = 2
     ntt = (NBasis*(NBasis+1))/2
     
     ANTOn = .FALSE.
          
     PRINT *
     PRINT *, " ****************************************************************** "
     PRINT *, " ***                                                            *** "
     PRINT *, " ***                      A l i c a n t e                       *** "
     PRINT *, " ***                      N a n o                               *** "
     PRINT *, " ***                      T r a n s p o r t                     *** "
     PRINT *, " ***                      G 0 9                                 *** "
     PRINT *, " ***                                                            *** "
     PRINT *, " ****************************************************************** "
     PRINT *, " ***                     Version: 2.3.4                         *** "
     PRINT *, " ****************************************************************** "
     PRINT *, " *  Copyright (c) by                                              * "
     PRINT *, " *                                                                * "
     PRINT *, " *  Juan Jose Palacios (1)                                        * "
     PRINT *, " *  David Jacob (2)                                               * "
     PRINT *, " *  Maria Soriano (1,5)                                           * "
     PRINT *, " *  Angel J. Perez-Jimenez (3)                                    * "
     PRINT *, " *  Emilio SanFabian (3)                                          * "
     PRINT *, " *  Jose Antonio Antonio Verges (4)                               * "
     PRINT *, " *  Enrique Louis (5)                                             * "
     PRINT *, " *                                                                * "
     PRINT *, " * (1) Departamento de Fisica de la Materia Condensada            * "
     PRINT *, " *     Universidad Autonoma de Madrid                             * "
     PRINT *, " *     28049 Madrid (SPAIN)                                       * "
     PRINT *, " * (2) Theory Department                                          * "
     PRINT *, " *     Max-Planck-Institute for Microstructure Physics            * "
     PRINT *, " *     Halle, 06120 (GERMANY)                                     * "
     PRINT *, " * (3) Departamento de Quimica Fisica                             * " 
     PRINT *, " *     Universidad de Alicante                                    * " 
     PRINT *, " *     03690 Alicante (SPAIN)                                     * "
     PRINT *, " * (4) Insto. de Ciencias de Materiales de Madrid (ICMM)          * "
     PRINT *, " *     Consejo Superior de Investigacion y Ciencia (CSIC)         * "
     PRINT *, " *     28049 Madrid (SPAIN)                                       * "
     PRINT *, " * (5) Departamento de Fisica Aplicada                            * "
     PRINT *, " *     Universidad de Alicante                                    * "
     PRINT *, " *     03690 Alicante (SPAIN)                                     * "
     PRINT *, " *                                                                * "
     PRINT *, " ****************************************************************** "
     PRINT *
     PRINT *, " ****************************************************************** "
     PRINT *, " *                                                                * "
     PRINT *, " *              Initializing ANT.G09 modules                      * "
     PRINT *, " *                                                                * "
     PRINT *, " ****************************************************************** "
     PRINT *
     
    !*********************************************************************
    !**** Before first cycle: Set all variables formerly in g09Common ****
    !*********************************************************************
    PRINT *, "Before first cycle: Set all variables formerly in g09Common"
    inputJAN=CinputJAN
    inputAOS=CinputAOS
    inputShellT=CinputShellT
    inputShellC=CinputShellC
    inputIAN=CinputIAN
    inputJAN=CinputJAN
    CALL SetNAtoms( inputNAtoms )
    CALL SetNShell( inputNShell )
    CALL SetAtm4Sh(inputNShell, inputJAN)
    CALL Set1stAO4Sh( inputNShell, inputAOs )
    CALL SetShellT( inputNShell, inputShellT )
    CALL SetShellC( inputNShell, inputShellC )
    CALL SetNAE( inputNAE )
    CALL SetNBE( inputNBE )
    CALL SetNE( inputNE )
    CALL SetNBasis( inputNBasis )
    CALL SetAN( inputNAtoms, inputIAN )
    CALL SetAtmChg( inputNAtoms, inputAtmChg )
    CALL SetAtmCo( inputNAtoms, inputAtmCo )

    do ia=1,inputNAtoms
      write(*, '(A,I2,A,10F16.10)')"inputAtmCo(1:3,",ia,") = ",( inputAtmCo(j,ia) ,j=1,3)
    end do
    PRINT *, "Before first cycle: End set all variables formerly in g09Common"

    !*********************************************************************
    ! Before first cycle: End setting all variables formerly in g09Common
    !*********************************************************************
     ! Reading parameter file

! THE COMMENTED CODE BELOW WAS USED TO READ FILE CALLED "NAME"
!     OPEN(UNIT=ifu_nam,FILE="name",IOSTAT=ios,STATUS='OLD')
!     IF( ios == 0 ) THEN
!        READ(ifu_nam,*) jobname
!        OPEN(UNIT=ifu_ini,FILE=trim(jobname)//'.ini',IOSTAT=ios,STATUS='OLD')
!        IF( ios == 0 ) THEN
!          CALL read_parameters( ifu_log, ifu_ini, ios, jobname )
!        ELSE
!          WRITE(ifu_log,*) " "
!          WRITE(ifu_log,*) "*** No parameter file found. Using default values."
!          WRITE(ifu_log,*) " "
!        END IF
!        CLOSE(ifu_ini)
!     ELSE
!        WRITE(ifu_log,*) " "
!        WRITE(ifu_log,*) "*** No name file found. Using default values."
!        WRITE(ifu_log,*) " "
!     END IF
!     CLOSE(UNIT=ifu_nam)
!     jobname = C_F_STRING(inputjobname)
     jobname = c_to_f_string(Cinputjobname)
     OPEN(UNIT=ifu_ini,FILE=trim(jobname)//'.ini',IOSTAT=ios,STATUS='OLD')
     IF( ios == 0 ) THEN
       CALL read_parameters( ifu_log, ifu_ini, ios, jobname )
     ELSE
       WRITE(ifu_log,*) " "
       WRITE(ifu_log,*) "*** No parameter file found. Using default values."
       WRITE(ifu_log,*) " "
     END IF
     CLOSE(ifu_ini)

     ant1dname = jobname
     call locase(ant1dname)

     CALL write_parameters( ifu_log )

     !Opening writting files
     OPEN(ifu_xyz,file=trim(jobname)//'_BL.xyz',status='unknown')
     IF (ElType(1) /= 'GHOST' .and. ElType(2) /= 'GHOST') OPEN(ifu_tra,file='T.'//trim(jobname)//'.dat',status='unknown')
     IF (RedTransmB < RedTransmE) OPEN(ifu_red,file='t.'//trim(jobname)//'.dat',status='unknown')
     IF (Hamilton) OPEN(ifu_ham,file='V.'//trim(jobname)//'.dat',status='unknown')
     IF (DiagFock) OPEN(ifu_diagham,file='DF.'//trim(jobname)//'.dat',status='unknown')
     IF (Mulliken) OPEN(ifu_mul,file='Q.'//trim(jobname)//'.dat',status='unknown')
     IF (LDOS_Beg <= LDOS_End) OPEN(ifu_dos,file='DOS.'//trim(jobname)//'.dat',status='unknown')

     ! Creating density matrix file name
     densitymatrix='P.' // trim(jobname) // '.dat'
     ! Creating Fock matrix file name
     fockmatrix='F.' // trim(jobname) // '.dat'
     !
     ! Allocate memory for dynamic arrays
     !
!     ALLOCATE( MT( NSpin, ntt ), STAT=AllocErr ) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!     IF( AllocErr /= 0 ) THEN ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!        PRINT *, "ESF/Error: could not allocate memory for MT(:)" ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!        STOP ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!     END IF ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
     PRINT *, "(NSpin,NBasis)=(",NSpin,",",NBasis,")"
     ALLOCATE( F(NSpin,NBasis,NBasis), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for F(:,:,:)"
        STOP
     END IF
!     DO iSpin=1,NSpin
!       DO i=1,NBasis
!          DO j=1,NBasis
!            F(iSpin,i,j) = pivOver(i,j) ! ADDED AND COMMENTED BY C.SALGADO ON 2017-03-25. FOCK MATRIX IS SET VIA F.name.dat.
!          END DO
!       END DO
!     END DO
      !FOCK MATRIX HD(iSpin,i,j) will be initialized in InitDevice by doing HD = F.

     if(.not. FMixing )then
!        ALLOCATE( D(NSpin,NBasis,NBasis), DD(npulay+1,NSpin,NBasis,NBasis),DDD(npulay,NSpin,NBasis*NBasis), STAT=AllocErr )
        ALLOCATE( DD(npulay+1,NSpin,NBasis,NBasis),DDD(npulay,NSpin,NBasis*NBasis), STAT=AllocErr ) ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
        IF( AllocErr /= 0 ) THEN
!           PRINT *, "ESF/Error: could not allocate memory for D(:,:,:), DD(:,:,:,:),DDD(:,:,:)"
           PRINT *, "ESF/Error: could not allocate memory for DD(:,:,:,:),DDD(:,:,:)" ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
           STOP
        END IF
        DD=0.0
        DDD=0.0
     end if


     !
     ! Obtain overlap matrix in lower triangular form
     ! and transform to regular form matrix
     !
     ALLOCATE( TS(ntt), S(NBasis,NBasis), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for S(:,:), TS(:)"
        STOP
     END IF

!     CALL FileIO(two,-IRwS1,ntt,TS,zero)
     acount = 1
     DO i=1,NBasis
        DO j=1,i
!           S(i,j) = TS(acount)
!           S(j,i) = TS(acount)
           S(i,j) = pivOver(i,j)
           S(j,i) = pivOver(j,i)
           acount = acount+1
        END DO
     END DO

     ! 
     ! Initialize module device 
     !
     CALL InitDevice( NBasis, UHF, S )

     CLOSE(ifu_xyz)

     if(FMixing)OPEN(ifu_fm,file=fockmatrix,iostat=ios,status='old')
     if(.not.FMixing)OPEN(ifu_dm,file=densitymatrix,iostat=ios,status='old')
     IF (ios /= 0) THEN
        FInit=.false.
     ELSE
        FInit=.true.
        if(FMixing)CLOSE(ifu_fm)
        if(.not.FMixing)CLOSE(ifu_dm)
     END IF

     IF( FInit )THEN
        !
        ! Initialization from file F.dat or P.dat
        ! File F.name.dat just never exists, so this is not the usual procedure to get the Hamiltonian.
        !
        if(FMixing)then
           CALL ReadFockMat(fockmatrix)
           DO i=1,NBasis
              DO j=1,NBasis
                 !F(1,i,j) = DevFockMat(1,i,j)
                 !if(UHF) F(2,i,j) = DevFockMat(2,i,j)
                 F(1,i,j) = pivFock(1,i,j)
                 if(UHF) F(2,i,j) = pivFock(2,i,j)
              end DO
           end DO
        else
           Write(*,'(A)')"DISPLAY D(iSpin,i,j) BEFORE ReadDensMat"
           DO ispin=1,NSpin
             if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY D(1,:,:)"
             if(iSpin==2)Write(*,'(A)')"BETA DENSITY D(2,:,:)"
             DO i=1,NBasis
!              DO j=1,NBasis
!                D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!                D(iSpin,j,i)=D(iSpin,i,j)
!              END DO
               Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
             END DO
           END DO
           Write(*,'(A)')"END D(iSpin,i,j) BEFORE ReadDensMat"
           CALL ReadDensMat(densitymatrix)        
           !
           ! Transform to triangular form
           !
           acount = 1
           DO i=1,NBasis
              DO j=1,i
                 D(1,i,j)=DevDensMat(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-25 BECAUSE D IS INOUT OF ANT.
                 !MT(1,acount)=D(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
                 IF(UHF) D(2,i,j) = DevDensMat(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-25 BECAUSE D IS INOUT OF ANT.
                 !IF(UHF) MT(2,acount) = D(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
                 acount = acount +1
              ENDDO
           ENDDO
           Write(*,'(A)')"DISPLAY D(iSpin,i,j) AFTER DevDensMat"
           DO ispin=1,NSpin
             if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY D(1,:,:)"
             if(iSpin==2)Write(*,'(A)')"BETA DENSITY D(2,:,:)"
             DO i=1,NBasis
!              DO j=1,NBasis
!                D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!                D(iSpin,j,i)=D(iSpin,i,j)
!              END DO
               Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
             END DO
           END DO
           Write(*,'(A)')"END D(iSpin,i,j) AFTER DevDensMat"
           !
           ! PASS DENSITY MATRIX TO GAUSSIAN VIA RWF
           !
!           CALL FileIO(one,-IRwPA,ntt,MT(1,:),zero)
!           IF (UHF) CALL FileIO(one,-IRwPB,ntt,MT(2,:),zero)
        end if
     ENDIF

     DEALLOCATE( S, TS, STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not DEAllocate memory for S(:,:), TS(:)"
        STOP
     END IF

     PRINT *, "RETURNING FROM ANT(...)"
     RETURN
  END IF

  IF( .NOT. LeadsOn() )THEN
     !
     ! Connect Leads in first cycle
     ! when reinitialized
     !
     IF( FInit )THEN
        ANTOn = .true.
        CALL SwitchOnLeads
        NCycLeadsOn = 0
        !
        ! Or if SL criterion is met:
        !
     ELSEIF( denerrj .LE. SL ) THEN
        ANTOn = .true.
        Call SwitchOnLeads
        NCycLeadsOn = 0
        if(.not.FMixing)then
           !
           ! Obtain density matrix from gaussian only in first step 
           ! with leads switched on and if not initialized from file
           !
!           CALL FileIO(two,-IRwPA,ntt,MT(1,:),zero)
!           IF(UHF) CALL FileIO(two,-IRwPB,ntt,MT(2,:),zero)

           Write(*,'(A)')"DISPLAY D(iSpin,i,j) BEFORE SetDevDensMat"
           DO ispin=1,NSpin
             if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY D(1,:,:)"
             if(iSpin==2)Write(*,'(A)')"BETA DENSITY D(2,:,:)"
             DO i=1,NBasis
!              DO j=1,NBasis
!                D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!                D(iSpin,j,i)=D(iSpin,i,j)
!              END DO
               Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
             END DO
           END DO
           Write(*,'(A)')"END D(iSpin,i,j) BEFORE SetDevDensMat"
           
           acount = 1
           DO i=1,NBasis
              DO j=1,i
!                 D(1,i,j)=MT(1,acount)
!                 D(1,j,i)=D(1,i,j)
                 !D(1,i,j)=pivDens(i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                 !D(1,j,i)=D(1,i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                 call SetDevDensMat( 1, i, j, D(1,i,j) )
                 call SetDevDensMat( 1, j, i, D(1,j,i) )
                 !call SetDevFockMat( 1, i, j, F(1,i,j) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                 !call SetDevFockMat( 1, i, j, F(1,j,i) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                 IF (UHF) THEN
!                    D(2,i,j)=MT(2,acount)
!                    D(2,j,i)=D(2,i,j)
                    !D(2,i,j)=pivDens(i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                    !D(2,j,i)=D(2,i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                    call SetDevDensMat( 2, i, j, D(2,i,j) )
                    call SetDevDensMat( 2, j, i, D(2,j,i) )
                    !call SetDevFockMat( 2, i, j, F(2,i,j) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                    !call SetDevFockMat( 2, i, j, F(2,j,i) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                 ENDIF
                 acount = acount+1
              ENDDO
           ENDDO

           Write(*,'(A)')"DISPLAY D(iSpin,i,j) AFTER SetDevDensMat"
           DO ispin=1,NSpin
             if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY D(1,:,:)"
             if(iSpin==2)Write(*,'(A)')"BETA DENSITY D(2,:,:)"
             DO i=1,NBasis
!              DO j=1,NBasis
!                D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!                D(iSpin,j,i)=D(iSpin,i,j)
!              END DO
               Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
             END DO
           END DO
           Write(*,'(A)')"END D(iSpin,i,j) AFTER SetDevDensMat"

        end if
     END IF
  END IF

  !IF( jcycle.EQ.1000 )THEN
  IF( jcycle.EQ.inputNCycles )THEN
     !NCycLeadsOn = 999
     NCycLeadsOn = inputNCycles-1
     CALL SwitchOnEvaluation()
     CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. .not. FullAcc .and. jcycle > 2) THEN
     IF( denerrj .LT. Crit*10) THEN
        ChargeAcc = ChargeA
        FermiAcc = FermiA
        PAcc = PA
     ELSE IF( denerrj .LT. Crit*100) THEN
        ChargeAcc = ChargeA*10.0
        FermiAcc = FermiA*10.0
        PAcc = PA*10.0
     ELSE IF( denerrj .LT. Crit*1000) THEN
        ChargeAcc = ChargeA*100.0
        FermiAcc = FermiA*100.0
        PAcc = PA*100.0
     ELSE
        ChargeAcc = ChargeA*1000.0
        FermiAcc = FermiA*1000.0
        PAcc = PA*1000.0
     END IF
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. FullAcc) THEN
     ChargeAcc = ChargeA
     FermiAcc = FermiA
     PAcc = PA
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. .not. FullAcc .and. jcycle <= 2) THEN
     ChargeAcc = ChargeA*10.0
     FermiAcc = FermiA*10.0
     PAcc = PA*10.0
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE
     !
     ! If leads are not switched on or if not in cycle 1000 return to Gaussian
     ! without calculation of density matrix with leads connected
     !
     RETURN
  END IF

  
  !*********************************************
  ! COMPUTE DENSITY MATRIX WITH LEADS CONNECTED
  !*********************************************

  WRITE(ifu_log,*) '-------------------------------------------------------------------------'
  WRITE(ifu_log,*) 'Computing the density matrix at cycle:', JCycle
  WRITE(ifu_log,*) '-------------------------------------------------------------------------'

  ! Obtain Fock matrix from Gaussian RWF
!  CALL FileIO(two,-IRwFA,ntt,MT(1,:),zero)
!  IF(UHF) CALL FileIO(two,-IRwFB,ntt,MT(2,:),zero)

  fmix=1.0d0
  ! For Fock matrix mixing 
  if(FMixing.and.NCycLeadsOn>0) fmix=alpha
  !if(NCycLeadsOn>0) fmix=alpha
  if(FMixing.and.FINIT.and.NCycLeadsOn.eq.0) fmix=0.0d0
  ! Transform from lower trinagular form to regular form
  !print *, "NCycLeadsOn=", NCycLeadsOn
  !print *, "fmix =", fmix
  acount = 1
  Write(*,*)"*********************************************************************************"
  Write(*,*)"***************************** UPDATING FOCK MATRIX ******************************"
  Write(*,*)"*********************************************************************************"
  DO i=1,NBasis
     DO j=1,i
        ! Mixing with old Fock matrix if fmix<1
!        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*MT(1,acount)
!        F(1,j,i)=F(1,i,j)
!        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*(pivHam(i,j)+pivFock(i,j)) ! COMMENTED BY C.SALGADO TO USE F(iSpin,:,:) AS J(Coulomb)+K(Exchange), removing pivHam.
!        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*(pivFock(i,j))
!        F(1,i,j)=Hart*pivFock(i,j)
!        F(1,i,j)=Hart*(pivHam(i,j)+pivFock(i,j))
        F(1,i,j)=(pivHam(i,j)+pivFock(1,i,j))
!        F(1,i,j)=pivHam(i,j)
        F(1,j,i)=F(1,i,j)
        IF (UHF) THEN
!           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*MT(2,acount)
!           F(2,j,i)=F(2,i,j)
!           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*(pivHam(i,j)+pivFock(i,j)) ! COMMENTED BY C.SALGADO TO USE F(iSpin,:,:) AS J(Coulomb)+K(Exchange), removing pivHam.
!           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*(pivFock(i,j))
!           F(2,i,j)=Hart*pivFock(i,j)
!           F(2,i,j)=Hart*(pivHam(i,j)+pivFock(i,j))
           F(2,i,j)=(pivHam(i,j)+pivFock(2,i,j))
!           F(2,i,j)=pivHam(i,j)
           F(2,j,i)=F(2,i,j)
        ENDIF
        acount = acount +1
     ENDDO
  ENDDO

!  write(strfmt,'(a,i0,a)') '(A,I1,A,',size(F,2),'(1x,F12.6))'
!  write(ifu_log,'(A)')strfmt
!  do iSpin=1,NSpin
!    do i=1,NBasis
!      write(*, '(10F12.6)')( F(iSpin,i,j) ,j=1,NBasis)
!    end do
!    !write(ifu_log, strfmt)"F(",iSpin,")",F(iSpin,:,:)
!  end do

  Write(*,*)"*********************************************************************************"
  Write(*,*)"*************************** END UPDATING FOCK MATRIX ****************************"
  Write(*,*)"*********************************************************************************"
  
  if(FMixing)then
     ntot=0
     DO i=1,GetNShell()-1
        IF (GetAtm4Sh(i).NE.GetAtm4Sh(i+1)) THEN
           NAO(GetAtm4Sh(i))=Get1stAO4Sh(i+1)-(ntot+1)
           ntot=ntot+NAO(GetAtm4Sh(i))
        ENDIF
     ENDDO
     NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot
     ! Write density matrix to file F.' // trim(jobname) // '.dat
!     if (NCycLeadsOn  < 1000) then
     if (NCycLeadsOn  < inputNCycles) then
        OPEN(ifu_fm,file=fockmatrix,status='unknown')
        WRITE(ifu_fm,*)DevShift()
        DO ispin=1,NSpin
           i=0
           do iAtom=1,GetNAtoms()
              do ii=1,NAO(iAtom)
                 i=i+1
                 j=0
                 do jAtom=1,GetNAtoms() 
                    do jj=1,NAO(jAtom)
                       j=j+1
                       if (i >= j) WRITE(ifu_fm,'(i2,2i6,e18.10,i5,i6,i5,i6)')ispin,i,j,F(ispin,i,j),iAtom,ii,jAtom,jj
                    end do
                 end do
              end do
           end do
        end do
        close(ifu_fm)
     end if
  end if
  
  ! Count cycles with leads connected
  NCycLeadsOn = NCycLeadsOn + 1
  !print*,NCycLeadsOn
  
  ! Turn on charge control every 5 steps in the first cycles
  IF(MOD(NCycLeadsOn-1,10) == 0 .and. NCycLeadsOn <= 21) CALL SwitchOnChargeCntr
  IF(MOD(NCycLeadsOn-1,20) == 0 .and. NCycLeadsOn > 21) CALL SwitchOnChargeCntr
  
  Write(*,'(A)')"DISPLAY D(iSpin,i,j) BEFORE TRANSPORT"
  DO ispin=1,NSpin
    if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY D(1,:,:)"
    if(iSpin==2)Write(*,'(A)')"BETA DENSITY D(2,:,:)"
    DO i=1,NBasis
!      DO j=1,NBasis
!        D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!        D(iSpin,j,i)=D(iSpin,i,j)
!      END DO
      Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
    END DO
  END DO
  Write(*,'(A)')"END D(iSpin,i,j) BEFORE TRANSPORT"
  ! Call subroutine that solves transport problem
  CALL Transport(F,ADDP)
  ! THE CALL TO DEVDENSMAT BELOW WAS ADDED BY C.SALGADO ON 2017-03-24.
!  Write(*,'(A)')"OVERWRITE D(iSpin,i,j) AFTER TRANSPORT"
  Write(*,'(A)')"DISPLAY D(iSpin,i,j) AFTER TRANSPORT"
  DO ispin=1,NSpin
    if(iSpin==1)Write(*,'(A)')"ALPHA DENSITY D(1,:,:)"
    if(iSpin==2)Write(*,'(A)')"BETA DENSITY D(2,:,:)"
    DO i=1,NBasis
!      DO j=1,NBasis
!        D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!        D(iSpin,j,i)=D(iSpin,i,j)
!      END DO
      Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
    END DO
  END DO
  Write(*,'(A)')"END D(iSpin,i,j) AFTER TRANSPORT"

  if( jcycle.EQ.inputNCycles )THEN
    if(CompFock)then
      Write(*,'(A)')"outHWFockA/B(i,j) AFTER TRANSPORT IN EVALUATION->COMPFOCK"
      Write(*,'(A)')"ALPHA HWFOCK outHWFockA(:,:)"
      do i=1,NBasis
        do j=1,NBasis
          outHWFockA(i,j)=DevHWFockMat(1,i,j)
        end do
        Write(*, '(10F8.4)')( outHWFockA(i,j) ,j=1,NBasis)
      end do
      if(NSpin==2)then
        Write(*,'(A)')"BETA HWFOCK outHWFockB(:,:)"
        do i=1,NBasis
          do j=1,NBasis
            outHWFockB(i,j)=DevHWFockMat(2,i,j)
          end do
          Write(*, '(10F8.4)')( outHWFockB(i,j) ,j=1,NBasis)
        end do
      end if
      Write(*,'(A)')"END outHWFockA/B(i,j) AFTER TRANSPORT IN EVALUATION->COMPFOCK"

      Write(*,'(A)')"outGibbsYA/B(i,j) AFTER TRANSPORT IN EVALUATION->COMPFOCK"
      Write(*,'(A)')"ALPHA GibbsY outGibbsYA(:,:)"
      do i=1,NBasis
        do j=1,NBasis
          outGibbsYA(i,j)=DevDGibbsYMat(1,i,j)
        end do
        Write(*, '(10F8.4)')( outGibbsYA(i,j) ,j=1,NBasis)
      end do
      if(NSpin==2)then
        Write(*,'(A)')"BETA GibbsY outGibbsYB(:,:)"
        do i=1,NBasis
          do j=1,NBasis
            outGibbsYB(i,j)=DevDGibbsYMat(2,i,j)
          end do
          Write(*, '(10F8.4)')( outGibbsYB(i,j) ,j=1,NBasis)
        end do
      end if
      Write(*,'(A)')"END outGibbsYA/B(i,j) AFTER TRANSPORT IN EVALUATION->COMPFOCK"
    end if
  end if

  CALL SwitchOffChargeCntr
  
  IF( SL <= 0.0d0 ) alpha = 1.0d0
 
  !
  ! Pulay accelaration for density matrix mixing (FMixing == false)
  !
  if(.not.FMixing)then 
  
     DD(1,:,:,:)=D(:,:,:)

     DO k=npulay,1,-1
     !DO k=1,npulay ! PUT HERE BY C.SALGADO WITHOUT IDEA.
        DD(k+1,:,:,:)=DD(k,:,:,:)
     END DO
  
     if (ADDP) then
        DO ispin=1,NSpin
           DO i=1,NBasis
              DO j=1,NBasis
                 DD(1,ispin,i,j)=(1.0d0-alpha)*DD(2,ispin,i,j)+alpha*DevDensMat(ispin,i,j)           
              END DO
           END DO
        END DO
     else
        DD(1,:,:,:)=0.5*DD(3,:,:,:)+0.5*DD(2,:,:,:)
     end if
  
!     write(*,*)"DD(kpulay,ispin,ibas,jbas)"
     DO k=1,npulay
        DO ispin=1,NSpin
           ic=0
           DO i=1,NBasis
              DO j=1,NBasis
                 ic=ic+1
                 DDD(k,ispin,ic)=DD(k,ispin,i,j)-DD(k+1,ispin,i,j)   
!                 write(*,'(A,I2,A,I2,A,I2,A,I2,A,F8.4)')"DD(",k,",",ispin,",",i,",",j,")=",DD(k,ispin,i,j)
              END DO
           END DO
        END DO
     END DO
     ! NEXT NESTED LOOP ONLY TO PRINT -DD(kpulay+1,ispin,ibas,jbas)
!     write(*,*)"-DD(kpulay+1,ispin,ibas,jbas)"
!     DO k=1,npulay
!        DO ispin=1,NSpin
!           DO i=1,NBasis
!              DO j=1,NBasis
!                 write(*,'(A,I2,A,I2,A,I2,A,I2,A,F8.4)')"-DD(",k+1,",",ispin,",",i,",",j,")=",DD(k+1,ispin,i,j)
!              END DO
!           END DO
!        END DO
!     END DO
     ! NEXT NESTED LOOP ONLY TO PRINT DDD(kpulay,ispin,ic)
!     write(*,*)"DDD(kpulay,ispin,ic)"
!     DO k=1,npulay
!        DO ispin=1,NSpin
!           DO ic=1,NBasis*NBasis
!             write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"DDD(",k,",",ispin,",",ic,")=",DDD(k,ispin,ic)
!           END DO
!        END DO
!     END DO

     ! Pulay acceleration procedure
     ndim=npulay+1
     if (mod(NCycLeadsOn,ndim)==npulay) then
        print *, "Pulay kick ............."
     
        D=0.0
        DO ispin=1,NSpin
           b=0.0
           c=0.0
           c(ndim)=-1.0
           !Matrix to obtain the coefficients
!           write(*,*)"ispin;b(k,j) MATRIX TO OBTAIN THE COEFFICIENTS FOR PULAY BEFORE DGETRF"
           DO k=1,ndim
              DO j=1,ndim
                 if (k.lt.ndim.and.j.lt.ndim) then
                    DO i=1,NBasis*NBasis
                       b(k,j)=b(k,j)+DDD(k,ispin,i)*DDD(j,ispin,i)
                    END DO
                 else if (k==ndim.and.j.lt.ndim) then
                    b(k,j)=-1.0
                 else if (j==ndim.and.k.lt.ndim) then
                    b(k,j)=-1.0
                 end if
!                 write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,b(",k,",",j,")=",b(k,j)
              END DO
           END DO
           ! EXTRACTED ZGETRF FROM NUMERIC FOR INTEL MKL.
           !call zgetrf(NDIM,NDIM,Sigma_aux,NDIM,ipiv,info) ! ORIGINAL FOR PGI IN BetheLattice.f90.
           !call zgetrf(NDim,NDim,Sigma_aux,NDim+1,ipiv,info) ! ADDING +1 TO PARAMETER 4 WORKS!!!
           Call DGETRF( ndim, ndim, b, ndim, IPIV, INFO) ! ORIGINAL CODE WORKS WITH GAUSSIAN. FAILS WITH INFO=2 WITH MKL.
           !Call DGETRF( ndim, ndim, b, ndim+1, IPIV, INFO)  ! ADDED +1 TO PARAMETER 4 BECAUSE IT WORKS. CHECK MAIL:  call zgemm
           !Call DGETRF( ndim, ndim, b, ndim, IPIV, INFO ) 
!           INFO=0 

           IF( INFO /= 0 ) THEN
             WRITE(ifu_log,*)'Pulay/Problems using DGETRF in ANT.f90'
             WRITE(ifu_log,*)'INFO=',INFO
!             STOP
           END IF


!           write(*,*)"ispin;b(k,j) MATRIX TO OBTAIN THE COEFFICIENTS FOR PULAY AFTER DGETRS"
!           DO k=1,ndim
!              DO j=1,ndim
!                 write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,b(",k,",",j,")=",b(k,j)
!              END DO
!           END DO
!           write(*,*)"ispin;IPIV(k) pivot indices from DGETRF to DGETRS"
!           DO k=1,ndim
!                 write(*,'(A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,IPIV(",k,")=",c(k)
!           END DO
           ! EXTRACTED ZGETRI FROM NUMERIC FOR INTEL MKL.
           !call zgetri(NDIM,Sigma_aux,NDIM,ipiv,work,4*NDIM,info) ! ORIGINAL FOR PGI IN BetheLattice.f90.
           !call zgetri(NDim,Sigma_aux,NDim+1,ipiv,work,4*NDim+1,info) ! ADDING +1 TO PARAMETERS 3,6 WORKS!!!
           !DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
           Call DGETRS( 'N', ndim, 1, b, ndim, IPIV, c, ndim, INFO) ! ORIGINAL CODE WORKS WITH GAUSSIAN.
           !Call DGETRS( 'N', ndim, 1, transpose(b), ndim, IPIV, c, ndim, INFO) ! ORIGINAL CODE WORKS WITH GAUSSIAN.
           !Call DGETRS( 'N', ndim, 1, b, ndim+1, IPIV, c, 4*ndim+1, INFO) ! ADDED +1 TO PARAMETERS 5,8 BECAUSE IT WORKS. CHECK MAIL:  call zgem


           IF( INFO /= 0 ) THEN
             WRITE(ifu_log,*)'Pulay/Problems using DGETRS in ANT.f90'
             WRITE(ifu_log,*)'INFO=',INFO
!             STOP
           END IF


!           write(*,*)"ispin;b(k,j) MATRIX TO OBTAIN THE COEFFICIENTS FOR PULAY AFTER DGETRS"
!           DO k=1,ndim
!              DO j=1,ndim
!                 write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,b(",k,",",j,")=",b(k,j)
!              END DO
!           END DO
!           write(*,*)"ispin;c(k) VECTOR OF THE COEFFICIENTS FOR PULAY AFTER DGETRS"
!           DO k=1,ndim
!                 write(*,'(A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,c(",k,")=",c(k)
!           END DO
!           write(*,*)"c(k)*DD(kpulay,ispin,ibas,jbas)"
           DO k=1,npulay
              DO i=1,NBasis
                 DO j=1,NBasis
!                    write(*,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,F8.4)')"c(",k,")*DD(",k,",",ispin,",",i,",",j,")=",c(k)*DD(k,ispin,i,j)
                    D(ispin,i,j)=D(ispin,i,j)+c(k)*DD(k,ispin,i,j)
                 END DO
              END DO
           END DO
           DO i=1,NBasis
             DO j=i,NBasis
               D(ispin,j,i)=D(ispin,i,j)
             END DO
           END DO

!           DO i=1,NBasis
!             Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
!           END DO

        END DO
     ELSE
       write(*,'(A)')"SUMMING CHARGE AFTER PULAY ACCELERATION"
       TraceCharge = (/ 0.0d0, 0.0d0 /)
       TraceChargeNoOv = (/ 0.0d0, 0.0d0 /)
        write(*,'(A,F8.4)')" ZERO TraceCharge(alpha+beta)=",SUM(TraceCharge(:))
       Allocate(PtimesS(NSpin,NBasis,NBasis),&
                SingleSpinDens(NBasis,NBasis),&
                S(NBasis,NBasis))
       S = pivOver

       write(*,'(A,I1,A,I2)')"NSpin=",NSpin,"; NBasis=",NBasis
       D=DD(1,:,:,:)
       DO ispin=1,NSpin
!         DO i=1,NBasis
!           DO j=i,NBasis
!             !write(*,'(I2,A,I2)'),i,",",j
!             D(ispin,j,i)=D(ispin,i,j) ! UNCOMMENTED ON 2017-04-14 BECAUSE OUTPUT DENSITIES WERE NOT SELF-TRANSPOSE.
!           END DO
!           Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
!         END DO
!         D(iSpin,:,:)=0.5*(D(iSpin,:,:)+transpose(D(iSpin,:,:))) ! SYMMETRIZE. ADDED ON 2017-04-14.

         SingleSpinDens(:,:)=D(iSpin,:,:)
         DO i=1,NBasis
           DO j=1,NBasis
             !write(*,'(A,I2,A,I2,A,F8.4)')" SingleSpinDens(",i,",",j,")=",SingleSpinDens(i,j)
             !write(*,'(A,I2,A,I2,A,F8.4)')" S(",i,",",j,")=",S(i,j)
           END DO
           TraceChargeNoOv(iSpin)=TraceChargeNoOv(iSpin)+SingleSpinDens(i,i)
         END DO
         write(*,'(A,I1,A,F8.4)')" TraceChargeNoOv(iSpin=",iSpin,")=",TraceChargeNoOv(iSpin)
         PtimesS(iSpin,:,:)=MATMUL(SingleSpinDens,S)
         DO i=1,NBasis
!           DO j=1,NBasis
!             write(*,'(A,I2,A,I2,A,I2,A,F8.4)')" D(",ispin,",",i,",",j,")=",D(ispin,i,j)
!           END DO
           TraceCharge(iSpin)=TraceCharge(iSpin)+PtimesS(iSpin,i,i)
         END DO
         write(*,'(A,I1,A,F8.4)')" TraceCharge(iSpin=",iSpin,")=",TraceCharge(iSpin)
       END DO
       write(*,'(A,F8.4)')" TraceChargeNoOv(alpha+beta)=",SUM(TraceChargeNoOv(:))
       write(*,'(A,F8.4)')" TraceCharge(alpha+beta)=",SUM(TraceCharge(:))
       Deallocate(PtimesS,SingleSpinDens)
       write(*,'(A)')"CORRECT DEALLOCATION OF PtimesS and SingleSpinDens"
         Write(*,'(A)')"*************************************************************"
         Write(*,'(A)')"******** NO NEED TO NORMALIZE D(A/B) TO inputN(A/B)E ********"
         Write(*,'(A)')"*************************************************************"
!       if(TraceCharge(1)/=inputNAE)then
!
!         Write(*,'(A)')"************************************************************"
!         Write(*,'(A)')"**************** NORMALIZING DA TO inputNAE ****************"
!         Write(*,'(A)')"************************************************************"
!         Write(*,'(A,I4)')"inputNAE = ",inputNAE
!         Write(*,'(A,F16.10)')"inputNAE/TraceCharge(1) = ",(inputNAE/TraceCharge(1))
!         Write(*,'(A)')"Before: "
!         do i=1,NBasis
!           Write(*, '(10F8.4)')( D(1,i,j) ,j=1,NBasis)
!         end do
!         if(TraceCharge(1)/=0.0d0)D(1,:,:)=D(1,:,:)*(inputNAE/TraceCharge(1))
!         Write(*,'(A)')"After: "
!         do i=1,NBasis
!           Write(*, '(10F8.4)')( D(1,i,j) ,j=1,NBasis)
!         end do
!       end if
!       if(NSpin==2)then
!         if(TraceCharge(2)/=inputNBE)then
!           Write(*,'(A)')"************************************************************"
!           Write(*,'(A)')"**************** NORMALIZING DB TO inputNBE ****************"
!           Write(*,'(A)')"************************************************************"
!           Write(*,'(A,I4)')"inputNBE = ",inputNBE
!           Write(*,'(A,F16.10)')"inputNBE/TraceCharge(2) = ",(inputNBE/TraceCharge(2))
!           Write(*,'(A)')"Before: "
!           do i=1,NBasis
!             Write(*, '(10F8.4)')( D(2,i,j) ,j=1,NBasis)
!           end do
!           if(TraceCharge(2)/=0.0d0)D(2,:,:)=D(2,:,:)*(inputNBE/TraceCharge(2))
!           Write(*,'(A)')"After: "
!           do i=1,NBasis
!             Write(*, '(10F8.4)')( D(2,i,j) ,j=1,NBasis)
!           end do
!         end if
!       end if


     END IF
     !write(*,'(A,F8.4)')"CHARGE AFTER PULAY ACCELERATION = ",Charge(:)
     ! End Pulay accelaration

     write(*,'(A)')"BEFORE DFTU"
     !     
     ! For DFT+U calculations the damped density matrix has
     ! to be fed back to device module in order to calculate
     ! the correct DFT+U potential in the next step 
     !  
     if( DFTU )then
        do ispin=1,NSpin
           do i=1,NBasis
              do j=1,NBasis
                 call SetDevDensMat( ispin, i, j, D(ispin,i,j) )
              end do
           end do
        end do
     end if
  
     write(*,'(A)')"BEFORE CALCULATE ntot"
     ntot=0
     DO i=1,GetNShell()-1
        IF (GetAtm4Sh(i).NE.GetAtm4Sh(i+1)) THEN
           NAO(GetAtm4Sh(i))=Get1stAO4Sh(i+1)-(ntot+1)
           ntot=ntot+NAO(GetAtm4Sh(i))
        ENDIF
     ENDDO
     NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot
     write(*,*)"CALCULATED NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot = ",NAO(GetAtm4Sh(GetNShell()))

     write(*,*)"BEFORE WRITING P.name.dat if (NCycLeadsOn  < 1000)"
     ! Write density matrix to file P.' // trim(jobname) // '.dat
!     if (NCycLeadsOn  < 1000) then
     if (NCycLeadsOn  < inputNCycles) then
        OPEN(ifu_dm,file=densitymatrix,status='unknown')
        WRITE(ifu_dm,*)DevShift()
        DO ispin=1,NSpin
           i=0
           do iAtom=1,GetNAtoms()
              do ii=1,NAO(iAtom)
                 i=i+1
                 j=0
                 do jAtom=1,GetNAtoms() 
                    do jj=1,NAO(jAtom)
                       j=j+1
                       if (i >= j) WRITE(ifu_dm,'(i2,2i6,e18.10,i5,i6,i5,i6)')ispin,i,j,D(ispin,i,j),iAtom,ii,jAtom,jj
                    end do
                 end do
              end do
           end do
        end do
        close(ifu_dm)
     end if
!     write(*,*)"AFTER WRITING P.name.dat if (NCycLeadsOn  < 1000)"
     write(*,*)"AFTER WRITING P.name.dat if (NCycLeadsOn  < inputNCycles)"
    
     if (NCycLeadsOn == 1 .and. PFix .and. .not. FInit) then
       write(*,*)"ENTERED if (NCycLeadsOn == 1 .and. PFix .and. .not. FInit)"
        CALL ReadDensMat(densitymatrix)        
        do ispin=1,NSpin
           do i=1,NBasis
              do j=1,NBasis
                 D(ispin,i,j)=DevDensMat(ispin,i,j)
              end do
              Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
           end do
        end do
     end if

  end if

  !*******************************************
  ! RETURN DENSITY MATRIX TO GAUSSIAN VIA RWF
  !*******************************************
  
  ! Transform density matrix D 
  ! to lower triangular form
  ! and return to Gaussian via RWF
  acount = 1
  DO i=1,NBasis
     DO j=1,i
        if(FMixing)then
           !MT(1,acount)=DevDensMat(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
           !IF(UHF) MT(2,acount) = DevDensMat(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
        else
           !MT(1,acount)=D(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
           !IF(UHF) MT(2,acount) = D(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
        end if
        acount = acount +1
     ENDDO
  ENDDO
!  CALL FileIO(one,-IRwPA,ntt,MT(1,:),zero)
!  IF (UHF) CALL FileIO(one,-IRwPB,ntt,MT(2,:),zero)
  
  ! In last cycle (1000): clean up and return to Gaussian 
  ! without returning density matrix
!  IF (jcycle.EQ.1000) THEN
  IF (jcycle.EQ.inputNCycles) THEN
     CLOSE(ifu_ham)
     CLOSE(ifu_tra)
     CLOSE(ifu_dos)
     CLOSE(ifu_mul)
     CLOSE(ifu_red)
     
     print *, " ----------------- End of ANT.G09 ------------------- "
     !DEALLOCATE( MT, F, STAT=AllocErr ) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
     !IF( AllocErr /= 0 ) THEN
     !   PRINT *, "ESF/Deallocation error for MT(:,:),F(:,:,:)"
     !   STOP
     !END IF
     ! I THINK NOT NECESSARY TO DEALLOCATE HERE, BECAUSE DEALLOCATION IN CleanUpDevice
     DEALLOCATE( F, STAT=AllocErr ) ! ADDED BY C.SALGADO ON 2017-03-22 TO DEALOCATE F WHILE REMOVING MT(NSpin, ntt) APPEARANCES.
     Write(*,*)"AllocErr = ",AllocErr
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Deallocation error for F(:,:,:)"
        STOP
     END IF
     Write(*,*) "DEALLOCATED F(:,:,:)"
!     if(.not.FMixing)then
     if(.not.FMixing)then
!       DEALLOCATE( D, DD, DDD, STAT=AllocErr ) ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
       DEALLOCATE( DD, DDD, STAT=AllocErr )
       Write(*,*)"AllocErr = ",AllocErr
       !IF( AllocErr /= 0 ) THEN
       IF(.NOT.(AllocErr == 0)) THEN
!         PRINT *, "ESF/Deallocation error for D(:,:,:), DD(:,:,:,:), DDD(:,:,:)" ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
         PRINT *, "ESF/Deallocation error for DD(:,:,:,:), DDD(:,:,:)"
         STOP
       END IF
       Write(*,*) "DEALLOCATED DD(:,:,:), DDD(:,:,:)"
     end if
     Write(*,*)"CALL CleanUpDevice"
     CALL CleanUpDevice
     Write(*,*)"RETURNING FROM ANT AFTER CleanUpDevice"
     RETURN
  END IF
  
  !*********
  ! BYE-BYE
  !*********
  !DA(:,:) = D(1,:,:)
  !if(NSpin==2)DB(:,:) = D(NSpin,:,:)
  do iSpin=1,NSpin
    do i=1,NBasis
      do j=1,NBasis
      if(iSpin.eq.1)DA(i,j) = D(iSpin,i,j)
      if(iSpin.eq.1)DB(i,j) = D(iSpin,i,j)
      end do
    end do
  end do
  Write(*,'(A)')"WRITE DA OUTPUT FROM ANTC"
  do i=1,NBasis
    Write(*, '(10F8.4)')( DA(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"END WRITE DA OUTPUT FROM ANTC"
  Write(*,'(A)')"WRITE DB OUTPUT FROM ANTC"
  do i=1,NBasis
    Write(*, '(10F8.4)')( DB(i,j) ,j=1,NBasis)
  end do
  Write(*,'(A)')"END WRITE DB OUTPUT FROM ANTC"
  PRINT *, "RETURNING FROM ANT BYE-BYE"
  RETURN
END SUBROUTINE antc


  SUBROUTINE ANT(UHF,JCycle,inputNCycles,inputjobname,inputjobname_len,D,pivHam,pivFock,pivCoulomb,pivExchange,pivOver,&
                 NBasis,inputNSpin,inputNAtoms,inputNShell,inputJAN,inputAOS,inputShellT,inputShellC,&
                 inputNAE,inputNBE,inputNE,inputIAN,inputAtmCo,inputAtmChg,&
                 denerrj,Crit,ANTOn)

!**********************************************************************************************************************
!* Interface subroutine with Gaussian                                                                                 *
!**********************************************************************************************************************
  USE Parameters, ONLY: SL, SwOffSPL, alpha, Read_Parameters, Write_Parameters, NSpinLock, npulay
  USE Parameters, ONLY: ChargeAcc,ChargeA,FermiAcc,FermiA,PAcc,PA,FullAcc,RedTransmB,RedTransmE,ElType,LDOS_Beg,LDOS_End
  USE Parameters, ONLY: Mulliken, Hamilton, PFix, DFTU, FMixing, IntEnergy, DiagFock, SpinMu
  USE constants, ONLY: Hart
  USE preproc
  USE device, ONLY: InitDevice, DevFockMat, DevDensMat, ReadDensMat, LeadsOn, DevShift, SwitchOnLeads, &
       SwitchOnEvaluation, SwitchOnSecant, SwitchOffSecant, SwitchOnSpinLock, SwitchOffSpinLock, &
       SwitchOnChargeCntr, SwitchOffChargeCntr, transport, CleanUpDevice, SetDevDensMat, SetDevFockMat, ReadFockMat
!  USE g09Common, ONLY: GetNShell, GetAtm4Sh, Get1stAO4Sh, GetNBasis, GetAN, GetAtmChg, GetAtmCo, GetNAtoms
  USE g09Common
  use ANTCommon
  use util
  use, intrinsic :: ISO_C_BINDING
  IMPLICIT NONE

  ! dummy arguments
  logical, INTENT(in)      :: UHF
  integer, INTENT(inout) :: JCycle, inputNCycles
  integer(C_INT), intent(IN), VALUE :: inputjobname_len
  character(inputjobname_len,kind=C_CHAR), intent(IN) :: inputjobname
  integer, INTENT(in)    :: NBasis,inputNSpin,inputNAtoms,inputNShell
  real,DIMENSION(inputNSpin,NBasis,NBasis),INTENT(inout) :: D !, pivDens
  real,DIMENSION(NBasis,NBasis),INTENT(in) :: pivHam, pivCoulomb, pivOver
  real,DIMENSION(inputNSpin,NBasis,NBasis),INTENT(in) :: pivFock, pivExchange
  integer, DIMENSION(inputNAtoms), INTENT(in):: inputIAN
  integer, DIMENSION(inputNShell), INTENT(in):: inputJAN, inputAOS, inputShellT, inputShellC
  integer, INTENT(in)    :: inputNAE,inputNBE,inputNE
  real, DIMENSION(inputNAtoms), INTENT(in) :: inputAtmChg
  real, DIMENSION(3,inputNAtoms), INTENT(in) :: inputAtmCo
  real,INTENT(in)        :: denerrj, Crit
!  integer, INTENT(in)    :: IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig!,NBasis
  integer    :: IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig!,NBasis
  integer :: inputNBasis ! LATER OVERWRITTEN BY NBasis, ADDED BY C.SALGADO ON 2017-03-23.

!  integer, DIMENSION(NBasis), INTENT(in)::

  logical,INTENT(inout)    :: ANTOn
  logical :: ADDP
  logical, SAVE :: FInit

  ! local variables
  real    :: val, density, fock, exspin, norma
  integer :: ic, ndim, i, j, info, ii, jj, ispin, acount, AllocErr, ios ,k, iatom, jatom, ntot
  integer, DIMENSION(MaxAtm) :: NAO

  CHARACTER(len=80), SAVE :: densitymatrix, fockmatrix

  integer, SAVE :: isw, NCycLeadsOn, NSpin
  integer, SAVE :: ntt, len
  integer, PARAMETER :: zero=0, two=2, one=1

  real, DIMENSION(:),ALLOCATABLE   :: TS
  real, DIMENSION(:,:),ALLOCATABLE :: S
  real, DIMENSION(:,:,:),ALLOCATABLE :: PtimesS

  !
  ! Matrices in lower trinagular form (alpha and beta)
  ! for communication with gaussian
  !
!  real, DIMENSION(:,:),ALLOCATABLE,SAVE ::  MT
  !
  ! Fock and density matrices in regular form
  !
  real, DIMENSION(:,:,:),ALLOCATABLE,SAVE :: F!, D ! COMMENTED MATRIX D ON 2017-03-22.
  real, DIMENSION(:,:,:,:),ALLOCATABLE,SAVE ::  DD
  real, DIMENSION(:,:,:),ALLOCATABLE,SAVE ::  DDD
  real, DIMENSION(npulay+1,npulay+1) ::  b
  real, DIMENSION(npulay+1) ::  c
  integer, DIMENSION(npulay+1) ::  ipiv

  real :: fmix
  real, DIMENSION(2) :: TraceCharge, TraceChargeNoOv
  real, DIMENSION(:,:),ALLOCATABLE ::  SingleSpinDens
  character(80) strfmt
  ! VERY IMPORTANT ADDED BY C.SALGADO ON 2017-03-23 TO OVERWRITE inputNBasis by NBasis.
  Write(*,*)"**********************************************"
  Write(*,*)"**************** ENTERED ANTC ****************"
  Write(*,*)"**********************************************"
  inputNBasis = NBasis
  NSpin = inputNSpin


!  do iSpin=1,NSpin
!    do i=1,NBasis
!      do j=1,NBasis
!        if(isnan(D(iSpin,i,j)))D(iSpin,i,j)=0.0d0
!        !write(*,'(A,I1,A,I2,A,I2,A,F14.6)')"D(",iSpin,",",i,",",j,") = ",D(iSpin,i,j)
!      end do
!    end do
!  end do
  !*******************************************************************
  ! Before first cycle: Initialize module device, allocate memory etc
  !*******************************************************************
  IF(JCycle.EQ.0) THEN

     NSpin = 1
     IF(UHF) NSpin = 2
     ntt = (NBasis*(NBasis+1))/2

     ANTOn = .FALSE.

     PRINT *
     PRINT *, " ****************************************************************** "
     PRINT *, " ***                                                            *** "
     PRINT *, " ***                      A l i c a n t e                       *** "
     PRINT *, " ***                      N a n o                               *** "
     PRINT *, " ***                      T r a n s p o r t                     *** "
     PRINT *, " ***                      G 0 9                                 *** "
     PRINT *, " ***                                                            *** "
     PRINT *, " ****************************************************************** "
     PRINT *, " ***                     Version: 2.3.4                         *** "
     PRINT *, " ****************************************************************** "
     PRINT *, " *  Copyright (c) by                                              * "
     PRINT *, " *                                                                * "
     PRINT *, " *  Juan Jose Palacios (1)                                        * "
     PRINT *, " *  David Jacob (2)                                               * "
     PRINT *, " *  Maria Soriano (1,5)                                           * "
     PRINT *, " *  Angel J. Perez-Jimenez (3)                                    * "
     PRINT *, " *  Emilio SanFabian (3)                                          * "
     PRINT *, " *  Jose Antonio Antonio Verges (4)                               * "
     PRINT *, " *  Enrique Louis (5)                                             * "
     PRINT *, " *                                                                * "
     PRINT *, " * (1) Departamento de Fisica de la Materia Condensada            * "
     PRINT *, " *     Universidad Autonoma de Madrid                             * "
     PRINT *, " *     28049 Madrid (SPAIN)                                       * "
     PRINT *, " * (2) Theory Department                                          * "
     PRINT *, " *     Max-Planck-Institute for Microstructure Physics            * "
     PRINT *, " *     Halle, 06120 (GERMANY)                                     * "
     PRINT *, " * (3) Departamento de Quimica Fisica                             * "
     PRINT *, " *     Universidad de Alicante                                    * "
     PRINT *, " *     03690 Alicante (SPAIN)                                     * "
     PRINT *, " * (4) Insto. de Ciencias de Materiales de Madrid (ICMM)          * "
     PRINT *, " *     Consejo Superior de Investigacion y Ciencia (CSIC)         * "
     PRINT *, " *     28049 Madrid (SPAIN)                                       * "
     PRINT *, " * (5) Departamento de Fisica Aplicada                            * "
     PRINT *, " *     Universidad de Alicante                                    * "
     PRINT *, " *     03690 Alicante (SPAIN)                                     * "
     PRINT *, " *                                                                * "
     PRINT *, " ****************************************************************** "
     PRINT *
     PRINT *, " ****************************************************************** "
     PRINT *, " *                                                                * "
     PRINT *, " *              Initializing ANT.G09 modules                      * "
     PRINT *, " *                                                                * "
     PRINT *, " ****************************************************************** "
     PRINT *

    !*********************************************************************
    !**** Before first cycle: Set all variables formerly in g09Common ****
    !*********************************************************************
    PRINT *, "Before first cycle: Set all variables formerly in g09Common"
    CALL SetNAtoms( inputNAtoms )
    CALL SetNShell( inputNShell )
    CALL SetAtm4Sh(inputNShell, inputJAN -1) ! IMPORTANT TO SUBSTRACT -1 BECAUSE OF THE C++ INDEX OF ARRAYS.
    CALL Set1stAO4Sh( inputNShell, inputAOs -1) ! IMPORTANT TO SUBSTRACT -1 BECAUSE OF THE C++ INDEX OF ARRAYS.
    CALL SetShellT( inputNShell, inputShellT )
    CALL SetShellC( inputNShell, inputShellC )
    CALL SetNAE( inputNAE )
    CALL SetNBE( inputNBE )
    CALL SetNE( inputNE )
    CALL SetNBasis( inputNBasis )
    CALL SetAN( inputNAtoms, inputIAN )
    CALL SetAtmChg( inputNAtoms, inputAtmChg )
    CALL SetAtmCo( inputNAtoms, inputAtmCo )

    do i=1,inputNAtoms
      write(*, '(A,I2,A,10F8.4)')"inputAtmCo(1:3,",i,") = ",( inputAtmCo(j,i) ,j=1,3)
    end do
    PRINT *, "Before first cycle: End set all variables formerly in g09Common"

    !*********************************************************************
    ! Before first cycle: End setting all variables formerly in g09Common
    !*********************************************************************
     ! Reading parameter file

! THE COMMENTED CODE BELOW WAS USED TO READ FILE CALLED "NAME"
!     OPEN(UNIT=ifu_nam,FILE="name",IOSTAT=ios,STATUS='OLD')
!     IF( ios == 0 ) THEN
!        READ(ifu_nam,*) jobname
!        OPEN(UNIT=ifu_ini,FILE=trim(jobname)//'.ini',IOSTAT=ios,STATUS='OLD')
!        IF( ios == 0 ) THEN
!          CALL read_parameters( ifu_log, ifu_ini, ios, jobname )
!        ELSE
!          WRITE(ifu_log,*) " "
!          WRITE(ifu_log,*) "*** No parameter file found. Using default values."
!          WRITE(ifu_log,*) " "
!        END IF
!        CLOSE(ifu_ini)
!     ELSE
!        WRITE(ifu_log,*) " "
!        WRITE(ifu_log,*) "*** No name file found. Using default values."
!        WRITE(ifu_log,*) " "
!     END IF
!     CLOSE(UNIT=ifu_nam)
!     jobname = C_F_STRING(inputjobname)
     jobname = c_to_f_string(inputjobname)
     OPEN(UNIT=ifu_ini,FILE=trim(jobname)//'.ini',IOSTAT=ios,STATUS='OLD')
     IF( ios == 0 ) THEN
       CALL read_parameters( ifu_log, ifu_ini, ios, jobname )
     ELSE
       WRITE(ifu_log,*) " "
       WRITE(ifu_log,*) "*** No parameter file found. Using default values."
       WRITE(ifu_log,*) " "
     END IF
     CLOSE(ifu_ini)

     ant1dname = jobname
     call locase(ant1dname)

     CALL write_parameters( ifu_log )

     !Opening writting files
     OPEN(ifu_xyz,file=trim(jobname)//'_BL.xyz',status='unknown')
     IF (ElType(1) /= 'GHOST' .and. ElType(2) /= 'GHOST') OPEN(ifu_tra,file='T.'//trim(jobname)//'.dat',status='unknown')
     IF (RedTransmB < RedTransmE) OPEN(ifu_red,file='t.'//trim(jobname)//'.dat',status='unknown')
     IF (Hamilton) OPEN(ifu_ham,file='V.'//trim(jobname)//'.dat',status='unknown')
     IF (DiagFock) OPEN(ifu_diagham,file='DF.'//trim(jobname)//'.dat',status='unknown')
     IF (Mulliken) OPEN(ifu_mul,file='Q.'//trim(jobname)//'.dat',status='unknown')
     IF (LDOS_Beg <= LDOS_End) OPEN(ifu_dos,file='DOS.'//trim(jobname)//'.dat',status='unknown')

     ! Creating density matrix file name
     densitymatrix='P.' // trim(jobname) // '.dat'
     ! Creating Fock matrix file name
     fockmatrix='F.' // trim(jobname) // '.dat'
     !
     ! Allocate memory for dynamic arrays
     !
!     ALLOCATE( MT( NSpin, ntt ), STAT=AllocErr ) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!     IF( AllocErr /= 0 ) THEN ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!        PRINT *, "ESF/Error: could not allocate memory for MT(:)" ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!        STOP ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
!     END IF ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
     PRINT *, "(NSpin,NBasis)=(",NSpin,",",NBasis,")"
     ALLOCATE( F(NSpin,NBasis,NBasis), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for F(:,:,:)"
        STOP
     END IF
!     DO iSpin=1,NSpin
!       DO i=1,NBasis
!          DO j=1,NBasis
!            F(iSpin,i,j) = pivOver(i,j) ! ADDED AND COMMENTED BY C.SALGADO ON 2017-03-25. FOCK MATRIX IS SET VIA F.name.dat.
!          END DO
!       END DO
!     END DO
      !FOCK MATRIX HD(iSpin,i,j) will be initialized in InitDevice by doing HD = F.

     if(.not. FMixing )then
        !ALLOCATE( D(NSpin,NBasis,NBasis), DD(npulay+1,NSpin,NBasis,NBasis),DDD(npulay,NSpin,NBasis*NBasis), STAT=AllocErr )
        ALLOCATE( DD(npulay+1,NSpin,NBasis,NBasis),DDD(npulay,NSpin,NBasis*NBasis), STAT=AllocErr ) ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
        IF( AllocErr /= 0 ) THEN
           !PRINT *, "ESF/Error: could not allocate memory for D(:,:,:), DD(:,:,:,:),DDD(:,:,:)"
           PRINT *, "ESF/Error: could not allocate memory for DD(:,:,:,:),DDD(:,:,:)" ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
           STOP
        END IF
        DD=0.0
        DDD=0.0
     end if

     !
     ! Obtain overlap matrix in lower triangular form
     ! and transform to regular form matrix
     !
     ALLOCATE( TS(ntt), S(NBasis,NBasis), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for S(:,:), TS(:)"
        STOP
     END IF

!     CALL FileIO(two,-IRwS1,ntt,TS,zero)
     acount = 1
     DO i=1,NBasis
        DO j=1,i
!           S(i,j) = TS(acount)
!           S(j,i) = TS(acount)
           S(i,j) = pivOver(i,j)
           S(j,i) = pivOver(j,i)
           acount = acount+1
        END DO
     END DO

     !
     ! Initialize module device
     !
     CALL InitDevice( NBasis, UHF, S )

     CLOSE(ifu_xyz)

     if(FMixing)OPEN(ifu_fm,file=fockmatrix,iostat=ios,status='old')
     if(.not.FMixing)OPEN(ifu_dm,file=densitymatrix,iostat=ios,status='old')
     IF (ios /= 0) THEN
        FInit=.false.
     ELSE
        FInit=.true.
        if(FMixing)CLOSE(ifu_fm)
        if(.not.FMixing)CLOSE(ifu_dm)
     END IF

     IF( FInit )THEN
        !
        ! Initialization from file F.dat or P.dat
        ! File F.name.dat just never exists, so this is not the usual procedure to get the Hamiltonian.
        !
        if(FMixing)then
           CALL ReadFockMat(fockmatrix)
           DO i=1,NBasis
              DO j=1,NBasis
                 !F(1,i,j) = DevFockMat(1,i,j)
                 !if(UHF) F(2,i,j) = DevFockMat(2,i,j)
                 F(1,i,j) = pivFock(1,i,j)
                 if(UHF) F(2,i,j) = pivFock(2,i,j)
              end DO
           end DO
        else
           CALL ReadDensMat(densitymatrix)
           !
           ! Transform to triangular form
           !
           acount = 1
           DO i=1,NBasis
              DO j=1,i
                 D(1,i,j)=DevDensMat(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-25 BECAUSE D IS INOUT OF ANT.
                 !MT(1,acount)=D(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
                 IF(UHF) D(2,i,j) = DevDensMat(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-25 BECAUSE D IS INOUT OF ANT.
                 !IF(UHF) MT(2,acount) = D(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
                 acount = acount +1
              ENDDO
           ENDDO
           !
           ! PASS DENSITY MATRIX TO GAUSSIAN VIA RWF
           !
!           CALL FileIO(one,-IRwPA,ntt,MT(1,:),zero)
!           IF (UHF) CALL FileIO(one,-IRwPB,ntt,MT(2,:),zero)
        end if
     ENDIF

     DEALLOCATE( S, TS, STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not DEAllocate memory for S(:,:), TS(:)"
        STOP
     END IF

     PRINT *, "RETURNING FROM ANT(...)"
     RETURN
  END IF

  IF( .NOT. LeadsOn() )THEN
     !
     ! Connect Leads in first cycle
     ! when reinitialized
     !
     IF( FInit )THEN
        ANTOn = .true.
        CALL SwitchOnLeads
        NCycLeadsOn = 0
        !
        ! Or if SL criterion is met:
        !
     ELSEIF( denerrj .LE. SL ) THEN
        ANTOn = .true.
        Call SwitchOnLeads
        NCycLeadsOn = 0
        if(.not.FMixing)then
           !
           ! Obtain density matrix from gaussian only in first step
           ! with leads switched on and if not initialized from file
           !
!           CALL FileIO(two,-IRwPA,ntt,MT(1,:),zero)
!           IF(UHF) CALL FileIO(two,-IRwPB,ntt,MT(2,:),zero)

           acount = 1
           DO i=1,NBasis
              DO j=1,i
!                 D(1,i,j)=MT(1,acount)
!                 D(1,j,i)=D(1,i,j)
                 !D(1,i,j)=pivDens(i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                 !D(1,j,i)=D(1,i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                 call SetDevDensMat( 1, i, j, D(1,i,j) )
                 call SetDevDensMat( 1, j, i, D(1,j,i) )
                 !call SetDevFockMat( 1, i, j, F(1,i,j) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                 !call SetDevFockMat( 1, i, j, F(1,j,i) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                 IF (UHF) THEN
!                    D(2,i,j)=MT(2,acount)
!                    D(2,j,i)=D(2,i,j)
                    !D(2,i,j)=pivDens(i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                    !D(2,j,i)=D(2,i,j) ! D HAS INOUT CHARACTER AFTER 2017-03-22. DO NOT LOAD FROM pivDens.
                    call SetDevDensMat( 2, i, j, D(2,i,j) )
                    call SetDevDensMat( 2, j, i, D(2,j,i) )
                    !call SetDevFockMat( 2, i, j, F(2,i,j) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                    !call SetDevFockMat( 2, i, j, F(2,j,i) ) ! ADDED BY C. SALGADO ON 2017-03-25 TO INITIALIZE FOCK MATRIX IN DEVICE HD(iSpin,i,j) WITHOUT BETHE LATTICE.
                 ENDIF
                 acount = acount+1
              ENDDO
           ENDDO
        end if
     END IF
  END IF

  !IF( jcycle.EQ.1000 )THEN
  IF( jcycle.EQ.inputNCycles )THEN
     !NCycLeadsOn = 999
     NCycLeadsOn = inputNCycles-1
     CALL SwitchOnEvaluation()
     CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. .not. FullAcc .and. jcycle > 2) THEN
     IF( denerrj .LT. Crit*10) THEN
        ChargeAcc = ChargeA
        FermiAcc = FermiA
        PAcc = PA
     ELSE IF( denerrj .LT. Crit*100) THEN
        ChargeAcc = ChargeA*10.0
        FermiAcc = FermiA*10.0
        PAcc = PA*10.0
     ELSE IF( denerrj .LT. Crit*1000) THEN
        ChargeAcc = ChargeA*100.0
        FermiAcc = FermiA*100.0
        PAcc = PA*100.0
     ELSE
        ChargeAcc = ChargeA*1000.0
        FermiAcc = FermiA*1000.0
        PAcc = PA*1000.0
     END IF
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. FullAcc) THEN
     ChargeAcc = ChargeA
     FermiAcc = FermiA
     PAcc = PA
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. .not. FullAcc .and. jcycle <= 2) THEN
     ChargeAcc = ChargeA*10.0
     FermiAcc = FermiA*10.0
     PAcc = PA*10.0
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE
     !
     ! If leads are not switched on or if not in cycle 1000 return to Gaussian
     ! without calculation of density matrix with leads connected
     !
     RETURN
  END IF


  !*********************************************
  ! COMPUTE DENSITY MATRIX WITH LEADS CONNECTED
  !*********************************************

  WRITE(ifu_log,*) '-------------------------------------------------------------------------'
  WRITE(ifu_log,*) 'Computing the density matrix at cycle:', JCycle
  WRITE(ifu_log,*) '-------------------------------------------------------------------------'

  ! Obtain Fock matrix from Gaussian RWF
!  CALL FileIO(two,-IRwFA,ntt,MT(1,:),zero)
!  IF(UHF) CALL FileIO(two,-IRwFB,ntt,MT(2,:),zero)

  fmix=1.0d0
  ! For Fock matrix mixing
  if(FMixing.and.NCycLeadsOn>0) fmix=alpha
  !if(NCycLeadsOn>0) fmix=alpha
  if(FMixing.and.FINIT.and.NCycLeadsOn.eq.0) fmix=0.0d0
  ! Transform from lower trinagular form to regular form
  !print *, "NCycLeadsOn=", NCycLeadsOn
  !print *, "fmix =", fmix
  acount = 1
  Write(*,*)"*********************************************************************************"
  Write(*,*)"***************************** UPDATING FOCK MATRIX ******************************"
  Write(*,*)"*********************************************************************************"
  DO i=1,NBasis
     DO j=1,i
        ! Mixing with old Fock matrix if fmix<1
!        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*MT(1,acount)
!        F(1,j,i)=F(1,i,j)
!        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*(pivHam(i,j)+pivFock(i,j)) ! COMMENTED BY C.SALGADO TO USE F(iSpin,:,:) AS J(Coulomb)+K(Exchange), removing pivHam.
!        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*(pivFock(i,j))
!        F(1,i,j)=Hart*pivFock(i,j)
!        F(1,i,j)=Hart*(pivHam(i,j)+pivFock(i,j))
        F(1,i,j)=(pivHam(i,j)+pivFock(1,i,j))
!        F(1,i,j)=pivHam(i,j)
        F(1,j,i)=F(1,i,j)
        IF (UHF) THEN
!           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*MT(2,acount)
!           F(2,j,i)=F(2,i,j)
!           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*(pivHam(i,j)+pivFock(i,j)) ! COMMENTED BY C.SALGADO TO USE F(iSpin,:,:) AS J(Coulomb)+K(Exchange), removing pivHam.
!           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*(pivFock(i,j))
!           F(2,i,j)=Hart*pivFock(i,j)
!           F(2,i,j)=Hart*(pivHam(i,j)+pivFock(i,j))
           F(2,i,j)=(pivHam(i,j)+pivFock(2,i,j))
!           F(2,i,j)=pivHam(i,j)
           F(2,j,i)=F(2,i,j)
        ENDIF
        acount = acount +1
     ENDDO
  ENDDO
  write(strfmt,'(a,i0,a)') '(A,I1,A,',size(F,2),'(1x,F12.6))'
  write(ifu_log,'(A)')strfmt
  do iSpin=1,NSpin
    do i=1,NBasis
      write(*, '(10F12.6)')( F(iSpin,i,j) ,j=1,NBasis)
    end do
    !write(ifu_log, strfmt)"F(",iSpin,")",F(iSpin,:,:)
  end do
  Write(*,*)"*********************************************************************************"
  Write(*,*)"*************************** END UPDATING FOCK MATRIX ****************************"
  Write(*,*)"*********************************************************************************"

  if(FMixing)then
     ntot=0
     DO i=1,GetNShell()-1
        IF (GetAtm4Sh(i).NE.GetAtm4Sh(i+1)) THEN
           NAO(GetAtm4Sh(i))=Get1stAO4Sh(i+1)-(ntot+1)
           ntot=ntot+NAO(GetAtm4Sh(i))
        ENDIF
     ENDDO
     NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot
     ! Write density matrix to file F.' // trim(jobname) // '.dat
!     if (NCycLeadsOn  < 1000) then
     if (NCycLeadsOn  < inputNCycles) then
        OPEN(ifu_fm,file=fockmatrix,status='unknown')
        WRITE(ifu_fm,*)DevShift()
        DO ispin=1,NSpin
           i=0
           do iAtom=1,GetNAtoms()
              do ii=1,NAO(iAtom)
                 i=i+1
                 j=0
                 do jAtom=1,GetNAtoms()
                    do jj=1,NAO(jAtom)
                       j=j+1
                       if (i >= j) WRITE(ifu_fm,'(i2,2i6,e18.10,i5,i6,i5,i6)')ispin,i,j,F(ispin,i,j),iAtom,ii,jAtom,jj
                    end do
                 end do
              end do
           end do
        end do
        close(ifu_fm)
     end if
  end if

  ! Count cycles with leads connected
  NCycLeadsOn = NCycLeadsOn + 1
  !print*,NCycLeadsOn

  ! Turn on charge control every 5 steps in the first cycles
  IF(MOD(NCycLeadsOn-1,10) == 0 .and. NCycLeadsOn <= 21) CALL SwitchOnChargeCntr
  IF(MOD(NCycLeadsOn-1,20) == 0 .and. NCycLeadsOn > 21) CALL SwitchOnChargeCntr

  ! Call subroutine that solves transport problem
  CALL Transport(F,ADDP)
  ! THE CALL TO DEVDENSMAT BELOW WAS ADDED BY C.SALGADO ON 2016-03-24.
!  Write(*,'(A)')"OVERWRITE D(iSpin,i,j) AFTER TRANSPORT"
  Write(*,'(A)')"DISPLAY D(iSpin,i,j) AFTER TRANSPORT"
  DO ispin=1,NSpin
    DO i=1,NBasis
!      DO j=1,NBasis
!        D(iSpin,i,j)=DevDensMat(iSpin,i,j)
!        D(iSpin,j,i)=D(iSpin,i,j)
!      END DO
      Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
    END DO
  END DO
  Write(*,'(A)')"END D(iSpin,i,j) AFTER TRANSPORT"

  CALL SwitchOffChargeCntr

  IF( SL <= 0.0d0 ) alpha = 1.0d0

  !
  ! Pulay accelaration for density matrix mixing (FMixing == false)
  !
  if(.not.FMixing)then

     DD(1,:,:,:)=D(:,:,:)

     DO k=npulay,1,-1
     !DO k=1,npulay ! PUT HERE BY C.SALGADO WITHOUT IDEA.
        DD(k+1,:,:,:)=DD(k,:,:,:)
     END DO

     if (ADDP) then
        DO ispin=1,NSpin
           DO i=1,NBasis
              DO j=1,NBasis
                 DD(1,ispin,i,j)=(1.0d0-alpha)*DD(2,ispin,i,j)+alpha*DevDensMat(ispin,i,j)
              END DO
           END DO
        END DO
     else
        DD(1,:,:,:)=0.5*DD(3,:,:,:)+0.5*DD(2,:,:,:)
     end if

!     write(*,*)"DD(kpulay,ispin,ibas,jbas)"
     DO k=1,npulay
        DO ispin=1,NSpin
           ic=0
           DO i=1,NBasis
              DO j=1,NBasis
                 ic=ic+1
                 DDD(k,ispin,ic)=DD(k,ispin,i,j)-DD(k+1,ispin,i,j)
!                 write(*,'(A,I2,A,I2,A,I2,A,I2,A,F8.4)')"DD(",k,",",ispin,",",i,",",j,")=",DD(k,ispin,i,j)
              END DO
           END DO
        END DO
     END DO
     ! NEXT NESTED LOOP ONLY TO PRINT -DD(kpulay+1,ispin,ibas,jbas)
!     write(*,*)"-DD(kpulay+1,ispin,ibas,jbas)"
!     DO k=1,npulay
!        DO ispin=1,NSpin
!           DO i=1,NBasis
!              DO j=1,NBasis
!                 write(*,'(A,I2,A,I2,A,I2,A,I2,A,F8.4)')"-DD(",k+1,",",ispin,",",i,",",j,")=",DD(k+1,ispin,i,j)
!              END DO
!           END DO
!        END DO
!     END DO
     ! NEXT NESTED LOOP ONLY TO PRINT DDD(kpulay,ispin,ic)
!     write(*,*)"DDD(kpulay,ispin,ic)"
!     DO k=1,npulay
!        DO ispin=1,NSpin
!           DO ic=1,NBasis*NBasis
!             write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"DDD(",k,",",ispin,",",ic,")=",DDD(k,ispin,ic)
!           END DO
!        END DO
!     END DO

     ! Pulay acceleration procedure
     ndim=npulay+1
     if (mod(NCycLeadsOn,ndim)==npulay) then
        print *, "Pulay kick ............."

        D=0.0
        DO ispin=1,NSpin
           b=0.0
           c=0.0
           c(ndim)=-1.0
           !Matrix to obtain the coefficients
!           write(*,*)"ispin;b(k,j) MATRIX TO OBTAIN THE COEFFICIENTS FOR PULAY BEFORE DGETRF"
           DO k=1,ndim
              DO j=1,ndim
                 if (k.lt.ndim.and.j.lt.ndim) then
                    DO i=1,NBasis*NBasis
                       b(k,j)=b(k,j)+DDD(k,ispin,i)*DDD(j,ispin,i)
                    END DO
                 else if (k==ndim.and.j.lt.ndim) then
                    b(k,j)=-1.0
                 else if (j==ndim.and.k.lt.ndim) then
                    b(k,j)=-1.0
                 end if
!                 write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,b(",k,",",j,")=",b(k,j)
              END DO
           END DO
           ! EXTRACTED ZGETRF FROM NUMERIC FOR INTEL MKL.
           !call zgetrf(NDIM,NDIM,Sigma_aux,NDIM,ipiv,info) ! ORIGINAL FOR PGI IN BetheLattice.f90.
           !call zgetrf(NDim,NDim,Sigma_aux,NDim+1,ipiv,info) ! ADDING +1 TO PARAMETER 4 WORKS!!!
           Call DGETRF( ndim, ndim, b, ndim, IPIV, INFO) ! ORIGINAL CODE WORKS WITH GAUSSIAN. FAILS WITH INFO=2 WITH MKL.
           !Call DGETRF( ndim, ndim, b, ndim+1, IPIV, INFO)  ! ADDED +1 TO PARAMETER 4 BECAUSE IT WORKS. CHECK MAIL:  call zgemm
           !Call DGETRF( ndim, ndim, b, ndim, IPIV, INFO )
!           INFO=0

           IF( INFO /= 0 ) THEN
             WRITE(ifu_log,*)'Pulay/Problems using DGETRF in ANT.f90'
             WRITE(ifu_log,*)'INFO=',INFO
!             STOP
           END IF


!           write(*,*)"ispin;b(k,j) MATRIX TO OBTAIN THE COEFFICIENTS FOR PULAY AFTER DGETRS"
!           DO k=1,ndim
!              DO j=1,ndim
!                 write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,b(",k,",",j,")=",b(k,j)
!              END DO
!           END DO
!           write(*,*)"ispin;IPIV(k) pivot indices from DGETRF to DGETRS"
!           DO k=1,ndim
!                 write(*,'(A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,IPIV(",k,")=",c(k)
!           END DO
           ! EXTRACTED ZGETRI FROM NUMERIC FOR INTEL MKL.
           !call zgetri(NDIM,Sigma_aux,NDIM,ipiv,work,4*NDIM,info) ! ORIGINAL FOR PGI IN BetheLattice.f90.
           !call zgetri(NDim,Sigma_aux,NDim+1,ipiv,work,4*NDim+1,info) ! ADDING +1 TO PARAMETERS 3,6 WORKS!!!
           !DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
           Call DGETRS( 'N', ndim, 1, b, ndim, IPIV, c, ndim, INFO) ! ORIGINAL CODE WORKS WITH GAUSSIAN.
           !Call DGETRS( 'N', ndim, 1, transpose(b), ndim, IPIV, c, ndim, INFO) ! ORIGINAL CODE WORKS WITH GAUSSIAN.
           !Call DGETRS( 'N', ndim, 1, b, ndim+1, IPIV, c, 4*ndim+1, INFO) ! ADDED +1 TO PARAMETERS 5,8 BECAUSE IT WORKS. CHECK MAIL:  call zgem


           IF( INFO /= 0 ) THEN
             WRITE(ifu_log,*)'Pulay/Problems using DGETRS in ANT.f90'
             WRITE(ifu_log,*)'INFO=',INFO
!             STOP
           END IF


!           write(*,*)"ispin;b(k,j) MATRIX TO OBTAIN THE COEFFICIENTS FOR PULAY AFTER DGETRS"
!           DO k=1,ndim
!              DO j=1,ndim
!                 write(*,'(A,I2,A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,b(",k,",",j,")=",b(k,j)
!              END DO
!           END DO
!           write(*,*)"ispin;c(k) VECTOR OF THE COEFFICIENTS FOR PULAY AFTER DGETRS"
!           DO k=1,ndim
!                 write(*,'(A,I2,A,I2,A,F8.4)')"iSpin=",iSpin,"; ,c(",k,")=",c(k)
!           END DO
!           write(*,*)"c(k)*DD(kpulay,ispin,ibas,jbas)"
           DO k=1,npulay
              DO i=1,NBasis
                 DO j=1,NBasis
!                    write(*,'(A,I2,A,I2,A,I2,A,I2,A,I2,A,F8.4)')"c(",k,")*DD(",k,",",ispin,",",i,",",j,")=",c(k)*DD(k,ispin,i,j)
                    D(ispin,i,j)=D(ispin,i,j)+c(k)*DD(k,ispin,i,j)
                 END DO
              END DO
           END DO
           DO i=1,NBasis
             DO j=i,NBasis
               D(ispin,j,i)=D(ispin,i,j)
             END DO
           END DO
           DO i=1,NBasis
             !DO j=1,NBasis
             !  !write(*,'(A,I2,A,I2,A,I2,A,F8.4)')" D(",ispin,",",i,",",j,")=",D(ispin,i,j)
             !END DO
             Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
           END DO
        END DO
     ELSE
       write(*,'(A)')"SUMMING CHARGE AFTER PULAY ACCELERATION"
       TraceCharge = (/ 0.0d0, 0.0d0 /)
       TraceChargeNoOv = (/ 0.0d0, 0.0d0 /)
        write(*,'(A,F8.4)')" ZERO TraceCharge(alpha+beta)=",SUM(TraceCharge(:))
       Allocate(PtimesS(NSpin,NBasis,NBasis),&
                SingleSpinDens(NBasis,NBasis),&
                S(NBasis,NBasis))
       S = pivOver

       write(*,'(A,I1,A,I2)')"NSpin=",NSpin,"; NBasis=",NBasis
       D=DD(1,:,:,:)
       DO ispin=1,NSpin
         DO i=1,NBasis
           DO j=i,NBasis
             !write(*,'(I2,A,I2)'),i,",",j
             D(ispin,j,i)=D(ispin,i,j)
           END DO
           Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
         END DO
         SingleSpinDens(:,:)=D(iSpin,:,:)
         DO i=1,NBasis
           DO j=1,NBasis
             !write(*,'(A,I2,A,I2,A,F8.4)')" SingleSpinDens(",i,",",j,")=",SingleSpinDens(i,j)
             !write(*,'(A,I2,A,I2,A,F8.4)')" S(",i,",",j,")=",S(i,j)
           END DO
           TraceChargeNoOv(iSpin)=TraceChargeNoOv(iSpin)+SingleSpinDens(i,i)
         END DO
         write(*,'(A,I1,A,F8.4)')" TraceChargeNoOv(iSpin=",iSpin,")=",TraceChargeNoOv(iSpin)
         PtimesS(iSpin,:,:)=MATMUL(SingleSpinDens,S)
         DO i=1,NBasis
           DO j=1,NBasis
             !write(*,'(A,I2,A,I2,A,I2,A,F8.4)')" D(",ispin,",",i,",",j,")=",D(ispin,i,j)
           END DO
           TraceCharge(iSpin)=TraceCharge(iSpin)+PtimesS(iSpin,i,i)
         END DO
         write(*,'(A,I1,A,F8.4)')" TraceCharge(iSpin=",iSpin,")=",TraceCharge(iSpin)
       END DO
       write(*,'(A,F8.4)')" TraceChargeNoOv(alpha+beta)=",SUM(TraceChargeNoOv(:))
       write(*,'(A,F8.4)')" TraceCharge(alpha+beta)=",SUM(TraceCharge(:))
       Deallocate(PtimesS,SingleSpinDens)
       write(*,'(A)')"CORRECT DEALLOCATION OF PtimesS and SingleSpinDens"
     END IF
     !write(*,'(A,F8.4)')"CHARGE AFTER PULAY ACCELERATION = ",Charge(:)
     ! End Pulay accelaration

     write(*,'(A)')"BEFORE DFTU"
     !
     ! For DFT+U calculations the damped density matrix has
     ! to be fed back to device module in order to calculate
     ! the correct DFT+U potential in the next step
     !
     if( DFTU )then
        do ispin=1,NSpin
           do i=1,NBasis
              do j=1,NBasis
                 call SetDevDensMat( ispin, i, j, D(ispin,i,j) )
              end do
           end do
        end do
     end if

     write(*,'(A)')"BEFORE CALCULATE ntot"
     ntot=0
     DO i=1,GetNShell()-1
        IF (GetAtm4Sh(i).NE.GetAtm4Sh(i+1)) THEN
           NAO(GetAtm4Sh(i))=Get1stAO4Sh(i+1)-(ntot+1)
           ntot=ntot+NAO(GetAtm4Sh(i))
        ENDIF
     ENDDO
     NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot
     write(*,*)"CALCULATED NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot = ",NAO(GetAtm4Sh(GetNShell()))

     write(*,*)"BEFORE WRITING P.name.dat if (NCycLeadsOn  < 1000)"
     ! Write density matrix to file P.' // trim(jobname) // '.dat
!     if (NCycLeadsOn  < 1000) then
     if (NCycLeadsOn  < inputNCycles) then
        OPEN(ifu_dm,file=densitymatrix,status='unknown')
        WRITE(ifu_dm,*)DevShift()
        DO ispin=1,NSpin
           i=0
           do iAtom=1,GetNAtoms()
              do ii=1,NAO(iAtom)
                 i=i+1
                 j=0
                 do jAtom=1,GetNAtoms()
                    do jj=1,NAO(jAtom)
                       j=j+1
                       if (i >= j) WRITE(ifu_dm,'(i2,2i6,e18.10,i5,i6,i5,i6)')ispin,i,j,D(ispin,i,j),iAtom,ii,jAtom,jj
                    end do
                 end do
              end do
           end do
        end do
        close(ifu_dm)
     end if
!     write(*,*)"AFTER WRITING P.name.dat if (NCycLeadsOn  < 1000)"
     write(*,*)"AFTER WRITING P.name.dat if (NCycLeadsOn  < inputNCycles)"

     if (NCycLeadsOn == 1 .and. PFix .and. .not. FInit) then
       write(*,*)"ENTERED if (NCycLeadsOn == 1 .and. PFix .and. .not. FInit)"
        CALL ReadDensMat(densitymatrix)
        do ispin=1,NSpin
           do i=1,NBasis
              do j=1,NBasis
                 D(ispin,i,j)=DevDensMat(ispin,i,j)
              end do
              Write(*, '(10F8.4)')( D(iSpin,i,j) ,j=1,NBasis)
           end do
        end do
     end if

  end if

  !*******************************************
  ! RETURN DENSITY MATRIX TO GAUSSIAN VIA RWF
  !*******************************************

  ! Transform density matrix D
  ! to lower triangular form
  ! and return to Gaussian via RWF
  acount = 1
  DO i=1,NBasis
     DO j=1,i
        if(FMixing)then
           !MT(1,acount)=DevDensMat(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
           !IF(UHF) MT(2,acount) = DevDensMat(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
        else
           !MT(1,acount)=D(1,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
           !IF(UHF) MT(2,acount) = D(2,i,j) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
        end if
        acount = acount +1
     ENDDO
  ENDDO
!  CALL FileIO(one,-IRwPA,ntt,MT(1,:),zero)
!  IF (UHF) CALL FileIO(one,-IRwPB,ntt,MT(2,:),zero)

  ! In last cycle (1000): clean up and return to Gaussian
  ! without returning density matrix
!  IF (jcycle.EQ.1000) THEN
  IF (jcycle.EQ.inputNCycles) THEN
     CLOSE(ifu_ham)
     CLOSE(ifu_tra)
     CLOSE(ifu_dos)
     CLOSE(ifu_mul)
     CLOSE(ifu_red)

     print *, " ----------------- End of ANT.G09 ------------------- "
     !DEALLOCATE( MT, F, STAT=AllocErr ) ! COMMENTED BY C.SALGADO ON 2017-03-22 TO REMOVE MT(NSpin, ntt) APPEARANCES.
     !IF( AllocErr /= 0 ) THEN
     !   PRINT *, "ESF/Deallocation error for MT(:,:),F(:,:,:)"
     !   STOP
     !END IF
     ! I THINK NOT NECESSARY TO DEALLOCATE HERE, BECAUSE DEALLOCATION IN CleanUpDevice
     DEALLOCATE( F, STAT=AllocErr ) ! ADDED BY C.SALGADO ON 2017-03-22 TO DEALOCATE F WHILE REMOVING MT(NSpin, ntt) APPEARANCES.
     Write(*,*)"AllocErr = ",AllocErr
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Deallocation error for F(:,:,:)"
        STOP
     END IF
     Write(*,*) "DEALLOCATED F(:,:,:)"
!     if(.not.FMixing)then
     if(.not.FMixing)then
       !DEALLOCATE( D, DD, DDD, STAT=AllocErr ) ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
       DEALLOCATE( DD, DDD, STAT=AllocErr )
       Write(*,*)"AllocErr = ",AllocErr
       !IF( AllocErr /= 0 ) THEN
       IF(.NOT.(AllocErr == 0)) THEN
         !PRINT *, "ESF/Deallocation error for D(:,:,:), DD(:,:,:,:), DDD(:,:,:)" ! D NOT ALLOCATABLE AFTER 2017-03-22 BECAUSE OF INOUT CHARACTER.
         PRINT *, "ESF/Deallocation error for DD(:,:,:,:), DDD(:,:,:)"
         STOP
       END IF
       Write(*,*) "DEALLOCATED DD(:,:,:), DDD(:,:,:)"
     end if
     Write(*,*)"CALL CleanUpDevice"
     CALL CleanUpDevice
     Write(*,*)"RETURNING FROM ANT AFTER CleanUpDevice"
     RETURN
  END IF

  !*********
  ! BYE-BYE
  !*********
  PRINT *, "RETURNING FROM ANT BYE-BYE"
  RETURN
END SUBROUTINE ANT

SUBROUTINE readmat(fileid,t)

real, allocatable, intent(out) :: t(:)
integer, intent(in) :: fileid
integer :: N
!open(1,file='temp.dat')
!open(1,file=filename)
read(fileid,*) N   ! your first line with 22
allocate( t(N-1) )  ! further on you only have 21 elements
read(fileid,*)t          ! so, read them in 
print*, t
!deallocate(t)
close(1)
END SUBROUTINE readmat

END MODULE AntMod

