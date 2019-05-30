!*********************************************************!
!*********************  ANT.G09-2.4.1  *******************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   Juan Jose Palacios (1)                                !
!   David Jacob (2)                                       !
!   Angel J. Perez-Jimenez (3)                            !
!   Emilio SanFabian (3)                                  !
!                                                         !
!  (1) Departamento de Fisica de la Materia Condensada    !
!      Universidad Autonoma de Madrid                     !      
!      28049 Madrid (SPAIN)                               !
!  (2) Theory Department                                  !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!  (3) Departamento de Quimica Fisica                     !
!      Universidad de Alicante                            !
!      03690 Alicante (SPAIN)                             !
!                                                         !
!*********************************************************!
  MODULE G09Common
!**********************************************************
!  Common blocks for communication with Gaussian09        !
!**********************************************************
  USE parameters, ONLY: Nalpha, Nbeta
  USE preproc, ONLY: DEFMAXSHL => MaxShl, DEFMAXATM => MaxAtm
  IMPLICIT NONE
  
  PRIVATE

  ! ***********************
  ! Gaussian Common block B 
  ! ***********************
   Integer MaxShl,MaxPrm,MaxSh1,MaxS21,JAN,ShellA,ShellN,ShellT,ShellC,ShlADF,AOS,JAnSav,NShell,MaxTyp,I5DB1,I7FB1
      Real*8 EXX,C1,C2,C3,C4,X,Y,Z,RLam,RLamSv
      Parameter (MaxShl=DEFMAXSHL,MaxPrm=(3*MaxShl),MaxSh1=(MaxShl+1),MaxS21=(2*MaxShl+1))
      Common/B/EXX(MaxPrm),C1(MaxPrm),C2(MaxPrm),C3(MaxPrm),X(MaxShl),&
        Y(MaxShl),Z(MaxShl),JAN(MaxShl),ShellA(MaxShl),ShellN(MaxShl),&
        ShellT(MaxShl),ShellC(MaxShl),AOS(MaxShl),JAnSav(MaxShl),&
        RLam(MaxShl),RLamSv(MaxShl),NShell,MaxTyp,I5DB1,I7FB1
      Dimension C4(MaxShl),ShlADF(MaxShl)
      Equivalence (C4(1),C3(MaxSh1)),(ShlADF(1),C3(MaxS21))

  PUBLIC :: GetAtm4Sh, Get1stAO4Sh, GetNShell, GetShellT, GetShellC
  PUBLIC :: SetAtm4Sh, Set1stAO4Sh, SetNShell, SetShellT, SetShellC
  PUBLIC :: SetNAtoms, SetNAE, SetNBE, SetNE, SetNBasis, SetAN, SetAtmChg, SetAtmCo
  ! *************************
  ! Gaussian Common block Mol
  ! *************************

      Integer MaxAtm,NAtoms,ICharg,Multip,NAE,NBE,NE,NBasis,IAn,NBsUse,&
        IAtWgt,IAtTpR,IAtFrg,IAtRes,NPtMol,NPDMol,NumTpS,MolDum,IAtSpn,&
        IAtTpS,MaxTpS,MicOpt,IAtTyp
      Real*8 AtmChg,C,AtChMM,AtmWgt,AtZEff,AtQMom,AtGFac
      Parameter (MaxAtm=DEFMAXATM,MaxTpS=MaxAtm)
      Common /Mol/ NAtoms,ICharg,Multip,NAE,NBE,NE,NBasis,IAn(MaxAtm),&
        NBsUse,AtmChg(MaxAtm),C(3,MaxAtm),IAtTyp(MaxAtm),AtChMM(MaxAtm),&
        AtmWgt(MaxAtm),IAtWgt(MaxAtm),IAtTpR(MaxAtm),IAtFrg(MaxAtm),&
        IAtRes(MaxAtm),NPtMol,NPDMol,NumTpS,MolDum,IAtSpn(MaxAtm),&
        AtZEff(MaxAtm),AtQMom(MaxAtm),AtGFac(MaxAtm),IAtTpS(2,MaxTpS),&
        MicOpt(MaxAtm)
  
  PUBLIC :: GetNAtoms, GetNAE, GetNBE, GetNE, GetNBasis, GetAN, GetAtmChg, GetAtmCo

  CONTAINS

    ! **********************************
    ! Get Atom number for shell number i
    ! **********************************
    integer FUNCTION GetAtm4Sh( i )
      IMPLICIT NONE
      integer, INTENT(in) :: i
      GetAtm4Sh = JAN(i)
    END FUNCTION GetAtm4Sh
    ! **********************************
    ! Set Atom number for shell number i
    ! **********************************
    SUBROUTINE SetAtm4Sh( inputNShell, inputJAN )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNShell
      integer, DIMENSION(inputNShell), INTENT(in):: inputJAN
      integer :: i
      DO i=1,inputNShell
        JAN(i) = inputJAN(i)+1
        Write(*,'(A,I2,A,I2)')"JAN(",i,")=",JAN(i)
      END DO
    END SUBROUTINE SetAtm4Sh

    ! *********************************************
    ! Get number of first atomic orbital on shell i
    ! *********************************************
    integer FUNCTION Get1stAO4Sh( i )
      IMPLICIT NONE
      integer, INTENT(in) :: i
      Get1stAO4Sh = AOS( i )
    END FUNCTION Get1stAO4Sh
    ! *********************************************
    ! Set number of first atomic orbital on shell i
    ! *********************************************
    SUBROUTINE Set1stAO4Sh( inputNShell, inputAOS )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNShell
      integer, DIMENSION(inputNShell), INTENT(in):: inputAOS
      integer :: i
      DO i=1,inputNShell
        AOS( i ) = inputAOS(i)+1 ! ADD +1 BECAUSE FORTRAN DOES NOT USE ELEMENT 0 IN ARRAY.
        Write(*,'(A,I2,A,I2)')"AOS(",i,")=",AOS( i )
      END DO
    END SUBROUTINE Set1stAO4Sh

    ! ********************************
    ! Get Number of shells in molecule
    ! ********************************
    integer FUNCTION GetNShell()
      IMPLICIT NONE
      GetNShell = NShell
    END FUNCTION GetNShell
    ! ********************************
    ! Set Number of shells in molecule
    ! ********************************
    SUBROUTINE SetNShell( inputNShell )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNShell
      NShell = inputNShell
    END SUBROUTINE SetNShell

    ! **************************
    ! Get shell type for shell i
    ! **************************
    ! s = 0, p = 1, d = 2 ,f = 3 ...
    ! sp = 1 also but shell constraint is different
    integer FUNCTION GetShellT( i )
      IMPLICIT NONE
      integer, INTENT(in) :: i
      GetShellT = ShellT( i )
    END FUNCTION GetShellT
    ! **************************
    ! Set shell type for shell i
    ! **************************
    ! s = 0, p = 1, d = 2 ,f = 3 ...
    ! sp = 1 also but shell constraint is different
    SUBROUTINE SetShellT( inputNShell, inputShellT )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNShell
      integer, DIMENSION(inputNShell), INTENT(in):: inputShellT
      integer :: i
      DO i=1,inputNShell
        ShellT( i ) = inputShellT( i )
      END DO
    END SUBROUTINE SetShellT

    ! ********************************
    ! Get shell constraint for shell i
    ! ********************************
    ! to distinguish between p and sp-shell
    ! sp = 2 , p = 1  
    integer FUNCTION GetShellC( i )
      IMPLICIT NONE
      integer, INTENT(in) :: i
      GetShellC = ShellC( i )
    END FUNCTION GetShellC
    ! ********************************
    ! Set shell constraint for shell i
    ! ********************************
    ! to distinguish between p and sp-shell
    ! sp = 2 , p = 1  
    SUBROUTINE SetShellC( inputNShell, inputShellC )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNShell
      integer, DIMENSION(inputNShell), INTENT(in):: inputShellC
      integer :: i
      DO i=1,inputNShell
        ShellC( i ) = inputShellC( i )
      END DO
    END SUBROUTINE SetShellC

    ! ***************************
    ! Number of atoms in molecule
    ! ***************************
    integer FUNCTION GetNAtoms()
      IMPLICIT NONE
      GetNAtoms = NAtoms
    END FUNCTION GetNAtoms
    ! ***************************
    ! Number of atoms in molecule
    ! ***************************
    SUBROUTINE SetNAtoms( inputNAtoms )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNAtoms
      NAtoms = inputNAtoms
    END SUBROUTINE SetNAtoms

    ! *****************************
    ! Get Number of Alpha electrons
    ! *****************************
    integer FUNCTION GetNAE()
      IMPLICIT NONE
      if (Nalpha  < 0) then
      GetNAE = NAE
      else
      GetNAE = Nalpha
      end if
    END FUNCTION GetNAE
    ! *****************************
    ! Set Number of Alpha electrons
    ! *****************************
    SUBROUTINE SetNAE( inputNAE )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNAE
      NAE = inputNAE
      Nalpha = -1 ! Important. Not used when negative. Use NAE always.
    END SUBROUTINE SetNAE
    
    ! ****************************
    ! Get Number of Beta electrons
    ! ****************************
    integer FUNCTION GetNBE()
      IMPLICIT NONE
      if (Nbeta  < 0) then
      GetNBE = NBE
      else
      GetNBE = Nbeta
      end if
    END FUNCTION GetNBE
    ! ****************************
    ! Set Number of Beta electrons
    ! ****************************
    SUBROUTINE SetNBE( inputNBE )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNBE
      NBE = inputNBE
      Nbeta = -1 ! Important. Not used when negative. Use NAE always.
    END SUBROUTINE SetNBE

    ! ***********************
    ! Get Number of electrons
    ! ***********************
    integer FUNCTION GetNE()
      IMPLICIT NONE
      if (Nbeta  >= 0 .and. Nalpha >= 0) then
      GetNE = Nalpha+Nbeta
      else
      GetNE = NE
      end if
    END FUNCTION GetNE
    ! ***********************
    ! Set Number of electrons
    ! ***********************
    SUBROUTINE SetNE( inputNE )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNE
      NE = inputNE
      Nalpha = -1
      Nbeta = -1
    END SUBROUTINE SetNE

    ! *****************************
    ! Get Number of basis functions
    ! *****************************
    integer FUNCTION GetNBasis()
      IMPLICIT NONE
      GetNBasis = NBasis
    END FUNCTION GetNBasis
    ! *****************************
    ! Set Number of basis functions
    ! *****************************
    SUBROUTINE SetNBasis( inputNBasis )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNBasis
      NBasis = inputNBasis
    END SUBROUTINE SetNBasis

    ! ****************************
    ! Get Atomic number of atom ia
    ! ****************************
    integer FUNCTION GetAN( ia )
      IMPLICIT NONE
      integer, INTENT(in) :: ia
      GetAN = IAN(ia)
    END FUNCTION GetAN
    ! ****************************
    ! Set Atomic number of atom ia
    ! ****************************
    SUBROUTINE SetAN( inputNAtoms, inputIAN )
      IMPLICIT NONE
      integer, INTENT(in) :: inputNAtoms
      integer, DIMENSION(inputNAtoms), INTENT(in):: inputIAN
      integer :: i
      DO i=1,inputNAtoms
        IAN(i) = inputIAN(i)
      END DO
    END SUBROUTINE SetAN

    ! ********************************************
    ! Get Atom core charge = nuke + core electrons
    ! ********************************************
    real FUNCTION GetAtmChg( ia )
      IMPLICIT NONE 
      integer, INTENT(in) :: ia
      GetAtmChg = AtmChg(ia)
    END FUNCTION GetAtmChg
    ! ********************************************
    ! Set Atom core charge = nuke + core electrons
    ! ********************************************
    SUBROUTINE SetAtmChg( inputNAtoms, inputAtmChg )
      IMPLICIT NONE 
      integer, INTENT(in) :: inputNAtoms
      real, DIMENSION(inputNAtoms), INTENT(in):: inputAtmChg
      integer :: i
      DO i=1,inputNAtoms
        AtmChg(i) = inputAtmChg(i)
      END DO
    END SUBROUTINE SetAtmChg

    ! ***********************
    ! Get Coordinates of atom
    ! ***********************
    real FUNCTION GetAtmCo( j, ia )
      IMPLICIT NONE
      integer, INTENT(in) :: j, ia
      GetAtmCo = C(j,ia)
    END FUNCTION GetAtmCo
    ! ***********************
    ! Set Coordinates of atom
    ! ***********************
    SUBROUTINE SetAtmCo( inputNAtoms, inputAtmCo )
!      use constants, only: au2ang
      !use iso_c_binding
      IMPLICIT NONE
      integer, INTENT(in) :: inputNAtoms
      real, DIMENSION(3,inputNAtoms), INTENT(in) :: inputAtmCo
      !real(c_long_double), INTENT(in) :: inputAtmCo(3,inputNAtoms)
      integer :: ia, j
      real :: au2ang = 0.52918d0
      write(*,'(A)')"ENTER SetAtmCo( inputNAtoms, inputAtmCo )"
      DO ia=1,inputNAtoms
        DO j=1,3
          C(j,ia) = inputAtmCo(j,ia)/au2ang
        END DO
        write(*,'(A,I2,A,3F8.4)')"C(1:3,",ia,") = ",( C(j,ia) ,j=1,3)
      END DO
      write(*,'(A)')"EXIT SetAtmCo( inputNAtoms, inputAtmCo )"
    END SUBROUTINE SetAtmCo

  END MODULE G09Common
