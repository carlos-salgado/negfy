!*********************************************************!
!*********************  ANT.G-2.4.1  *********************!
!*********************************************************!
!                                                         !
!  Copyright (c) by                                       !
!                                                         !
!  Juan Jose Palacios (1)                                 !
!  David Jacob (2)                                        !
!  Maria Soriano (1)                                      !
!  Angel J. Perez-Jimenez (3)                             !
!  Emilio SanFabian (3)                                   !
!  Jose Antonio Antonio Verges (4)                        !
!  Enrique Louis (5)                                      !
!                                                         !
! (1) Departamento de Fisica de la Materia Condensada     !
!     Universidad Autonoma de Madrid                      !    
!     28049 Madrid (SPAIN)                                !
! (2) Theory Department                                   !
!     Max-Planck-Institute for Microstructure Physics     !
!     Halle, 06120 (GERMANY)                              !
! (3) Departamento de Quimica Fisica                      !
!     Universidad de Alicante                             !
!     03690 Alicante (SPAIN)                              !
! (4) Insto. de Ciencias de Materiales de Madrid (ICMM)   !
!     Consejo Superior de Investigacion y Ciencia (CSIC)  !
!     28049 Madrid (SPAIN)                                !
! (5) Departamento de Fisica Aplicada                     !
!     Universidad de Alicante                             !    
!     03690 Alicante (SPAIN)                              !
!                                                         !
!*********************************************************!
  MODULE device
!*********************************************************!
!  Main module for computation of device Green's function !
!  and transport                                          !
!*********************************************************!
  use ANTCommon
  use parameters, only: Debug, DebugDev
  !use util, only: PrintCMatrix
  implicit none
  save

  private

  public :: DevNAOrbs, DevNSpin, DevShift, DevFockMat, DevDensMat, SetDevDensMat, SetDevFockMat
  public :: DevHWFockMat, SetDevHWFockMat, DevDGibbsYMat, SetDevDGibbsYMat
  public :: DevDGibbsYKernel1Mat, SetDevDGibbsYKernel1Mat, DevDGibbsYKernel2Mat, SetDevDGibbsYKernel2Mat
  public :: BuildLiouvillian, CompCGibbsY
  public :: LeadsOn, SwitchOnLeads, SecantOn, SwitchOnSecant, SwitchOffSecant
  public :: EvaluationOn, SwitchOnEvaluation
  public :: SwitchOnChargeCntr, SwitchOffChargeCntr, SwitchOffSpinLock, SwitchOnSpinLock
  public :: InitDevice, ReadDensMat, ReadFockMat, CleanUpDevice, InitElectrodes, Transport
  !--------------------------------
  public :: IntDDOStimesE, DDOStimesE0, WorkFock, WorkEnergy
  !--------------------------------
  !*****************************
  ! Module's internal variables
  !*****************************

  ! *** Number of atomic orbitals, number of non-degenerate spin-channels ***
  integer :: NAOrbs, NSpin, DNAOrbs

  ! *** Total number of electrons in central device region
  integer :: NCDEl, NCDAO1, NCDAO2

  ! *** Actual electron charge for certain Fermi energy ***
  real :: QAlpha, QBeta, Q_SOC

  ! *** Overlap matrix S of device ***
  real, dimension(:,:),allocatable :: SD, InvSD
  real, dimension(:,:), allocatable :: S_SOC

  ! *** Complex S^+1/2 matrix ***
  complex*16, dimension(:,:),allocatable :: SPH, SNH
  complex*16, dimension(:,:,:),allocatable :: PDOUT
  complex*16, dimension(:,:,:),allocatable :: PDOUTGIBBS
  complex*16, dimension(:,:,:),allocatable :: HWOUT
  complex, dimension(:,:,:),allocatable :: LiouvSOp
  complex, dimension(:,:,:),allocatable :: CGibbsY, CGibbsYKernel1, CGibbsYKernel2
  complex*16, dimension(:,:),allocatable :: PDOUT_SOC
  
  ! *** Hamiltonian and density matrix of device
  real, dimension(:,:,:),allocatable :: HD
  real, dimension(:,:,:),allocatable :: PD
  real, dimension(:,:,:),allocatable :: PDGIBBS
  real, dimension(:,:,:),allocatable :: HW
  real, dimension(:,:,:),allocatable :: DGibbsY, DGibbsYKernel1, DGibbsYKernel2
  real, dimension(:,:),allocatable :: PD_SOC
  complex*16, dimension(:,:), allocatable :: H_SOC

  ! *** Orthognalization matrix for device ***
  real, dimension(:,:),allocatable :: OD

  ! *** Energy shifts ***
  real :: shift, ShiftUp, ShiftDown

  ! *** spin excess charge ***
  real :: exspin

  ! *** internal spin variable 1=up,2=down ***
  integer  :: ispin
  
  ! *** Lowest and highest eigen value of Fock matrix ***
  real :: LEV,HEV
 
  ! *** lower and upper band edge of electrode ***
  real :: EMinEc, EMaxEc

  ! *** lower and upper energy bound ***
  real :: EMin, EMax

  ! *** Density of states projected on atoms at the Fermi energy
  real, dimension(:,:), allocatable :: AtomDOSEF

  ! *** Control Switches ***
  logical :: Leads       = .false.
  logical :: Secant      = .false.
  logical :: Evaluation  = .false.
  logical :: UDTrans     = .false.
  logical :: ChargeCntr  = .false.
  logical :: SpinLock    = .true.

  ! Whether device Hamiltonain has been orthogonalized
  logical :: HDOrtho = .false.
  
  ! *** In which atom to calculate the spinmu ***
  integer :: spinatom
  ! *** Number of electrons projected on atoms at the Fermi energy
  ! *** Actual electron charge for certain Fermi energy ***
  real :: LocalQAlpha, LocalQBeta
  real, dimension(:), allocatable :: alphaelec, betaelec, alphalocalshift, betalocalshift
  real, dimension(:,:,:),allocatable :: LocPD
  complex*16, dimension(:,:,:),allocatable :: LocPDOUT

!!$OMP THREADPRIVATE(shift)
  contains

  !********************************
  !*** Public module procedures ***
  !********************************

  !**************************************
  ! Access functions to private entities
  !**************************************


  ! *** Total number of atomic orbitals in device Hilbert space ***
  integer function DevNAOrbs()
    implicit none
    DevNAOrbs=NAOrbs
  end function DevNAOrbs

  ! *** Number of non-degenerate spin bands ***
  integer function DevNSpin() 
    implicit none
    DevNSpin=NSpin
  end function DevNSpin

  ! *** Spin excess charge to equilibrate Fermi level up and down ***
  real function DevShift()
    implicit none
    DevShift=-shift
  end function DevShift

  ! *** Get matrix Element of Density matrix ***
  real function DevFockMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevFockMat = HD(is, i, j)
  end function DevFockMat

  ! *** Set matrix Element of Fock matrix ***
  subroutine SetDevFockMat( is, i, j, fij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: fij
    HD(is, i, j) = fij
  end subroutine SetDevFockMat

  ! *** Get matrix Element of Density matrix ***
  real function DevDensMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDensMat = PD(is, i, j)
  end function DevDensMat

  ! *** Set matrix Element of Density matrix ***
  subroutine SetDevDensMat( is, i, j, pij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: pij
    PD(is, i, j) = pij
  end subroutine SetDevDensMat

  ! *** Get matrix Element of HWFock matrix ***
  real function DevHWFockMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevHWFockMat = HW(is, i, j)
  end function DevHWFockMat

  ! *** Set matrix Element of HWFock matrix ***
  subroutine SetDevHWFockMat( is, i, j, hwij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: hwij
    HW(is, i, j) = hwij
  end subroutine SetDevHWFockMat

  ! *** Get matrix Element of CGibbsY matrix ***
  real function DevDGibbsYMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDGibbsYMat = DGibbsY(is, i, j)
  end function DevDGibbsYMat

  ! *** Set matrix Element of CGibbsY matrix ***
  subroutine SetDevDGibbsYMat( is, i, j, DGibbsYij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: DGibbsYij
    DGibbsY(is, i, j) = DGibbsYij
  end subroutine SetDevDGibbsYMat

  ! *** Get matrix Element of CGibbsYKernel1 matrix ***
  real function DevDGibbsYKernel1Mat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDGibbsYKernel1Mat = DGibbsYKernel1(is, i, j)
  end function DevDGibbsYKernel1Mat

  ! *** Set matrix Element of CGibbsYKernel1 matrix ***
  subroutine SetDevDGibbsYKernel1Mat( is, i, j, DGibbsYKernel1ij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: DGibbsYKernel1ij
    DGibbsYKernel1(is, i, j) = DGibbsYKernel1ij
  end subroutine SetDevDGibbsYKernel1Mat

  ! *** Get matrix Element of CGibbsYKernel2 matrix ***
  real function DevDGibbsYKernel2Mat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDGibbsYKernel2Mat = DGibbsYKernel2(is, i, j)
  end function DevDGibbsYKernel2Mat

  ! *** Set matrix Element of CGibbsYKernel2 matrix ***
  subroutine SetDevDGibbsYKernel2Mat( is, i, j, DGibbsYKernel2ij )
    implicit none
    integer, intent(in) :: is, i, j
    real, intent(in) :: DGibbsYKernel2ij
    DGibbsYKernel2(is, i, j) = DGibbsYKernel2ij
  end subroutine SetDevDGibbsYKernel2Mat

  ! ***
  logical function LeadsOn()
    implicit none
    LeadsOn = Leads
  end function LeadsOn

  ! *** 
  subroutine SwitchOnLeads
    implicit none
    Leads = .true.
  end subroutine SwitchOnLeads

  ! ***
  logical function SecantOn()
    implicit none
    SecantOn = Secant
  end function SecantOn

  ! ***
  subroutine SwitchOnSecant()
    implicit none
    Secant = .true.
  end subroutine SwitchOnSecant

  ! ***
  subroutine SwitchOffSecant()
    implicit none
    Secant = .false.
  end subroutine SwitchOffSecant

  ! ***
  logical function EvaluationOn()
    implicit none
    EvaluationOn = Evaluation
  end function EvaluationOn

  ! ***
  subroutine SwitchOnEvaluation()
    implicit none
    Evaluation = .true.
  end subroutine SwitchOnEvaluation

  subroutine SwitchOnChargeCntr()
    implicit none
    ChargeCntr = .true.
  end subroutine SwitchOnChargeCntr

  ! ***
  subroutine SwitchOffChargeCntr()
    implicit none
    ChargeCntr = .false.
  end subroutine SwitchOffChargeCntr

  ! ***
  subroutine SwitchOffSpinLock()
    implicit none
    SpinLock = .false.
  end subroutine SwitchOffSpinLock

  ! ***
  subroutine SwitchOnSpinLock()
    implicit none
    SpinLock = .true.
  end subroutine SwitchOnSpinLock
  

  !***********************************************
  !* Initialize device for transport calculation *
  !***********************************************
  subroutine InitDevice( NBasis, UHF, S )
    use constants, only: d_zero
    use numeric, only: RMatPow
    use parameters, only: ElType, FermiStart, Overlap, HybFunc, CompFock!, boolComputeCGibbsY exists not yet.
    use cluster, only: AnalyseCluster, AnalyseClusterElectrodeOne, AnalyseClusterElectrodeTwo, NAOAtom, NEmbedBL
    use g09Common, only: GetNAtoms, GetAtmChg
    use correlation
    use orthogonalization
    use ANTCommon
    
    implicit none

    integer, intent(in) ::  NBasis!, NSpin
    logical, intent(in) :: UHF
    real, dimension(NBasis,NBasis),intent(in) :: S
!    real, dimension(NSpin,NBasis,NBasis),intent(in) :: F

    integer :: AllocErr, ios, iatom, NEmbed1, NEmbed2

    real, dimension(NBasis,NBasis) :: RSPH 

    write(ifu_log,*) "-------------------"
    write(ifu_log,*) "Initializing device"
    write(ifu_log,*) "-------------------"

    NAOrbs = NBasis
    DNAOrbs = 2*NBasis
    NSpin=1 ! COMMENTED BY C.SALGADO ON 2017-03-25 BECAUSE NSpin IS NOW INPUT.
    if( UHF ) NSpin=2 ! COMMENTED BY C.SALGADO ON 2017-03-25 BECAUSE NSpin IS NOW INPUT.
    exspin=d_zero

    ! Dynamic arrays 
    allocate( SD(NAOrbs,NAOrbs), InvSD(NAOrbs,NAOrbs), &
         HD(NSpin,NAOrbs,NAOrbs), &
         PD(NSpin,NAOrbs,NAOrbs), &
         PDOUT(NSpin,NAOrbs,NAOrbs), &
         STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Allocation error for SD, InvSD, SMH, SPH, H, P"
       stop
    end if
    if(CompFock)then
      allocate( HW(NSpin,NAOrbs,NAOrbs), &
         HWOUT(NSpin,NAOrbs,NAOrbs), &
         STAT=AllocErr )
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Allocation error for HW, HWOUT"
        stop
      end if
      allocate( CGibbsY(NSpin,NAOrbs,NAOrbs), &
                DGibbsY(NSpin,NAOrbs,NAOrbs), &
                CGibbsYKernel1(NSpin,NAOrbs,NAOrbs), &
                DGibbsYKernel1(NSpin,NAOrbs,NAOrbs), &
                CGibbsYKernel2(NSpin,NAOrbs,NAOrbs), &
                DGibbsYKernel2(NSpin,NAOrbs,NAOrbs), &
                STAT=AllocErr )
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Allocation error for CGibbsY"
        stop
      end if
      allocate( PDGIBBS(NSpin,NAOrbs,NAOrbs), &
                PDOUTGIBBS(NSpin,NAOrbs,NAOrbs), &
                STAT=AllocErr )
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Allocation error for PDGIBBS"
        stop
      end if
    end if

    SD = S
    call RMatPow( SD, -1.0d0, InvSD )

    !HD = F ! ADDED BY C.SALGADO ON 2017-03-25 TO INITIALIZE HAMILTONIAN WITHOUT SELF-ENERGY.

    if( HybFunc ) call InitCorrelation(NAOrbs,NSpin)

    allocate( SPH(NAorbs,NAOrbs), STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/InitDevice: Allocation error for SPH(:,:)"
       stop
    end if
    ! Computing transformation matrix S^+1/2
    call RMatPow( SD,  0.5d0, RSPH )
    SPH = RSPH

    shift = -FermiStart
    shiftup = shift
    shiftdown = shift

    IF( ElType(1) == "BETHE" .and. ElType(2) == "BETHE" ) THEN 
      call AnalyseCluster
    ELSE IF  (ElType(1) == "BETHE" .and. ElType(2) == "GHOST" ) THEN
      call AnalyseClusterElectrodeOne
    ELSE IF  (ElType(1) == "GHOST" .and. ElType(2) == "BETHE" ) THEN
      call AnalyseClusterElectrodeTwo
    ELSE IF  (ElType(1) == "GHOST" .and. ElType(2) == "GHOST" ) THEN
      continue                           
    ELSE 
      print *, 'These electrodes are not implemented yet !!!'
      stop
    END IF

    call InitElectrodes

    EMin = 0.0d0
    EMax = 0.0d0

    LEV = 0.0d0
    HEV = 0.0d0

    ! Compute number of electrons 
    ! in central device region
    

    IF (Overlap < 0.01) THEN
       NEmbed1=0
       NEmbed2=0
    ELSE
       NEmbed1=NEmbedBL(1)
       NEmbed2=NEmbedBL(2)
    END IF

   
    print *, "---------------------------------------------------"
    print *, " Details on device and contacting atoms -----------"
    print *, "---------------------------------------------------"

    print *, "NEmbed(1) =", NEmbedBL(1)
    print *, "NEmbed(2) =", NEmbedBL(2)

    NCDEl = 0
    do iatom=NEmbed1+1,GetNAtoms()-NEmbed2
      NCDEl = NCDEl + GetAtmChg(iatom)
    end do
    
    print *, "Number of electrons in neutral reduced device"
    print *, "NCDEl = ", NCDEl
    
    print *, "First and last orbital in reduced device"
    IF (Overlap < 0.01) THEN
       NCDAO1 = 1
    ELSE
       NCDAO1 = 0
    END IF
    do iatom=1,NEmbed1
       NCDAO1 = NCDAO1 + NAOAtom(iatom) 
    end do

    print *, "NCDAO1 = ", NCDAO1

    NCDAO2 = NCDAO1-1
    do iatom = NEmbed1+1,GetNAtoms()-NEmbed2
       NCDAO2 = NCDAO2 + NAOAtom(iatom) 
    end do

    IF  (ElType(1) == "GHOST" .and. ElType(2) == "GHOST" ) NCDAO2 = NAOrbs
    print *, "NCDAO2 = ", NCDAO2
    print *, "---------------------------------------------------"

  end subroutine InitDevice
  
  !***********************************************
  !* Read initial density matrix from file P.dat *
  !***********************************************
  subroutine ReadDensMat(densitymatrix)
    use parameters, only: NSpinEdit, SpinEdit, MRStart, SpinDel, PFix, NFix, IFix, densitymatrixx
    use constants, only: d_zero
    use numeric, only: RMatPow
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
    use g09Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
    use ANTCommon
    implicit none
    
    integer :: norb, ni, nj, isp, n, iatom, jatom, is, i, j, ii, jj, AOStart, AO_BEG, AO_END, iorb, jorb, iato, jato
    real :: density, TrP, xxx
    real, dimension(NSpin,NAOrbs,NAOrbs) :: PDMod
    real, dimension(NAOrbs,NAOrbs) :: PDAux
    character (len=50) :: densitymatrix
  
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"
    write(ifu_log,*) "Starting calculation with density matrix from file  ", densitymatrix
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"

    PD=0.0d0
    PDMod=0.0d0

    open(ifu_dm,file=densitymatrix,status='old')
    read(ifu_dm,*,end=21) shift
    shift = -shift
    do i=1,2*NAOrbs*NAOrbs
       read(ifu_dm,*,end=21)is,ii,jj,density
       PDMod(is,ii,jj)=density
       PDMod(is,jj,ii)=density
    end do
 21 close (ifu_dm)

    if (PFix) then
       write(ifu_log,*) "-------------------------------------------------------------------------------------------"
       write(ifu_log,*) "  ... and using supplementary density matrix from file  ", densitymatrixx
       write(ifu_log,*) "-------------------------------------------------------------------------------------------"
       if (NFix == 0) print *,'Warning ... NFix = 0'
       open(ifu_dmx,file=densitymatrixx,status='old')
       read(ifu_dmx,*) xxx   
       do norb=1,8*NAOrbs*NAOrbs
          read(ifu_dmx,*,end=22)is,ni,nj,density,iato,iorb,jato,jorb
          do isp=1,NSpin
             i=0
             do iAtom=1,GetNAtoms()
                do n=1,NFix
                   if (iAtom == IFix(n)) then 
                      i=i+NAOAtom(iAtom)
                      goto 11
                   end if
                end do
                do ii=1,NAOAtom(iAtom)
                   i=i+1
                   j=0
                   do jAtom=1,GetNAtoms()
                      do n=1,NFix
                         if (jAtom == IFix(n)) then 
                            j=j+NAOAtom(jAtom)
                            goto 12
                         end if
                      end do
                      do jj=1,NAOAtom(jAtom)
                         j=j+1
                         if (is == isp .and. iato == iAtom .and. jato == jAtom .and. iorb == ii .and. jorb == jj) then
                            PDMod(isp,i,j)=density
                            PDMod(isp,j,i)=density
                         end if
                      end do
12                 end do
                end do
11           end do
          end do
       end do
22     close(ifu_dmx)
       
    end if
  
    !open(111,file='readdensitymatrix',status='unknown')
    !do is=1,NSpin
    !i=0
    !do iAtom=1,GetNAtoms()
    ! do ii=1,NAOAtom(iAtom)
    !    i=i+1
    !    j=0
    !    do jAtom=1,GetNAtoms()
    !       do jj=1,NAOAtom(jAtom)
    !          j=j+1
    !          write(111,'(i3,4i5,e18.10)')is,iAtom,ii,jAtom,jj,PDMod(is,i,j)
    !end do
    !end do
    !end do
    !end do
    !end do
    !close(111)

    !
    ! Manipulate atomic spins of initial guess 
    ! if SpinEdit or MRStart set
    !
    if( NSpinEdit > 0 .or. MRStart > 0 .or. SpinDel )then

       if( MRStart > 0 )then
          do iAtom=MRStart,GetNAtoms()
             SpinEdit( iAtom ) = -1
          end do
       end if

       if( SpinDel )then
          do iAtom=1,GetNAtoms()
             SpinEdit( iAtom ) = 0
          end do
       end if

       PD = PDMod
       PDMod = d_zero

       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             if( SpinEdit(iAtom) ==  1 .and. SpinEdit(jAtom) ==  1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)
                      PDMod(1,i,j) = PD(1,i,j)
                      PDMod(2,i,j) = PD(2,i,j)
                   end do
                end do
             else if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      PDMod(2,i,j) = PD(1,i,j)
                      PDMod(1,i,j) = PD(2,i,j)
                   end do
                end do
             else
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      PDMod(1,i,j) = 0.5d0*(PD(1,i,j)+PD(2,i,j))            
                      PDMod(2,i,j) = 0.5d0*(PD(1,i,j)+PD(2,i,j))            
                   end do
                end do
             end if
          end do
       end do

    end if

    ! Normalize to correct number of electrons when recommended
     
    TrP = d_zero
    do is=1,NSpin
       PDAux=MATMUL(PDMod(is,:,:),SD)
       do i=1,NAOrbs
          TrP = TrP + PDAux(i,i)
       end do
    end do
    if (NSpin ==1) then
       PRINT *, "Tr[P*S] of initial guess =", TrP*2.0
       if( NSpinEdit > 0 .or. MRStart > 0 .or. PFIX) then
         PD = PDMod*GetNE()/(TrP*2.0)
       else
         PD = PDMod
       end if
    else
       PRINT *, "Tr[P*S] of initial guess =", TrP
       if( NSpinEdit > 0 .or. MRStart > 0 .or. PFIX) then
          PD = PDMod*GetNE()/(TrP)
       else
          PD = PDMod
       end if
    end if
    PRINT *, "--------------------------------"
       
  end subroutine ReadDensMat

  !********************************************
  !* Read initial Fock matrix from file F.dat *
  !********************************************
  subroutine ReadFockMat(fockmatrix)
    use parameters, only: NSpinEdit, SpinEdit, MRStart, SpinDel !!!, PFix, NFix, IFix, densitymatrixx
    use constants, only: d_zero
    use numeric, only: RMatPow
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
    use g09Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
    use ANTCommon
    implicit none
    
    integer :: norb, ni, nj, isp, n, iatom, jatom, is, i, j, ii, jj, AOStart, AO_BEG, AO_END, iorb, jorb, iato, jato
    real :: fock !, TrP, xxx

    real, dimension(NSpin,NAOrbs,NAOrbs) :: HDMod
    character (len=50) :: fockmatrix
  
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"
    write(ifu_log,*) "Starting calculation with fock matrix from file  ", fockmatrix
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"

    HD=0.0d0
    HDMod=0.0d0

    open(ifu_fm,file=fockmatrix,status='old')
    read(ifu_fm,*,end=21) shift
    shift = -shift
    do i=1,2*NAOrbs*NAOrbs
       read(ifu_fm,*,end=21)is,ii,jj,fock
       HD(is,ii,jj)=fock
       HD(is,jj,ii)=fock
    end do
    Write(ifu_log,'(A)')"HAMILTONIAN READ FROM F.name.dat is: "
    do iSpin=1,NSpin
      do i=1,NAOrbs
        write(ifu_log, '(10F8.4)')(HD(iSpin,i,j) ,j=1,NAOrbs)
      end do
    end do
    Write(ifu_log,'(A)')"END DISPLAY HAMILTONIAN"
21 close (ifu_fm)

    !
    ! Manipulate atomic spins of initial guess 
    ! if SpinEdit or MRStart set
    !
    if( NSpinEdit > 0 .or. MRStart > 0 .or. SpinDel )then

       if( MRStart > 0 )then
          do iAtom=MRStart,GetNAtoms()
             SpinEdit( iAtom ) = -1
          end do
       end if

       if( SpinDel )then
          do iAtom=1,GetNAtoms()
             SpinEdit( iAtom ) = 0
          end do
       end if

       HDMod = d_zero

       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             if( SpinEdit(iAtom) ==  1 .and. SpinEdit(jAtom) ==  1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)
                      HDMod(1,i,j) = HD(1,i,j)
                      HDMod(2,i,j) = HD(2,i,j)
                   end do
                end do
             else if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      HDMod(2,i,j) = HD(1,i,j)
                      HDMod(1,i,j) = HD(2,i,j)
                   end do
                end do
             else
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      HDMod(1,i,j) = 0.5d0*(HD(1,i,j)+HD(2,i,j))            
                      HDMod(2,i,j) = 0.5d0*(HD(1,i,j)+HD(2,i,j))            
                   end do
                end do
             end if
          end do
       end do

       HD = HDMod

    end if

  end subroutine ReadFockMat


  !**************************
  ! Deallocate dynamic arrays
  !**************************
  subroutine CleanUpDevice
    use BetheLattice, only: CleanUpBL, LeadBL
    use parameters, only: ElType, BiasVoltage, CompFock
    implicit none
    integer :: AllocErr, LeadNo
    Write(*,*)"ENTERED CleanUpDevice"
    print *, "DEALLOCATING SD"
    deallocate( SD, STAT=AllocErr )
    Write(*,*)"AllocErr = ",AllocErr
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Deallocation error for SD"
       stop
    end if
    print *, "CORRECT DEALLOCATION OF SD"
!    deallocate( InvSD, STAT=AllocErr )
!    Write(*,*)"AllocErr = ",AllocErr
!    if( AllocErr /= 0 ) then
!       print *, "DEVICE/Deallocation error for InvSD"
!       stop
!    end if
!    print *, "CORRECT DEALLOCATION OF InvSD"
!    print *, "DEALLOCATING HD"
!    deallocate( HD, STAT=AllocErr )
!    Write(*,*)"AllocErr = ",AllocErr
!    if( AllocErr /= 0 ) then
!       print *, "DEVICE/Deallocation error for HD"
!       stop
!    end if
!    print *, "CORRECT DEALLOCATION OF HD"
    print *, "DEALLOCATING PD"
    deallocate( PD, STAT=AllocErr )
    Write(*,*)"AllocErr = ",AllocErr
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Deallocation error for PD"
       stop
    end if
    print *, "CORRECT DEALLOCATION OF PD"
    print *, "DEALLOCATING PDOUT"
    deallocate( PDOUT, STAT=AllocErr )
    Write(*,*)"AllocErr = ",AllocErr
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Deallocation error for PDOUT"
       stop
    end if
    print *, "CORRECT DEALLOCATION OF PDOUT"
    print *, "DEALLOCATING SPH"
    deallocate( SPH, STAT=AllocErr )
    Write(*,*)"AllocErr = ",AllocErr
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Deallocation error for SPH"
       stop
    end if
    print *, "CORRECT DEALLOCATION OF SPH"

    if(CompFock)then
      print *, "DEALLOCATING HW"
      deallocate( HW, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for HW"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF HW"
      print *, "DEALLOCATING HWOUT"
      deallocate( HWOUT, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for HWOUT"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF HWOUT"

      print *, "DEALLOCATING DGibbsY"
      deallocate( DGibbsY, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for DGibbsY"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF DGibbsY"
      print *, "DEALLOCATING CGibbsY"
      deallocate( CGibbsY, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for CGibbsY"
        stop
      end if

      print *, "DEALLOCATING DGibbsYKernel1"
      deallocate( DGibbsYKernel1, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for DGibbsYKernel1"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF DGibbsYKernel1"
      print *, "DEALLOCATING CGibbsYKernel1"
      deallocate( CGibbsYKernel1, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for CGibbsYKernel1"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF CGibbsYKernel1"

      print *, "DEALLOCATING DGibbsYKernel2"
      deallocate( DGibbsYKernel2, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for DGibbsYKernel2"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF DGibbsYKernel2"
      print *, "DEALLOCATING CGibbsYKernel2"
      deallocate( CGibbsYKernel2, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for CGibbsYKernel2"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF CGibbsYKernel2"

      print *, "DEALLOCATING PDGIBBS"
      deallocate( PDGIBBS, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for PDGIBBS"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF PDOUTGIBBS"
      print *, "DEALLOCATING PDOUTGIBBS"
      deallocate( PDOUTGIBBS, STAT=AllocErr )
      Write(*,*)"AllocErr = ",AllocErr
      if( AllocErr /= 0 ) then
        print *, "DEVICE/Deallocation error for PDOUTGIBBS"
        stop
      end if
      print *, "CORRECT DEALLOCATION OF PDOUTGIBBS"
    end if

!    IF(.NOT.(BiasVoltage==0.0d0))THEN
!      deallocate( PDOUT, STAT=AllocErr )
!      if( AllocErr /= 0 ) then
!        print *, "DEVICE/Deallocation error for PDOUT"
!        stop
!      end if
!    END IF
!    print *, "CORRECT DEALLOCATION OF PDOUT"
    do LeadNo=1,2
       select case( ElType(LeadNo) )
       case( "BETHE" )
          call CleanUpBL( LeadBL(LeadNo) ) 
       end select
    end do
    print *, "LEAVING CleanUpDevice"
  end subroutine CleanUpDevice


  !*************************
  !* Initialize electrodes *
  !*************************
  subroutine InitElectrodes
    use BetheLattice, only: InitBetheLattice, LeadBL, BL_EMin, BL_EMax
    use OneDLead, only: Init1DLead, Lead1d, L1D_EMin, L1D_EMax
    use parameters, only: ElType
    implicit none 
    integer :: LeadNo
    real, dimension(2) :: EMin, EMax

    do LeadNo=1,2
       select case( ElType(LeadNo) )
       case( "BETHE" )
          call InitBetheLattice( LeadBL(LeadNo), LeadNo )
          EMin(LeadNo) = BL_EMin( LeadBL(LeadNo) )
          EMax(LeadNo) = BL_EMax( LeadBL(LeadNo) )
       case( "1DLEAD" )
          call Init1DLead( Lead1d(LeadNo), LeadNo )
          EMin(LeadNo) = L1D_EMin( Lead1d(LeadNo) )
          EMax(LeadNo) = L1D_EMax( Lead1d(LeadNo) )
       case( "GHOST" )
          EMin(LeadNo) = -100.0                       
          EMax(LeadNo) =  100.0                       
       case DEFAULT
          print *, "Unknown option for ElType:", ElType(LeadNo)
          stop
       end select
    end do
    
    EMinEc = min( Emin(1),EMin(2) )
    EMaxEc = max( EMax(1),EMax(2) )
  end subroutine InitElectrodes

  
  !***************************
  !* Solve transport problem *
  !***************************
  subroutine Transport(F,ADDP) 
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, &
                          IntEnergy, BiasEnergy, DiagFock, SPINMU, BiasVoltage, CompFock, CompGibbsY
    use numeric, only: RMatPow, RSDiag
    use correlation
    use orthogonalization
    use constants, only: Hart
    implicit none

    logical,intent(out) :: ADDP
    real, dimension(NSpin,NAOrbs,NAOrbs),intent(in) :: F

    real, dimension(NAOrbs) :: evals
    real, dimension(:,:),allocatable :: SPM

    real :: diff !!,TrP,QD
    integer :: i,j,is, info, AllocErr
    character(80) strfmt

    HD = F

!    write(strfmt,'(a,i0,a)') '(A,I1,A,',size(HD,2),'(1x,F12.6))'
!    do iSpin=1,NSpin
!      do i=1,NAOrbs
!        write(*, '(A,I1,A,I2,A,10F10.5)')"HD(",iSpin,",",i,",:) = ",(HD(iSpin,i,j) ,j=1,NAOrbs)
!      end do
!      !write(ifu_log, strfmt)"F(",iSpin,")",HD(iSpin,:,:)/Hart
!    end do

    !
    ! Estimate upper bound for maximal eigenvalue of HD and use it for upper and lower energy boundaries
    !
    if( NSpin == 1 ) EMax = maxval(sum(abs(HD(1,:,:)),1))
    if( NSpin == 2 ) EMax = max(maxval(sum(abs(HD(1,:,:)),1)),maxval(sum(abs(HD(2,:,:)),1)))
    EMin = -EMax
    print *, "EMin=", EMin
    print *, "EMax=", EMax

    if( .not. DMImag .and. ChargeCntr )then
       print *, "--------------"
       print *, "Charge Control"
       print *, "--------------"
       ! Find upper and lower energy bound 
       call FindEnergyBounds
    endif

    if( DFTU ) call Add_DFT_plus_U_Pot( PD, HD )

    ! 1st FALSE MEANS NOT COMPUTE FOCK WITH QHWTot. 2nd FALSE MEANS NOT COMPUTE CGibbsY.
    if(.not.DMImag) call CompDensMat(ADDP,.false.,.false.)
    if(DMImag) call CompDensMat2(ADDP)
! Next if(Evaluation) piece was added by C. Salgado to integrate E*DOS between Emin and Efermi.
    if( Evaluation )then
      print *
      print *, "****************************************** "
      print *, "*                                        * "
      print *, "*        ANT.G09 final analysis          * "
      print *, "*                                        * "
      print *, "****************************************** "
      print *
      if( SPINMU )then
        !call Hamiltonian
        !call LDOS
        !call MullPop
        if(.not.DMImag) call CompLocalMu(ADDP)
        if(DMImag) call CompLocalMu(ADDP)
      end if
      !if( IntEnergy ) call WorkEnergy
      if( DiagFock ) call WorkFock
      if( IntEnergy ) call WorkEnergy
      if( BiasEnergy ) call WorkBiasEnergy
      if( CompFock )then
      !if( CompGibbsY )then
        if(BiasVoltage==0.0d0)then
          ! 1st TRUE MEANS COMPUTE FOCK WITH QHWTot. 2nd FALSE MEANS NOT COMPUTE CGibbsY.
          !call CompDensMat(ADDP,.true.,.false.)
          call BuildLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
          call CompDensMat(ADDP,.true.,.true.)
          call DeallocateLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
        elseif(BiasVoltage/=0.0d0)then
          call BuildLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
          ! 1st TRUE MEANS COMPUTE FOCK WITH QHWTot. 2nd TRUE MEANS COMPUTE CGibbsY.
          call CompDensMat(ADDP,.true.,.true.)
          call DeallocateLiouvillian ! COMMENTED ON 2018-04-24 BECAUSE LIOUVILIAN IS NOT USED ANYMORE.
        end if
      end if
!    end if
!    if( Evaluation )then
!       print *
!       print *, "****************************************** "
!       print *, "*                                        * "
!       print *, "*        ANT.G09 final analysis          * "
!       print *, "*                                        * "
!       print *, "****************************************** "
!       print *

       IF( ANT1DInp ) call WriteANT1DInput

       if( POrtho )then
          allocate( OD(NAorbs,NAOrbs), STAT=AllocErr )
          if( AllocErr /= 0 ) then
             print *, "DEVICE/InitDevice: Allocation error for OD(:,:)"
             stop
          end if
          do ispin=1,NSpin
             PD(ispin,:,:) = matmul( SD, matmul(PD(ispin,:,:), SD ) )
          end do
          HDOrtho = .true.
          call ProjOrtho(cix, SD, OD )
          ! Othogonalize density matrix density matrix and Hamiltonian
          do ispin=1,NSpin
             HD(ispin,:,:) = matmul( transpose(OD), matmul(F(ispin,:,:), OD) )
             PD(ispin,:,:) = matmul( OD, matmul(PD(ispin,:,:), transpose(OD) ) )
          end do
          HDOrtho = .true.
       end if
       if( DiagCorrbl ) call DiagCorrBlocks( HD, SD )
       call Hamiltonian
       IF ( HybFunc ) call CompHybFunc
       IF ((ElType(1) == "GHOST" .or. ElType(2) == "GHOST") .and. LDOS_Beg <= LDOS_End) CALL LDOS
       IF (ElType(1) /= "GHOST" .and. ElType(2) /= "GHOST") THEN            
          IF( RedTransmE >= RedTransmB  ) call EigenChannelAnalysis
          call transmission
       END IF
       call MullPop
       return
    end if
    !ADDP=.FALSE. ! ADDED BY C.SALGADO TO SKIP PROBLEMS ON 2017-04-03
  end subroutine Transport

  subroutine BuildLiouvillian
!    use g09Common, only: GetNAE, GetNBE
    use parameters, only: eta
    use constants, only: d_zero, d_one, c_zero, c_one, ui, eleccharge, hbar
    use numeric, only: CInv
    implicit none

!    complex, dimension(:,:,:),allocatable :: LiouvSOp!, LiouvSOpL, LiouvSOpR
    complex, dimension(NAOrbs,NAOrbs) :: IDD
    complex, dimension(NAOrbs*NAOrbs,NAOrbs*NAOrbs) :: LiouvSOpPiv
!    real ::
!    real ::
    integer :: iSpin, i, j, k, l
    integer :: info, allocerr
    integer :: n, ipiv(NAOrbs*NAOrbs)
    complex*16, DIMENSION( 4*NAOrbs*NAOrbs ) :: work

    if(DebugDev)then
    Write(*,'(A)')"***********************************"
    Write(*,'(A)')"**** ENTER BuildLiouvillian!!! ****"
    Write(*,'(A)')"***********************************"
    end if

    allocate( LiouvSOp(NSpin,NAorbs*NAorbs,NAorbs*NAorbs),&
              STAT=AllocErr )
    if( AllocErr /= 0 ) then
      print *, "DEVICE/BuildLiouvillian: Allocation error for LiouvSOp(:,:,:)"
      stop
    end if
    LiouvSOpPiv = c_zero
    LiouvSOp = c_zero

    IDD = c_zero                           ! Initialize the array.
    forall(j = 1:NAOrbs) IDD(j,j) = c_one     ! Set the diagonal.

    Write(*,'(A)')"DON'T FORGET TO MULTIPLY BY THE ELECTRON CHARGE!!!"
    do iSpin=1,NSpin
      LiouvSOpPiv = c_zero
      do i=1,NAOrbs
        do j=1,NAOrbs
          do k=1,NAOrbs
            do l=1,NAOrbs
              !LiouvSOpL((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)=HD(i,j)*IDD(k,l)
              !LiouvSOpR((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)=IDD(i,j)*HD(k,l)
              ! FOR EACH SUPEROP ELEMENT, 1st SUMANDO BELONGS TO LiouvSOpL & 2nd TO LiouvSOpR, NEGATIVE.
              ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.
              !LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)= -(COMPLEX(HD(iSpin,i,j),0.0D0)*IDD(k,l) - IDD(i,j)*COMPLEX(HD(iSpin,k,l),0.0D0))
              LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l)= - HD(iSpin,i,j)*IDD(k,l) + IDD(i,j)*HD(iSpin,k,l) ! STRICTLY -L.
              ! ADD COMPLEX ETA TO THE DIAGONAL.
              if((i==j) .and. (k==l))then
                LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l) = LiouvSOpPiv((i-1)*NAOrbs+k,(j-1)*NAOrbs+l) + ui*eta ! BECAUSE -L+i*eta
              end if
              !Write(*,'(I4,I4,I4,F12.8,F12.8)')iSpin,i,j,LiouvSOp(iSpin,i,j)
            end do
          end do
        end do
      end do

      if(DebugDev)then
      if(iSpin==1)Write(*,'(A)')"ALPHA LIOUVILLIAN BEFORE INVERSION"
      if(iSpin==2)Write(*,'(A)')"BETA LIOUVILLIAN BEFORE INVERSION"
      call PrintCMatrix(LiouvSOpPiv)
      end if

      !info = Cinv(LiouvSOpPiv)
      n = SIZE( LiouvSOpPiv, 1)
      CALL zgetrf(n,n,LiouvSOpPiv,n,ipiv,info)
      CALL zgetri(n,LiouvSOpPiv,n,ipiv,work,4*n,info)
      if( info /= 0 ) THEN
        WRITE(ifu_log,*)'Device/BuildLiouvillian using CInv (zgetrf,zgetri) in device.f90'
        WRITE(ifu_log,*)'INFO=',info
!       STOP
      end if

      do i=1,NAOrbs*NAOrbs
        do j=1,NAOrbs*NAOrbs
          LiouvSOp(iSpin,i,j)=LiouvSOpPiv(i,j)
        end do
      end do

    end do

    if(DebugDev)then
    do iSpin=1,NSpin
      if(iSpin==1)Write(*,'(A)')"ALPHA INVERSE LIOUVILLIAN"
      if(iSpin==2)Write(*,'(A)')"BETA INVERSE LIOUVILLIAN"
      !do i=1,NAOrbs*NAOrbs
      !  do j=1,NAOrbs*NAOrbs
      !    Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)')iSpin,i,j,"(",DREAL(LiouvSOp(iSpin,i,j)),") + i*(",DIMAG(LiouvSOp(iSpin,i,j)),")"
      !  end do
      !end do
      if(iSpin==1)Write(*,'(A)')"ALPHA LIOUVILLIAN AFTER INVERSION"
      if(iSpin==2)Write(*,'(A)')"BETA LIOUVILLIAN AFTER INVERSION"
      call PrintCMatrix(LiouvSOpPiv)
    end do
    end if

    if(DebugDev)then
    Write(*,'(A)')"***********************************"
    Write(*,'(A)')"**** EXIT BuildLiouvillian!!! *****"
    Write(*,'(A)')"***********************************"
    end if

  end subroutine BuildLiouvillian

  subroutine DeallocateLiouvillian
    implicit none
    integer :: allocerr
    !deallocate(LiouvSOp,stat=allocerr)
    if (allocerr /= 0 ) then
      print*,"Problems deallocating LiouvSOp"
      !stop
    end if
  end subroutine DeallocateLiouvillian

  subroutine WorkEnergy
    use g09Common, only: GetNAE, GetNBE
    use constants, only: d_zero, d_pi
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, DiagFock, IntEnergy, QEXCESS
    implicit none
    ! *** Eigenvalues of the Fock-Matrix ***
    real, dimension(:,:,:),allocatable :: OHD
    real :: upIntegerDOSE, downIntegerDOSE, IntegerDOSE, E1
    real :: hartreeupIntegerDOSE, hartreedownIntegerDOSE, hartreeIntegerDOSE,evperhartree
    integer :: i,j,is, info, AllocErr, NROT
    
    evperhartree = 2.721138D1

        !if(DMImag)then
        if(NSpin==2 .and. SPINLOCK )then
        !if(NSpin==2)then
          print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
          do ispin=1,NSpin
            !if(shiftup/=shiftdown)then
              if (ispin.eq.1)then
                E1=shiftup
                write(ifu_log,*)'ispin',ispin
                write(ifu_log,*)'- SHIFTUP',-shiftup
              end if
              if (ispin.eq.2)then
                E1=shiftdown
                write(ifu_log,*)'ispin',ispin
                write(ifu_log,*)'- SHIFTDOWN',-shiftdown
              end if
              if (ispin.eq.1)upIntegerDOSE = FullEnergy(d_zero)
              if (ispin.eq.2)downIntegerDOSE = FullEnergy(d_zero)
            !end if
            write(ifu_log,*)'--------------------------------------------------------'
            if (ispin.eq.1) then
              write(ifu_log,'(A,F9.5)') ' Fermi energy for alpha electrons= ', -shiftup
              write(ifu_log,*)
              write(ifu_log,'(A,F16.8)') ' Energy of the alpha electrons: ', upIntegerDOSE
            end if
            if (ispin.eq.2) then
              write(ifu_log,'(A,F9.5)') ' Fermi energy for beta electrons=  ', -shiftdown
              write(ifu_log,*)
              write(ifu_log,'(A,F16.8)') ' Energy of the beta electrons: ', downIntegerDOSE
            end if
            write(ifu_log,*)'--------------------------------------------------------'
          end do
        else if(NSpin==2 .and. .not. SPINLOCK )then
          print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
          ispin = 1
            E1=shift
            write(ifu_log,*)'ispin',ispin
            write(ifu_log,*)'- SHIFT',-shift
            IntegerDOSE = FullEnergy(d_zero)
            write(ifu_log,*)'--------------------------------------------------------'
            write(ifu_log,'(A,F9.5)') ' Fermi energy for alpha/beta electrons= ', -shift
            write(ifu_log,*)
            write(ifu_log,'(A,F16.8)') ' Energy of the alpha/beta electrons: ', IntegerDOSE
            write(ifu_log,*)'--------------------------------------------------------'

 !if(.not.DMImag) then
        !if(NSpin==1)then
        else
          print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
          do ispin=1,NSpin
            E1=shift
            write(ifu_log,*)'ispin',ispin
            write(ifu_log,*)'- SHIFT',-shift
            IntegerDOSE = FullEnergy(d_zero)
            write(ifu_log,*)'--------------------------------------------------------'
            write(ifu_log,'(A,F9.5)') ' Fermi energy for alpha/beta electrons= ', -shift
            write(ifu_log,*)
            write(ifu_log,'(A,F16.8)') ' Energy of the alpha/beta electrons: ', IntegerDOSE
            write(ifu_log,*)'--------------------------------------------------------'
          end do
        end if
        if(NSpin==2 .and. SPINLOCK )then
          IntegerDOSE = upIntegerDOSE + downIntegerDOSE
          hartreeIntegerDOSE = IntegerDOSE/evperhartree
        else if(NSpin==2 .and. .not. SPINLOCK )then
          IntegerDOSE = 2.0d0*IntegerDOSE
          hartreeIntegerDOSE = IntegerDOSE/evperhartree
        else if(NSpin==1)then
          IntegerDOSE = 2.0d0*IntegerDOSE
          hartreeIntegerDOSE = IntegerDOSE/evperhartree
        end if
        write(ifu_log,*)
        write(ifu_log,*)'--------------------------------------------------------'
        !write(ifu_log,'(A,F10.5)') ' Total energy of electrons:  ', IntegerDOSE
        write(ifu_log,'(A,F16.8)') ' Total energy of electrons (eV):  ', IntegerDOSE
        write(ifu_log,'(A,F10.5)') ' Total energy of electrons (Hartree):  ', hartreeIntegerDOSE
        write(ifu_log,*)'--------------------------------------------------------'
  end subroutine WorkEnergy

  subroutine WorkBiasEnergy
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, DiagFock, BiasEnergy
    implicit none
    ! *** Eigenvalues of the Fock-Matrix ***
    real, dimension(:,:,:),allocatable :: OHD
    real :: upIntegerDOSE, downIntegerDOSE, IntegerDOSE
    real :: hartreeupIntegerDOSE, hartreedownIntegerDOSE, hartreeIntegerDOSE,evperhartree
    integer :: i,j,is, info, AllocErr, NROT
    
    evperhartree = 2.721138D1

 !if(.not.DMImag) then
        if(NSpin==1)then
            write(ifu_log,*)''
            write(ifu_log,*)''
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE DOS(E)*E'
            write(ifu_log,*)'- SHIFT',-shift
            !write(ifu_log,*)'SHIFT',shift
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'' 

            IntegerDOSE = QXTotEnergy(shift) ! THE - SIGN DOESN'T WORK.
            !IntegerDOSE = CompEnergy(shift) ! THE BEST ONE. INTEGRATES TO MU.
            !IntegerDOSE = CompEnergy(-shift) ! THE BEST ONE. INTEGRATES TO MU.

            hartreeIntegerDOSE = IntegerDOSE/evperhartree
            !print *, "Integral of DOS*E =", upIntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (eV) =',IntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (Hartree) =',hartreeIntegerDOSE
            write(ifu_log,*)''
        end if
        !if(DMImag)then
        if(NSpin==2)then
          write(ifu_log,*)'- SHIFTUP',-shiftup
          write(ifu_log,*)'- SHIFTDOWN',-shiftdown

          if(shiftup/=shiftdown)then
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE ALPHA DOS(E)*E'
            write(ifu_log,*)'- SHIFTUP',-shiftup
            !write(ifu_log,*)'SHIFTUP',shiftup
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE BETA DOS(E)*E'
            write(ifu_log,*)'- SHIFTDOWN',-shiftdown
            !write(ifu_log,*)'SHIFTDOWN',shiftdown
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)''

            upIntegerDOSE = QXTotEnergy(shiftup) ! THE - SIGN DOESN'T WORK.
            !upIntegerDOSE = CompEnergy(shiftup) ! THE BEST ONE. INTEGRATES TO MU.
            !upIntegerDOSE = CompEnergy(-shiftup) ! THE BEST ONE. INTEGRATES TO MU.

            hartreeupIntegerDOSE = upIntegerDOSE/evperhartree
            !print *, "Integral of spin-up DOS*E =", upIntegerDOSE
     write(ifu_log,*)'Integral of spin-up DOS*E (eV) =',upIntegerDOSE
     write(ifu_log,*)'Integral of spin-up DOS*E (Hartree) =',hartreeupIntegerDOSE
            write(ifu_log,*)''

            downIntegerDOSE = QXTotEnergy(shiftdown) ! THE - SIGN DOESN'T WORK.
            !downIntegerDOSE = CompEnergy(shiftdown) ! THE BEST ONE. INTEGRATES TO MU.
            !downIntegerDOSE = CompEnergy(-shiftdown) ! THE BEST ONE. INTEGRATES TO MU.

     hartreedownIntegerDOSE = downIntegerDOSE/evperhartree
            !print *, "Integral of spin-down DOS*E =", downIntegerDOSE
     write(ifu_log,*)'Integral of spin-down DOS*E (eV) =',downIntegerDOSE
     write(ifu_log,*)'Integral of spin-down DOS*E (Hartree) =',hartreedownIntegerDOSE
            write(ifu_log,*)''

            IntegerDOSE = upIntegerDOSE + downIntegerDOSE
            hartreeIntegerDOSE = IntegerDOSE/evperhartree
          else
            write(ifu_log,*)''
            write(ifu_log,*)''
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'ENERGY UP TO WHICH WE INTEGRATE DOS(E)*E'
            write(ifu_log,*)'- SHIFT',-shift
            !write(ifu_log,*)'SHIFT',shift
            write(ifu_log,*)'*************************************************************'
            write(ifu_log,*)'' 

            IntegerDOSE = QXTotEnergy(shift) ! THE - SIGN DOESN'T WORK.
            !IntegerDOSE = CompEnergy(shift) ! THE BEST ONE. INTEGRATES TO MU.
            !IntegerDOSE = CompEnergy(-shift) ! THE BEST ONE. INTEGRATES TO MU.

            hartreeIntegerDOSE = IntegerDOSE/evperhartree
            !print *, "Integral of DOS*E =", upIntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (eV) =',IntegerDOSE
     write(ifu_log,*)'Integral of DOS*E (Hartree) =',hartreeIntegerDOSE
            write(ifu_log,*)''
          end if
        end if

        write(ifu_log,*)'*************************************************************'
        write(ifu_log,*)''
 write(ifu_log,*)'Integral of total DOS*E (eV) =',IntegerDOSE
 write(ifu_log,*)'Integral of total DOS*E (Hartree) =',hartreeIntegerDOSE 
        write(ifu_log,*)''
        write(ifu_log,*)'*************************************************************'

  end subroutine WorkBiasEnergy

  subroutine WorkFock
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, DiagFock, IntEnergy, glue
    use numeric, only: RMatPow, RSDiag, CDiag, Jacobi, eigsrt, ordks, balanc, RInv, GAUSSJ
    use correlation
    use util
    use orthogonalization
    !use lapack_blas, only: zgetri, zgetrf
    implicit none
    ! *** Eigenvalues of the Fock-Matrix ***
    real, dimension(:,:,:),allocatable :: OHD
    complex*16, dimension(:,:,:),allocatable :: COHD
    real, dimension(:,:),allocatable :: kseigenvalues
    real, dimension(:),allocatable :: upkseigenvalues, DupOHD
    real, dimension(:),allocatable :: downkseigenvalues, DdownOHD
    real, dimension(:),allocatable :: allkseigenvalues, ordkseigenvalues, hartreekseigenvalues
    real, dimension(:,:),allocatable :: upHD, upOHD, invupOHD, pivupOHD, VupOHD
    real, dimension(:,:),allocatable :: downHD, downOHD, invdownOHD, pivdownOHD, VdownOHD
    complex*16, dimension(:),allocatable :: DupCOHD
    complex*16, dimension(:,:),allocatable :: upCOHD, VupCOHD
    complex*16, dimension(:),allocatable :: DdownCOHD
    complex*16, dimension(:,:),allocatable :: downCOHD, VdownCOHD
    real :: kohnshamenergy, upkohnshamenergy, downkohnshamenergy
    real :: hartreeks, hartreeupks, hartreedownks, evperhartree
    integer :: i,j,is, info, AllocErr, NROT
    integer :: upocckscount, downocckscount, norbks
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr
    complex*16 :: cshiftup, cshiftdown


! Dynamic arrays 
    allocate( OHD(NSpin,NAOrbs,NAOrbs), &
         kseigenvalues(NAOrbs,NAOrbs), &
         allkseigenvalues(2*NAOrbs), &
         hartreekseigenvalues(2*NAOrbs), &
         ordkseigenvalues(2*NAOrbs), &
         upkseigenvalues(NAOrbs), &
         downkseigenvalues(NAOrbs), &
         upHD(NAOrbs,NAOrbs), &
         downHD(NAOrbs,NAOrbs), &
         upOHD(NAOrbs,NAOrbs), &
         downOHD(NAOrbs,NAOrbs), &
         invupOHD(NAOrbs,NAOrbs), &
         invdownOHD(NAOrbs,NAOrbs), &
         pivupOHD(NAOrbs,NAOrbs), &
         pivdownOHD(NAOrbs,NAOrbs), &
         DupOHD(NAOrbs), &
         DdownOHD(NAOrbs), &
         VupOHD(NAOrbs,NAOrbs), &
         VdownOHD(NAOrbs,NAOrbs), &
         COHD(NSpin,NAOrbs,NAOrbs), &
         upCOHD(NAOrbs,NAOrbs), &
         downCOHD(NAOrbs,NAOrbs), &
         DupCOHD(NAOrbs), &
         DdownCOHD(NAOrbs), &
         VupCOHD(NAOrbs,NAOrbs), &
         VdownCOHD(NAOrbs,NAOrbs), &
         STAT=AllocErr )
    
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Allocation error for SD, InvSD, SMH, SPH, H, P"
       stop
    end if
    
    evperhartree = 2.721138D1   
    upkohnshamenergy = 0.d0
    downkohnshamenergy = 0.d0
    !write(ifu_log,*)''
    !write(ifu_log,*)'******************************************************************'
    !write(ifu_log,*)'Before orthogonalization and diagonalization (looks yet symmetric)'
      !NAOrbs = SIZE(HD(1,:,:),2)
      !NBasis = SIZE(HD(1,:,:),2)
10111 format(a6,i4)
10112 format(f12.6)
!10111 format(a6,i4,f12.6)

        !write(ifu_log,*)'SNH'
        !call PrintCMatrix( SNH(1:4,1:4) )
        !write(ifu_log,*)'SPH'
        !call PrintCMatrix( SPH(1:4,1:4) )

        do ispin=1,NSpin
            !write(ifu_log,*)'HD'
            !call PrintRMatrix( HD(ispin,1:4,1:4) )
        end do

        do ispin=1,NSpin
           !OHD(ispin,:,:) = HD(ispin,:,:)
           ! Looks like this is the correct way.
           !OHD(ispin,:,:) = matmul( SNH, matmul( HD(ispin,:,:), SPH ) )
           ! Orthogonalization by SNH*HD*SNH
           OHD(ispin,:,:) = matmul( SNH, matmul( HD(ispin,:,:), SNH ) ) ! This seems to be the correct expression.
           !OHD(ispin,:,:) = HD(ispin,:,:)
           !OHD(ispin,:,:) = matmul( SNH, matmul( F(ispin,:,:), SPH ) )

           !OHD(ispin,:,:) = matmul( transpose(SD), matmul(HD(ispin,:,:), SD) )
           !OHD(ispin,:,:) = matmul( transpose(OD), matmul(F(ispin,:,:), OD) )
           !PD(ispin,:,:) = matmul( OD, matmul(PD(ispin,:,:), transpose(OD) ) )
        end do

        if(.not.DMIMag)then
          cshiftup = dcmplx(shift)
        else
          cshiftup = dcmplx(shiftup)
          cshiftdown = dcmplx(shiftdown)
        end if
        do ispin=1,NSpin
            !if ((ispin == 1)) upOHD=OHD(ispin,:,:)+sigl+sigr
            if ((ispin == 1))then
              call CompSelfEnergies( ispin, cshiftup, sigl, sigr )
              sigr = glue*sigr
              sigl = glue*sigl
              upCOHD=dcmplx(OHD(ispin,:,:))+sigl+sigr
            end if
            !if(DMImag)
            !if ((ispin == 2)) downOHD=OHD(ispin,:,:)+sigl+sigr
            if ((ispin == 2))then
              call CompSelfEnergies( ispin, cshiftdown, sigl, sigr )
              sigr = glue*sigr
              sigl = glue*sigl
              downCOHD=dcmplx(OHD(ispin,:,:))+sigl+sigr
            end if
        end do
    !write(ifu_log,*)'***************************************************************'
    !write(ifu_log,*)''
    !	write(ifu_log,*)'After orthogonalization'
    !write(ifu_log,*)'***************************************************************'
    !write(ifu_log,*)''

    do ispin=1,NSpin
      if ((ispin == 1))then
        !write(ifu_log,*)'Spin-up Hamiltonian'
        !call PrintRMatrix( upOHD(1:4,1:4) )
        !call PrintCMatrix( upCOHD(1:4,1:4) )
        !write(ifu_log,*)'***************************************************************'
      end if
      !if(DMImag)
      if ((ispin == 2))then
        !write(ifu_log,*)'Spin-down Hamiltonian'
        !call PrintRMatrix( downOHD(1:4,1:4) )
        !call PrintCMatrix( downCOHD(1:4,1:4) )
        !write(ifu_log,*)'***************************************************************'
      end if
    end do
    do ispin=1,NSpin
      if ((ispin == 1)) then
        !call balanc(upOHD)
        !call zgetrf(NAOrbs,NAOrbs,upCOHD,NAOrbs,ipiv,info)
        !write(ifu_log,*)'After balancing'
        !write(ifu_log,*)'***************************************************************'
        !write(ifu_log,*)''
        !write(ifu_log,*)'Spin-up Hamiltonian'
        !call PrintRMatrix( upOHD(1:4,1:4) )
        !call PrintCMatrix( upCOHD(1:4,1:4) )
      end if
      !if(DMImag)
      if ((Nspin == 2)) then
        !call balanc(downOHD)
        !call zgetrf(NAOrbs,NAOrbs,downCOHD,NAOrbs,ipiv,info)
        !write(ifu_log,*)'After balancing'
        !write(ifu_log,*)'***************************************************************'
        !write(ifu_log,*)''
        !write(ifu_log,*)'Spin-down Hamiltonian'
        !call PrintRMatrix( downOHD(1:4,1:4) )
        !call PrintCMatrix( downCOHD(1:4,1:4) )
      end if
    end do

    !write(ifu_log,*)'***************************************************************'
    !write(ifu_log,*)''

    !write(ifu_log,*)'*************************************************************'
    !write(ifu_log,*)''
    !write(ifu_log,*)'Before diagonalization'
    do ispin=1,NSpin
      if ((ispin == 1)) then
        call CDiag(upCOHD,DupCOHD,info)
        ! Sorting eigenvalues and eigenvectors in ascending order of eigenvalues.
        upOHD = real(upCOHD,8)
        DupOHD = real(DupCOHD,8)
        DupOHD = -DupOHD
        call eigsrt(DupOHD,VupOHD)
        DupOHD = -DupOHD
        !print *, ' DupOHD(j)'
        !print *, ( DupOHD(j), j=1,NAOrbs )
      end if
      !if(DMImag)
      if ((ispin == 2)) then
        call CDiag(downCOHD,DdownCOHD,info)
        ! Sorting eigenvalues and eigenvectors in ascending order of eigenvalues.
        downOHD = real(downCOHD,8)
        DdownOHD = real(DdownCOHD,8)
        DdownOHD = -DdownOHD
        call eigsrt(DdownOHD,VdownOHD)
        DdownOHD = -DdownOHD
        !print *, ' DdownOHD(j)'
        !print *, ( DdownOHD(j), j=1,NAOrbs )
      end if
    end do
    
    !if(.not.DMImag)
    if ((Nspin == 1)) DdownOHD = DupOHD !because beta and alpha eigenvalues are equal, overwrite the betas with the alphas. YOU DON'T NEED TO MULTIPLY LATER BY TWO.
    if ((Nspin == 1)) norbks = nint(2*QAlpha)
    if ((Nspin == 2)) norbks = nint(QAlpha+QBeta)

    do i=1,NAOrbs
      do j=i,NAOrbs
        do ispin=1,NSpin
          if ((ispin == 1)) upOHD(j,i)=upOHD(i,j)
          if ((ispin == 2)) downOHD(j,i)=downOHD(i,j)
        end do
      end do
    end do

    !write(ifu_log,*)'************************************************************'
    !write(ifu_log,*)''
    !write(ifu_log,*)'After diagonalization'
        
    allkseigenvalues(1:NAOrbs) = DupOHD(:)
    ordkseigenvalues(1:NAOrbs) = 1.0
    allkseigenvalues(NAOrbs+1:2*NAOrbs) = DdownOHD(:)
    ordkseigenvalues(NAOrbs+1:2*NAOrbs) = 2.0

    allkseigenvalues = -allkseigenvalues
    call ordks(allkseigenvalues,ordkseigenvalues)
    allkseigenvalues = -allkseigenvalues

    !do i=1,NAOrbs+1 !I first wrote this, don't remember why.

!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S occ. eigen. (eV)',(allkseigenvalues(i),i=1,norbks)
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S virt. eigen. (eV)',(allkseigenvalues(i),i=norbks+1,2*NAOrbs)
!    write(ifu_log,*)''
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
!    !do i=1,2*NAOrbs
!    !   shiftallkseigenvalues(i) = allkseigenvalues(i) - shift
!    !end do
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)'SHIFT = ', -SHIFT
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S occ. eigen. (eV) (shifted)',(allkseigenvalues(i)-shift,i=1,norbks)
!    write(ifu_log,*)''
!    write(ifu_log,*)'K-S virt. eigen. (eV) (shifted)',(allkseigenvalues(i)-shift,i=norbks+1,2*NAOrbs)
!    write(ifu_log,*)''
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)''
    do i=1,2*NAOrbs
      hartreekseigenvalues(i) = allkseigenvalues(i)/27.21184
      !shifthartreekseigenvalues(i) = shiftallkseigenvalues(i)/27.21184
    end do

    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
    write(ifu_log,*)''
    write(ifu_log,*)'K-S occ. eigen. (Hartree)',(hartreekseigenvalues(i),i=1,norbks)
    write(ifu_log,*)''
    !write(ifu_log,*)'K-S virt. eigen. (Hartree)',(hartreekseigenvalues(i),i=norbks+1,2*NAOrbs)
    !write(ifu_log,*)''
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)'SHIFT = ', -SHIFT
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Total num of K-S occ. orbitals',norbks
    write(ifu_log,*)''
    write(ifu_log,*)'K-S occ. eigen. (Hartree) (shifted)',(hartreekseigenvalues(i)-shift/evperhartree,i=1,norbks)
    write(ifu_log,*)''
    !write(ifu_log,*)'K-S virt. eigen. (Hartree) (shifted)',(hartreekseigenvalues(i)-shift/evperhartree,i=norbks+1,2*NAOrbs)
    !write(ifu_log,*)''
    write(ifu_log,*)'************************************************************'
    write(ifu_log,*)''

    upocckscount = 0
    downocckscount = 0
    do i=1,norbks
      if(ordkseigenvalues(i).eq.1.0) upocckscount = upocckscount + 1
      if(ordkseigenvalues(i).eq.2.0) downocckscount = downocckscount + 1
      IF (upocckscount + downocckscount .GE. norbks) EXIT
    end do

      !do i=1,NAOrbs
      !  !if(DupOHD(i).le.0.0d0)then
      !  if(i.le.NAOrbs/2)then
      !    upkohnshamenergy=upkohnshamenergy+DupOHD(i)
      !  end if
      !  !if(DdownOHD(i).le.0.0d0)then
      !  if(i.le.NAOrbs/2)then
      !    downkohnshamenergy=downkohnshamenergy+DdownOHD(i)
      !  end if
      !end do

    !do i=1,2*NAOrbs
    !  shiftDupOHD(i) = DupOHD(i)
    !  shiftDdownOHD(i) = DdownOHD(i)
    !end do

    upkohnshamenergy = sum(DupOHD(1:upocckscount))
    downkohnshamenergy = sum(DdownOHD(1:downocckscount))

    !shiftupkohnshamenergy = upkohnshamenergy - upocckscount*shift
    !shiftdownkohnshamenergy = downkohnshamenergy - downocckscount*shift

    !kohnshamenergy = upkohnshamenergy + downkohnshamenergy
    kohnshamenergy = sum(allkseigenvalues(1:norbks))
    hartreeupks = upkohnshamenergy/evperhartree
    hartreedownks = downkohnshamenergy/evperhartree
    hartreeks = kohnshamenergy/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied alpha eigenvalues',upocckscount
!    !print '(1000(ES14.4))', ( DupOHD(j), j=1,upocckscount )
!    print *, ( DupOHD(j), j=1,upocckscount )
!    write(ifu_log,*)'up Kohn-Sham energy sum',upkohnshamenergy
!    write(ifu_log,*)'upKS in Hartrees',hartreeupks
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied beta eigenvalues',downocckscount
!    !print '(1000(ES14.4))', ( DdownOHD(j), j=1,downocckscount )
!    print *, ( DdownOHD(j), j=1,downocckscount )
!    write(ifu_log,*)'down Kohn-Sham energy sum',downkohnshamenergy
!    write(ifu_log,*)'downKS in Hartrees',hartreedownks
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
!    write(ifu_log,*)'Total Kohn-Sham energy sum',kohnshamenergy
!    write(ifu_log,*)'KS in Hartrees',hartreeks
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!
!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)'SHIFT = ', -SHIFT
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied alpha eigenvalues (shifted)',upocckscount
!    !print '(1000(ES14.4))', ( DupOHD(j)-shift, j=1,upocckscount )
!    print *, ( DupOHD(j)-shift, j=1,upocckscount )
!    write(ifu_log,*)'up Kohn-Sham energy sum (shifted)',upkohnshamenergy-upocckscount*shift
!    write(ifu_log,*)'upKS in Hartrees (shifted)',hartreeupks-upocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied beta eigenvalues (shifted)',downocckscount
!    !print '(1000(ES14.4))', ( DdownOHD(j)-shift, j=1,downocckscount )
!    print *, ( DdownOHD(j)-shift, j=1,downocckscount )
!    write(ifu_log,*)'down Kohn-Sham energy sum (shifted)',downkohnshamenergy-downocckscount*shift
!    write(ifu_log,*)'downKS in Hartrees (shifted)',hartreedownks-downocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
!    write(ifu_log,*)'Total Kohn-Sham energy sum (shifted)',kohnshamenergy-norbks*shift
!    write(ifu_log,*)'KS in Hartrees (shifted)',hartreeks-norbks*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''

    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Occupied alpha eigenvalues in Hartrees',upocckscount
    !print '(1000(ES14.4))', ( DupOHD(j), j=1,upocckscount )
    print *, ( DupOHD(j)/evperhartree, j=1,upocckscount )
    write(ifu_log,*)'up Kohn-Sham energy sum',upkohnshamenergy
    write(ifu_log,*)'upKS in Hartrees',hartreeupks
    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Occupied beta eigenvalues in Hartrees',downocckscount
    !print '(1000(ES14.4))', ( DdownOHD(j), j=1,downocckscount )
    print *, ( DdownOHD(j)/evperhartree, j=1,downocckscount )
    write(ifu_log,*)'down Kohn-Sham energy sum',downkohnshamenergy
    write(ifu_log,*)'downKS in Hartrees',hartreedownks
    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''
    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
    write(ifu_log,*)'Total Kohn-Sham energy sum',kohnshamenergy
    write(ifu_log,*)'KS in Hartrees',hartreeks
    write(ifu_log,*)'****************************************************************'
    write(ifu_log,*)''

!    write(ifu_log,*)'************************************************************'
!    write(ifu_log,*)'SHIFT = ', -SHIFT
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied alpha eigenvalues in Hartrees (shifted)',upocckscount
!    !print '(1000(ES14.4))', ( DupOHD(j)-shift, j=1,upocckscount )
!    print *, ( (DupOHD(j)-shift)/evperhartree, j=1,upocckscount )
!    write(ifu_log,*)'up Kohn-Sham energy sum (shifted)',upkohnshamenergy-upocckscount*shift
!    write(ifu_log,*)'upKS in Hartrees (shifted)',hartreeupks-upocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Occupied beta eigenvalues in Hartrees (shifted)',downocckscount
!    !print '(1000(ES14.4))', ( DdownOHD(j)-shift, j=1,downocckscount )
!    print *, ( (DdownOHD(j)-shift)/evperhartree, j=1,downocckscount )
!    write(ifu_log,*)'down Kohn-Sham energy sum (shifted)',downkohnshamenergy-downocckscount*shift
!    write(ifu_log,*)'downKS in Hartrees (shifted)',hartreedownks-downocckscount*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''
!    write(ifu_log,*)'Total Kohn-Sham eigenvalues',norbks
!    write(ifu_log,*)'Total Kohn-Sham energy sum (shifted)',kohnshamenergy-norbks*shift
!    write(ifu_log,*)'KS in Hartrees (shifted)',hartreeks-norbks*shift/evperhartree
!    write(ifu_log,*)'****************************************************************'
!    write(ifu_log,*)''



      !do ispin=1,NSpin
      !  do j=1,6
      !    if(NSpin==1) write(ifu_log,10111)'Elem:',j
      !    if(NSpin==1) write(ifu_log,10112)(upOHD(j,i),i=1,3)
      !    if(NSpin==2) write(ifu_log,10111)'Elem:',j
      !    if(NSpin==2) write(ifu_log,10112)(downOHD(j,i),i=1,3)
      !  end do
      !end do

      !call DiagFullHamiltonian( HD, SD)
      !call DiagHamiltonian
      deallocate( OHD, kseigenvalues, allkseigenvalues, upkseigenvalues, downkseigenvalues, upHD, downHD, upOHD, downOHD, DupOHD, DdownOHD, VupOHD, VdownOHD )
      !deallocate( shiftallkseigenvalues, shiftupkseigenvalues, shiftdownkseigenvalues, shiftDupOHD, shiftDdownOHD )
  end subroutine WorkFock

  subroutine PrintCMatrix( A )
    implicit none
    complex*16,dimension(:,:),intent(in) :: A
    integer :: i,j,dim1,dim2
    dim1 = size(A,1)
    dim2 = size(A,2)
    do i=1,dim1
       !print '(1000(ES14.4))', ( A(i,j), j=1,dim2 )
       print '(100(g15.5,g15.5,2x))', ( real(A(i,j)),&
                                        AIMAG(A(i,j)), j=1,dim2 )
    end do
  end subroutine PrintCMatrix
!
!  subroutine PrintRMatrix( A )
!    implicit none
!    real*8,dimension(:,:),intent(in) :: A
!    integer :: i,j,dim1,dim2
!    dim1 = size(A,1)
!    dim2 = size(A,2)
!    do i=1,dim1
!       !print '(1000(ES14.4))', ( A(i,j), j=1,dim2 )
!       print '(100(g15.5))', ( A(i,j), j=1,dim2 )
!    end do
!  end subroutine PrintRMatrix

  subroutine WriteANT1DInput
    use parameters, only: eta
    use AntCommon
    implicit none

    real :: dsmall
    integer :: i,j

    CHARACTER(len=55) :: fname

    dsmall = eta
    
    fname='dev.'//trim(ant1dname)//'.dat'
    open(ifu_ant,file=fname,status='unknown')

    write(ifu_ant,'(A)'),       "&DevParams"
    write(ifu_ant,'(A,I1)'),    "NDSpin = ", NSpin
    write(ifu_ant,'(A,I4)'),    "NDAO = ", NAOrbs
    write(ifu_ant,'(A,I4)'),    "NDEl = ", NCDEl
    write(ifu_ant,'(A,F12.8)'), "EFermi = ", -shift
    write(ifu_ant,'(A)'),       "sparse = .true."
    write(ifu_ant,'(A)'),       "/"
    write(ifu_ant,*)
    write(ifu_ant,'(A)'),       "! Hamiltonian"
    do ispin=1,NSpin
       if( NSpin == 2 .and. ispin == 1 ) write(ifu_ant,'(A)'),       "! Spin-up"
       if( NSpin == 2 .and. ispin == 2 ) write(ifu_ant,'(A)'),       "! Spin-down"
       do i=1,NAOrbs
          do j=1,NAOrbs
             if(abs(HD(ispin,i,j))>=dsmall) write(ifu_ant,'(I6,I6,ES20.8)'), i, j, HD(ispin,i,j)
          enddo
       enddo
       write(ifu_ant,'(I6,I6,ES20.8)'), 0, 0, 0.0d0
       write(ifu_ant,*)
    end do
    write(ifu_ant,'(A)'),       "! Overlap"
    do i=1,NAOrbs
       do j=1,NAOrbs
          if(abs(SD(i,j))>=dsmall) write(ifu_ant,'(I6,I6,ES20.8)'), i, j, SD(i,j)
       enddo
    enddo
    write(ifu_ant,'(I6,I6,ES20.8)'), 0, 0, 0.0d0
    write(ifu_ant,*)
    close(ifu_ant)
  end subroutine WriteANT1DInput
  
  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensMat(ADDP, boolComputeFock, boolComputeCGibbsY)
    use numeric, only: Secant, Muller, Muller_OMP, BISEC
    use g09Common, only: GetNAE, GetNBE
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    logical,intent(in) :: boolComputeFock, boolComputeCGibbsY
    integer :: i,j, k,cond
    real :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    real :: TotChargeFromGlesser
    logical :: root_fail
    
    Write(*,'(A)')"ENTERED Device/CompDensMat"

    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    if( NSpin == 2 .and. SPINLOCK )then
       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
       do ispin=1,NSpin
          root_fail = .true.
          if (ispin.eq.1) E1=shiftup
          if (ispin.eq.2) E1=shiftdown
          E0=E1-Z
          E2=E1+Z
          if( root_fail )then
             print*,'MULLER method'
             call MULLER(F,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             !call MULLER_OMP(F,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
                print *, 'Warning: MULLER method failed to find root. Using SECANT.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if( root_fail )then
             print*,'SECANT method'
             call SECANT(F,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
                print *, 'Warning: SECANT method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if (root_fail) then
             print *, 'BISEC method'
             if (ispin.eq.1) shiftup = BISEC(F,EMin,EMax,Delta,5*Max,K)
             if (ispin.eq.2) shiftdown = BISEC(F,EMin,EMax,Delta,5*Max,K)
             DE=Delta
             if(k.lt.5*Max) root_fail = .false.
             if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
          end if
          write(ifu_log,*)'--------------------------------------------------------'
          if (ispin.eq.1) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
          end if
          if (ispin.eq.2) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
          end if
          write(ifu_log,*)'--------------------------------------------------------'
       end do
    else
       root_fail = .true.
       E0=shift-Z 
       E1=shift
       E2=shift+Z
       if (root_fail) then
          print*,'MULLER method'
          if(boolComputeFock)then
            call MULLER(QHWTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          else
            call MULLER(QXTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          end if
          !call MULLER_OMP(QXTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
             print *, 'Warning: MULLER method failed to find root. Using SECANT.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print*,'SECANT method'
          if(boolComputeFock)then
            call SECANT(QHWTot,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          else
            call SECANT(QXTot,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          end if
          if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: SECANT method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print *, 'BISEC method'
          if(boolComputeFock)then
            shift = BISEC(QHWTot,EMin,EMax,Delta,5*Max,K)
          else
            shift = BISEC(QXTot,EMin,EMax,Delta,5*Max,K)
          end if
          DE=Delta
          if(k.lt.5*Max) root_fail = .false.
          if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
       end if

       if(boolComputeCGibbsY)then
         ! MUST NOT CALL MULLER METHOD RELYING ON  QCGibbsYTot BECAUSE IT USES GLESSER IN THE FULL ENERGY RANGE.
         !call MULLER(QCGibbsYTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
         write(ifu_log,*)''
         write(ifu_log,*)''
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'--- BEFORE QCGibbsYTot(shift) ---'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'--- boolComputeCGibbsY IS true ---'
         write(ifu_log,*)'--- ENTERING QCGibbsYTot(shift) ---'
         TotChargeFromGlesser = QCGibbsYTot(shift)
         write(ifu_log,'(A,F10.5)') ' TotChargeFromGlesser:  ', TotChargeFromGlesser
         write(ifu_log,*)'--- AFTER QCGibbsYTot(shift) ---'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)''
         write(ifu_log,*)''
       else
         write(ifu_log,*)''
         write(ifu_log,*)''
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'--- boolComputeCGibbsY IS false ---'
         write(ifu_log,*)'--- NOT ENTERING QCGibbsYTot(shift) ---'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)'-----------------------------------------------'
         write(ifu_log,*)''
         write(ifu_log,*)''
         TotChargeFromGlesser = 0.0d0
       end if

       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
       write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
       write(ifu_log,*)'-----------------------------------------------'
    end if
    ADDP = .not. root_fail
  end subroutine CompDensMat


subroutine CompLocalMu(ADDP)
    !use numeric, only: LocalSecant, LocalMuller, LocalBISEC
    use numeric, only: Secant, Muller, BISEC
    use g09Common, only: GetNAE, GetNBE
    use cluster, only: LoAOrbNo, HiAOrbNo, LoAOrbNo, HiAOrbNo
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess,BiasVoltage,SPIN_Beg,SPIN_End
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail

!-------------------------------------------------------
    integer :: n,l
    real :: sdeg, ro_a, ro_b, chargeregion, spinregion
    real, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b, tmp
!-------------------------------------------------------
    
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'------- I am in CompLocalMu ---------'
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Mulliken population analysis ---'
    write(ifu_log,*)'-------------------------------------'


    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    rho_a = matmul( PD(1,:,:), SD )
    if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )

    write(ifu_log,*)'---------------------------------------------------------------------'
    write(ifu_log,*)'Charges in selected region to compute local electrochemical potential'
    write(ifu_log,*)'---------------------------------------------------------------------'

    chargeregion=0.0d0
    spinregion=0.0d0
    do j=SPIN_Beg,SPIN_End
       ro_a=0.0d0
       ro_b=0.0d0
       do i=LoAOrbNo(j),HiAOrbNo(j)
          ro_a=ro_a+rho_a(i,i)
          !if(NSpin==1) ro_b=ro_a
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       alphaelec(j-SPIN_Beg+1)=ro_a
       if(NSpin==1) betaelec(j-SPIN_Beg+1)=ro_a
       if(NSpin==2) betaelec(j-SPIN_Beg+1)=ro_b
       if(NSpin ==1 ) write(ifu_log,'(A,I5,A,f9.5)')'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin ==2 ) write(ifu_log,'(A,I5,A,f9.5,A,f9.5)')'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)

       chargeregion=chargeregion+(ro_a+ro_b)
       if (NSpin == 2) spinregion=spinregion+(ro_a-ro_b)
    end do
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in selected region:',chargeregion*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in selected region:',spinregion
   write(ifu_log,*)'----------------------------------------------------'

!-------------------------------------------------------------------------------------------------
    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    !Epsilon=ChargeAcc*(NCDEl+QExcess)
    Epsilon=ChargeAcc*(chargeregion)

!-------------------------------------------------------------------------------------------------

    !--- ADDED BY CARLOS -----------------------------------------------
    print*,'BiasVoltage = ', BiasVoltage
    !-------------------------------------------------------------------
!if( NSpin == 2 .and. SPINLOCK )then
if( NSpin == 2 )then
!  print*,'SPINLOCK is true'
  print*,'NSpin = 2'
  print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
  do spinatom=SPIN_Beg,SPIN_End
    alphalocalshift(spinatom) = shiftup
    betalocalshift(spinatom) = shiftdown
    do ispin=1,NSpin
       root_fail = .true.
       !if (ispin.eq.1) E1=shiftup
       !if (ispin.eq.2) E1=shiftdown
       if (ispin.eq.1) E1=alphalocalshift(spinatom)
       if (ispin.eq.2) E1=betalocalshift(spinatom)
       E0=E1-Z
       E2=E1+Z
       if( root_fail )then
          print*,'MULLER method'
          call MULLER(FPart,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: MULLER method failed to find root. Using SECANT.'
             root_fail = .true.
          else
             !if (ispin.eq.1)shiftup=E3
             !if (ispin.eq.2)shiftdown=E3
             if (ispin.eq.1)alphalocalshift(spinatom)=E3
             if (ispin.eq.2)betalocalshift(spinatom)=E3
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print*,'SECANT method'
          call SECANT(FPart,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: SECANT method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             !if (ispin.eq.1)shiftup=E3
             !if (ispin.eq.2)shiftdown=E3
             if (ispin.eq.1)alphalocalshift(spinatom)=E3
             if (ispin.eq.2)betalocalshift(spinatom)=E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print *, 'BISEC method'
          !if (ispin.eq.1) shiftup = LocalBISEC(FPart,EMin,EMax,Delta,5*Max,K)
          !if (ispin.eq.2) shiftdown = LocalBISEC(FPart,EMin,EMax,Delta,5*Max,K)
          if (ispin.eq.1) alphalocalshift(spinatom) = BISEC(FPart,EMin,EMax,Delta,5*Max,K)
          if (ispin.eq.2) betalocalshift(spinatom) = BISEC(FPart,EMin,EMax,Delta,5*Max,K)          
          DE=Delta
          if(k.lt.5*Max) root_fail = .false.
          if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
       end if
       write(ifu_log,*)'--------------------------------------------------------'
       if (ispin.eq.1) then
          write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -alphalocalshift(spinatom),'  +/-',dabs(DE)
          write(ifu_log,*)
          write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', LocalQAlpha
       end if
       if (ispin.eq.2) then
          write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -betalocalshift(spinatom),'  +/-',dabs(DE)
          write(ifu_log,*)
          write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', LocalQBeta
          write(ifu_log,*)
          write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', LocalQAlpha+LocalQBeta
       end if
       write(ifu_log,*)'--------------------------------------------------------'
    end do
  end do
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'-----------------------------------------------'
  do ispin=1,NSpin 
    do spinatom=SPIN_Beg,SPIN_End
      !write(ifu_log,'(A,F9.5,A,f9.5)') ' Local Mu= ',-alphalocalshift(spinatom),'  +/-',dabs(DE)
      if (ispin.eq.1)write(ifu_log,'(A,I5,A,F9.5)') 'Atom',spinatom,'	Alpha Local Mu:  ',-alphalocalshift(spinatom)
      if (ispin.eq.2)write(ifu_log,'(A,I5,A,F9.5)') 'Atom',spinatom,'	Beta Local Mu:  ',-betalocalshift(spinatom)
    end do
    write(ifu_log,*)'-----------------------------------------------'
  end do
else
  do spinatom=SPIN_Beg,SPIN_End
    !print*,'SPINLOCK is false'
    print*,'NSpin = 1'
    !alphalocalshift(spinatom) = shiftup
    !betalocalshift(spinatom) = shiftdown
    alphalocalshift(spinatom) = shift
    root_fail = .true.
    E0=alphalocalshift(spinatom)-Z 
    E1=alphalocalshift(spinatom)
    E2=alphalocalshift(spinatom)+Z
    if (root_fail) then
      print*,'MULLER method'
      !call MULLER(QXPart,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
      call MULLER(QXPart,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
      if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
         print *, 'Warning: MULLER method failed to find root. Using SECANT.'
         root_fail = .true.
      else
         alphalocalshift(spinatom) = E3
         root_fail = .false.
      end if
    end if
    if (root_fail) then
      print*,'SECANT method'
      !call SECANT(QXPart,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
      call SECANT(QXPart,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
      if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
         print *, 'Warning: SECANT method failed to find root. Using BISEC.'
         root_fail = .true.
      else
         alphalocalshift(spinatom) = E3
         root_fail = .false.
      end if
    end if
    if (root_fail) then
      print *, 'BISEC method'
      !Delta = 10*Delta !Added by Carlos Salgado
      !Max = 10*Max !Added by Carlos Salgado
      !alphalocalshift(spinatom) = BISEC(QXPart,EMin,EMax,Delta,5*Max,K)
      alphalocalshift(spinatom) = BISEC(QXPart,EMin,EMax,Delta,5*Max,K)
      !alphalocalshift(spinatom) = BISEC(QXPart,EMin,EMax,1.0D2*Delta,500*Max,K)
      DE=Delta
      if(k.lt.5*Max) root_fail = .false.
      if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
    end if

    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,'(A,F9.5,A,f9.5)') ' Alpha Local Mu:  ',-alphalocalshift(spinatom),'  +/-',dabs(DE)
    write(ifu_log,*)
    write(ifu_log,'(A,F10.5)') ' Number of alpha electrons:  ', LocalQAlpha
    write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', LocalQBeta
    write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', LocalQAlpha+LocalQBeta
    write(ifu_log,*)'-----------------------------------------------'

  end do
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'-----------------------------------------------'
  do spinatom=SPIN_Beg,SPIN_End
    !write(ifu_log,'(A,F9.5,A,f9.5)') ' Local Mu= ',-alphalocalshift(spinatom),'  +/-',dabs(DE)
    write(ifu_log,'(A,I5,A,F9.5)') 'Atom',spinatom,'	Local Mu:  ',-alphalocalshift(spinatom)
  end do
    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'-----------------------------------------------'

end if
  ADDP = .not. root_fail
end subroutine CompLocalMu
  
  ! 
  ! Computes the density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real function CompPD( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
    use g09Common, only: GetNAE, GetNBE
     implicit none

     ! chemical potential
     real, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     real, dimension(2*npmax) :: x, w
     real :: Ei, dEdx, Q, QQ, DPD, E0 !, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     
     real, parameter :: x0 = 0.5d0
     real :: aa, bb, cc, EM

     EM = EMax
     aa = EM/x0
     bb = aa*(1.0d0-x0)**2
     cc = EM - bb/(1-x0)

     shift = mu

     Q = d_zero

     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD=d_zero
        QAlpha = d_zero; QBeta = d_zero
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

        do ispin=1,NSpin
!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMP DO
           do i=1,np
              Ei = 2.0d0*EMax*x(i)
              if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
              dEdx = 2.0d0*EMax
              if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
              call GPlus0( ui*Ei, GD )
!$OMP CRITICAL
              do k=1,NAOrbs
                 do l=1,NAOrbs
                    DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                    PD(ispin,k,l) = PD(ispin,k,l) + DPD
                    if(ispin.eq.1) QAlpha = QAlpha + DPD*SD(l,k)
                    if(ispin.eq.2) QBeta = QBeta + DPD*SD(l,k)
                 end do
              end do
!$OMP END CRITICAL
           end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        end do
        if(NSpin.eq.1) QBeta=QAlpha
        Q = QAlpha + QBeta        
        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDEl ) exit
     end do
        
     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     CompPD = Q - dble(NCDEl)

  end function CompPD


  ! 
  ! Computes the density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real function CompEnergy( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
    use g09Common, only: GetNAE, GetNBE
     implicit none

     ! chemical potential
     real, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     !integer, parameter :: nmin=1, nmax=16, npmax=2**nmax-1
     real, dimension(2*npmax) :: x, w
     real :: Ei, dEdx, Q, QQ, DPD, E0 !, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     
     real, parameter :: x0 = 0.5d0
     real :: aa, bb, cc, EM

     real :: IntDOSE, IntDOSEAlpha, IntDOSEBeta, IntDOSEE, CompPD

     EM = EMax
     aa = EM/x0
     bb = aa*(1.0d0-x0)**2
     cc = EM - bb/(1-x0)
     
     print *, " - SHIFT:", -shift
     shift = mu ! ADDED BECAUSE IT IS IN OTHER FUNCTIONS
     print *, " - mu:", -mu
     !shift = mu ! I HAVE TO CLEAR THIS TO AVOID OVERWRITING SHIFT
     print *, " - SHIFT:", -shift
     Q = d_zero
     IntDOSE = d_zero

     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD=d_zero
        QAlpha = d_zero; QBeta = d_zero
        IntDOSEE = IntDOSE
        IntDOSEAlpha = d_zero; IntDOSEBeta = d_zero
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

        do ispin=1,NSpin
!$OMBLABLABLAP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMBLABLABLAP DO
           do i=1,np
              Ei = 2.0d0*EMax*x(i)
              if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
              dEdx = 2.0d0*EMax
              if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
              call GPlus0( ui*Ei, GD )
!$OMBLABLABLAP CRITICAL
              do k=1,NAOrbs
                 do l=1,NAOrbs
                    DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                    PD(ispin,k,l) = PD(ispin,k,l) + DPD
                    if(ispin.eq.1) QAlpha = QAlpha + DPD*SD(l,k)
                    if(ispin.eq.2) QBeta = QBeta + DPD*SD(l,k)
                    if(ispin.eq.1) IntDOSEAlpha = IntDOSEAlpha + (x(i)-mu)*DPD*SD(l,k)
                    if(ispin.eq.2) IntDOSEBeta = IntDOSEBeta + (x(i)-mu)*DPD*SD(l,k)
                 end do
              end do
!$OMBLABLABLAP END CRITICAL
           end do
           print *, ' Grid points:',np
           !print *, ' Charge alpha ', QAlpha
           !print *, ' Charge beta ', QBeta
           print *, ' Charge alpha + beta ', QAlpha + QBeta
           !print *, ' IntDOSE alpha ', IntDOSEAlpha
           !print *, ' IntDOSE beta ', IntDOSEBeta
           print *, ' IntDOSE alpha + beta ', IntDOSEAlpha + IntDOSEBeta
!$OMBLABLABLAP END DO
!$OMBLABLABLAP BARRIER
!$OMBLABLABLAP END PARALLEL
        end do
        if(NSpin.eq.1) QBeta=QAlpha
        Q = QAlpha + QBeta 
        if(NSpin.eq.1) IntDOSEBeta=IntDOSEAlpha  
        IntDOSE = IntDOSEAlpha + IntDOSEBeta     
        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDEl ) exit
        !if( n > nmin .and. abs(Q-QQ) < (1/16)*PAcc*NCDEl ) exit
     end do
        
     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     print *, ' Computed Charge alpha + beta ', Q

     CompPD = Q - dble(NCDEl)
     CompEnergy = IntDOSE

  end function CompEnergy
  
  ! 
  ! Computes the spin-resolved density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real function CompSpinPD( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
    use g09Common, only: GetNAE, GetNBE
     implicit none

     ! chemical potential
     real, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     real, dimension(2*npmax) :: x, w
     real :: Ei, dEdx, Q, QQ, DPD !, E0, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     real :: NCDAB

     ! Number of alpha OR beta electrons in central region C of device
     if(ispin.eq.1) NCDAB = dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE())
     if(ispin.eq.2) NCDAB = dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE())
     
     shift = mu

     Q=0.0d0
     do n=nmin,nmax
        np=2**n-1
        
        QQ = Q
        Q = 0.0d0
        PD(ispin,:,:) = 0.0d0
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMP DO
        do i=1,np
           Ei = 2.0d0*EMax*x(i)
           if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
           dEdx = 2.0d0*EMax
           if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
           call GPlus0( ui*Ei, GD )
!$OMP CRITICAL
           do k=1,NAOrbs
              do l=1,NAOrbs
                 DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                 PD(ispin,k,l) = PD(ispin,k,l) + DPD
                 Q = Q + DPD*SD(l,k)
              end do
           end do
!$OMP END CRITICAL
        end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDAB ) exit
     end do
        
     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     if(ispin.eq.1) QAlpha = Q
     if(ispin.eq.2) QBeta = Q

     CompSpinPD = Q - NCDAB

   end function CompSpinPD


  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensMat2(ADDP)
    use numeric, only: Secant, Muller, Muller_OMP, BISEC
    use g09Common, only: GetNAE, GetNBE
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail
    
    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    root_fail = .true.
    if( NSpin == 2 .and. SPINLOCK )then
       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
       do ispin=1,NSpin
          if (ispin.eq.1) E0=shiftup
          if (ispin.eq.2) E0=shiftdown
          E1=E0-Z
          E2=E0+Z 
          if( root_fail ) then
             print*,'Secant method'
             call SECANT(CompSpinPD,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
             if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
                print *, 'Warning: Secant method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E2
                if (ispin.eq.2)shiftdown=E2
                root_fail = .false.
             end if
          end if
          if( root_fail ) then
             print*,'Muller method'
             call MULLER(CompSpinPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             !call MULLER_OMP(CompSpinPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
                print *, 'Warning: Muller method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if( root_fail )then
             print *, 'BISEC method'
             if( ispin.eq.1) shiftup = BISEC(CompSpinPD,EMin,EMax,Delta,Max,K)
             if( ispin.eq.2) shiftdown = BISEC(CompSpinPD,EMin,EMax,Delta,Max,K)
             DE=Delta
             if(k.lt.Max) root_fail = .false.
          end if
          write(ifu_log,*)'--------------------------------------------------------'
          if (ispin.eq.1) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
          end if
          if (ispin.eq.2) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
          end if
          write(ifu_log,*)'--------------------------------------------------------'
       end do
    else
       E0=shift
       E1=E0-Z 
       E2=E0+Z 
       if( root_fail )then
          print*,'Secant method'
          call SECANT(CompPD,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
          if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
             print *, 'Warning: Secant method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E2
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print*,'Muller method'
          call MULLER(CompPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          !call MULLER_OMP(CompPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
             print *, 'Warning: Muller method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print *, 'BISEC method'
          print *, EMin, EMax
          shift = BISEC(CompPD,EMin,EMax,Delta,Max,K)
          DE=Delta
          if(k.lt.Max) root_fail = .false.
       end if
       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
       write(ifu_log,*)'-----------------------------------------------'
    end if
    ADDP = .not. root_fail
  end subroutine CompDensMat2

!------------------------------------------------------------------------------------------------
!--------------------- MULLSTORE TO SORT AND STORE MULLPOP --------------------------------------
!------------------------------------------------------------------------------------------------
  subroutine MullStore
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol
    USE parameters, only: Mulliken, LDOS_Beg, LDOS_End
    use g09Common, only: GetAtmCo
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n, l
    real :: sdeg, ro_a, ro_b, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b, tmp

    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Mulliken population analysis ---'
    write(ifu_log,*)'-------------------------------------'


    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    rho_a = matmul( PD(1,:,:), SD )
    if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 1'
    write(ifu_log,*)'----------------------'
    I1=0
    chargelead1=0.0
    spinlead1=0.0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin ==2 ) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead1=chargelead1+(ro_a+ro_b)
       if (NSpin == 2) spinlead1=spinlead1+(ro_a-ro_b)
       I1 = I1 + NAOAtom(j)
    end do
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 1:',chargelead1*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 1:',spinlead1
   write(ifu_log,*)'----------------------------------------------------'

    chargemol=0.0d0
    spinmol=0.0d0
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'-------------------'
       write(ifu_log,*)'Charges in molecule'
       write(ifu_log,*)'-------------------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+rho_a(i,i)
             if(NSpin==2) ro_b=ro_b+rho_b(i,i)
          end do
          if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',sdeg*ro_a
          if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
          IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
          END IF
          IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
          END IF
          chargemol=chargemol+ro_a+ro_b
          if (NSpin == 2) spinmol=spinmol+(ro_a-ro_b)
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in molecule:',chargemol*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in molecule:',spinmol
   write(ifu_log,*)'----------------------------------------------------'

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 2'
    write(ifu_log,*)'----------------------'
    chargelead2=0.0
    spinlead2=0.0
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead2=chargelead2+ro_a+ro_b
       if (NSpin == 2) spinlead2=spinlead2+(ro_a-ro_b)
    end do

   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 2:',chargelead2*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 2:',spinlead2
   write(ifu_log,*)'----------------------------------------------------'

1011 format(a6,i4,a10,f8.4)
1012 format(a6,i4,a10,f8.4,a10,2f9.4)
1013 format(7f9.4)
  end subroutine MullStore

  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensLocal2(ADDP)
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, LoAOrbNo, HiAOrbNo
    use numeric, only: Secant, Muller, BISEC
    use g09Common, only: GetNAE, GetNBE, GetNAtoms, GetAtmCo

    use constants, only: bohr
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess, NSpinMuBl, SpinMuBlBeg, SpinMuBlEnd, Mulliken, Ldos_beg, Ldos_end
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail

    !---------------------------------------------------------------------------------------------
    !---------- ADDED BY CARLOS TO SORT AND STORE THE MULLPOP ------------------------------------
    !---------------------------------------------------------------------------------------------
    integer :: I1, is ,n, l
    real :: sdeg, ro_a, ro_b, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b, tmp

    integer :: iblock, iao, jao, ispin
    integer :: iatbl
    !---------------------------------------------------------------------------------------------

    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'- Charges in Block to calculate the spatially -'
    write(ifu_log,*)'---- and spin resolved chemical potential -----'
    write(ifu_log,*)'-----------------------------------------------'


    !--- MULLIKEN POPULATION -----------------------
    rho_a = matmul( PD(1,:,:), SD )
    if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )
    !-----------------------------------------------

    I1=0

    do iblock = 1,NSpinMuBl
      do iatbl = SpinMuBlBeg(iblock),SpinMuBlEnd(iblock)
        I1=LoAOrbNo(iatbl)
        ro_a=0.0d0
        ro_b=0.0d0
        !do i=I1,I1+NAOAtom(j)-1
        do i=LoAOrbNo(iatbl),HiAOrbNo(iatbl)
           ro_a=ro_a+rho_a(i,i)
           if(NSpin==2) ro_b=ro_b+rho_b(i,i)
        end do
        if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
        if(NSpin ==2 ) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)

        

        IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
          if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
          if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
        END IF
        IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
          if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
          if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
        END IF

        !chargelead1=chargelead1+(ro_a+ro_b)
        !if (NSpin == 2) spinlead1=spinlead1+(ro_a-ro_b)

        !I1 = I1 + NAOAtom(j)




      end do
    end do


1011 format(a6,i4,a10,f8.4)
1012 format(a6,i4,a10,f8.4,a10,2f9.4)
1013 format(7f9.4) 

    ADDP = .not. root_fail
  end subroutine CompDensLocal2

  ! 
  ! Computes the spin-resolved density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
!  real function CompLocalSpinPD( iatom, iatommu, iatomQ )
!     use constants
!     use util
!     use numeric, only: CHDiag, gauleg, RTrace
!    !-------------------------------------
!    use Cluster, only : hiaorbno, loaorbno
!    !-------------------------------------
!     use parameters, only: eta, PAcc
!#ifdef G03ROOT
!    use g03Common, only: GetNAE, GetNBE
!#endif
!#ifdef 1
!    use g09Common, only: GetNAE, GetNBE
!#endif
!     implicit none
!
!     ! chemical potential
!     real, intent(in) :: iatommu
!     integer, intent(in) :: iatom
!
!     complex*16, dimension(NAOrbs,NAOrbs) :: GD
!
!     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
!     real, dimension(2*npmax) :: x, w
!     real :: Ei, dEdx, Q, QQ, DPD !, E0, Qi
!     integer :: n, np, i, j, k, l !, info, ierr
!     real :: NCDAB
!
!     ! Number of alpha OR beta electrons in central region C of device
!     !if(ispin.eq.1) NCDAB = dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE())
!     !if(ispin.eq.2) NCDAB = dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE())
!     ! WE DON'T INTEGRATE ANY MORE TO THE FIXED CHARGE OF THE ELEMENT
!     ! NOW WE INTEGRATE TO THE LOCAL KNOWN CHARGE CALCULATED BEFORE iatomq
!     ! Is shift we use locally defined by iatommu? In principle, yes.
!     
!     !shift = mu
!     iatomshift = iatommu
!
!     !Q=0.0d0
!     iatomQ=0.0d0
!     do n=nmin,nmax
!        np=2**n-1
!        
!        !QQ = Q
!        !Q = 0.0d0
!        !PD(ispin,:,:) = 0.0d0
!        iatomQQ = iatomQ
!        iatomQ = 0.0d0
!        iatomdens(ispin,:,:) = 0.0d0
!        
!        ! Compute Gauss-Legendre abcsissas and weights
!        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)
!
!!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!!$OMP DO
!        do i=1,np
!           Ei = 2.0d0*EMax*x(i)
!           if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
!           dEdx = 2.0d0*EMax
!           if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
!           call GPlus0( ui*Ei, GD )
!!$OMP CRITICAL
!           do k=1,NAOrbs
!              do l=1,NAOrbs
!                 DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
!                 PD(ispin,k,l) = PD(ispin,k,l) + DPD
!                 Q = Q + DPD*SD(l,k)
!              end do
!           end do
!!$OMP END CRITICAL
!        end do
!!$OMP END DO
!!$OMP BARRIER
!!$OMP END PARALLEL
!        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDAB ) exit
!     end do
!        
!     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
!     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 
!
!     if(ispin.eq.1) QAlpha = Q
!     if(ispin.eq.2) QBeta = Q
!
!     CompSpinPD = Q - NCDAB
!
!!-------------------------------------------------------------------------
!SG = matmul( SD, green )
!          ! computing total DOS
!          DOS=d_zero
!          AtomDOS=d_zero
!          do j=1,GetNAtoms()
!          do i=LoAOrbNo(j),HiAOrbNo(j)
!             AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
!             DOS=DOS-dimag(SG(i,i))/d_pi
!          end do
!          end do
!!-------------------------------------------------------------------------
!
!   end function CompLocalSpinPD


  !*************************************
  !* Compute retarded Green's function *
  !*************************************
!  subroutine gplus0Subspace(z,partgreen,iniOrb,endOrb)
!    use PARAMETERS, only: eta,glue
!    use constants, only: c_zero, ui
!    use lapack_blas, only: zgetri, zgetrf
!
!    implicit none
!
!    integer :: i, j, info, omp_get_thread_num
!    integer, dimension(NAOrbs) :: ipiv
!    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
!    complex*16 :: work(4*NAOrbs) 
!
!    complex*16, intent(in) :: z 
!
!    integer, intent(in) :: iniOrb, endOrb
!    
!    integer :: partdim
!    partdim = endOrb - iniOrb
!    
!    !complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: partgreen
!    complex*16, dimension(partdim,partdim), intent(out) :: partgreen
!
!    ! Initilization 
!    partgreen=c_zero
!    partsigr=-ui*eta*SD(iniorb:endOrb,iniOrb:endOrb) 
!    partsigl=-ui*eta*SD(iniorb:endOrb,iniOrb:endOrb) 
!
!    call CompSelfEnergies( ispin, z, partsigl, partsigr )
!    partsigr=glue*partsigr
!    partsigl=glue*partsigl
!    
!    !************************************************************************
!    !c Retarded "Green" function
!    !************************************************************************
!    do i=1,NAOrbs
!       do j=1,NAOrbs
!          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
!       enddo
!    enddo
!
!    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
!    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)
!
!  end subroutine gplus0Subspace

  
  !*************************************
  !* Compute retarded Green's function *
  !*************************************
  subroutine gplus0(z,green)
    use PARAMETERS, only: eta,glue
    use constants, only: c_zero, ui
   !use lapack_blas, only: zgetri, zgetrf
   !use lapack95, only: zgetri, zgetrf
   !use lapack95
   !use blas95

    implicit none
    external               zgetri, zgetrf

    integer :: i, j, info, omp_get_thread_num
    integer, dimension(NAOrbs) :: ipiv
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
    complex*16 :: work(4*NAOrbs) 

    complex*16, intent(in) :: z 

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: green
    

    ! Initilization 
    green=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

  end subroutine gplus0


  !*************************************
  !* Compute retarded Green's function *
  !*************************************
  subroutine gplus(z,green,gammar,gammal)
    use PARAMETERS, only: eta, glue
    use constants, only: c_zero, ui
   !use lapack_blas, only: zgetri, zgetrf
   !use lapack95, only: zgetri, zgetrf
   !use lapack95
   !use blas95

    implicit none
    external              zgetri, zgetrf
!    external              cgetri, cgetrf
!    external              LAPACKE_zgetrf

    integer :: i, j, omp_get_thread_num
    integer*8 :: info
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: green
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gammar,gammal
    

    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo
           !call zgetrf(NDIM,NDIM,Sigma_aux,NDIM,ipiv,info) ! ORIGINAL FOR PGI IN BetheLattice.f90.
           !call zgetrf(NDim,NDim,Sigma_aux,NDim+1,ipiv,info) ! ADDING +1 TO PARAMETER 4 WORKS!!!
    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info) ! ORIGINAL CODE FOR ANT.G IN device.f90.
!    IF( INFO /= 0 ) THEN
!      WRITE(ifu_log,*)'Problems using zgetrf in device.f90'
!      WRITE(ifu_log,*)'INFO=',INFO
!      IF( INFO<0 )THEN
!        WRITE(ifu_log,'(I2,A)')-INFO,"th parameter had an illegal value."
!      ELSE
!        WRITE(ifu_log,'(I8,A)')INFO,"th parameter is exactly 0. The matrix is singular."
!      END IF
!      !STOP
!    END IF
           !call zgetri(NDIM,Sigma_aux,NDIM,ipiv,work,4*NDIM,info) ! ORIGINAL FOR PGI IN BetheLattice.f90.
           !call zgetri(NDim,Sigma_aux,NDim+1,ipiv,work,4*NDim+1,info) ! ADDING +1 TO PARAMETERS 3,6 WORKS!!!
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info) ! ORIGINAL CODE FOR ANT.G IN device.f90.
!    IF( INFO /= 0 ) THEN
!      WRITE(ifu_log,*)'Problems using zgetri in device.f90'
!      WRITE(ifu_log,*)'INFO=',INFO
!      IF( INFO<0 )THEN
!        WRITE(ifu_log,'(I2,A)')-INFO,"th parameter had an illegal value."
!      ELSE
!        WRITE(ifu_log,'(I2,A)')INFO,"th parameter is exactly 0. The matrix is singularand the inversion could not be completed."
!      END IF
!      !STOP
!    END IF

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus

  !*************************************
  !* Compute lesser Green's function *
  !*************************************
  subroutine glesser(z,gless)
    use parameters, only: eta, biasvoltage, glue
    use constants, only: c_zero, ui, c_one
    use util, only: PrintCMatrix
   !use lapack_blas, only: zgetri, zgetrf, zgemm
   !use lapack95, only: zgetri, zgetrf, zgemm
   !use lapack95
   !use blas95

    implicit none
    external               zgetri, zgetrf, zgemm

    integer :: i, j, info, error
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)

    complex*16, intent(in) :: z

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gless
    complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

    if(DebugDev)Write(*,'(A)')"ENTERED Device/glesser"

     allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then 
        print*,"Problems allocating" 
        stop
     end if
     allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then 
        print*,"Problems allocating" 
        stop
     end if

    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     ! argument 'C' means zgemm uses the conjugate hermitian of green (retarded), which is the advanced.
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammal,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessL,NAOrbs)
     ! argument 'C' means zgemm uses the conjugate hermitian of green (retarded), which is the advanced.
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammar,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessR,NAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

     deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
     if (error /= 0 ) then 
        print*,"Problems deallocating" 
        stop
     end if

     if(DebugDev)then
       Write(*,'(A)')"Device/glesser output gless(:,:) = "
       call PrintCMatrix(gless)
     end if

  end subroutine glesser

  !*************************************
  !* Compute lesser Green's function *
  !*************************************
  subroutine CompCGibbsY(z,gless,GibbsY,GibbsYKernel1,GibbsYKernel2) !NoHanWorks20180425
    use parameters, only: eta, biasvoltage, glue
    use constants, only: c_zero, ui, c_one
    use util, only: PrintCMatrix
   !use lapack_blas, only: zgetri, zgetrf, zgemm
   !use lapack95, only: zgetri, zgetrf, zgemm
   !use lapack95
   !use blas95

    implicit none
    external               zgetri, zgetrf, zgemm

    integer :: i, j, info, error
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)

    complex*16, intent(in) :: z

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gless
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: GibbsY, GibbsYKernel1, GibbsYKernel2
    complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

    complex*16, dimension(:,:), allocatable :: auxPD, aux1Kernel1, aux2Kernel1, aux1Kernel2, aux2Kernel2!, glessKernelL, glessKernelR

     if(DebugDev)then
    Write(*,*)"********************************************************************"
    Write(*,*)"*********************** ENTERED CompCGibbsY ************************"
    Write(*,*)"********************************************************************"
     end if

    if(DebugDev)Write(*,'(A)')"ENTERED Device/CompCGibbsY"

     allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then
        print*,"Problems allocating"
        stop
     end if
     allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then
        print*,"Problems allocating"
        stop
     end if
     allocate(auxPD(NAOrbs,NAOrbs), aux1Kernel1(NAOrbs,NAOrbs), aux2Kernel1(NAOrbs,NAOrbs), aux1Kernel2(NAOrbs,NAOrbs), aux2Kernel2(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then
        print*,"Problems allocating"
        stop
     end if

    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD
    sigl=-ui*eta*SD

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl

    !************************************************************************
    !c Coupling matrices
    !************************************************************************
    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     ! argument 'C' means zgemm uses the conjugate hermitian of green (retarded), which is the advanced.
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammal,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessL,NAOrbs)
     ! argument 'C' means zgemm uses the conjugate hermitian of green (retarded), which is the advanced.
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammar,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessR,NAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

     !GibbsYL=(Fermi-biasvoltage/2)*ui*(glessL*theta(dble(z)-(biasvoltage/2)))
     !GibbsYR=(Fermi+biasvoltage/2)*ui*(glessR*theta(dble(z)-(-biasvoltage/2)))
     ! REMOVE Fermi because everything is shifted to ZERO.
     !GibbsYL=(-biasvoltage/2)*ui*(glessL*theta(dble(z)-(biasvoltage/2)))
     !GibbsYR=(+biasvoltage/2)*ui*(glessR*theta(dble(z)-(-biasvoltage/2)))

     ! ADD OR SUBSTRACT shift TO PUT HERE THE CORRECT FERMI LEVEL
     !GibbsY = - ( (-shift-biasvoltage/2)*(glessL*theta(dble(z)-(biasvoltage/2))) &
     !         & + (-shift+biasvoltage/2)*(glessR*theta(dble(z)-(-biasvoltage/2))) )
     glessL = (glessL*theta(dble(z)-(biasvoltage/2)))
     glessR = (glessR*theta(dble(z)-(-biasvoltage/2)))
     GibbsY = - ( (-shift-biasvoltage/2)*glessL + (-shift+biasvoltage/2)*glessR )


     do i=1,NAOrbs
       do j=1,NAOrbs
         auxPD(i,j)=PD(ispin,i,j)
         aux1Kernel1(i,j) = (0.0,0.0)
         aux2Kernel1(i,j) = (0.0,0.0)
       end do
     end do

     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,auxPD,NAOrbs,glessL,NAOrbs, &
          &           c_zero,aux1Kernel1,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,glessL,NAOrbs,aux1Kernel1,NAOrbs, &
          &           c_zero,aux2Kernel1,NAOrbs)

     GibbsYKernel1 = (-shift-biasvoltage/2)*(z-shift)*aux2Kernel1;
     GibbsYKernel2 = (-shift-biasvoltage/2)*( - aux2Kernel1 );

     do i=1,NAOrbs
       do j=1,NAOrbs
         auxPD(i,j)=PD(ispin,i,j)
         aux1Kernel1(i,j) = (0.0,0.0)
         aux2Kernel1(i,j) = (0.0,0.0)
       end do
     end do

     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,auxPD,NAOrbs,glessR,NAOrbs, &
          &           c_zero,aux1Kernel1,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,glessR,NAOrbs,aux1Kernel1,NAOrbs, &
          &           c_zero,aux2Kernel1,NAOrbs)

     GibbsYKernel1 = GibbsYKernel1 + (-shift+biasvoltage/2)*(z-shift)*aux2Kernel1;
     GibbsYKernel2 = GibbsYKernel2 + (-shift+biasvoltage/2)*( - aux2Kernel1 );


     deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
     if (error /= 0 ) then
        print*,"Problems deallocating"
        stop
     end if
     deallocate(auxPD, aux1Kernel1, aux2Kernel1, aux1Kernel2, aux2Kernel2,stat=error)
     if (error /= 0 ) then
        print*,"Problems deallocating"
        stop
     end if

     if(DebugDev)then
       Write(*,'(A)')"Device/CompCGibbsY output gless(:,:) = "
       call PrintCMatrix(gless)

!       Write(*,'(A)')"Device/CompCGibbsY output GibbsY(:,:) = "
!       call PrintCMatrix(GibbsY)
!
!       Write(*,'(A)')"Device/CompCGibbsY output GibbsYKernel1(:,:) = "
!       call PrintCMatrix(GibbsYKernel1)
!
!       Write(*,'(A)')"Device/CompCGibbsY output GibbsYKernel2(:,:) = "
!       call PrintCMatrix(GibbsYKernel2)

    Write(*,*)"********************************************************************"
    Write(*,*)"************************* EXIT CompCGibbsY *************************"
    Write(*,*)"********************************************************************"
    end if

  end subroutine CompCGibbsY !NoHanWorks20180425

  !*************************************
  !* Compute lesser Green's function *
  !*************************************
  subroutine CompCGibbsYHan(z,gless,GibbsY)
    use parameters, only: eta, biasvoltage, glue
    use constants, only: c_zero, ui, c_one, eleccharge, hbar
    use util, only: PrintCMatrix, PrintRMatrix
   !use lapack_blas, only: zgetri, zgetrf, zgemm
   !use lapack95, only: zgetri, zgetrf, zgemm
   !use lapack95
   !use blas95

    implicit none
    external               zgetri, zgetrf, zgemm

    integer :: i, j, info, error
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)

    complex*16, intent(in) :: z

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gless
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: GibbsY
    complex*16, dimension(NAOrbs,NAOrbs) :: CGibbsYHan
    complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

    ! 2018-04-25 CompCGibbsYHan DECLARATIONS
    complex*16, dimension(NAOrbs,NAOrbs) :: HGOp, OUT1, OUT2
    complex*16, dimension(NAOrbs,NAOrbs) :: HDPiv
    complex*16, dimension(NAOrbs*NAOrbs) :: HGOpCol, CGibbsYHanCol
    complex*16, dimension(NAOrbs*NAOrbs,NAOrbs*NAOrbs) :: LiouvSOpPiv
    !complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL
    !INTEGER MKL_GET_MAX_THREADS
    !INTEGER MAX_THREADS
    integer :: N, N2
    !DOUBLE PRECISION ALPHA, BETA
    real :: ALPHA, BETA



    if(DebugDev)then
    Write(*,*)"********************************************************************"
    Write(*,*)"********************** ENTERED CompCGibbsYHan **********************"
    Write(*,*)"********************************************************************"
    end if

    if(DebugDev)Write(*,'(A)')"ENTERED Device/CompCGibbsY"

     allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then
        print*,"Problems allocating"
        stop
     end if
     allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then
        print*,"Problems allocating"
        stop
     end if

    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD
    sigl=-ui*eta*SD

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl

    !************************************************************************
    !c Coupling matrices
    !************************************************************************
    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammal,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessL,NAOrbs)
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammar,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessR,NAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

     !GibbsYL=(Fermi-biasvoltage/2)*ui*(glessL*theta(dble(z)-(biasvoltage/2)))
     !GibbsYR=(Fermi+biasvoltage/2)*ui*(glessR*theta(dble(z)-(-biasvoltage/2)))
     ! REMOVE Fermi because everything is shifted to ZERO.
     !GibbsYL=(-biasvoltage/2)*ui*(glessL*theta(dble(z)-(biasvoltage/2)))
     !GibbsYR=(+biasvoltage/2)*ui*(glessR*theta(dble(z)-(-biasvoltage/2)))
     GibbsY= - ( (-biasvoltage/2)*(glessL*theta(dble(z)-(biasvoltage/2))) + (+biasvoltage/2)*(glessR*theta(dble(z)-(-biasvoltage/2))) )

     deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
     if (error /= 0 ) then
        print*,"Problems deallocating"
        stop
     end if

     if(DebugDev)then
       Write(*,'(A)')"Device/CompCGibbsYHan output gless(:,:) = "
       call PrintCMatrix(gless)
     end if



    ALPHA = 1.0d0
    BETA = 0.0d0
    N = NAOrbs

    do i=1,NAOrbs
      do j=1,NAOrbs
        HDPiv(i,j)=HD(iSpin,i,j)
      end do
    end do

    if(DebugDev)then
      Write(*,'(A)')"Device/CompCGibbsYHan output HDPiv(:,:) = "
      call PrintCMatrix(HDPiv)
    end if

    !MAX_THREADS = MKL_GET_MAX_THREADS()
    !CALL MKL_SET_NUM_THREADS(MAX_THREADS)

    ! pERFORM COMMUTATOR HD*glesser - glesser*HD
    !CALL ZGEMM('N','N',M,N,P,ALPHA,A,M,B,P,BETA,C,M)
    CALL ZGEMM ('N','N',N,N,N,c_one,HDPiv,N,gless,N,c_zero,OUT1,N ) ! 2018-04-25 REPLACED glesserBY gless BECAUSE ALL NaN WHEN PRINTING HGOp BELOW.
    CALL ZGEMM ('N','N',N,N,N,c_one,gless,N,HDPiv,N,c_zero,OUT2,N ) ! 2018-04-25 REPLACED glesserBY gless BECAUSE ALL NaN WHEN PRINTING HGOp BELOW.

    !HGOp = -(ui/hbar)*(OUT1 - OUT2)  ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.
    do i=1,NAOrbs
      do j=1,NAOrbs
        !HGOp(i,j) = -(ui/hbar)*(OUT1(i,j) - OUT2(i,j))  ! SEEMS TO BE THE CORRECT, BUT NON-HERMITIAN, MAYBE ANTIHERMITIAN.
        HGOp(i,j) = -(c_one/hbar)*(OUT1(i,j) - OUT2(i,j))  ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.
      end do
    end do

    !************************************************************************
    !c CONVERT TO SUPEROPERATOR FORM.
    !************************************************************************
    if(DebugDev)then
    if(iSpin==1)Write(*,'(A)')"ALPHA CURRENT OPERATOR"
    if(iSpin==2)Write(*,'(A)')"BETA CURRENT OPERATOR"
    Write(*,*)DREAL(HGOp(1,1))
    Write(*,*)DIMAG(HGOp(1,1))
    Write(*,*)DREAL(HGOp(1,2))
    Write(*,*)DIMAG(HGOp(1,2))
    call PrintCMatrix( HGOp )
    end if
    do i=1,NAOrbs
      do j=1,NAOrbs
        HGOpCol((i-1)*NAOrbs+j)=HGOp(i,j)
        !Write(*,'(I2,I2,A,F12.8,A,F12.8,A)'),i,j,"(",DREAL(HGOp(i,j)),")+i(",DIMAG(HGOp(i,j)),")"
      end do
      !Write(*, '(10(A,F12.8,A,F12.8,A))')( ("(",DREAL(HGOp(i,j)),") +i*(",DIMAG(HGOp(i,j)),")") ,j=1,NAOrbs)
    end do

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    N2 = NAOrbs*NAOrbs
    if(DebugDev)then
    if(iSpin==1)Write(*,'(A)')"ALPHA INVERSE LIOUVILLIAN"
    if(iSpin==2)Write(*,'(A)')"BETA INVERSE LIOUVILLIAN"
    end if

    if(.not.ALLOCATED(LiouvSop))then
      call BuildLiouvillian
    end if

    if(DebugDev)then
    Write(*,*)DREAL(LiouvSop(1,1,1))
    Write(*,*)DIMAG(LiouvSop(1,1,1))
    end if

    do i=1,N2
      do j=1,N2
        !Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)')iSpin,i,j,"(",DREAL(LiouvSOp(iSpin,i,j)),") + i*(",DIMAG(LiouvSOp(iSpin,i,j)),")"
        LiouvSOpPiv(i,j)=LiouvSOp(iSpin,i,j)
        !Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)')iSpin,i,j,"(",DREAL(LiouvSOpPiv(i,j)),")+i(",DIMAG(LiouvSOpPiv(i,j)),")"
      end do
    end do
    !alpha and beta are scalars,
    !A, B and C are matrices:
    !op(A) is an m-by-k matrix,
    !op(B) is a k-by-n matrix,
    !C is an m-by-n matrix.
    !call zgemm(transa, transb, m, n, k, alpha, A, lda, b, ldb, beta, c, ldc)
    ! THE 1 FOR THE 4th PARAMETER n IS THE NUMBER OF COLUMNS OF B AND C, 1 IN THIS CASE.
    call zgemm('N','N',N2,1,N2,c_one,LiouvSOpPiv,N2,HGOpCol,N2,c_zero,CGibbsYHanCol,N2)



    do i=1,NAOrbs
      do j=1,NAOrbs
        CGibbsYHan(i,j)=CGibbsYHanCol((i-1)*NAOrbs+j)
      enddo
      !Write(*, '(10(A,F12.8,A,F12.8,A))')( ("(",DREAL(CGibbsYHan(i,j)),") +i*(",DIMAG(CGibbsYHan(i,j)),")") ,j=1,NAOrbs)
    enddo


    if(DebugDev)then
    if(iSpin==1)Write(*,'(A)')"ALPHA CGibbsYHan(:,:)"
    if(iSpin==2)Write(*,'(A)')"BETA CGibbsYHan(:,:)"
    Write(*,*)DREAL(CGibbsYHan(1,1))
    Write(*,*)DIMAG(CGibbsYHan(1,1))
    Write(*,*)DREAL(CGibbsYHan(1,2))
    Write(*,*)DIMAG(CGibbsYHan(1,2))
    call PrintCMatrix( CGibbsYHan )
    end if

    if(DebugDev)then
    Write(*,*)"********************************************************************"
    Write(*,*)"************************ EXIT CompCGibbsYHan ***********************"
    Write(*,*)"********************************************************************"
    end if

  end subroutine CompCGibbsYHan

  !*************************************
  !* Compute lesser Green's function *
  !*************************************
  subroutine CompCGibbsYHanOld(glesser,CGibbsYHan)
    use parameters, only: eta, biasvoltage, glue
    use constants, only: c_zero, ui, c_one, eleccharge, hbar
   !use lapack_blas, only: zgetri, zgetrf, zgemm
   !use lapack95, only: zgetri, zgetrf, zgemm
   !use lapack95
   !use blas95

    implicit none
    external               zgemm

    integer :: i, j!, info, error
    integer, dimension(NAOrbs) :: ipiv
!    complex*16 :: work(4*NAOrbs)
!    complex*16, intent(in) :: z
    complex*16, dimension(NAOrbs,NAOrbs), intent(in) :: glesser
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: CGibbsYHan
    complex*16, dimension(NAOrbs,NAOrbs) :: HGOp, OUT1, OUT2
    complex*16, dimension(NAOrbs,NAOrbs) :: HDPiv
    complex*16, dimension(NAOrbs*NAOrbs) :: HGOpCol, CGibbsYHanCol
    complex*16, dimension(NAOrbs*NAOrbs,NAOrbs*NAOrbs) :: LiouvSOpPiv
    !complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL
    !INTEGER MKL_GET_MAX_THREADS
    !INTEGER MAX_THREADS
    integer :: N, N2
    !DOUBLE PRECISION ALPHA, BETA
    real :: ALPHA, BETA

    Write(*,*)"********************************************************************"
    Write(*,*)"********************** ENTERED CompCGibbsYHan **********************"
    Write(*,*)"********************************************************************"
    ALPHA = 1.0d0
    BETA = 0.0d0
    N = NAOrbs

    do i=1,NAOrbs
      do j=1,NAOrbs
        HDPiv(i,j)=HD(iSpin,i,j)
      end do
    end do

    !MAX_THREADS = MKL_GET_MAX_THREADS()
    !CALL MKL_SET_NUM_THREADS(MAX_THREADS)

    ! pERFORM COMMUTATOR HD*glesser - glesser*HD
    !CALL ZGEMM('N','N',M,N,P,ALPHA,A,M,B,P,BETA,C,M)
    CALL ZGEMM ('N','N',N,N,N,c_one,HDPiv,N,glesser,N,c_zero,OUT1,N )
    CALL ZGEMM ('N','N',N,N,N,c_one,glesser,N,HDPiv,N,c_zero,OUT2,N )

    !HGOp = -(ui/hbar)*(OUT1 - OUT2)  ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.
    do i=1,NAOrbs
      do j=1,NAOrbs
        !HGOp(i,j) = -(ui/hbar)*(OUT1(i,j) - OUT2(i,j))  ! SEEMS TO BE THE CORRECT, BUT NON-HERMITIAN, MAYBE ANTIHERMITIAN.
        HGOp(i,j) = -(c_one/hbar)*(OUT1(i,j) - OUT2(i,j))  ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.
      end do
    end do

    !************************************************************************
    !c CONVERT TO SUPEROPERATOR FORM.
    !************************************************************************
    if(iSpin==1)Write(*,'(A)')"ALPHA CURRENT OPERATOR"
    if(iSpin==2)Write(*,'(A)')"BETA CURRENT OPERATOR"
    Write(*,*)DREAL(HGOp(1,1))
    Write(*,*)DIMAG(HGOp(1,1))
    Write(*,*)DREAL(HGOp(1,2))
    Write(*,*)DIMAG(HGOp(1,2))
    call PrintCMatrix( HGOp )
    do i=1,NAOrbs
      do j=1,NAOrbs
        HGOpCol((i-1)*NAOrbs+j)=HGOp(i,j)
        !Write(*,'(I2,I2,A,F12.8,A,F12.8,A)'),i,j,"(",DREAL(HGOp(i,j)),")+i(",DIMAG(HGOp(i,j)),")"
      end do
      !Write(*, '(10(A,F12.8,A,F12.8,A))')( ("(",DREAL(HGOp(i,j)),") +i*(",DIMAG(HGOp(i,j)),")") ,j=1,NAOrbs)
    end do

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    N2 = NAOrbs*NAOrbs
    if(iSpin==1)Write(*,'(A)')"ALPHA INVERSE LIOUVILLIAN"
    if(iSpin==2)Write(*,'(A)')"BETA INVERSE LIOUVILLIAN"
    if(.not.ALLOCATED(LiouvSop))then
      call BuildLiouvillian
    end if
    Write(*,*)DREAL(LiouvSop(1,1,1))
    Write(*,*)DIMAG(LiouvSop(1,1,1))
    do i=1,N2
      do j=1,N2
        !Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)')iSpin,i,j,"(",DREAL(LiouvSOp(iSpin,i,j)),") + i*(",DIMAG(LiouvSOp(iSpin,i,j)),")"
        LiouvSOpPiv(i,j)=LiouvSOp(iSpin,i,j)
        !Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)')iSpin,i,j,"(",DREAL(LiouvSOpPiv(i,j)),")+i(",DIMAG(LiouvSOpPiv(i,j)),")"
      end do
    end do
    !alpha and beta are scalars,
    !A, B and C are matrices:
    !op(A) is an m-by-k matrix,
    !op(B) is a k-by-n matrix,
    !C is an m-by-n matrix.
    !call zgemm(transa, transb, m, n, k, alpha, A, lda, b, ldb, beta, c, ldc)
    ! THE 1 FOR THE 4th PARAMETER n IS THE NUMBER OF COLUMNS OF B AND C, 1 IN THIS CASE.
    call zgemm('N','N',N2,1,N2,c_one,LiouvSOpPiv,N2,HGOpCol,N2,c_zero,CGibbsYHanCol,N2)


    if(iSpin==1)Write(*,'(A)')"ALPHA CGibbsYHan(:,:)"
    if(iSpin==2)Write(*,'(A)')"BETA CGibbsYHan(:,:)"
    Write(*,*)DREAL(CGibbsYHan(1,1))
    Write(*,*)DIMAG(CGibbsYHan(1,1))
    Write(*,*)DREAL(CGibbsYHan(1,2))
    Write(*,*)DIMAG(CGibbsYHan(1,2))
    call PrintCMatrix( CGibbsYHan )
    do i=1,NAOrbs
      do j=1,NAOrbs
        CGibbsYHan(i,j)=CGibbsYHanCol((i-1)*NAOrbs+j)
      enddo
      !Write(*, '(10(A,F12.8,A,F12.8,A))')( ("(",DREAL(CGibbsYHan(i,j)),") +i*(",DIMAG(CGibbsYHan(i,j)),")") ,j=1,NAOrbs)
    enddo
    Write(*,*)"********************************************************************"
    Write(*,*)"************************ EXIT CompCGibbsYHan ***********************"
    Write(*,*)"********************************************************************"

  end subroutine CompCGibbsYHanOld


  ! *************
  ! Step function
  ! *************
  function theta(x)
    implicit none
    real,intent(in) :: x
    real :: theta
    theta=0.0d0
    if(x.lt.0.0d0) theta=1.0
    return
  end function theta

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  !****************************************
  double precision function f(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: chargeup, chargedown, rrr, a, b, Q
    integer :: i,j,M

    shift=x
    
    ! Radius of complex contour integration
    ! add 10eV just in case 
    rrr = 0.5*abs(EMin)+10.0d0;

    !c c Integral limits ... (a,b)
    a = 0.d0
    b = d_pi
    M=1000
    call intpj(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))
      
    ! Density matrix out of equilibirum
    if (biasvoltage /= 0.0) then
       M=100
       call intch(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
    !  print *,PDOUT
       do i=1,NAOrbs
          do j=1,NAOrbs
             PD(ispin,i,j)=PD(ispin,i,j)+PDOUT(ispin,i,j)
          end do
       end do
    end if
    ! Density matrix out of equilibirum
      
    Q=d_zero
    !do i=1,NAOrbs
    do i=NCDAO1, NCDAO2
       do j=1,NAOrbs
          Q=Q+PD(ispin,i,j)*SD(j,i)
       end do
    end do
    if (ispin.eq.1) QAlpha = Q
    if (ispin.eq.2) QBeta  = Q
 !if (ispin.eq.1) f = QAlpha - dble(GetNAE())
 !if (ispin.eq.2) f = QBeta  - dble(GetNBE())
    if (ispin.eq.1) f = QAlpha - dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
    if (ispin.eq.2) f = QBeta  - dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
    return
  end function f


  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  !****************************************
  double precision function QXTot(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    use cluster, only: NembedBL
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    if(DebugDev)Write(*,'(A)')"ENTERED Device/QxTot"
    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       call intpj(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=100
          if(NembedBL(1)==1)M=1000
          call intch(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          do i=1,NAOrbs
             do j=1,NAOrbs
                PD(ispin,i,j)=PD(ispin,i,j)+PDOUT(ispin,i,j)
             end do
          end do
       end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PD(ispin,i,j)*SD(j,i)
          end do
       end do
       if( ispin == 1 ) QAlpha = Q
       if( ispin == 2 ) QBeta  = Q
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QXTot = QAlpha + QBeta - dble(NCDEl) - QExcess
    if(DebugDev)Write(*,'(A,g16.6,A,g16.6)')"QXTot(",x,") = ",QXTot
    return
  end function QXTot

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  ! and integrates FOCK HW
  ! in parallel to PD for Pulay Forces.
  !****************************************
  double precision function QHWTot(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       call intpjHW(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=100
          call intchHW(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          do i=1,NAOrbs
             do j=1,NAOrbs
                PD(ispin,i,j)=PD(ispin,i,j)+PDOUT(ispin,i,j)
                HW(ispin,i,j)=HW(ispin,i,j)+HWOUT(ispin,i,j)
             end do
          end do
       end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PD(ispin,i,j)*SD(j,i)
          end do
       end do
       if( ispin == 1 ) QAlpha = Q
       if( ispin == 2 ) QBeta  = Q
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QHWTot = QAlpha + QBeta - dble(NCDEl) - QExcess
    return
  end function QHWTot

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  ! and integrates FOCK HW
  ! in parallel to PD for Pulay Forces.
  ! ADDED OM 2018-04-23 IN REPLACEMENT OF QCGibbsYTotOld to integrate glesserL/R from bottom instead bias window.
  !****************************************
  double precision function QCGibbsYTot(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    write(ifu_log,'(A51,i4)')' ENTERING QCGibbsYTot(x)'

    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case
       !rrr = 0.5*abs(EMin)+10.0d0; ! ORIGINAL CODE COMMENTED ON 2018-04.23
       rrr = 0.5*abs(EMin)+4.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       !call intpjHW(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0)) ! EL 2019-04-25 ME LO HE ENCONTRADO COMENTADO, PERO ESO NO ES COMPATIBLE CON EL RESULTADO TEÓRICO.

       ! Density matrix out of equilibirum
       !if (biasvoltage /= 0.0) then
          !M=100
          !call intchCGibbsY(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          !M=1000
          M=10000
          call intchCGibbsY(d_zero-dabs(biasvoltage/2.0)-rrr,dabs(biasvoltage/2.0),M)

!          THE DENSITY MATRIX PD WAS ALREADY FILLED, THIS IS ONLY QCGibbsYTot POSTPROCESSING.
          do i=1,NAOrbs
             do j=1,NAOrbs
                PDGIBBS(ispin,i,j)=PDGIBBS(ispin,i,j)+PDOUTGIBBS(ispin,i,j)
                !HW(ispin,i,j)=HW(ispin,i,j)+HWOUT(ispin,i,j)
             end do
          end do
          if(DebugDev)then
          if(iSpin==1)Write(*,'(A)')"WRITE ALPHA COMPLEX CGibbsY(1,:,:)"
          if(iSpin==2)Write(*,'(A)')"WRITE BETA COMPLEX CGibbsY(2,:,:)"
          Write(*,*)DREAL(CGibbsY(iSpin,1,1))
          Write(*,*)DIMAG(CGibbsY(iSpin,1,1))
          Write(*,*)DREAL(CGibbsY(iSpin,1,2))
          Write(*,*)DIMAG(CGibbsY(iSpin,1,2))
          call PrintCMatrix(CGibbsY(iSpin,:,:))

          if(iSpin==1)Write(*,'(A)')"WRITE ALPHA COMPLEX CGibbsYKernel1(1,:,:)"
          if(iSpin==2)Write(*,'(A)')"WRITE BETA COMPLEX CGibbsYKernel1(2,:,:)"
          Write(*,*)DREAL(CGibbsYKernel1(iSpin,1,1))
          Write(*,*)DIMAG(CGibbsYKernel1(iSpin,1,1))
          Write(*,*)DREAL(CGibbsYKernel1(iSpin,1,2))
          Write(*,*)DIMAG(CGibbsYKernel1(iSpin,1,2))
          call PrintCMatrix(CGibbsYKernel1(iSpin,:,:))

          if(iSpin==1)Write(*,'(A)')"WRITE ALPHA COMPLEX CGibbsYKernel2(1,:,:)"
          if(iSpin==2)Write(*,'(A)')"WRITE BETA COMPLEX CGibbsYKernel2(2,:,:)"
          Write(*,*)DREAL(CGibbsYKernel2(iSpin,1,1))
          Write(*,*)DIMAG(CGibbsYKernel2(iSpin,1,1))
          Write(*,*)DREAL(CGibbsYKernel2(iSpin,1,2))
          Write(*,*)DIMAG(CGibbsYKernel2(iSpin,1,2))
          call PrintCMatrix(CGibbsYKernel2(iSpin,:,:))

          end if
          do i=1,NAOrbs
             do j=1,NAOrbs
             !   Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)'),iSpin,i,j,"(",DREAL(CGibbsY(iSpin,i,j)),") + i*(",DIMAG(CGibbsY(iSpin,i,j)),")"
             end do
            !Write(*, '(10(A,F12.8,A,F12.8,A))')( ("(",DREAL(CGibbsY(iSpin,i,j)),") +i*(",DIMAG(CGibbsY(iSpin,i,j)),")") ,j=1,NAOrbs)
          end do
          if(MAXVAL(ABS(DIMAG(CGibbsY(iSpin,:,:)))) > 1.0D-3)then
            Write(*,'(A)')"***************************************************************************"
            Write(*,'(A,I1,A)')"******** CGibbsY(",iSpin,",:,:) CONTAINS NON-VANISHING IMAGINARY PART ********"
            Write(*,'(A)')"***************************************************************************"
          end if
          do i=1,NAOrbs
             do j=1,NAOrbs
                DGibbsY(ispin,i,j)=CGibbsY(ispin,i,j)
             end do
          end do
          if(MAXVAL(ABS(DIMAG(CGibbsYKernel1(iSpin,:,:)))) > 1.0D-3)then
            Write(*,'(A)')"***************************************************************************"
            Write(*,'(A,I1,A)')"***** CGibbsYKernel1(",iSpin,",:,:) CONTAINS NON-VANISHING IMAGINARY PART ****"
            Write(*,'(A)')"***************************************************************************"
          end if
          do i=1,NAOrbs
             do j=1,NAOrbs
                DGibbsYKernel1(ispin,i,j)=CGibbsYKernel1(ispin,i,j)
             end do
          end do
          if(MAXVAL(ABS(DIMAG(CGibbsYKernel2(iSpin,:,:)))) > 1.0D-3)then
            Write(*,'(A)')"***************************************************************************"
            Write(*,'(A,I1,A)')"***** CGibbsYKernel2(",iSpin,",:,:) CONTAINS NON-VANISHING IMAGINARY PART ****"
            Write(*,'(A)')"***************************************************************************"
          end if
          do i=1,NAOrbs
             do j=1,NAOrbs
                DGibbsYKernel2(ispin,i,j)=CGibbsYKernel2(ispin,i,j)
             end do
          end do
       !else
       !   do i=1,NAOrbs
       !      do j=1,NAOrbs
       !         DGibbsY(ispin,i,j)=(-shift)*PD(ispin,i,j)
       !      end do
       !   end do
       !end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PDGIBBS(ispin,i,j)*SD(j,i)
          end do
       end do

       ! COMMENT THIS BECAUSE DON'T WANT TO REPLACE THE CORRECT QAlpha AND  Qbeta.
       !if( ispin == 1 ) QAlpha = Q
       !if( ispin == 2 ) QBeta  = Q
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QCGibbsYTot = QAlpha + QBeta - dble(NCDEl) - QExcess
    return
  end function QCGibbsYTot

    !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  ! and integrates FOCK HW
  ! in parallel to PD for Pulay Forces.
  !****************************************
  double precision function QCGibbsYTotOld(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       call intpjHW(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=100
          call intchCGibbsY(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          do i=1,NAOrbs
             do j=1,NAOrbs
                PD(ispin,i,j)=PD(ispin,i,j)+PDOUT(ispin,i,j)
                HW(ispin,i,j)=HW(ispin,i,j)+HWOUT(ispin,i,j)
             end do
          end do
          if(DebugDev)then
          if(iSpin==1)Write(*,'(A)')"WRITE ALPHA COMPLEX CGibbsY(1,:,:)"
          if(iSpin==2)Write(*,'(A)')"WRITE BETA COMPLEX CGibbsY(2,:,:)"
          Write(*,*)DREAL(CGibbsY(iSpin,1,1))
          Write(*,*)DIMAG(CGibbsY(iSpin,1,1))
          Write(*,*)DREAL(CGibbsY(iSpin,1,2))
          Write(*,*)DIMAG(CGibbsY(iSpin,1,2))
          call PrintCMatrix(CGibbsY(iSpin,:,:))
          end if
          do i=1,NAOrbs
             do j=1,NAOrbs
             !   Write(*,'(I4,I4,I4,A,F12.8,A,F12.8,A)'),iSpin,i,j,"(",DREAL(CGibbsY(iSpin,i,j)),") + i*(",DIMAG(CGibbsY(iSpin,i,j)),")"
             end do
            !Write(*, '(10(A,F12.8,A,F12.8,A))')( ("(",DREAL(CGibbsY(iSpin,i,j)),") +i*(",DIMAG(CGibbsY(iSpin,i,j)),")") ,j=1,NAOrbs)
          end do
          if(MAXVAL(ABS(DIMAG(CGibbsY(iSpin,:,:)))) > 1.0D-3)then
            Write(*,'(A)')"***************************************************************************"
            Write(*,'(A,I1,A)')"******** CGibbsY(",iSpin,",:,:) CONTAINS NON-VANISHING IMAGINARY PART ********"
            Write(*,'(A)')"***************************************************************************"
          end if
          do i=1,NAOrbs
             do j=1,NAOrbs
                DGibbsY(ispin,i,j)=CGibbsY(ispin,i,j)
             end do
          end do
       end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PD(ispin,i,j)*SD(j,i)
          end do
       end do
       if( ispin == 1 ) QAlpha = Q
       if( ispin == 2 ) QBeta  = Q
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QCGibbsYTotOld = QAlpha + QBeta - dble(NCDEl) - QExcess
    return
  end function QCGibbsYTotOld

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  !****************************************
  double precision function QXTotEnergy(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    real :: QXTot, IntDOSE, intpjIntDOSE, intchIntDOSE, IntDOSEAlpha, IntDOSEBeta

    print*,"I am in QXTotEnergy" 

    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       !IntDOSE = d_zero
       !intchIntDOSE = d_zero
       !intpjIntDOSE = d_zero
       LocPD = 0.d0
       LocPDOUT = 0.d0

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       print*,"-shift = ", -shift
       !call intpjenergy(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))
       intpjIntDOSE=intpjenergy(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))
       !print*,"After call intpjIntDOSE = ", intpjIntDOSE
       ! Density matrix out of equilibirum
       !print*,"biasvoltage = ", biasvoltage
       if (biasvoltage /= 0.d0) then
          !print*,"Calculating bias window energy"
          M=100
          !call intchenergy(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          intchIntDOSE=intchenergy(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          do i=1,NAOrbs
             do j=1,NAOrbs
                LocPD(ispin,i,j)=LocPD(ispin,i,j)+LocPDOUT(ispin,i,j)
             end do
          end do
       end if
       ! Density matrix out of equilibirum

       !intchIntDOSE = 0.d0
       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+LocPD(ispin,i,j)*SD(j,i)
             ! WITH THIS, I TAKE INTO ACCOUNT INTPJ (EQ.) AND INTCH (NONEQ.)
             !intchIntDOSE=intchIntDOSE+(d_zero-x)*PD(ispin,i,j)*SD(j,i)
             !intchIntDOSE=intchIntDOSE+(x-shift)*PD(ispin,i,j)*SD(j,i)
             !intchIntDOSE=intchIntDOSE+x*PD(ispin,i,j)*SD(j,i)
          end do
       end do
       print*,"intpjenergy IntDOSE = ", intpjIntDOSE
       print*,"intchenergy IntDOSE = ", intchIntDOSE
       IntDOSE = intpjIntDOSE + intchIntDOSE
       print*,"IntDOSE = ", IntDOSE
       if( ispin == 1 ) QAlpha = Q
       if( ispin == 2 ) QBeta  = Q
       if( ispin == 1 ) IntDOSEAlpha = IntDOSE
       if( ispin == 2 ) IntDOSEBeta  = IntDOSE
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QXTot = QAlpha + QBeta - dble(NCDEl) - QExcess
    if( NSpin == 1 ) IntDOSEBeta = IntDOSEAlpha
    QXTotEnergy = IntDOSEAlpha + IntDOSEBeta
    return
  end function QXTotEnergy

  !****************************************
  ! Function gives excess charge of selected atom
  ! in dependence of Fermi energy x
  !****************************************
  double precision function FPart(x)
    use cluster, only: NAOAtom, LoAOrbNo, HiAOrbNo
    use parameters, only: biasvoltage,QExcess, SPIN_Beg, SPIN_End
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: chargeup, chargedown, rrr, a, b, Q
    integer :: i,j,M

    !write(ifu_log,*)'-------------------------------------'
    !write(ifu_log,*)'------- I am in FPart ---------------'
    !write(ifu_log,*)'-------------------------------------'

    shift=x
    
    ! Radius of complex contour integration
    ! add 10eV just in case 
    rrr = 0.5*abs(EMin)+10.0d0;

    !c c Integral limits ... (a,b)
    a = 0.d0
    b = d_pi
    M=1000
    call localintpj(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))
      
    ! Density matrix out of equilibirum
    if (biasvoltage /= 0.0) then
       M=100
       call localintch(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
    !  print *,PDOUT
       !do i=1,NAOrbs
       do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom)
          do j=1,NAOrbs
             LocPD(ispin,i,j)=LocPD(ispin,i,j)+LocPDOUT(ispin,i,j)
          end do
       end do
    end if
    ! Density matrix out of equilibirum
      
    Q=d_zero
    !do i=1,NAOrbs
    !do i=NCDAO1, NCDAO2
    do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom)
       do j=1,NAOrbs
          Q=Q+LocPD(ispin,i,j)*SD(j,i)
       end do
    end do
    if (ispin.eq.1) LocalQAlpha = Q
    if (ispin.eq.2) LocalQBeta  = Q
!    !if (ispin.eq.1) f = QAlpha - dble(GetNAE())
!    !if (ispin.eq.2) f = QBeta  - dble(GetNBE())
!    !if (ispin.eq.1) f = QAlpha - dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
!    !if (ispin.eq.2) f = QBeta  - dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
!    if (ispin.eq.1) FPart = LocalQAlpha + alphaelec(spinatom-SPIN_Beg+1) !- GetAtmChg(spinatom)
!    if (ispin.eq.1) FPart = LocalQBeta + betaelec(spinatom-SPIN_Beg+1) !- GetAtmChg(spinatom)


    if (ispin.eq.1) FPart = LocalQAlpha - alphaelec(spinatom-SPIN_Beg+1) ! - GetAtmChg(spinatom)
    if (ispin.eq.2) FPart = LocalQBeta - betaelec(spinatom-SPIN_Beg+1) ! - GetAtmChg(spinatom)
    !QXPart = LocalQAlpha + LocalQBeta
    write(ifu_log,'(A,I4,A,f9.5)')' Curr. spinatom:',spinatom,'		FPart: ',FPart
    return
  end function FPart

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  !****************************************
  double precision function QXPart(x)
    use cluster, only: NAOAtom, LoAOrbNo, HiAOrbNo
    use parameters, only: biasvoltage,QExcess, SPIN_Beg, SPIN_End
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE, GetNAtoms, GetAtmChg
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    !real :: LocalQAlpha, LocalQBeta
    !integer, intent(in) :: iatom

    !print*,"I am in QXTot" 
    !write(ifu_log,*)'-------------------------------------'
    !write(ifu_log,*)'------- I am in QXPart ---------------'
    !write(ifu_log,*)'-------------------------------------'
    
    !do ispin=1,NSpin
    !  do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom)
    !    do j=1,NAOrbs
    !      LocPD(ispin,i,j)=PD(ispin,i,j)
    !    end do
    !  end do
    !end do

    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       call localintpj(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.d0) then
          M=100
          call localintch(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          do i=1,NAOrbs
          !do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom)
             do j=1,NAOrbs
                LocPD(ispin,i,j)=LocPD(ispin,i,j)+LocPDOUT(ispin,i,j)
             end do
          end do
       end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       !do i=NCDAO1, NCDAO2
       !do i=LoAOrbNo(SPIN_Beg), HiAOrbNo(SPIN_End)
       do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom)
          do j=1,NAOrbs
             Q=Q+LocPD(ispin,i,j)*SD(j,i)
          end do
       end do
       !if( ispin == 1 ) QAlpha = Q
       !if( ispin == 2 ) QBeta  = Q
       if( ispin == 1 ) LocalQAlpha = Q
       if( ispin == 2 ) LocalQBeta  = Q
    end do

    !if( NSpin == 1 ) QBeta = QAlpha
    if( NSpin == 1 ) LocalQBeta = LocalQAlpha
!    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
!    !QXTot = QAlpha + QBeta - dble(NCDEl) - QExcess
!    write(ifu_log,*)'Current spinatom: ', spinatom
!    write(ifu_log,*)' Alpha Local Shift= ',x
!    write(ifu_log,*)'LocalQAlpha', LocalQAlpha
!    write(ifu_log,*)'LocalQBeta', LocalQBeta
!    !write(ifu_log,*)'LocalQAlpha + LocalQBeta', LocalQAlpha + LocalQBeta
!    write(ifu_log,*)'alphaelec', alphaelec(spinatom-SPIN_Beg+1)
!    write(ifu_log,*)'betaelec', betaelec(spinatom-SPIN_Beg+1)
!    !write(ifu_log,*)'alphaelec + betaelec', alphaelec + betaelec
!    write(ifu_log,*)'QXPart = ', LocalQAlpha + LocalQBeta - alphaelec(spinatom-SPIN_Beg+1) - betaelec(spinatom-SPIN_Beg+1)! - GetAtmChg(spinatom)

    QXPart = LocalQAlpha + LocalQBeta - alphaelec(spinatom-SPIN_Beg+1) - betaelec(spinatom-SPIN_Beg+1)! - GetAtmChg(spinatom)
    !QXPart = LocalQAlpha + LocalQBeta
    write(ifu_log,'(A,I4,A,f9.5)')' Curr. spinatom:',spinatom,'		QXPart: ',QXPart
    return
  end function QXPart


  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x with SOC
  !****************************************
  double precision function QXTot_SOC(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real, intent(in) :: x

    real :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000
       call intpj_SOC(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
      !if (biasvoltage /= 0.0) then
      !   M=100
      !   call intch_SOC(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
      !   PD_SOC=PD_SOC+PDOUT_SOC
      !end if
       ! Density matrix out of equilibirum

       do i=NCDAO1, NCDAO2
          do j=1,DNAOrbs
             Q=Q+PD_SOC(i,j)*S_SOC(j,i)
          end do
       end do
       do i=NCDAO1+NAOrbs, NCDAO2+NAOrbs
          do j=1,DNAOrbs
             Q=Q+PD_SOC(i,j)*S_SOC(j,i)
          end do
       end do

    Q_SOC = Q
    QXTot_SOC = Q - dble(NCDEl) - QExcess
    return
  end function QXTot_SOC


  !ccccccccccccccccccccccccccccccc
  !c                                                                              c
  !c     Change of variable no. 3 for the numerical integration:                  c
  !c                                                                              c
  !c     int_{-infty}^{Eq} DOS(E)dE = int_{-1}^{1} DOS(E(x))*(dE/dx)dx            c
  !c                                                                              c
  !c     E = Em*(1-bx)/(1+x)                                                      c
  !c                                                                              c
  !ccccccccccccccccccccccccccccccc
  real function edex3(Em,b,x)
    implicit none
    real :: Em, b ,x
    edex3 = 0.5d0*((Em-b)*x + (Em+b))
    return
  end function edex3
  
  !******************************!
  ! Writing the Hamiltonian      !
  !******************************!
  subroutine Hamiltonian
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, NEmbedBL
    USE parameters, only: Hamilton
    use g09Common, only: GetAtmCo
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n
    real :: ro_a, ro_b

    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Self-consistent Hamiltonian  ---'
    write(ifu_log,*)'-------------------------------------'

    
    !open(ifu_ham,file='V.'//trim(xxx)//'.dat',status='unknown')
    write(ifu_log,*)'---------'
    write(ifu_log,*)'Left lead'
    write(ifu_log,*)'---------'
    I1=0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+hd(1,i,i)
          if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       end do
       ro_a=ro_a/NAOAtom(j)
       ro_b=ro_b/NAOAtom(j)
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j, ro_a
       if(NSpin ==2 ) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
       if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
       END IF
       I1 = I1 + NAOAtom(j)
    end do
    
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'--------'
       write(ifu_log,*)'Molecule'
       write(ifu_log,*)'--------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+hd(1,i,i)
             if(NSpin==2) ro_b=ro_b+hd(2,i,i)
          end do
          ro_a=ro_a/NAOAtom(j)
          ro_b=ro_b/NAOAtom(j)
          if(NSpin==1) write(ifu_log,1011)'Atom:',j, ro_a
          if(NSpin==2) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
          IF (Hamilton) THEN
          if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
          if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
          END IF
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
    write(ifu_log,*)'----------'
    write(ifu_log,*)'Right lead'
    write(ifu_log,*)'----------'
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+hd(1,i,i)
          if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       end do
       ro_a=ro_a/NAOAtom(j)
       ro_b=ro_b/NAOAtom(j)
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j, ro_a
       if(NSpin==2) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
       if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
       END IF
    end do
1012 format(4f10.5)
1011 format(a6,i4,f12.6)
     !close(ifu_ham)
  end subroutine Hamiltonian

  !*******************************************!
  ! Writing the Diagonalized Hamiltonian      !
  !*******************************************!
  subroutine DiagHamiltonian
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, NEmbedBL
    USE parameters, only: Hamilton
    use g09Common, only: GetAtmCo
    use constants, only: Bohr
    implicit none

    real, dimension(:,:),allocatable :: kseigenvalues
    integer :: i,j, I1, is ,n
    real :: ro_a, ro_b

    write(ifu_log,*)'--------------------------------------------------'
    write(ifu_log,*)'---  Self-consistent Diagonalized Hamiltonian  ---'
    write(ifu_log,*)'--------------------------------------------------'

    
    !open(ifu_diagham,file='V.'//trim(xxx)//'.dat',status='unknown')
    !write(ifu_log,*)'---------'
    !write(ifu_log,*)'Left lead'
    !write(ifu_log,*)'---------'
    !I1=0
    do j=1,NAOrbs
       !ro_a=0.0d0
       !ro_b=0.0d0
       !do i=I1+1,I1+NAOAtom(j)
       !   ro_a=ro_a+hd(1,i,i)
       !   if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       !end do
       !ro_a=ro_a/NAOAtom(j)
       !ro_b=ro_b/NAOAtom(j)
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j, kseigenvalues(1,j)
       if(NSpin ==2 ) write(ifu_log,1011)'Atom:',j, kseigenvalues(2,j)
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_diagham,1012)j,kseigenvalues(1,j),(hd(1,i,j),i=1,NAOrbs)
       if(NSpin ==2 ) write(ifu_diagham,1012)j,kseigenvalues(2,j),(hd(2,i,j),i=1,NAOrbs)
       END IF
       !I1 = I1 + NAOAtom(j)
    end do

1012 format(4f10.5)
1011 format(a6,i4,f12.6)
     !close(ifu_diagham)
  end subroutine DiagHamiltonian


  !******************************!
  ! Mulliken population analysis !
  !******************************!
  subroutine MullPop
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol
    USE parameters, only: Mulliken, LDOS_Beg, LDOS_End
    use g09Common, only: GetAtmCo
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n, l
    real :: sdeg, ro_a, ro_b, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b, tmp

    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Mulliken population analysis ---'
    write(ifu_log,*)'-------------------------------------'


    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    rho_a = matmul( PD(1,:,:), SD )
    if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 1'
    write(ifu_log,*)'----------------------'
    I1=0
    chargelead1=0.0
    spinlead1=0.0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin ==2 ) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead1=chargelead1+(ro_a+ro_b)
       if (NSpin == 2) spinlead1=spinlead1+(ro_a-ro_b)
       I1 = I1 + NAOAtom(j)
    end do
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 1:',chargelead1*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 1:',spinlead1
   write(ifu_log,*)'----------------------------------------------------'

    chargemol=0.0d0
    spinmol=0.0d0
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'-------------------'
       write(ifu_log,*)'Charges in molecule'
       write(ifu_log,*)'-------------------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+rho_a(i,i)
             if(NSpin==2) ro_b=ro_b+rho_b(i,i)
          end do
          if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',sdeg*ro_a
          if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
          IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
          END IF
          IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
          END IF
          chargemol=chargemol+ro_a+ro_b
          if (NSpin == 2) spinmol=spinmol+(ro_a-ro_b)
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in molecule:',chargemol*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in molecule:',spinmol
   write(ifu_log,*)'----------------------------------------------------'

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 2'
    write(ifu_log,*)'----------------------'
    chargelead2=0.0
    spinlead2=0.0
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead2=chargelead2+ro_a+ro_b
       if (NSpin == 2) spinlead2=spinlead2+(ro_a-ro_b)
    end do

   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 2:',chargelead2*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 2:',spinlead2
   write(ifu_log,*)'----------------------------------------------------'

1011 format(a6,i4,a10,f8.4)
1012 format(a6,i4,a10,f8.4,a10,2f9.4)
1013 format(7f9.4)
  end subroutine MullPop

!*******************************************************************************
!* Subroutine to compute LDOS(E) when it is not computed by Transmission        *
!*******************************************************************************
  SUBROUTINE LDOS
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: c_one, c_zero, d_zero, d_pi
    use parameters, only: LDOS_Beg, LDOS_End, EW1, EW2, EStep, DOSEnergy
    use preproc, only: MaxAtm
    use omp_lib
    use g09Common, only: GetNAtoms

    implicit none
    real :: energy, trans, DOS, energ
    real, dimension(10001) :: xxx
    real, dimension(MaxAtm) :: AtomDOS
    complex*16 :: cenergy,ctrans
    integer :: n, nsteps, i, imin, imax, info ,j
    
    complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG

    print *
    print *, "-------------------------"
    print *, "--- Calculating  LDOS ---"
    print *, "-------------------------"
    print *

    allocate ( AtomDOSEF(2,MaxAtm) )

    nsteps = (EW2-EW1)/EStep + 1
    do ispin=1,NSpin

       open(333,file='tempDOS',status='unknown')
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,green,gammar,gammal)
!!$OMP DO SCHEDULE(DYNAMIC)
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(cenergy,Green,GammaR,GammaL)

          ! Mulliken DOS 
!!$OMP CRITICAL
          SG = matmul( SD, green )
          ! computing total DOS
          DOS=d_zero
          AtomDOS=d_zero
          do j=1,GetNAtoms()
          do i=LoAOrbNo(j),HiAOrbNo(j)
             AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
             DOS=DOS-dimag(SG(i,i))/d_pi
          end do
          end do

          if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

          ! print out DOS and atomic orbital resolved DOS ***
          imin = LoAOrbNo(LDOS_Beg)
          if( imin < 1 ) imin = 1
          imax = HiAOrbNo(LDOS_End)
          if( imax > NAOrbs ) imax = NAOrbs
          call flush(333)
          write(333,3333) energy,DOS*(-1)**(ispin+1),(AtomDOS(j)*(-1)**(ispin+1),j=LDOS_Beg,LDOS_End),(-dimag(SG(i,i))*(-1)**(ispin+1)/d_pi,i=imin,imax)

!!$OMP END CRITICAL

       end do ! End of energy loop
      
!!$OMP END DO
!!$OMP END PARALLEL

  ! Reordering in energy for nice output
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
 
      close(333,status='delete')

        write(ifu_dos,*) '    '                    
      end do ! End of spin loop


3333 format(f10.5,10000E14.5)

  END SUBROUTINE LDOS

!*******************************
!* Subroutine to evaluate T(E) *
!*******************************
  subroutine transmission
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: c_one, c_zero, d_zero, d_pi, ui
    use parameters, only: NChannels,HTransm,EW1,EW2,EStep,LDOS_Beg,LDOS_End, DOSEnergy, soc, FermiAcc, QExcess, ChargeAcc, eta
    use numeric, only: CMatPow, CHDiag, CDiag, sort, MULLER, MULLER_OMP
   !use lapack_blas, only: zgemm
   !use lapack95
   !use blas95
    use preproc, only: MaxAtm
    use omp_lib
    use g09Common, only: GetNAtoms

    implicit none
    external               zgemm
    real :: energy, trans, DOS, energ, E0, E1, E2, E3, Delta, Epsilon, DE
    real, dimension(MaxAtm) :: AtomDOS
    real, dimension(10001) :: xxx
    complex*16 :: cenergy,ctrans
    integer :: n, nsteps, i, imin, imax, info, j, AllocErr, cond, k
    integer :: Max = 20
    complex*16, dimension(:,:), allocatable :: GammaL, GammaR, Green, T, temp, SG
    complex*16, dimension(:,:), allocatable :: DGammaL, DGammaR, DGreen, DT, Dtemp, DSG
    complex*16, dimension(:,:),allocatable :: dummy
    real, dimension(:), allocatable :: tn,Dtn
    complex*16, dimension(:), allocatable :: ctn,Dctn
    real, dimension(:),allocatable   :: tchan1,tchan2

    print *
    print *, "--------------------------------"
    print *, "--- Calculating Transmission ---"
    print *, "---  (and DOS if required)   ---"
    print *, "--------------------------------"
    print *

    if (soc /= 0.0) then
       write(ifu_log,*)' Adding spin-orbit coupling ...'
       write(ifu_log,*)'... and finding new Fermi level'
       allocate(DGammaL(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DGammaR(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DGreen(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DT(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Dtemp(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DSG(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(S_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(H_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(PD_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Dtn(DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Dctn(DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       DSG=c_zero
       DT=c_zero
       Dtemp=c_zero
       S_SOC=d_zero
       do i=1,NAOrbs
       do j=1,NAOrbs
          S_SOC(i,j)=SD(i,j)
          S_SOC(i+NAOrbs,j+NAOrbs)=SD(i,j)
       end do
       end do
       call spin_orbit
       
    else
       allocate(GammaL(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(GammaR(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Green(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(T(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(temp(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(SG(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(tn(NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(ctn(NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
    end if

    if( NChannels > 0 ) then
       print *, "Number of eigen-channels to print out: ", NChannels 
       allocate( tchan1(NChannels), tchan2(NChannels), dummy(NChannels,NChannels), STAT = AllocErr )
       if( AllocErr /= 0 ) stop
    end if

    allocate( AtomDOSEF(2,MaxAtm), STAT=AllocErr)
    if( AllocErr /= 0 ) stop

    ! finding new Fermi energy with SOC
    if (soc /= 0.0) then
      
       E0=shift-10.0d0*FermiAcc
       E1=shift
       E2=shift+10.0d0*FermiAcc
       Delta=FermiAcc
       Epsilon=ChargeAcc*(NCDEl+QExcess)
       print*,'MULLER method'
       call MULLER(QXTot_SOC,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
       !call MULLER_OMP(QXTot_SOC,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
       if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
         write(ifu_log,*) 'I could not accurately find the new Fermi level ...'
         write(ifu_log,*) ' ...using the best approximation'
       end if
       shift = E3
       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of electrons:  ', Q_SOC
       write(ifu_log,*)'-----------------------------------------------'
    end if
         
    nsteps = (EW2-EW1)/EStep + 1

    do ispin=1,NSpin

      if (ispin == 2 .and. soc /= 0.0) exit

      if (LDOS_Beg <= LDOS_End ) open(333,file='tempDOS',status='unknown')
      open(334,file='tempT',status='unknown')

      if (soc == 0.0d0) then
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,green,gammar,gammal,T,temp)
!!$OMP DO SCHEDULE(static,1)
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          !cenergy=dcmplx(energy) ! ORIGINAL CODE BUT SEEMS TO FAIL IN LAST CYCLE FOR CERTAIN SYSTEMS: Li2 DIMER.
          cenergy=dcmplx(energy+ui*eta)
          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          !Write(ifu_log,'(A,F8.4,A,F8.4,A)')"call gplus(",DREAL(cenergy),"+i",DIMAG(cenergy),",Green,GammaR,GammaL)"
          call gplus(cenergy,Green,GammaR,GammaL)

          if( .not. HTransm )then
             !*************************************************************
             !* Here we use the following non-Hermitian expression  for T *
             !* [Gamma_L G^a Gamma_R G^r]                                 *
             !* It works better for large clusters                        *
             !*************************************************************
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL,NAorbs, Green,  NAOrbs, c_zero, T,    NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green,  NAOrbs, c_zero, T,    NAOrbs)
          else
             !********************************************************
             !* Here we use the following Hermitian expression for T *
             !* [Gamma_L^1/2 G^a Gamma_R G^r Gamma_L^1/2]            *
             !********************************************************
             call CMatPow(GammaL,0.5d0,temp)
             GammaL=temp
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,GammaL,NAOrbs,Green, NAOrbs,c_zero,temp,NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,temp,  NAOrbs,GammaR,NAOrbs,c_zero,T,   NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,T,     NAOrbs,Green, NAOrbs,c_zero,temp,NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,temp,  NAOrbs,GammaL,NAOrbs,c_zero,T,   NAOrbs)
          end if

!!$OMP CRITICAL
          ! Mulliken DOS 
          if (LDOS_Beg <= LDOS_End ) then
            SG = matmul( SD, green )
            ! computing total DOS
            DOS=d_zero
            AtomDOS=d_zero
            do j=1,GetNAtoms()
            do i=LoAOrbNo(j),HiAOrbNo(j)
               AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
               DOS=DOS-dimag(SG(i,i))/d_pi
            end do
            end do

            if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

            ! print out DOS and atomic orbital resolved DOS ***
            imin = LoAOrbNo(LDOS_Beg)
            if( imin < 1 ) imin = 1
            imax = HiAOrbNo(LDOS_End)
            if( imax > NAOrbs ) imax = NAOrbs
            call flush(333)
            write(333,3333) energy,DOS*(-1)**(ispin+1),(AtomDOS(j)*(-1)**(ispin+1),j=LDOS_Beg,LDOS_End),(-dimag(SG(i,i))*(-1)**(ispin+1)/d_pi,i=imin,imax)
          end if

          ! computing transmission T
          ctrans=c_zero
          do i=1,NAOrbs
             ctrans=ctrans + T(i,i)
          end do

          if (dimag(ctrans).gt.1.0d-5) then
             write(ifu_log,*)'Transmission not real !!!'
             stop
          end if
          trans=ctrans

          ! Diagonalize the T matrix 
          ! to get eigen channels
          if( NChannels > 0 )then
             if( HTransm ) then 
                call CHDiag( T, tn, info )
             else
                call CDiag( T, ctn, info )
                do i=1,NAOrbs
                  tn(i) = dble( ctn(i) )
                end do
                ! sort eigenvalues smallest to biggest
                call sort(NAOrbs,tn)
             end if
             if( n > 3 ) call SeparateSpaghettis( tchan1, tchan2, tn(NAOrbs-NChannels+1:NAOrbs), dummy, NChannels)
             tchan1=tchan2
             tchan2=tn(NAOrbs-NChannels+1:NAOrbs)
          end if
          call flush(334)
          write(334,1002)energy,trans,(tn(i),i=NAOrbs,NAOrbs-NChannels+1,-1)
          
!!$OMP END CRITICAL
       end do ! End of energy loop
!!$OMP END DO
!!$OMP END PARALLEL

       else !SOC case

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,Dgreen,Dgammar,Dgammal,DT,Dtemp) 
!!$OMP DO SCHEDULE(STATIC,10)
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
             call gplus_SOC(cenergy,DGreen,DGammaR,DGammaL)
             call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one, DGammaL,DNAOrbs, DGreen,  DNAOrbs, c_zero, DT,    DNAOrbs)
             call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one, DT,     DNAOrbs, DGammaR, DNAOrbs, c_zero, Dtemp, DNAOrbs)
             call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one, Dtemp,  DNAOrbs, DGreen,  DNAOrbs, c_zero, DT,    DNAOrbs)

!!$OMP CRITICAL
      ! Mulliken DOS 
           if (LDOS_Beg <= LDOS_End ) then
             DSG = matmul( S_SOC, DGreen )
             DOS=d_zero
             AtomDOS=d_zero
             do j=1,GetNAtoms()
             do i=LoAOrbNo(j),HiAOrbNo(j)
                AtomDOS(j)=AtomDOS(j)-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs))/(2*d_pi)
                DOS=DOS-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs))/(2*d_pi)
             end do
             end do

             if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

     ! print out DOS and atomic orbital resolved DOS ***
             imin = LoAOrbNo(LDOS_Beg)
             if( imin < 1 ) imin = 1
             imax = HiAOrbNo(LDOS_End)
             if( imax > NAOrbs ) imax = NAOrbs
             call flush(333)
             write(333,3333) energy,DOS,(AtomDOS(j),j=LDOS_Beg,LDOS_End),((-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs)))/(2*d_pi),i=imin,imax)
           end if

          ! computing transmission T
          ctrans=c_zero
          do i=1,DNAOrbs
             ctrans=ctrans + DT(i,i)
          end do

          if (dimag(ctrans).gt.1.0d-5) then
             write(ifu_log,*)'Transmission not real !!!'
             stop
          end if
          trans=ctrans/2.0   ! in units of 2e^2/h

          ! Diagonalize the T matrix 
          ! to get eigen channels
          if( NChannels > 0 )then
             if( HTransm ) then 
                call CHDiag( DT, Dtn, info )
             else
                call CDiag( DT, Dctn, info )
                do i=1,DNAOrbs
                  Dtn(i) = dble( Dctn(i) )
                end do
                ! sort eigenvalues smallest to biggest
                call sort(DNAOrbs,Dtn)
             end if
             if( n > 3 ) call SeparateSpaghettis( tchan1, tchan2, Dtn(DNAOrbs-NChannels+1:DNAOrbs), dummy, NChannels)
             tchan1=tchan2
             tchan2=tn(DNAOrbs-NChannels+1:DNAOrbs)
          end if

          call flush(334)
          !write(334,1002)energy,trans,(Dtn(i),i=DNAOrbs,DNAOrbs-NChannels+1,-1)
          write(334,1002)energy,trans

!!$OMP END CRITICAL
       end do ! End of energy loop
!!$OMP END DO
!!$OMP END PARALLEL

       end if

       if (LDOS_Beg <= LDOS_End ) then
  ! Reordering in energy for nice DOS output
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
      write(ifu_dos,*)'   '
      close(333,status='delete')
      end if

  ! Reordering in energy for nice T output
      do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(334)
          do i=1,10000000000
          read(334,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(334)
             read(334,1002) (xxx(j),j=1,2+NChannels)
             write(ifu_tra,1002) (xxx(j),j=1,2+NChannels)
             exit
          end if
          end do
       end do
      close(334,status='delete')

      write(ifu_tra,*)'   '

      end do ! End of spin loop

      if (soc /= 0.0) then
         deallocate(DGammaL)
         deallocate(DGammaR)
         deallocate(DGreen)
         deallocate(DT)
         deallocate(Dtemp)
         deallocate(DSG)
         deallocate(Dtn)
         deallocate(Dctn)
         deallocate(S_SOC)
         deallocate(H_SOC)
      else 
         deallocate(GammaL)
         deallocate(GammaR)
         deallocate(Green)
         deallocate(T)
         deallocate(temp)
         deallocate(SG)
         deallocate(tn)
         deallocate(ctn)
      end if

    if( NChannels > 0 ) then
       deallocate( tchan1, tchan2, dummy )
    end if

1002 format(f10.5,10000E14.5)
3333 format(f10.5,10000E14.5)

  end subroutine transmission

  !
  ! Computes Hybridization functions for all correlated subspaces
  ! for all energy points defined in mesh.dat and writes them to 
  ! file
  !
  subroutine CompHybFunc
    use constants
    use parameters, only: NCorrBl, CorrBeg, CorrEnd, eta
    use util
    use correlation
    use antcommon
    use omp_lib
    implicit none

    integer :: ios, n, nmax, iblock, i, j, iao,jao
    real :: En, Ep
    complex*16 :: zn
    character(len=3) :: istr
    character(len=100) :: fname

    real, dimension(:), allocatable :: EMesh
    
    ! Device Green's function
    complex*16, dimension(NAorbs,NAOrbs) :: GD

    complex*16, dimension(:,:,:,:,:), allocatable :: delta

    print *, "-----------------------------------------"
    print *, "--- Computing Hybridization functions ---"
    print *, "-----------------------------------------"

    call SetHamOvl( HD, SD )
    !
    ! Read mesh file mesh.dat
    ! 
    print *, "Reading energy mesh from file mesh.dat" 
    ! open mesh file
    open(unit=ifu_msh,file='mesh.dat',status='old',iostat=ios)
    if( ios /= 0 )then
       print *, "Device/CompHybFunc/Error: Could not open energy mesh file mesh.dat. Abort."
       STOP
    end if
    ios = 0; nmax=0; Ep=-1e+10
    do while( ios == 0 )
       read(unit=ifu_msh,fmt=*,iostat=ios), En
       if( ios /= 0 ) exit
       if( En <= Ep ) exit
       Ep = En
       nmax=nmax+1
       !print *, nmax, En
    end do
    print *, "Mesh file has", nmax, " data points."
    allocate( EMesh( nmax ) )
    rewind(ifu_msh)
    do n=1,nmax
       read(unit=ifu_msh,fmt=*,iostat=ios), EMesh(n)
       if( ios /= 0 ) exit
    end do
    close(ifu_msh)
    !
    ! calculate hybridization function for all mesh points
    ! 
    allocate( delta(nmax,NSpin,NCorrBl,NMaxCorr,NMaxCorr) ) 
    delta = c_zero

    do ispin=1,NSpin
!$OMP PARALLEL PRIVATE(En,zn,GD)
!$OMP DO
       do n=1,nmax
          En = EMesh(n)
          zn = En+ui*eta
          call gplus0(zn,GD)
!$OMP CRITICAL
          call CompDelta( ispin, En, -shift, GD, delta(n,ispin,:,:,:) )
!$OMP END CRITICAL
       end do
!$OMP END DO
!$OMP END PARALLEL
    end do
    !
    ! write hybridization functions to files
    !
    do iblock=1,NCorrBl
       print *, "Correlated block", iblock
       call int2str( iblock, istr )
       fname='delta.out.'//istr
       print *, "Output file for Hybridization:", fname
       open(unit=ifu_hyb,file=fname,status='unknown',iostat=ios)
       fname='Ac.out.'//istr 
       print *, "Output file for Bath Sepctral function:", fname
       open(unit=ifu_ac,file=fname,status='unknown',iostat=ios)
       do n=1,nmax
          En=EMesh(n)
          !
          ! write hybridization function Delta
          !
          write(ifu_hyb,fmt='(E20.10,1000E14.6)'),& 
               En, ( (delta(n,ispin,iblock,i,i),ispin=1,NDSpin),i=1,ncorrao(iblock) )
          call flush(ifu_hyb)
          !
          ! write bath spectral function = -Im(Delta)/pi
          !
          write(ifu_ac,fmt='(E20.10,1000E14.6)'), &
               En, ( (-AIMAG(delta(n,ispin,iblock,i,i))/d_pi,ispin=1,NDSpin),i=1,ncorrao(iblock) )
          call flush(ifu_ac)
       end do
       close(ifu_hyb)
       close(ifu_ac)
    end do
    print *, "done." 

    deallocate( EMesh, delta )
    
  end subroutine CompHybFunc


  !*******************************************!
  ! Routine for orbital eigen-channel         !
  ! analysis with reduced transmission matrix !
  !*******************************************!
  subroutine EigenChannelAnalysis
    use parameters, only:  RedTransmB,RedTransmE, eta, EW1, EW2,EStep
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: ui, d_pi
    use numeric, only: CInv, CMatPow, CHDiag
    implicit none
    
    real :: rho, phi,DE,energy,delta,mindelta,tsave
    real, parameter :: dsmall = 1.0d-10
    integer :: NMaxEV = 10 ! At how many points to print out eign vectors: 2*NMaxEV+1

    integer :: NSD, NLD, NRD, NS1, NS2, N,NMax, info

    complex*16,dimension(NAOrbs,NAOrbs) :: SigmaL, SigmaR
    complex*16,dimension(:,:),allocatable :: GSD, SigmaLD, SigmaRD, GammaLD, GammaRD, GammaLDph, TS, TS0
    complex*16,dimension(:,:),allocatable :: gLD, gRD
    complex*16,dimension(:,:),allocatable :: VLS, VRS

    complex*16 :: cenergy, ci
    integer :: is, i, nchan, j, jmin, k

    real,dimension(:),allocatable :: tchan,tchan1,tchan2

    ! Dimension of scattering region inside device
    NSD = HiAOrbNo(RedTransmE)-LoAOrbNo(RedTransmB)+1
    ! Dimension of left region of device
    NLD = LoAOrbNo(RedTransmB)-1
    ! Dimension of right region of device
    NRD = NAOrbs-HiAOrbNo(RedTransmE)
    
    NS1 = LoAOrbNo(RedTransmB)
    NS2 = HiAOrbNo(RedTransmE)

    allocate( GSD(NSD,NSD), SigmaLD(NSD,NSD), SigmaRD(NSD,NSD), &
         GammaLD(NSD,NSD), GammaRD(NSD,NSD), GammaLDph(NSD,NSD), TS(NSD,NSD), TS0(NSD,NSD), &
         gLD(NLD,NLD), gRD(NRD,NRD), VLS(NLD,NSD), VRS(NRD,NSD), tchan(NSD), tchan1(NSD), tchan2(NSD) )

    !open(ifu_red,file='t.dat',status='unknown')

    print *
    print *, "-----------------------------------------------------------------------"
    print *, "--- Orbital Eigen-channel Analysis with reduced Transmission matrix ---"
    print *, "-----------------------------------------------------------------------"
    print *
    print *, "Begin of scattering region: AO ", NS1
    print *, "End of scattering region:   AO ", NS2
    print *

    do is=1,NSpin
       if(NSpin ==2 .and. is==1)print*,"Spin UP"
       if(NSpin ==2 .and. is==2)print*,"Spin DOWN"       

       NMax = int(EW2/EStep)

       !DE = d_zero
       !IF(NMaxE>0) DE = EWindow/DBLE(NMaxE)

       do N = -NMax,NMax!-NMaxE,NMaxE  
          energy = N*EStep
          cenergy = energy

          call CompSelfEnergies(is,cenergy,SigmaL,SigmaR)
       
          do i=1,NLD
             do j=1,NLD
                gLD(i,j) = (energy-shift+ui*eta)*SD(i,j) - HD(is,i,j) - SigmaL(i,j)
             end do
          end do

          do i=1,NRD
             do j=1,NRD
                gRD(i,j) = (energy-shift+ui*eta)*SD(NS2+i,NS2+j) - HD(is,NS2+i,NS2+j) - SigmaR(NS2+i,NS2+j)
             end do
          end do

          info = CInv( gLD )
          info = CInv( gRD )
          
          VLS = HD(is,1:NLD, NS1:NS2) - (energy-shift)*SD(1:NLD, NS1:NS2)
          VRS = HD(is,NS2+1:NAOrbs,NS1:NS2) - (energy-shift)*SD(NS2+1:NAOrbs,NS1:NS2)
          
          SigmaLD = matmul( conjg(transpose(VLS)), matmul( gLD, VLS ) )
          SigmaRD = matmul( conjg(transpose(VRS)), matmul( gRD, VRS ) )
          
          GSD = (energy-shift+ui*eta)*SD(NS1:NS2,NS1:NS2) & 
               - HD(is,NS1:NS2,NS1:NS2) - SigmaLD - SigmaRD

          info = CInv( GSD )

          GammaLD = ui*(SigmaLD-conjg(transpose(SigmaLD)))
          GammaRD = ui*(SigmaRD-conjg(transpose(SigmaRD)))

          call CMatPow( GammaLD, 0.5d0, GammaLDph )
          
          ![Gamma_L^1/2 G^a Gamma_R G^r Gamma_L^1/2] 
          TS = matmul( GammaLDph, conjg(transpose(GSD)) )
          TS = matmul( TS, GammaRD )
          TS = matmul( TS, GSD )
          TS = matmul( TS, GammaLDph )

          call CHDiag( TS, tchan, info )
          if( info /= 0 )then
             print*, "Error diagonalizing Reduced transmission matrix: info=", info
             return
          end if

          ! *** Ordering of eigenchannels ***
          if( N >= -NMax+2 ) call SeparateSpaghettis( tchan1, tchan2, tchan, TS, NSD )
          tchan1=tchan2
          tchan2=tchan

          write(ifu_red,'(F10.5,100E14.5)'), energy,(tchan(i),i=NSD,1,-1)

          if( N ==0 )then 
             print *, "Eigenchannel composition at Fermi level:"
             print *, "----------------------------------------"
             do nchan=NSD,1,-1
                print '(A,I2,A,F9.3)',"Channel ", NSD-nchan+1, ": Transmission = ", tchan(nchan)
                do i=1,NSD
                   ci = TS(i,nchan)
                   rho = abs(ci)
                   if( abs(real(ci)) < dsmall * abs(aimag(ci)) )then
                      phi = sign(0.5*d_pi,aimag(ci))
                   else
                      phi = atan( aimag(ci)/real(ci) )
                      if( real(ci) < 0 .and. aimag(ci) > 0 ) phi = phi + 0.5*d_pi
                      if( real(ci) < 0 .and. aimag(ci) < 0 ) phi = phi - 0.5*d_pi
                   end if
                   print '(I4,A,F8.4,F8.4)', i, " : ", rho, phi
                end do
             end do
          end if

       end do
       write(ifu_red,*)
       print *, " "
    end do

    deallocate( GSD, SigmaLD, SigmaRD, GammaLD, GammaRD, GammaLDph, TS, TS0, gLD, gRD, VLS, VRS , tchan, tchan1, tchan2 )

    !close(ifu_red)

  end subroutine EIGENCHANNELANALYSIS


  !
  ! *** Subroutine to separate individual eigen channel     ***
  ! *** transmissions (= spaghettis) using first derivative ***
  !
  subroutine SeparateSpaghettis( evals1, evals2, evals, evecs, N )
    implicit none
    
    ! eigen values at before last, last and actual energy
    real, dimension(N) :: evals1, evals2, evals
    ! eigenvectors at actual energy
    complex*16, dimension(N,N) :: evecs
    ! number of eigenvectors
    integer :: N
    
    integer :: i, j, jmin
    real :: yex, delta, mindelta, valsave

    complex*16, dimension(N) :: vecsave

    do i=1,N

       ! extrapolate actual eigenvalue from last 2 eigenvalues
       yex = 2.0d0*evals2(i)-evals1(i)

       ! Find actual eigenvalue which deviates minimally 
       ! from extrapolated value 
       mindelta = abs(yex-evals(i))
       jmin = i
       do j=i+1,N
          delta =  abs(yex-evals(j))
          if( delta < mindelta )then 
             jmin = j 
             mindelta = delta
          end if
       end do

       ! change eigenvector i with eigenvector jmin
       if( jmin /= i )then 
          vecsave = evecs(:,jmin)
          evecs(:,jmin)=evecs(:,i)
          evecs(:,i)=vecsave
          valsave = evals(jmin)
          evals(jmin) = evals(i)
          evals(i) = valsave
       end if
    end do
    
  end subroutine SeparateSpaghettis
  

  ! ***********************************
  ! Compute Self energies of electrodes
  ! ***********************************
  subroutine CompSelfEnergies( spin, cenergy, Sigma1, Sigma2 )
    use parameters, only: ElType, DD, UD, DU, Overlap
    use BetheLattice, only: CompSelfEnergyBL, LeadBL 
    use OneDLead, only: CompSelfEnergy1D, Lead1D
    use constants
    use util, only: PrintCMatrix, PrintRMatrix
   !use lapack_blas, only: zgemm
   !use blas95
   !use lapack95
    implicit none

    external               zgemm
    integer, intent(in) :: spin
    complex*16, intent(in) :: cenergy
    complex*16, dimension(:,:),intent(inout) :: Sigma1
    complex*16, dimension(:,:),intent(inout) :: Sigma2
    complex*16, dimension(NAOrbs,NAOrbs) :: temp
    integer :: is,omp_get_thread_num
    
!    if(DebugDev)Write(*,'(A)')"ENTERED Device/CompSelfEnergies"

    !Sigma1=(0.0,0.0)
    !Sigma2=(0.0,0.0)

    is = spin
    if( DD .and. spin == 1 ) is=2
    if( DD .and. spin == 2 ) is=1
    if( UD .and. spin == 1 ) is=1
    if( UD .and. spin == 2 ) is=2
    if( DU .and. spin == 1 ) is=2
    if( DU .and. spin == 2 ) is=1
    !write(ifu_log,*)omp_get_thread_num(),'in CompSelfEnergies',shift,cenergy,'lead 1'
    select case( ElType(1) )
    case( "BETHE" ) 
       call CompSelfEnergyBL( LeadBL(1), is, cenergy, Sigma1 )
       ! transform to non-orthogonal basis
       if( Overlap < -0.01 .and. .not. HDOrtho )then
          ! Sigma1 -> S^1/2 * Sigma1 * S^1/2
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, Sigma1,NAorbs, SPH,  NAOrbs, c_zero, temp,   NAOrbs)
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, SPH,   NAorbs, temp, NAOrbs, c_zero, Sigma1, NAOrbs)
       endif
    case( "1DLEAD" )
       call CompSelfEnergy1D( Lead1D(1), is, cenergy, Sigma1 )
    case( "GHOST" )
        continue
    end select

    if( DD .and. spin == 1 ) is=2
    if( DD .and. spin == 2 ) is=1
    if( UD .and. spin == 1 ) is=2
    if( UD .and. spin == 2 ) is=1
    if( DU .and. spin == 1 ) is=1
    if( DU .and. spin == 2 ) is=2
    
    !write(ifu_log,*)omp_get_thread_num(),'in CompSelfEnergies',shift,cenergy,'lead 2'
    select case( ElType(2) )
    case( "BETHE" )
       call CompSelfEnergyBL( LeadBL(2), is, cenergy, Sigma2 )
       ! transform to non-orthogonal basis 
       if( Overlap < -0.01 .and. .not. HDOrtho )then
          ! Sigma2 -> S^1/2 * Sigma2 * S^1/2
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, Sigma2,NAorbs, SPH,  NAOrbs, c_zero, temp,   NAOrbs)
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, SPH,   NAorbs, temp, NAOrbs, c_zero, Sigma2, NAOrbs)
       endif
    case( "1DLEAD" )
       call CompSelfEnergy1D( Lead1D(2), is, cenergy, Sigma2 )
    case( "GHOST" )
        continue
    end select

    if(DebugDev)then
!      Write(*,'(A)')"Device/CompSelfEnergies Sigma1"
!      call PrintCMatrix(Sigma1)
!      Write(*,'(A)')"Device/CompSelfEnergies Sigma2"
!      call PrintCMatrix(Sigma2)
    end if

  end subroutine CompSelfEnergies


  !
  ! *** Integrand for charge integration on complex contour ***
  !
  real function DDOS( E0, R, phi )
    use constants, only: c_zero, ui, d_pi

    implicit none
    real, intent(in) :: phi, E0, R
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal
    complex*16 :: TrGS, z

    z = E0 - R*(cos(phi) - ui*sin(phi)) 

    TrGS=c_zero
    do ispin=1,NSpin
       call gplus(z,green,gammar,gammal)
       !do i=1,NAOrbs
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             TrGS = TrGS + green(i,j)*SD(j,i)
          end do
       end do
    end do
    ! Account for spin degeneracy
    if(NSpin==1) TrGS = TrGS * 2.0d0
    DDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
  end function DDOS

  !
  ! *** Integrand for charge integration on complex contour ***
  !
  real function CDOS( E0, R, phi )
    use constants, only: c_zero, ui, d_pi

    real, intent(in) :: phi, E0, R
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal
    complex*16 :: TrGS, z

    z = E0 - R*(cos(phi) - ui*sin(phi)) 

    TrGS=c_zero
    !do ispin=1,NSpin
       call gplus(z,green,gammar,gammal) ! Depends on ispin.
       !do i=1,NAOrbs
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             TrGS = TrGS + green(i,j)*SD(j,i)
          end do
       end do
    !end do ! The ispin loop is in WorkEnergy.
    ! Account for spin degeneracy
    !if(NSpin==1) TrGS = TrGS * 2.0d0 ! I will call this for each spin component.
    CDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
    !print *, "E0=", E0, " CDOS", CDOS
  end function CDOS

  !
  ! *** Integrand for charge integration on real contour ***
  !
  real function RDOS( E0 )
    use Cluster, only : hiaorbno, loaorbno
    use g09Common, only: GetNAtoms
    !complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG
    use constants, only: c_zero, ui, d_pi, d_zero

    real :: DOS
    real, intent(in) :: E0
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal,SG
    complex*16 :: TrGS, z

    z=dcmplx(E0)
    DOS=d_zero
    !AtomDOS=d_zero
    do ispin=1,NSpin
          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          !call gplus(z,green,GammaR,GammaL)
          call glesser(z,green)
          ! Mulliken DOS 
          SG = matmul( SD, green )
          ! computing total DOS
          DOS=d_zero
          !AtomDOS=d_zero
          do j=1,GetNAtoms()
            do i=LoAOrbNo(j),HiAOrbNo(j)
              !AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
              !DOS=DOS-dimag(SG(i,i))/d_pi
              DOS=DOS+dimag(SG(i,i))/d_pi
            end do
          end do
    end do
    ! Account for spin degeneracy
    !if(NSpin==1) TrGS = TrGS * 2.0d0
    !if(NSpin==1) DOS = DOS * 2.0d0 ! I do this later.
    RDOS=DOS
    !print *, "E0=", E0, " RDOS", RDOS
  end function RDOS

  !
  ! *** Integrand for charge integration on complex contour ***
  !
  real function CDOSE( E0, R, phi )
    use constants, only: c_zero, ui, d_pi

    real, intent(in) :: phi, E0, R
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal
    complex*16 :: TrGS, z

    real :: CDOS

    z = E0 - R*(cos(phi) - ui*sin(phi)) 

    TrGS=c_zero
    !do ispin=1,NSpin
       call gplus(z,green,gammar,gammal) ! Depends on ispin.
       !call glesser(z,green) ! As seen in Stefanucci.
       !do i=1,NAOrbs
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             TrGS = TrGS + green(i,j)*SD(j,i)
          end do
       end do
    !end do ! The ispin loop is in WorkEnergy.
    ! Account for spin degeneracy
    !if(NSpin==1) TrGS = TrGS * 2.0d0 ! I will call this for each spin component.
    CDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
!    CDOSE = CDOS*E0
    !print *, "E0=", E0, " CDOS", CDOS
    if( NSpin == 2 .and. SPINLOCK )then
      if( ispin == 1 )CDOSE = CDOS*(E0-shiftup)
      if( ispin == 2 )CDOSE = CDOS*(E0-shiftdown)
    else if( NSpin == 2 .and. .not. SPINLOCK )then
      CDOSE = CDOS*(E0-shift)
    else if ( NSpin == 1 )then
      CDOSE = CDOS*(E0-shift)
    end if
  end function CDOSE

  !
  ! *** Integrand for charge integration on real contour ***
  !
  real function RDOSE( E0 )
    use Cluster, only : hiaorbno, loaorbno
    use g09Common, only: GetNAtoms
    !complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG
    use constants, only: c_zero, ui, d_pi, d_zero

    real :: DOS
    real, intent(in) :: E0
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal,SG
    complex*16 :: TrGS, z

    z=dcmplx(E0)
    DOS=d_zero
    !AtomDOS=d_zero
    do ispin=1,NSpin
          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(z,green,GammaR,GammaL)
          !call glesser(z,green) ! As seen in Stefanucci.
          ! Mulliken DOS 
          SG = matmul( SD, green )
          ! computing total DOS
          DOS=d_zero
          !AtomDOS=d_zero
          do j=1,GetNAtoms()
            do i=LoAOrbNo(j),HiAOrbNo(j)
              !AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
              !DOS=DOS-dimag(SG(i,i))/d_pi
              DOS=DOS+dimag(SG(i,i))/d_pi
            end do
          end do
    end do
    ! Account for spin degeneracy
    !if(NSpin==1) TrGS = TrGS * 2.0d0
    !if(NSpin==1) DOS = DOS * 2.0d0 ! I do this later.
!    RDOSE=DOS*E0
    !print *, "E0=", E0, " RDOS", RDOS
    if( NSpin == 2 .and. SPINLOCK )then
      if( ispin == 1 )RDOSE = DOS*(E0-shiftup)
      if( ispin == 2 )RDOSE = DOS*(E0-shiftdown)
    else if( NSpin == 2 .and. .not. SPINLOCK )then
      RDOSE = DOS*(E0-shift)
    else if ( NSpin == 1 )then
      RDOSE = DOS*(E0-shift)
    end if
  end function RDOSE

  !
  ! *** Integrand for charge integration on real contour ***
  !
  real function RDOSOld( E0 )
    use Cluster, only : hiaorbno, loaorbno
    use g09Common, only: GetNAtoms
    !complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG
    use constants, only: c_zero, ui, d_pi, d_zero

    real :: DOS
    real, intent(in) :: E0
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal,SG
    complex*16 :: TrGS, z

    !z = E0 - R*(cos(phi) - ui*sin(phi)) 
    !z = E0 +ui*d_zero
    !z = E0
    !z = E0 - Ex
    !cenergy=dcmplx(energy)
    !z=dcmplx(E0 - Ex)
    z=dcmplx(E0)
    DOS=d_zero
    do ispin=1,NSpin
          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          !call gplus(z,green,gammar,gammal)
          call glesser(z,green)
          ! Mulliken DOS 
          SG = matmul( SD, green )
          ! computing total DOS

          !AtomDOS=d_zero
          !do j=1,GetNAtoms()
          !  do i=LoAOrbNo(j),HiAOrbNo(j)
          !   !AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
          !   DOS=DOS-dimag(SG(i,i))/d_pi
          !  end do
          !end do
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             DOS = DOS-dimag(green(i,j)*SD(j,i))/d_pi
          end do
       end do
    end do
    !DDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
    !print *, "E0 = ", E0,"RDOS =", -DIMAG(E0*TrGS)/d_pi
    !RDOS = -DIMAG(E0*TrGS)/d_pi
    ! Account for spin degeneracy
    !if(NSpin==1) TrGS = TrGS * 2.0d0
    if(NSpin==1) DOS = DOS * 2.0d0 ! I do this later.
    RDOSOld=DOS
    !print *, "E0=", E0, " RDOS", RDOS
  end function RDOSOld


!
  ! *** Integrand for charge integration on real contour ***
  !
  real function RDOStimesE( E0 )
    use Cluster, only : hiaorbno, loaorbno
    use g09Common, only: GetNAtoms
    !complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG
    use constants, only: c_zero, ui, d_pi, d_zero

    real :: DOS, RDOS
    real, intent(in) :: E0
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal,SG
    complex*16 :: TrGS, z

    !z = E0 - R*(cos(phi) - ui*sin(phi)) 
    !z = E0 +ui*d_zero
    !z = E0
    
    !cenergy=dcmplx(energy)
    z=dcmplx(E0)
          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(z,green,gammar,gammal)
          !call glesser(z,green)

          ! Mulliken DOS 
          SG = matmul( SD, green )
          ! computing total DOS
          DOS=d_zero
          !AtomDOS=d_zero
          do j=1,GetNAtoms()
            do i=LoAOrbNo(j),HiAOrbNo(j)
             !AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
             DOS=DOS-dimag(SG(i,i))/d_pi
            end do
          end do

    !DDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
    !print *, "E0 = ", E0,"RDOS =", -DIMAG(E0*TrGS)/d_pi
    !RDOS = -DIMAG(E0*TrGS)/d_pi
    ! Account for spin degeneracy
    !if(NSpin==1) TrGS = TrGS * 2.0d0
    if(NSpin==1) DOS = DOS * 2.0d0 ! DO IT HERE
    !print *, "E0=", E0, " DOS", DOS
    RDOS=DOS
    ! I have to substract the Fermi level.
    RDOStimesE = RDOS*(E0 - shift)
    !RDOStimesE = RDOS*(E0)
  end function RDOStimesE
  
  ! 
  ! *** Total charge up to energy E ***
  ! 
  ! Integrates Green's function along imaginary
  ! axis so that no lower energy bound is needed.
  ! Uses routine qromo of Num. Rec. with midpnt
  ! rule up to some point on the path and midinf 
  ! rule for the tail.
  !
  real function TotEnergyOld( E )
    use constants, only: d_pi
   !use parameters, only: ChargeAcc
    use numeric, only: midpnt, midinf, qromoc
    implicit none
    
    real, intent(in) :: E
    integer omp_get_thread_num
    real :: s, q, y0, y1
    real :: E0, intDOSE, R
    integer :: i, iwindow, nwindow
    real :: windowwidth

    print *, "I am in TotEnergy"
    !    R  = 0.5*(E - EMin)
    !R  = (E - EMin)
    !y0 = EMin
    E0 = E
    !write(ifu_log,*)omp_get_thread_num(),'in TotCharge',E0

    IntDOSE = 0.0d0
    ! integration from [0:y0] with midpnt rule
    ! and from [y0:inf] with midinf rule 
    ! (assumes integrand decaying approx. ~1/x)
    !y0=-20.0d0
    !y0=-25.0d0
    !y0 = (EMin-E0)/2.0d0

    !y1=E0


    print *, "[EMin,EMax] [", EMin,",", EMax,"]"

    windowwidth = 5.0
    nwindow = floor(abs(EMin/windowwidth))

    do iwindow=1,nwindow
      !y0=-10.0d0
      !y1=0.0d0
      y0 = -windowwidth*(iwindow)
      y1 = -windowwidth*(iwindow-1)
      !y1 = E0
      print *, "Integrating interval nr.", iwindow
      print *, "Integrating interval [y0,y1] [", y0,",", y1,"]"
      s = 0.0d0
      !call qromoc( RDOStimesE,y0,1.d30,s,midinf)
      !call qromoc( RDOStimesE,y0,0.0d0,s,midpnt)
      call qromoc( RDOStimesE,y0,y1,s,midpnt)
      !do i=1,SIZE(s)
      !  Write(ifu_diagham,*)i,s(i),
      print *, "Partial Integral =",s
      IntDOSE = IntDOSE + s
    print *, "Accumulated Integral =",IntDOSE
    end do

    y0=EMin
    y1 = -windowwidth*(nwindow)
    !y1=-60.0d0
      print *, "Integrating interval nr.", iwindow
    print *, "Integrating interval [y0,y1] [", y0,",", y1,"]"
    s = 0.0d0
    !call qromoc( RDOStimesE,y0,1.d30,s,midinf)
    !call qromoc( RDOStimesE,y0,0.0d0,s,midpnt)
    call qromoc( RDOStimesE,y0,y1,s,midpnt)
    print *, "Partial Integral =",s
    IntDOSE = IntDOSE + s
    print *, "Accumulated Integral =",IntDOSE

    !q = q + LeadBL(WhichLead)%NAOrbs

    !!print *, "Whichlead = ", WhichLead
    !!PRINT *, " E=", E, "    TotCharge=", q

    !TotCharge = q-ChargeOffSet
    TotEnergyOld = IntDOSE
  end function TotEnergyOld

  ! 
  ! *** Total charge up to energy E ***
  ! 
  ! Integrates Green's function along imaginary
  ! axis so that no lower energy bound is needed.
  ! Uses routine qromo of Num. Rec. with midpnt
  ! rule up to some point on the path and midinf 
  ! rule for the tail.
  !
real function FullEnergy( E )
    use constants, only: d_zero, d_pi
    use parameters, only: ChargeAcc, biasvoltage
    use numeric, only: midpnt, midinf, qromoc, gauleg
    implicit none
    
    real, intent(in) :: E
    integer omp_get_thread_num
    real :: s, q, y0, y1
    real :: E0, intDOSE, R
    integer :: i
    real :: Chrg

!-------------------------------------------------------------------------------
    integer, parameter :: nmax = 1023
    real :: qq, w_j, phi_j
    integer :: n, n1, n2, j
    
    real, dimension(nmax) :: x, w
!-------------------------------------------------------------------------------

    print *, "I am in FullEnergy"

    Chrg = 0.0d0
    IntDOSE = 0.0d0
    q = 0.0d0


  !do ispin=1,NSpin
    !    R  = 0.5*(E - EMin)
    !R  = (E - EMin)
    !y0 = EMin
    !E0 = E
    !write(ifu_log,*)omp_get_thread_num(),'in TotCharge',E0
    IntDOSE = 0.0d0
    ! integration from [0:y0] with midpnt rule
    ! and from [y0:inf] with midinf rule 
    ! (assumes integrand decaying approx. ~1/x)
    !y0=-20.0d0
    !y0=-25.0d0
    !y0 = (EMin-E0)/2.0d0

    !y1=E0
    print *, "ispin", ispin
    print *, "[EMin,EMax] [", EMin,",", EMax,"]"

    !---------------------------------------------------------------------
    !----------- BELOW THE BIAS WINDOW -----------------------------------
    !---------------------------------------------------------------------
    ! Integration contour parameters:
    !E0 = 0.5*(E + EMin) ! Non-modified code.
    !R  = 0.5*(E - EMin) ! Non-modified code.
    E0 = 0.5*(E-dabs(biasvoltage/2.0) + EMin)
    R  = 0.5*(E-dabs(biasvoltage/2.0) - EMin)
    !E0 = 0.5*(d_zero-dabs(biasvoltage/2.0) + EMin)
    !R  = 0.5*(d_zero-dabs(biasvoltage/2.0) - EMin)
    ! Computing integral of DOS over 
    ! complex contour using Gauss-Legendre 
    ! quadrature
    n=1
    do  
       q = d_zero
       do j=1,n
          call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
          !q = q + w(j)*CDOS( E0, R, x(j) )
          q = q + w(j)*CDOSE( E0, R, x(j) )
       end do
       !print *, "j:",j,"TotCharge", q
       if( n > 1 .and. (q == d_zero .or. abs(q-qq) < ChargeAcc*NCDEl ) ) exit  
       n=2*n+1
       if( n > nmax )then
          print *, "TotCharge/gaussian quadrature has not converged after", nmax, " steps."
          Chrg = 2.0d0*(NCDAO2-NCDAO1+1) - 10.0d0*ChargeAcc*NCDEl
          return
       end if
       qq = q
    end do
    !print *, "gaussian quadrature converged after", n, " steps. Error:", abs(q-qq)

    Chrg = q
    IntDOSE = q
    print *, "Partial Integral =",IntDOSE
    !---------------------------------------------------------------------
    !----------- WITHIN THE BIAS WINDOW ----------------------------------
    !---------------------------------------------------------------------
    if (biasvoltage /= 0.d0) then
      q = d_zero
      y0 = E-dabs(biasvoltage/2.0)
      y1 = E+dabs(biasvoltage/2.0)
      !y0 = d_zero-dabs(biasvoltage/2.0)
      !y1 = d_zero+dabs(biasvoltage/2.0)
      !y1 = E0
      print *, "Integrating interval [y0,y1] [", y0,",", y1,"]"
      !call qromoc( RDOS,y0,y1,s,midpnt)
      call qromoc( RDOSE,y0,y1,s,midpnt)
      print *, "Partial Integral =",s
      q = q + s
      print *, "Accumulated Integral =",q
      
      Chrg = Chrg + q
      IntDOSE = IntDOSE + q

    end if
  !end do ! ispin loop.
    !q = q + LeadBL(WhichLead)%NAOrbs

    !!print *, "Whichlead = ", WhichLead
    !!PRINT *, " E=", E, "    TotCharge=", q

    PRINT *, " ispin=", ispin, "    TotCharge=", Chrg
    PRINT *, " ispin=", ispin, "    TotEnergy=", IntDOSE
    !PRINT *, " E=", d_zero, "    TotCharge=", q
    !PRINT *, " E=", d_zero, "    TotEnergy=", IntDOSE

    !TotCharge = q-Charg
    FullEnergy = IntDOSE
end function FullEnergy

  ! 
  ! *** Total charge up to energy E, lower bound is EMin ***
  ! 
  real function TotCharge( E )
    use constants, only: d_zero, d_pi
    use parameters, only: ChargeAcc
    use numeric, only: gauleg 
    implicit none
    
    real, intent(in) :: E

    integer, parameter :: nmax = 2047
    
    real :: q,qq, E0, R, w_j, phi_j
    integer :: n, n1, n2, j, i
    
    real, dimension(nmax) :: x, w

!    integer :: c1
!    real :: t1

    !c1 = 0.0d0
    
    ! Integration contour parameters:
    E0 = 0.5*(E + EMin)
    R  = 0.5*(E - EMin)
    ! Computing integral of DOS over 
    ! complex contour using Gauss-Legendre 
    ! quadrature
    n=1
    do  
       q = d_zero
       do j=1,n
          call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
          q = q + w(j)*DDOS( E0, R, x(j) )
       end do
       if( n > 1 .and. (q == d_zero .or. abs(q-qq) < ChargeAcc*NCDEl ) ) exit  
       n=2*n+1
     
!       if( mod(n,3) == 0)then
!          WRITE(*,*) "STEPS: ",n
!          CALL CPU_TIME(t1)
!          WRITE(*,*) "CPU_TIME: ",t1
!          WRITE(*,*) "TotCharge = ",q
!          !CALL SYSTEM_CLOCK(c1)
!          !WRITE(*,*) "SYSTEM_CLOCK: ",c1
!       end if

       if( n > nmax )then
          !WRITE(*,*) "STEPS: ",n
          !CALL CPU_TIME(t1)
          !WRITE(*,*) "CPU_TIME: ",t1
          !WRITE(*,*) "TotCharge = ",q
          !CALL SYSTEM_CLOCK(c1)
          !WRITE(*,*) "SYSTEM_CLOCK: ",c1
          print *, "TotCharge/gaussian quadrature has not converged after", nmax, " steps."
          TotCharge = 2.0d0*(NCDAO2-NCDAO1+1) - 10.0d0*ChargeAcc*NCDEl
          return
       end if
       qq = q
    end do
    !print *, "gaussian quadrature converged after", n, " steps. Error:", abs(q-qq)

    TotCharge = q
  end function TotCharge

!-------------------------------------------------------------------------------
  ! This function based on DDOS is the integrand for E*DOS. Written by C. Salgado.
  real function DDOStimesE0( E0, R, phi )
!    use Cluster, only : hiaorbno, loaorbno
!#ifdef G03ROOT
!    use g03Common, only: GetNAtoms
!#endif
!#ifdef 1
!    use g09Common, only: GetNAtoms
!#endif
    use constants, only: c_zero, ui, d_pi

    real, intent(in) :: phi, E0, R
    real :: DDOS
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal
    complex*16 :: TrGS, z
    complex*16,dimension(NAOrbs,NAOrbs) :: SG
    real :: DOS

    !z = E0 - R*(cos(phi) - ui*sin(phi)) 
    z = dcmplx(E0)
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!------ ALTERNATIVE 1 -----------------------------------------------
    TrGS=c_zero
    do ispin=1,NSpin
       call gplus(z,green,gammar,gammal)
       !do i=1,NAOrbs
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             TrGS = TrGS + green(i,j)*SD(j,i)
          end do
       end do
    end do
    ! Account for spin degeneracy
    if(NSpin==1) TrGS = TrGS * 2.0d0
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!------ ALTERNATIVE 2 -----------------------------------------------
 !   SG = matmul( SD, green )
 !   ! computing total DOS
 !   DOS=d_zero
 !   !AtomDOS=d_zero
 !   do j=1,GetNAtoms()
 !     do i=LoAOrbNo(j),HiAOrbNo(j)
 !       !AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
 !       DOS=DOS-dimag(SG(i,i))/d_pi
 !     end do
 !   end do
 !   ! Account for spin degeneracy
 !   !if(NSpin==1) DOS = DOS * 2.0d0
!--------------------------------------------------------------------
!--------------------------------------------------------------------
    DDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
    !DDOStimesE0 = DDOS*E0
    !DDOStimesE0 = DDOS*(E0 - R*(cos(phi)))
    DDOStimesE0 = DDOS*(E0 - shift)
  end function DDOStimesE0

  ! 
  ! *** Total charge up to energy E, lower bound is EMin ***
  ! 
  real function IntDDOStimesE( E )
    use constants, only: d_zero, d_pi
    use parameters, only: ChargeAcc
    use numeric, only: gauleg 
    implicit none
    
    real, intent(in) :: E

    integer, parameter :: nmax = 2047
    
    real :: q,qq, E0, R, w_j, phi_j
    real :: intDOSE, intDOSEE
    integer :: n, n1, n2, j, i
    
    real, dimension(nmax) :: x, w
    
    print *, "I am in IntDDOStimesE" 

    ! Integration contour parameters:
    E0 = 0.5*(E + EMin)
    R  = 0.5*(E - EMin)
    ! Computing integral of DOS*E over 
    ! complex contour using Gauss-Legendre 
    ! quadrature
    n=1
    do  
       !q = d_zero
       intDOSE = d_zero
       do j=1,n
          !call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
   call gauleg(EMin,E0,x(1:n),w(1:n),n)
          !q = q + w(j)*DDOS( E0, R, x(j) )
          intDOSE = intDOSE + w(j)*DDOStimesE0( E0, R, x(j) )
   print *, "j",j,"intDOSE", intDOSE
          
       end do
       !if( n > 1 .and. (q == d_zero .or. abs(q-qq) < ChargeAcc*NCDEl ) ) exit  
       if( n > 1 .and. (intDOSE == d_zero .or. abs(intDOSE-intDOSEE) < 1.0d-6 ) ) exit 
       n=2*n+1
       if( n > nmax )then
          print *, "IntDDOStimesE/gaussian quadrature has not converged after", nmax, " steps."
          !TotCharge = 2.0d0*(NCDAO2-NCDAO1+1) - 10.0d0*ChargeAcc*NCDEl
          IntDDOStimesE = intDOSE
          return
       end if
       !qq = q
       intDOSEE = intDOSE
    end do
    !print *, "gaussian quadrature converged after", n, " steps. Error:", abs(q-qq)

    !TotCharge = q
    IntDDOStimesE = intDOSE
  end function IntDDOStimesE
  !-----------------------------------------------------------

  ! *************************************
  !  Estimates upper/lower energy 
  !  boundary EMin/EMax,
  !  above/below which DOS is gauranteed 
  !  to be zero
  ! *************************************
  subroutine FindEnergyBounds
    use parameters, only: ChargeAcc,eta

    implicit none
    integer :: i, cond,k 
    real :: EStep, Q
    
    print *, "Searching energy boundaries [EMin, EMax] such that"
    print '(A,I4)', " Total Integrated Spectral Weight (TISW)=", 2*(NCDAO2-NCDAO1+1)

    EStep = 10.0d0 +10000.0 *eta

    do
       Q = TotCharge( EMax )
       print '(A,F12.5,A,F12.5,A,F10.5)', " EMin=", EMin, "  EMax=", EMax , "  TISW=", Q
       if( abs(Q - 2.0d0*(NCDAO2-NCDAO1+1)) < ChargeAcc*NCDEl*10.0 ) then
          exit
       end if
       EMin = EMin - EStep
       EMax = EMax + EStep
    end do
    print *, "--------------------------------------------------"

  end subroutine FindEnergyBounds

!----------------------------------------------------------------------------------------------------------------
!-------------- ADAPTED ROUTINES TO INTEGRATE ENERGY ------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------

!ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!        IntDOSE:  The value of the energy integral. Interval [-1,1]           c
!ccccccccccccccccccccccccccccccc
  !subroutine intchenergy(Er,El,M)
  real function intchenergy(Er,El,M)
    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc 
   
    real,intent(in) :: Er,El
    integer,intent(inout) :: M
    real, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm 
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

    complex*16, dimension(NAOrbs,NAOrbs) :: gammar, gammal
    real :: DOS

    if(DebugDev)Write(*,'(A)')"ENTERED Device/intchenergy"
      intchenergy = d_zero

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      !call glesser(E0,green)
      call gplus(E0,green,gammal,gammar)
      do i=1,NAOrbs
       do j=1,NAOrbs
        LocPDOUT(ispin,i,j) = -ui*green(i,j)/(2*pi)
        p(i,j) = LocPDOUT(ispin,i,j)
        DOS = DOS-ui*green(i,j)/(2*pi)
       enddo
      enddo
      intchenergy = intchenergy + DOS*E0
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,NAOrbs
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*LocPDOUT(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       !call glesser(Em,greenm)
       call gplus(Em,greenm,gammal,gammar)
!!$OMP  SECTION
       !call glesser(Ep,greenp)
       call gplus(Ep,greenp,gammal,gammar)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,NAOrbs
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
            DOS=DOS-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
          intchenergy=intchenergy+DOS*(0.5d0*Em+0.5d0*Ep)
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             LocPDOUT(ispin,i,j)=LocPDOUT(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,NAOrbs
       do j=1,NAOrbs
        rCHp=dble(LocPDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(LocPDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(LocPDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(LocPDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,NAOrbs
      do j=1,NAOrbs
        rCHp=dble(LocPDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(LocPDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(LocPDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(LocPDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        LocPDOUT(ispin,i,j) = 16*LocPDOUT(ispin,i,j)/(3*(n+1))
        LocPDOUT(ispin,i,j) = LocPDOUT(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      DOS = 16*DOS/(3*(n+1))
      DOS = DOS*(El-Er)/2
      intchenergy = 16*intchenergy/(3*(n+1))
      intchenergy = intchenergy*(El-Er)/2
      write(ifu_log,'(A51,i4)')' Integration of POUT has needed a max no. of loops=',(((n-1)/2)+1)/2

      return
    !end subroutine intchenergy
    end function intchenergy

  !ccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        IntDOSE:  On output, the value of the integral (energy density).      c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !ccccccccccccccccccccccccccccccc
  !subroutine intpjenergy(rrr,bi,Emi,M,Eq,intpjIntDOSE)
  real function intpjenergy(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui

    implicit none

    real,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real, dimension(NAOrbs,NAOrbs) :: PDP

    real :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk,omp_get_thread_num,omp_get_num_threads
    real, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn 

    logical :: omp_get_nested

    complex*16, dimension(NAOrbs,NAOrbs) :: SG,gammar, gammal
    real :: DOS

    !real,intent(out) :: intpjIntDOSE

    if(DebugDev)Write(*,'(A)')"ENTERED Device/intpjenergy"
 
    a = 1.d0/d_pi
    b = bi   ! This is input 0.0
    Em = Emi   ! This is input pi
    LocPD(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
    !call gplus0(E0,green)
    call gplus(E0,green,gammar,gammal)
    CH = 0.d0
    intpjenergy = 0.d0
    !---------------------------------------------------------------------------
    !-------------- OLD CODE ---------------------------------------------------
    !---------------------------------------------------------------------------
    !do i = 1,NAOrbs
    !   do j =1,NAOrbs
    !      LocPD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
    !      CH = CH + LocPD(ispin,i,j)*SD(j,i)
    !      !intpjIntDOSE = intpjIntDOSE + (real(E0)-shift)*PD(ispin,i,j)*SD(j,i)
    !      !intpjIntDOSE = intpjIntDOSE + (real(E0))*PD(ispin,i,j)*SD(j,i)
    !   enddo
    !enddo
    !---------------------------------------------------------------------------
    !-------------- NEW CODE ---------------------------------------------------
    !---------------------------------------------------------------------------
    !SG = matmul( SD, green ) ! Not necessary.
    ! computing total DOS
    DOS=d_zero
    !AtomDOS=d_zero
    do i = 1,NAOrbs
      do j =1,NAOrbs
        LocPD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
        CH = CH + LocPD(ispin,i,j)*SD(j,i)
        !AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
        !DOS=DOS-dimag(SD(i,j)*green(j,i))/d_pi
        DOS=DOS-a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
        
      end do
    end do
    intpjenergy=DOS*E0

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_nested(.true.)
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------' 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
    !chunk=max(((n+1)/2)/omp_get_num_threads(),1)
     chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!!$OMP  DO SCHEDULE(DYNAMIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          !call gplus0(EE(k),greenn(k,:,:))
          call gplus(EE(k),greenn(k,:,:),gammar,gammal)
       end do
!!$OMP  END DO
!!$OMP  END PARALLEL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
             !intpjIntDOSE = intpjIntDOSE + (real((EE(1)+EE(2))/2.0)-shift)*PDP(i,j)*SD(j,i)
             DOS=DOS-a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
       intpjenergy = intpjenergy + DOS*(0.5d0*real(EE(1)) + 0.5d0*real(EE(2)))
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             LocPD(ispin,i,j)=LocPD(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    !intpjIntDOSE = 0.d0
    do k=1,NAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]   
       do k1=1,NAOrbs
          CH = CH + LocPD(ispin,k,k1)*SD(k1,k)
          !intpjIntDOSE = intpjIntDOSE + (real((EE(1)+EE(2))/2.0)-shift)*PD(ispin,k,k1)*SD(k1,k)
          !intpjIntDOSE = intpjIntDOSE + (real((EE(1)+EE(2))/2.0))*PD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    !if ((CH-xp)*(CH-xp)*16.gt.(3*(n+1)*abs(CH-q)*PAcc*NCDEl)/1000000.and.n.le.M*1000000) goto 1   ! TOLERANCE REDUCED BY A FACTOR OF 10, TO GET A BETTER ACCURACY.
    ! Test for successfullness and integral final value
    M = 0
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    !if ((CH-xp)*(CH-xp)*16.gt.(3*(n+1)*abs(CH-q)*PAcc*NCDEl)/1000000) M = 1
    CH = 16*CH/(3*(n+1))
    !intpjIntDOSE = 16*intpjIntDOSE/(3*(n+1))
    DOS = 16*DOS/(3*(n+1))
    intpjenergy = 16*intpjenergy/(3*(n+1))
    do k=1,NAOrbs
       do l=1,NAOrbs
          LocPD(ispin,k,l) = 16*LocPD(ispin,k,l)/(3*(n+1))
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P has needed a max no. of loops=',(((n-1)/2)+1)/2
    !print*,"intpjenergy IntDOSE = ", intpjIntDOSE
    return
  !end subroutine intpjenergy
  end function intpjenergy


!-----------------------------------------------------------------------------------------------
!----------- LOCALINTCH AND LOCALINTPJ ---------------------------------------------------------
!-----------------------------------------------------------------------------------------------

!ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccc
  subroutine localintch(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc 
   
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
    use g09Common, only: GetNAtoms, GetAtmChg

    real,intent(in) :: Er,El
    integer,intent(inout) :: M
    real, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

    complex*16, dimension(NAOrbs,NAOrbs) :: gammar, gammal

    !integer,intent(in) :: iatom

    !print*,"I am in intch" 
    !write(ifu_log,*)'-------------------------------------'
    !write(ifu_log,*)'------- I am in localintch ----------'
    !write(ifu_log,*)'-------------------------------------'

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      !call glesser(E0,green)
      call gplus(E0,green,gammar,gammal)
      !do i=1,NAOrbs
      do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
       do j=1,NAOrbs
        LocPDOUT(ispin,i,j) = -ui*green(i,j)/(2*pi)
        p(i,j) = LocPDOUT(ispin,i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      !do i=1,NAOrbs
      do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*LocPDOUT(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       !call glesser(Em,greenm)
       call gplus(Em,green,gammar,gammal)
!!$OMP  SECTION
       !call glesser(Ep,greenp)
       call gplus(Ep,green,gammar,gammal)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          !do i=1,NAOrbs
          do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       !do i = 1,NAOrbs
       do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
          do j = 1,NAOrbs
             LocPDOUT(ispin,i,j)=LocPDOUT(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      !do i=1,NAOrbs
      do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
       do j=1,NAOrbs
        rCHp=dble(LocPDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(LocPDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(LocPDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(LocPDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        !if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*(GetAtmChg(spinatom)).and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        !if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*(GetAtmChg(spinatom)).and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      !do i=1,NAOrbs
      do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
      do j=1,NAOrbs
        rCHp=dble(LocPDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(LocPDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(LocPDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(LocPDOUT(ispin,i,j)-q(i,j))
        !if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*(GetAtmChg(spinatom))) M = 1
        !if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*(GetAtmChg(spinatom))) M = 1
        LocPDOUT(ispin,i,j) = 16*LocPDOUT(ispin,i,j)/(3*(n+1))
        LocPDOUT(ispin,i,j) = LocPDOUT(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      write(ifu_log,'(A51,i4)')' Integration of POUT has needed a max no. of loops=',(((n-1)/2)+1)/2

      return
    end subroutine localintch

  !ccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !ccccccccccccccccccccccccccccccc
  subroutine localintpj(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui

    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
    use g09Common, only: GetNAtoms, GetAtmChg


    implicit none
    real,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real, dimension(NAOrbs,NAOrbs) :: PDP

    real :: a,b,Em,S0,c0,x0,er0,der0,CH,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk,omp_get_thread_num,omp_get_num_threads
    real, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn 

    complex*16, dimension(NAOrbs,NAOrbs) :: gammar, gammal

    logical :: omp_get_nested

    !integer,intent(in) :: iatom

    !print*,"I am in intpj" 
    !write(ifu_log,*)'-------------------------------------'
    !write(ifu_log,*)'------- I am in localintpj ----------'
    !write(ifu_log,*)'-------------------------------------'
 
    a = 1.d0/d_pi
    b = bi
    Em = Emi
    LocPD(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
    

    !--------------------------------------------------------------------------------
    !------------- ORIGINAL CODE ----------------------------------------------------
    !--------------------------------------------------------------------------------
    !call gplus0(E0,green)  
    !CH = 0.d0
    !do i = 1,NAOrbs
    !   do j =1,NAOrbs
    !      PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
    !      CH = CH + PD(ispin,i,j)*SD(j,i)
    !   enddo
    !enddo
    !do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
    !   do j =1,NAOrbs
    !      PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
    !      CH = CH + PD(ispin,i,j)*SD(j,i)
    !   enddo
    !enddo
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !------------ MODIFIED CODE ----------------------------------------------------
    !-------------------------------------------------------------------------------
    !cenergy=dcmplx(energy)

    !*********************************************************************
    !* Evaluation of the retarded "Green" function and coupling matrices *
    !*********************************************************************
    !call gplus(cenergy,Green,GammaR,GammaL)
    call gplus(E0,green,gammar,gammal)
    !call gplus0(E0,green)
    ! Mulliken DOS 
!$O!MP CRITICAL
    !SG = matmul( SD, green ) !Not necessary
    ! computing total DOS
    !DOS=d_zero
    !AtomDOS=d_zero
    !do j=1,GetNAtoms()
    !  do i=LoAOrbNo(j),HiAOrbNo(j)
    !    AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
    !    DOS=DOS-dimag(SG(i,i))/d_pi
    !  end do
    !end do
    CH = 0.d0
    do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom)
       do j =1,NAOrbs 
          LocPD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          CH = CH + LocPD(ispin,i,j)*SD(j,i)
      end do
    end do
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_nested(.true.)
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------' 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
    !chunk=max(((n+1)/2)/omp_get_num_threads(),1)
     chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!!$OMP  DO SCHEDULE(DYNAMIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          !call gplus0(EE(k),greenn(k,:,:))
          call gplus(EE(k),greenn(k,:,:),gammar,gammal)
       end do
!!$OMP  END DO
!!$OMP  END PARALLEL
       !do i = 1,NAOrbs
       do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do

    end do
!$OMP END DO
!$OMP CRITICAL
       !do i = 1,NAOrbs
       do i=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
          do j = 1,NAOrbs
             LocPD(ispin,i,j)=LocPD(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    !do k=1,NAOrbs
    do k=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
       ! Non-orthogonal basis: Ch = Tr[P*SD]   
       do k1=1,NAOrbs
          CH = CH + LocPD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    !if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*(GetAtmChg(spinatom)).and.n.le.M) goto 1 
    ! Test for successfullness and integral final value
    M = 0
    !if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*(GetAtmChg(spinatom))) M = 1 
    CH = 16*CH/(3*(n+1))
    !do k=1,NAOrbs
    do k=LoAOrbNo(spinatom), HiAOrbNo(spinatom) 
       do l=1,NAOrbs
          LocPD(ispin,k,l) = 16*LocPD(ispin,k,l)/(3*(n+1))
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P has needed a max no. of loops=',(((n-1)/2)+1)/2

    return
  end subroutine localintpj

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------------------------------

!ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccc
  subroutine intch(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc 
    use omp_lib
   
    implicit none
    real,intent(in) :: Er,El
    integer,intent(inout) :: M
    real, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm 
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

      if(DebugDev)Write(*,'(A)')"ENTERED Device/intch"

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      call glesser(E0,green)
      do i=1,NAOrbs
       do j=1,NAOrbs
        PDOUT(ispin,i,j) = -ui*green(i,j)/(2*pi)
        p(i,j) = PDOUT(ispin,i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,NAOrbs
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*PDOUT(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        if(DebugDev)Write(*,'(A,I4)')"Device/intch l = ",l
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       call glesser(Em,greenm)
!!$OMP  SECTION
       call glesser(Ep,greenp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,NAOrbs
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDOUT(ispin,i,j)=PDOUT(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,NAOrbs
       do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,NAOrbs
      do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        PDOUT(ispin,i,j) = 16*PDOUT(ispin,i,j)/(3*(n+1))
        PDOUT(ispin,i,j) = PDOUT(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      write(ifu_log,'(A51,i4)')' Integration of POUT has needed a max no. of loops=',(((n-1)/2)+1)/2

      return
    end subroutine intch

!ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccc
  subroutine intchHW(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc
    use omp_lib

    implicit none
    real,intent(in) :: Er,El
    integer,intent(inout) :: M
    real, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP, HWP ! HWP ADDED FOR FOCK PULAY FORCES INTEGRATION.
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      call glesser(E0,green)
      do i=1,NAOrbs
       do j=1,NAOrbs
        PDOUT(ispin,i,j) = -ui*green(i,j)/(2*pi)
        HWOUT(ispin,i,j) = -DREAL(E0)*ui*green(i,j)/(2*pi)
        p(i,j) = PDOUT(ispin,i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,NAOrbs
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*PDOUT(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
       HWP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       call glesser(Em,greenm)
!!$OMP  SECTION
       call glesser(Ep,greenp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,NAOrbs
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
            HWP(i,j) = HWP(i,j)-ui*(DREAL(Em)*greenm(i,j)+DREAL(Ep)*greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDOUT(ispin,i,j)=PDOUT(ispin,i,j)+PDP(i,j)
             HWOUT(ispin,i,j)=HWOUT(ispin,i,j)+HWP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,NAOrbs
       do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,NAOrbs
      do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        PDOUT(ispin,i,j) = 16*PDOUT(ispin,i,j)/(3*(n+1))
        HWOUT(ispin,i,j) = 16*HWOUT(ispin,i,j)/(3*(n+1))
        PDOUT(ispin,i,j) = PDOUT(ispin,i,j)*(El-Er)/2
        HWOUT(ispin,i,j) = HWOUT(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      write(ifu_log,'(A51,i4)')' Integration of POUT/HWOUT has needed a max no. of loops=',(((n-1)/2)+1)/2

      return
    end subroutine intchHW

  !ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccc
  subroutine intchCGibbsY(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc
    use omp_lib

    implicit none
    real,intent(in) :: Er,El
    integer,intent(inout) :: M
    real, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP, HWP, CGibbsYP, CGibbsYPm, CGibbsYPp, CGibbsYPr  ! HWP ADDED FOR FOCK PULAY FORCES INTEGRATION.
    complex*16, dimension(NAOrbs,NAOrbs) :: CGibbsYKernel1P, CGibbsYKernel1Pm, CGibbsYKernel1Pp, CGibbsYKernel1Pr
    complex*16, dimension(NAOrbs,NAOrbs) :: CGibbsYKernel2P, CGibbsYKernel2Pm, CGibbsYKernel2Pp, CGibbsYKernel2Pr
!    complex*16, dimension(NAOrbs,NAOrbs) :: GibbsYL, GibbsYR
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

    write(ifu_log,'(A51,i4)')' ENTERING intchCGibbsY(Er,El,M)'

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
!      call glesser(E0,green) ! COMMENTED ON 2018-04-23 BECAUSE WITH NEW CompCGibbsY(E0,green,GibbsYL,GibbsYR) I ALREADY COMPUTE green.
!      call CompCGibbsY(greenm,CGibbsYPr)
      call CompCGibbsY(E0,green,CGibbsYPr,CGibbsYKernel1Pr,CGibbsYKernel2Pr) ! ADDED ON 2018-04-23 IT ALSO PROVIDES green BUT WITH THE ADDITIONAL GibbsYL, GibbsYR.
      do i=1,NAOrbs
       do j=1,NAOrbs
        PDOUTGIBBS(ispin,i,j) = -ui*green(i,j)/(2*pi)
!        YOUT(ispin,i,j) = -DREAL(E0)*ui*green(i,j)/(2*pi)
!        CGibbsY(ispin,i,j) = CGibbsYPr(i,j) ! COMMENTED ON 2018-04-23 BECAUSE WITH NEW CompCGibbsY(E0,green,GibbsYL,GibbsYR) I ALREADY COMPUTE green.
        CGibbsY(ispin,i,j) = CGibbsYPr(i,j) !GibbsYL(i,j)+GibbsYR(i,j)
        CGibbsYKernel1(ispin,i,j) = CGibbsYKernel1Pr(i,j)
        CGibbsYKernel2(ispin,i,j) = CGibbsYKernel2Pr(i,j)
        p(i,j) = PDOUTGIBBS(ispin,i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,NAOrbs
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*PDOUTGIBBS(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
       HWP=d_zero
       CGibbsYP=d_zero
       CGibbsYKernel1P=d_zero
       CGibbsYKernel2P=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
!       call glesser(Em,greenm) ! COMMENTED ON 2018-04-23 BECAUSE WITH NEW CompCGibbsY(E0,green,GibbsYL,GibbsYR) I ALREADY COMPUTE green.
       call CompCGibbsY(Em,greenm,CGibbsYPm,CGibbsYKernel1Pm,CGibbsYKernel2Pm)
!!$OMP  SECTION
!       call glesser(Ep,greenp) ! COMMENTED ON 2018-04-23 BECAUSE WITH NEW CompCGibbsY(E0,green,GibbsYL,GibbsYR) I ALREADY COMPUTE green.
       call CompCGibbsY(Ep,greenp,CGibbsYPp,CGibbsYKernel1Pp,CGibbsYKernel2Pp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
!       call CompCGibbsY(greenm,CGibbsYPm)
!       call CompCGibbsY(greenp,CGibbsYPp)
          do i=1,NAOrbs
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
!            HWP(i,j) = HWP(i,j)-ui*(DREAL(Em)*greenm(i,j)+DREAL(Ep)*greenp(i,j))*xs(l)**4/(2*pi)
!            CGibbsYP(i,j) = CGibbsYP(i,j) + CGibbsYPm(i,j) + CGibbsYPp(i,j)
            CGibbsYP(i,j) = CGibbsYP(i,j)-(CGibbsYPm(i,j)+CGibbsYPp(i,j))*xs(l)**4/(2*pi)
            CGibbsYKernel1P(i,j) = CGibbsYKernel1P(i,j)-(CGibbsYKernel1Pm(i,j)+CGibbsYKernel1Pp(i,j))*xs(l)**4/(2*pi)
            CGibbsYKernel2P(i,j) = CGibbsYKernel2P(i,j)-(CGibbsYKernel2Pm(i,j)+CGibbsYKernel2Pp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo

      enddo
!$OMP END DO
!$OMP CRITICAL
       ! FILLING PDOUT AND HWOUT WAS ALREADY DONE IN intchHW BUT NOW WE ARE FILLING PDOUTGIBBS AND CGibbsY.
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDOUTGIBBS(ispin,i,j)=PDOUTGIBBS(ispin,i,j)+PDP(i,j)
!             HWOUT(ispin,i,j)=HWOUT(ispin,i,j)+HWP(i,j)
             CGibbsY(ispin,i,j)=CGibbsY(ispin,i,j)+CGibbsYP(i,j)
             CGibbsYKernel1(ispin,i,j)=CGibbsYKernel1(ispin,i,j)+CGibbsYKernel1P(i,j)
             CGibbsYKernel2(ispin,i,j)=CGibbsYKernel2(ispin,i,j)+CGibbsYKernel2P(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,NAOrbs
       do j=1,NAOrbs
        rCHp=dble(PDOUTGIBBS(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUTGIBBS(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUTGIBBS(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUTGIBBS(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,NAOrbs
      do j=1,NAOrbs
        rCHp=dble(PDOUTGIBBS(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUTGIBBS(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUTGIBBS(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUTGIBBS(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        PDOUTGIBBS(ispin,i,j) = 16*PDOUTGIBBS(ispin,i,j)/(3*(n+1))
!        HWOUT(ispin,i,j) = 16*HWOUT(ispin,i,j)/(3*(n+1))
        CGibbsY(ispin,i,j) = 16*CGibbsY(ispin,i,j)/(3*(n+1))
        CGibbsYKernel1(ispin,i,j) = 16*CGibbsYKernel1(ispin,i,j)/(3*(n+1))
        CGibbsYKernel2(ispin,i,j) = 16*CGibbsYKernel2(ispin,i,j)/(3*(n+1))
        PDOUTGIBBS(ispin,i,j) = PDOUTGIBBS(ispin,i,j)*(El-Er)/2
!        HWOUT(ispin,i,j) = HWOUT(ispin,i,j)*(El-Er)/2
        CGibbsY(ispin,i,j) = CGibbsY(ispin,i,j)*(El-Er)/2
        CGibbsYKernel1(ispin,i,j) = CGibbsYKernel1(ispin,i,j)*(El-Er)/2
        CGibbsYKernel2(ispin,i,j) = CGibbsYKernel2(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      write(ifu_log,'(A51,i4)')' Integration of POUTGIBBS/CGibbsY has needed a max no. of loops=',(((n-1)/2)+1)/2

      return
    end subroutine intchCGibbsY

  !ccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !ccccccccccccccccccccccccccccccc
  subroutine intpj(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui
    use omp_lib

    implicit none

    real,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real, dimension(NAOrbs,NAOrbs) :: PDP

    real :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk
    real, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn 

    if(DebugDev)Write(*,'(A)')"ENTERED Device/intpj"

    !logical :: omp_get_nested
 
    !call omp_set_nested(.true.)
    !call mkl_set_dynamic(.false.)
    a = 1.d0/d_pi
    b = bi
    Em = Emi
    PD(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
    call gplus0(E0,green)
    CH = 0.d0
    do i = 1,NAOrbs
       do j =1,NAOrbs
          PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          CH = CH + PD(ispin,i,j)*SD(j,i)
       enddo
    enddo

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------' 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
    !chunk=max(((n+1)/2)/omp_get_num_threads(),1)
     chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!$OMP  DO SCHEDULE(STATIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          call gplus0(EE(k),greenn(k,:,:))
       end do
!$OMP  END DO
!$OMP  END PARALLEL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PD(ispin,i,j)=PD(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do k=1,NAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]   
       do k1=1,NAOrbs
          CH = CH + PD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    CH = 16*CH/(3*(n+1))
    do k=1,NAOrbs
       do l=1,NAOrbs
          PD(ispin,k,l) = 16*PD(ispin,k,l)/(3*(n+1))
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P has needed a max no. of loops=',(((n-1)/2)+1)/2

    !call mkl_set_dynamic(.true.)
    return
  end subroutine intpj

  subroutine intpj_SOC(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui
    use omp_lib

    implicit none

    real,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real, dimension(DNAOrbs,DNAOrbs) :: PDP

    real :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk !,omp_get_thread_num,omp_get_num_threads
    real, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(DNAOrbs,DNAOrbs) :: green
    complex*16, dimension(2,DNAOrbs,DNAOrbs) :: greenn 

    if(DebugDev)Write(*,'(A)')"ENTERED Device/intpj_SOC"
    !logical :: omp_get_nested
 
    a = 1.d0/d_pi
    b = bi
    Em = Emi
    PD_SOC = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
    call gplus0_SOC(E0,green)
    CH = 0.d0
    do i = 1,DNAOrbs
       do j =1,DNAOrbs
          PD_SOC(i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          CH = CH + PD_SOC(i,j)*S_SOC(j,i)
       enddo
    enddo

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_nested(.true.)
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------' 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
    !chunk=max(((n+1)/2)/omp_get_num_threads(),1)
     chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!$OMP  DO SCHEDULE(DYNAMIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          call gplus0_SOC(EE(k),greenn(k,:,:))
       end do
!$OMP  END DO
!$OMP  END PARALLEL
       do i = 1,DNAOrbs
          do j = 1,DNAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,DNAOrbs
          do j = 1,DNAOrbs
             PD_SOC(i,j)=PD_SOC(i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do k=1,DNAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]   
       do k1=1,DNAOrbs
          CH = CH + PD_SOC(k,k1)*S_SOC(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    CH = 16*CH/(3*(n+1))
    do k=1,DNAOrbs
       do l=1,DNAOrbs
          PD_SOC(k,l) = 16*PD_SOC(k,l)/(3*(n+1))
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P has needed a max no. of loops=',(((n-1)/2)+1)/2

    return
  end subroutine intpj_SOC


  !ccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !ccccccccccccccccccccccccccccccc
  subroutine intpjHW(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc
    use constants, only: d_pi, d_zero, ui
    use omp_lib

    implicit none

    real,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real, dimension(NAOrbs,NAOrbs) :: PDP, HWP ! HWP ADDED FOR FOCK PULAY FORCES INTEGRATION.

    real :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk
    real, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn

    if(DebugDev)Write(*,'(A)')"ENTERED Device/intpjHW"

    !logical :: omp_get_nested

    !call omp_set_nested(.true.)
    !call mkl_set_dynamic(.false.)
    a = 1.d0/d_pi
    b = bi ! bi is 0.0d0
    Em = Emi ! Emi is d_pi
    PD(ispin,:,:) = d_zero
    HW(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0) ! Em is d_pi inherited from Emi and b is 0.0d0 inherited from bi.
    ! er0 seems to be d_pi/2 because: edex3 = 0.5d0*((Em-b)*x + (Em+b)) = 0.5d0*((Em-0.0)*0.0 + (Em+0.0))
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq ! rrr = 0.5*abs(EMin)+10.0d0;
    call gplus0(E0,green)
    CH = 0.d0
    do i = 1,NAOrbs
       do j =1,NAOrbs
          PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          HW(ispin,i,j)=DREAL(E0)*PD(ispin,i,j)
          CH = CH + PD(ispin,i,j)*SD(j,i)
       enddo
    enddo

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
    HWP=d_zero
    !chunk=max(((n+1)/2)/omp_get_num_threads(),1)
     chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!$OMP  DO SCHEDULE(STATIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          call gplus0(EE(k),greenn(k,:,:))
       end do
!$OMP  END DO
!$OMP  END PARALLEL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
             HWP(i,j) = HWP(i,j) + a*(DREAL(EE(1))*dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +DREAL(EE(2))*dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PD(ispin,i,j)=PD(ispin,i,j)+PDP(i,j)
             HW(ispin,i,j)=HW(ispin,i,j)+HWP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do k=1,NAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]
       do k1=1,NAOrbs
          CH = CH + PD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    CH = 16*CH/(3*(n+1))
    do k=1,NAOrbs
       do l=1,NAOrbs
          PD(ispin,k,l) = 16*PD(ispin,k,l)/(3*(n+1))
          HW(ispin,k,l) = DREAL(E0)*PD(ispin,k,l)
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P/HW has needed a max no. of loops=',(((n-1)/2)+1)/2

    !call mkl_set_dynamic(.true.)
    return
  end subroutine intpjHW


  !ccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !ccccccccccccccccccccccccccccccc
  subroutine intpjCGibbsY(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc
    use constants, only: d_pi, d_zero, ui
    use omp_lib

    implicit none

    real,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real, dimension(NAOrbs,NAOrbs) :: PDP, HWP ! HWP ADDED FOR FOCK PULAY FORCES INTEGRATION.

    real :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk
    real, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn

    if(DebugDev)Write(*,'(A)')"ENTERED Device/intpjCGibbsY"

    !logical :: omp_get_nested

    !call omp_set_nested(.true.)
    !call mkl_set_dynamic(.false.)
    a = 1.d0/d_pi
    b = bi ! bi is 0.0d0
    Em = Emi ! Emi is d_pi
    PD(ispin,:,:) = d_zero
    HW(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0) ! Em is d_pi inherited from Emi and b is 0.0d0 inherited from bi.
    ! er0 seems to be d_pi/2 because: edex3 = 0.5d0*((Em-b)*x + (Em+b)) = 0.5d0*((Em-0.0)*0.0 + (Em+0.0))
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq ! rrr = 0.5*abs(EMin)+10.0d0;
    call gplus0(E0,green)
    CH = 0.d0
    do i = 1,NAOrbs
       do j =1,NAOrbs
          PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          HW(ispin,i,j)=DREAL(E0)*PD(ispin,i,j)
          !HW(ispin,i,j)=DREAL(-shift)*PD(ispin,i,j)
          CH = CH + PD(ispin,i,j)*SD(j,i)
       enddo
    enddo

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------'
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
    HWP=d_zero
    !chunk=max(((n+1)/2)/omp_get_num_threads(),1)
     chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!$OMP  DO SCHEDULE(STATIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          call gplus0(EE(k),greenn(k,:,:))
       end do
!$OMP  END DO
!$OMP  END PARALLEL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
             HWP(i,j) = HWP(i,j) + a*(DREAL(EE(1))*dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +DREAL(EE(2))*dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
             !HWP(i,j) = HWP(i,j) + a*(DREAL(-shift)*dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
             !     &   +DREAL(-shift)*dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PD(ispin,i,j)=PD(ispin,i,j)+PDP(i,j)
             HW(ispin,i,j)=HW(ispin,i,j)+HWP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do k=1,NAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]
       do k1=1,NAOrbs
          CH = CH + PD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    CH = 16*CH/(3*(n+1))
    do k=1,NAOrbs
       do l=1,NAOrbs
          PD(ispin,k,l) = 16*PD(ispin,k,l)/(3*(n+1))
          HW(ispin,k,l) = DREAL(E0)*PD(ispin,k,l)
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P/HW has needed a max no. of loops=',(((n-1)/2)+1)/2

    !call mkl_set_dynamic(.true.)
    return
  end subroutine intpjCGibbsY



  !****************************************************
  !* Compute retarded Green's function with SOC*
  !****************************************************
  subroutine gplus_SOC(z,green,gammar,gammal)
    use parameters, only: eta, glue
    use constants, only: c_zero, ui
   !use lapack_blas, only: zgetri, zgetrf
   !use lapack95, only: zgetri, zgetrf
   !use lapack95
   !use blas95

    implicit none
    external               zgetri, zgetrf

    integer :: i, j, info, omp_get_thread_num
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: gammar,gammal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigmar=c_zero
    sigmal=c_zero
    sigr1=-ui*eta*SD 
    sigl1=-ui*eta*SD 
    sigr2=-ui*eta*SD 
    sigl2=-ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1 )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2 )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    gammar=ui*(sigmar-conjg(transpose(sigmar)))
    gammal=ui*(sigmal-conjg(transpose(sigmal)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus_SOC

  !****************************************************
  !* Compute retarded Green's function with SOC*
  !****************************************************
  subroutine gplus0_SOC(z,green)
    use parameters, only: eta, glue
    use constants, only: c_zero, ui
   !use lapack_blas, only: zgetri, zgetrf
   !use lapack95, only: zgetri, zgetrf
   !use lapack95
   !use blas95

    implicit none
    external               zgetri, zgetrf

    integer :: i, j, info, omp_get_thread_num
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    
    ! Initilization 
    green=c_zero
    sigmar=c_zero
    sigmal=c_zero
    sigr1=-ui*eta*SD 
    sigl1=-ui*eta*SD 
    sigr2=-ui*eta*SD 
    sigl2=-ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1 )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2 )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus0_SOC


!ccccccccccccccccccccccccccccccc
! SOC subroutine
!ccccccccccccccccccccccccccccccc

 subroutine spin_orbit

    use parameters, only: soc, SOCEdit
    use cluster, only: NAOAtom
    use g09Common, only: GetNAtoms, GetAtmChg
    implicit none
    complex*16, dimension(DNAOrbs,DNAOrbs) :: hamil,hamil_SO
    integer :: i,j,aux,aux1,aux2,natoms,totdim
    integer, dimension(DNAOrbs) :: identity,orbital
    complex*16, dimension(8,8) :: reference

 totdim=NAOrbs   
 natoms=GetNAtoms()
 hamil = dcmplx(0.0d0,0.0d0)
 hamil_SO = dcmplx(0.0d0,0.0d0)
 reference=dcmplx(0.d0,0.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create the intra atomic SOC matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 reference(2,3)=dcmplx(0.d0,-soc/2.d0)
 reference(2,8)=dcmplx(soc/2.d0,0.d0)

 reference(3,2)=dcmplx(0.d0,soc/2.d0)
 reference(3,8)=dcmplx(0.d0,-soc/2.d0)
 
 reference(4,6)=dcmplx(-soc/2.d0,0.d0)
 reference(4,7)=dcmplx(0.d0,soc/2.d0)
 
 reference(6,4)=dcmplx(-soc/2.d0,0.d0)
 reference(6,7)=dcmplx(0.d0,soc/2.d0)

 reference(7,4)=dcmplx(0.d0,-soc/2.d0)
 reference(7,6)=dcmplx(0.d0,-soc/2.d0)

 reference(8,2)=dcmplx(soc/2.d0,0.d0)
 reference(8,3)=dcmplx(0.d0,soc/2.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Duplicate the size of the Hamiltonian to include up and down
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (NSpin == 2) then
    do i=1,NAOrbs
    do j=1,NAOrbs
       hamil(i,j)=dcmplx(HD(1,i,j),0.0d0)
       hamil(i+NAOrbs,j+NAOrbs)=dcmplx(HD(2,i,j),0.0d0)
    end do
    end do
 else 
    do i=1,NAOrbs
    do j=1,NAOrbs
       hamil(i,j)=dcmplx(HD(1,i,j),0.0d0)
       hamil(i+NAOrbs,j+NAOrbs)=dcmplx(HD(1,i,j),0.0d0)
    end do
    end do
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the ordering of the basis set  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 aux = 0
 do i=1, natoms
    do j=1,NAOAtom(i)
     aux=aux+1
     orbital(aux)=j
     orbital(aux+totdim)=j
     identity(aux)=i
     identity(aux+totdim)=i
  end do
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create the H_soc for the whole system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do i=1, 2*totdim
  do j=1,2*totdim
   aux1=orbital(i)
   aux2=orbital(j)
   if(identity(i).eq.identity(j) .and. SOCEdit(identity(i)) == 1) then
      if ( aux1 > 4  .or. aux2 > 4) print *, 'something wrong in SOC subroutine'
      hamil_SO(i,j)=reference(aux1,aux2)
   end if
  end do
 end do
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return the input matrix with double size plus the soc intection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do i=1, totdim*2
    do j=1, totdim*2
       H_SOC(i,j)=hamil(i,j)+hamil_SO(i,j)
    end do
 end do

 return
 end subroutine spin_orbit

  END MODULE device
