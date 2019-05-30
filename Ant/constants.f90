!*********************************************************!
!*********************  ANT.G-2.4.0  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   Juan Jose Palacios (1)                                !
!   David Jacob (2)                                       !
!                                                         !
!  (1) Departamento de Fisica de la Materia Condensada    !
!      Universidad Autonoma de Madrid                     !      
!      28049 Madrid (SPAIN)                               !
!  (2) Theory Department                                  !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!                                                         !
!*********************************************************!
  MODULE constants
!*********************************************************!
!  Module containing mathematical and physical constants  !
!*********************************************************!

  IMPLICIT NONE

  ! Mathematical constants
  real, PARAMETER :: d_pi   = 3.14159265358979323846d0
  real, PARAMETER :: d_zero = 0.0d0
  real, PARAMETER :: d_one  = 1.0d0
  
  complex*16, PARAMETER :: c_zero = (0.0d0,0.0d0)
  complex*16, PARAMETER :: c_one  = (1.0d0,0.0d0)
  complex*16, PARAMETER :: ui   = (0.0d0,1.0d0)

  ! Physical constants
  real, PARAMETER :: Hart = 27.2113834d0
  real, PARAMETER :: Ryd  = 0.5d0*Hart
  real, PARAMETER :: Bohr  = 0.5291772108d0

  real, PARAMETER :: eleccharge = 1.6d-19
!  real, PARAMETER :: hbar = 1.05457d-34 ! J*s
  real, PARAMETER :: hbar = 6.582d-16 ! eV*S ! WITH THIS VALUE THE GibbsY OPERATOR RESULTS HERMITIAN.

  ! Computational constants
  
  END MODULE constants
