!*********************************************************!
!*********************  ANT.G-2.4.1  *********************!
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
      MODULE preproc
!*********************************************************!
      IMPLICIT NONE
!DDEFMAXSHL=250000
!DDEFMAXATM=250000
!DDEFMAXTIT=100
!DDEFMAXRTE=4000
      integer, PARAMETER :: DEFMAXATM=250000
      integer, PARAMETER :: DEFMAXSHL=250000
      integer, PARAMETER :: DEFMAXTIT=100
      integer, PARAMETER :: DEFMAXRTE=4000

      integer, PARAMETER :: MAXATM=DEFMAXATM
      integer, PARAMETER :: MAXSHL=DEFMAXSHL
      integer, PARAMETER :: MAXTIT=DEFMAXTIT
      integer, PARAMETER :: MAXRTE=DEFMAXRTE

      END MODULE preproc
