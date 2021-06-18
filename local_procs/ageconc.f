!     !     Testing locally defined WAQ processes
! other modules have no underscores, so i'm sticking with that.
      
      subroutine ageconc ( pmsa   , fl     , ipoint , increm , noseg  ,
     &                      noflux , iexpnt , iknmrk , noq1   , noq2   ,
     &                      noq3   , noq4   )
!>\file
!>   Age concentration in the sense of CART

! Age concentration description:
! Name    T   L I/O   Description                                   Units
! ----    --- -  -    -------------------                            ----
! CONC    R*4 1 I concentration                                     [g/m3]
!     FL (1)  R*4 1 O age concentration flux                         [gN/m3/d]

! plus DELT

!     Old text from nitrif.f
!     Description of the module :
!
!        ----- old version -----
! Name    T   L I/O   Description                                   Units
! ----    --- -  -    -------------------                            ----
! CFL     R*4 1 I constant in O2FUNC                                   [-]
! CONC    R*4 1 I ammonium concentration                            [g/m3]
! COX     R*4 1 I critical oxygen concentratio for nitrification    [g/m3]
! CRTEMP  R*4 1 I critical temperature for nitrification              [oC]
! FL (1)  R*4 1 O nitrification flux                             [gN/m3/d]
! RC      R*4 1 I first order nitrification rate                     [1/d]
! O2FUNC  R*4 1 I function for OXY effect on the nitrification rate    [-]
! OOX     R*4 1 I critical concentr. dissolved oxygen               [g/m3]
! OXY     R*4 1 I concentration of dissolved oxygen                 [g/m3]
! POROS   R*4 1 L porosity                                             [-]
! SKEWN   R*4 1 I constant in O2FUNC                                   [-]
! TC      R*4 1 I temperature coefficient for nitrification            [-]
! TEMP    R*4 1 I ambient temperature                                 [oC]
! TEMP20  R*4 1 L ambient temperature - stand. temp (20)              [oC]
! ZERO    R*4 1 I zeroth order nitrification rate                [gN/m3/d]
!
!        ----- new version -----
! Name    T   L I/O   Description                                   Units
! ----    --- -  -    -------------------                            ----
! AMFUNC  R*4 1 I function for NH4 effect on the nitrification rate    [-]
! CRTEMP  R*4 1 I critical temperature for nitrification              [oC]
! CROXY   R*4 1 I critical oxygen concentratio for nitrification    [g/m3]
! FL (1)  R*4 1 O nitrification flux                             [gN/m3/d]
! K0NIT   R*4 1 I zeroth order nitrification rate                [gN/m3/d]
! K0TEMP  R*4 1 I zeroth order nitrification rate below CRTEMP   [gN/m3/d]
! K0NOX   R*4 1 I zeroth order nitrification rate below CROXY    [gN/m3/d]
! KNIT    R*4 1 I MM nitrification rate                          [gN/m3/d]
! KSAM    R*4 1 I half saturation constant for ammonium            [gN/m3]
! KSOX    R*4 1 I half saturation constant for oxygen               [g/m3]
! NH4     R*4 1 I ammonium concentration                            [g/m3]
! OXFUNC  R*4 1 I function for OXY effect on the nitrification rate    [-]
! OXY     R*4 1 I concentration of dissolved oxygen                 [g/m3]
! POROS   R*4 1 L porosity                                             [-]
! TC      R*4 1 I temperature coefficient for nitrification            [-]
! TEMP    R*4 1 I ambient temperature                                 [oC]
! TEMP20  R*4 1 L ambient temperature - stand. temp (20)              [oC]
!
!     Logical Units : -


!     original code has 26 items listed as input above.
!     dupe CRTEMP, NH4, OXY, TC, TEMP, TEMP20
!     Modules called : -

!     Name     Type   Library
!     ------   -----  ------------
!
      IMPLICIT NONE
!
      REAL     PMSA  ( * ) , FL    (*)
      INTEGER  IPOINT( * ) , INCREM(*) , NOSEG , NOFLUX,
     +         IEXPNT(4,*) , IKNMRK(*) , NOQ1, NOQ2, NOQ3, NOQ4
!
      INTEGER  IP1
      INTEGER  IN1
      INTEGER  IFLUX, ISEG
      REAL     CONC
      REAL     DELT
      REAL     FLNIT
      
! original code had 19 of these
      IN1  = INCREM( 1) ! CONC
      IP1  = IPOINT( 1) ! CONC
!
      IFLUX = 0
      DO 9000 ISEG = 1 , NOSEG
         !! CALL DHKMRK(1,IKNMRK(ISEG),IKMRK1)
         !! IF ( IKMRK1 .GT. 0) THEN

         ! How do I get BTEST?
         IF (BTEST(IKNMRK(ISEG),0)) THEN
            CONC   = MAX ( 0.0, PMSA(IP1 ) )
            DELT   = PMSA(IP2) 
            FL( 1 + IFLUX ) = CONC
            
         ENDIF

         IFLUX = IFLUX + NOFLUX
         IP1   = IP1   + IN1

 9000 CONTINUE

      RETURN
!
      END
