cas-------------------------------------------------------------------------
cas SUBROUTINES AND FUNCTIONS OF VPSC USED BY THE ADAPTATIVE SAMPLING METHOD
cas
cas                      Ricardo Lebensohn (Jun 2009)
cas
cas * THIS FILE CONTAINS:
cas
cas 1) SUBROUTINE vpsc_input_as, CALLED BY FUNCTION init_vpsc (OF MODULE 
cas    vpscVpETxt_interface),
cas 2) SUBROUTINE vpsc_evol_as, CALLED BY FUNCTION EvalEvolution_vpsc (OF 
cas    MODULE vpscVpETxt_interface),
cas 3) OTHER SUBROUTINES AND FUNCTIONS (NOT RELATED TO SECOND-ORDER METHOD), 
cas    IN ALPHABETIC ORDER
cas 4) SECOND-ORDER SUBROUTINES
cas
cas  * THE SUBROUTINES WHOSE NAMES END IN "_as" ARE EITHER NEW OR HAVE 
cas  SUBSTANTIAL CHANGES WRT THE STANDARD VPSC7 CODE
cas
cas  * MARKS "cas" INDICATE CHANGES WRT TO STANDARD VPSC7 CODE
cas
cas----------------------------------------------------------------------------
c
C     *******************************
C     SUBROUTINES AND FUNCTIONS LIST:
C     *******************************
C
C     SUBROUTINE VPSC_INPUT_AS
C     SUBROUTINE VPSC_EVOL_AS
C
C     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
C     SUBROUTINE CRYSTAL_SYMMETRY   --->   VERSION 04/SEP/06
c     SUBROUTINE DATA_CRYSTAL_AS
c     SUBROUTINE DATA_CRYSTAL_VOCE      --->      VERSION 12/JUN/2001
c     SUBROUTINE DATA_GRAIN_AS
c     FUNCTION DET
C     SUBROUTINE ESHELBY      --->      VERSION 07/OCT/05
C     SUBROUTINE ESH_INV3_VOIGT   --->   VERSION 23/JUL/01
C     SUBROUTINE ESH_INV4_VOIGT   --->   VERSION 20/JUL/01
C     SUBROUTINE ESH_MULT_VOIGT
C     SUBROUTINE EULER
C     SUBROUTINE GET_TWINOR_AS
C     SUBROUTINE GRAIN_RATE_AND_MODULI (ex MICRO) -->  VERSION 17/DEC/02
C     SUBROUTINE GRAIN_STRESS (ex VISC)   --->   version of 19/MAR/2003
C     SUBROUTINE INIT_CRSS_VOCE_AS
C     SUBROUTINE INITIAL_STATE_GUESS     --->     VERSION 07/DEC/05
c     SUBROUTINE LOAD_CONDITIONS_AS
C     SUBROUTINE LU_INVERSE
C     SUBROUTINE LU_EQSYSTEM
C     SUBROUTINE NEWTON_RAPHSON (ex SNLNR)   --->   VERSION APR/2003
C     SUBROUTINE SCALE_3   ---->   VERSION 26/MAY/2000
C     SUBROUTINE STATE_5x5      --->      VERSION 26/jan/2000
C     SUBROUTINE STATE_6x6      --->      VERSION 16/NOV/97
C     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
C     FUNCTION TMISMATCH5   ---->   VERSION OF 02/JAN/03
C     FUNCTION TNORM5   ---->   VERSION OF 02/JAN/03
C     SUBROUTINE TWIN_ORIENTATION    --->    VERSION OF 04/dec/02
C     SUBROUTINE UPDATE_CRSS_VOCE_AS
C     SUBROUTINE UPDATE_F2IJDOT_AS
C     SUBROUTINE UPDATE_ORIENTATION_AS
C     SUBROUTINE UPDATE_SCHMID
C     SUBROUTINE UPDATE_SHAPE_AS
C     SUBROUTINE UPDATE_TWINNING_AS
C     SUBROUTINE VOIGT   ---->   VERSION OF 09/02/98
C     SUBROUTINE VPSC     --->      VERSION 14/SEP/04
C
C     FLUCTUATIONS AND SECOND-ORDER SUBROUTINES:
C
C     SUBROUTINE FLUCTUATIONS
C     SUBROUTINE LINSOLVER25
C     SUBROUTINE GET_THEFLU
C     SUBROUTINE SOP
C     SUBROUTINE SOMOD
C     SUBROUTINE GRAIN_STRESS_ALT
C     SUBROUTINE EXTRAPOLSO
C     SUBROUTINE VOIGT10
C     SUBROUTINE GET_GAMDOT
C
C------------------------------------------------------------------------------
c
cas
cas   SUBROUTINE VPSC_INPUT
cas
      subroutine vpsc_input_as 
     #   (vpscin, orients, itnphase, ngrains, volfrac, ! INPUT
!RL
     #    nstrengthsMax,itngrtot,                      ! INPUT

!RL     #    strengths, shapes, ntwinsys, porosity)       ! OUTPUT
     #    strengths, nstrengths, shapes, ntwinsys, porosity,        ! OUTPUT
     &      c_scaling)                                ! input 
cas
      INCLUDE 'vpsc_as.dim'
      DIMENSION AUX5(5),AUX55(5,5),AUX3(3),AUX33(3,3),AUX3333(3,3,3,3)
cas
      CHARACTER*255 vpscin

cRL      dimension strengths(:), shapes(:,:), orients(:,:),
cRL     #          ngrains(:),volfrac(:),ntwinsys(:)
cass      dimension strengths(:,:), shapes(:,:), orients(:,:),
cass     #          ngrains(:),volfrac(:),ntwinsys(:),nstrengths(:)
cass
      dimension strengths(nstrengthsMax,itnphase)
      dimension shapes(6,itnphase+1), orients(3,itngrtot)
      dimension ngrains(itnphase),volfrac(itnphase)
      dimension ntwinsys(itnphase),nstrengths(itnphase)

      real*8 c_scaling

      character(len=255) :: dataDir

      call get_environment_variable("VPSC_INPUT_PATH", dataDir)
      if (dataDir == "") then 
         dataDir = '../../CoEVP/CM/src/fine_scale_models/tantalum/'
      endif

cass
      UR0= 20     ! VPSCIN 
      UR1= 21     ! FILECRYS

      OPEN(UR0,FILE=vpscin,STATUS='OLD')
      
cass      OPEN(UR0,FILE='vpsc_as_try.in',STATUS='OLD')
cas

C *********   INITIALIZATION BLOCK   ***************************

      PI=4.*ATAN(1.)

      DO I=1,5
      DO J=1,5
        XID5(I,J)=0.
        IF(I.EQ.J) XID5(I,J)=1.
      ENDDO
      ENDDO

      DO I=1,3
      DO J=1,3
        ZERO33(I,J)=0.
        XID3(I,J)=0.
        IF(I.EQ.J) XID3(I,J)=1.
      ENDDO
      ENDDO

C     CALCULATES TENSORS OF THE SYMMETRIC BASIS 'B(3,3,6)'
      CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,0,6)

C     SEED FOR RANDOM NUMBER GENERATOR (RAN2) (USED FOR TWINNING AND RX)
      JRAN=-1

C     READS # OF ELEMENTS, # OF PHASES AND RELATIVE VOLUME OF PHASE IN ELEM.
C     PHASES ARE LABELED SEQUENTIALLY FOR A MULTIELEMENT RUN (i.e. IF ELEMENT
C     #1 CONSISTS OF phase1=fcc AND phase2=bcc, THEN phase3 AND phase4 ARE
C     THE fcc AND bcc PHASES IN ELEMENT #2...AND SO ON)

cas      READ(UR0,*) NELEM
      nelem=1
      READ(UR0,*) NPH
      READ(UR0,*) (WPH(I),I=1,NPH)
cas
cas   check consistency between nph, wph with itphase, volfrac (from AS)
cas
      if(nph.ne.itnphase) then
      write(*,*) 'INCONSISTENCY, NPH'
      stop
      endif
c
      do iph=1,nph
      if(wph(iph).ne.volfrac(iph)) then
      write(*,*) 'INCONSISTENCY, WPH - PH#',iph
      write(*,*) 'WPH=',wph(iph)
      write(*,*) 'VOLFRAC=',volfrac(iph)
      stop
      endif
      enddo
cas
      if(nph .gt. nphmx) then
        write(*,'('' number of phases exceeds multiphase dimens !!'')')
        write(*,'('' --> increase parameter NPHMX to'',i5)') nph
        stop
      endif

C ***************************************************************************
C     LOOP OVER PHASES.
C     FOR EACH PHASE READS CRYSTAL SYSTEMS, TEXTURE AND GRAIN SHAPES.
C ***************************************************************************

      IPHBOT=1
      IPHTOP=NPH
      NGR(0)=0
      ISHAPE(0)=0

C *** INITIALIZE DEFORMATION GRADIENT (AND ELLIPSOID) OF THE PHASE

      DO I=1,3
      DO J=1,3
ca        FIJPH(I,J,0)=XID3(I,J)
        F2IJPH(I,J,0)=XID3(I,J)
      ENDDO
      ENDDO
cas      CALL UPDATE_SHAPE (0)
      CALL UPDATE_SHAPE_AS (0)

      DO IPH=1,NPH

        READ(UR0,'(a)') PROSA
        READ(UR0,*)     ISHAPE(IPH)
cas
        if(ishape(iph).ne.0) then
          write(*,*) 'ISHAPE(',IPH,').NE.0'
          stop
        endif
cas
        READ(UR0,*)     (AXISPH(0,I,IPH),I=1,3)
        READ(UR0,*)     (EULERPH(I,IPH) ,I=1,3)

        READ(UR0,'(a)') PROSA

        READ(UR0,'(a)') FILECRYS

        filecrys = trim(dataDir)//filecrys


C *** READS SLIP AND TWINNING MODES FOR THE PHASE

        OPEN (UR1,FILE=FILECRYS,STATUS='OLD')
cas
cas          CALL DATA_CRYSTAL(IPH)
cas
          CALL DATA_CRYSTAL_AS (IPH,nstrengthsMax)  
cas
      CLOSE(UR1)     ! closes SX file after reading hardening parameters
cas
casC *** READS INITIAL TEXTURE (BUNGE.ROE,KOCKS CONVENTIONS) FROM 'FILETEXT'
casC *** INITIALIZES DEFORMATION GRADIENT AND ELLIPSOID FOR EACH PHASE.
casC *** IF ISHAPE>2 READS INITIAL GRAIN AXES AND ORIENTATION FROM FILEAXES
casC     AND INITIALIZES DEF GRADIENT & ELLIPSOID FOR INDIVIDUAL GRAINS.
cas
cas          CALL DATA_GRAIN(IPH)
cas
cass          CALL DATA_GRAIN_AS (IPH,NGRAINS(IPH),ORIENTS)
cass
          CALL DATA_GRAIN_AS (IPH,NGRAINS(IPH),ITNGRTOT,ORIENTS)
cass
cas
        IFLAT(IPH)=0      ! WILL UPDATE THE SHAPE OF THE PHASE

cas        CALL UPDATE_SHAPE (IPH)
        CALL UPDATE_SHAPE_AS (IPH)

      ENDDO     ! END OF DATA INPUT LOOP OVER ALL PHASES

C ***************************************************************************

      NGTOT=NELEM*NGR(NPH)
!      WRITE(*,'('' --> TOTAL NUMBER OF GRAINS IS'',I6)') NGTOT
      IF(NGTOT.GT.NGRMX) THEN
        WRITE(*,'('' --> INCREASE PARAMETER NGRMX IN vpsc_as.dim'')')
        STOP
      ENDIF

C ****************************************************************************
C *** READS SETTINGS FOR CONVERGENCE PROCEDURES.

      READ(UR0,'(A)') PROSA
      READ(UR0,*) ERRS,ERRD,ERRM,ERRSO
      READ(UR0,*) ITMAXEXT,ITMAXINT,ITMAXSO
      READ(UR0,*) IRSVAR,JXRSINI,JXRSFIN,JXRSTEP
      READ(UR0,*) IBCINV

!     modify the convergence criterion 
      errm = errm*c_scaling

C *** INITIALIZE GAUSS-LEGENDRE COORDINATES AND ARRAYS USED IN THE DOUBLE
C *** INTEGRALS GIVING THE ESHELBY TENSORS

      CALL ESHELBY(AUX3,AUX3333,0.0+0,AUX3333,AUX3333,AUX33,AUX33,PDIL,
     #                  AUX3333,AUX3333,0)


C ****************************************************************************
C *** READS MODELING CONDITIONS.

      READ(UR0,'(A)') PROSA
cas      READ(UR0,*) IHARDLAW
      ihardlaw=0
cas
      READ(UR0,*) IRATESENS
      READ(UR0,*) INTERACTION
      READ(UR0,*) IUPDORI, IUPDSHP, IUPDHAR

C *** CHECKS IF RATE SENSITIVITY IS THE SAME FOR ALL SYSTEMS.
C *** SEARCH FOR NRSMIN (NEEDED TO GET TAUMAX INSIDE NR SUBROUTINE)

      NUNIQUE=1
      NCOMPA=NRS(1,1)
      NRSMIN=NRS(1,1)
      DO IPH=1,NPH
        DO IS=1,NSYST(IPH)
          IF(NRS(IS,IPH).NE.NCOMPA) NUNIQUE=0
          IF(NRS(IS,IPH).LT.NRSMIN) NRSMIN=NRS(IS,IPH)
        ENDDO
      ENDDO
      IF(NUNIQUE.EQ.0) THEN
        WRITE(*,'('' --> THIS RUN CONTAINS MIXED nrs EXPONENTS !!'')')
        IF(INTERACTION.NE.1) THEN
          WRITE(*,'('' --> WILL FORCE AFFINE (INTERACTION=1)'')')
          STOP
        ENDIF
      ENDIF

      iflu=0
      if(interaction.eq.5) iflu=1

      IF(IUPDORI.EQ.0 .AND. IUPDSHP.EQ.1) THEN
        WRITE(*,'('' --> IUPDSHP=1 AND IUPDORI=0'')')
        WRITE(*,'('' --> CANNOT UPDATE SHAPE WITHOUT ALSO UPDATING'',
     #            '' ORIENTATION'')')
        STOP
      ENDIF

C *******************************************************************
C *** INITIALIZE ARRAYS ASSOCIATED WITH GRAINS:
C     HARDENING, ACCUMULATED SHEAR, POWER, TWINNING PARAMETERS, etc

      DO IPH=1,NPH

        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          GTOTGR(KKK)=0.
          GAMD0G(KKK)=1.0
          KTWSMX(KKK)=0
          NTWEVENTS(KKK)=0
          DO ITS=1,NTWSYS(IPH)
            TWFRSY(ITS,KKK)=0.
          ENDDO
        ENDDO

        DO ITM=1,NTWMOD(IPH)
          EFTWFR(ITM,IPH)=0.
          TWFRPH(ITM,IPH)=0.
        ENDDO

      ENDDO

C *** INITIALIZE CRSS IN EACH SYSTEM OF EACH GRAIN.
C     FOR 'VOCE' THE PARAMETERS HAVE BEEN READ FROM SX FILE.
C     FOR 'MTS' READ HARDENING PARAMETERS HERE (OPTION=0). WILL ALSO
C     INITIALIZE THE STRUCTURE/TEMP/RATE INDEPENDENT ARRAY 'TAUE'.

cas      IF(IHARDLAW.EQ.0) CALL UPDATE_CRSS_VOCE (1)
      IF(IHARDLAW.EQ.0) CALL INIT_CRSS_VOCE_AS
cas
cas   output to init_vpsc
cas
      do iph=1,nph
cas
      nstrengths(iph)=nmodes(iph)
cas
cas      strengths(iph)=tau(1,0,iph)
cas
      is=1
      do im=1,nmodes(iph)
        if(im.gt.1) is=is+nsm(im-1,iph)
        strengths(im,iph)=tau(is,0,iph)
      enddo
cas
ca      do i=1,3
ca      do j=1,3
ca      aux33(i,j)=0.
ca      do k=1,3
ca      aux33(i,j)=aux33(i,j)+fijph(i,k,iph)*fijph(j,k,iph)
ca      enddo
ca      enddo
ca      enddo
c
      shapes(1,iph)=f2ijph(1,1,iph)
      shapes(2,iph)=f2ijph(2,2,iph)
      shapes(3,iph)=f2ijph(3,3,iph)
      shapes(4,iph)=f2ijph(2,3,iph)
      shapes(5,iph)=f2ijph(3,1,iph)
      shapes(6,iph)=f2ijph(1,2,iph)
c
      ntwinsys(iph)=ntwsys(iph)
c
      porosity=0.d0  ! for the time being
c
      enddo
      close(UR0)
cas
      RETURN
      END
c
c------------------------------------------------------------------------------
c
      subroutine vpsc_evol_as
     # (stressSvecP,wgts,strengths,gammas,shapes,porosity,    ! input
cass
     # itnPhase, itnGrTot, itnStrngthMax, itnTwinSysMax,      ! input
cass
     #  velgrad,ratesens,strengthRates,gammarates,shapeRates, ! output
     #  porosityRate,reorRates,twVFRates)                     ! output

      INCLUDE 'vpsc_as.dim'

cass      dimension stressSvecP(:), velgrad(:,:), wgts(:)
cRL      dimension strengths(:), gammas(:),shapes(:,:)
cRL      dimension strengthRates(:), gammarates(:), shapeRates(:,:)
cass      dimension strengths(:,:), gammas(:),shapes(:,:)
cass      dimension strengthRates(:,:), gammarates(:), shapeRates(:,:)
cass
      dimension stressSvecP(7), velgrad(3,3), wgts(itnGrTot)
      dimension strengths(itnStrngthMax,itnphase)
      dimension gammas(itnphase),shapes(6,itnphase+1)
      dimension strengthRates(itnStrngthMax,itnphase)
      dimension gammarates(itnphase), shapeRates(6,itnphase)
      dimension reorRates(3,itnGrTot)
      dimension twVFRates(itnTwinSysMax,itnGrTot)
cass
      DIMENSION AUX5(5),AUX33(3,3),AUX55(5,5),AUX3333(3,3,3,3)
cas
cass      write(*,*) 'itnStrngthMax', itnStrngthMax
cass      write(*,*) 'strengths(1,1)',strengths(1,1)
cas
cass      write(*,*) strengths
cass      pause
cass
cass  provisory, improve later
cass
      do i=1,itngrtot
      wgt(i)=wgts(i)
      enddo
cass
      scauchy(1,1)=stressSvecP(1)+stressSvecP(7)
      scauchy(2,2)=stressSvecP(2)+stressSvecP(7)
      scauchy(3,3)=stressSvecP(3)+stressSvecP(7)
      scauchy(2,3)=stressSvecP(4)
      scauchy(3,2)=scauchy(2,3)
      scauchy(3,1)=stressSvecP(5)
      scauchy(1,3)=scauchy(3,1)
      scauchy(1,2)=stressSvecP(6)
      scauchy(2,1)=scauchy(1,2)
c
      CALL CHG_BASIS (SBAR,SCAUCHY,AUX55,AUX3333,2,5)
cas
      CALL UPDATE_CRSS_VOCE_AS (1,STRENGTHS,GAMMAS,
     #                       STRENGTHRATES,GAMMARATES
cass
     #                       ,itnStrngthMax,itnphase)
cass
cas
      do iph=1,nph
      f2ijph(1,1,iph)=shapes(1,iph)
      f2ijph(2,2,iph)=shapes(2,iph)
      f2ijph(3,3,iph)=shapes(3,iph)
      f2ijph(2,3,iph)=shapes(4,iph)
      f2ijph(3,1,iph)=shapes(5,iph)
      f2ijph(1,2,iph)=shapes(6,iph)
      call update_shape_as(iph)
      enddo
cas     
cass        WRITE(*,'(''   ITSGR    SIGGR      SIGAV        DAV    ITTAN'',
cass     #               ''     MTAN    ITSEC     MSEC'',/)')
cass
!        WRITE(*,'(''   ITSGR    SIGGR      SIGAV        DAV'')')
cass

C     IRATESENS=0 ALLOWS TO ELIMINATE RATE SENSITIVITY.

          IF(IRATESENS.EQ.0) CALL SCALE_3 (ISTEP)   ! rate-insensitive response

C     IMPOSING STRESS -> MAKES A SACHS GUESS FOR EVERY STEP.

          CALL INITIAL_STATE_GUESS

cas        IF(ISTEP.EQ.1.AND.IRECOVER.EQ.1) THEN
cas          KGX=1
cas          DO IPH=IPHBOT,IPHTOP
cas            IPHEL=IPH-IPHBOT+1
cas          DO KKK=NGR(IPH-1)+1,NGR(IPH)
cas            CALL GRAIN_RATE_AND_MODULI (1,1,KGX,KKK,IPHEL,IPH)
cas            KGX=KGX+1
cas          ENDDO
cas          ENDDO
cas        ENDIF

C *************************************************************************
C     VPSC CALCULATION (OR TAYLOR CALCULATION WHEN INTERACTION=0) OF
C     STRESS AND STRAIN-RATE FOR EVERY GRAIN.

cas        CALL VPSC (ISTEP)      !!!     THIS IS THE CORE OF THE CODE
        CALL VPSC      !!!     THIS IS THE CORE OF THE CODE

casC ***************************************************************************
casC     VPSC PROVIDES THE MACROSCOPIC STRAIN RATE 'DSIM'.
casC     CALCULATES MACROSCOPIC ROTATION RATE AND MACROSCOPIC VELOCITY GRADIENT
casC     DEPENDING ON THE BOUNDARY CONDITIONS.
cas
cas      DO I=1,3
cas        ROTBAR(I,I)=0.
cas        DO J=1,3
cas          IF(IUDOT(I,J).EQ.1.AND.IUDOT(J,I).EQ.1) THEN
cas            ROTBAR(I,J)=(UDOT(I,J)-UDOT(J,I))/2.
cas          ELSE IF(IUDOT(I,J).EQ.1) THEN
cas            ROTBAR(I,J)=UDOT(I,J)-DSIM(I,J)
cas            ROTBAR(J,I)=-ROTBAR(I,J)
cas          ENDIF
cas          UDOT(I,J)=DSIM(I,J)+ROTBAR(I,J)
cas        ENDDO
cas      ENDDO

C     VON MISES STRAIN-RATE & STRESS.
C     IF INTERACT=0: DBAR IS IMPOSED AND WE DEFINE SBAR=SAV.
C     IF INTERACT>0: DBAR & SBAR FOLLOW FROM THE CALL TO SUBR. STATE6x6.

cas      SVM=0.
cas      DVM=0.
cas      DO I=1,5
cas        SVM=SVM+SBAR(I)*SBAR(I)
cas        DVM=DVM+SBAR(I)*DBAR(I)
cas      ENDDO
cas      SVM=SQRT(SVM*3./2.)
cas      DVM=DVM/SVM

cas      CALL STAT_SHEAR_ACTIVITY
cas      IF(IHARDLAW.NE.2) CALL STAT_GRAIN_SHAPE
cas      IF(IHARDLAW.NE.2) CALL STAT_STRESS_STRAIN

cas        CALL UPDATE_FIJ(0)        ! UPDATES DEFORM TENSOR OF ELEMENT
cas
           CALL UPDATE_F2IJDOT_AS (0) 
cas
cas        CALL UPDATE_SHAPE(0)      ! UPDATES SHAPE OF ELEMENT
cas
cas          DO IPH=IPHBOT,IPHTOP
cas            IPHEL=IPH-IPHBOT+1
cas            IF(NTWMOD(IPHEL).NE.0) CALL UPDATE_TWINNING (IPH)
cas          ENDDO
cas
         IF(itnTwinSysMax.gt.0) then
cass           CALL UPDATE_TWINNING_AS (twVFRates)
           CALL UPDATE_TWINNING_AS (twVFRates,itnTwinSysMax,itnGrTot)
         ENDIF
cas
cas          IF(IUPDORI.EQ.1) CALL UPDATE_ORIENTATION
cas
      IF(IUPDORI.EQ.1) THEN
cass         CALL UPDATE_ORIENTATION_AS (REORRATES)
         CALL UPDATE_ORIENTATION_AS (REORRATES,itnGrTot)
       ELSE
          do I=1,3
             do J=1,kgrtot
                REORRATES(I,J)=0.
             enddo
          enddo
          do I=1,3
             do J=1,3
                ROTBAR_AS(I,J)=0.
             enddo
          enddo
       ENDIF 

      DO I=1,3
        DO J=1,3
          VELGRAD(I,J)=DSIM(I,J)+ROTBAR_AS(I,J)
        ENDDO
      ENDDO
cas
          DO IPH=IPHBOT,IPHTOP
cas
cas            CALL UPDATE_FIJ(IPH)    ! UPDATES DEF TENSOR OF PHASE & GRAINS
cas
            CALL UPDATE_F2IJDOT_AS (IPH)
cas
      shaperates(1,iph)=f2ijdotph(1,1,iph)
      shaperates(2,iph)=f2ijdotph(2,2,iph)
      shaperates(3,iph)=f2ijdotph(3,3,iph)
      shaperates(4,iph)=f2ijdotph(2,3,iph)
      shaperates(5,iph)=f2ijdotph(3,1,iph)
      shaperates(6,iph)=f2ijdotph(1,2,iph)
cas
cas            IF(IUPDSHP.EQ.1) CALL UPDATE_SHAPE (IPH)
cas
          ENDDO

          IF(IUPDHAR.EQ.1) THEN
casC *** VOCE HARDENING PLUS PREDOMINANT TWIN REORIENTATION SCHEME
cas            IF(IHARDLAW.EQ.0)  CALL UPDATE_CRSS_VOCE (2)
cas
      IF(IHARDLAW.EQ.0)  THEN
      CALL UPDATE_CRSS_VOCE_AS (2,STRENGTHS,GAMMAS,
     #                       STRENGTHRATES,GAMMARATES
cass
     #                       ,itnStrngthMax,itnphase)
cass
      ENDIF
cas
          ENDIF
cas
      ratesens=1./nrs(1,1)
cas
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C **************************************************************************
C     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
C
C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
C     (modif. 10/MAY/01 - KDIM version - R.L.)
C
C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
C     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
C     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
C             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
C     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
C             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
C **************************************************************************

      SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      PARAMETER (SQR2=1.41421356237309   )
      PARAMETER (RSQ2=0.70710678118654744)
      PARAMETER (RSQ3=0.57735026918962584)
      PARAMETER (RSQ6=0.40824829046386304)

      DIMENSION CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3)

C     DIMENSION B(3,3,6)
C     DATA B /RSQ6,0,   0,   0,   RSQ6,0,   0,   0,  -2*RSQ6,
C    #        RSQ2,0,   0,   0,  -RSQ2,0,   0,   0,   0,
C    #        0,   0,   0,   0,   0,   RSQ2,0,   RSQ2,0,
C    #        0,   0,   RSQ2,0,   0,   0,   RSQ2,0,   0,
C    #        0,   RSQ2,0,   RSQ2,0,   0,   0,   0,   0,
C    #        RSQ3,0,   0,   0,   RSQ3,0,   0,   0,   RSQ3/

      COMMON/BASIS/ B(3,3,6)

      IF(IOPT.EQ.0) THEN
C *** CALCULATES BASIS TENSORS B(N)

        DO I=1,3
          DO J=1,3
            DO N=1,6
              B(I,J,N)=0.0
            ENDDO
          ENDDO
        ENDDO

        B(1,1,2)=-RSQ6
        B(2,2,2)=-RSQ6
        B(3,3,2)= 2.D0*RSQ6

        B(1,1,1)=-RSQ2
        B(2,2,1)= RSQ2

        B(2,3,3)=RSQ2
        B(3,2,3)=RSQ2

        B(1,3,4)=RSQ2
        B(3,1,4)=RSQ2

        B(1,2,5)=RSQ2
        B(2,1,5)=RSQ2

        B(1,1,6)=RSQ3
        B(2,2,6)=RSQ3
        B(3,3,6)=RSQ3

      ENDIF

C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
      IF(IOPT.EQ.1) THEN
        DO 40 I=1,3
        DO 40 J=1,3
        C2(I,J)=0.0
        DO 40 N=1,KDIM
   40   C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
      ENDIF

C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
      IF(IOPT.EQ.2) THEN
        DO 50 N=1,KDIM
        CE2(N)=0.0
        DO 50 I=1,3
        DO 50 J=1,3
   50   CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
      ENDIF

C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
      IF(IOPT.EQ.3) THEN
        DO 20 I=1,3
        DO 20 J=1,3
        DO 20 K=1,3
        DO 20 L=1,3
        C4(I,J,K,L)=0.0
        DO 20 N=1,KDIM
        DO 20 M=1,KDIM
   20   C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
      ENDIF

C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
      IF(IOPT.EQ.4) THEN
        DO 30 N=1,KDIM
        DO 30 M=1,KDIM
        CE4(N,M)=0.0
        DO 30 I=1,3
        DO 30 J=1,3
        DO 30 K=1,3
        DO 30 L=1,3
   30   CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
      ENDIF

      RETURN
      END
c
c------------------------------------------------------------------------------
c
c ***********************************************************************
c     SUBROUTINE CRYSTAL_SYMMETRY   --->   version 04/SEP/06
c
c *** If IOPTION=1:
c     Reads crystal symmetry 'icrysym' and unit cell parameters.
c     Generates vectors 'cvec(i,n)' of the unit cell.
c     Generates symmetry operators 'h(i,j,nsymop)' for crystal symmetry.
c *** If IOPTION=2:
c     Reads Miller indices of systems in 3 or 4-index notation 'isn(i)'
c     & 'isb(i)'. Calculates normal & burgers vectors 'sn(i)' & 'sb(i)'
c *** If IOPTION=3:
c     Reads Miller indices of diffraction planes 'isn(i)' and spherical
c     angles 'chi , eta'of diffraction direction.
c     Generates crystallographically equivalent orientations sneq(i,n) of
c     a sn(i) by applying all the symmetry operations to it.
c     Discards multiplicity and defines 'npol'
c *** Simmetry parameter ICRYSYM:
c        1: CUBIC
c        2: HEXAGONAL
c        3: TRIGONAL
c        4: TETRAGONAL
c        5: ORTHORHOMBIC
c        6: MONOCLINIC
c        7: TRICLINIC
c ***********************************************************************

      subroutine crystal_symmetry (ioption,ur1,icrysym,sn,sneq,sb,npol)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      dimension h(3,3,24),hx(3,3,6),itag(24)
      dimension isn(4),sn(3),sneq(3,24)
      dimension isb(4),sb(3)
      dimension cdim(3),cang(3),cvec(3,3)
      integer ur1
      character*5 crysym
      save h,nsymop,cvec
      pi=3.1415926535898

      !write(*,*) 'In crystal_summetry ur1 = ', ur1
c ****************************************************************************

      if(ioption.eq.1) then

        read(ur1,*)
        read(ur1,'(a)') crysym
        icrysym=0
        if(crysym.eq.'cubic' .or. crysym.eq.'CUBIC') icrysym=1
        if(crysym.eq.'hexag' .or. crysym.eq.'HEXAG') icrysym=2
        if(crysym.eq.'trigo' .or. crysym.eq.'TRIGO') icrysym=3
        if(crysym.eq.'tetra' .or. crysym.eq.'TETRA') icrysym=4
        if(crysym.eq.'ortho' .or. crysym.eq.'ORTHO') icrysym=5
        if(crysym.eq.'monoc' .or. crysym.eq.'MONOC') icrysym=6
        if(crysym.eq.'tricl' .or. crysym.eq.'TRICL') icrysym=7
        if(icrysym.eq.0) then
          write(*,*) ' *** CANNOT RECOGNIZE THE CRYSTAL SYMMETRY !!'
          stop
        endif

        READ(UR1,*) (CDIM(i),i=1,3),(CANG(i),i=1,3)
        DO I=1,3
          CANG(I)=CANG(I)*PI/180.
        ENDDO
        CVEC(1,1)=1.
        CVEC(2,1)=0.
        CVEC(3,1)=0.
        CVEC(1,2)=COS(CANG(3))
        CVEC(2,2)=SIN(CANG(3))
        CVEC(3,2)=0.
        CVEC(1,3)=COS(CANG(2))
        CVEC(2,3)=(COS(CANG(1))-COS(CANG(2))*COS(CANG(3)))/SIN(CANG(3))
        CVEC(3,3)=SQRT(1.-CVEC(1,3)**2-CVEC(2,3)**2)
        DO J=1,3
        DO I=1,3
          CVEC(I,J)=CDIM(J)*CVEC(I,J)
        ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          DO M=1,6
            HX(I,J,M)=0.d0
          ENDDO
          DO N=1,24
            H(I,J,N)=0.d0
          ENDDO
        ENDDO
        ENDDO

c *** identity operation ---> triclinic & all symmetries
      do i=1,3
        h(i,i,1)=1.d0
      enddo
      nsymop=1

c *** 180 deg rotation around (001) ---> orthorhombic, monoclinic
      if(icrysym.eq.5 .or. icrysym.eq.6) then
        h(1,1,2)= cos(pi)
        h(2,2,2)= cos(pi)
        h(3,3,2)= 1.d0
        h(1,2,2)=-sin(pi)
        h(2,1,2)= sin(pi)
        nsymop=2
      endif

c *** x-mirror & y-mirror ---> orthorhombic
      if(icrysym.eq.5) then
        h(1,1,3)=-1.d0
        h(2,2,3)= 1.d0
        h(3,3,3)= 1.d0

        h(1,1,4)= 1.d0
        h(2,2,4)=-1.d0
        h(3,3,4)= 1.d0
        nsymop=4
      endif

c *** cubic symmetry
      if(icrysym.eq.1) then

c *** rotations of (pi/3) & (2*pi/3) around <111>
        hx(1,3,1)= 1.d0
        hx(2,1,1)= 1.d0
        hx(3,2,1)= 1.d0

        hx(1,2,2)= 1.d0
        hx(2,3,2)= 1.d0
        hx(3,1,2)= 1.d0

        do m=1,2
          do n=1,nsymop
            mn=m*nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
            enddo
            enddo
            enddo
          enddo
        enddo
        nsymop=mn

c *** mirror across the plane (110)
        hx(1,2,3)= 1.d0
        hx(2,1,3)= 1.d0
        hx(3,3,3)= 1.d0

        do n=1,nsymop
          mn=nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              h(i,j,mn)=h(i,j,mn)+hx(i,k,3)*h(k,j,n)
            enddo
            enddo
            enddo
        enddo
        nsymop=mn

c *** rotations of 90, 180, 270 around x3

        do m=1,3
          ang=pi/2.*float(m)
          hx(1,1,m)= cos(ang)
          hx(2,2,m)= cos(ang)
          hx(3,3,m)= 1.0
          hx(1,2,m)=-sin(ang)
          hx(2,1,m)= sin(ang)
          hx(1,3,m)= 0.0
          hx(3,1,m)= 0.0
          hx(2,3,m)= 0.0
          hx(3,2,m)= 0.0
        enddo

        do m=1,3
          do n=1,nsymop
            mn=m*nsymop+n
              do i=1,3
              do j=1,3
              do k=1,3
                h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
              enddo
              enddo
              enddo
          enddo
        enddo
        nsymop=mn

      endif                    !end of condition for icrysym=1

c *** hexagonal, trigonal and tetragonal symmetry

      if(icrysym.ge.2 .and. icrysym.le.4) then
        if(icrysym.eq.2) nrot=6
        if(icrysym.eq.3) nrot=3
        if(icrysym.eq.4) nrot=4

c *** mirror plane at 30 deg or 60 deg or 45 deg with respect to x1
        ang=pi/float(nrot)
        h(1,1,2)= cos(ang)**2-sin(ang)**2
        h(2,2,2)=-h(1,1,2)
        h(3,3,2)= 1.d0
        h(1,2,2)= 2.*cos(ang)*sin(ang)
        h(2,1,2)= h(1,2,2)
        nsymop=2

c *** rotations of 2*pi/6 around axis <001> for hexagonals.
c *** rotations of 2*pi/3 around axis <001> for trigonals.
c *** rotations of 2*pi/8 around axis <001> for trigonals.
        do nr=1,nrot-1
          ang=nr*2.*pi/nrot
          hx(1,1,nr)= cos(ang)
          hx(2,2,nr)= cos(ang)
          hx(3,3,nr)= 1.d0
          hx(1,2,nr)=-sin(ang)
          hx(2,1,nr)= sin(ang)
        enddo

        do m=1,nrot-1
          do n=1,nsymop
            mn=m*nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
            enddo
            enddo
            enddo
          enddo
        enddo
        nsymop=mn

      endif               !end of condition for icrysym= 2,3,4

c     write(4,*)
c     write(4,'(''  # of symmetry operations='',i4)') nsymop
c     write(4,'(''  symmetry matrices'')')
c     write(4,'(i3,9f7.3)') (n,((h(i,j,n),j=1,3),i=1,3),n=1,nsymop)

      endif               !end of condition for ioption=1

c **************************************************************************
c   Reads Miller-Bravais indices for cubic (1), tetragonal (4), ortho-
c   rhombic (5), monoclinic (6) & triclinic (7) systems in 3-index notation.
c   For hexagonal (2) & trigonal (3) systems reads 4-index notation.
c   Converts indices of plane normal and slip direction into normalized
c   vectors sn(i) and sb(i), respectively.
c **************************************************************************

      if (ioption.eq.2 .or. ioption.eq.3) then

        nind=3
        if (icrysym.eq.2 .or. icrysym.eq.3) nind=4
        if (ioption.eq.2) then
          read(ur1,*) (isn(i),i=1,nind),(isb(i),i=1,nind)
        else if (ioption.eq.3) then
          read(ur1,*) (isn(i),i=1,nind),chi,eta
cas          if(nind.eq.3) write(10,'(3i4,2f10.1)') (isn(i),i=1,3),chi,eta
cas          if(nind.eq.4) write(10,'(4i4,2f10.1)') (isn(i),i=1,4),chi,eta
          eta=eta*pi/180.
          chi=chi*pi/180.
          sb(1)=cos(eta)*sin(chi)
          sb(2)=sin(eta)*sin(chi)
          sb(3)=         cos(chi)
        endif
        if(nind.eq.4) isn(3)=isn(4)

        sn(1)= isn(1)/cvec(1,1)
        sn(2)=(isn(2)-cvec(1,2)*sn(1))/cvec(2,2)
        sn(3)=(isn(3)-cvec(1,3)*sn(1)-cvec(2,3)*sn(2))/cvec(3,3)
        snnor=sqrt(sn(1)**2+sn(2)**2+sn(3)**2)
        do j=1,3
          sn(j)=sn(j)/snnor
          if(abs(sn(j)).lt.1.e-03) sn(j)=0.
        enddo

        if(ioption.eq.2) then

          IF (ICRYSYM.EQ.2 .OR. ICRYSYM.EQ.3) THEN
            ISB(1)=ISB(1)-ISB(3)
            ISB(2)=ISB(2)-ISB(3)
            ISB(3)=ISB(4)
          ENDIF
          do i=1,3
            sb(i)=isb(1)*cvec(i,1)+isb(2)*cvec(i,2)+isb(3)*cvec(i,3)
          enddo
          sbnor=sqrt(sb(1)**2+sb(2)**2+sb(3)**2)
          do j=1,3
            sb(j)=sb(j)/sbnor
            if(abs(sb(j)).lt.1.e-03) sb(j)=0.
          enddo

          prod=sn(1)*sb(1)+sn(2)*sb(2)+sn(3)*sb(3)
          IF(PROD.GE.1.E-3) THEN
            WRITE(*,'('' SYSTEM IS NOT ORTHOGONAL !!'')')
            WRITE(*,'('' ISN='',3I7)') (ISN(J),J=1,3)
            WRITE(*,'('' ISB='',3I7)') (ISB(J),J=1,3)
            WRITE(*,'(''   N='',3F7.3)') (SN(J),J=1,3)
            WRITE(*,'(''   B='',3F7.3)') (SB(J),J=1,3)
            STOP
          ENDIF

        endif      ! end of if(ioption.eq.2)

      endif

c **************************************************************************
c *** generates all symmetry related vectors sneq(i,n) with z>0.
c *** eliminates redundant poles: coincidents and opposites
c **************************************************************************

      if(ioption.eq.3) then

        do n=1,nsymop
          itag(n)=0
          do i=1,3
          sneq(i,n)=0.d0
            do j=1,3
              sneq(i,n)=sneq(i,n)+h(i,j,n)*sn(j)
            enddo
          enddo
        enddo

        do m=1,nsymop-1
          if(itag(m).eq.0) then
            do n=m+1,nsymop
              snpro=sneq(1,m)*sneq(1,n)+sneq(2,m)*sneq(2,n)
     #             +sneq(3,m)*sneq(3,n)
              if(snpro.le. 1.001 .and. snpro.ge. 0.999) itag(n)=1
              if(snpro.ge.-1.001 .and. snpro.le.-0.999) itag(n)=1
            enddo
          endif
        enddo

        npol=0
        do n=1,nsymop
          if(itag(n).eq.0) then
            npol=npol+1
            isign=1
            if(sneq(3,n).lt.0.) isign=-1
            sneq(1,npol)=isign*sneq(1,n)
            sneq(2,npol)=isign*sneq(2,n)
            sneq(3,npol)=isign*sneq(3,n)
          endif
        enddo

      endif            !end of ioption=3
c **************************************************************************

      return
      end
c
c------------------------------------------------------------------------------
c
cas
cas      SUBROUTINE DATA_CRYSTAL (IPH)
cas
      SUBROUTINE DATA_CRYSTAL_AS (IPH, nstrengthsMax)
cas
      INCLUDE 'vpsc_as.dim'

      DIMENSION SN(3),SB(3)
      DIMENSION AUX3(3,24),AUX5(5),AUX33(3,3)
      DIMENSION AUX55(5,5),AUX3333(3,3,3,3)
      PARAMETER (NMFILE=12)
      DIMENSION HSELFX(NMFILE),HLATEX(NMFILE,NMFILE),MODE(NMFILE)

casC
casC     WRITES CRYSTAL DATA FILE INTO 'RUN_LOG.OUT' FILE
cas      WRITE(10,*)
cas      WRITE(10,'('' **** CRYSTAL DATA FILE ****'')')
cas      DO IDUM=1,100
cas        READ(UNIT=UR1,END=10,FMT='(A)') PROSA
cas        WRITE(10,'(A)') PROSA
cas      ENDDO
cas   10 REWIND UR1
cas      WRITE(10,*)
C
C *** Reads crystal symmetry and unit cell parameters.
C *** Generates all symmetry operations associated with CRYSYM.

      !write(*,*) 'ur1 = ', ur1
      call crystal_symmetry (1,ur1,icrysym,sn,aux3,sb,npoles)

C *** READS SINGLE CRYSTAL ELASTIC STIFFNESS
      READ(UR1,'(A)') PROSA
      READ(UR1,*)     ((C2CA(I,J,IPH),J=1,6),I=1,6)
C *** READS SINGLE CRYSTAL THERMAL EXPANSION COEFFICIENTS
      READ(UR1,'(A)') PROSA
      READ(UR1,'(A)') PROSA

C *** ESTIMATES ELASTIC BULK AND COMPRESSIBILITY MODULI (PROVISORY)
C     ALF=(C2CA(1,1,IPH)+C2CA(2,2,IPH)+C2CA(3,3,IPH))/3.
C     BET=(C2CA(1,2,IPH)+C2CA(1,3,IPH)+C2CA(2,3,IPH))/3.
C     BULK=(3.*ALF+6.*BET)/9.

C *** READS INFORMATION ABOUT SLIP AND TWINNING SYSTEMS
      READ(UR1,'(A)') PROSA
      READ(UR1,*)     NMODESX
      READ(UR1,*)     NMODES(IPH)
      READ(UR1,*)     (MODE(I),I=1,NMODES(IPH))

      IF(NMODMX.GT.NMFILE) THEN
        WRITE(*,'('' NMODMX IN SX FILE IS'',I3)') NMODMX
        WRITE(*,'('' CHANGE PARAMETER NMFILE INSIDE DATA_CRYSTAL'')')
        STOP
      ENDIF

      IF(NMODES(IPH).GT.NMODMX) THEN
        WRITE(*,'('' NMODES IN PHASE'',I3,'' IS'',I3)') IPH,NMODES(IPH)
        WRITE(*,'('' CHANGE PARAMETER NMODMX IN VPSC6.DIM'')')
        STOP
      ENDIF
cas
      IF(NMODES(IPH).GT.nstrengthsMax) THEN
        WRITE(*,'('' NMODES IN PHASE'',I3,'' IS'',I3)') IPH,NMODES(IPH)
        WRITE(*,
     #'('' CHANGE PARAM maxStrengthsPerPhase IN FUNCTION init_vpsc'')')
        STOP
      ENDIF
cas
      ICS=1              ! ICS=1 --> CENTRO-SYMMETRIC SCYS
      NTWMOD(iph)=0
      NSYST(iph) =0
      NTWSYS(iph)=0
      KOUNT=1

C *** READS DEFORMATION MODES AND ASSOCIATED PARAMETERS FROM FILECRYS

      DO 100 NM=1,NMODESX

        READ(UR1,'(a)') PROSA
        READ(UR1,*)     MODEX,NSMX,NRSX,ISENSEX
        IF(MODEX.NE.NM) THEN
          WRITE(*,*) ' WARNING !!!'
          WRITE(*,*) ' MODE NUMBERS MUST BE SEQUENTIAL IN CRYSTAL FILE'
          STOP
        ENDIF

C *** SKIPS NON-ACTIVE MODE IF IT IS NOT IN THE LIST.
        IF(MODEX.NE.MODE(KOUNT)) THEN
          READ(UR1,*)
          READ(UR1,*)
          READ(UR1,*)
          DO IDUMMY=1,NSMX
            READ(UR1,*)
          ENDDO
          GO TO 100
        ENDIF

        IF(ISENSEX.EQ.0) ICS=0      ! ICS=0 --> NON-CENTRO-SYMM SCYS

        READ(UR1,*)     TWSHX,ISECTWX,THRES1X,THRES2X
c$$$        READ(UR1,*)     TAU0X,TAU1X,THET0X,THET1X,HPFACX,GNDFACX
        READ(UR1,*)     TAU0X,TAU1X,THET0X,THET1X
        HSELFX(KOUNT)=1.
        READ(UR1,*)     (HLATEX(KOUNT,JM),JM=1,NMODES(IPH))

C *** CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C       TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C       TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C       THETA0=0 FORCES NO-HARDENING.
C *** IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.

        IF(IHARDLAW.NE.1)
     #    CALL DATA_CRYSTAL_VOCE(KOUNT,IPH,TAU0X,TAU1X,THET0X,THET1X)

      NSM(KOUNT,IPH)=NSMX

C *** VERIFICATION OF TWINNING DATA TO BE SURE PROGRAM WILL RUN PROPERLY

      IF(TWSHX.EQ.0. .AND. NTWMOD(IPH).NE.0) THEN
        WRITE(*,'('' TWIN SHEAR='',E10.3,''  MODE'',I4,''  PHASE'',I4)')
     #                        TWSHX,NM,IPH
        WRITE(*,*) ' TWINNING MODES MUST FOLLOW SLIP MODES'
        WRITE(*,*) ' -->   REORDER CRYSTAL FILE'
        STOP
      ENDIF

      IF(TWSHX.NE.0.) THEN
        NTWMOD(IPH)=NTWMOD(IPH)+1
        IF(NTWMOD(IPH).GT.NTWMMX) THEN
         WRITE(*,'('' NTWMOD IN PHASE'',I3,'' IS'',I3)') IPH,NTWMOD(IPH)
         WRITE(*,'('' CHANGE PARAMETER NTWMMX IN VPSC6.DIM'')')
         STOP
        ENDIF
        TWSH(NTWMOD(IPH),IPH)     =TWSHX
        TWTHRES(1,NTWMOD(IPH),IPH)=THRES1X
        TWTHRES(2,NTWMOD(IPH),IPH)=THRES2X
      ENDIF

      DO 205 JS=1,NSM(kount,iph)

        NSYST(iph)=NSYST(iph)+1
        NSYSX=NSYST(iph)

        IF(NSYST(IPH).GT.NSYSMX) THEN
          WRITE(*,'('' NSYST IN PHASE'',I3,'' IS'',I3)') IPH,NSYST(IPH)
          WRITE(*,'('' CHANGE PARAMETER NSYSMX IN vpsc_as.dim'')')
          STOP
        ENDIF

        IF(TWSHX.NE.0.) THEN
         NTWSYS(iph)=NTWSYS(iph)+1
         IF(NTWSYS(IPH).GT.NTWSMX) THEN
         WRITE(*,'('' NTWSYS IN PHASE'',I3,'' IS'',I3)') IPH,NTWSYS(IPH)
         WRITE(*,'('' CHANGE PARAMETER NTWSMX IN VPSC6.DIM'')')
         STOP
         ENDIF
        ENDIF

C   INITIALIZES PARAMETERS ASSOCIATED WITH EACH SYSTEM IN THE MODE.

        NRS(NSYSX,iph)   =NRSX
        ISENSE(NSYSX,iph)=ISENSEX
        ISECTW(NSYSX,iph)=ISECTWX
        TAU(NSYSX,0,iph) =TAU0X
        TAU(NSYSX,1,iph) =TAU1X
        THET(NSYSX,0,iph)=THET0X
        THET(NSYSX,1,iph)=THET1X
cas        HPFAC(NSYSX,iph) =HPFACX
cas        GNDFAC(NSYSX,iph)=GNDFACX          ! added mar/25/2005

C *** CALCULATES CARTESIAN COMPONENTS OF SLIP AND NORMAL VECTORS

        CALL crystal_symmetry(2,UR1,ICRYSYM,SN,AUX3,SB,NPOLES)
        DO J=1,3
          DNCA(J,NSYSX,IPH)=SN(J)
          DBCA(J,NSYSX,IPH)=SB(J)
        ENDDO

C ***  DEFINES SCHMID VECTOR IN CRYSTAL AXES FOR EACH SYSTEM

        DO I=1,3
        DO J=1,3
          AUX33(I,J)=(DNCA(I,NSYSX,IPH)*DBCA(J,NSYSX,IPH)+
     #                DNCA(J,NSYSX,IPH)*DBCA(I,NSYSX,IPH))/2.
        ENDDO
        ENDDO
        CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,2,5)
        DO I=1,5
          SCHCA(I,NSYSX,IPH)=AUX5(I)
        ENDDO

  205 CONTINUE    ! END OF LOOP OVER DEFORMATION MODES

      KOUNT=KOUNT+1

  100 CONTINUE    ! END OF LOOP OVER ALL MODES IN PHASE 'IPH'

C *** CHECKS WHETHER THE SINGLE CRYSTAL YIELD SURFACE IS OPEN

      NSLSYS=NSYST(IPH)-NTWSYS(IPH)
      DO ICOMP=1,5
        ICLOSEPOS=0
        ICLOSENEG=0
        DO NS=1,NSLSYS
          IF(ABS(SCHCA(ICOMP,NS,IPH)).GT.1.E-3) THEN
            ICLOSEPOS=1
            ICLOSENEG=1
          ENDIF
        ENDDO
        IF(NTWSYS(IPH).NE.0) THEN
          DO NS=NSLSYS+1,NSYST(IPH)
            IF(SCHCA(ICOMP,NS,IPH).GT. 1.E-3) ICLOSEPOS=1
            IF(SCHCA(ICOMP,NS,IPH).LT.-1.E-3) ICLOSENEG=1
          ENDDO
        ENDIF
        IF(ICLOSEPOS.NE.1 .OR. ICLOSENEG.NE.1) THEN
          WRITE(*,'('' WARNING ! THE SCYS IS OPEN FOR PHASE'',I5,
     #              '' ALONG DIRECTION'',I5)') IPH,ICOMP
          STOP
        ENDIF
      ENDDO

C *** INITIALIZE SELF & LATENT HARDENING COEFS FOR EACH SYSTEM OF THE PHASE.
C     ABSOLUTE UNITS ARE ACCOUNTED FOR BY MODULATING FACTOR IN HARDENING LAW.

      I=0
      DO IM=1,NMODES(IPH)
      DO IS=1,NSM(IM,IPH)
        I=I+1
        J=0
        DO JM=1,NMODES(IPH)
        DO JS=1,NSM(JM,IPH)
          J=J+1
          HARD(I,J,IPH)=HLATEX(IM,JM)
        ENDDO
        ENDDO
        HARD(I,I,IPH)=HSELFX(IM)
      ENDDO
      ENDDO

cas      I=0
cas      DO IM=1,NMODES(IPH)
cas        WRITE(10,*)
cas        WRITE(10,'(''  NORMAL & BURGERS FOR MODE'',I3,'' IN PHASE'',
cas     #             I3)') IM,IPH
cas        DO IS=1,NSM(IM,IPH)
cas          I=I+1
cas          WRITE(10,'(3F10.3,3X,3F10.3)') (DNCA(J,I,IPH),J=1,3),
cas     #                                   (DBCA(J,I,IPH),J=1,3)
cas        ENDDO
cas      ENDDO

C     WRITE(10,*)
C     WRITE(10,'(''  HARDENING MATRIX FOR PHASE'',I3)') IPH
C     DO IS=1,NSYST(IPH)
C       WRITE(10,'(24F5.1)') (HARD(IS,JS,IPH),JS=1,NSYST(IPH))
C     ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *********************************************************************
C     SUBROUTINE DATA_CRYSTAL_VOCE      --->      VERSION 12/JUN/2001
C *********************************************************************
C
C     CHECKS WHETHER VOCE PARAMETERS ARE KOSHER:
C        TAU0>0 , TAU1 >= 0 , THET0 >= THET1 >= 0
C        TAU1=0   CORRESPONDS TO LINEAR HARDENING.
C        THETA0=0 FORCES NO-HARDENING.
C     IF VOCE PARAMETERS ARE NON-KOSHER CHECKS FOR ILL-POSED HARDENING.
C *********************************************************************

      SUBROUTINE DATA_CRYSTAL_VOCE(KOUNT,IPH,TAU0X,TAU1X,THET0X,THET1X)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

   70 FORMAT(' *** MODE',I3,'   IN PHASE',I3)
   71 FORMAT('     TAU0.LE.0 --> ILL-POSED HARDENING  !!!')
   72 FORMAT('     THETA1<0 --> NON-KOSHER HARDENING MAY GIVE TAU<0')
   73 FORMAT('     TAU1<0  --> NON-KOSHER HARDENING MAY GIVE TAU<0')
   75 FORMAT('     THETA0=0 --> WILL RESET TAU1=THETA1=0')
   76 FORMAT('     |THETA1|.GE.|THETA0| --> NON-KOSHER HARDENING')
   77 FORMAT('     THETA0<0 --> HARDENING DEPENDS ON |THETA0| ANYWAY')

      TINY=1.E-4*TAU0X
      IF(TAU0X.LE.0.0) THEN
        WRITE(*,70) KOUNT,IPH
        WRITE(*,71)
        STOP
      ENDIF
      IF(ABS(THET1X).LE.TINY) THEN
        THET1X=0.
      ELSE IF(THET1X.LT.0.0) THEN
        WRITE(*,70) KOUNT,IPH
        WRITE(*,72)
        STOP
      ENDIF
      IF(TAU1X.LT.0) THEN
        WRITE(*,70) KOUNT,IPH
        WRITE(*,73)
        STOP
      ENDIF
      IF(ABS(THET0X).LE.TINY) THEN
        IF(ABS(TAU1X).LE.TINY) THEN
          TAU1X =0.
          THET0X=THET1X
        ENDIF
        IF(ABS(TAU1X).GT.TINY) THEN
          WRITE(*,70) KOUNT,IPH
          WRITE(*,75)
          STOP
          TAU1X =0.
          THET0X=0.
          THET1X=0.
        ENDIF
      ENDIF
      IF(ABS(TAU1X).LE.TINY) THEN
        TAU1X =0.0
        THET0X=THET1X
      ENDIF
      IF(THET0X.LT.0.0) THEN
        WRITE(*,70) KOUNT,IPH
        WRITE(*,77)
        STOP
        THET0X=ABS(THET0X)
      ENDIF
      IF(TAU1X.NE.0.) THEN
        IF(ABS(THET1X).GE.ABS(THET0X)) THEN
          WRITE(*,70) KOUNT,IPH
          WRITE(*,76)
          STOP
        ENDIF
      ENDIF

      RETURN
      END
c
c------------------------------------------------------------------------------
c
cas
cas      SUBROUTINE DATA_GRAIN (IPH)
cas
      SUBROUTINE DATA_GRAIN_AS (IPH,NGRAIN,ITNGRTOT,EULANGBIG)
cas
      INCLUDE 'vpsc_as.dim'

      DIMENSION AA(3,3),FIJX(3,3),FNEW(3,3),EULANG(3)
      CHARACTER NOMEN*1
      DIMENSION AX(3)
cas
cass      dimension eulangbig(:,:)
      dimension eulangbig(3,ITNGRTOT)
cass
      nomen='B'
cass
      NGR(IPH)=NGR(IPH-1)+NGRAIN

      if(ngr(iph).gt.ngrmx) then
        write(*,'('' number of grains exceeds dimension !!'')')
        write(*,'('' --> increase parameter NGRMX to'',i7)') ngr(iph)
        stop
      endif

C ***************************************************************************
C     READS EULER ANGLES, CONVERTS TO BUNGE NOTATION, CALCULATES ROT MATRIX

      TOTWGT=0.
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
cas
cas        READ(UR2,*) (EULANG(I),I=1,3),WGT(KKK)
cas
         eulang(1)=eulangbig(1,kkk)
         eulang(2)=eulangbig(2,kkk)
         eulang(3)=eulangbig(3,kkk)
cas
cas        TOTWGT=TOTWGT+WGT(KKK)
cas
        if(nomen.eq.'B' .or. nomen.eq.'b') then
          eul1=  eulang(1)
          eul2=  eulang(2)
          eul3=  eulang(3)
        else if(nomen.eq.'K' .or. nomen.eq.'k') then
          eul1= (eulang(1)-90.)
          eul2= -eulang(2)
          eul3=(-eulang(3)-90.)
        else if(nomen.eq.'R' .or. nomen.eq.'r') then
          eul1= (eulang(1)+90.)
          eul2=  eulang(2)           ! fixed 29/apr/02
          eul3= (eulang(3)-90.)
        else
          write(*,'(/,'' CANNOT IDENTIFY EULER ANGLE CONVENTION  !!'')')
          stop
        endif

C     CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
C     SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.

        CALL EULER(2,EUL1,EUL2,EUL3,AA)
        DO J=1,3
        DO K=1,3
          AG(J,K,KKK)=AA(K,J)
        ENDDO
        ENDDO

      ENDDO     ! END OF LOOP OVER GRAINS IN EACH PHASE

casC     DOUBLE RENORMALIZATION OF WEIGHTS: FIRST NORMALIZE THE WEIGHTS WITHIN
casC     EACH PHASE, NEXT RENORMALIZE TO THE VOLUME FRACTION OF THE PHASE.
cas
cas      DO KKK=NGR(IPH-1)+1,NGR(IPH)
cas        WGT(KKK)=WGT(KKK)/TOTWGT*WPH(IPH)
cas      ENDDO

C ***************************************************************************
C     INITIAL F TENSOR AND EIGENVECTORS:

C     * IF ISHAPE=0 ASSUMES SAME INITIAL SHAPE (axisph) AND ORIENTATION
C       (eulerph) FOR ALL THE GRAINS IN THE PHASE.
C       CALCULATES ESHELBY TENSOR WITH THE AVERAGE GRAIN SHAPE AND DOES
C       NOT KEEP TRACK OF LOCAL GRAIN SHAPE.
C
C       'eulerph' ANGLES OF (G) WRT (S).
C       'aa'    TRANSFORMS FROM (S) TO (G)
C       'fijx'  COLUMNS ARE GRAIN AXES EXPRESSED IN GRAIN SYSTEM
C       'fnew'  COLUMNS ARE GRAIN AXES EXPRESSED IN SAMPLE SYSTEM

      da=eulerph(1,iph)
      db=eulerph(2,iph)
      dc=eulerph(3,iph)

      call euler(2,da,db,dc,aa)

      do i=1,3
      do j=1,3
        fijx(i,j)=(i/j)*(j/i)*AXISPH(0,I,IPH)
      enddo
      enddo

      do j=1,3
      do i=1,3
        fnew(i,j)=0.
        do m=1,3
          fnew(i,j)=fnew(i,j)+aa(m,i)*fijx(m,j)
        enddo
      enddo
      enddo

cas      do i=1,3
cas      do j=1,3
cas        fijph(i,j,iph)=fnew(i,j)
cas      enddo
cas      enddo

      do i=1,3
      do j=1,3
       aux=0.
       do k=1,3
        aux=aux+fnew(i,k)*fnew(j,k)
       enddo
       f2ijph(i,j,iph)=aux
      enddo
      enddo

      RETURN
      END
c
c------------------------------------------------------------------------------
c
      function det(a)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      dimension a(3,3)

      det=a(1,1)*a(2,2)*a(3,3)
     #   +a(1,2)*a(2,3)*a(3,1)
     #   +a(1,3)*a(2,1)*a(3,2)
     #   -a(1,3)*a(2,2)*a(3,1)
     #   -a(2,3)*a(3,2)*a(1,1)
     #   -a(1,2)*a(2,1)*a(3,3)
      return
      end
c
c------------------------------------------------------------------------------
c
C
C ***********************************************************************
C     SUBROUTINE ESHELBY      --->      VERSION 07/OCT/05
C
C     IOPTION=0: Initialize arrays assoc. with Gauss integration points.
C     IOPTION=1: Calculate elastic Eshelby tensor for elastic inclusion.
C     IOPTION=2: Calculate incompressible Eshelby tensors ESIM (strain-
C                rate) & ESCR (spin-rate) associated with the visco-
C                plastic inclusion.
C     IOPTION=3: Calculate incompressible and hydrostatic Eshelby tensors
C                PESH (deviatoric pressure), PDIL (spherical pressure) &
C                DESH (dilatation) for visco-plastic inclusion.
C     IOPTION=4: Calculates d(S)/d(M) (term 1)
C     IOPTION=5: Calculates d(S)/d(M) (term 2)
C
C     Options 2-3-4-5 admit a non-zero visco-plastic bulk modulus KEFF.
C
C     Algorithms are based in Lebensohn et al, MSMSE 6 (1998) p.447.
C     Uses explicit matrix inversion and explicit Voigt notation (when
C     possible) to optimize computing time.
C
C     Modified oct/2005 to adapt number of integration points to the shape
C     of the ellipsoid in order to keep Eshelby tensor within a certain
C     tolerance (based on analysis done by Gwenaelle Proust).
C     Aspect ratio criterion was adopted for the case when AXIS(2) is
C     largest and AXIS(3) is smallest.
C ***********************************************************************

      SUBROUTINE ESHELBY (axis,c4,keff,esim,escr,
     #                    desh,pesh,pdil,dldm,dsddm,ioption)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION c4(3,3,3,3),esim(3,3,3,3),escr(3,3,3,3)
      DIMENSION p(3,3),pesh(3,3),desh(3,3)
      DIMENSION c2(6,6),gamma2(6,6),gamma4(3,3,3,3)
      DIMENSION axis(3),x1(10),a1(10),a1inv(10)
      DIMENSION aa1x(6),aa2x(3,3),aaww1x(6),aaww2x(3,3)

      dimension dldm(3,3,3,3),dsddm(3,3,3,3),dldm2(6,6)
CFEB
      dimension da1(6),da(4,4),dainv(4,4),ainv(4,4)
      dimension aux33(3,3),aux66(6,6),aux3333(3,3,3,3),aux44(4,4)
CFEE
      PARAMETER (ngaumx=32,ngaumx2=1024)
      DIMENSION xph(ngaumx),xth(ngaumx),wph(ngaumx),wth(ngaumx)
      COMMON/ESHELBY1/ngaussph(3),ngaussth(3)
      COMMON/ESHELBY2/alpha(3,3,ngaumx2),aa1(3,6,ngaumx2),
     #                aww(3,3,ngaumx2),aaww1(3,6,ngaumx2),ww(3,ngaumx2)

      REAL*8    keff
      INTEGER   case

      pi=4.*atan(1.0)

C ***********************************************************************
C     INITIALIZATION RUN
C     Calculates Gauss-Legendre integration points and weights in the
c     interval [0,pi].
C     Initializes arrays associated with each point to avoid repeating
C     its calculation at every call.
C ***********************************************************************

      if(ioption.eq.0) then

        ngaussph(1)=16
        ngaussth(1)=16
        ngaussph(2)=32
        ngaussth(2)=32
        ngaussph(3)=16
        ngaussth(3)=32

        do case=1,3
          if(ngaussph(case).gt.ngaumx.or.ngaussth(case).gt.ngaumx) then
            write(*,*) ' DIMENSION ngaumx EXCEEDED IN SUBR ESHELBY !!'
            stop
          endif

          call gauleg(0.0,pi,xph,wph,ngaussph(case))
          call gauleg(0.0,pi,xth,wth,ngaussth(case))

c *** integration [0,pi][0,pi] adds a factor 2 in Eqs. B11 & B14.

          do ith=1,ngaussth(case)
            sinth=sin(xth(ith))
            costh=cos(xth(ith))
            simbtet=wth(ith)*sinth/(2.0*pi)

            do iph=1,ngaussph(case)
              ny=iph+(ith-1)*ngaussph(case)
              ww(case,ny)=simbtet*wph(iph)
              alpha(case,1,ny)=sinth*cos(xph(iph))
              alpha(case,2,ny)=sinth*sin(xph(iph))
              alpha(case,3,ny)=costh

              do i=1,3
              do j=1,3
                aa2x(i,j)  =alpha(case,i,ny)*alpha(case,j,ny)
                aaww2x(i,j)=aa2x(i,j)*ww(case,ny)
              enddo
              enddo
              call voigt(aa1x  ,aa2x  ,c2,c4,2)
              call voigt(aaww1x,aaww2x,c2,c4,2)
              do i=1,6
                aa1(case,i,ny)  =aa1x(i)
                aaww1(case,i,ny)=aaww1x(i)
              enddo

c *** Array AWW is used only if ICAUCHY=1.
              do i=1,3
                aww(case,i,ny)=alpha(case,i,ny)*ww(case,ny)
              enddo
            enddo
          enddo
        enddo      ! end of do case=1,3

        return

      endif      ! ENDIF FOR IOPTION=0

C ***********************************************************************
C     CALCULATION OF ESHELBY TENSORS
C     For given stiffness 'C4' and ellipsoid axes 'AXIS'
C ***********************************************************************

      if(ioption.ge.1) then

        abc=axis(1)*axis(2)*axis(3)
        ratio1=axis(2)/axis(3)
        ratio2=axis(1)/axis(3)
        if(ratio1.lt.25) case=1
        if(ratio1.ge.25 .and. ratio2.le.5) case=2
        if(ratio1.ge.25 .and. ratio2.gt.5) case=3
        npoints=ngaussph(case)*ngaussth(case)

        pdil=0.
        do j=1,3
        do i=1,3
          p(i,j)=0.
        enddo
        enddo
        do j=1,6
        do i=1,6
          gamma2(i,j)=0.
        enddo
        enddo

        call voigt(aa1x,aa2x,c2,c4,4)
        IF(IOPTION.EQ.5) CALL VOIGT(AA1X,AA2X,DLDM2,DLDM,4)

        do ny=1,npoints

c   Compute Voigt components A1(1)-A(6) of tensor A(3,3) defined by Eq.B3:
c   --->  A(i,j)=L(i,j,k,l)*a(j)*a(l)

          do i=1,6
            aa1x(i)=aa1(case,i,ny)
          enddo
          call esh_mult_voigt(c2,aa1x,a1)
cw
      IF(IOPTION.EQ.1) THEN

c   If solving an elastic inclusion invert the system
c   --> A(3,3) x X(3,3) = C(3,3)
c   Inverts A(3,3) using explicit Voigt notation.
c   Uses explicit form of C(3,3) to calculate solution in Voigt notation.

            call esh_inv3_voigt(a1,a1inv)
            do i=1,6
              x1(i)=a1inv(i)
            enddo

      ELSE IF(IOPTION.GE.2) THEN

c   If solving a visco-plastic inclusion defines components A1(7) to A1(10).
c   Solves the system given by Eq.B4 --> A(4,4) x X(4,4) = C(4,4)
c   Inverts A(4,4) using explicit Voigt notation.
c   Uses explicit form of C(4,4) to calculate solution in Voigt notation.
c   The solution is symmetric. Numerical deviation from symmetry is averaged.

            a1(7) = alpha(case,1,ny)
            a1(8) = alpha(case,2,ny)
            a1(9) = alpha(case,3,ny)
            a1(10)= 0.
            if(keff.gt.0.) a1(10)=-1./keff

            call esh_inv4_voigt(a1,a1inv)

            do i=1,10
              x1(i)=a1inv(i)
            enddo

      ENDIF
CFEB
      IF(IOPTION.EQ.5) THEN
c
        CALL VOIGT10(A1INV,AINV,1)
c
        call esh_mult_voigt(dldm2,aa1x,da1)
        call voigt(da1,aux33,aux66,aux3333,1)
c
        do i=1,3
        do j=1,3
          da(i,j)=aux33(i,j)
        enddo
        enddo
c
        do i=1,3
        da(i,4)=0.
        da(4,i)=0.
        enddo
        da(4,4)=0.
c
c           dAinv/dp = - Ainv  : dA/dp : Ainv
c                x1  = - a1inv :  da1  : a1inv
c
        do i=1,4
        do j=1,4
        dummy=0.
        do k=1,4
        dummy=dummy+da(i,k)*ainv(k,j)
        enddo
        aux44(i,j)=dummy
        enddo
        enddo

        do i=1,4
        do j=1,4
        dummy=0.
        do k=1,4
          dummy=dummy+ainv(i,k)*aux44(k,j)
        enddo
        dainv(i,j)=-dummy
        enddo
        enddo

        CALL VOIGT10(X1,DAINV,2)

      ENDIF
CFEE
          ro3=((alpha(case,1,ny)*axis(1))**2+
     #         (alpha(case,2,ny)*axis(2))**2+
     #         (alpha(case,3,ny)*axis(3))**2)**1.5
          abcoro3=abc/ro3

c   Compute the Eshelby integral Eq.B11 defining:
c         Gamma(m,j,n,i)=T(m,n,i,j)=a(m)*a(j)*G(n,i)
c   with the property:
c         Gamma(m,j,n,i)=Gamma(j,m,n,i)=Gamma(m,j,i,n)

          do i=1,6
          do j=1,6
            gamma2(i,j)=gamma2(i,j)+aaww1(case,i,ny)*x1(j)*abcoro3
          enddo
          enddo

c   Compute the pressure related Eshelby integral Eq.B14
          if(ioption.eq.3) then
            do j=1,3
            do i=1,3
              p(i,j)=p(i,j)+aww(case,j,ny)*x1(i+6)*abcoro3
            enddo
            enddo
            pdil=pdil+ww(case,ny)*x1(10)*abcoro3
          endif

        end do   ! end of loop over double integration

c ********************************************************************
c   Go back to the 3*3*3*3 notation
        call voigt(aa1x,aa2x,gamma2,gamma4,3)

c   Compute symmetric (distortion) Eshelby tensor from Eq.B9.
c       esim(n,m,k,l)=0.5*(gamma(m,j,n,i)+gamma(n,j,m,i))*c4(i,j,k,l)
c   Compute anti-symmetric (rotation) Eshelby tensor from Eq.B9.
c       escr(n,m,k,l)=0.5*(gamma(m,j,n,i)-gamma(n,j,m,i))*c4(i,j,k,l)

        do l=1,3
        do k=1,3
        do m=1,3
        do n=1,3
c
          dumsim=0.
          dumscr=0.
c
        do j=1,3
        do i=1,3
c
        IF(IOPTION.NE.4) THEN
          dumsim=dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))*c4(i,j,k,l)
          dumscr=dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))*c4(i,j,k,l)
        ELSE
          dumsim=dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))*dldm(i,j,k,l)
          dumscr=dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))*dldm(i,j,k,l)
        ENDIF
c
        enddo
        enddo
c
          IF(IOPTION.LT.4) THEN
            esim(n,m,k,l)=0.5*dumsim
            escr(n,m,k,l)=0.5*dumscr
          ELSE
            dsddm(n,m,k,l)=0.5*dumsim
          ENDIF
c
        enddo
        enddo
        enddo
        enddo

c   Compute pressure & dilatation related Eshelby tensors (Eq.B13)

        if(ioption.eq.3) then
          do l=1,3
          do k=1,3
            pesh(k,l)=0.
            do j=1,3
            do i=1,3
              pesh(k,l)=pesh(k,l)+p(i,j)*c4(i,j,k,l)
            end do
            end do
          end do
          end do
          do j=1,3
          do i=1,3
            desh(i,j)=(p(i,j)+p(j,i))/2.
          end do
          end do
        endif

      endif      !  endif for IOPTION.GE.1

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *************************************************************************
C     SUBROUTINE ESH_INV3_VOIGT   --->   version 23/jul/01
C
C     Inverts the 3x3 symmetric matrix 'A' using explicit Voigt notation:
C     11->1, 22->2, 33->3, 23=32->4, 31=13->5, 12=21->6
C *************************************************************************

      SUBROUTINE ESH_INV3_VOIGT (A,AINV)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION A(10),AINV(10)

      DET = A(1)*A(2)*A(3) + 2*A(4)*A(5)*A(6) - A(1)*A(4)*A(4)
     #     - A(2)*A(5)*A(5) - A(3)*A(6)*A(6)

      AINV(1) = ( A(2)*A(3) - A(4)*A(4))/DET
      AINV(2) = ( A(1)*A(3) - A(5)*A(5))/DET
      AINV(3) = ( A(1)*A(2) - A(6)*A(6))/DET
      AINV(4) = (-A(1)*A(4) + A(5)*A(6))/DET
      AINV(5) = ( A(4)*A(6) - A(2)*A(5))/DET
      AINV(6) = (-A(3)*A(6) + A(4)*A(5))/DET

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C **********************************************************************
C     SUBROUTINE ESH_INV4_VOIGT   --->   VERSION 20/JUL/01

C     Inverts the 4*4 symmetric matrix 'A' using explicit Voigt notation:
C     11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6
C     14-->7, 24-->8, 34-->9, 44-->10.
C **********************************************************************

      SUBROUTINE ESH_INV4_VOIGT (A,AINV)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION A(10),AINV(10)


      ainv(1) = a(2)*a(3)*a(10)+2*a(4)*a(8)*a(9) -
     #          a(2)*a(9)*a(9)-a(3)*a(8)*a(8)-a(4)*a(4)*a(10)

      ainv(2) = a(1)*a(3)*a(10)+2*a(5)*a(7)*a(9) -
     #          a(1)*a(9)*a(9)-a(3)*a(7)*a(7)-a(5)*a(5)*a(10)

      ainv(3) = a(1)*a(2)*a(10)+2*a(6)*a(7)*a(8) -
     #          a(1)*a(8)*a(8)-a(2)*a(7)*a(7)-a(6)*a(6)*a(10)

      ainv(4) = a(1)*a(4)*a(10)+a(5)*a(7)*a(8)+a(6)*a(7)*a(9) -
     #          a(1)*a(8)*a(9)-a(4)*a(7)*a(7)-a(5)*a(6)*a(10)
      ainv(4) =-ainv(4)

      ainv(5) = a(4)*a(6)*a(10)+a(2)*a(7)*a(9)+a(5)*a(8)*a(8) -
     #          a(4)*a(7)*a(8)-a(6)*a(8)*a(9)-a(2)*a(5)*a(10)

      ainv(6) = a(3)*a(6)*a(10)+a(5)*a(8)*a(9)+a(4)*a(7)*a(9) -
     #          a(3)*a(7)*a(8)-a(6)*a(9)*a(9)-a(4)*a(5)*a(10)
      ainv(6) =-ainv(6)

      ainv(7) = a(4)*a(6)*a(9)+a(4)*a(5)*a(8)+a(2)*a(3)*a(7) -
     #          a(4)*a(4)*a(7)-a(2)*a(5)*a(9)-a(3)*a(6)*a(8)
      ainv(7) =-ainv(7)

      ainv(8) = a(1)*a(4)*a(9)+a(5)*a(5)*a(8)+a(3)*a(6)*a(7) -
     #          a(4)*a(5)*a(7)-a(5)*a(6)*a(9)-a(1)*a(3)*a(8)

      ainv(9) = a(1)*a(2)*a(9)+a(5)*a(6)*a(8)+a(4)*a(6)*a(7) -
     #          a(2)*a(5)*a(7)-a(6)*a(6)*a(9)-a(1)*a(4)*a(8)
      ainv(9) =-ainv(9)

      ainv(10)=a(1)*a(2)*a(3)+2*a(4)*a(5)*a(6) -
     #         a(1)*a(4)*a(4)-a(2)*a(5)*a(5)-a(3)*a(6)*a(6)

      det=   a(1)*ainv(1)+   a(2)*ainv(2)+   a(3)*ainv(3)+
     #    2.*a(4)*ainv(4)+2.*a(5)*ainv(5)+2.*a(6)*ainv(6)+
     #    2.*a(7)*ainv(7)+2.*a(8)*ainv(8)+2.*a(9)*ainv(9)+
     #       a(10)*ainv(10)
      det=   det/4.

      do i=1,10
        ainv(i)=ainv(i)/det
      enddo

      return
      end
c
c------------------------------------------------------------------------------
c
C ***********************************************************
      SUBROUTINE ESH_MULT_VOIGT(B,C,A)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

C     Performs the multiplication:
C        A(i,k)=B(i,j,k,l)*C(j,l) using Voigt's notation
C        B is a 6*6 symmetric matrix
C        C is a 3*3 symmetric tensor
C        A will be a 3*3 symmetric tensor

      DIMENSION B(6,6),C(6),A(6)

      A(1)=B(1,1)*C(1)+B(6,6)*C(2)+B(5,5)*C(3)
     #    +2*(B(5,6)*C(4)+B(1,5)*C(5)+B(1,6)*C(6))

      A(2)=B(6,6)*C(1)+B(2,2)*C(2)+B(4,4)*C(3)
     #    +2*(B(2,4)*C(4)+B(4,6)*C(5)+B(2,6)*C(6))

      A(3)=B(5,5)*C(1)+B(4,4)*C(2)+B(3,3)*C(3)
     #    +2*(B(3,4)*C(4)+B(3,5)*C(5)+B(4,5)*C(6))

      A(4)=B(5,6)*C(1)+B(2,4)*C(2)+B(3,4)*C(3)
     #      +(B(2,3)+B(4,4))*C(4)
     #      +(B(3,6)+B(4,5))*C(5)
     #      +(B(4,6)+B(2,5))*C(6)

      A(5)=B(1,5)*C(1)+B(4,6)*C(2)+B(3,5)*C(3)
     #      +(B(3,6)+B(4,5))*C(4)
     #      +(B(1,3)+B(5,5))*C(5)
     #      +(B(1,4)+B(5,6))*C(6)

      A(6)=B(1,6)*C(1)+B(2,6)*C(2)+B(4,5)*C(3)
     #      +(B(4,6)+B(2,5))*C(4)
     #      +(B(1,4)+B(5,6))*C(5)
     #      +(B(1,2)+B(6,6))*C(6)

      RETURN
      END
c
c------------------------------------------------------------------------------
c
c *****************************************************************************
      subroutine euler (iopt,ph,th,tm,a)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES (in degrees) OF ca REFERRED TO sa.
c *****************************************************************************

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      dimension a(3,3)
      pi=4.*atan(1.d0)

      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=sin(ph*pi/180.)
        cph=cos(ph*pi/180.)
        sth=sin(th*pi/180.)
        cth=cos(th*pi/180.)
        stm=sin(tm*pi/180.)
        ctm=cos(tm*pi/180.)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif

      return
      end
c
c------------------------------------------------------------------------------
c
cass      subroutine get_twinor_as (iph, twinor)
      subroutine get_twinor_as (iph, twinor, ntws)
c
      INCLUDE 'vpsc_as.dim'
c
cass      dimension twinor(:,:,:)
cass
      dimension twinor(3,3,ntws)
cass
      dimension atwin(3,3),btwin(3)
c
        ITS=NSYST(IPH)-NTWSYS(IPH)
        KTS=0
        DO ITM=1,NTWMOD(IPH)
          KK=ITM+NMODES(IPH)-NTWMOD(IPH)
          DO NSMX=1,NSM(KK,IPH)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
c
            btwin(1)=dbca(1,ITS,iph)
            btwin(2)=dbca(2,ITS,iph)
            btwin(3)=dbca(3,ITS,iph)
c     
            call twin_orientation (btwin,atwin)
            
C     NOTE ABOUT SUBROUTINE TWIN_ORIENTATION (IN VPSC.SUB):
C     GIVEN THE BURGERS VECTOR 'BTWIN' OF THE TWIN SYSTEM (IN CRYSTAL 
C     AXES) MAKES A ROTATION OF 180 DEG AROUND IT TO DEFINE THE 
C     ORIENTATION OF THE TWIN-RELATED CRYSTAL. THE MATRIX 'ATWIN' 
C     TRANSFORMS FROM TWIN TO CRYSTAL AXES
            
            do I=1,3
               do J=1,3
                  twinor(I,J,KTS)=atwin(I,J)
               enddo
            enddo

          ENDDO
        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE
c
      return
      end
c
c------------------------------------------------------------------------------
c
C ****************************************************************************
C     SUBROUTINE GRAIN_RATE_AND_MODULI (ex MICRO) -->  VERSION 17/DEC/02
C
C     GIVEN THE STRESS 'X' IN GRAIN 'KKK', CALCULATES STRAIN-RATE AND
C     VISCO-PLASTIC MODULI USING THE RATE SENSITIVITY KINEMATIC LAW.
C ****************************************************************************

      SUBROUTINE GRAIN_RATE_AND_MODULI (JSC,KCALMOD,KGX,KKK,IPHEL,IPH)

      INCLUDE 'vpsc_as.dim'

      DIMENSION SX(5),DX(5),RSS(NSYSMX)
      REAL*8 RSSX

      COMMON/NLNR /DB(5),XMASTX(5,5),SCX(5,NSYSMX),TAUX(NSYSMX),
     #             GAMD0X(NSYSMX),NSYSTX,NRSX(NSYSMX),ISENSEX(NSYSMX)

C     COPY MAIN ARRAYS INTO AUXILIAR ARRAYS FOR COMPUTATIONAL EFFICIENCY

      NSYSTX=NSYST(IPHEL)
      DO IS=1,NSYSTX
        ISENSEX(IS)=ISENSE(IS,IPHEL)
        NRSX(IS)   =NRS(IS,IPHEL)
        if(jsc.eq.1.and.irsvar.eq.1) NRSX(IS)=JXRS
cw
        TAUX(IS)  =CRSS(IS,KKK)
        GAMD0X(IS)=GAMD0G(KGX)
        DO J=1,5
          SCX(J,IS)=SCH(J,IS,KGX)
        ENDDO
      ENDDO

c     write(*,*)
c     write(*,*) 'SG inside grain_rate for KKK=',kkk
c     write(*,'(5e12.3)') (sg(i,kkk),i=1,5)
c     pause

      DO I=1,5
        SX(I)=SG(I,KKK)
      ENDDO

C     GETS RESOLVED SHEAR STRESSES 'rssx' AND SHEAR RATES 'gamdot'.
C       SIGN(GAMDOT)=SIGN(RSSX).
C       NRS CAN BE EVEN OR ODD.
C       RSS IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

      DO IS=1,NSYSTX
        RSSX=SX(1)*SCX(1,IS)+SX(2)*SCX(2,IS)+SX(3)*SCX(3,IS)+
     #       SX(4)*SCX(4,IS)+SX(5)*SCX(5,IS)
        IF(.NOT.(RSSX.GT.0 .OR. ISENSEX(IS).EQ.1)) RSSX=0.
        RSSX=RSSX/TAUX(IS)

c     write(*,*)
c     write(*,'('' RSSX and TAUX for sys'',i3,2e12.3)') is,rssx,taux(is)

        RSS(IS)       =GAMD0X(IS)*ABS(RSSX**(NRSX(IS)-1))/TAUX(IS)
        GAMDOT(IS,KGX)=GAMD0X(IS)*ABS(RSSX**NRSX(IS))*SIGN(1.0,RSSX)
      ENDDO

C     CALCULATE STRAIN-RATE IN GRAIN 'dg'

      DO I=1,5
        DG(I,KGX)=0.
        DO IS=1,NSYSTX
          DG(I,KGX)=DG(I,KGX)+SCX(I,IS)*GAMDOT(IS,KGX)
        ENDDO
        DX(I)=DG(I,KGX)
      ENDDO

c     write(*,*)
c     write(*,*) 'DG inside grain_rate for KGX=',kgx
c     write(*,'(5e12.3)') (dg(i,kgx),i=1,5)
c     pause

      IF(KCALMOD.EQ.1) THEN

C     CALCULATE CRYSTAL COMPLIANCE --> explain next lines!

      DO I=1,5
      DO J=1,5
        XMCTG (I,J,KGX)=0.
        DO IS=1,NSYSTX
          if(interaction.eq.2.or.interaction.eq.3.or.
     #      interaction.eq.4) then
            XMCTG(I,J,KGX)=XMCTG(I,J,KGX)+SCX(I,IS)*SCX(J,IS)*RSS(IS)
          else
            XMCTG(I,J,KGX)=XMCTG(I,J,KGX)+NRSX(IS)*
     #                    SCX(I,IS)*SCX(J,IS)*RSS(IS)
          endif
        ENDDO
      ENDDO
      ENDDO
CFEB
      if(interaction.eq.5) then
         write(*,*) "Not implemented"
c$$$        DO IS=1,NSYSTX
c$$$          aso(is,kgx)=nrsx(is)*rss(is)
c$$$          eso(is,kgx)=(1-nrsx(is))*gamdot(is,kgx)
c$$$        ENDDO
      endif
CFEE
      IF(INTERACTION.EQ.1) THEN
        do i=1,5
          dczero(i,kgx)=dx(i)
          do j=1,5
            dczero(i,kgx)=dczero(i,kgx)-xmctg(i,j,kgx)*sx(j)
          enddo
        enddo
cc      WRITE (*,'(i5,5F10.3)') KGX,(DCZERO(I,KGX),I=1,5)
cc      pause
      ELSE
        do i=1,5
          dczero(i,kgx)=0.
        enddo
      ENDIF        ! INTERACTION ENDIF

      ENDIF        ! KCALMOD ENDIF

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *****************************************************************************
C     SUBROUTINE GRAIN_STRESS (ex VISC)   --->   version of 19/MAR/2003
C
C     GIVEN A GUESS STRESS 'X(I)' AND GRAIN INDEX 'KKK', SOLVES INTERACTION
C     EQUATION TO FIND THE STRESS 'X' COMPATIBLE WITH MACROSCOPIC COMPLIANCE.
C     IF INTERX=-1 SOLVES ONLY THE POWER LAW FOR A RELAXED CONSTRAINTS CASE.
C     IF INTERX= 0 SOLVES ONLY THE POWER LAW FOR A TAYLOR CASE.
C     IF INTERX> 0 SOLVES INTERACTION EQUATION FOR A SELF-CONSISTENT CASE.
C *****************************************************************************

      SUBROUTINE GRAIN_STRESS (INTERX,KGX,KKK,IPHEL,IPH)

      INCLUDE 'vpsc_as.dim'

      DIMENSION X(5),XORI(5)
cw      COMMON/ASTER/SASTAV(5),SASTBAR(5),DAST(5)
      COMMON/NLNR /DB(5),XMASTX(5,5),SCX(5,NSYSMX),TAUX(NSYSMX),
     #             GAMD0X(NSYSMX),NSYSTX,NRSX(NSYSMX),ISENSEX(NSYSMX)


C     EMPIRIC ALGORITHM TO GET TAULIM FOR NR SUBROUTINE (RL: 01/FEB/00)

      taulim=2.*(tnorm(dbar,5,1)/gamd0g(1))**(1./nrsmin)
      if(taulim.lt.2.) taulim=2.

      iprint=0   ! controls diagnostic print-out
cas      ngg=699
cas      iprx=iprint*(ngg/kkk)*(kkk/ngg)+
cas     #     iprint*((ngg+ngr(1))/kkk)*(kkk/(ngg+ngr(1)))
cas      if(iprx.eq.1) then
cas        write(10,'('' INSIDE GRAIN_STRESS: GRAIN'',i5 )') KKK
cas      endif

C     COPY MAIN ARRAYS INTO AUXILIAR ARRAYS FOR COMPUTATIONAL EFFICIENCY
C     AND TO MAKE 'NEWTON_RAPHSON' A 'STAND-ALONE' SUBROUTINE

      NSYSTX=NSYST(IPHEL)
      DO IS=1,NSYSTX
        ISENSEX(IS)=ISENSE(IS,IPHEL)
        NRSX(IS)   =NRS(IS,IPHEL)
        if(interx.eq.1.and.irsvar.eq.1) NRSX(IS)=JXRS
cw
        TAUX(IS)  =CRSS(IS,KKK)
        GAMD0X(IS)=GAMD0G(KGX)
        DO J=1,5
          SCX(J,IS)=SCH(J,IS,KGX)
        ENDDO
      ENDDO

      DO I=1,5
        X(I)=STRY(I,KGX)
      ENDDO

      IRC=0

C *** CORRECTS STRESS 'X' IF IT EXCEEDS THE YIELD SURFACE.

      TAUMAX=0.
      DO IS=1,NSYSTX
        RSSX=X(1)*SCX(1,IS)+X(2)*SCX(2,IS)+X(3)*SCX(3,IS)+
     #       X(4)*SCX(4,IS)+X(5)*SCX(5,IS)
        IF(.NOT.(RSSX.GT.0 .OR. ISENSEX(IS).EQ.1)) RSSX=0.
        RSSX=RSSX/TAUX(IS)
        IF(ABS(RSSX).GT.TAUMAX) TAUMAX=ABS(RSSX)
      ENDDO

cas      if(iprx.eq.1) then
cas        write(10,'('' GRAIN'',i5,''  TAUMAX='',F10.3)') KKK,TAUMAX
cas      endif

      IF(TAUMAX. LT. 1.E-12) THEN
        WRITE(*,'('' TAUMAX<1e-12 inside subroutine GRAIN_STRESS'')')
        WRITE(*,'('' GRAIN #'',i5)') KKK
        WRITE(*,'('' STRY  ='',5e12.4)') X
        WRITE(*,'('' CRSS  ='',6E12.4)') (TAUX(IS),IS=1,NSYST(IPHEL))
        STOP
      ENDIF

      IF(TAUMAX.GT.TAULIM .OR. INTERX.LE.0) THEN
        DO I=1,5
          X(I)=X(I)/TAUMAX
        ENDDO
      ENDIF

      IF(INTERX.LE.0) THEN
        DO I=1,5
          DB(I)=DBAR(I)
          DO J=1,5
            XMASTX(I,J)=0.
          ENDDO
        ENDDO
      ELSE IF(INTERX.GT.0) THEN
        IF(ISHAPE(IPH).LE.1) THEN
          DO I=1,5
          DO J=1,5
            XMASTX(I,J)=XMASTPH(I,J,IPH)
          ENDDO
          ENDDO
CFEB
cas        ELSE IF(ISHAPE(IPH).GT.1) THEN
cas          DO I=1,5
cas          DO J=1,5
cas            XMASTX(I,J)=XMASTGR(I,J,KKK)
cas          ENDDO
cas          ENDDO
CFEE
        ENDIF

        DO I=1,5
cw          DB(I)=DAST(I)
          DB(I)=DAST(I)+DZERO(I)
          DO J=1,5
            DB(I)=DB(I)+XMASTX(I,J)*SASTAV(J)
          ENDDO
        ENDDO

      ENDIF

C *** CALLS Newton-Raphson SUBROUTINE TO CALCULATE GRAIN STRESS
C *** Internally it does a 5D (IRC=0) or a 3D (IRC=1) convergence.

      DO I=1,5
        XORI(I)=X(I)
      ENDDO
      ITMX=1000
      EPS=5.e-04

        iprx=0
c     if(kgx.eq.1) then
c       write(10,*) 'stress in grain #1 before entering N-R'
c       write(10,'(5e12.3)') x
c       iprx=1
c     endif

cpp      write(*,*) 'KGX=',kgx
cpp      write(*,*) 'DZERO=',dzero
cpp      write(*,*) 'DB=',db
cpp      write(*,*) 'DAST=',dast
cpp      write(*,*) 'SASTAV=',sastav
cpp      write(*,*) 'XMASTX=',xmastx
cpp      pause
c
      CALL NEWTON_RAPHSON (X,IRC,ITMX,EPS,TAULIM,IERROR ,iprx)

      IF(IERROR.GT.0) THEN
        IF(IERROR.EQ.1) WRITE(*,'('' SINGULAR SYSTEM IN NEWTRAPH -->'',
     #       '' CANNOT SOLVE GRAIN'',I6,'' IN PHASE'',I6)') KKK,IPH
        IF(IERROR.EQ.2) WRITE(*,'('' ITMAX WAS REACH IN NEWTRAPH -->'',
     #       '' CANNOT SOLVE GRAIN'',I6,'' IN PHASE'',I6)') KKK,IPH
        WRITE(*,'('' THE INPUT STRESS IS RETAINED'')')
        DO I=1,5
          X(I)=XORI(I)
        ENDDO
      ENDIF

      DO I=1,5
        SG(I,KKK)=X(I)
      ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
      SUBROUTINE INIT_CRSS_VOCE_AS

      INCLUDE 'vpsc_as.dim'

        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO IS=1,NSYST(IPHEL)
              CRSS(IS,KKK)=TAU(IS,0,IPHEL)
            ENDDO
          ENDDO
        ENDDO
c
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C **************************************************************************
C     SUBROUTINE INITIAL_STATE_GUESS     --->     VERSION 07/DEC/05

C     IF STRAIN IS IMPOSED, USES A FIRST STRESS GUESS COLINEAR WITH THE
C     STRAIN RATE AND CALCULATES GRAIN STRESS GIVEN BY TAYLOR.
C     IF STRESS IS IMPOSED, SETS THE GRAIN STRESS EQUAL TO THE MACROSCOPIC.
C     CALCULATES GRAIN STRAIN RATES 'DG' FOR EITHER CASE.
C     CALCULATES AVERAGE STRESS 'SAV' AND STRAIN-RATE 'DAV'.
C     CALCULATES VISCO-PLASTIC MODULI 'XMCTG' FOR EVERY GRAIN.
C     CALCULATES AN INITIAL GUESS FOR THE MACROSCOPIC VISCO-PLASTIC MODULUS.
C **************************************************************************

      SUBROUTINE INITIAL_STATE_GUESS

      INCLUDE 'vpsc_as.dim'

      integer flag
      integer diagnostics

      DIMENSION AUXTAN(5,5)

      flag = 0
      diagnostics = 0
!as      DBARNORM=TNORM(DBAR,5,1)
      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
!as        IF(STRAIN_CONTROL.EQ.1) THEN
!as            DO J=1,5
!as              STRY(J,KGX)=DBAR(J)/DBARNORM
!as            ENDDO
!as            CALL GRAIN_STRESS (0,KGX,KKK,IPHEL,IPH)
!as          ELSE IF(STRAIN_CONTROL.EQ.0) THEN
            DO J=1,5
              SG(J,KKK)=SBAR(J)
            ENDDO
!as          ENDIF
          CALL GRAIN_RATE_AND_MODULI (0,1,KGX,KKK,IPHEL,IPH)
          KGX=KGX+1
        ENDDO
      ENDDO

      DO I=1,5
        DAV(I)=0.
        SAV(I)=0.
        DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)
          DAV(I)=DAV(I)+DG(I,KKK)*WGT(KKK)
          SAV(I)=SAV(I)+SG(I,KKK)*WGT(KKK)
        ENDDO
        ENDDO
      ENDDO

!     CALCULATE INITIAL GUESS FOR MACROSCOPIC MODULI 'Mtg' AS THE INVERSE
!     OF THE AVERAGE OF THE GRAIN'S STIFFNESSES

      DO I=1,5
        dzero(i)=0.
        DO J=1,5
          XMTG(I,J)=0.
        ENDDO
      ENDDO

      KGX=1
      DO IPH=IPHBOT,IPHTOP
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        DO I=1,5
        DO J=1,5
          AUXTAN(I,J)=XMCTG(I,J,KGX)
        ENDDO
        ENDDO
        KGX=KGX+1

         flag = 0
         CALL LU_INVERSE(AUXTAN,5,flag,diagnostics)
         if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'AUXTAN is singular'
               DO I=1,5
                  write(*,*) (AUXTAN(I,J),j=1,5)
               ENDDO
            endif
            auxtan=0.0
         endif

         DO I=1,5
            dzero(i)=dzero(i)+dczero(i,kgx)*wgt(kkk)
            DO J=1,5
               XMTG(I,J)=XMTG(I,J)+AUXTAN(I,J)*WGT(KKK)
            ENDDO
         ENDDO
      ENDDO
      ENDDO
cw
      DO I=1,5
         DO J=1,5
            XLTG(I,J)=XMTG(I,J)
         ENDDO
      ENDDO

      flag = 0
      CALL LU_INVERSE(XMTG,5,flag,diagnostics)
        if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'XMTG is singular'
               DO I=1,5
                  write(*,*) (xmtg(I,J),j=1,5)
               ENDDO
            endif
        endif


c      write( *,*)
c      write( *,'('' inside initial_state_guess'')')
c      write( *,'('' sav='',5e12.3)') sav
c      write( *,'('' dav='',5e12.3)') dav
c      write( *,'('' xmtg ='',5e12.3)') xmtg
c      write( *,*)
c      pause

      RETURN
      END
c
c------------------------------------------------------------------------------
c
      SUBROUTINE LOAD_CONDITIONS_AS
c
      INCLUDE 'vpsc_as.dim'
c
      do I=1,6
         idsim(I)=0
         iscau(I)=1
      enddo
      strain_control=0
c
      RETURN
      END
c
c------------------------------------------------------------------------------
c
      SUBROUTINE LU_INVERSE (A,N,FLAG,diagnostics)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)
      
      integer flag
      integer diagnostics

C *** INVERTS A MATRIX USING LU DECOMPOSITION
      INTEGER N,ISINGULAR
      REAL*8 A,Y,INDX
C      DIMENSION A(N,N),Y(N,N),INDX(N)      ! MAY CHOKE SOME COMPILERS
      DIMENSION A(5,5),Y(5,5),INDX(5)
      
      flag=0

c     write(*,*) 'A(i,j) matrix inside lu_inverse'
c     write(*,'(5e12.3)') ((a(i,j),j=1,n),i=1,n)
c     pause

C **************************************************************
C *** BLOCK ADDED 03/DEC/05 TO AVOID NUMERICALLY SINGULAR MATRIX
      AMAX=0.D0
      DO I=1,N
         DO J=1,N
            DUM=ABS(A(I,J))
            IF(DUM .GT. AMAX) AMAX=DUM
         ENDDO
      ENDDO
      DO I=1,N
         DO J=1,N
            A(I,J)=A(I,J)/AMAX      ! normalize the matrix
         ENDDO
      ENDDO
C **************************************************************

      DO I=1,N
         DO J=1,N
            Y(I,J)=0.
         ENDDO
         Y(I,I)=1.
      ENDDO

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
      IF(ISINGULAR.EQ.1) THEN
         if (diagnostics.eq.1) then
            WRITE(*,*) ' *** SINGULAR MATRIX IN LU_INVERSE !!'
            write(*,*) 'A(i,j) matrix inside lu_inverse'
            write(*,'(5e12.3)') ((a(i,j),j=1,n),i=1,n)
         endif
         flag=1
         return
        !STOP
      ENDIF

      DO J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
      ENDDO

      DO I=1,N
         DO J=1,N
            A(I,J)=Y(I,J) /AMAX      ! renormalize the inverse
         ENDDO
      ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
      SUBROUTINE LU_EQSYSTEM(A,B,N,ISINGULAR)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)
      INTEGER N,ISINGULAR
      REAL*8 A,INDX
C *** SOLVES A*X=B USING LU DECOMPOSITION

      DIMENSION A(N,N),B(N),INDX(N)      ! MAY CHOKE SOME COMPILERS
!     DIMENSION A(5,5),B(5),INDX(5)

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)

      IF(ISINGULAR.EQ.1) RETURN

      CALL LUBKSB(A,N,N,INDX,B)

      RETURN
      END
c
c------------------------------------------------------------------------------
c
      SUBROUTINE LU_EQSYS25(A,B,N,ISINGULAR)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

C *** SOLVES A*X=B USING LU DECOMPOSITION

c      DIMENSION A(N,N),B(N),INDX(N)      ! MAY CHOKE SOME COMPILERS
      DIMENSION A(25,25),B(25),INDX(25)

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
c
      do j=1,25
         d=d*a(j,j)
      enddo

c      write(*,*) 'DET25=',d
c      pause

      IF(ISINGULAR.EQ.1) RETURN

      CALL LUBKSB(A,N,N,INDX,B)

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C ***************************************************************************
C
C     SUBROUTINE NEWTON_RAPHSON (ex SNLNR)   --->   VERSION APR/2003
C
C     GIVEN AN INPUT GRAIN STRESS 'X' SOLVES THE VISCOPLASTIC EQUATION
C     USING NEWTON-RAPHSON LINEARIZATION AND ITERATING UNTIL CONVERGENCE.
C     - Fi(X)   : ARE THE FUNCTIONS TO MINIMIZE.
C     - FGRADij : ARE THE DERIVATIVES OF Fi WITH RESPECT TO Xj.
C     CORRECTIONS TO 'X' AND RELAXATION OF TOLERANCE ELIMINATED (FEB/2000).
C     PASSING ALL RELEVANT CRYSTAL ARRAYS THROUGH COMMON (13/DEC/02).
C     ADDED RELAXED CONSTRAINTS OPTION (IRC=1) (APR/03)
C
C ***************************************************************************

      SUBROUTINE NEWTON_RAPHSON (X,IRC,KMAX,EPS,TAULIM,IERROR ,iprx)

      INCLUDE 'vpsc_as.dim'

      DIMENSION FGRAD(5,5),F(5),FGRADX(3,3),FX(3),X(5),XOLD(5)
      DIMENSION RSS(NSYSMX),GD(NSYSMX)
      COMMON/NLNR /DB(5),XMASTX(5,5),SCX(5,NSYSMX),TAUX(NSYSMX),
     #             GAMD0X(NSYSMX),NSYSTX,NRSX(NSYSMX),ISENSEX(NSYSMX)

      COEF=0.2
      IERROR=0

      DO 1000 K=1,KMAX

C *** RESOLVED SHEARS CALCULATION & OUTSIDE Y.S. ERROR MANAGING BLOCK.
C *** NRS MAY BE EVEN OR ODD:
C     GD ALWAYS > 0 --> GD IS USED TO GET DERIVATIVES TO BUILD THE
C     COEFFICIENT MATRIX FOR N-R METHOD IN F(I) CALCULATION (INDEPENDENT
C     TERM FOR N-R METHOD)   RSS*GD=GAMDOT -> SIGN(RSS*GD)=SIGN(RSS)

        DO IS=1,NSYSTX
          RSS(IS)=SCX(1,IS)*X(1)+SCX(2,IS)*X(2)+SCX(3,IS)*X(3)
     #           +SCX(4,IS)*X(4)+SCX(5,IS)*X(5)
          IF(.NOT.(RSS(IS).GT.0. .OR. ISENSEX(IS).EQ.1)) RSS(IS)=0.
          RSS(IS)=RSS(IS)/TAUX(IS)
          IF(ABS(RSS(IS)).LT.1.E-10) RSS(IS)=0.
c
cpp          write(*,*) 'X=',x
cpp          write(*,*) 'K,IS,RSS(IS)=',k,is,rss(is)
cpp          pause
c
          GD(IS)=GAMD0X(IS)
          IF(NRSX(IS).NE.1) GD(IS)=GD(IS)*ABS(RSS(IS)**(NRSX(IS)-1))

          IF(ABS(RSS(IS)).GT.TAULIM) THEN
            DO I=1,5
              X(I)=XOLD(I)+COEF*(X(I)-XOLD(I))
            ENDDO
            GO TO 1000
          ENDIF
        ENDDO

cw        if(iprx.eq.1) then
cpp        write(*,'(''db & rss='',5f10.4)') db
cpp        write(*,'(12f10.4)') (rss(is),is=1,nsystx)
cpp        write(*,'(''nrsx & xmast='',30i3)') (nrsx(is),is=1,nsystx)
cpp        write(*,'(5f13.7)') ((xmastx(i,j),j=1,5),i=1,5)
cw        endif

        DO I=1,5
          F(I)=-DB(I)
          DO J=1,5
            F(I)=F(I)+XMASTX(I,J)*X(J)
          ENDDO
        ENDDO
        DO I=1,5
          DO IS=1,NSYSTX
            F(I)=F(I)+SCX(I,IS)*RSS(IS)*GD(IS)
          ENDDO
        ENDDO
        DO I=1,5
          DO J=1,5
            FGRAD(I,J)=-XMASTX(I,J)
          ENDDO
        ENDDO

        DO IS=1,NSYSTX
          SFACTOR=NRSX(IS)*GD(IS)/TAUX(IS)
          DO I=1,5
          DO J=1,5
            FGRAD(I,J)=FGRAD(I,J)-SFACTOR*SCX(I,IS)*SCX(J,IS)
          ENDDO
          ENDDO
        ENDDO

        DO I=1,5
          XOLD(I)=X(I)
        ENDDO

C *** SOLVES LINEAR SYSTEM
        IF(IRC.EQ.0) THEN

cw          if(iprx.eq.1) then
cpp          write(*,'(''f    ='',5e12.3)') f
cpp          write(*,'(''fgrad='',5e12.3,/,(6x,5e12.3))') fgrad
cw          endif

          CALL LU_EQSYSTEM (FGRAD,F,5,ISINGULAR)
          IF(ISINGULAR.EQ.1) THEN
            IERROR=1
            write(*,*) 'Failed here with fGrad ='
            do i=1,5
            write(*,'(5f10.5)') (fgrad(i,j),j=1,5)
            enddo
            RETURN
          ENDIF

C *** BOUNDS THE STRESS CORRECTION TO AVOID LARGE OSCILATIONS IN CONVERGENCE
cas          if(iprx.eq.1) then
cas          rcorr=TNORM(F,5,1)/TNORM(XOLD,5,1)
cas          write(10,'(''*** NR correction'',5f9.3,f12.5)')
cas     #                 (f(i),i=1,5),rcorr
cas          endif

          DO I=1,5
            X(I)=XOLD(I)+F(I)
          ENDDO
        ENDIF

        IF(IRC.EQ.1) THEN          ! RELAXED CONSTRAINTS CASE
          FX(1)=F(1)
          FX(2)=F(2)
          FX(3)=F(5)
          FGRADX(1,1)=FGRAD(1,1)
          FGRADX(1,2)=FGRAD(1,2)
          FGRADX(1,3)=FGRAD(1,5)
          FGRADX(2,1)=FGRAD(2,1)
          FGRADX(2,2)=FGRAD(2,2)
          FGRADX(2,3)=FGRAD(2,5)
          FGRADX(3,1)=FGRAD(5,1)
          FGRADX(3,2)=FGRAD(5,2)
          FGRADX(3,3)=FGRAD(5,5)
          CALL LU_EQSYSTEM (FGRADX,FX,3,ISINGULAR)
          IF(ISINGULAR.EQ.1) THEN
            IERROR=1
            RETURN
          ENDIF
          X(1)=XOLD(1)+FX(1)
          X(2)=XOLD(2)+FX(2)
          X(3)=0.
          X(4)=0.
          X(5)=XOLD(5)+FX(3)
        ENDIF

        RERROR=TMISMATCH(X,XOLD,5,1)
        IF(RERROR.LT.EPS) RETURN

1000  CONTINUE      ! END OF MASTER DO
C *******************************************************************

      IERROR=2

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *************************************************************************
C     SUBROUTINE SCALE_3   ---->   VERSION 26/MAY/2000
C
C     THIS SUBROUTINE IS MEANT TO BE USED IN COMBINATION WITH MTS MODEL
C     WHERE RATE EFFECTS ARE ACCOUNTED FOR IN THE FUNCTIONAL FORM OF CRSS.
C
C     BY REDEFINING THE FACTOR 'GAMD0' OF THE ORDER OF THE STRAIN-RATE
C     THE RATIO (RSS/CRSS) WILL BE OF ORDER ONE AND THERE ARE NO RATE
C     SENSITIVITY EFFECTS ASSOCIATED WITH THE N'th POWER.
C     THE N'th POWER IS ONLY A WAY TO AVOID AMBIGUITIES AND IS UNIQUE.
C *************************************************************************

      SUBROUTINE SCALE_3 (ISTEP)

cas      INCLUDE 'vpsc7.dim'
      INCLUDE 'vpsc_as.dim'

C     DIMENSION DGX(5)

      REFRATE=TNORM(DBAR,5,1)
C     REFRATE=DVM      ! REFERENCE USED BY KOK & BEAUDOIN

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        DO KKK=NGR(IPH-1)+1,NGR(IPH)

C          IF(ISTEP.GT.1) THEN      ! DG(I,KGX) NOT STORED FOR MULTIELEM CALC
C            DO I=1,5
C              DGX(I)=DG(I,KGX)
C            ENDDO
C            REFRATE=TNORM(DGX,5,1)
C          ENDIF

          GAMD0G(KGX)=REFRATE
          KGX=KGX+1
        ENDDO
      ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *****************************************************************************
C     SUBROUTINE STATE_5x5      --->      VERSION 26/jan/2000
C
C     BASED ON SECANT EQUATION:   D(i)=MS(i,j)*S(j)
C     SHIFTS KNOWN COMPONENTS TO INDEPENDENT TERM, CALCULATES COEFFICIENTS
C     FOR UNKNOWN COMPONENTS AND INVERTS THE LINEAR SYSTEM.
C *****************************************************************************

      subroutine state_5x5 (ibc_d,BC_D,ibc_s,BC_S,XMS)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      dimension aux11(5),aux21(5,5)
      dimension XMS(5,5)
      dimension ibc_d(5),bc_d(5),ibc_s(5),bc_s(5)

      do i=1,5
        aux11(i)=-1.d0*ibc_d(i)*bc_d(i)
        do j=1,5
          aux11(i)=aux11(i)+xms(i,j)*ibc_s(j)*bc_s(j)
          aux21(i,j)=ibc_s(j)*(i/j)*(j/i)-ibc_d(j)*xms(i,j)
        enddo
      enddo

      CALL LU_EQSYSTEM(AUX21,AUX11,5,IER)

      if(ier.eq.1) then
      write(*,*) 'SINGULAR SYSTEM IN STATE_5X5'
      stop
      endif

      do i=1,5
        bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux11(i)
        bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux11(i)
      enddo

      return
      end
c
c------------------------------------------------------------------------------
c
C
C *****************************************************************************
C     SUBROUTINE STATE_6x6      --->      VERSION 16/NOV/97
C
C     BASED ON SECANT EQUATION:   D(i)=MS(i,j)*S(j)
C     SHIFTS KNOWN COMPONENTS TO INDEPENDENT TERM, CALCULATES COEFFICIENTS
C     FOR UNKNOWN COMPONENTS AND INVERTS THE LINEAR SYSTEM.
C *****************************************************************************
c
      subroutine state_6x6 (ibc_d,AUX_D,ibc_s,AUX_S,AUX_MS,flag)
c
      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      dimension aux4(3,3,3,3),aux2(6,6),aux11(6),aux21(6,6)
      dimension AUX_D(3,3),AUX_S(3,3),AUX_MS(5,5)
      dimension ibc_d(6),bc_d(6),ibc_s(6),bc_s(6)
      dimension profac(6)
      dimension aux5(5),aux33(3,3)

      integer flag

      flag = 0
c
      do i=1,6
        profac(i)=1.d0+(i/4)
      enddo
C
C *** PUTS CARTESIAN 'D' & 'S' IN VOIGT.
        CALL VOIGT (BC_D,AUX_D,AUX2,AUX4,2)
        CALL VOIGT (BC_S,AUX_S,AUX2,AUX4,2)
C *** PUTS b-BASIS 'MS' IN VOIGT.
        CALL CHG_BASIS(AUX5,AUX33,AUX_MS,AUX4,3,5)
        CALL VOIGT (AUX11,AUX33,AUX2,AUX4,4)

      do i=1,6
        aux11(i)=-1.d0*ibc_d(i)*bc_d(i)
        do j=1,6
          aux11(i)=aux11(i)+aux2(i,j)*ibc_s(j)*bc_s(j)*profac(j)
          aux21(i,j)=ibc_s(j)*(i/j)*(j/i)-ibc_d(j)*aux2(i,j)*profac(j)
        enddo
      enddo

      CALL LU_EQSYSTEM(AUX21,AUX11,6,IER)

      if(ier.eq.1) then
         write(*,*) 'SINGULAR SYSTEM IN STATE_6X6'
         flag = 1
         write(*,*) 'aux11 ='
         write(*,'(6e12.3)') (aux11(j),j=1,6)
         write(*,*) 'aux21 ='
         write(*,'(6e12.3)') ((aux21(i,j),j=1,6),i=1,6)
         return
         !stop
      endif

      do i=1,6
        bc_d(i)=ibc_d(i)*bc_d(i)+ibc_s(i)*aux11(i)
        bc_s(i)=ibc_s(i)*bc_s(i)+ibc_d(i)*aux11(i)
      enddo
      CALL VOIGT (BC_D,AUX_D,AUX2,AUX4,1)
      CALL VOIGT (BC_S,AUX_S,AUX2,AUX4,1)

      return
      end
c
c------------------------------------------------------------------------------
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NRowsxNCols MATRICES
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
C     OF BOTH DATA.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION TMISMATCH (v1,v2,NR,NC)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION v1(nr*nc),v2(nr*nc)
      DIMENSION v_dif(nr*nc),v_ave(nr*nc)
        do i=1,NR*NC
          v_dif(i)=v1(i)-v2(i)
          v_ave(i)=0.5*(v1(i)+v2(i))
        enddo
        tmismatch =tnorm (v_dif,NR,NC)/tnorm (v_ave,NR,NC)
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NRows,NCols =< 6)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION TNORM (V,NR,NC)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION v(nr*nc)
        tnorm =0.
        do i=1,NR*NC
          tnorm =tnorm +v(i)*v(i)
        enddo
        tnorm =sqrt(tnorm)
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TMISMATCH5   ---->   VERSION OF 02/JAN/03
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NRowsxNCols MATRICES
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
C     OF BOTH DATA.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION TMISMATCH5(v1,v2,NR,NC)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION v1(25),v2(25)
      DIMENSION v_dif(25),v_ave(25)
        do i=1,NR*NC
          v_dif(i)=v1(i)-v2(i)
          v_ave(i)=0.5*(v1(i)+v2(i))
        enddo
        tmismatch5=tnorm5(v_dif,NR,NC)/tnorm5(v_ave,NR,NC)
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     FUNCTION TNORM5   ---->   VERSION OF 02/JAN/03
C
C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NRows,NCols =< 5)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION TNORM5(V,NR,NC)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION v(25)
        tnorm5=0.
        do i=1,NR*NC
          tnorm5=tnorm5+v(i)*v(i)
        enddo
        tnorm5=sqrt(tnorm5)
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *****************************************************************************
C     SUBROUTINE TWIN_ORIENTATION    --->    VERSION OF 04/dec/02
C
C     GIVEN THE BURGERS VECTOR 'BUR' OF THE TWIN SYSTEM (IN CRYSTAL AXES)
C     MAKES A ROTATION OF 180 DEG AROUND IT TO DEFINE THE ORIENTATION OF
C     THE TWIN RELATED CRYSTAL.
C     THE MATRIX 'ATWIN' TRANSFORMS FROM TWIN TO CRYSTAL AXES.
C *****************************************************************************

      SUBROUTINE TWIN_ORIENTATION (BUR,ATWIN)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION BUR(3),HPI(3,3),ATWIN(3,3),AUX(3,3)

      DATA HPI/-1.,0.,0.,0.,-1.,0.,0.,0.,1./
      PI=4.*ATAN(1.D0)

      ANG1=ATAN2(BUR(2),BUR(1))*180./PI+90.
      ANG2=SQRT(BUR(1)**2+BUR(2)**2)
      ANG2=ATAN2(ANG2,BUR(3))*180./PI
      CALL EULER (2,ANG1,ANG2,0.0+0,AUX)

      DO I=1,3
      DO J=1,3
        ATWIN(I,J)=0.
        DO K1=1,3
        DO K2=1,3
          ATWIN(I,J)=ATWIN(I,J)+AUX(K1,I)*HPI(K1,K2)*AUX(K2,J)
        ENDDO
        ENDDO
      ENDDO
      ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_VOCE     --->      VERSION OF 28/dec/00
C
C     A VOCE LAW FUNCTION OF THE ACCUMULATED SHEAR IN EACH GRAIN IS ADDED
C     AS A MULTIPLICATIVE FACTOR THAT MODULATES THE ORIGINAL LINEAR HARDENING.
C     THE UNITS (AND THE STRENGTH) OF THE HARDENING ARE CARRIED BY THE
C     MULTIPLICATIVE FACTOR 'VOCE'.
C     THE SELF & LATENT COUPLING COEFFICIENTS 'HARD' ARE DIMENSIONLESS
C     CONSTANTS RELATIVE TO THE FACTOR 'VOCE'.
C     THE INCREMENT OF CRSS IS GIVEN BY AN ANALYTIC INTEGRAL EXPRESSION OF
C     THE VOCE MODULATING FUNCTION, RATHER THAN USING FORWARD EXTRAPOLATION
C     AS WAS DONE BEFORE 28/SET/00.
C     THE PARAMETERS IN VOCE EXPRESSION MAY ADOPT 'NON-KOSHER' VALUES (DEC/00)
*******************************************************************************

cas      SUBROUTINE UPDATE_CRSS_VOCE (IOPTION)
      SUBROUTINE UPDATE_CRSS_VOCE_AS (IOPTION,STRENGTHS,GAMMAS,
     #                            STRENGTHRATES,GAMMARATES
     #                       ,itnStrngthMax,itnphase)
cass
cas
      INCLUDE 'vpsc_as.dim'
cas
cRL      dimension strengths(:,:),gammas(:),strengthrates(:,:),gammarates(:)
cass      dimension strengths(:,:),gammas(:)
cass      dimension strengthrates(:,:),gammarates(:)
cass
      dimension strengths(itnStrngthMax,itnphase)
      dimension gammas(itnphase)
      dimension strengthRates(itnStrngthMax,itnphase)
      dimension gammarates(itnphase)
cass
cas
cas   this is now done in init_crss_voce_as
cas
cas      IF (IOPTION.EQ.1) THEN
cas
cas        DO IPH=IPHBOT,IPHTOP
cas          IPHEL=IPH-IPHBOT+1
cas          DO KKK=NGR(IPH-1)+1,NGR(IPH)
cas            DO IS=1,NSYST(IPHEL)
cas              CRSS(IS,KKK)=TAU(IS,0,IPHEL)
cas            ENDDO
cas          ENDDO
cas        ENDDO
cas
      IF (IOPTION.EQ.1) THEN

        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
c
cRL
cRL            DO IS=1,NSYST(IPHEL)
cRL              CRSS(IS,KKK)=
cRL     #        STRENGTHS(IPHEL)*TAU(IS,0,IPHEL)/TAU(1,0,IPHEL)
cRL            ENDDO
cRL
          js=0
           do im=1,nmodes(iphel)
            do is=1,nsm(im,iphel)
             js=js+1
cass             write (*,*)'IM,IS,JS,KKK,IPHEL=',im,is,js,kkk,iphel
             CRSS(JS,KKK)=STRENGTHS(IM,IPHEL)
            enddo
           enddo
c
          ENDDO
        ENDDO

      ELSE IF (IOPTION.EQ.2) THEN

        KGX=1
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
cas
          GAMTOTX=GAMMAS(IPHEL)
cas
          DO KKK=NGR(IPH-1)+1,NGR(IPH)

cas          GAMTOTX=GTOTGR(KKK)
c     as          DELTGAM=0.0
             DO IS=1,NSYST(IPHEL)
cas            DELTGAM=DELTGAM+ABS(GAMDOT(IS,KGX))*TINCR
                GAMMARATES(IPHEL)=
     #GAMMARATES(IPHEL)+ABS(GAMDOT(IS,KGX))*WGT(KKK)
cas
             ENDDO
cRL
cRL          DO IS=1,NSYST(IPHEL)
cRL
             is=0
             do im=1,nmodes(iphel)
c
                STRENGTHRATES(IM,IPHEL)=0.
c
                do ks=1,nsm(im,iphel)
                   is=is+1
c
                   DTAU=0.
                   DO JS=1,NSYST(IPHEL)
c     as              DTAU=DTAU+HARD(IS,JS,IPHEL)*ABS(GAMDOT(JS,KGX))*TINCR
                      DTAU=DTAU+HARD(IS,JS,IPHEL)*ABS(GAMDOT(JS,KGX))
                   ENDDO
                   TAU0 =TAU (IS,0,IPHEL)
                   TAU1 =TAU (IS,1,IPHEL)
                   THET0=THET(IS,0,IPHEL)
                   THET1=THET(IS,1,IPHEL)
                   TINY=1.E-4*TAU0
                   
                   VOCE=0.0
                   IF(ABS(THET0).GT.TINY) THEN
c     as              VOCE=THET1*DELTGAM
                      VOCE=THET1
                      IF(ABS(TAU1).GT.TINY) THEN
                         FACT=ABS(THET0/TAU1)
                         EXPINI=EXP(-GAMTOTX*FACT)
c     as                EXPDEL=EXP(-DELTGAM*FACT)
c     as                VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*
c     as     #            (EXPDEL-1.)-THET1/FACT*EXPINI*
c     as     #            (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
                         
                         VOCE=VOCE+
     #(FACT*TAU1-THET1)*EXPINI+FACT*THET1*GAMTOTX*EXPINI
                      ENDIF
                   ENDIF
c     as            CRSS(IS,KKK)=CRSS(IS,KKK)+DTAU*VOCE/DELTGAM
c     as
c     as      The expression below has to change if we consider 
c     as      as many strengths per as modes per phase
c     as
c     as            STRENGTHRATES(IPHEL)=
c     as     #      STRENGTHRATES(IPHEL)+DTAU*VOCE*WGT(KKK)/NSYST(IPHEL)
c     as
                   STRENGTHRATES(IM,IPHEL)=
     #STRENGTHRATES(IM, IPHEL)+DTAU*VOCE*WGT(KKK)/NSM(IM,IPHEL)
c     as
                enddo           ! ks enddo
             enddo              ! im enddo
c     
             KGX=KGX+1
          ENDDO
       ENDDO
       
      ENDIF
      
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C *****************************************************************************
C     SUBROUTINE UPDATE_FIJ      --->      VERSION OF NOV/28/2005
C
C     USES THE VELOCITY GRADIENT (AVERAGE, PHASE or GRAIN) IN THE STEP
C     TO UPDATE INCREMENTALLY THE CORRESPONDING DEFORMATION TENSOR 'FIJ'
C
C     REPLACES PREVIOUS SUBR. UPDFIJ_AVERAGE & SUBR. UPDFIJ_LOCAL  (SEP/04)
C *****************************************************************************

      SUBROUTINE UPDATE_F2IJDOT_AS (IPH)

      INCLUDE 'vpsc_as.dim'

      DIMENSION FNEW(3,3)


C *** UPDATES THE DEFORM GRAD IN THE ELEMENT 'FIJPH(i,j,0)' USING THE
C     MACROSCOPIC VELOCITY GRADIENT 'UDOT'
C *** UPDATES THE DEFORM GRAD IN EACH PHASE 'FIJPH(i,j,IPH)' USING THE
C     AVERAGE VELOCITY GRADIENT FOR THE PHASE 'LIJPH' CALCULATED INSIDE
C     SUBROUTINE UPDATE_ORIENTATION.
C *** LIJPH ACCOUNTS FOR ROTATIONS BUT NOT FOR STRETCH WHEN IFLAT(iph)=1.
C *** FIJPH COINCIDES WITH FIJ OF ELEMENT IF NPH=1 AND IFLAT(1)=0.

      DO I=1,3
         DO J=1,3
cas        FNEW(I,J)=0.0
cas
      aux1=0.
      aux2=0.
cas
        IF(IPH.EQ.0) THEN
          DO K=1,3
cas            FNEW(I,J)=FNEW(I,J)+(TINCR*UDOT(I,K)+XID3(I,K))*FIJPH(K,J,0)
cas
            aux1=aux1+UDOT(I,K)*F2IJPH(K,J,0)
            aux2=aux2+F2IJPH(I,K,0)*UDOT(J,K)
cas
          ENDDO
        ELSE IF(IPH.GT.0) THEN
          DO K=1,3
cas            FNEW(I,J)=FNEW(I,J)+(TINCR*LIJPH(I,K,IPH)+XID3(I,K))
cas     #                         *FIJPH(K,J,IPH)
cas
      if(interaction.ne.5) then
            aux1=aux1+LIJPH(I,K,IPH)*F2IJPH(K,J,0)
            aux2=aux2+F2IJPH(I,K,0)*LIJPH(J,K,IPH)
      else
            aux1=aux1+UDOT(I,K)*F2IJPH(K,J,0)
            aux2=aux2+F2IJPH(I,K,0)*UDOT(J,K)
      endif
cas
          ENDDO
        ENDIF
cas
      FNEW(I,J)=aux1+aux2
cas
      ENDDO
      ENDDO
      DO I=1,3
      DO J=1,3
        F2IJDOTPH(I,J,IPH)=FNEW(I,J)
      ENDDO
      ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C **************************************************************************
C     SUBROUTINE UPDATE_ORIENTATION    --->   VERSION 28/MAR/2007
C
C     UPDATES GRAIN ORIENTATIONS DUE TO CRYSTALLOGRAPHIC SHEAR, BUT DOES
C     NOT PERFORM TWIN REORIENTATION.
C     UPDATES GRAIN AND PHASE DISTORTION TENSORS REQUIRED FOR UPDATING THE
C     PHASE AND GRAIN ELLIPSOID IN SUBROUTINES UPDATE_FIJ & UPDATE_SHAPE.
C
C     CNT: modified to use only phase-associated tensors when ISHAPE.LE.1
C     RAL: split in 3 DO loops to deal with co-rotations (17/02/00)
C     CNT: moved last loop to UPDATE_TWINNING (05/06/02)
C **************************************************************************

cas      SUBROUTINE UPDATE_ORIENTATION
cas
cass      SUBROUTINE UPDATE_ORIENTATION_AS (REORRATES)
      SUBROUTINE UPDATE_ORIENTATION_AS (REORRATES,kgrtot)
cas

      INCLUDE 'vpsc_as.dim'
cas
      dimension reorrates(3,kgrtot)
cas
      DIMENSION as(3,3,3,3),east(3,3),aa(3,3)
      dimension dnsa(3),dbsa(3),rotslip(3,3),rotloc(3,3)
      dimension aux5(5),aux55(5,5),aux3333(3,3,3,3)
cas      dimension rot(3,3,NGRMX),arot(3,3),corot(3,3)
      real*8      LIJGRX(3,3),LIJGR0(3,3)

cas      RSLBAR=0.
cas      RLCBAR=0.

      DO I=1,3
        DO J=1,3
           ROTBAR_AS(I,J)=0
        ENDDO
      ENDDO
      

      KGX=1
      DO 2000 IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1

        DO I=1,3
        DO J=1,3
          LIJPH(I,J,IPH)=0.
        ENDDO
        ENDDO

        IF(ISHAPE(IPH).LE.1) THEN
          DO I=1,3
          DO J=1,3
          DO K=1,3
          DO L=1,3
            AS(I,J,K,L)=ASPH(I,J,K,L,IPH)
          ENDDO
          ENDDO
          ENDDO
          ENDDO
        ENDIF

casC *** THE FOLLOWING TO AVOID EMPTY 'CHILD' PHASE IN COMP GRAIN MODEL
cas        IF(WPH(IPH).LT.1.E-6) GO TO 2000

      DO 1000 KKK=NGR(IPH-1)+1,NGR(IPH)

C *** CALCULATES LOCAL ROTATION ROTLOC=PI*S**(-1)*(DG-DAV) FOR EVERY GRAIN.
C *** ROTLOC IS ZERO FOR TAYLOR CALCULATION.

        IF(INTERACTION.LE.0) THEN
          DO I=1,3
          DO J=1,3
            ROTLOC(I,J)=0.
          ENDDO
          ENDDO
        ELSE IF(INTERACTION.GT.0) THEN
CFEB
cas          IF(ISHAPE(IPH).GE.2) THEN
cas            DO I=1,3
cas            DO J=1,3
cas            DO K=1,3
cas            DO L=1,3
cas              AS(I,J,K,L)=ASGR(I,J,K,L,KKK)
cas            ENDDO
cas            ENDDO
cas            ENDDO
cas            ENDDO
cas          ENDIF
CFEE
          DO I=1,5
            AUX5(I)=DG(I,KGX)-DAV(I)
          ENDDO

          CALL CHG_BASIS(AUX5,EAST,AUX55,AUX3333,1,5)

          DO I=1,3
          DO J=1,3
          ROTLOC(I,J)=0.
            DO K=1,3
            DO L=1,3
              ROTLOC(I,J)=ROTLOC(I,J)+AS(I,J,K,L)*EAST(K,L)
            ENDDO
            ENDDO
          ENDDO
          ENDDO

        ENDIF

C *** CALCULATES VELOCITY GRADIENT IN EACH PHASE AND EACH GRAIN

        do i=1,3
        do j=1,3
          aa(i,j)=ag(i,j,KKK)
          LIJGR0(I,J)=0.
        enddo
        enddo

        do is=1,nsyst(IPHEL)
          do i=1,3
            dnsa(i)=0.
            dbsa(i)=0.
            do j=1,3
              dnsa(i)=dnsa(i)+aa(i,j)*dnca(j,is,IPHEL)
              dbsa(i)=dbsa(i)+aa(i,j)*dbca(j,is,IPHEL)
            enddo
          enddo
          do i=1,3
          do j=1,3
            LIJGR0(i,j)=LIJGR0(i,j)+dbsa(i)*dnsa(j)*gamdot(is,KGX)
          enddo
          enddo
        enddo

        DO I=1,3
        DO J=1,3
cas          LIJGRX(I,J)=ROTBAR(I,J)+ROTLOC(I,J)
          LIJGRX(I,J)=ROTLOC(I,J)
          IF(IFLAT(IPH).EQ.0) THEN
            LIJGRX(I,J)=LIJGRX(I,J)+(LIJGR0(I,J)+LIJGR0(J,I))/2.
          ENDIF
          LIJPH(I,J,IPH)=LIJPH(I,J,IPH)+LIJGRX(I,J)*WGT(KKK)/WPH(IPH)
        ENDDO
        ENDDO
CFEB
cas        IF(ISHAPE(IPH).GT.0) THEN
cas         DO I=1,3
cas          DO J=1,3
cas            LIJGR(I,J,KKK)=LIJGRX(I,J)
cas          ENDDO
cas          ENDDO
cas        ENDIF
CFEE
cas
cas     Below we are considering that shears in BOTH slip and twinning
cas     systems contribute to the reorientation rate of the grain. See
cas     also comment in update_twinning_as
cas
        DO I=1,3
        DO J=1,3
          ROTSLIP(I,J)=(LIJGR0(I,J)-LIJGR0(J,I))/2.
        ENDDO
        ENDDO

cas        DO I=1,3
cas        DO J=1,3
cas          ROT(I,J,KGX)=(ROTBAR(I,J)+ROTLOC(I,J)-ROTSLIP(I,J))*TINCR
cas        ENDDO
cas        ENDDO

          REORRATES(1,KGX)=ROTLOC(2,3)-ROTSLIP(2,3)
          REORRATES(2,KGX)=ROTLOC(3,1)-ROTSLIP(3,1)
          REORRATES(3,KGX)=ROTLOC(1,2)-ROTSLIP(1,2)

casC *** AVERAGE SLIP AND LOCAL ROTATION (FOR STATISTICAL PURPOSES ONLY)
cas        RSLBAR=RSLBAR+SQRT(ROTSLIP(3,2)**2+ROTSLIP(1,3)**2+
cas     #         ROTSLIP(2,1)**2)*WGT(KKK)
cas        RLCBAR=RLCBAR+SQRT(ROTLOC(3,2)**2+ROTLOC(1,3)**2+
cas     #         ROTLOC(2,1)**2)*WGT(KKK)

        DO I=1,3
        DO J=1,3
        ROTBAR_AS(I,J)=ROTBAR_AS(I,J)+ROTSLIP(I,J)*WGT(KKK)
        ENDDO
        ENDDO
cas
        KGX=KGX+1
1000  CONTINUE      ! END OF DO LOOP OVER GRAINS
2000  CONTINUE      ! END OF DO LOOP OVER PHASES

cas      KGX=1
cas      DO IPH=IPHBOT,IPHTOP
cas        IPHEL=IPH-IPHBOT+1
cas        DO KKK=NGR(IPH-1)+1,NGR(IPH)
cas
cas        DO I=1,3
cas        DO J=1,3
cas          AA(I,J)=AG(I,J,KKK)
cas          COROT(I,J)=0.
cas          do in=0,nneigh
cas            corot(i,j)=corot(i,j)+rot(i,j,neigh(in,kgx))*wneigh(in,kgx)
cas          enddo
cas        ENDDO
cas        ENDDO
cas
casC *** CALCULATE THE NEW TRASFORMATION MATRIX AND UPDATE
cas
cas        CALL REORIENT_GRAIN (AROT,COROT)
cas
cas        DO I=1,3
cas        DO J=1,3
cas          AG(I,J,KKK)=0.
cas          DO K=1,3
cas            AG(I,J,KKK)=AG(I,J,KKK)+AROT(I,K)*AA(K,J)
cas          ENDDO
cas        ENDDO
cas        ENDDO
cas
cas        KGX=KGX+1
cas        ENDDO      ! END OF DO LOOP OVER GRAINS
cas      ENDDO      ! END OF DO LOOP OVER PHASES


      RETURN
      END
c
c------------------------------------------------------------------------------
c
C **************************************************************************
C     SUBROUTINE UPDATE_SCHMID
C
C     ROTATES SCHMID TENSORS OF EACH GRAIN FROM CRYSTAL TO SAMPLE AXES
C **************************************************************************
      SUBROUTINE UPDATE_SCHMID

      INCLUDE 'vpsc_as.dim'

      DIMENSION aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3)
      DIMENSION aux33r(3,3)

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        DO IS=1,NSYST(IPHEL)

          DO J=1,5
            AUX5(J)=SCHCA(J,IS,IPHEL)
          ENDDO

          CALL CHG_BASIS(AUX5,AUX33,AUX55,AUX3333,1,5)

          DO I=1,3
          DO J=1,3
            AUX33R(I,J)=0.
            DO I1=1,3
            DO J1=1,3
              AUX33R(I,J)=AUX33R(I,J)+AG(I,I1,KKK)*AG(J,J1,KKK)*
     #                                AUX33(I1,J1)
            ENDDO
            ENDDO
          ENDDO
          ENDDO

          CALL CHG_BASIS(AUX5,AUX33R,AUX55,AUX3333,2,5)

          DO J=1,5
            SCH(J,IS,KGX)=AUX5(J)
          ENDDO
        ENDDO

      KGX=KGX+1
      ENDDO
      ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C **************************************************************************
C     SUBROUTINE UPDATE_SHAPE      --->      VERSION 28/NOV/2005
C
C     CALCULATES THE DIRECTION AND LENGTH (EIGENVECTORS AND EIGENVALUES OF
C     DEFORMATION TENSOR FIJ) OF THE AXES OF THE ELLIPSOID ASSOCIATED WITH
C     AVERAGE, PHASE AND/OR GRAIN ACCUMULATED DEFORMATION.
C     CALCULATES THE EULER ANGLES OF ELLIPSOID AXES WRT SAMPLE AXES.
C
C     REPLACES PREVIOUS SUBR. GRAXES_AVERAGE & SUBR. GRAXES_LOCAL  (SEP/04)
C **************************************************************************

cas      SUBROUTINE UPDATE_SHAPE (IPH)
      SUBROUTINE UPDATE_SHAPE_AS (IPH)

      INCLUDE 'vpsc_as.dim'

      DIMENSION W(3),BX(3,3),B(3,3),BT(3,3)


C *** IF IPH=0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN ELEMENT
C *** IF IPH>0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN PHASE 'IPH'
C *** CALCULATES EIGENVALUES, EIGENVECTORS & EULER ANGLES OF ELEMENT GRAIN,
C     PHASE GRAIN, OR INDIVIDUAL GRAIN
C *** 'AXISPH' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.

      DO I=1,3
      DO J=1,3
cas
cas        BX(I,J)=0.
cas        DO K=1,3
cas
cwx
cwx          BX(I,J)=BX(I,J)+FIJPH(I,K,IPH)*FIJPH(J,K,IPH)
cwx
cwx     uses UDOT forall phases to prevent
cwx     numerical instability in SO procedure
cwx
cas        if(interaction.eq.5) then
cas          BX(I,J)=BX(I,J)+FIJPH(I,K,0)*FIJPH(J,K,0)
cas          BX(I,J)=F2IJPH(I,K,0)*F2IJPH(J,K,0)
cas        else
cas          BX(I,J)=BX(I,J)+FIJPH(I,K,IPH)*FIJPH(J,K,IPH)
cass          BX(I,J)=F2IJPH(I,K,IPH)*F2IJPH(J,K,IPH)
          BX(I,J)=F2IJPH(I,J,IPH)
cas        endif
cwx
cas
cas        ENDDO
cas
      ENDDO
      ENDDO

      CALL JACOBI(BX,3,3,W,B,NROT,IER)
      CALL EIGSRT(W,B,3,3)
      IF (IER.EQ.1) THEN
        WRITE(*,*) 'ERROR IN UPDATE_SHAPE FOR PHASE ELLIPSOID',IPH
        STOP
      ENDIF

C *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
C *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR.
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING 1 AND 2.

      SIGN=-1.
      IF(DET(B).LE.0.) SIGN=1.
      DO I=1,3
        EXCHANGE=B(I,1)
        B(I,1)=B(I,2)
        B(I,2)=EXCHANGE*SIGN
      ENDDO
      EXCHANGE=W(1)
      W(1)=W(2)
      W(2)=EXCHANGE

      DO I=1,3
        AXISPH(0,I,IPH)=SQRT(W(I))
        DO J=1,3
          AXISPH(I,J,IPH)=B(I,J)
          BT(I,J)        =B(J,I)
        ENDDO
      ENDDO

      CALL EULER(1,ANG1,ANG2,ANG3,BT)
      EULERPH(1,IPH)=ANG1
      EULERPH(2,IPH)=ANG2
      EULERPH(3,IPH)=ANG3

C *** STOPS UPDATING THE ELLIPSOID FOR A LIMIT ASPECT RATIO

      IF(SQRT(W(2)/W(3)).GT.50. .AND.IFLAT(IPH).EQ.0.AND.IPH.NE.0) THEN
CFEB
cas        write(10,*) '***********************************************'
cas        write(10,*) 'NO FURTHER ELLIPSOID SHAPE UPDATING - PHASE ',iph
cas        write(10,*) '***********************************************'
        write(*, *) '***********************************************'
        write(*, *) 'NO FURTHER ELLIPSOID SHAPE UPDATING - PHASE ',iph
        write(*, *) '***********************************************'
CFEE
        IFLAT(IPH)=1
      ENDIF
CFEB
cas      write(10,'(''IPH='',I3,''  FIJPH  '',9F9.4)')
cas     #           IPH,((FIJPH(I,J,IPH),J=1,3),I=1,3)
cas      write(10,'(9X,''EIGNVAL'',3F9.4)')  (AXISPH(0,J,IPH),J=1,3)
cas      write(10,'(9X,''EIGNVEC'',9F9.4)') ((AXISPH(I,J,IPH),J=1,3),I=1,3)

C ********************************************************************
C     CALCULATES ELLIPSOID REPRESENTING DEFORMATION OF EACH GRAIN.
C     EIGENVALUES AND EIGENVECTORS OF EACH GRAIN 'KKK' IN PHASE 'IPH'.
C     'AXISGR(I,J)' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.

cas      IF(ISHAPE(IPH).GT.0) THEN
cas
cas      DO KKK=NGR(IPH-1)+1,NGR(IPH)
cas
cas      DO I=1,3
cas      DO J=1,3
cas        BX(I,J)=0.
cas        DO K=1,3
cas          BX(I,J)=BX(I,J)+FIJGR(I,K,KKK)*FIJGR(J,K,KKK)
cas        ENDDO
cas      ENDDO
cas      ENDDO
cas
cas      CALL JACOBI(BX,3,3,W,B,NROT,IER)
cas      CALL EIGSRT(W,B,3,3)
cas      IF (IER.EQ.1) THEN
cas        WRITE(*,'(''ERROR IN UPDATE_SHAPE FOR GRAIN'',I5,
cas     #            ''  IN PHASE'',I3)') KKK,IPH
cas        STOP
cas      ENDIF
cas
casC *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
casC *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
casC     CALCULATION OF THE ESHELBY TENSOR.
casC *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
casC     HANDED BY EXCHANGING 1 AND 2.
cas
cas      SIGN=-1.
cas      IF(DET(B).LE.0.) SIGN=1.
cas      DO I=1,3
cas        EXCHANGE=B(I,1)
cas        B(I,1)=B(I,2)
cas        B(I,2)=EXCHANGE*SIGN
cas      ENDDO
cas      EXCHANGE=W(1)
cas      W(1)=W(2)
cas      W(2)=EXCHANGE
cas
cas      DO I=1,3
cas        AXISGR(0,I,KKK)=SQRT(W(I))
cas        DO J=1,3
cas           AXISGR(I,J,KKK)=B(I,J)
cas        ENDDO
cas      ENDDO
cas
casC     IF(KKK.LT.5)
casC    # WRITE(10,'(''KKK='',I3,'' EIGVAL^2'',3F7.2,'' EIGVEC'',9F7.2)')
casC    #           KKK,W,((AXISGR(I,J,KKK),J=1,3),I=1,3)
cas
casC *** STOPS UPDATING THE ELLIPSOID WHEN THE 'WORST OFFENDER' GRAIN REACHES
casC     A LIMIT ASPECT RATIO
cas
cas      IF(SQRT(W(2)/W(3)).GT.50. .AND.IFLAT(IPH).EQ.0.AND.IPH.NE.0) THEN
cas        write(10,*) '***********************************************'
cas        write(10,*) 'NO FURTHER ELLIPSOID SHAPE UPDATING - PHASE ',iph
cas        write(10,*) '***********************************************'
cas        write(*, *) '***********************************************'
cas        write(*, *) 'NO FURTHER ELLIPSOID SHAPE UPDATING - PHASE ',iph
cas        write(*, *) '***********************************************'
cas        IFLAT(IPH)=1
cas      ENDIF
cas
cas      ENDDO      ! END OF DO KKK
cas      ENDIF      ! END OF IF ISHAPE>0
casC ********************************************************************
CFEE
      RETURN
      END
c
c------------------------------------------------------------------------------
c
cass      SUBROUTINE UPDATE_TWINNING_AS (twVFRates)
      SUBROUTINE UPDATE_TWINNING_AS (twVFRates,ktwsmax,kgrtot)
c
      INCLUDE 'vpsc_as.dim'
c
      dimension twVFRates(ktwsmax,kgrtot)
c
      do I=1,ktwsmax
         do J=1,kgrtot
            twVFRates(I,J)=0.
         enddo
      enddo

      do iph=1,nph
c
      if(ntwmod(iph).gt.0) then
c
      KGX=1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
c
        ITS=NSYST(IPH)-NTWSYS(IPH)
        KTS=0
        DO ITM=1,NTWMOD(IPH)
          KK=ITM+NMODES(IPH)-NTWMOD(IPH)
          TWSHX=TWSH(ITM,IPH)
          DO NSMX=1,NSM(KK,IPH)
            KTS=KTS+1             ! shifted counter for twin systems
            ITS=ITS+1             ! absolute counter over all systems
            IF(GAMDOT(ITS,KGX).GT.0.) THEN
            TWVFRATES(KTS,KKK)=ABS(GAMDOT(ITS,KGX))/TWSHX
cas
cas      To prevent twinning gamdots from contributing to 
cas      reorientation rate of the grain one has to execute:
cas
cas      GAMDOT(ITS,KGX)=0
cas
            ENDIF
          ENDDO
        ENDDO      ! END OF LOOP OVER TWIN MODES IN THE PHASE
c
        KGX=KGX+1
      ENDDO      ! END OF LOOP OVER GRAINS IN THE PHASE
c
      endif
c
      enddo   ! nph endif

      RETURN
      END
c
c------------------------------------------------------------------------------
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINE VOIGT   ---->   VERSION OF 09/02/98
C
C     TRANSFORMS 6X1 MATRIX T1 INTO SECOND ORDER TENSOR T2 IF IOPT=1
C     AND VICEVERSA IF IOPT=2.
C     TRANSFORMS 6X6 MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF IOPT=3
C     AND VICEVERSA IF IOPT=4.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE VOIGT(T1,T2,C2,C4,IOPT)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
      DIMENSION IJV(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

      IF(IOPT.EQ.1) THEN
      DO 30 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      T2(I1,I2)=T1(I)
   30 T2(I2,I1)=T1(I)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 40 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
   40 T1(I)=T2(I1,I2)
      ENDIF
C
      IF (IOPT.EQ.3) THEN
      DO 10 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 10 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
   10 C4(I2,I1,J2,J1)=C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.4) THEN
      DO 20 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 20 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
   20 C2(I,J)=C4(I1,I2,J1,J2)
      ENDIF
C
      RETURN
      END
c
c------------------------------------------------------------------------------
c
C********************************************************************
C     SUBROUTINE VPSC     --->      VERSION 29/nov/05
C********************************************************************
cas
cas      SUBROUTINE VPSC (ISTEP)
cas
      SUBROUTINE VPSC

      INCLUDE 'vpsc_as.dim'

      DIMENSION EGA(3,3,3,3),RGA(3,3,3,3),
     #          C4SA(3,3,3,3),ESA(3,3,3,3),RSA(3,3,3,3),
     #          EINVSA(3,3,3,3),PGA(3,3)
cw      DIMENSION PGA(3,3),,PSA(3,3),P5(5)
      DIMENSION E5(5,5),E5INV(5,5)
      DIMENSION XMTNEW(5,5)
      DIMENSION XIMSINV(5,5),FS(5,5),XMAST(5,5)
      DIMENSION XMCBCTG(5,5),XMTBAVE(5,5)
      DIMENSION BC1(5,5),BC2(5,5),BCINV(5,5),BCAVE(5,5)
      DIMENSION AUX5(5),AUX33(3,3),AUX55(5,5),AUX3333(3,3,3,3)
      DIMENSION DBAUX(5)
      DIMENSION EIGB(3,3),AXB(3)
c
      DIMENSION IDLEQ(5),ISLEQ(5)
cw      DIMENSION DAUX(6),SAUX(6),AUX66(6,6)
c
      dimension xmcphc(5),xmtphave(5),dsimaux(3,3)
      dimension phave(5),dzero1(5),dznew(5)

      integer diagnostics
      integer flag

      diagnostics = 0
C *** FOR A TAYLOR CASE (INTERACTION=0) SOLVES VP EQUATION FOR GRAIN STRESS.
C *** CALCULATE STRAIN-RATE 'DG' AND MODULI 'XMCTG' FOR EVERY GRAIN.

      IF(INTERACTION.EQ.0) THEN
        KGX=1
        DO IPH=IPHBOT,IPHTOP
          IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            DO I=1,5
              STRY(I,KGX)=SG(I,KKK)
            ENDDO
              CALL GRAIN_STRESS (INTERACTION,KGX,KKK,IPHEL,IPH)
              CALL GRAIN_RATE_AND_MODULI (0,0,KGX,KKK,IPHEL,IPH)
           KGX=KGX+1
          ENDDO
        ENDDO

C *** CALCULATE AVERAGE STRESS AND STRAIN-RATE

        DO I=1,5
          SAV(I)=0.
          DAV(I)=0.
          KGX=1
          DO IPH=IPHBOT,IPHTOP
            DO KKK=NGR(IPH-1)+1,NGR(IPH)
              SAV(I)=SAV(I)+SG(I,KKK)*WGT(KGX)
              DAV(I)=DAV(I)+DG(I,KGX)*WGT(KGX)
              KGX=KGX+1
            ENDDO
          ENDDO
          SBAR(I)=SAV(I)      ! ONLY REQUIRED TO CALCULATE DVM & SVM
        ENDDO
        CALL CHG_BASIS(SBAR,SDEVIAT,AUX55,AUX3333,1,5)

        RETURN
      ENDIF

C *** CASE OF SELF-CONSISTENT CALCULATION (INTERACTION>0)
c
      if(irsvar.eq.0) then
        jxrsini=nrs(1,1)
        jxrsfin=nrs(1,1)
        jxrstep=nrs(1,1)
      endif
cwx
cwx   comment next line to perform
cwx   rs loop at every step
cwx
cas
cas      if(istep.gt.1) jxrsini=jxrsfin
cas
c
      IRS=0
      DO JXRS=JXRSINI,JXRSFIN,JXRSTEP
         IRS=IRS+1
CFEB
C *** QUADRATIC EXTRAPOLATION GUESS FOR ASO and ESO
      IF(INTERACTION.EQ.5 .AND. IRS.GE.4) THEN
         write (*,*) "Second order not implemented"
         STOP
      ENDIF
CFEE
      if(jxrsini.ne.jxrsfin) then
        write(*,*)
        write(*,*) 'NRS ITERATION', IRS
        write(*,*)
      endif
c
      ITSO=0
      ERRESO=2.*ERRSO
      ERRASO=2.*ERRSO
      KSO=1
c
      if (diagnostics.eq.1) then
         write(*,*)
         IF(INTERACTION.EQ.1) write(*,'(a)') ' AFFINE CALCULATION'
         IF(INTERACTION.EQ.2) write(*,'(a)') ' SECANT CALCULATION'
         IF(INTERACTION.EQ.3) write(*,'(a)') ' NEFF=10 CALCULATION'
         IF(INTERACTION.EQ.4) write(*,'(a)') ' TANGENT CALCULATION'
         write(*,*)
      endif

      DO WHILE( KSO.EQ.1 .AND.
     #  (ERRESO.GT.ERRSO.OR.ERRASO.GT.ERRSO) .AND. ITSO.LT.ITMAXSO )

      RELSGR=2*ERRS
      RELS=2*ERRS
      RELD=2*ERRD

C *** OUTER LOOP: VARIES STRESS AND COMPLIANCE IN THE GRAINS

      IT2=0
c
      DO WHILE( (RELSGR.GT.ERRS.OR.RELD.GT.ERRD.OR.RELS.GT.ERRS)
     # .AND. IT2.LT.ITMAXEXT)
c
      IT2=IT2+1

C *** INNER LOOP: VARIES OVERALL TANGENT COMPLIANCE

      IT1TG=0
      RER1TG=2*ERRM
c
      DO WHILE(RER1TG.GT.ERRM.AND.IT1TG.LT.ITMAXINT)
c
cpp      write(*,*) 'XLTG='
cpp      write(*,*) xltg
cpp      write(*,*) 'DZERO='
cpp      write(*,*) dzero
cpp      pause
c
        IT1TG=IT1TG+1
c
        DO I=1,5
          XMTPHAVE(I)=0.
          PHAVE(I)=0.
          DO J=1,5
            BCAVE(I,J)=0.
            XMTBAVE(I,J)=0.
          ENDDO
        ENDDO

C *************************************************************************
C     CALCULATE ESHELBY TENSOR S=S(Mtg) AND FS=[1/(I-S)]S.
C     S(Mtg) WILL BE THE SAME FOR EVERY GRAIN IN A GIVEN PHASE WHEN
C     ISHAPE.LE.1.
C     IT WILL BE DIFFERENT FOR EVERY GRAIN WHEN ISHAPE.GE.2 BECAUSE Mtg OR
C     Msec HAVE TO BE ROTATED TO ELLIPSOID AXES AND, IN ADDITION, THE
C     SIZE OF THE ELLIPSOID AXES IS DIFFERENT FOR EVERY GRAIN.
C     WHEN Mtg=n*MseC THEN S(Mtg)=S(Msec), INDEPENDENT OF ISHAPE VALUE.
C *************************************************************************

C *** SKIP EVERY OTHER CALCULATION OF MSTAR
      ISKIP=1
      IF(ISKIP.EQ.1) THEN
c
      CALL CHG_BASIS(AUX5,AUX33,XLTG,C4SA,3,5)

C **************************************************************************
C *** LOOP #1 OVER PHASES AND GRAINS
C **************************************************************************

      KGX=1
      DO IPH=IPHBOT,IPHTOP

        IF(ISHAPE(IPH).LE.1) THEN
          DO I=1,3
            AXB(I)=AXISPH(0,I,IPH)
            DO J=1,3
              EIGB(I,J)=AXISPH(I,J,IPH)
            ENDDO
          ENDDO
        ENDIF

      DO KKK=NGR(IPH-1)+1,NGR(IPH)
CFEB
cas        if(ishape(iph).ge.2) then      ! INDIVIDUAL GRAIN SHAPE
cas          do i=1,3
cas            axb(i)=AXISGR(0,I,KKK)
cas            do j=1,3
cas              eigb(i,j)=AXISGR(i,j,kkk)
cas            enddo
cas          enddo
cas        endif
CFEE
C *** ESHELBY CALCULATION FOR EVERY PHASE OR FOR EVERY GRAIN

      IF(ISHAPE(IPH).GE.2 .OR. KKK.EQ.NGR(IPH-1)+1) THEN

C     ROTATION OF STIFFNESS 'C4SA' TO ELLIPSOID PRINCIPAL AXES
      DO 95 I=1,3
      DO 95 J=1,3
      DO 95 M=1,3
      DO 95 N=1,3
        DUMMY=0.
        DO 90 I1=1,3
        DO 90 J1=1,3
        DO 90 M1=1,3
        DO 90 N1=1,3
          DUMMY=DUMMY+EIGB(I1,I)*EIGB(J1,J)*EIGB(M1,M)
     #           *EIGB(N1,N)*C4SA(I1,J1,M1,N1)
   90   CONTINUE
        C4GA(I,J,M,N)=DUMMY
   95 CONTINUE

      IF(ICAUCHY.EQ.0) IOPTION=2
      IF(ICAUCHY.EQ.1) IOPTION=3
      CALL ESHELBY(AXB,C4GA,0.0+0,EGA,RGA,AUX33,PGA,PDIL,
     #             AUX3333,AUX3333,IOPTION)

c     write(10,'(''  distortion Eshelby tensor'')')
c     write(10,'(9f9.5)') ((((ega(i,j,k,l),l=1,3),k=1,3),j=1,3),i=1,3)
c     write(10,'(''  rotation Eshelby tensor'')')
c     write(10,'(9f9.5)') ((((rga(i,j,k,l),l=1,3),k=1,3),j=1,3),i=1,3)
c     write(10,'(''  hydrostatic Eshelby tensor'')')
c     write(10,'(3f9.5)') ((pga(i,j),j=1,3),i=1,3)
c     write(10,'(''  hydrostatic Eshelby factor'')')
c     write(10,'(1f9.5)') pdil

C     ROTATES THE DISTORTION, ROTATION AND PRESSURE ESHELBY TENSORS
C     FOR THE PHASE OR FOR EACH GRAIN BACK INTO SAMPLE AXES.

cw      DO 110 I=1,3
cw      DO 110 J=1,3
cw        DUMMYP=0.
cw        DO 100 I1=1,3
cw        DO 100 J1=1,3
cw          DUMMYP=DUMMYP+EIGB(I,I1)*EIGB(J,J1)*PGA(I1,J1)
cw  100   CONTINUE
cw        PSA(I,J)=DUMMYP
cw  110 CONTINUE

      DO 130 I=1,3
      DO 130 J=1,3
      DO 130 M=1,3
      DO 130 N=1,3
        DUMMYE=0.
        DUMMYR=0.
        DO 120 I1=1,3
        DO 120 J1=1,3
        DO 120 M1=1,3
        DO 120 N1=1,3
          DUMMYE=DUMMYE+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     #           *EIGB(N,N1)*EGA(I1,J1,M1,N1)
          DUMMYR=DUMMYR+EIGB(I,I1)*EIGB(J,J1)*EIGB(M,M1)
     #           *EIGB(N,N1)*RGA(I1,J1,M1,N1)
  120   CONTINUE
        ESA(I,J,M,N)=DUMMYE
        RSA(I,J,M,N)=DUMMYR
  130 CONTINUE

      CALL CHG_BASIS(AUX5,AUX33,E5,ESA,4,5)
cw      CALL CHG_BASIS(P5,PSA,AUX55,AUX3333,2,5)

      DO I=1,5
      DO J=1,5
        E5INV(I,J)=E5(I,J)
        XIMSINV(I,J) =XID5(I,J)-E5(I,J)
      ENDDO
      ENDDO

      flag = 0
      CALL LU_INVERSE(XIMSINV,5,flag,diagnostics)
        if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'XIMSINV is singular'
               DO I=1,5
                  write(*,*) (ximsinv(I,J),j=1,5)
               ENDDO
            endif
        endif

      flag = 0
      CALL LU_INVERSE(E5INV,5,flag,diagnostics)
        if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'E5INV is singular'
               DO I=1,5
                  write(*,*) (e5inv(I,J),j=1,5)
               ENDDO
            endif
        endif


      CALL CHG_BASIS(AUX5,AUX33,E5INV,EINVSA,3,5)

      DO I=1,5
      DO J=1,5
        FS(I,J)=0.
        DO K1=1,5
          FS(I,J)=FS(I,J)+XIMSINV(I,K1)*E5(K1,J)
        ENDDO
      ENDDO
      ENDDO

      DO I=1,5
      DO J=1,5
      XIMSINVPH(I,J,IPH)=XIMSINV(I,J)
      FSPH(I,J,IPH)=FS(I,J)
      ENDDO
      ENDDO
c
c
c     CALCULATE M* = inv(I-S) * S * Mtg
c
      DO I=1,5
      DO J=1,5
cw
        DUMMY=0.
        DO K=1,5
        DUMMY=DUMMY+FS(I,K)*XMTG(K,J)
        ENDDO
c
        if(interaction.eq.3) then      ! neff=10 case
c
        XMAST(I,J)=10*DUMMY
c
        else if(interaction.eq.4) then      ! tangent case
cw
cw        XMAST(I,J)=NRS(1,1)*DUMMY
cw
          if(irsvar.eq.1) then
            XMAST(I,J)=jxrs*DUMMY
          else
            XMAST(I,J)=NRS(1,1)*DUMMY
          endif
c
        else
c
          XMAST(I,J)=DUMMY
c
        endif
cw
      ENDDO
      ENDDO
C
C *** COPIES TENSORS INTO PHASE ARRAYS OR INTO GRAIN ARRAYS.
C *** 'ASPH' OR 'ASGR' ARE USED IN SUBR. UPDATE_ORIENTATION TO CALCULATE
C     SPIN-RATE DEVIATIONS.
C *** 'D5PH' OR 'D5GR' ARE USED IN SUBR. CALC_CAUCHY TO CALCULATE
C     PRESSURE DEVIATIONS.

      IF(ISHAPE(IPH).LE.1) THEN
        DO 140 I=1,3
        DO 140 J=1,3
        DO 140 K=1,3
        DO 140 L=1,3
          DUMMY=0.
          DO K1=1,3
          DO L1=1,3
            DUMMY=DUMMY+RSA(I,J,K1,L1)*EINVSA(K1,L1,K,L)
          ENDDO
          ENDDO
          ASPH(I,J,K,L,IPH)=DUMMY
  140   CONTINUE
      ENDIF
CFEB
cas      IF(ISHAPE(IPH).GE.2) THEN
cas        DO 150 I=1,3
cas        DO 150 J=1,3
cas        DO 150 K=1,3
cas        DO 150 L=1,3
cas          DUMMY=0.
cas          DO K1=1,3
cas          DO L1=1,3
cas            DUMMY=DUMMY+RSA(I,J,K1,L1)*EINVSA(K1,L1,K,L)
cas          ENDDO
cas          ENDDO
cas          ASGR(I,J,K,L,KKK)=DUMMY
cas  150   CONTINUE
cas      ENDIF
cas
cw      IF(ICAUCHY.EQ.1) THEN
cw        DO I=1,5
cw          DUMMY=0.
cw          DO J=1,5
cw            DUMMY=DUMMY+P5(J)*E5INV(J,I)
cw          ENDDO
cw          IF(ISHAPE(IPH).LE.1) D5PH(I,IPH)=DUMMY
cw          IF(ISHAPE(IPH).GE.2) D5GR(I,KKK)=DUMMY
cw        ENDDO
cw      ENDIF
CFEE
      IF(ISHAPE(IPH).LE.1) THEN
        DO I=1,5
        DO J=1,5
          XMASTPH(I,J,IPH)=XMAST(I,J)
        ENDDO
        ENDDO
      ENDIF
CFEB
cas      IF(ISHAPE(IPH).GE.2) THEN
cas        DO I=1,5
cas        DO J=1,5
cas          XMASTGR(I,J,KKK)=XMAST(I,J)
cas        ENDDO
cas        ENDDO
cas    ENDIF
CFEE
      ENDIF           !  END OF IF(ISHAPE(IPH.GE.2 .OR. KKK.EQ.NGR(IPH-1)+1)

      KGX=KGX+1
      ENDDO           !  END OF DO OVER GRAINS   (LOOP #1)
      ENDDO           !  END OF DO OVER PHASES   (LOOP #1)

      ENDIF           !  END OF ESHELBY CALCULATION BIG IF

C **************************************************************************
C *** LOOP #2 OVER PHASES AND GRAINS
C **************************************************************************

C *** SOLVES ITERATIVELY SC EQUATION FOR GENERAL CASE OF DIFFERENT SHAPES.
C *** REQUIRED FOR MULTI-PHASE WHEN GRAINS OF EACH PHASE HAVE DIFFERENT SHAPE,
C     OR FOR ISHAPE>1, WHEN EVERY GRAIN HAS DIFFERENT SHAPE, OR BOTH.

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        IF(ISHAPE(IPH).LE.1) THEN    ! shouldn't this be done outside DO KKK ?
          DO I=1,5
          DO J=1,5
            XMAST(I,J)=XMASTPH(I,J,IPH)
            BC2(I,J)  =XMAST(I,J)+XMTG(I,J)
          ENDDO
          ENDDO
        ENDIF
CFEB
cas        IF(ISHAPE(IPH).GE.2) THEN
cas          DO I=1,5
cas          DO J=1,5
cas            XMAST(I,J)=XMASTGR(I,J,KKK)
cas            BC2(I,J)  =XMAST(I,J)+XMTG(I,J)
cas          ENDDO
cas          ENDDO
cas        ENDIF
CFEE
        DO I=1,5
        DO J=1,5
          BC1(I,J)=XMAST(I,J)+XMCTG(I,J,KGX)
        ENDDO
        ENDDO

      flag = 0
        CALL LU_INVERSE(BC1,5,flag,diagnostics)
        if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'BC1 is singular'
               DO I=1,5
                  write(*,*) (bc1(I,J),j=1,5)
               ENDDO
            endif
        endif

cc
        DO I=1,5
        DO J=1,5
c
          CHIFLU(I,J,KGX)=BC1(I,J)
c
          DUMMY=0.
          DUM2=0.
          DO K1=1,5
            DUMMY=DUMMY+BC1(I,K1)*BC2(K1,J)
            DUM2=DUM2+XMCTG(I,K1,KGX)*BC1(K1,J)
          ENDDO
          BC(I,J,KGX)=DUMMY
          BETFLU(I,J,KGX)=DUM2
        ENDDO
        ENDDO
cc
        DO I=1,5
          DUMMY=0.
          DO J=1,5
            DUMMY=DUMMY-
     #            BC1(I,J)*(DCZERO(J,KGX)-DZERO(J))
          ENDDO
          PHIC(I,KGX)=DUMMY
        ENDDO
cc
        DO I=1,5
        DO J=1,5
          XMCBCTG(I,J) =0.
          DO K1=1,5
            XMCBCTG(I,J)=XMCBCTG(I,J)+XMCTG(I,K1,KGX)*BC(K1,J,KGX)
          ENDDO
        ENDDO
        ENDDO
cc
        DO I=1,5
        XMCPHC(I)=DCZERO(I,KGX)
        DO J=1,5
        XMCPHC(I)=XMCPHC(I)+XMCTG(I,J,KGX)*PHIC(J,KGX)
        ENDDO
        ENDDO
cc
        DO I=1,5
cc
          XMTPHAVE(I)=XMTPHAVE(I)+XMCPHC(I)*WGT(KKK)
          IF(IBCINV.EQ.1) THEN
            PHAVE(I)=PHAVE(I)+PHIC(I,KGX)*WGT(KKK)
          ELSE
            PHAVE(I)=0.
          ENDIF
cc
        DO J=1,5
          BCAVE(I,J)  =BCAVE(I,J)  +BC(I,J,KGX)  *WGT(KGX)
          XMTBAVE(I,J)=XMTBAVE(I,J)+XMCBCTG(I,J) *WGT(KGX)
        ENDDO
        ENDDO
        KGX=KGX+1

      ENDDO
      ENDDO

C ****************************************************************
C *** END OF LOOP #2 OVER GRAINS AND PHASES
C ****************************************************************

      DO I=1,5
      DO J=1,5
c
        IF(IBCINV.EQ.1) THEN
          BCINV(I,J)=BCAVE(I,J)
        ELSE
          BCINV(I,J)=XID5(I,J)
         ENDIF
c
      ENDDO
      ENDDO

      flag = 0
      CALL LU_INVERSE(BCINV,5,flag,diagnostics)
        if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'BCINV is singular'
               DO I=1,5
                  write(*,*) (bcinv(I,J),j=1,5)
               ENDDO
            endif
        endif


      DO I=1,5
      DO J=1,5
        XMTNEW(I,J)=0.
        DO K1=1,5
          XMTNEW(I,J)=XMTNEW(I,J)+XMTBAVE(I,K1)*BCINV(K1,J)
        ENDDO
      ENDDO
      ENDDO


C *************************************************************************
C     Mold-Mnew (tangent) COMPARISON

        RER1TG=TMISMATCH5(XMTG,XMTNEW,5,5)
CFEB
c       WRITE(*,'(1H+,5X,''IT1TAN='',I3,'' --> TAN MOD REL ERROR'',
c    #          2E12.3)') IT1TG,RER1TG
c       WRITE(*,'(1H+,I5,3E11.3,I7,E11.3,I7,E11.3)')
c    #        IT2,RELSGR,RELS,RELD, IT1TG,RER1TG
CFEE
      DO I=1,5
      DO J=1,5
        XMTG(I,J)=XMTNEW(I,J)
        XLTG(I,J)=XMTNEW(I,J)
      ENDDO
      ENDDO

      flag = 0
      CALL LU_INVERSE(XLTG,5,flag,diagnostics)
        if (flag.eq.1)  then
            if (diagnostics.eq.1) then 
               write(*,*) 'XLTG is singular'
               DO I=1,5
                  write(*,*) (xltg(I,J),j=1,5)
               ENDDO
            endif
        endif


      if(interaction.eq.1.or.interaction.eq.5) then       ! affine or SO
        DO I=1,5
          DZERO1(I)=0.
          DO K1=1,5
            DZERO1(I)=DZERO1(I)+XMTNEW(I,K1)*PHAVE(K1)
          ENDDO
        ENDDO
        DO I=1,5
          DZNEW(I)=XMTPHAVE(I)-DZERO1(I)
        ENDDO
        RERDZERO=TMISMATCH(DZERO,DZNEW,5,1)
        DO I=1,5
           DZERO(I)=DZNEW(I)
        ENDDO
cc      write(*,*)
cc      write(*,*) 'DZERO=',dzero
      else
        DO I=1,5
          DZERO(I)=0.
        ENDDO
        RERDZERO=0.
      endif
      RERDZERO=RERDZERO    ! to fool the compiler
c
cw       WRITE(*,'(1H+,5X,''IT1TAN='',I3,'' --> TG,D0 REL ERRORS'',
cw     #          2E12.3)') IT1TG,RER1TG,RERDZERO
c
      ENDDO         !  END OF (DO..WHILE) FOR INNER ITERATION (IT1TG)

C *************************************************************************
C     BOUNDARY CONDITIONS:
C     * IF DIAGONAL IMPOSED CALLS SUBROUTINE 'STATE_5x5' (RL: 26/1/00)
C       AND SOLVES STRESS-STRAIN COMPONENTS IN DEVIATORIC SPACE
C     * ELSE, CALLS SUBROUT 'STATE_6x6' WITH DSIM(3,3) & SCAUCHY(3,3)
C       PLUS INDICES OF THE KNOWN COMPONENTS. SOLVES FOR UNKNOWN COMPONENTS
C       IN 6-DIM CAUCHY SPACE.
c
      DO I=1,5
        DBAUX(I)=DBAR(I)-DZERO(I)
      ENDDO
c
      IF(IDSIM(1)*IDSIM(2)*IDSIM(3).EQ.1) THEN

        idleq(1)=1
        idleq(2)=1
        idleq(3)=idsim(4)
        idleq(4)=idsim(5)
        idleq(5)=idsim(6)

        isleq(1)=0
        isleq(2)=0
        isleq(3)=iscau(4)
        isleq(4)=iscau(5)
        isleq(5)=iscau(6)

        CALL STATE_5x5 (IDLEQ,DBAUX,ISLEQ,SBAR,XMTG)
c
        DO I=1,5
          DBAR(I)=DBAUX(I)+DZERO(I)
        ENDDO
c
        CALL CHG_BASIS (DBAR,DSIM   ,AUX55,AUX3333,1,5)
        CALL CHG_BASIS (SBAR,SDEVIAT,AUX55,AUX3333,1,5)
c
cw      write(*,*) dzero
cw      write(*,*) dbar
cw      write(*,*) sbar
cw      pause

      ELSE

        CALL CHG_BASIS (DBAUX,DSIMAUX,AUX55,AUX3333,1,5)
        CALL STATE_6x6 (IDSIM,DSIMAUX,ISCAU,SCAUCHY,XMTG, flag)
         if (flag.eq.1) then
            write(*,*) 'idsim = '
            write(*,*) (idsim(J),j=1,6)
            
            write(*,*) 'dsimaux = '
            DO I=1,3
               write(*,*) (dsimaux(I,J),j=1,3)
            ENDDO
            
            write(*,*) 'iscau = '
            write(*,*) (iscau(J),j=1,6)
            
            write(*,*) 'scauchy = '
            DO I=1,3
               write(*,*) (scauchy(I,J),j=1,3)
            ENDDO

            write(*,*) 'xmtg = '
            DO I=1,5
               write(*,*) (xmtg(I,J),j=1,5)
            ENDDO

         endif
         
c
        PREMAC=1./3.*(SCAUCHY(1,1)+SCAUCHY(2,2)+SCAUCHY(3,3))
c
        CALL CHG_BASIS (DBAUX,DSIMAUX,AUX55,AUX3333,2,5)
c
        DO I=1,5
          DBAR(I)=DBAUX(I)+DZERO(I)
        ENDDO
c
        CALL CHG_BASIS (DBAR,DSIM,AUX55,AUX3333,1,5)
c
        CALL CHG_BASIS (SBAR,SCAUCHY,AUX55,AUX3333,2,5)
        CALL CHG_BASIS (SBAR,SDEVIAT,AUX55,AUX3333,1,5)
c
      ENDIF
c
        SVM=SQRT(3./2.)*TNORM(SBAR,5,1)

C *************************************************************************
C     S* & D* REQUIRED FOR DIFFERENT GRAIN SHAPES
C     SASTAV(I) IS USED INSIDE GRAIN_STRESS
C     SASTBAR(I) IS USED BELOW TO CALCULATE DAST(I)

      SASTMIX=0.
      if(interaction.ne.3.and.interaction.ne.4) SASTMIX=0.75
      DO I=1,5
        SASTAV(I) =0.
        SASTBAR(I)=0.
        DO J=1,5
          SASTAV(I) =SASTAV(I) +BCINV(I,J)*(SAV(J) -PHAVE(J))
          SASTBAR(I)=SASTBAR(I)+BCINV(I,J)*(SBAR(J)-PHAVE(J))
        ENDDO
        SAVEX=SASTAV(I)
        SBARX=SASTBAR(I)
        SASTAV(I) =(1.-SASTMIX)*SAVEX+SASTMIX*SBARX
        SASTBAR(I)=SASTMIX*SAVEX+(1.-SASTMIX)*SBARX
      ENDDO

      DO I=1,5
        DAST(I)=0.
        DO J=1,5
         DAST(I)=DAST(I)+XMTG(I,J)*SASTBAR(J)
        ENDDO
      ENDDO

C *** SAVE CURRENT STRESS IN GRAIN.
C *** STARTING FROM CONVERGED 'XMTG' AND 'STRY' RECALCULATES 'SG' USING THE
C     INTERACTION EQUATION.
C *** USES NEW 'SG' TO CALCULATE STRAIN-RATE 'DG', SHEAR RATES 'GAMDOT' AND
C     MODULI 'XMCTG' FOR EVERY GRAIN.

      KGX=1
      DO IPH=IPHBOT,IPHTOP
        IPHEL=IPH-IPHBOT+1
      DO KKK=NGR(IPH-1)+1,NGR(IPH)
        DO I=1,5
          STRY(I,KGX)=SG(I,KKK)
        ENDDO
CFEB
cas        if((irs.eq.1.and.itso.le.2.and.istep.eq.1)
cas     #                       .or.interaction.ne.5) then
cas
        if(irs.eq.1.and.itso.le.2.or.interaction.ne.5) then
cas
cwx        if((irs.eq.1.and.itso.eq.1.and.istep.eq.1)
cwx     #                       .or.interaction.ne.5) then
CFEE
          CALL GRAIN_STRESS (1,KGX,KKK,IPHEL,IPH)
          CALL GRAIN_RATE_AND_MODULI (1,1,KGX,KKK,IPHEL,IPH)
CFEB
        else
           write(*,*) "second order is not implemented"
           stop
        endif
CFEE
c
      KGX=KGX+1
      ENDDO
      ENDDO

C *** UPDATE AVERAGE STRESS AND AVERAGE STRAIN-RATE.
C *** CALCULATES STRESS CONVERGENCE IN GRAINS.

      SDEV=0.
      DO I=1,5
        SAV(I)=0.
        DAV(I)=0.
        KGX=1
        DO IPH=IPHBOT,IPHTOP
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            SAV(I)=SAV(I)+SG(I,KKK)*WGT(KGX)
            DAV(I)=DAV(I)+DG(I,KGX)*WGT(KGX)
            SDEV=SDEV+(SG(I,KKK)-STRY(I,KGX))**2*WGT(KGX)
            KGX=KGX+1
          ENDDO
        ENDDO
      ENDDO

C *** CHECK CONSISTENCY OF LOCAL, MACROSCOPIC AND AVERAGE MAGNITUDES.
C      a) <SG-STRY> .LT.ERR
C      b) /DBAR-DAV/.LT.ERR
C      c) /SBAR-SAV/.LT.ERR

      RELSGR=SQRT(SDEV)/TNORM(SBAR,5,1)
      RELD  =TMISMATCH(DBAR,DAV,5,1)
      RELS  =TMISMATCH(SBAR,SAV,5,1)

C *** DEFINES EMPIRICAL MIX COEF FOR SASTAV & SASTBAR BASED ON CONVERGENCE
C     CONVNEW=RELSGR+RELD+RELS
C     CONVDIF=CONVNEW-CONVOLD
C     SASTMIX=0.
C     IF(IT2.GT.10 .AND. CONVDIF.GT.0.) SASTMIX=0.25
C     IF(IT2.GT.20 .AND. CONVDIF.GT.0.) SASTMIX=0.50
C     IF(IT2.GT.30 .AND. CONVDIF.GT.0.) SASTMIX=0.75
C     CONVOLD=CONVNEW
CFEB
!      WRITE(*,'(1H+,I5,3E11.3,I7,E11.3,I7,E11.3)')
!     #      IT2,RELSGR,RELS,RELD
ccas      WRITE(UW1,'(  I5,3E12.3,I7,E12.3,I7,E12.3)')
ccas     #      IT2,RELSGR,RELS,RELD
CFEE

      ENDDO     ! END OF (DO..WHILE) FOR OUTER ITERATION (IT2) ON 'SG'
c
      SVM=0.
      DVM=0.
      DO I=1,5
        SVM=SVM+SBAR(I)*SBAR(I)
        DVM=DVM+SBAR(I)*DBAR(I)
      ENDDO
      SVM=SQRT(SVM*3./2.)
      DVM=DVM/SVM
CFEB
      if(iflu.eq.1) then
         write(*,*) "second order is not implemented"
         stop
      endif
c
      if(interaction.eq.5.and.itso.ge.2) then
c     as      write(97,'(2i5,2x,3e12.3)') irs,itso,erraso,erreso
         write(*,*) "second order is not implemented"
         stop
      endif
CFEE
      ITSO = ITSO + 1
      ENDDO     ! END OF (DO..WHILE) FOR SO ITERATION (ITSO)
CFEB
cas      if(iflu.eq.1) then
cas       write(83,'(a,i7)') ' NRS = ',jxrs
cas       if(icubcom.eq.1) write(84,'(a,i7)') ' NRS = ',jxrs
cas       call sdpx
cas       if(interaction.eq.5) then
cas        write(97,'(a,i7,4e11.4)')
cas     #   ' NRS,SVM,SDS,SDD,UTILDE = '
cas     #   ,jxrs,svm,sdseqintra,sddeqintra,utilde
cas        write(97,*)
cas       endif
cas      endif

C *** UPDATES VALUES FOR PERFORMING EXTRAPOLATION ON ESO and ASO
      IF(INTERACTION.EQ.5) THEN
         write(*,*) "Second order not implemented"
         STOP
      ENDIF
CFEE

      ENDDO     ! END OF DO JXRS FOR NRS ITERATION
CFEB
      IF(INTERACTION.EQ.5) then
         write(*,*) "Second order not implemented"
         stop
      endif
CFEE
      RETURN
      END

c------------------------------------------------------------------------------
c
C ************************************************************
C     SUBROUTINE VOIGT10     --->     VERSION 4/JAN/07
C ************************************************************

      SUBROUTINE VOIGT10(T1,T2,IOPT)

      IMPLICIT INTEGER (I-N), REAL*8 (A-H, O-Z)

      DIMENSION T1(10),T2(4,4)
      DIMENSION IJV(10,2)
      DATA ((IJV(N,M),M=1,2),N=1,10)
     #  /1,1,2,2,3,3,2,3,1,3,1,2,1,4,2,4,3,4,4,4/

      IF(IOPT.EQ.1) THEN
      DO 30 I=1,10
      I1=IJV(I,1)
      I2=IJV(I,2)
      T2(I1,I2)=T1(I)
   30 T2(I2,I1)=T1(I)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 40 I=1,10
      I1=IJV(I,1)
      I2=IJV(I,2)
   40 T1(I)=T2(I1,I2)
      ENDIF
C
      RETURN
      END
c
