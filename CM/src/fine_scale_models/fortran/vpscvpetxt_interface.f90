! -*-f90-*-
! DO-NOT-DELETE revisionify.begin() 
!
!   Copyright (c) 2007-2009 Lawrence Livermore National Security,
!   LLC. Produced at the Lawrence Livermore National Laboratory (Nathan
!   Barton <barton22@llnl.gov>) CODE-OCEC-08-104.
!   
!   Please also read the file NOTICES.
!   
!   This file is part of the mdef package (version 0.2) and is
!   free software: you can redistribute it and/or modify it under the
!   terms of the GNU Lesser General Public License as published by the
!   Free Software Foundation, either version 3 of the License, or (at your
!   option) any later version.
!   
!   This program is distributed in the hope that it will be useful, but
!   WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Lesser General Public License for more details.
!   
!   A copy of the GNU Lesser General Public License may be found in the
!   file NOTICES. If this file is missing, see
!   <http://www.gnu.org/licenses/>.
!
! DO-NOT-DELETE revisionify.end() 

! Assumptions, at least for now:
!
!       (*) Fully stress driven
!
!       (*) All state evolution is handled outside of VPSC, VPSC does
!       evalutions at fixed state; note that dt is not passed as an
!       argument
!
! Other considerations:
!
!       (*) It would be nice to have different rate sensitivities for
!       twinning and slip, especially for capturing slip/twinning
!       transitions. But then computing a meaningful aggregate average
!       rate sensitivity becomes perhaps more complicated. As
!       necessary, we can live with using the same rate sensitivity
!       for both twinning and slip for now. Actually, capturing
!       slip/twinning transitions may ultimately require a more
!       sophisticated slip and dislocation multiplication kinetics
!       which include drag effects.

MODULE vpscVpETxt_interface
  !
  USE vpETxt_interface_base

CONTAINS

  FUNCTION init_vpsc( &
       & it, &
       & fnameIn, & 
       & phaseDataList, &
       & MPICommF, c_scaling) RESULT(status)
    !
    !USE rot_params_mod, ONLY : BUNGE_p
    !
    IMPLICIT NONE
    !
    ! pass back non-zero on error:
    INTEGER :: status
    !
    ! type for carrying any information we want to keep handy:
    TYPE (vpETxt_interface_type), INTENT(INOUT)  :: it

    !
    ! file from which to read VPSC information:
    CHARACTER(len=*), INTENT(IN) :: fnameIn
    !
    ! data for this material, some of which gets set up in here
    TYPE (vpETxt_phase_data_type), INTENT(INOUT), TARGET :: phaseDataList(:)
    !
    ! fortran style MPI communicator number;
    ! not needed yet but included for completeness:
    INTEGER, INTENT(IN) :: MPICommF 
    real(idp), intent(in) :: c_scaling
    !
    ! ------------------------------------------------------------
    ! local variables
    !
    INTEGER :: iPhase, ngrtot
    REAL(idp) :: porosity
    REAL(idp), POINTER :: &
         & strengths(:,:), &  ! (nPhase)
         & shapes(:,:), &   ! (SVEC,nPhase)
         & orients(:,:), &  ! (nOrParam,ngrtot)
         & volFracPhase(:)  ! (nPhase)
    INTEGER, PARAMETER :: maxStrengthsPerPhase = 10 !NB: not elegant, but convenient
    INTEGER, POINTER :: &
         & nGrPhase(:), &   ! (nPhase)
!RL         & nTwinSys(:), &   ! (nPhase)
!RL         & nStrngths(:)     ! (nPhase)
         & nTwSysPhase(:), &  ! (nPhase)
         & nStrPhase(:)       ! (nPhase)
    !
    ! ------------------------------------------------------------

    status = 0

    !it%ISCR = BUNGE_p

    it%phaseDataList => phaseDataList
    it%nPhase = SIZE(it%phaseDataList, DIM=1)

    ngrtot = 0
    DO iPhase = 1, it%nPhase
      ngrtot = ngrtot + phaseDataList(iPhase)%nGrains
    END DO 
    it%nGrTot = ngrtot !NB: for convenience later, used in automatic allocations
    !
    ALLOCATE(&
         & orients(EULERDIM_p,ngrtot), &
         & strengths(maxStrengthsPerPhase, it%nPhase), &
         & shapes(SVEC,it%nPhase+1), & ! +1 in case of porosity
         & nGrPhase(it%nPhase), &
!RL         & nTwinSys(it%nPhase), &
!RL         & nStrngths(it%nPhase), &
         & nTwSysPhase(it%nPhase), &
         & nStrPhase(it%nPhase), &
         & volFracPhase(it%nPhase), &
         & STAT=status)
    IF (status /= 0) RETURN

    ngrtot = 0
    DO iPhase = 1, it%nPhase
      !
      nGrPhase(iPhase)     = phaseDataList(iPhase)%nGrains
      volFracPhase(iPhase) = phaseDataList(iPhase)%volfrac
      
      IF (SIZE(phaseDataList(iPhase)%orientations,DIM=1) /= EULERDIM_p) THEN
        status = -1; RETURN
      END IF
      !
      orients(:,ngrtot+1:ngrtot+phaseDataList(iPhase)%nGrains) = phaseDataList(iPhase)%orientations(:,:)

      ngrtot = ngrtot + phaseDataList(iPhase)%nGrains
      !
    END DO
    
    ! parse input file

    CALL vpsc_input_as (&
         & fnameIn, &           ! input
         & orients, &           ! input
         & it%nPhase, &         ! input
         & nGrPhase, &          ! input
         & volFracPhase, &      ! input
!RL
         & maxStrengthsPerPhase, & ! input
!
         & it%nGrTot, &         ! input
!RL
         & strengths, &         ! output
!RL         & nStrngths, &         ! output !NB: new output
         & nStrPhase, &         ! output
         & shapes, &            ! output
!RL         & nTwinSys, &       ! output
         & nTwSysPhase, &       ! output
         & porosity, &           ! output
!RL         & maxStrengthsPerPhase & ! parameter, for dimensioning !NB: new input
         & c_scaling &
         & )
    CALL load_conditions_as
    CALL update_Schmid
    !
    it%withPorosity = porosity > zero

    it%nH = 0
    it%nStrngthMax = 0
    it%nTwinSysMax = 0
    !
    phase_loop : DO iPhase = 1, it%nPhase
      !
!RL      it%nStrngthMax = MAX(it%nStrngthMax, nStrngths(iPhase)) !NB: for convenience later
!RL      it%nTwinSysMax = MAX(it%nTwinSysMax, nTwinSys(iPhase))  !NB: for convenience later
      it%nStrngthMax = MAX(it%nStrngthMax, nStrPhase(iPhase))
      it%nTwinSysMax = MAX(it%nTwinSysMax, nTwSysPhase(iPhase))

      phaseDataList(iPhase)%nTwinSystems = nTwSysPhase(iPhase)

      IF (phaseDataList(iPhase)%nTwinSystems > 0) THEN
        ALLOCATE( &
             & phaseDataList(iPhase)%twinOR(DIMS, DIMS, phaseDataList(iPhase)%nTwinSystems), &
             & STAT=status)
        IF (status /= 0) RETURN
!RL        CALL get_twinor_as (iPhase, phaseDataList(iPhase)%twinOR)
        CALL get_twinor_as (iPhase, phaseDataList(iPhase)%twinOR, phaseDataList(iPhase)%nTwinSystems)
      ELSE
        NULLIFY(phaseDataList(iPhase)%twinOR)
      END IF

      ! want storage for strength and average grain shape for this phase
      phaseDataList(iPhase)%iHVLB = it%nH + 1
!RL      phaseDataList(iPhase)%iHVUB = it%nH + nStrngths(iPhase) + 1 + SVEC
!RL      it%nH = it%nH + nStrngths(iPhase) + 1 + SVEC
      phaseDataList(iPhase)%iHVUB = it%nH + nStrPhase(iPhase) + 1 + SVEC
      it%nH = it%nH + nStrPhase(iPhase) + 1 + SVEC

!RL      phaseDataList(iPhase)%nStrngth = nStrngths(iPhase)
      phaseDataList(iPhase)%nStrengths = nStrPhase(iPhase)

    END DO phase_loop
    !
    IF (it%withPorosity) THEN
      it%porosityIHVLB = it%nH + 1
      it%porosityIHVUB = it%nH + 1 + SVEC
      it%nH = it%nH + 1 + SVEC
    ELSE
      it%porosityIHVLB = -1
      it%porosityIHVUB = -1
    END IF

!    WRITE(*,*) 'GALEN about to allocate hVecInit with nH = ', it%nH
    ALLOCATE( &
         & it%hVecInit(it%nH), &
         & STAT=status)
    IF (status /= 0) RETURN
    !
    phase_loop_hVecInit : DO iPhase = 1, it%nPhase

      ! set default initial strength for the phase
      it%hVecInit(phaseDataList(iPhase)%iHVLB:&
!RL           & phaseDataList(iPhase)%iHVLB+phaseDataList(iPhase)%nStrngth-1) = &
!RL           & strengths(1:phaseDataList(iPhase)%nStrngth, iPhase)
           & phaseDataList(iPhase)%iHVLB+phaseDataList(iPhase)%nStrengths-1) = &
           & strengths(1:phaseDataList(iPhase)%nStrengths, iPhase)
      !
      ! set to zero initial gamma for the phase
!RL      it%hVecInit(phaseDataList(iPhase)%iHVLB+phaseDataList(iPhase)%nStrngth) = 0.
      it%hVecInit(phaseDataList(iPhase)%iHVLB+phaseDataList(iPhase)%nStrengths) = 0.  
      !
      !
      ! set default initial grain shape for the phase
!RL      it%hVecInit(phaseDataList(iPhase)%iHVLB+phaseDataList(iPhase)%nStrngth+1:&
      it%hVecInit(phaseDataList(iPhase)%iHVLB+phaseDataList(iPhase)%nStrengths+1:&
           & phaseDataList(iPhase)%iHVUB) = &
           & shapes(:,iPhase)

    END DO phase_loop_hVecInit
    !
    IF (it%withPorosity) THEN

      ! set default initial porosity parameters for the material
      !
      !NB:
      ! You had shapes(:,it%nPhase) below, but it%nPhase is the number of solid phases --
      ! there is a lot of code which would have to change elsewhere if it%nPhase were changed
      ! to include a count of the porosity as a phase; I have therefore changed the index on 
      ! shapes below to it%nPhase+1 hoping that is correct. It is possible that the use of nPhase 
      ! needs to change elsewhere too, but I do not know enough about how vpsc_input is supposed to
      ! work to know for sure 
      !
      !RL:
      !OK
      it%hVecInit(it%porosityIHVLB) = porosity
      it%hVecInit(it%porosityIHVLB+1:it%porosityIHVUB) = &
           & shapes(:,it%nPhase+1)

    END IF

    DEALLOCATE( &
         & strengths, &
         & shapes, &
         & orients, &
         & volFracPhase, &
         & nGrPhase, &
!RL         & nTwinSys, &
!RL         & nStrngths)
         & nTwSysPhase, &
         & nStrPhase)

  END FUNCTION init_vpsc

  FUNCTION evalEvolution_vpsc( &
       & it, &
       & tK, &
       & stressSvecP, &
       & hVec, &
       & L_o, m, g, & 
       & hVecDot, &
       & phaseResponseList &
       & ) RESULT(status)
    !
    IMPLICIT NONE
    !
    ! pass back non-zero on error:
    INTEGER                                 :: status
    !
    TYPE (vpETxt_interface_type), INTENT(IN)  :: it
    !
    ! absolute temperature:
    REAL(idp), INTENT(IN)                   :: tK
    !
    ! applied stress;
    ! stress components: 11', 22', 33', 23', 31', 12', p=-tr(sigma)/3
    ! deviatoric components are of a Kirchhoff type
    REAL(idp), INTENT(IN)                   :: stressSvecP(SVEC+1)
    !
    ! vector of state variables
    REAL(idp), INTENT(IN), TARGET           :: hVec(it%nH) 
    !
    ! vector of rate of change of state variables
    REAL(idp), INTENT(OUT), TARGET          :: hVecDot(it%nH)
    ! for each phase:
    !   hVecDot(phaseDataList%iHVLB:phaseDataList%iHVUB)
    ! if (it%withPorosity), need to have set:
    !   hVecDot(it%porosityIHVLB:it%porosityIHVUB)
    !
    ! data types for evolving quantities, per phase, contains both 
    ! inputs and outputs
    TYPE(vpETxt_phase_response_type), INTENT(INOUT), TARGET :: phaseResponseList(it%nPhase)
    ! outputs for each iPhase:
    !   phaseResponseList(iPhase)%twinVFRates(:,:)
    !   phaseResponseList(iPhase)%reorientationRates(:,:)
    !
    ! homogenized plastic flow response parameters,
    ! see Equation 33 in http://dx.doi.org/10.1016/j.ijplas.2007.03.004
    !     reference velocity gradient
    REAL(idp), INTENT(OUT)                  :: L_o(DIMS,DIMS)
    !     homogenized rate sensitivity
    REAL(idp), INTENT(OUT)                  :: m
    !     homogenized scalar flow strength
    REAL(idp), INTENT(OUT)                  :: g
    !
    ! ------------------------------------------------------------
    ! locals
    !
    REAL(idp) :: solidVolFrac, porosity,porosityRate 
    REAL(idp), POINTER :: poreShape(:) ! (SVEC)
    !
    INTEGER   :: iPhase
    INTEGER   :: ngrtot
    !
    ! these are automatically allocated:
    REAL(idp) :: &
         & strengths(it%nStrngthMax, it%nPhase), & 
         & gammas(it%nPhase), &
         & shapes(SVEC,it%nPhase+1), & ! +1 in case of porosity
         & wgts(it%nGrTot), &
         & strengthRates(it%nStrngthMax, it%nPhase), &
         & gammaRates(it%nPhase), &
         & shapeRates(SVEC,it%nPhase+1), & ! +1 in case of porosity
         & reorRates(WVEC,it%nGrTot), &
         & twVFRates(it%nTwinSysMax, it%nGrTot) !NB: switched to nTwinSysMax
    !
    TYPE(vpETxt_phase_response_type), POINTER :: phaseResponse
    TYPE(vpETxt_phase_data_type), POINTER :: phaseData
    !
    REAL(idp) :: L(DIMS, DIMS), D(DIMS,DIMS), dMag, sMag

    status = 0

    NULLIFY(phaseResponse, phaseData)

    IF (it%withPorosity) THEN
      porosity     =  hVec(it%porosityIHVLB)
      solidVolFrac =  one - porosity
      poreShape    => hVec(it%porosityIHVLB+1:it%porosityIHVUB)
      ! poreShape(:) is now a 6-vector for the shape of the pore
    ELSE
      porosity = zero !NB: I added this just for completeness
      solidVolFrac = one
      NULLIFY(poreShape)
    END IF

    ngrtot = 0
    DO iPhase = 1, it%nPhase
      !
      wgts(ngrtot+1:ngrtot+it%phaseDataList(iPhase)%nGrains) = &
           &  phaseResponseList(iPhase)%weights(:) * &
           &   ( solidVolFrac * it%phaseDataList(iPhase)%volFrac )
      !
      ngrtot = ngrtot + it%phaseDataList(iPhase)%nGrains
      !
    END DO

    phase_loop : DO iPhase = 1, it%nPhase

      ! assigments just for convenience
      !
      phaseResponse => phaseResponseList(iPhase)
      phaseData     => it%phaseDataList(iPhase)

!RL      strengths(1:phaseData%nStrngth,iPhase) =  hVec(phaseData%iHVLB:phaseData%iHVLB+phaseData%nStrngth-1)
!RL      gammas(iPhase) =  hVec(phaseData%iHVLB+phaseData%nStrngth)
!RL      shapes(:,iPhase)  = hVec(phaseData%iHVLB+phaseData%nStrngth+1:phaseData%iHVUB)
      strengths(1:phaseData%nStrengths,iPhase) =  &
   &  hVec(phaseData%iHVLB:phaseData%iHVLB+phaseData%nStrengths-1)
      gammas(iPhase) =  hVec(phaseData%iHVLB+phaseData%nStrengths)
      shapes(:,iPhase)  = hVec(phaseData%iHVLB+phaseData%nStrengths+1:phaseData%iHVUB)
    END DO phase_loop
    !
    NULLIFY(phaseResponse, phaseData)

    CALL vpsc_evol_as (&
         & stressSvecP, &       ! input
         & wgts, &              ! input
         & strengths, &         ! input
         & gammas, &            ! input
         & shapes, &            ! input
         & porosity, &          ! input
!RL
         & it%nPhase, &         ! input
         & it%nGrTot, &         ! input
         & it%nStrngthMax, &    ! input
         & it%nTwinSysMax, &    ! input
!RL
         & L,m, &               ! output
         & strengthRates, &     ! output
         & gammaRates, &        ! output
         & shapeRates, &        ! output
         & porosityRate, &      ! output
         & reorRates, &         ! output
         & twVFRates &         ! output
!RL         & it%nStrngthMax &     ! input !NB: new input for dimensioning
         & )

    ngrtot = 0
    !
    phase_loop_2 : DO iPhase = 1, it%nPhase

      phaseData  => it%phaseDataList(iPhase)

      phaseResponseList(iPhase)%twinVFRates &
           & (1:phaseData%nTwinSystems,1:phaseData%nGrains)= &
           & twVFRates(1:phaseData%nTwinSystems, & !NB: indexing changed
           &           ngrtot+1:ngrtot+phaseData%nGrains)

      phaseResponseList(iPhase)%reorientationRates &
           & (:,1:phaseData%nGrains)= &
           & reorRates(:,ngrtot+1:ngrtot+phaseData%nGrains)
           
      ngrtot = ngrtot + phaseData%nGrains

!RL      hVecDot(phaseData%iHVLB:phaseData%iHVLB+phaseData%nStrngth-1) = & !NB: indexing changed
!RL           & strengthRates(1:phaseData%nStrngth,iPhase) 
!RL      hVecDot(phaseData%iHVLB+phaseData%nStrngth)=gammaRates(iPhase) !NB: indexing changed
!RL      hVecDot(phaseData%iHVLB+phaseData%nStrngth+1:phaseData%iHVUB)=shapeRates(:,iPhase) !NB: indexing changed

      hVecDot(phaseData%iHVLB:phaseData%iHVLB+phaseData%nStrengths-1) = &
           & strengthRates(1:phaseData%nStrengths,iPhase) 
      hVecDot(phaseData%iHVLB+phaseData%nStrengths)=gammaRates(iPhase)
      hVecDot(phaseData%iHVLB+phaseData%nStrengths+1:phaseData%iHVUB)=shapeRates(:,iPhase)

    END DO phase_loop_2
    !
    IF (it%withPorosity) THEN
      hVecDot(it%porosityIHVLB) = porosityRate
      hVecDot(it%porosityIHVLB+1:it%porosityIHVUB) = &
           & shapeRates(:,it%nPhase+1)
    END IF

    ! have L(:,:) from vpsc calculation;
    ! have computed some meaningful effective m;
    !
    D = onehalf * (L + TRANSPOSE(L))
    !
    ! note that these magnitudes do not, on purpose, have funny
    ! von-Mises TYPE factors in them
    dMag = SQRT( &
         &        D(1,1)*D(1,1) + D(2,2)*D(2,2) + D(3,3)*D(3,3) + &
         & two * (D(2,3)*D(2,3) + D(3,1)*D(3,1) + D(1,2)*D(1,2))  &
         & )
    sMag = SQRT( &
         &       SUM(stressSvecP(1:3)*stressSvecP(1:3)) + &
         & two * SUM(stressSvecP(4:6)*stressSvecP(4:6))   &
         & )
    !
    ! set reference velocity gradient such that magnitude of deformation rate
    ! part comes out to dRef:
    !L_o(:,:) = L(:,:) * it%dRef / dMag
    ! try returning only the symmetric part to conform to the spec : Jamal

    ! this is the problem line JAMAL
    ! L_o(:,:) = D(:,:) * it%dRef / dMag
    L_o(:,:) = D(:,:) 

    !
    ! have 
    !   dMag / dRef = (sMag / g) ** (1/m)
    ! so that: 
    !   sMag / g = (dMag / dRef) ** m
    ! so that: 
    g = sMag * (it%dRef / dMag) ** m

  END FUNCTION evalEvolution_vpsc

!RL  SUBROUTINE getTauStarParams_vpsc(&
!RL       & it, &
!RL       & hVec, tK, tK_ref, minH, xm_eff)
!RL    !
!RL    IMPLICIT NONE
!RL    !
!RL    TYPE (vpETxt_interface_type), INTENT(IN)  :: it
!RL    !
!RL    ! current history variable values:
!RL    REAL(idp), INTENT(IN)  :: hVec(it%nH)
!RL    !
!RL    ! current and reference temperatures (on an absolute scale):
!RL    REAL(idp), INTENT(IN)  :: tK, tK_ref
!RL    !
!RL    ! strength of the softest deformation mode:
!RL    REAL(idp), INTENT(OUT) :: minH
!RL    !
!RL    ! estimate of the effective rate sensitivity at the current temperature:
!RL    REAL(idp), INTENT(OUT) :: xm_eff
!RL    !
!RL    ! ------------------------------------------------------------
!RL
!RL    ! ... implementation here!
!RL
!RL  END SUBROUTINE getTauStarParams_vpsc

END MODULE vpscVpETxt_interface

