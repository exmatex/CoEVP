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

MODULE vpETxt_interface_base
  !
  ! things shared by all aggregate level model implementations (VPSC and Taylor so far)

  IMPLICIT NONE

  INTEGER, PARAMETER, PUBLIC :: idp = SELECTED_REAL_KIND(15, 60)
  INTEGER, PARAMETER :: DIMS  = 3
  INTEGER, PARAMETER :: SVEC  = 6
  INTEGER, PARAMETER :: TVEC  = SVEC-1 !5
  INTEGER, PARAMETER :: WVEC  = DIMS !3
  !
  INTEGER, PARAMETER :: EULERDIM_p = 3
  !
  REAL(idp), PARAMETER :: &
       &   zero       = 0.0d0, &
       &   one        = 1.0d0, &
       &   two        = 2.0d0, &
       &   onehalf    = 0.5d0

  CHARACTER(21), PARAMETER, PRIVATE :: location='vpETxt_interface_base'

  ! this data type is meant for quantities which are fixed at the beginning
  ! of the calculations; whether set up by the aggregate level model or at a higher level
  !
  TYPE vpETxt_phase_data_type

    ! number of grains, for convenience
    INTEGER :: nGrains
    
    ! These first few things are inputs to the aggregate level model;
    !
    ! note that orientations are effectively fixed -- from the
    ! perspective of aggregate level model, texture evolution happens by changes
    ! in the weights
    REAL(idp), POINTER :: &
         & orientations(:,:), & ! (nOrParam, nGrains)
         & Cbin(:,:,:)          ! (3, 3, nGrains)
    !
    ! this is the volume fraction of the solid part of the material occupied
    ! by this phase (does not change if porosity is included)
    REAL(idp) :: volFrac
    !
    ! symmetry group code coming in can be checked against what
    ! aggregate level model reads from its input file as a sanity
    ! check
    INTEGER :: iSymmCode
    
    ! And these things need to be set by the aggregate level model;
    !
    ! this is the total number of systems, not just the number of twin types
    INTEGER :: nTwinSystems
    !
    ! twin orientation relationships: rotation matrices taking components
    ! in the parent crystal frame to components in the twin crystal frame;
    ! these are needed at the next level up to compute the non-local 
    ! contributions across orientation space when twinning happens;
    ! these _must_ be ordered to correspond with the order used for 
    ! twinVFRates in vpETxt_phase_response_type
    REAL(idp), POINTER :: twinOR(:,:,:) ! (DIMS, DIMS, nTwinSystems)
    !
    ! convenience variables for indexing into hVec and hVecDot
    INTEGER :: iHVLB, iHVUB

    ! And these may only be used by some aggregate level models
    !TYPE(crystal), POINTER :: crys
!RL    INTEGER :: nStrngth !NB: for vpsc
    INTEGER :: nStrengths ! RL name fixed
    
  END TYPE vpETxt_phase_data_type

  ! this data type is meant for quantities which evolve through the 
  ! calculations, whether set by aggregate level model or not
  !
  TYPE vpETxt_phase_response_type
    
    ! Inputs to aggregate level evaluation:
    !
    ! sum of each phase adds to 1; separate volume fraction for each phase
    REAL(idp), POINTER :: weights(:) ! (nGrains)
    !
    REAL(idp), POINTER :: dw_dOdfDof(:,:) ! (nGrain, nOdfDof)
    INTEGER :: nOdfDof ! for convenience

    ! Outputs from aggregate level model evaluation: 
    !
    ! volume fraction rate of twinning:
    REAL(idp), POINTER :: twinVFRates(:,:) ! (nTwinSystems, nGrains)
    !
    ! rates of crystal lattice spin in global coordinates (23, 31, 12):
    REAL(idp), POINTER :: reorientationRates(:,:) ! (WVEC, nGrains)

  END TYPE vpETxt_phase_response_type

  TYPE vpETxt_interface_type

    ! chosen orientation parameterization; see rot_params_mod
    INTEGER :: iSCR

    ! number of hardness state variables for the collection of phases/grains
    INTEGER :: nPhase, nH

    ! reference deformation rate
    REAL(idp) :: dRef

    ! phase data
    TYPE(vpETxt_phase_data_type), POINTER :: phaseDataList(:)

    ! information related to porosity:
    LOGICAL :: withPorosity
    ! convenience variables for indexing into hVec and hVecDot
    INTEGER :: porosityIHVLB, porosityIHVUB

    ! default initial values for hVec
    REAL(idp), POINTER :: hVecInit(:) ! nH

    ! And this may only be used by some aggregate level models
    !TYPE (nlsolve_tr_dogleg_nJ_control), POINTER :: nlControl

    LOGICAL :: enableUnitStress

    !NB: not necessarily used by all aggregate level models, here for convenience
    INTEGER :: nTwinSysMax, nStrngthMax, nGrTot

  END TYPE vpETxt_interface_type

CONTAINS

  SUBROUTINE ctor_init_vpETxt(this)
    IMPLICIT NONE
    TYPE (vpETxt_interface_type), INTENT(INOUT)  :: this

    this%nPhase = -1
    this%nH     = -1
    NULLIFY(this%phaseDataList)

    this%withPorosity = .FALSE.
    this%porosityIHVLB = -1
    this%porosityIHVUB = -1

    NULLIFY(this%hVecInit) ! this%nlControl

    this%enableUnitStress = .FALSE.

    this%nTwinSysMax = -1
    this%nStrngthMax = -1
    this%nGrTot      = -1
    
  END SUBROUTINE ctor_init_vpETxt
  !
  SUBROUTINE dtor_vpETxt(this)
    IMPLICIT NONE
    TYPE (vpETxt_interface_type), INTENT(INOUT)  :: this
    !
    INTEGER :: iPhase

    IF (ASSOCIATED(this%phaseDataList)) THEN
      DO iPhase = LBOUND(this%phaseDataList,DIM=1), UBOUND(this%phaseDataList,DIM=1)
        CALL dtor_vpETxtPD(this%phaseDataList(iPhase))
      END DO
      DEALLOCATE(this%phaseDataList)
    END IF

    IF (ASSOCIATED(this%hVecInit)) THEN
      DEALLOCATE(this%hVecInit)
    END IF

    ! IF (ASSOCIATED(this%nlControl)) THEN
    !   CALL dealloc_nl(this%nlControl)
    !   DEALLOCATE(this%nlControl)
    ! END IF

    CALL ctor_init_vpEtxt(this)
    
  END SUBROUTINE dtor_vpETxt

  SUBROUTINE ctor_init_vpETxtPD(this)
    IMPLICIT NONE
    TYPE (vpETxt_phase_data_type), INTENT(INOUT) :: this

    this%nGrains    = -1
    NULLIFY(this%orientations, this%Cbin)
    this%volFrac    = -1.0
    this%iSymmCode  = -1
    this%nTwinSystems = -1
    NULLIFY(this%twinOR)
    this%iHVLB      = -1
    this%iHVUB      = -1
!RL    this%nStrngth   = -1
    this%nStrengths   = -1
    !
    !NULLIFY(this%crys)

  END SUBROUTINE ctor_init_vpETxtPD
  !
  SUBROUTINE dtor_vpETxtPD(this)
    IMPLICIT NONE
    TYPE (vpETxt_phase_data_type), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%orientations)) THEN
      DEALLOCATE(this%orientations)
    END IF
    IF (ASSOCIATED(this%Cbin)) THEN
      DEALLOCATE(this%Cbin)
    END IF
    IF (ASSOCIATED(this%twinOR)) THEN
      DEALLOCATE(this%twinOR)
    END IF

    ! IF (ASSOCIATED(this%crys)) THEN
    !   ! ... would dealloc crys here if had code to do it
    !   NULLIFY(this%crys)
    ! END IF

    CALL ctor_init_vpETxtPD(this)

  END SUBROUTINE dtor_vpETxtPD

  SUBROUTINE ctor_init_vpETxtPR(this)
    IMPLICIT NONE
    TYPE (vpETxt_phase_response_type), INTENT(INOUT) :: this

    NULLIFY(this%weights, this%dw_dOdfDof, this%twinVFRates, this%reorientationRates)
    this%nOdfDof = 0

  END SUBROUTINE ctor_init_vpETxtPR
  !
  SUBROUTINE dtor_vpETxtPR(this)
    IMPLICIT NONE
    TYPE (vpETxt_phase_response_type), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%weights)) THEN
      DEALLOCATE(this%weights)
    END IF

    IF (ASSOCIATED(this%dw_dOdfDof)) THEN
      DEALLOCATE(this%dw_dOdfDof)
    END IF

    IF (ASSOCIATED(this%twinVFRates)) THEN
      DEALLOCATE(this%twinVFRates)
    END IF

    IF (ASSOCIATED(this%reorientationRates)) THEN
      DEALLOCATE(this%reorientationRates)
    END IF

    CALL ctor_init_vpETxtPR(this)

  END SUBROUTINE dtor_vpETxtPR

END MODULE vpETxt_interface_base

