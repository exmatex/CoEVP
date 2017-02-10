subroutine printIt (it, phaseDataList, phaseResponseList)

   !
   USE vpETxt_interface_base
   USE vpscVpETxt_interface
   !
   IMPLICIT NONE
   !
   type (vpETxt_interface_type) :: it
   type (vpETxt_phase_data_type) :: phaseDataList(it%nPhase)
   type (vpETxt_phase_response_type) :: phaseResponseList(it%nPhase)

   integer :: iPhase, i, j

   ! print the diagnostics here to make sure we got it all
      write(*,*) 'Interface data'
      write(*,*) 'dRef = ', it%dRef
      write(*,*) 'nPhase = ', it%nPhase
      write(*,*) 'nH = ', it%nH
      write(*,*) 'nGrTot= ', it%nGrTot
      write(*,*) 'nTwinSysMax = ', it%nTwinSysMax
      write(*,*) 'nStrngthMax = ', it%nStrngthMax
      write(*,*) 'withPorosity= ', it%withPorosity
      write(*,*) 'porosityIHVLB = ', it%porosityIHVLB
      write(*,*) 'porosityIHVUB = ', it%porosityIHVUB
      write(*,*) '   hVecInit'
      do i=1,it%nH
         write(*,*) it%hVecInit(i)
      enddo

      write(*,*) 'Phase Data'
      write(*,*) 'volFrac = ',(phaseDataList(i)%volFrac,i=1,it%nPhase)
      write(*,*) 'nGrains = ',(phaseDataList(i)%nGrains,i=1,it%nPhase)
      write(*,*) 'iHVLB = ',(phaseDataList(i)%iHVLB,i=1,it%nPhase)
      write(*,*) 'iHVUB = ',(phaseDataList(i)%iHVUB,i=1,it%nPhase)
      write(*,*) 'nStrPhase = ',(phaseDataList(i)%nStrengths,i=1,it%nPhase)
      write(*,*) 'nTwSysPhase = ',(phaseDataList(i)%nTwinSystems,i=1,it%nPhase)

      write(*,*) '   Grain orientstions, volume fractions:'
      do iPhase = 1, it%nPhase
         write(*,*) 'Phase ',iPhase
         do i = 1,phaseDataList(iPhase)%nGrains
            write(*,*) i, (phaseDataList(iPhase)%orientations(j,i), j=1,3), phaseResponseList(iPhase)%weights(i)
         enddo
      enddo

end subroutine printIt


      SUBROUTINE vpsc_init(&
!     interface data
      dRef,                &
      nPhase,              &
      nH,                  &
      hVecInit,            &
      porosityIHVLB,       &
      porosityIHVUB,       &
      nTwinSysMax,         &
      nStrngthMax,         &
      nGrTot,              &
      intPorosity,         &
!     phase data
      nGrains,             &
      orients,             &
      volFrac,             &
      nStrPhase,           &
      nTwSysPhase,         &
      iHVLB,               &
      iHVUB,               &
      strengths,           &
!     phase response data
      weights,             &
!     dw_dOdfDof,          &  ! apparently never used
!     nOdfDof,             &  ! apparently never used
      twinVFRates,         &
      reorientationRates,  &
!     extra
      c_scaling,           &
      diagnostics          &
      )


!     
      USE vpETxt_interface_base
      USE vpscVpETxt_interface
!
 
      IMPLICIT NONE
!     
!     !  integer :: init_vpsc, evalEvolution_vpsc
!     
!     init_vpsc and evalEvolution_vpsc arguments
!     
      INTEGER, PARAMETER :: maxStrengthsPerPhase = 10 !NB: not elegant, but convenient
      integer, parameter :: nOdfDof = 300
!     interface data
      real(idp):: dRef
      integer :: nPhase, nH
      real(idp) :: hVecInit(nH)
      integer :: porosityIHVLB
      integer :: porosityIHVUB
      integer :: nTwinSysMax, nStrngthMax, nGrTot
      logical :: withPorosity
!     phase data
      integer :: nGrains(nPhase)
      real(idp) :: orients(EULERDIM_p*ngrtot)
      real(idp) :: volFrac(nPhase)
      integer :: nStrPhase(nPhase)
      integer :: nTwSysPhase(nPhase)
      integer :: iHVLB(nPhase)
      integer :: iHVUB(nPhase)
      real(idp) :: strengths(nPhase*maxStrengthsPerPhase)
!     phase response data
      real(idp) :: weights(nGrTot)
      real(idp) :: dw_dOdfDof(nGrTot*nOdfDof)
!     integer :: nOdfDof
      real(idp) :: twinVFRates(nTwinSysMax*nGrTot)
      real(idp) :: reorientationRates(WVEC*nGrTot)

      integer :: intPorosity
      integer :: status
      type (vpETxt_interface_type) :: it
!     
!     init_vpsc arguments
!     
      character(len=255) :: fnameIn
      type (vpETxt_phase_data_type),allocatable :: phaseDataList(:)
      integer :: MPICommF
!     
!     evalEvolution_vpsc arguments
!     
      type (vpETxt_phase_response_type),allocatable :: phaseResponseList(:)
      real(idp) tk,stressSvecP(SVEC+1), L_o(DIMS,DIMS),m,g
      real(idp) D(DIMS,DIMS)
!     
      integer :: iPhase, iWeight
!     
      integer :: kgr,i,j,f_eval, iGlobal
      integer :: diagnostics
      real(idp) :: c_scaling
      real(idp) :: wgtTot(nPhase)
      character(len=255) :: dataDir

      call get_environment_variable("VPSC_INPUT_PATH", dataDir)
      if (dataDir == "") then 
         dataDir = '../../CoEVP/CM/src/fine_scale_models/tantalum/'
      endif
      
      
      
!     
!     
      if (diagnostics == 1) then
         write(*,*) 'nPhase passed in = ', nPhase
      endif

      fnameIn=trim(adjustl(dataDir))//'vpsc_as_try.in'

      open(62,FILE=fnameIn,STATUS='OLD')

      read(62,*) iPhase
      
      if (diagnostics == 1) then
         write(*,*) 'nPhase read in = ', iPhase
      endif

      it%nH=nH 
      it%nPhase=nPhase
      it%dRef=1.
      dRef = 1.

      if (nPhase > 2) then
         write(*,*) 'No more than 2 phases allowed right now!'
         stop
      endif

      allocate(phaseDataList(nPhase))
      allocate(phaseResponseList(nPhase))
!     
!     
      read(62,*) (phaseDataList(iPhase)%volfrac, iPhase=1,nPhase)
      close(62)
      

!     worst piece of code ever
      if (nPhase == 1) then 
         open(55,FILE=trim(adjustl(dataDir))//'rand419.tex',status='old')
         read(55,*) nGrains(1)
      endif
      if (diagnostics == 1) then
         write(*,*) (nGrains(i), i=1,nPhase)
      endif
      nGrTot = 0
      
      if (nPhase == 1) then 
         wgtTot(1) = 0.0
                  
         do iPhase = 1, nPhase
            nGrTot = nGrTot + nGrains(iPhase)
            phaseDataList(iPhase)%ngrains=nGrains(iPhase)
            allocate(phaseDataList(iPhase)%orientations(3,nGrains(iPhase)))
            allocate(phaseResponseList(iPhase)%weights(nGrains(iPhase)))
            allocate(phaseResponseList(iPhase)%reorientationRates(3,nGrains(iPhase)))
            allocate(phaseResponseList(iPhase)%twinVFrates(12,nGrains(iPhase)))
         enddo
         
         do kgr=1,nGrains(1)
            read(55,*) (phaseDataList(1)%orientations(i,kgr),i=1,3), &
            &           phaseResponseList(1)%weights(kgr)
            wgtTot(1)=wgtTot(1)+phaseResponseList(1)%weights(kgr)
         enddo
         
         close(55)
         
         
      else if (nPhase == 2) then
         wgtTot(1) = 0.0
         wgtTot(2) = 0.0
         
         open(55,file=trim(dataDir)//'orient_ph1.in',status='old')
         open(56,file=trim(dataDir)//'orient_ph2.in',status='old')
!     
         read(55,*) nGrains(1)
         read(56,*) nGrains(2)
         
         if (diagnostics == 1) then
            write(*,*) (nGrains(i), i=1,nPhase)
         endif
         nGrTot = 0

         do iPhase = 1, nPhase
            nGrTot = nGrTot + nGrains(iPhase)
            phaseDataList(iPhase)%ngrains=nGrains(iPhase)
            allocate(phaseDataList(iPhase)%orientations(3,nGrains(iPhase)))
            allocate(phaseResponseList(iPhase)%weights(nGrains(iPhase)))
            allocate(phaseResponseList(iPhase)%reorientationRates(3,nGrains(iPhase)))
            allocate(phaseResponseList(iPhase)%twinVFrates(12,nGrains(iPhase)))
         enddo
         
         do kgr=1,nGrains(1)
            read(55,*) (phaseDataList(1)%orientations(i,kgr),i=1,3), &
            &           phaseResponseList(1)%weights(kgr)
            wgtTot(1)=wgtTot(1)+phaseResponseList(1)%weights(kgr)
         enddo
         close(55)
         do kgr=1,nGrains(2)
            read(56,*) (phaseDataList(2)%orientations(i,kgr),i=1,3), &
            &           phaseResponseList(2)%weights(kgr)
            wgtTot(2)=wgtTot(2)+phaseResponseList(2)%weights(kgr)
         enddo
         close(56)
      endif
!     
      if (diagnostics == 1) then
         write(*,*) 'Total of weights per phase:'
         write(*,*) (wgtTot(i), i=1,2)
      endif
!     normalize the weights if needed
      do iPhase = 1, nPhase
         do kgr=1,nGrains(iPhase)
            phaseResponseList(iPhase)%weights(kgr)=phaseResponseList(iPhase)%weights(kgr)/wgtTot(iPhase)
         enddo
      enddo
      
      f_eval=init_vpsc(it,fnameIn,phaseDataList,MPICommF, c_scaling)

      if (f_eval /= 0) then
         write(*,*) 'PROBLEM IN init_vpsc'
         stop
      endif
      
      if (diagnostics == 1) then
         write(*,*) 'init_vpsc completed'
!     
!     generate the flattened arrays to be passed up to the C code
!     interface type
         write(*,*) 'dRef = ', dRef
         write(*,*) 'nPhase = ', nPhase
         write(*,*) 'nH = ', nH
         write(*,*) 'nGrTot= ', nGrTot
         write(*,*) 'nTwinSysMax = ', nTwinSysMax
         write(*,*) 'nStrngthMax = ', nStrngthMax

         write(*,*) 'After running init_vpsc, interface values are:'
         write(*,*) 'dRef = ', it%dRef
         write(*,*) 'nPhase = ', it%nPhase
         write(*,*) 'nH = ', it%nH
         write(*,*) 'nGrTot= ', it%nGrTot
         write(*,*) 'nTwinSysMax = ', it%nTwinSysMax
         write(*,*) 'nStrngthMax = ', it%nStrngthMax
         write(*,*) 'withPorosity= ', it%withPorosity
         write(*,*) 'porosityIHVLB = ', porosityIHVLB
         write(*,*) 'porosityIHVUB = ', porosityIHVUB
      endif
      
!     convert fortran logical to C int
      intPorosity = 0;
!     if (withPorosity) then intPorosity = 1;
      
      dRef = it%dRef 
      nH = it%nH
      nGrTot = it%nGrTot
      nTwinSysMax = it%nTwinSysMax
      nStrngthMax = it%nStrngthMax
      withPorosity = it%withPorosity
      porosityIHVLB = it%porosityIHVLB
      porosityIHVUB = it%porosityIHVUB

!      write(*,*) 'GALEN about to loop through hVecInit with nH = ', nH
!      write(*,*) '      and it%nH = ', it%nH
      do i = 1, nH
         hVecInit(i) = it%hVecInit(i)
      enddo
      
      iGlobal=1
      iWeight=1
!     phase date list
      do iPhase = 1, nPhase
         volFrac(iPhase) = phaseDataList(iPhase)%volfrac
         nGrains(iPhase) = phaseDataList(iPhase)%nGrains
         iHVLB(iPhase) = phaseDataList(iPhase)%iHVLB
         iHVUB(iPhase) = phaseDataList(iPhase)%iHVUB
         nStrPhase(iPhase) = phaseDataList(iPhase)%nStrengths
         nTwSysPhase(iPhase) = phaseDataList(iPhase)%nTwinSystems
         
!     write orientatations to flat array
         do i=1, nGrains(iPhase)
            do j=1, EULERDIM_p
               orients(iGlobal) = phaseDataList(iPhase)%orientations(j,i)
               iGlobal = iGlobal+1
            enddo
            weights(iWeight) = phaseResponseList(iPhase)%weights(i)
            iWeight = iWeight + 1
         enddo
         
      enddo
      
      if (diagnostics == 1) then
         
!     print the diagnostics here to make sure we got it all
         write(*,*) 'At the end of the fortran function vpsc_init;'
         write(*,*) 'Interface data'
         write(*,*) 'dRef = ', dRef
         write(*,*) 'nPhase = ', nPhase
         write(*,*) 'nH = ', nH
         write(*,*) 'nGrTot= ', nGrTot
         write(*,*) 'nTwinSysMax = ', nTwinSysMax
         write(*,*) 'nStrngthMax = ', nStrngthMax
         write(*,*) 'withPorosity= ', withPorosity
         write(*,*) '   hVecInit'
         do i=1,nH
            write(*,*) hVecInit(i)
         enddo
         
         write(*,*) 'Phase Data'
         write(*,*) 'volFrac = ',(volFrac(i),i=1,nPhase)
         write(*,*) 'nGrains = ',(nGrains(i),i=1,nPhase)
         write(*,*) 'iHVLB = ',(iHVLB(i),i=1,nPhase)
         write(*,*) 'iHVUB = ',(iHVUB(i),i=1,nPhase)
         write(*,*) 'nStrPhase = ',(nStrPhase(i),i=1,nPhase)
         write(*,*) 'nTwSysPhase = ',(nTwSysPhase(i),i=1,nPhase)
         
         write(*,*) '   Grain orientstions, weights:'
         do i =1,nGrTot
            iGlobal = EULERDIM_p*(i-1)
            write(*,*) i, orients(iGlobal+1), orients(iGlobal+2), orients(iGlobal+3), weights(i)
         enddo
         

         write(*,*) 'Got to end of vpsc_init'
         
      endif
      
      end subroutine vpsc_init
      
subroutine vpsc_run(    &
      stressIn,         &
      strainOut,        &
      ! interface variables
      m, g, &
      dRef,    &
      nPhase, nH, &
      hVecInit,            &
      porosityIHVLB,       &
      porosityIHVUB,       &
      nTwinSysMax,   &
      nStrngthMax,   &
      nGrTot,        &
      withPorosity,  &
      ! phase data
      nGrains,             &
      orients,       &
      volFrac,       &
      !shapes,        &
      nStrPhase,     &
      nTwSysPhase,    &
      iHVLB, iHVUB,        &
      strengths,     &
      ! phase response data
      weights,             &
      twinVFRates,         &
      reorientationRates,  &
      diagnostics          &
      )
   !
   USE vpETxt_interface_base
   USE vpscVpETxt_interface
   !
   IMPLICIT NONE
   !
   INTEGER, PARAMETER :: maxStrengthsPerPhase = 10 !NB: not elegant, but convenient

   ! input temp variables
   ! multidimensional arrays unrolled flat
   real(idp):: dRef
   integer :: nPhase, nH, nTwinSysMax, nStrngthMax, nGrTot
   logical :: withPorosity
   real(idp) :: hVecInit(nH)
   integer :: porosityIHVLB
   integer :: porosityIHVUB
   integer :: nStrPhase(nPhase)
   integer :: nTwSysPhase(nPhase)
   real(idp) :: volFrac(nPhase)
   integer :: iHVLB(nPhase)
   integer :: iHVUB(nPhase)
   real(idp) :: strengths(nPhase*maxStrengthsPerPhase)
   !real(idp) :: shapes(nPhase*SVEC)
   integer :: nGrains(nPhase)
   real(idp) :: orients(EULERDIM_p*ngrtot)
   real(idp) :: weights(nGrTot)
   real(idp) :: twinVFRates(nTwinSysMax*nGrTot)
   real(idp) :: reorientationRates(WVEC*nGrTot)

   integer:: i,j, f_eval, status
   integer :: iPhase
   real :: tStart, tEnd
   real(idp) tk,stressSvecP(SVEC+1), L_o(DIMS,DIMS),m,g
   real(idp), pointer :: hVec(:), hVecDot(:)
   real(idp) stressIn(6), strainOut(6)
   real(idp) :: scaling
   type (vpETxt_interface_type) :: it
   type (vpETxt_phase_data_type),allocatable, target :: phaseDataList(:)
   type (vpETxt_phase_response_type),allocatable :: phaseResponseList(:)

   integer :: iGlobal, iWeight
   integer :: diagnostics

   withPorosity = .false.

   if (diagnostics == 1) then 
      write(*,*) 'Starting vpsc_run'
   endif

      call cpu_time(tStart)
   ! scale the stress
      scaling = 1.0
   ! index rotation to conform to VPSC internal notation
      stressSvecP=0.
      ! old rotation
      !stressSvecP(1)=stressIn(1)*scaling
      !stressSvecP(2)=stressIn(3)*scaling
      !stressSvecP(3)=stressIn(6)*scaling
      !stressSvecP(4)=stressIn(5)*scaling
      !stressSvecP(5)=stressIn(4)*scaling
      !stressSvecP(6)=stressIn(2)*scaling

      ! new rotation as per Milo
      stressSvecP(1)=stressIn(1)*scaling
      stressSvecP(2)=stressIn(2)*scaling
      stressSvecP(3)=stressIn(3)*scaling
      stressSvecP(4)=stressIn(5)*scaling
      stressSvecP(5)=stressIn(6)*scaling
      stressSvecP(6)=stressIn(4)*scaling

      stressSvecP(7)=0.*(1./3.)*scaling

   ! allocate all the structs used by the interface structs

   allocate(phaseDataList(nPhase))
   allocate(phaseResponseList(nPhase))

   ! populate the structs from the data passed in
   ! phaseDataList

   iGlobal=1
   iWeight=1
   do iPhase = 1, nPhase
      phaseDataList(iPhase)%ngrains=nGrains(iPhase)
      allocate(phaseDataList(iPhase)%orientations(3,nGrains(iPhase)))
      allocate(phaseResponseList(iPhase)%reorientationRates(3,nGrains(iPhase)))
      allocate(phaseResponseList(iPhase)%twinVFrates(12,nGrains(iPhase)))
      allocate(phaseResponseList(iPhase)%weights(nGrains(iPhase)))
      do i = 1,nGrains(iPhase)
      do j = 1,3
         phaseDataList(iPhase)%orientations(j,i) = orients(iGlobal)
         iGlobal = iGlobal + 1
      enddo
      phaseResponseList(iPhase)%weights(i) = weights(iWeight)
      iWeight = iWeight + 1
      enddo
   enddo

   it%nStrngthMax = nStrngthMax
   it%nTwinSysMax = nTwinSysMax
   it%nGrTot = nGrTot
   it%dRef = dRef
   it%nPhase = nPhase
   it%withPorosity = withPorosity
   it%porosityIHVLB = porosityIHVLB
   it%porosityIHVUB = porosityIHVUB

   it%phaseDataList => phaseDataList

   allocate(it%hVecInit(nH))
   do i=1,nH
      it%hVecInit(i) = hVecInit(i)
   enddo

   do iPhase = 1, nPhase
      
      phaseDataList(iPhase)%volFrac = volFrac(iPhase)
      phaseDataList(iPhase)%iHVLB = iHVLB(iPhase)
      phaseDataList(iPhase)%iHVUB = iHVUB(iPhase)

      phaseDataList(iPhase)%nGrains = nGrains(iPhase)
      phaseDataList(iPhase)%nStrengths = nStrPhase(iPhase)
      phaseDataList(iPhase)%nTwinSystems = nTwSysPhase(iPhase)


      IF (phaseDataList(iPhase)%nTwinSystems > 0) THEN
        ALLOCATE( &
             & phaseDataList(iPhase)%twinOR(DIMS, DIMS, phaseDataList(iPhase)%nTwinSystems), &
             & STAT=status)
        IF (status /= 0) RETURN
        CALL get_twinor_as (iPhase, phaseDataList(iPhase)%twinOR, phaseDataList(iPhase)%nTwinSystems)
      ELSE
        NULLIFY(phaseDataList(iPhase)%twinOR)
      END IF


   enddo


   if (diagnostics == 1) then 

      write(*,*) 'In vpsc_run, before calling vpsc'
   ! print the diagnostics here to make sure we got it all
      !call printIt(it, phaseDataList, phaseResponseList)

      write(*,*) 'stressIn, stressSvecP: '
      do i=1,7
         write(*,*) stressIn(i), stressSvecP(i)
      enddo

   !
      write(*,*)
      write(*,'('' --> TOTAL NUMBER OF STATE VARIABLES IS'',I6)') it%nH
      write(*,*)   
   !  pause
   !
   endif

      allocate(hVec(it%nH))
      allocate(hVecDot(it%nH))
   !
      hvec=it%hVecInit

      f_eval=evalEvolution_vpsc(it,tK,stressSvecP,hVec, &
         & L_o,m,g,hVecDot,phaseResponseList)
   !
      if (f_eval /= 0) then
         write(*,*) 'PROBLEM IN evalEvolution_vpsc'
         stop
      endif
   !
   if (diagnostics == 1) then 
      write(*,*)
      write(*,*) 'L_o='
      write(*,*)  L_o
      write(*,*)
      do i=1,dims
      do j=1,i
         write(*,*) (L_o(i,j)+L_o(j,i))*0.5
         enddo
         enddo
   
   !
      call cpu_time(tEnd)
      write(*,*) tEnd-tStart
   endif

    ! old rotation
   !strainOut(1) = L_o(1,1)
   !strainOut(2) = (L_o(1,2)+L_o(2,1))*0.5
   !strainOut(3) = L_o(2,2)
   !strainOut(4) = (L_o(1,3)+L_o(3,1))*0.5
   !strainOut(5) = (L_o(2,3)+L_o(3,2))*0.5
   !strainOut(6) = L_o(3,3)

    ! new rotation
   strainOut(1) = L_o(1,1)
   strainOut(2) = L_o(2,2)
   strainOut(3) = L_o(3,3)
   strainOut(4) = (L_o(1,2)+L_o(2,1))*0.5
   strainOut(5) = (L_o(2,3)+L_o(3,2))*0.5
   strainOut(6) = (L_o(1,3)+L_o(3,1))*0.5

END subroutine vpsc_run

