PROGRAM caller
   !
   USE vpETxt_interface_base
   USE vpscVpETxt_interface
   !
   IMPLICIT NONE
   !
   !!  integer :: init_vpsc, evalEvolution_vpsc
   !
   ! init_vpsc and evalEvolution_vpsc arguments
   !
   integer :: status
   type (vpETxt_interface_type) :: it
   !
   ! init_vpsc arguments
   !
   character(len=20) :: fnameIn
   type (vpETxt_phase_data_type),allocatable :: phaseDataList(:)
   integer :: MPICommF
   !
   ! evalEvolution_vpsc arguments
   !
   type (vpETxt_phase_response_type),allocatable :: phaseResponseList(:)
   real(idp) tk,stressSvecP(SVEC+1), L_o(DIMS,DIMS),m,g
   real(idp) D(DIMS,DIMS)
   real(idp), pointer :: hVec(:), hVecDot(:)
   real(idp) stressIn(6)
   !
   ! Outer iteration loop variables
   !
   real :: tStart, tEnd
   integer :: numCases, iCase, iLine
   integer :: numPhases, iPhase
   integer, allocatable :: nGrains(:)
   !
   integer :: kgr,i,j,f_eval
   real(idp) :: scaling
   integer :: ng1, ng2 ! number of grains in phases 
   !
   !
   !  List of stress inputs from tracefiles
   open(60, file='../results/stressIn.in', status='old')
   read(60,*) numCases
   !  output file for iteration and timing
   open(61, file='timing.out', status='unknown')

   fnameIn='vpsc_as_try.in'
   open(62,FILE=fnameIn,STATUS='OLD')
   read(62,*) numPhases

   it%nPhase=numPhases
   it%dRef=1.

   if (numPhases > 2) then
      write(*,*) 'No more than 2 phases allowed right now!'
      stop
   endif

   allocate(nGrains(numPhases))
   allocate(phaseDataList(numPhases))
   allocate(phaseResponseList(numPhases))

   !nGrains(1) = 200
   !nGrains(2) = 100

   !
   !
   read(62,*) (phaseDataList(iPhase)%volfrac, iPhase=1,numPhases)
   close(62)
   !
   open(55,file='orient_ph1.in',status='old')
   open(56,file='orient_ph2.in',status='old')
   !
   read(55,*) nGrains(1)
   read(56,*) nGrains(2)

   do iPhase = 1, numPhases
      phaseDataList(iPhase)%ngrains=nGrains(iPhase)
      allocate(phaseDataList(iPhase)%orientations(3,nGrains(iPhase)))
      allocate(phaseResponseList(iPhase)%weights(nGrains(iPhase)))
      allocate(phaseResponseList(iPhase)%reorientationRates(3,nGrains(iPhase)))
      allocate(phaseResponseList(iPhase)%twinVFrates(12,nGrains(iPhase)))
   enddo

   do kgr=1,nGrains(1)
      read(55,*) (phaseDataList(1)%orientations(i,kgr),i=1,3), &
&           phaseResponseList(1)%weights(kgr)
   enddo
   close(55)
   do kgr=1,nGrains(2)
      read(56,*) (phaseDataList(2)%orientations(i,kgr),i=1,3), &
&           phaseResponseList(2)%weights(kgr)
   enddo
   close(56)
   !

   iLine = 0
      f_eval=init_vpsc(it,fnameIn,phaseDataList,MPICommF)
   !
      if (f_eval /= 0) then
         write(*,*) 'PROBLEM IN init_vpsc'
         stop
      endif
   !
   write(*,*) 'Size of interface type struct = ',sizeof(it), ' bytes'
   write(*,*) 'Size of phase data list struct = ',sizeof(phaseDataList), ' bytes'
   write(*,*) 'Size of phase response list struct = ',sizeof(phaseResponseList), ' bytes'

   do iCase = 1, numCases

      call cpu_time(tStart)

      read(60,*) (stressIn(i), i=1,6)
      iLine = iLine + 1

   !   skip lines of data
      do i=1,1000
         read(60,*) 
         iLine = iLine + 1
      enddo

   ! scale the stress
      scaling = 1.0
   ! index rotation to conform to VPSC internal notation
      stressSvecP=0.
      stressSvecP(1)=stressIn(1)*scaling
      stressSvecP(2)=stressIn(3)*scaling
      stressSvecP(3)=stressIn(6)*scaling
      stressSvecP(4)=stressIn(5)*scaling
      stressSvecP(5)=stressIn(4)*scaling
      stressSvecP(6)=stressIn(2)*scaling
      stressSvecP(7)=20.*(1./3.)*scaling

   !stressSvecP(1)=20.*(-1./3.)
   !stressSvecP(2)=20.*(-1./3.)
   !stressSvecP(3)=20.*(2./3.)
   !stressSvecP = stressSvecP*1.0
   !
      write(*,*)
      write(*,'('' --> TOTAL NUMBER OF STATE VARIABLES IS'',I6)') it%nH
      write(*,*)   
   !  pause
   !
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
      write(61,*) iCase, tEnd-tStart

   enddo ! loop over test cases
   END PROGRAM caller
