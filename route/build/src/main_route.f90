module main_route_module

! variable types
USE nrtype                                                  ! variable types, etc.

! mapping HRU runoff to reach
USE remapping,           only : basin2reach
! subroutines: basin routing
USE basinUH_module,      only : IRF_route_basin             ! perform UH convolution for basin routing

! subroutines: river routing
USE accum_runoff_module, only : accum_runoff                ! upstream flow accumulation
USE kwt_route_module,    only : kwt_route                  ! kinematic wave routing method
USE irf_route_module,    only : irf_route                  ! unit hydrograph (impulse response function) routing method
!USE kwt_route_module,    only : kwt_route => kwt_route_orig ! kinematic wave routing method
!USE irf_route_module,    only : irf_route => irf_route_orig ! river network unit hydrograph method

! subroutines: general utility
USE nr_utility_module,   only : findIndex                   ! find index within a vector
USE nr_utility_module,   only : arth

implicit none

private

public::main_route

contains

 ! ******
 ! public subroutine: main HRU/reach routing routines
 ! ************************
 subroutine main_route(iens,           &  ! input:  ensemble index
                       ierr, message)     ! output: error control
   ! Details:
   ! Given HRU (basin) runoff, perform hru routing (optional) to get reach runoff, and then channel routing ! Restriction:
   ! 1. Reach order in NETOPO, RPARAM, RCHFLX, KROUTE must be in the same orders
   ! 2. Process a list of reach indices (in terms of NETOPO etc.) given by ixRchProcessed
   ! 3. basinRunoff_in is given in the order of NETOPO(:)%HRUIX.

   ! shared data
   USE public_var, only : routOpt
   USE public_var, only : doesBasinRoute
   USE public_var, only : doesAccumRunoff
   USE public_var, only : allRoutingMethods
   USE public_var, only : kinematicWave
   USE public_var, only : impulseResponseFunc
   USE globalData, only : TSEC                    ! beginning/ending of simulation time step [sec]
   USE globalData, only : ixPrint                 ! desired reach index to be on-screen print

   USE globalData, only : NETOPO           ! entire river reach netowrk topology structure
   USE globalData, only : RPARAM           ! entire river reach parameter structure
   USE globalData, only : RCHFLX           ! entire reach flux structure
   USE globalData, only : KROUTE           ! entire river reach kwt sate structure
   USE globalData, only : runoff_data      ! runoff data structure
   USE globalData, only : river_basin      ! OMP basin decomposition
   USE globalData, only : nRch             ! number of reaches in the whoel river network

   implicit none

   ! input
   integer(i4b),               intent(in)    :: iens                 ! ensemble member
   ! output
   integer(i4b),               intent(out)   :: ierr                 ! error code
   character(len=strLen),      intent(out)   :: message              ! error message
   ! local variables
   character(len=strLen)                     :: cmessage             ! error message of downwind routine
   real(dp)                                  :: T0,T1                ! beginning/ending of simulation time step [sec]
   real(dp),      allocatable                :: reachRunoff_local(:) ! reach runoff (m/s)
   integer(i4b)                              :: mainstem=2           ! basin indicator (mainstem = 2 )
   integer(i4b)                              :: iSeg                 ! index of reach
   ! timing variables
   integer*8                                 :: cr, startTime, endTime ! rate, start and end time stamps
   real(dp)                                  :: elapsedTime            ! elapsed time for the process

  CALL system_clock(count_rate=cr)

  ! initialize errors
  ierr=0; message = "main_routing/"

  ! define the start and end of the time step
  T0=TSEC(0); T1=TSEC(1)

  allocate(reachRunoff_local(nRch), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating arrays for [reachRunoff_local]'; return; endif

  ! 1. subroutine: map basin runoff to river network HRUs
  ! map the basin runoff to the stream network...
  call system_clock(startTime)
  call basin2reach(&
                  ! input
                  runoff_data%basinRunoff, & ! basin runoff (m/s)
                  NETOPO,                  & ! reach topology
                  RPARAM,                  & ! reach parameter
                  ! output
                  reachRunoff_local,       & ! intent(out): reach runoff (m3/s)
                  ierr, cmessage)            ! intent(out): error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  call system_clock(endTime)
  elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
  write(*,"(A,1PG15.7,A)") '      elapsed-time [basin2reach] = ', elapsedTime, ' s'

  ! 2. subroutine: basin route
  if (doesBasinRoute == 1) then
    call system_clock(startTime)
    ! instantaneous runoff volume (m3/s) to data structure
    do iSeg = 1,nRch
     RCHFLX(iens,iSeg)%BASIN_QI = reachRunoff_local(iSeg)
    enddo
    ! perform Basin routing
    call IRF_route_basin(iens,              &  ! input:  ensemble index
                         RCHFLX,            &  ! inout:  reach flux data structure
                         ierr, cmessage)       ! output: error controls
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [IRF_route_basin] = ', elapsedTime, ' s'
  else
    ! no basin routing required (handled outside mizuRoute))
    do iSeg = 1,nRch
    RCHFLX(iens,iSeg)%BASIN_QR(0) = RCHFLX(iens,iSeg)%BASIN_QR(1)   ! streamflow from previous step
    RCHFLX(iens,iSeg)%BASIN_QR(1) = reachRunoff_local(iSeg)         ! streamflow (m3/s)
    end do
  end if

  ! 3. subroutine: river reach routing
   ! perform upstream flow accumulation
   if (doesAccumRunoff == 1) then
     call system_clock(startTime)
     call accum_runoff(iens,              &  ! input: ensemble index
                       ixPrint,           &  ! input: index of verbose reach
                       NETOPO,            &  ! input: reach topology data structure
                       RCHFLX,            &  ! inout: reach flux data structure
                       ierr, cmessage)       ! output: error controls
     if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
     call system_clock(endTime)
     elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
     write(*,"(A,1PG15.7,A)") '      elapsed-time [accum_runoff] = ', elapsedTime, ' s'
   endif

   ! perform KWT routing
   if (routOpt==allRoutingMethods .or. routOpt==kinematicWave) then
    call system_clock(startTime)
    call kwt_route(iens,                 & ! input: ensemble index
                   river_basin,          & ! input: river basin data type
                   T0,T1,                & ! input: start and end of the time step
                   mainstem,             & ! input: basinType (1-> tributary, 2->mainstem)
                   ixPrint,              & ! input: index of the desired reach
                   NETOPO,               & ! input: reach topology data structure
                   RPARAM,               & ! input: reach parameter data structure
                   KROUTE,               & ! inout: reach state data structure
                   RCHFLX,               & ! inout: reach flux data structure
                   ierr,cmessage)          ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [kwt_route] = ', elapsedTime, ' s'
   endif

   ! perform IRF routing
   if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
    call system_clock(startTime)
    call irf_route(iens,                & ! input: ensemble index
                   river_basin,         & ! input: river basin data type
                   ixPrint,             & ! input: index of the desired reach
                   NETOPO,              & ! input: reach topology data structure
                   RCHFLX,              & ! inout: reach flux data structure
                   ierr,cmessage)         ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    call system_clock(endTime)
    elapsedTime = real(endTime-startTime, kind(dp))/real(cr)
    write(*,"(A,1PG15.7,A)") '      elapsed-time [irf_route] = ', elapsedTime, ' s'
   endif

 end subroutine main_route

end module main_route_module
