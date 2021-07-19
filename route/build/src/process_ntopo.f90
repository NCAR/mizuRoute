module process_ntopo

! data types
USE nrtype,    ONLY: i4b,dp,lgt          ! variable types, etc.
USE nrtype,    ONLY: strLen              ! length of characters
USE dataTypes, ONLY: var_ilength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength,dlength ! double precision type: var(:)%dat, or dat

! global vars
USE public_var, ONLY: iulog              ! i/o logical unit number
USE public_var, ONLY: idSegOut           ! ID for stream segment at the bottom of the subset

! options
USE public_var, ONLY: topoNetworkOption  ! option to compute network topology
USE public_var, ONLY: computeReachList   ! option to compute reach list
USE public_var, ONLY: hydGeometryOption  ! option to obtain routing parameters
USE public_var, ONLY: routOpt            ! option for desired routing method
USE public_var, ONLY: allRoutingMethods  ! option for routing methods - all the methods
USE public_var, ONLY: kinematicWave      ! option for routing methods - Lagrangian kinematic wave only
USE public_var, ONLY: kinematicWaveEuler ! option for routing methods - Euler kinematic wave only
USE public_var, ONLY: impulseResponseFunc! option for routing methods - IRF only

! named variables
USE public_var, ONLY: true,false         ! named integers for true/false

! named variables
USE var_lookup, ONLY: ixSEG              ! index of variables for the stream segments
USE var_lookup, ONLY: ixNTOPO            ! index of variables for the network topology
USE var_lookup, ONLY: ixPFAF             ! index of variables for the pfafstetter code

! common variables
USE public_var, ONLY: compute            ! compute given variable
USE public_var, ONLY: doNotCompute       ! do not compute given variable
USE public_var, ONLY: readFromFile       ! read given variable from a file
USE public_var, ONLY: realMissing        ! missing value for real
USE public_var, ONLY: integerMissing     ! missing value for integers

implicit none

! privacy -- everything private unless declared explicitly
private
public::check_river_properties
public::augment_ntopo
public::put_data_struct

contains

 ! *********************************************************************
 ! public subroutine: augment river network data
 ! *********************************************************************
 subroutine augment_ntopo(&
                  ! input: model control
                  nHRU,             & ! number of HRUs
                  nSeg,             & ! number of stream segments
                  ! inout: populate data structures
                  structHRU,        & ! ancillary data for HRUs
                  structSEG,        & ! ancillary data for stream segments
                  structHRU2seg,    & ! ancillary data for mapping hru2basin
                  structNTOPO,      & ! ancillary data for network toopology
                  ! output:
                  ierr, message,    & ! error control
                  ! optional output:
                  tot_hru,          & ! total number of all the upstream hrus for all stream segments
                  tot_upseg,        & ! total number of immediate upstream segments for all  stream segments
                  tot_upstream,     & ! total number of all the upstream segments for all stream segments
                  tot_uh,           & ! total number of unit hydrograph from all the stream segments
                  ixHRU_desired,    & ! indices of desired hrus
                  ixSeg_desired     ) ! indices of desired reaches

 ! external subroutines/data
 ! network topology routine
 USE network_topo,     only:hru2segment           ! get the mapping between HRUs and segments
 USE network_topo,     only:up2downSegment        ! get the mapping between upstream and downstream segments
 USE network_topo,     only:reachOrder            ! define the processing order
 USE network_topo,     only:streamOrdering        ! define stream order (Strahler)
 USE network_topo,     only:reach_list            ! reach list
 USE network_topo,     only:reach_mask            ! identify all reaches upstream of a given reach
 ! Routing parameter estimation routine
 USE routing_param,    only:make_uh               ! construct reach unit hydrograph
 ! routing spatial constant parameters
 USE globalData,       only:mann_n, wscale        ! KWT routing parameters (Transfer function parameters)
 USE globalData,       only:velo, diff            ! IRF routing parameters (Transfer function parameters)

 USE public_var, ONLY: dt                        ! simulation time step [sec]

 ! This subroutine populate river network topology data strucutres
 implicit none
 ! output: model control
 integer(i4b),       intent(in)                    :: nHRU             ! number of HRUs
 integer(i4b),       intent(in)                    :: nSeg             ! number of stream segments
 ! inout: populate data structures
 type(var_dlength), intent(inout), allocatable     :: structHRU(:)     ! HRU properties
 type(var_dlength), intent(inout), allocatable     :: structSEG(:)     ! stream segment properties
 type(var_ilength), intent(inout), allocatable     :: structHRU2seg(:) ! HRU-to-segment mapping
 type(var_ilength), intent(inout), allocatable     :: structNTOPO(:)   ! network topology
 ! output: error control
 integer(i4b)      , intent(out)                   :: ierr             ! error code
 character(*)      , intent(out)                   :: message          ! error message
 ! optional output:
 integer(i4b), optional, intent(out)               :: tot_upstream     ! total number of all of the upstream stream segments for all stream segments
 integer(i4b), optional, intent(out)               :: tot_upseg        ! total number of immediate upstream segments for all  stream segments
 integer(i4b), optional, intent(out)               :: tot_hru          ! total number of all the upstream hrus for all stream segments
 integer(i4b), optional, intent(out)               :: tot_uh           ! total number of unit hydrograph from all the stream segments
 integer(i4b), optional, intent(out),  allocatable :: ixHRU_desired(:) ! indices of desired hrus
 integer(i4b), optional, intent(out),  allocatable :: ixSeg_desired(:) ! indices of desired reaches
 ! --------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=strLen)                             :: cmessage             ! error message of downwind routine
 integer(i4b)                                      :: tot_upstream_tmp     ! temporal storage for tot_upstream
 integer(i4b)                                      :: tot_upseg_tmp        ! temporal storage tot_upseg_tmp
 integer(i4b)                                      :: tot_hru_tmp          ! temporal storage tot_hru_tmp
 integer(i4b)                                      :: tot_uh_tmp           ! temporal storage tot_uh_tmp
 integer(i4b), allocatable                         :: ixHRU_desired_tmp(:) ! temporal storage ixHRU_desired_tmp
 integer(i4b), allocatable                         :: ixSeg_desired_tmp(:) ! temporal storage ixSeg_desired_tmp
 integer(i4b)                                      :: iSeg                 ! indices for stream segment
 real(dp)     , allocatable                        :: seg_length(:)        ! temporal array for segment length
 type(dlength), allocatable                        :: temp_dat(:)          ! temporal storage for dlength data structure
 integer(i4b), parameter                           :: maxUpstreamFile=10000000 ! 10 million: maximum number of upstream reaches to enable writing
 integer*8                                         :: time0,time1,cr       ! for timing

 ierr=0; message='augment_ntopo/'

 call system_clock(count_rate=cr)
 call system_clock(time0)

 ! ---------- get the mapping between HRUs and segments ------------------------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! get the mapping between HRUs and basins
  call hru2segment(&
                   ! input
                   nHRU,          & ! input: number of HRUs
                   nSeg,          & ! input: number of stream segments
                   ! input-output: data structures
                   structHRU,     & ! ancillary data for HRUs
                   structSEG,     & ! ancillary data for stream segments
                   structHRU2seg, & ! ancillary data for mapping hru2basin
                   structNTOPO,   & ! ancillary data for network toopology
                   ! output
                   tot_hru_tmp,   & ! output: total number of all the upstream hrus for all stream segments
                   ierr, cmessage)  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,1PG15.7,A)') 'after hru2segment: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 endif  ! if need to compute network topology

 ! ---------- get the mapping between upstream and downstream segments ---------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! get the mapping between upstream and downstream segments
  call up2downSegment(&
                      ! input
                      nSeg,          & ! input: number of stream segments
                      ! input-output: data structures
                      structNTOPO,   & ! ancillary data for network toopology
                      ! output
                      tot_upseg_tmp, & ! output: sum of immediate upstream segments
                      ierr, cmessage)  ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,1PG15.7,A)') 'after up2downSegment: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 endif  ! if need to compute network topology

 ! ---------- get the processing order -----------------------------------------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! defines the processing order for the individual stream segments in the river network
  call REACHORDER(nSeg,         &   ! input:        number of reaches
                  structNTOPO,  &   ! input:output: network topology
                  ierr, cmessage)   ! output:       error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,1PG15.7,A)') 'after reachOrder: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 endif  ! if need to compute network topology

 ! ---------- compute Strahler stream order for each reach ---------------------------------------------------

 ! check the need to compute network topology
 if(topoNetworkOption==compute)then

  ! defines the processing order for the individual stream segments in the river network
  call streamOrdering(nSeg,         &   ! input:        number of reaches
                      structNTOPO,  &   ! input:output: network topology
                      ierr, cmessage)   ! output:       error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,1PG15.7,A)') 'after streamOrder: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 endif  ! if need to compute network topology

 ! ---------- get the list of all upstream reaches above a given reach ---------------------------------------

 ! get the list of all upstream reaches above a given reach
 call reach_list(&
                 ! input
                 nSeg,                        & ! Number of reaches
                 (computeReachList==compute), & ! flag to compute the reach list
                 structNTOPO,                 & ! Network topology
                 ! output
                 structSEG,                   & ! input: ancillary data for stream segments
                 tot_upstream_tmp,            & ! Total number of upstream reaches for all reaches
                 ierr, cmessage)                ! Error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 !write(*,'(a,1x,1PG15.7,A)') 'after reach_list: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 ! ---------- Compute routing parameters  --------------------------------------------------------------------

 ! compute hydraulic geometry (width and Manning's "n")
 if(hydGeometryOption==compute)then

  ! (hydraulic geometry only needed for the kinematic wave method)
  if (routOpt==allRoutingMethods .or. routOpt==kinematicWave .or. routOpt==kinematicWaveEuler) then
   do iSeg=1,nSeg
    structSEG(iSeg)%var(ixSEG%width)%dat(1) = wscale * sqrt(structSEG(iSeg)%var(ixSEG%totalArea)%dat(1))  ! channel width (m)
    structSEG(iSeg)%var(ixSEG%man_n)%dat(1) = mann_n                                                      ! Manning's "n" paramater (unitless)
   end do
  end if

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,1PG15.7,A)') 'after river geometry: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 endif  ! computing hydraulic geometry

 ! get the channel unit hydrograph
 if(topoNetworkOption==compute)then

  ! (channel unit hydrograph is only needed for the impulse response function)
  if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then

   ! extract the length information from the structure and place it in a vector
   allocate(seg_length(nSeg), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': seg_length'; return; endif
   do iSeg = 1,nSeg
     seg_length(iSeg) = structSEG(iSeg)%var(ixSEG%length)%dat(1)
   end do

   ! compute lag times in the channel unit hydrograph
   call make_uh(seg_length, dt, velo, diff, temp_dat, ierr, cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! put the lag times in the data structures
   tot_uh_tmp = 0
   do iSeg=1,nSeg
    allocate(structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat(size(temp_dat(iSeg)%dat)), stat=ierr, errmsg=cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage)//': structSEG%var(ixSEG%uh)%dat'; return; endif
    structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat(:) = temp_dat(iSeg)%dat(:)
    tot_uh_tmp = tot_uh_tmp+size(temp_dat(iSeg)%dat)
   enddo

  endif ! if using the impulse response function

  ! get timing
  call system_clock(time1)
  !write(*,'(a,1x,1PG15.7,A)') 'after reach parameters: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 endif ! if there is a need to compute the channel unit hydrograph

 ! ---------- get the mask of all upstream reaches above a given reach ---------------------------------------

 ! get the mask of all upstream reaches above a given reach
 call reach_mask(&
                 ! input
                 idSegOut,          &  ! input: reach index
                 structNTOPO,       &  ! input: network topology structures
                 structSeg,         &  ! input: river reach properties
                 nHRU,              &  ! input: number of HRUs
                 nSeg,              &  ! input: number of reaches
                 ! output: updated dimensions
                 tot_hru_tmp,       &  ! input+output: total number of all the upstream hrus for all stream segments
                 tot_upseg_tmp,     &  ! input+output: sum of immediate upstream segments
                 tot_upstream_tmp,  &  ! input+output: total number of upstream reaches for all reaches
                 tot_uh_tmp,        &  ! input+output: total number of unit hydrograph dimensions
                 ! output: dimension masks
                 ixHRU_desired_tmp, &  ! output: indices of desired hrus
                 ixSeg_desired_tmp, &  ! output: indices of desired reaches
                 ! output: error control
                 ierr, cmessage )  ! output: error control
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get timing
 call system_clock(time1)
 !write(*,'(a,1x,1PG15.7,A)') 'after reach_mask: time = ', real(time1-time0,kind(dp))/real(cr), ' s'

 ! for optional output
 if (present(tot_hru))       tot_hru=tot_hru_tmp
 if (present(tot_upseg))     tot_upseg=tot_upseg_tmp
 if (present(tot_upstream))  tot_upstream=tot_upstream_tmp
 if (present(tot_uh))        tot_uh=tot_uh_tmp
 if (present(ixSeg_desired)) then
   allocate(ixSeg_desired(size(ixSeg_desired_tmp)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem in allocating [isSeg_desire]'; return; endif
   ixSeg_desired=ixSeg_desired_tmp
 endif
 if (present(ixHRU_desired)) then
   allocate(ixHRU_desired(size(ixHRU_desired_tmp)), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'problem in allocating [ixHRU_desired]'; return; endif
   ixHRU_desired=ixHRU_desired_tmp
 endif

end subroutine augment_ntopo

 ! *********************************************************************
 ! public subroutine: check network data is physically valid
 ! *********************************************************************
 subroutine check_river_properties(structNTOPO, structHRU, structSEG, &  ! input: data structure for physical river network data
                                   ierr, message)
  ! saved global data
  USE public_var, ONLY: min_slope          ! minimum slope
  implicit none
  ! input
  type(var_ilength)           , intent(in)       :: structNTOPO(:)   ! network topology
  type(var_dlength)           , intent(inout)    :: structHRU(:)     ! HRU properties
  type(var_dlength)           , intent(inout)    :: structSEG(:)     ! stream reach properties
  ! output: error control
  integer(i4b)                , intent(out)      :: ierr             ! error code
  character(*)                , intent(out)      :: message          ! error message
  ! local varialbles
  integer(i4b)                                   :: nSeg,nHru        ! number of stream reaches and HRUs
  integer(i4b)                                   :: iSeg,iHru        ! loop indices

  ierr=0; message='check_river_properties/'

  nSeg = size(structSEG)
  nHru = size(structHRU)

  ! check reach
  do iSeg = 1,nSeg
    associate(segId => structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1))
    ! Check reach length
    if (structSEG(iSeg)%var(ixSEG%length)%dat(1) <= 0 )then
     ierr=10
     write(message,'(a,i0,a,1PG15.7)') trim(message)//'reach length for reach id ', segId, ' is negative:', structSEG(iSeg)%var(ixSEG%length)%dat(1)
     return
    endif
    end associate
  enddo

  do iHru = 1,nHru
  ! check somehting for hru properties
  enddo

  end subroutine check_river_properties

 ! *********************************************************************
 ! public subroutine: populate old data strucutures
 ! *********************************************************************
 ! ---------- temporary code: populate old data structures --------------------------------------------------
 subroutine put_data_struct(nSeg, structSEG, structNTOPO, &
                            RPARAM_in, NETOPO_in , ierr, message)
  ! saved global data
  USE dataTypes,     ONLY: RCHPRP             ! Reach parameters
  USE dataTypes,     ONLY: RCHTOPO            ! Network topology
  USE globalData,    ONLY: fshape, tscale     ! basin IRF routing parameters (Transfer function parameters)
  USE public_var,    ONLY: min_slope          ! minimum slope
  USE public_var,    ONLY: dt                 ! simulation time step [sec]
  USE public_var,    ONLY: is_lake_sim        ! lake simulation option
  ! external subroutines
  USE routing_param, ONLY: basinUH            ! construct basin unit hydrograph
  implicit none
  ! input
  integer(i4b)                , intent(in)       :: nSeg             ! number of stream segments
  type(var_dlength)           , intent(in)       :: structSEG(:)     ! stream segment properties
  type(var_ilength)           , intent(in)       :: structNTOPO(:)   ! network topology
  ! output data structure
  type(RCHPRP)  , allocatable , intent(out)      :: RPARAM_in(:)     ! Reach Parameters
  type(RCHTOPO) , allocatable , intent(out)      :: NETOPO_in(:)     ! River Network topology
  ! output: error control
  integer(i4b)                , intent(out)      :: ierr             ! error code
  character(*)                , intent(out)      :: message          ! error message
  ! local varialbles
  character(len=strLen)                          :: cmessage         ! error message of downwind routine
  integer(i4b)                                   :: nUps             ! number of upstream segments for a segment
  integer(i4b)                                   :: iSeg,iUps        ! loop indices

  ierr=0; message='put_data_struct/'

  ! get lag times in the basin unit hydrograph (not sure this is right place...)
  call basinUH(dt, fshape, tscale, ierr, cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! allocate space
  allocate(RPARAM_in(nSeg), NETOPO_in(nSeg), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'unable to allocate space for old data structures'; return; endif

  ! loop through stream segments
  do iSeg=1,nSeg

   ! ----- reach parameters -----

   ! copy data into the reach parameter structure
   RPARAM_in(iSeg)%RLENGTH         =     structSEG(iSeg)%var(ixSEG%length)%dat(1)
   RPARAM_in(iSeg)%R_SLOPE         = max(structSEG(iSeg)%var(ixSEG%slope)%dat(1), min_slope)
   RPARAM_in(iSeg)%R_MAN_N         =     structSEG(iSeg)%var(ixSEG%man_n)%dat(1)
   RPARAM_in(iSeg)%R_WIDTH         =     structSEG(iSeg)%var(ixSEG%width)%dat(1)

   if (is_lake_sim) then
     RPARAM_in(iSeg)%D03_MaxStorage  =     structSEG(iSeg)%var(ixSEG%D03_MaxStorage)%dat(1)
     RPARAM_in(iSeg)%D03_Coefficient =     structSEG(iSeg)%var(ixSEG%D03_Coefficient)%dat(1)
     RPARAM_in(iSeg)%D03_Power       =     structSEG(iSeg)%var(ixSEG%D03_Power)%dat(1)
     RPARAM_in(iSeg)%H06_Smax        =     structSEG(iSeg)%var(ixSEG%H06_Smax)%dat(1)
     RPARAM_in(iSeg)%H06_alpha       =     structSEG(iSeg)%var(ixSEG%H06_alpha)%dat(1)
     RPARAM_in(iSeg)%H06_envfact     =     structSEG(iSeg)%var(ixSEG%H06_envfact)%dat(1)
     RPARAM_in(iSeg)%H06_S_ini       =     structSEG(iSeg)%var(ixSEG%H06_S_ini)%dat(1)
     RPARAM_in(iSeg)%H06_c1          =     structSEG(iSeg)%var(ixSEG%H06_c1)%dat(1)
     RPARAM_in(iSeg)%H06_c2          =     structSEG(iSeg)%var(ixSEG%H06_c2)%dat(1)
     RPARAM_in(iSeg)%H06_exponent    =     structSEG(iSeg)%var(ixSEG%H06_exponent)%dat(1)
     RPARAM_in(iSeg)%H06_denominator =     structSEG(iSeg)%var(ixSEG%H06_denominator)%dat(1)
     RPARAM_in(iSeg)%H06_c_compare   =     structSEG(iSeg)%var(ixSEG%H06_c_compare)%dat(1)
     RPARAM_in(iSeg)%H06_frac_Sdead  =     structSEG(iSeg)%var(ixSEG%H06_frac_Sdead)%dat(1)
     RPARAM_in(iSeg)%H06_E_rel_ini   =     structSEG(iSeg)%var(ixSEG%H06_E_rel_ini)%dat(1)
     RPARAM_in(iSeg)%H06_I_Jan       =     structSEG(iSeg)%var(ixSEG%H06_I_Jan)%dat(1)
     RPARAM_in(iSeg)%H06_I_Feb       =     structSEG(iSeg)%var(ixSEG%H06_I_Feb)%dat(1)
     RPARAM_in(iSeg)%H06_I_Mar       =     structSEG(iSeg)%var(ixSEG%H06_I_Mar)%dat(1)
     RPARAM_in(iSeg)%H06_I_Apr       =     structSEG(iSeg)%var(ixSEG%H06_I_Apr)%dat(1)
     RPARAM_in(iSeg)%H06_I_May       =     structSEG(iSeg)%var(ixSEG%H06_I_May)%dat(1)
     RPARAM_in(iSeg)%H06_I_Jun       =     structSEG(iSeg)%var(ixSEG%H06_I_Jun)%dat(1)
     RPARAM_in(iSeg)%H06_I_Jul       =     structSEG(iSeg)%var(ixSEG%H06_I_Jul)%dat(1)
     RPARAM_in(iSeg)%H06_I_Aug       =     structSEG(iSeg)%var(ixSEG%H06_I_Aug)%dat(1)
     RPARAM_in(iSeg)%H06_I_Sep       =     structSEG(iSeg)%var(ixSEG%H06_I_Sep)%dat(1)
     RPARAM_in(iSeg)%H06_I_Oct       =     structSEG(iSeg)%var(ixSEG%H06_I_Oct)%dat(1)
     RPARAM_in(iSeg)%H06_I_Nov       =     structSEG(iSeg)%var(ixSEG%H06_I_Nov)%dat(1)
     RPARAM_in(iSeg)%H06_I_Dec       =     structSEG(iSeg)%var(ixSEG%H06_I_Dec)%dat(1)
     RPARAM_in(iSeg)%H06_D_Jan       =     structSEG(iSeg)%var(ixSEG%H06_D_Jan)%dat(1)
     RPARAM_in(iSeg)%H06_D_Feb       =     structSEG(iSeg)%var(ixSEG%H06_D_Feb)%dat(1)
     RPARAM_in(iSeg)%H06_D_Mar       =     structSEG(iSeg)%var(ixSEG%H06_D_Mar)%dat(1)
     RPARAM_in(iSeg)%H06_D_Apr       =     structSEG(iSeg)%var(ixSEG%H06_D_Apr)%dat(1)
     RPARAM_in(iSeg)%H06_D_May       =     structSEG(iSeg)%var(ixSEG%H06_D_May)%dat(1)
     RPARAM_in(iSeg)%H06_D_Jun       =     structSEG(iSeg)%var(ixSEG%H06_D_Jun)%dat(1)
     RPARAM_in(iSeg)%H06_D_Jul       =     structSEG(iSeg)%var(ixSEG%H06_D_Jul)%dat(1)
     RPARAM_in(iSeg)%H06_D_Aug       =     structSEG(iSeg)%var(ixSEG%H06_D_Aug)%dat(1)
     RPARAM_in(iSeg)%H06_D_Sep       =     structSEG(iSeg)%var(ixSEG%H06_D_Sep)%dat(1)
     RPARAM_in(iSeg)%H06_D_Oct       =     structSEG(iSeg)%var(ixSEG%H06_D_Oct)%dat(1)
     RPARAM_in(iSeg)%H06_D_Nov       =     structSEG(iSeg)%var(ixSEG%H06_D_Nov)%dat(1)
     RPARAM_in(iSeg)%H06_D_Dec       =     structSEG(iSeg)%var(ixSEG%H06_D_Dec)%dat(1)
     RPARAM_in(iSeg)%H06_purpose     =     structSEG(iSeg)%var(ixSEG%H06_purpose )%dat(1)
     RPARAM_in(iSeg)%H06_I_mem_F     =     (structSEG(iSeg)%var(ixSEG%H06_I_mem_F)%dat(1)==1)
     RPARAM_in(iSeg)%H06_D_mem_F     =     (structSEG(iSeg)%var(ixSEG%H06_D_mem_F)%dat(1)==1)
     RPARAM_in(iSeg)%H06_I_mem_L     =     structSEG(iSeg)%var(ixSEG%H06_I_mem_L )%dat(1)
     RPARAM_in(iSeg)%H06_D_mem_L     =     structSEG(iSeg)%var(ixSEG%H06_D_mem_L )%dat(1)
   end if

   ! compute variables
   RPARAM_in(iSeg)%BASAREA = structSEG(iSeg)%var(ixSEG%basArea)%dat(1)
   RPARAM_in(iSeg)%UPSAREA = structSEG(iSeg)%var(ixSEG%upsArea)%dat(1)
   RPARAM_in(iSeg)%TOTAREA = structSEG(iSeg)%var(ixSEG%totalArea)%dat(1)

   ! NOT USED: MINFLOW -- minimum environmental flow
   RPARAM_in(iSeg)%MINFLOW = structSEG(iSeg)%var(ixSEG%minFlow)%dat(1)

   ! ----- network topology -----

   ! reach indices
   NETOPO_in(iSeg)%REACHIX = structNTOPO(iSeg)%var(ixNTOPO%segIndex)%dat(1)     ! reach index (1, 2, 3, ..., nSeg)
   NETOPO_in(iSeg)%REACHID = structNTOPO(iSeg)%var(ixNTOPO%segId)%dat(1)        ! reach ID (unique reach identifier)

   ! downstream reach indices
   NETOPO_in(iSeg)%DREACHI = structNTOPO(iSeg)%var(ixNTOPO%downSegIndex)%dat(1) ! Immediate Downstream reach index
   NETOPO_in(iSeg)%DREACHK = structNTOPO(iSeg)%var(ixNTOPO%downSegId)%dat(1)    ! Immediate Downstream reach ID

   ! allocate space for immediate upstream reach indices
   nUps = size(structNTOPO(iSeg)%var(ixNTOPO%upSegIds)%dat)
   allocate(NETOPO_in(iSeg)%UREACHI(nUps), NETOPO_in(iSeg)%UREACHK(nUps), NETOPO_in(iSeg)%goodBas(nUps), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for upstream structures'; return; endif

   ! populate immediate upstream data structures
   if(nUps>0)then
    do iUps=1,nUps   ! looping through upstream reaches
     NETOPO_in(iSeg)%UREACHI(iUps) = structNTOPO(iSeg)%var(ixNTOPO%upSegIndices)%dat(iUps)      ! Immediate Upstream reach indices
     NETOPO_in(iSeg)%UREACHK(iUps) = structNTOPO(iSeg)%var(ixNTOPO%upSegIds    )%dat(iUps)      ! Immediate Upstream reach Ids
     NETOPO_in(iSeg)%goodBas(iUps) = (structNTOPO(iSeg)%var(ixNTOPO%goodBasin)%dat(iUps)==true) ! "good" basin
    end do  ! Loop through upstream reaches
   endif

   ! define the reach order
   NETOPO_in(iSeg)%RHORDER = structNTOPO(iSeg)%var(ixNTOPO%rchOrder)%dat(1)  ! Processing sequence

   ! allocate space for contributing HRUs
   nUps = structNTOPO(iSeg)%var(ixNTOPO%nHRU)%dat(1)
   allocate(NETOPO_in(iSeg)%HRUID(nUps), NETOPO_in(iSeg)%HRUIX(nUps), NETOPO_in(iSeg)%HRUWGT(nUps), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for contributing HRUs'; return; endif

   ! HRU2SEG topology
   if(nUps>0)then
     do iUps = 1, nUps
       NETOPO_in(iSeg)%HRUID(iUps) = structNTOPO(iSeg)%var(ixNTOPO%hruContribId)%dat(iUps)
       NETOPO_in(iSeg)%HRUIX(iUps) = structNTOPO(iSeg)%var(ixNTOPO%hruContribIx)%dat(iUps)
       NETOPO_in(iSeg)%HRUWGT(iUps) = structSEG(iSeg)%var(ixSEG%weight)%dat(iUps)
     end do  ! Loop through contributing HRU loop
   end if

   ! lake parameters
   if (is_lake_sim) then
     NETOPO_in(iSeg)%isLake       = (structNTOPO(iSeg)%var(ixNTOPO%isLake)%dat(1)==true)
     NETOPO_in(iSeg)%LakeTargVol  = (structNTOPO(iSeg)%var(ixNTOPO%LakeTargVol)%dat(1)==true)
     NETOPO_in(iSeg)%LakeModelType= structNTOPO(iSeg)%var(ixNTOPO%LakeModelType)%dat(1) ! type of the parameteric lake
     ! NOT USED: lake parameters
     NETOPO_in(iSeg)%LAKE_IX = integerMissing  ! Lake index (0,1,2,...,nlak-1)
     NETOPO_in(iSeg)%LAKE_ID = integerMissing  ! Lake ID (REC code?)
     NETOPO_in(iSeg)%BASULAK = realMissing     ! Area of basin under lake
     NETOPO_in(iSeg)%RCHULAK = realMissing     ! Length of reach under lake
     NETOPO_in(iSeg)%LAKINLT = .false.         ! .TRUE. if reach is lake inlet, .FALSE. otherwise
     NETOPO_in(iSeg)%USRTAKE = .false.         ! .TRUE. if user takes from reach, .FALSE. otherwise
   end if

   ! NOT USED: Location (available in the input files)
   NETOPO_in(iSeg)%RCHLAT1 = realMissing     ! Start latitude
   NETOPO_in(iSeg)%RCHLAT2 = realMissing     ! End latitude
   NETOPO_in(iSeg)%RCHLON1 = realMissing     ! Start longitude
   NETOPO_in(iSeg)%RCHLON2 = realMissing     ! End longitude

   ! reach unit hydrograph
   if (routOpt==allRoutingMethods .or. routOpt==impulseResponseFunc) then
     allocate(NETOPO_in(iSeg)%UH(size(structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat)), stat=ierr, errmsg=cmessage)
     if(ierr/=0)then; message=trim(message)//trim(cmessage)//': NETOPO_in(iSeg)%UH'; return; endif
     NETOPO_in(iSeg)%UH(:) =  structSEG(iSeg)%var(ixSEG%timeDelayHist)%dat(:)
   end if

   ! upstream reach list
   allocate(NETOPO_in(iSeg)%RCHLIST(size(structNTOPO(iSeg)%var(ixNTOPO%allUpSegIndices)%dat)), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//': NETOPO_in(iSeg)%RCHLIST'; return; endif
   NETOPO_in(iSeg)%RCHLIST(:) =  structNTOPO(iSeg)%var(ixNTOPO%allUpSegIndices)%dat(:)

  end do  ! looping through stream segments

 end subroutine put_data_struct

end module process_ntopo
