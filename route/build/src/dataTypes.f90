MODULE dataTypes

! used to create specific data types

USE nrtype,     ONLY: i4b,dp,lgt
USE nrtype,     ONLY: strLen
USE public_var, ONLY: realMissing
USE public_var, ONLY: integerMissing

implicit none

 ! everything private unless specified otherwise
 private

 ! ---------- information on the data structures -----------------------------------------------------------

 ! data structure information
 type,public :: struct_info
  character(len=strLen)  :: structName = 'empty'          ! name of the data structure
  integer(i4b)           :: varType    = integerMissing   ! type of variable in each data structure
  integer(i4b)           :: nSpace     = integerMissing   ! length of the spatial dimension in each data structure
  integer(i4b)           :: nVars      = integerMissing   ! number of variables in each data structure
 end type struct_info

 ! ---------- information on the dimensions -----------------------------------------------------------------

 ! dimension information
 type,public :: dim_info
  character(len=strLen)  :: dimName   = 'empty'          ! name of the data structure
  integer(i4b)           :: dimId     = integerMissing   ! dimension id
  integer(i4b)           :: dimLength = integerMissing   ! dim length
 end type dim_info

 ! ---------- metadata structures --------------------------------------------------------------------------

 ! define derived type for model variables, including name, description, and units
 type,public :: var_info
  character(len=strLen)  :: varName  = 'empty'         ! variable name
  character(len=strLen)  :: varDesc  = 'empty'         ! variable description
  character(len=strLen)  :: varUnit  = 'empty'         ! variable units
  integer(i4b)           :: varType  = integerMissing  ! variable type (vectors of different size)
  logical(lgt)           :: varFile  = .true.          ! .true. if the variable should be read from a file
 endtype var_info

 ! ---------- time structures ------------------------------------------------------------------------------
 type,public :: time
  integer(i4b)           :: iy       = integerMissing  ! year
  integer(i4b)           :: im       = integerMissing  ! month
  integer(i4b)           :: id       = integerMissing  ! day
  integer(i4b)           :: ih       = integerMissing  ! hour
  integer(i4b)           :: imin     = integerMissing  ! minute
  real(dp)               :: dsec     = realMissing     ! second
 endtype time

 ! ---------- states structure --------------------------------------------------------------------------
 !
 type,public :: var
  integer(i4b),  allocatable  :: array_2d_int(:,:)
  integer(i4b),  allocatable  :: array_3d_int(:,:,:)
  real(dp),      allocatable  :: array_2d_dp(:,:)
  real(dp),      allocatable  :: array_3d_dp(:,:,:)
  logical(lgt),  allocatable  :: array_2d_lgt(:,:)
  logical(lgt),  allocatable  :: array_3d_lgt(:,:,:)
 end type var

 type,public :: states
  type(var),     allocatable :: var(:)
 end type states

 ! ---------- basin data structures ----------------------------------------------------------------------
 ! segIndex points to the segment in the entire river network data
 ! segOrder is order within subset of mainstem segments or tributary segments
 type,public :: subbasin_mpi
  character(32)              :: pfaf                  ! subbasin pfaf code - mainstem starting "-"
  integer(i4b)               :: basinType             ! basin type identifier: tributary->1, mainstem->2
  integer(i4b)               :: idNode                ! node ID
  integer(i4b)               :: outletIndex           ! reach index of a domain outlet
  integer(i4b),  allocatable :: segIndex(:)           ! reach indices within a subbasin
  integer(i4b),  allocatable :: hruIndex(:)           ! hru indices within a subbasin
 end type subbasin_mpi

 ! Data structures to reach indices within each leg of stream order separately
 ! Used for openMP
 type,public :: reach
  integer(i4b), allocatable :: segIndex(:)           ! index of segment index
  integer(i4b)              :: nRch                  ! number of reach
 end type reach

 type,public :: subbasin_omp
   type(reach), allocatable :: branch(:)
 end type subbasin_omp

 ! Data structures to hold mainstem and independent tributary reaches separately
 ! Used for openMP
 type,public :: subbasin_omp_tmp
  integer(i4b)               :: outIndex             ! index of outlet segment based on segment array
  type(reach), allocatable   :: mainstem(:)          ! mainstem reach
  type(reach), allocatable   :: tributary(:)         ! tributary reach
 end type subbasin_omp_tmp

 ! ---------- general data structures ----------------------------------------------------------------------

 ! ** double precision type
 ! dat vector
 type, public :: dlength
  real(dp),allocatable                :: dat(:)    ! dat(:)
 endtype dlength
 ! var vector
 type, public :: var_dlength
  type(dlength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_dlength

 ! ** integer type
 ! dat vector
 type, public :: ilength
  integer(i4b),allocatable            :: dat(:)    ! dat(:)
 endtype ilength
 ! var vector
 type, public :: var_ilength
  type(ilength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_ilength

 ! ** character type
 ! dat vector
 type, public :: clength
  character(len=32),allocatable       :: dat(:)    ! dat(:)
 endtype clength
 ! var vector
 type, public :: var_clength
  type(clength),allocatable           :: var(:)    ! var(:)%dat
 endtype var_clength

 ! ---------- forcing input file strcuture -----------------------------------------------------------------
  ! input file name and strcuture for nc files
  type, public ::  infileinfo
   integer(i4b)                            :: nTime            ! number of time step in a nc file
   integer(i4b)                            :: iTimebound(1:2)  ! time index of start and end of the
   real(dp)                 , allocatable  :: timevar(:)       ! the time varibale from the netcdf file
   real(dp)                                :: convTime2Days    ! the time varibale from the netcdf file
   real(dp)                                :: ncrefjulday      ! the julian day for the reference of the nc file
   character(len=strLen)                   :: infilename       ! the name of the input file name
   character(len=strLen)                   :: calendar         ! the calendar
   character(len=strLen)                   :: unit             ! the unit of time
  end type infileinfo


 ! ---------- mapping data structures ----------------------------------------------------------------------

 ! data to remap runoff hru to river network hrus
 type, public :: remap
   ! information in the mapping file
   integer(i4b)             , allocatable  :: hru_id(:)    ! Id of hrus associated with river network (="hru")
   integer(i4b)             , allocatable  :: qhru_id(:)   ! Id of hrus associated with runoff simulation (="qhru")
   integer(i4b)             , allocatable  :: num_qhru(:)  ! number of "qhru" within "hru"
   integer(i4b)             , allocatable  :: i_index(:)   ! Index in the y dimension of the runoff grid (starting with 1,1 in LL corner)
   integer(i4b)             , allocatable  :: j_index(:)   ! Index in the x dimension of the runoff grid (starting with 1,1 in LL corner)
   real(dp)                 , allocatable  :: weight(:)    ! area weight of "qhru" intersecting "hru"
   ! ancillary index vectors
   integer(i4b)             , allocatable  :: hru_ix(:)    ! Index of hrus associated with river network (="hru")
   integer(i4b)             , allocatable  :: qhru_ix(:)   ! Index of hrus associated with runoff simulation (="qhru")
 end type remap

 ! simulated runoff data
 type, public :: runoff
   integer(i4b)                            :: nTime         ! number of time steps
   integer(i4b)                            :: nSpace(1:2)   ! number of spatial dimension
   real(dp)                                :: time          ! time variable at one time step
   real(dp)                 , allocatable  :: sim(:)        ! flux simulation (HM_HRU) at one time step (size: nSpace(1))
   real(dp)                 , allocatable  :: sim2D(:,:)    ! flux simulation (x,y) at one time step (size: /nSpace(1),nSpace(2)/)
   integer(i4b)             , allocatable  :: hru_id(:)     ! id of HM_HRUs or RN_HRUs at which runoff is stored (size: nSpace(1))
   integer(i4b)             , allocatable  :: hru_ix(:)     ! Index of RN_HRUs associated with river network (used only if HM_HRUs = RN_HRUs)
   real(dp)                 , allocatable  :: basinRunoff(:)! remapped river network catchment runoff (size: number of nHRU)
   real(dp)                 , allocatable  :: basinEvapo(:) ! remapped river network catchment runoff (size: number of nHRU)
   real(dp)                 , allocatable  :: basinPrecip(:)! remapped river network catchment runoff (size: number of nHRU)
 end type runoff

 ! water management data; fluxes to/from reaches or target volume
 type, public :: wm
   integer(i4b)                            :: nTime           ! number of time steps
   integer(i4b)                            :: nSpace(1:2)     ! number of spatial dimension, in this case only one dimentonal
   real(dp)                                :: time            ! time variable at one time step
   real(dp)                 , allocatable  :: sim(:)          ! user specified flux add/subtract, or volume at one time step (size: nSpace)
   real(dp)                 , allocatable  :: sim2D(:,:)      ! to provide modularity for reading data
   integer(i4b)             , allocatable  :: seg_id(:)       ! id of reach in data (size: nSpace)
   integer(i4b)             , allocatable  :: seg_ix(:)       ! Index of river network reach IDs corresponding reach ID in data
   real(dp)                 , allocatable  :: flux_wm(:)      ! allocated flux to existing river network using sort_flux (size: number of nRCH)
   real(dp)                 , allocatable  :: vol_wm(:)       ! allocated target vol to existing river network using sort_flux (size: number of nRCH)
 end type

 ! ---------- reach parameters ----------------------------------------------------------------------------

 ! Reach Parameters
 type, public ::  RCHPRP
  real(dp)                                   :: R_SLOPE
  real(dp)                                   :: R_MAN_N
  real(dp)                                   :: R_WIDTH
  real(dp)                                   :: RLENGTH
  real(dp)                                   :: UPSAREA        ! upstream area (zero if headwater basin)
  real(dp)                                   :: BASAREA        ! local basin area
  real(dp)                                   :: TOTAREA        ! UPSAREA + BASAREA
  real(dp)                                   :: MINFLOW        ! minimum environmental flow
  real(dp)                                   :: D03_MaxStorage ! Doll 2003; maximume storage [m3]
  real(dp)                                   :: D03_Coefficient! Doll 2003; Coefficient [?]
  real(dp)                                   :: D03_Power      ! Doll 2003; Power [-]
  real(dp)                                   :: H06_Smax       ! Hanasaki 2006; maximume reservoir storage [m3]
  real(dp)                                   :: H06_alpha      ! Hanasaki 2006; fraction of active storage compared to total storage [-]
  real(dp)                                   :: H06_envfact    ! Hanasaki 2006; fraction of inflow that can be used to meet demand [-]
  real(dp)                                   :: H06_S_ini      ! Hanasaki 2006; initial storage used for initial estimation of release coefficient [m3]
  real(dp)                                   :: H06_c1         ! Hanasaki 2006; coefficient 1 for target release for irrigation reseroir [-]
  real(dp)                                   :: H06_c2         ! Hanasaki 2006; coefficient 2 for target release for irrigation reseroir [-]
  real(dp)                                   :: H06_exponent   ! Hanasaki 2006; Exponenet of actual release for "within-a-year" reservoir [-]
  real(dp)                                   :: H06_denominator! Hanasaki 2006; Denominator of actual release for "within-a-year" reservoir [-]
  real(dp)                                   :: H06_c_compare  ! Hanasaki 2006; Criterion for distinguish of "within-a-year" or "multi-year" reservoir [-]
  real(dp)                                   :: H06_frac_Sdead ! Hanasaki 2006; Fraction of dead storage to maximume storage [-]
  real(dp)                                   :: H06_E_rel_ini  ! Hanasaki 2006; Initial release coefficient [-]
  real(dp)                                   :: H06_I_Jan      ! Hanasaki 2006; Average January   inflow [m3/s]
  real(dp)                                   :: H06_I_Feb      ! Hanasaki 2006; Average Februrary inflow [m3/s]
  real(dp)                                   :: H06_I_Mar      ! Hanasaki 2006; Average March     inflow [m3/s]
  real(dp)                                   :: H06_I_Apr      ! Hanasaki 2006; Average April     inflow [m3/s]
  real(dp)                                   :: H06_I_May      ! Hanasaki 2006; Average May       inflow [m3/s]
  real(dp)                                   :: H06_I_Jun      ! Hanasaki 2006; Average June      inflow [m3/s]
  real(dp)                                   :: H06_I_Jul      ! Hanasaki 2006; Average July      inflow [m3/s]
  real(dp)                                   :: H06_I_Aug      ! Hanasaki 2006; Average August    inflow [m3/s]
  real(dp)                                   :: H06_I_Sep      ! Hanasaki 2006; Average September inflow [m3/s]
  real(dp)                                   :: H06_I_Oct      ! Hanasaki 2006; Average October   inflow [m3/s]
  real(dp)                                   :: H06_I_Nov      ! Hanasaki 2006; Average November  inflow [m3/s]
  real(dp)                                   :: H06_I_Dec      ! Hanasaki 2006; Average December  inflow [m3/s]
  real(dp)                                   :: H06_D_Jan      ! Hanasaki 2006; Average January   demand [m3/s]
  real(dp)                                   :: H06_D_Feb      ! Hanasaki 2006; Average Februrary demand [m3/s]
  real(dp)                                   :: H06_D_Mar      ! Hanasaki 2006; Average March     demand [m3/s]
  real(dp)                                   :: H06_D_Apr      ! Hanasaki 2006; Average April     demand [m3/s]
  real(dp)                                   :: H06_D_May      ! Hanasaki 2006; Average May       demand [m3/s]
  real(dp)                                   :: H06_D_Jun      ! Hanasaki 2006; Average June      demand [m3/s]
  real(dp)                                   :: H06_D_Jul      ! Hanasaki 2006; Average July      demand [m3/s]
  real(dp)                                   :: H06_D_Aug      ! Hanasaki 2006; Average Agust     demand [m3/s]
  real(dp)                                   :: H06_D_Sep      ! Hanasaki 2006; Average September demand [m3/s]
  real(dp)                                   :: H06_D_Oct      ! Hanasaki 2006; Average October   demand [m3/s]
  real(dp)                                   :: H06_D_Nov      ! Hanasaki 2006; Average November  demand [m3/s]
  real(dp)                                   :: H06_D_Dec      ! Hanasaki 2006; Average December  demand [m3/s]
  integer(i4b)                               :: H06_purpose    ! Hanasaki 2006; reservoir purpose; (0= non-irrigation, 1=irrigation) [-]
  logical(lgt)                               :: H06_I_mem_F    ! Hanasaki 2006; Flag to transition to modelled inflow [-]
  logical(lgt)                               :: H06_D_mem_F    ! Hanasaki 2006; Flag to transition to modelled/provided demand [-]
  integer(i4b)                               :: H06_I_mem_L    ! Hanasaki 2006; Memory length in years for inflow [year]
  integer(i4b)                               :: H06_D_mem_L    ! Hanasaki 2006; Memory length in years for demand [year]
 end type RCHPRP

 ! River Network topology
 type, public :: RCHTOPO
  integer(i4b)                               :: REACHIX      ! Reach index (1,2,...,nrch)
  integer(i4b)                               :: REACHID      ! Reach ID (REC code)
  real(dp)                                   :: RCHLAT1      ! Start latitude
  real(dp)                                   :: RCHLAT2      ! End latitude
  real(dp)                                   :: RCHLON1      ! Start longitude
  real(dp)                                   :: RCHLON2      ! End longitude
  integer(i4b)                               :: DREACHI      ! Immediate Downstream reach index
  integer(i4b)                               :: DREACHK      ! Immediate Downstream reach ID
  integer(i4b),dimension(:),allocatable      :: UREACHI      ! Immediate Upstream reach indices
  integer(i4b),dimension(:),allocatable      :: UREACHK      ! Immediate Upstream reach IDs
  integer(i4b),dimension(:),allocatable      :: RCHLIST      ! all upstream reach indices
  integer(i4b),dimension(:),allocatable      :: HRUID        ! all contributing HRU IDs
  integer(i4b),dimension(:),allocatable      :: HRUIX        ! all contributing HRU indices
  real(dp),    dimension(:),allocatable      :: HRUWGT       ! areal weight for contributing HRUs
  logical(lgt),dimension(:),allocatable      :: goodBas      ! Flag to denote a good basin
  character(len=32),dimension(:),allocatable :: pfafCode     ! pfafstetter code
  integer(i4b)                               :: RHORDER      ! Processing sequence
  real(dp)    ,dimension(:),allocatable      :: UH           ! Unit hydrograph for upstream
  integer(i4b)                               :: LAKE_IX      ! Lake index (1,2,...,nlak)
  integer(i4b)                               :: LAKE_ID      ! Lake ID (REC code)
  real(dp)                                   :: BASULAK      ! Area of basin under lake
  real(dp)                                   :: RCHULAK      ! Length of reach under lake
  logical(lgt)                               :: LAKINLT      ! .TRUE. if reach is lake inlet, .FALSE. otherwise
  logical(lgt)                               :: USRTAKE      ! .TRUE. if user takes from reach, .FALSE. otherwise
  logical(lgt)                               :: ISLAKE       ! .TRUE. if the object is a lake
  logical(lgt)                               :: LAKETARGVOL  ! .TRUE. if the lake follow a given target volume
  integer(i4b)                               :: LAKEMODELTYPE! 1=Doll, 2=Hanasaki, else=non-parameteric
 end type RCHTOPO

 ! ---------- reach states --------------------------------------------------------------------

 !---------- Lagrangian kinematic wave states (collection of particles) ---------------------------------
 ! Individual flow particles
 TYPE, public :: FPOINT
  real(dp)                                   :: QF           ! Flow
  real(dp)                                   :: QM           ! Modified flow
  real(dp)                                   :: TI           ! initial time of point in reach
  real(dp)                                   :: TR           ! time point expected to exit reach
  logical(lgt)                               :: RF           ! routing flag (T if point has exited)
 END TYPE FPOINT

 ! Collection of flow points within a given reach
 TYPE, public :: LKWRCH
  type(FPOINT),allocatable             :: KWAVE(:)
 END TYPE LKWRCH

 ! ---------- Eulerian kinematic wave ---------------------------------
 type, public :: EKWRCH
   real(dp)    :: Q(1:4)          ! Discharge at upstream and downstream of reach at current and previous time step(m3/s)
   real(dp)    :: A(1:4)          ! Flow area at upstream and downstream of reach at current and previous time step(m3/s)
 end type EKWRCH

 ! ---------- irf states (future flow series ) ---------------------------------
 ! Future flow series
 type, public :: IRFRCH
  real(dp), allocatable                :: qfuture(:)    ! runoff volume in future time steps for IRF routing (m3/s)
 END TYPE IRFRCH

 type, public :: STRSTA
  type(IRFRCH)       :: IRF_ROUTE
  type(LKWRCH)       :: LKW_ROUTE
  type(EKWRCH)       :: EKW_ROUTE
 end type STRSTA


 ! ---------- reach fluxes --------------------------------------------------------------------

 ! fluxes and states in each reach
 TYPE, public :: strflx
  real(dp), allocatable                :: QFUTURE(:)         ! runoff volume in future time steps (m3/s)
  real(dp), allocatable                :: QFUTURE_IRF(:)     ! runoff volume in future time steps for IRF routing (m3/s)
  real(dp), allocatable                :: QPASTUP_IRF(:,:)   ! runoff volume in the past time steps for lake upstream (m3/s)
  real(dp), allocatable                :: DEMANDPAST_IRF(:,:)! demand volume for lake (m3/s)
  real(dp)                             :: BASIN_QI           ! instantaneous runoff volume from the local basin (m3/s)
  real(dp)                             :: BASIN_QR(0:1)      ! routed runoff volume from the local basin (m3/s)
  real(dp)                             :: REACH_Q            ! time-step average streamflow (m3/s)
  real(dp)                             :: REACH_Q_IRF        ! time-step average streamflow (m3/s) from IRF routing
  real(dp)                             :: UPSTREAM_QI        ! sum of upstream streamflow (m3/s)
  real(dp)                             :: REACH_VOL(0:1)     ! volume of water at previous and current time step [m3]
  real(dp)                             :: REACH_WM_FLUX      ! water management fluxes to and from each reach
  real(dp)                             :: REACH_WM_VOL       ! target volume from the second water management file (m3)
  real(dp)                             :: TAKE               ! average take
  logical(lgt)                         :: isRoute            ! .true. if the reach is routed
  real(dp)                             :: basinEvapo         ! remapped river network catchment Evaporation (size: number of nHRU)
  real(dp)                             :: basinPrecip        ! remapped river network catchment Precipitation (size: number of nHRU)
 END TYPE strflx

 ! ---------- lake data types -----------------------------------------------------------------

 ! Lake Parameters
 TYPE, public :: LAKPRP
  real(dp)                             :: AREAREF            ! lake area
  real(dp)                             :: LAKREFLEV          ! lake elevation
  real(dp)                             :: LAKAVGLEV          ! lake average level (for initialization)
  real(dp)                             :: HE2AR_C            ! water height-surface area parameter
  real(dp)                             :: HE2AR_D            ! water height-surface area parameter
  real(dp)                             :: HGHTLOW            ! minimum water height for discharge
  real(dp)                             :: HGHTECO            ! minimum height for ecological concerns
  real(dp)                             :: HGHTSPL            ! spillway height
  real(dp)                             :: DSCHECO            ! discharge at "ecological" height
  real(dp)                             :: DSCHSPL            ! discharge at spillway height
  real(dp)                             :: RATECVA            ! discharge rating curve parameter
  real(dp)                             :: RATECVB            ! discharge rating curve parameter
 END TYPE LAKPRP

 ! Lake topology
 TYPE, public :: LAKTOPO
  integer(i4b)                         :: LAKE_IX           ! Lake index (1,2,...,nlak)
  integer(i4b)                         :: LAKE_ID           ! Lake ID (REC code)
  real(dp)                             :: LAKLAT1           ! Centroid latitude
  real(dp)                             :: LAKLAT2           ! Outlet latitude
  real(dp)                             :: LAKLON1           ! Centroid longitude
  real(dp)                             :: LAKLON2           ! Outlet longitude
  integer(i4b)                         :: DREACHI           ! Downstream reach index
  integer(i4b)                         :: DREACHK           ! Downstream reach ID
  integer(i4b)                         :: DLAKE_I           ! Downstream lake index
  integer(i4b)                         :: DLAKE_K           ! Downstream lake ID
 END TYPE LAKTOPO

 ! Lake fluxes
 TYPE, public :: LKFLX
  real(dp)                             :: LAKE_Qav          ! lake discharge (average over time step) (m3 s-1)
  real(dp)                             :: LAKE_Q            ! lake discharge (instantaneous) (m3 s-1)
  real(dp)                             :: LAKE_P            ! lake precipitation (m3)
  real(dp)                             :: LAKE_E            ! lake evaporation (m3)
  real(dp)                             :: LAKE_I            ! inflow to lake (m3 s-1)
 END TYPE LKFLX

END MODULE dataTypes

MODULE objTypes

 USE nrtype,     only: i4b,dp,lgt
 USE nrtype,     only: strLen   ! string length
 USE public_var, only: realMissing
 USE public_var, only: integerMissing

 ! define derived type for model variables, including name, description, and units
 type, public :: var_info_new
   character(len=strLen)    :: varName  = 'empty'         ! variable name
   character(len=strLen)    :: varDesc  = 'empty'         ! variable description
   character(len=strLen)    :: varUnit  = 'empty'         ! variable units
   integer(i4b)             :: varType  = integerMissing  ! variable type
   integer(i4b),allocatable :: varDim(:)                  ! dimension ID associated with variable
   logical(lgt)             :: varFile  = .true.          ! .true. if the variable should be read from a file
 CONTAINS
   procedure, pass :: init
 end type var_info_new

 CONTAINS

  SUBROUTINE init(this, vName, vDesc, vUnit, vType, vDim, vFile)
    implicit none
    class(var_info_new)                 :: this
    character(*),            intent(in) :: vName    ! variable name
    character(*),            intent(in) :: vDesc    ! variable description
    character(*),            intent(in) :: vUnit    ! variable units
    integer(i4b),            intent(in) :: vType    ! variable type
    integer(i4b),            intent(in) :: vDim(:)  ! dimension ID
    logical(lgt),            intent(in) :: vFile    ! .true. if the variable should be read from a file
    integer(i4b)                        :: n        ! size of dimension

    n = size(vDim)
    allocate(this%varDim(n))
    this%varName      = vName
    this%varDesc      = vDesc
    this%varUnit      = vUnit
    this%varType      = vType
    this%varDim(1:n) = vDim(1:n)
    this%varFile      = vFile
  END SUBROUTINE init

END MODULE objTypes
