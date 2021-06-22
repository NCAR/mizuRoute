module routing_param

!numeric type
USE nrtype
! global parameters
USE public_var,         only : realMissing    ! missing value for real number
USE public_var,         only : integerMissing ! missing value for integer number

! privary
implicit none
private

public::basinUH
public::make_uh

contains

  SUBROUTINE basinUH(dt, fshape, tscale, IERR, MESSAGE)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007 (FUSE)
  ! Martyn Clark, 2012 (modified to use allocatable arrays, and to remove option of no routing)
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the fraction of runoff in future time steps
  ! ---------------------------------------------------------------------------------------
  ! Modules Modified:
  ! -----------------
  ! MODULE reach_flux -- runoff fractions stored in FRAC_FUTURE(:)
  ! ---------------------------------------------------------------------------------------
  USE gamma_func_module, ONLY : gammp                   ! interface for the incomplete gamma function
  USE globalData,        ONLY : FRAC_FUTURE             ! fraction of runoff in future time steps
  IMPLICIT NONE
  ! input
  REAL(DP), INTENT(IN)                   :: dt          ! model time step
  REAL(DP), INTENT(IN)                   :: fshape      ! shape parameter in gamma distribution
  REAL(DP), INTENT(IN)                   :: tscale      ! time scale parameter in gamma distribution
  ! output
  INTEGER(I4B), INTENT(OUT)              :: IERR        ! error code
  CHARACTER(*), INTENT(OUT)              :: MESSAGE     ! error message
  ! locals
  REAL(DP)                               :: ntdh_min    ! minimum number of time delay points
  REAL(DP)                               :: ntdh_max    ! maximum number of time delay points
  REAL(DP)                               :: ntdh_try    ! trial number of time delay points
  INTEGER(I4B)                           :: itry        ! index of trial value
  INTEGER(I4B), PARAMETER                :: MAXTRY=100  ! maximum number of trial values
  INTEGER(I4B)                           :: NTDH        ! number of values on the time delay histogram
  INTEGER(I4B)                           :: JTIM        ! (loop through future time steps)
  REAL(DP)                               :: TFUTURE     ! future time (units of dt)
  REAL(DP)                               :: X_VALUE     ! xvalue to evaluate using gammp
  REAL(DP)                               :: CUMPROB     ! cumulative probability at JTIM
  REAL(DP)                               :: PSAVE       ! cumulative probability at JTIM-1
  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  ierr=0; message='basinUH/'
  ! use a Gamma distribution with shape parameter, fshape, and time parameter, tscale, input
  ! find the desired number of future time steps
  ! check if the cummulative Gamma distribution is close to 1.00 for given model time step, tscale and fsahpe.
  X_VALUE = dt/tscale 
  cumprob = gammp(fshape, X_VALUE)
  if(cumprob > 0.999_dp) then ! in case if the cumprob is close to 1 in one model time step
   ntdh_try = 1.999_dp
  else
   ntdh_min = 1._dp
   ntdh_max = 1000._dp
   ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
   do itry=1,maxtry
    x_value = dt*ntdh_try/tscale
    cumprob = gammp(fshape, x_value)
    !print*, tscale, ntdh_try, cumprob, x_value, itry
    if(cumprob < 0.99_dp)  ntdh_min = ntdh_try
    if(cumprob > 0.999_dp) ntdh_max = ntdh_try
    if(cumprob > 0.99_dp .and. cumprob < 0.999_dp) exit
    ntdh_try = 0.5_dp*(ntdh_min + ntdh_max)
    if(itry==maxtry)then; ierr=20; message=trim(message)//'cannot identify the maximum number of bins for the tdh'; return; endif
   end do
  endif
  ntdh = ceiling(ntdh_try)
  ! allocate space for the time-delay histogram
  if (.not.allocated(FRAC_FUTURE)) then
    allocate(FRAC_FUTURE(ntdh), stat=ierr)
    if(ierr/=0)then; message=trim(message)//'unable to allocate space for the time delay histogram'; return; endif
  endif
  ! loop through time steps and compute the fraction of runoff in future time steps
  PSAVE = 0.                                                 ! cumulative probability at JTIM-1
  DO JTIM=1,NTDH
   TFUTURE            = REAL(JTIM, kind(dp))*DT       ! future time
   CUMPROB            = gammp(fshape,TFUTURE/tscale)   ! cumulative probability at JTIM
   FRAC_FUTURE(JTIM)  = MAX(0._DP, CUMPROB-PSAVE)     ! probability between JTIM-1 and JTIM
   PSAVE              = CUMPROB                       ! cumulative probability at JTIM-1
   !WRITE(*,'(I5,1X,F20.5,1X,2(F11.5))') JTIM, TFUTURE, FRAC_FUTURE(JTIM), CUMPROB
  END DO
  ! ensure that the fractions sum to 1.0 (account for rounding errors, and not enough bins)
  FRAC_FUTURE(:) = FRAC_FUTURE(:) / SUM(FRAC_FUTURE(:))
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE basinUH


! *********************************************************************
! subroutine: compute normalized UH from Saint-Venant Eq. at sim. time step
!             for all the upstream segment
! *********************************************************************
 subroutine make_uh(&
                    ! Input
                    length,     &       ! input: river segment array [meter]
                    dt,         &       ! input: time step interval [sec]
                    velo,       &       ! input: IRF parameter 1 - celerity C for each stream segment [m/s]
                    diff,       &       ! input: IRF parameter 2 - diffusivity D for each stream segment [m^2/s]
                    ! Output
                    seg_uh,     &       ! output: unit hydrograph ordinates for a given segment length array
                    ierr, message)   ! output: error control
 ! ----------------------------------------------------------------------------------------
 ! Purpose:
 !
 !   Calculate UH at given simulation time step, C, and D for each stream segment
 !   based on Saint-Venant equation - see equations (13)-(15) in Lohmann, et al. (1996) Tellus
 !   IRF and UH are used interchangebly.
 !
 ! ----------------------------------------------------------------------------------------
  ! global variables
  USE public_var, only : pi      ! pi
  USE dataTypes,  only : dlength
  implicit none
  ! input variables
  real(dp),                      intent(in)   :: length(:)     ! river segment length
  real(dp),                      intent(in)   :: dt            ! Time step interval [sec]
  real(dp),                      intent(in)   :: velo          ! Wave velocity C for each segment [m/sec]
  real(dp),                      intent(in)   :: diff          ! Diffusivity D for each segment [m2/sec]
  ! output variables
  type(dlength),allocatable,     intent(out)  :: seg_uh(:)     ! unit hydrograph ordinates for each segment
  integer(i4b),                  intent(out)  :: ierr          ! error code
  character(*),                  intent(out)  :: message       ! error message
  ! local variables
  real(dp), allocatable                       :: UHM(:)        ! Unit hydrograph derived from S-V equation
  real(dp), allocatable                       :: UHQ(:)        !
  real(dp), allocatable                       :: fr(:)         ! Unit runoff depth evenly distributed over the simulation duration at hourly step
  real(dp)                                    :: seg_length    ! length of rech segment
  real(dp)                                    :: POT           ! Inside of exponential function of IRF
  real(dp)                                    :: H             ! IRF, h(x,t) function h function in eq.15 Lohmann et al.1996
  real(dp)                                    :: UHQ0          !
  real(dp)                                    :: INTE          ! Accumulation (integration) of UH
  real(dp)                                    :: sec           ! hr at each time step  [hr]
  real(dp),parameter                          :: dTUH=3600.0   ! UH time step [sec]
  integer(i4b)                                :: nSeg          ! Numer of segment
  integer(i4b)                                :: iHrStrt       ! index of UH time step where rising limb of UH start
  integer(i4b)                                :: iHrLast       ! index of UH time step where recession limb of UH become zero
  integer(i4b)                                :: nTSub         ! number of time steps where 1/nTsub [m] of runoff is inserted
  integer(i4b)                                :: iSeg          ! Loop index
  integer(i4b)                                :: iHr,jHr       ! Loop index of hour
  integer(i4b)                                :: iTagg         ! index for aggregated (i.e. simulation) time step
  integer(i4b),parameter                      :: nTMAX=240     ! Maximum hour of UH [hr] - 10 days times 24hrs
  integer(i4b),parameter                      :: nHr=240       ! Maximum hour of UH [hr] - 10 days times 24hrs
  character(len=strLen)                       :: cmessage      ! error message from subroutine

 ! initialize error control
 ierr=0; message='make_uh/'

 ! Dynamically assigned parameters
 nTsub=ceiling(dt/dTUH)
 !nTsub=floor(dt/dTUH)
 nSeg = size(length)

 ! Memory allocation
 allocate(seg_uh(nSeg), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': seg_uh'; return; endif
 allocate(fr(nTMAX), stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': fr'; return; endif
 allocate(UHQ(nTMAX),stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': UHQ'; return; endif
 allocate(UHM(nHr),stat=ierr, errmsg=cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage)//': UHM'; return; endif

 ! Generate Unit Runoff depth [1/nTsub, 1/nTsub, ..., 1/nTsub, 0, 0, ...0]
 fr=0._dp
 fr(1:nTsub) = 1.0_dp / nTsub

 ! starting segment loop
 do iSeg=1,nSeg

  !Compute S-V based UH, UHM
  seg_length = length(iSeg)
  !seg_length = RPARAM(iSeg)%RLENGTH

  INTE = 0._dp
  sec = 0._dp
  UHM(:) = 0._dp
  do iHr=1,nHr
    sec = sec + dTUH
    if (velo .GT. 0.0) then
      POT = ((velo*sec-seg_length)**2.0)/(4.0_dp*diff*sec)
      if (POT .GT. 69.0_dp) then
        H = 0.0_dp
      else
        H = 1.0_dp/(2.0_dp*sqrt(pi*diff*sec))*seg_length*exp(-POT)
      endif
    else
      H = 0.0_dp
    endif
    UHM(iHr) = H
    INTE = INTE + H
  enddo

  !Normalize ordinates of S-V UH by its sum
  if (INTE > 0.0) then
    UHM= UHM/INTE
  end if

  !Get hour indices of start of rising part of UHM and end of recession part of UHM
  INTE = 0._dp
  do iHr=1,nTMAX
    INTE= INTE+UHM(iHr)
    iHrLast=iHr
    if (INTE > 0.99999) exit
  enddo
  INTE = 0._dp
  do iHr=nTMAX,1,-1
    INTE= INTE+UHM(iHr)
    iHrStrt=iHr
    if (INTE > 0.99999) exit
  enddo

  !Convolute with unit runoff depth, fr
  INTE = 0._dp
  UHQ(:) = 0._dp
  do jHr = 1,nTMAX
    UHQ0=0._dp
    do iHr = iHrStrt,iHrLast
      if ((jHr-iHr) > 0) then
        if (jHr-iHr <= nTsub ) then
          UHQ0 = UHQ0 + fr(jHr-iHr)*UHM(iHr)
        endif
      else
        exit
      endif
    enddo
    UHQ(jHr) = UHQ0
    INTE = INTE+UHQ0
  end do

  ! Normalize ordinates of UHQ by its sum
  if (INTE > 0.0_dp) then
    UHQ = UHQ/INTE
  end if

  !Get hour indices of end of recession part of UH
  INTE = 0._dp
  do iHr=1,nTMAX
    INTE= INTE+UHQ(iHr)
    iHrLast=iHr
    if (INTE > 0.9999) exit
  enddo

  ! Re-normalize the UHQ by its sum
  UHQ = UHQ/INTE

  !Aggregate hourly unit hydrograph to simulation time step
  allocate(seg_uh(iSeg)%dat((iHrLast+nTsub-1)/nTsub),stat=ierr,errmsg=cmessage)
  if(ierr/=0)then; message=trim(message)//trim(cmessage)//': seg_uh%dat'; return; endif

  seg_uh(iSeg)%dat(:)=0._dp
  do jHr = 1,iHrLast
    iTagg = (jHr+nTsub-1)/nTsub
    seg_uh(iSeg)%dat(iTagg) = seg_uh(iSeg)%dat(iTagg) + UHQ(jHr)
  end do

 end do ! sSeg loop

 end subroutine make_uh


end module routing_param

