!> Declaration of variables that control the configuration of a run, mostly set by commands.
! Variables related to output files are in the module "declare_output".

module declare_configuration

   use declare_limits
   use mloclib_messages
   
   implicit none
   save

   ! Logging and verbosity
   logical :: verbose_screen
   logical :: verbose_log
   logical :: debug
   
   ! Flags
   integer :: ibias
   integer :: istcor
   integer :: ittabl
   integer :: iwght
   integer :: longitude_range

   ! Windowing residuals
   real :: wind1
   real :: wind2
   real :: windloclim

   ! Epicentral distance ranges for hypocentroid and cluster vectors
   integer :: nclim
   integer :: nhlim
   real, dimension(3,2) :: clim
   real, dimension(3,2) :: hlim
   
   ! Use P arrivals only for hypocentroid
   logical :: ponly
   
   ! Data weighting and bias correction (in the Jordan & Sverdrup sense)
   logical :: bias_corr ! Bias correction of hypocentroid
   logical :: data_weight ! Weighting data by inverse of error
   logical :: pttt ! perfect theoretical travel times assumed for hypocentroid

    ! Test mode
    logical :: test_mode
   
   ! Selection of a direct calibration subcluster
   integer :: subc_nmin
   integer :: subc_nconnect
   real :: subc_delt
   logical :: subc_set
   
   ! Convergence limits
   integer :: convergence_test_index
   real :: cl_epi_h
   real :: cl_epi_c
   real :: cl_dep_h
   real :: cl_dep_c
   real :: cl_ot_h
   real :: cl_ot_c
   
   ! Kill events
   integer :: kill_count
   logical :: kill_all
   logical :: kill_one
   character(len=76) :: kill_reason
   
   ! Skip specified stations, phases and authors
   integer :: n_skip
   logical :: skip
   character(len=8), dimension(n_skip_max,3) :: skip_params
      
   ! Take starting locations from HDF file of a previous run
   integer, allocatable, dimension(:) :: hdf_day
   integer, allocatable, dimension(:) :: hdf_mon
   integer, allocatable, dimension(:) :: hdf_yr4
   integer :: nrhdf
   real, allocatable, dimension(:) :: hdf_dep
   real, allocatable, dimension(:) :: hdf_dep_minus
   real, allocatable, dimension(:) :: hdf_dep_plus
   real, allocatable, dimension(:) :: hdf_lat
   real, allocatable, dimension(:) :: hdf_lon
   real, allocatable, dimension(:) :: hdf_time
   character(len=2), allocatable, dimension(:) :: hdf_dep_code
   character(len=10), allocatable, dimension(:) :: hdf_evid
   character(len=60) :: rhdf_filnam
   logical :: read_hdf
      
contains

!*****************************************************************************************
   subroutine allocate_declare_configuration ()

      integer :: error
      integer :: n
      character(len=32) :: p
      
      n = nhdfmax
      p = 'allocate_declare_configuration'

      allocate (hdf_day(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_day', error)
      allocate (hdf_dep_code(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_dep_code', error)
      allocate (hdf_dep_minus(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_dep_minus', error)
      allocate (hdf_dep_plus(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_dep_plus', error)
      allocate (hdf_dep(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_dep', error)
      allocate (hdf_evid(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_evid', error)
      allocate (hdf_lat(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_lat', error)
      allocate (hdf_lon(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_lon', error)
      allocate (hdf_mon(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_mon', error)
      allocate (hdf_time(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_time', error)
      allocate (hdf_yr4(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdf_yr4', error)
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)

   end subroutine allocate_declare_configuration
   
   
!*****************************************************************************************
end module declare_configuration

