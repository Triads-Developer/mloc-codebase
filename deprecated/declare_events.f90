!> Declaration and allocation of variables related to events
module declare_events

   use declare_limits ! to reference nevmax
   use mloclib_messages
   
   implicit none
   save

   integer, allocatable, dimension(:) :: mtiev
   integer, allocatable, dimension(:,:) :: mindx
   character(len=76), allocatable, dimension(:) :: event_comment
   character(len=10), allocatable, dimension(:) :: evid
   character(len=6), allocatable, dimension(:) :: evid_source
   character(len=30), allocatable, dimension(:) :: evtnam
   character(len=100), allocatable, dimension(:) :: infile ! event data file pathname (nev)
   character(len=20), allocatable, dimension(:) :: infile20 ! event data filename (nev)

   ! Depth-related
   integer :: n_cfd
   real, allocatable, dimension(:,:) :: depth_cf
   real :: depth_default
   real :: depth_default_minus
   real :: depth_default_plus
   real, allocatable, dimension(:,:) :: depth_hdf
   real, allocatable, dimension(:,:) :: depth_inp
   real, allocatable, dimension(:) :: depthp_minus
   real, allocatable, dimension(:) :: depthp_plus
   real :: median_constrained_depths
   real :: rescfd
   character(len=1), allocatable, dimension(:) :: depth_inp_c
   character(len=1), allocatable, dimension(:) :: depth_hdf_c
   character(len=1), allocatable, dimension(:) :: depth_cf_c
   character(len=1) :: depth_default_c
   
   ! Ellipse-related
   real, allocatable, dimension(:) :: alic
   real, allocatable, dimension(:) :: alphac
   real, allocatable, dimension(:) :: alphacg
   real :: alphah
   real, allocatable, dimension(:) :: blic
   real, allocatable, dimension(:,:,:) :: ccv
   real, allocatable, dimension(:,:) :: hcv
   real, allocatable, dimension(:) :: kcritc
   real, allocatable, dimension(:) :: xl1c
   real, allocatable, dimension(:) :: xl1cg
   real :: xl1h
   real, allocatable, dimension(:) :: xl2c
   real, allocatable, dimension(:) :: xl2cg
   real :: xl2h

   integer, allocatable, dimension(:) :: hour_inp
   integer, allocatable, dimension(:) :: idye
   integer, allocatable, dimension(:) :: iyre
   integer, allocatable, dimension(:) :: min_inp
   integer, allocatable, dimension(:) :: mone
   integer, allocatable, dimension(:) :: nevt
   real, allocatable, dimension(:) :: lat_inp
   real, allocatable, dimension(:) :: lon_inp
   real, allocatable, dimension(:) :: rmag
   real, allocatable, dimension(:) :: sec_inp
   real, allocatable, dimension(:) :: time_inp
   character(len=9), allocatable, dimension(:) :: hypo_author
   character(len=9), allocatable, dimension(:) :: magnitude_author
   character(len=2), allocatable, dimension(:) :: mmag
   
contains

!*****************************************************************************************
   subroutine events_allocate ()
   
   ! Allocate variable arrays related to events
   
      integer :: error
      integer :: n
      character(len=32) :: p
      
      n = nevmax
      p = 'events_allocate' ! procedure
      
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (depth_cf_c(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depth_cf_c', error)
      allocate (depth_cf(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'depth_cf', error)
      allocate (depth_hdf_c(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depth_hdf_c', error)
      allocate (depth_hdf(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'depth_hdf', error)
      allocate (depth_inp_c(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depth_inp_c', error)
      allocate (depth_inp(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'depth_inp', error)
      allocate (depthp_minus(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depthp_minus', error)
      allocate (depthp_plus(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depthp_plus', error)
      allocate (event_comment(n), stat=error); if (error .gt. 0) call allocation_error (p, 'event_comment', error)
      allocate (evid_source(n), stat=error); if (error .gt. 0) call allocation_error (p, 'evid_source', error)
      allocate (evid(n), stat=error); if (error .gt. 0) call allocation_error (p, 'evid', error)
      allocate (evtnam(n), stat=error); if (error .gt. 0) call allocation_error (p, 'evtnam', error)
      allocate (hour_inp(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hour_inp', error)
      allocate (hypo_author(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hypo_author', error)
      allocate (idye(n), stat=error); if (error .gt. 0) call allocation_error (p, 'idye', error)
      allocate (infile(n), stat=error); if (error .gt. 0) call allocation_error (p, 'infile', error)
      allocate (infile20(n), stat=error); if (error .gt. 0) call allocation_error (p, 'infile20', error)
      allocate (iyre(n), stat=error); if (error .gt. 0) call allocation_error (p, 'iyre', error)
      allocate (lat_inp(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lat_inp', error)
      allocate (lon_inp(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lon_inp', error)
      allocate (magnitude_author(n), stat=error); if (error .gt. 0) call allocation_error (p, 'magnitude_author', error)
      allocate (min_inp(n), stat=error); if (error .gt. 0) call allocation_error (p, 'min_inp', error)
      allocate (mindx(n,4), stat=error); if (error .gt. 0) call allocation_error (p, 'mindx', error)
      allocate (mmag(n), stat=error); if (error .gt. 0) call allocation_error (p, 'mmag', error)
      allocate (mone(n), stat=error); if (error .gt. 0) call allocation_error (p, 'mone', error)
      allocate (mtiev(n), stat=error); if (error .gt. 0) call allocation_error (p, 'mtiev', error)
      allocate (nevt(n), stat=error); if (error .gt. 0) call allocation_error (p, 'nevt', error)
      allocate (rmag(n), stat=error); if (error .gt. 0) call allocation_error (p, 'rmag', error)
      allocate (sec_inp(n), stat=error); if (error .gt. 0) call allocation_error (p, 'sec_inp', error)
      allocate (time_inp(n), stat=error); if (error .gt. 0) call allocation_error (p, 'time_inp', error)
      
      return
   
   end subroutine events_allocate


!*****************************************************************************************
end module declare_events
