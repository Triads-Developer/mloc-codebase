!> Declaration of global variables related to seismic phases
module declare_phases

   use declare_limits
   use mloclib_messages
   
   implicit none
   save

   ! Junk Phases
   integer :: n_junk
   character(len=8), dimension(n_junk_max) :: junk_phase
   
   ! Phase re-identification
   integer :: n_no_phreid
   character(len=1), allocatable, dimension(:,:) :: no_phid ! Flag to prevent phase reidentification for particular readings
   character(len=8), dimension(n_no_phreid_max) :: no_phreid
   logical :: phid ! .true. if phase identification is to be done
   logical :: plogout ! ~.plog output file
   logical, allocatable, dimension(:,:) :: phidird ! Individual reading flag to allow change of phase ID if .true.
   
   ! Spread and baseline offset of specific phases
   integer :: nsprd
   real, dimension(nsprdmax) :: sprd
   real, dimension(nsprdmax) :: ttoffset
   logical :: read_ttsprd
   character(len=8), dimension(nsprdmax) :: sprdph
   character(len=100) :: ttsprdfname
   
    ! Crustal heterogeneity factor for differential-time data
   real :: cscale
   real :: vmr
   
contains

!*****************************************************************************************
   subroutine phases_allocate ()
   
   ! Allocate variable arrays related to seismic phases
   
      integer :: error
      integer :: m
      integer :: n
      character(len=32) :: p
      
      m = npmax
      n = nevmax
      p = 'phases_allocate' ! procedure
      
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (no_phid(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'no_phid', error)
      allocate (phidird(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'phidird', error)
      
      return
   
   end subroutine phases_allocate
  
   
!*****************************************************************************************
end module declare_phases

