!> Declaration and allocation of variables related to the hypocentroid

module declare_hypocentroid

   use declare_limits
   use mloclib_messages
   
   implicit none
   save

   integer, dimension(0:itmax1) :: hourh
   integer, dimension(0:itmax1) :: minh
   integer, dimension(0:itmax1) :: nqmth
   integer :: nstep
   
   real, dimension(4) :: bcorr
   real, dimension(4,0:itmax1) :: delx0
   real, dimension(0:itmax1) :: depthh
   real, dimension(12) :: dimph
   real, allocatable, dimension(:,:) :: dtmph
   real, dimension(0:itmax1) :: ehatsqh
   real, allocatable, dimension(:,:) :: ehiev
   real :: hlatshift
   real :: hlonshift
   real :: hdepthshift
   real :: htimeshift
   real, dimension(0:itmax1) :: lath
   real, dimension(0:itmax1) :: lathgc ! hypocentroid latitude in geocentric coordinates.
   real, dimension(0:itmax1) :: lonh
   real, dimension(4) :: sdxhath
   real, dimension(0:itmax1) :: sech
   real, dimension(0:itmax1) :: shath
   
   character(len=5), dimension(4) :: fixprh
   
   logical :: depthfh
   logical :: latfh
   logical :: lonfh
   logical :: timefh
   
contains

!*****************************************************************************************
   subroutine hypocentroid_allocate ()
   
   ! Allocate variable arrays related to the hypocentroid
   
      integer :: error
      integer :: it1
      integer :: n
      integer :: m
      character(len=32) :: p
      
      n = nevmax
      m = npmax
      it1 = itmax1
      p = 'hypocentroid_allocate'
      
      allocate (dtmph(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'dtmph', error)
      allocate (ehiev(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'ehiev', error)
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)

      return
   
   end subroutine hypocentroid_allocate


!*****************************************************************************************
end module declare_hypocentroid

