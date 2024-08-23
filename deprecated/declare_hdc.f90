!> Declaration of global variables related to the hypocentroidal decomposition algorithm

module declare_hdc

   use declare_limits
   use mloclib_messages
   
   implicit none
   save
   

   integer, allocatable, dimension(:) :: jb
   integer, allocatable, dimension(:,:) :: ntq
   integer, allocatable, dimension(:) :: ntqi
   real, allocatable, dimension(:,:,:) :: a
   real, dimension(0:itmax1) :: otsh
   real, allocatable, dimension(:,:) :: otsp
   real, allocatable, dimension(:) :: sighatj
   real, allocatable, dimension(:) :: vnhat
   real, allocatable, dimension(:) :: wq
   real, allocatable, dimension(:) :: wq2
   
contains

!*****************************************************************************************
   subroutine hdc_allocate ()
   
   ! Allocate variable arrays related to HDC
   
      integer :: error
      integer :: it1
      integer :: j
      integer :: k
      integer :: l
      integer :: m
      integer :: n
      character(len=32) :: p
      
      it1 = itmax1
      j = ntmax1
      k = nqmax
      l = ntqmax
      m = npmax
      n = nevmax
      p = 'hdc_allocate' ! procedure
      
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (a(n,m,4), stat=error); if (error .gt. 0) call allocation_error (p, 'a', error)
      allocate (jb(j), stat=error); if (error .gt. 0) call allocation_error (p, 'jb', error)
      allocate (ntq(l,k), stat=error); if (error .gt. 0) call allocation_error (p, 'ntq', error)
      allocate (ntqi(k), stat=error); if (error .gt. 0) call allocation_error (p, 'ntqi', error)
      allocate (otsp(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'otsp', error)
      allocate (sighatj(j), stat=error); if (error .gt. 0) call allocation_error (p, 'sighatj', error)
      allocate (vnhat(j), stat=error); if (error .gt. 0) call allocation_error (p, 'vnhat', error)
      allocate (wq(k), stat=error); if (error .gt. 0) call allocation_error (p, 'wq', error)
      allocate (wq2(k), stat=error); if (error .gt. 0) call allocation_error (p, 'wq2', error)

      return
      
   end subroutine hdc_allocate
   
   
!*****************************************************************************************
end module declare_hdc

