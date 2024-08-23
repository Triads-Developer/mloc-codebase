!> Declaration and allocation of variables related to the cluster vectors

module declare_cluster_vectors

   use declare_limits
   use mloclib_messages
   
   implicit none
   save

   integer, allocatable, dimension(:,:) :: hourp
   integer, allocatable, dimension(:,:) :: minp
   integer :: ntc
   integer :: nqc
   
   real, allocatable, dimension(:,:) :: depthp
   real, dimension(12) :: dimpc
   real, allocatable, dimension(:,:) :: dimpciev
   real, allocatable, dimension(:,:) :: dtmpc
   real, allocatable, dimension(:,:,:) :: dxp
   real, allocatable, dimension(:,:) :: eci
   real, allocatable, dimension(:,:) :: eciev
   real, dimension(0:itmax1) :: ehatsqc
   real, allocatable, dimension(:) :: lat_cf
   real, allocatable, dimension(:) :: lat_hdf
   real, allocatable, dimension(:) :: lon_cf
   real, allocatable, dimension(:) :: lon_hdf
   real :: radius_cvff
   real, allocatable, dimension(:,:) :: sdxhatc
   real, dimension(0:itmax1) :: shatc
   real, allocatable, dimension(:) :: shatsqci
   real :: tikhonov_factor
   real, allocatable, dimension(:) :: time_cf
   real, allocatable, dimension(:) :: time_hdf
   real, allocatable, dimension(:,:) :: latp
   real, allocatable, dimension(:,:) :: latpgc ! event latitude in geocentric coordinates.
   real, allocatable, dimension(:,:) :: lonp
   real, allocatable, dimension(:,:) :: secp
   
   double precision :: shatsqc
   double precision, allocatable, dimension(:,:) :: vhatc

   character(len=1), allocatable, dimension(:) :: depset_pr
   character(len=5), allocatable, dimension(:,:) :: fixpr
   
   logical :: damping
   logical, allocatable, dimension(:) :: depthf
   logical, allocatable, dimension(:) :: latf
   logical, allocatable, dimension(:) :: lonf
   logical, allocatable, dimension(:) :: timef
   
contains
   
!*****************************************************************************************
   subroutine cluster_vectors_allocate ()
   
   ! Allocate variable arrays related to cluster vectors
   
      integer :: error
      integer :: it1
      integer :: it2
      integer :: m
      integer :: n
      character(len=32) :: p
      
      n = nevmax
      m = npmax
      it1 = itmax1
      it2 = itmax2
      p = 'cluster_vectors_allocate' ! procedure
      
      allocate (depset_pr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depset_pr', error)
      allocate (depthf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depthf', error)
      allocate (depthp(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'depthp', error)
      allocate (dimpciev(n,12), stat=error); if (error .gt. 0) call allocation_error (p, 'dimpciev', error)
      allocate (dtmpc(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'dtmpc', error)
      allocate (dxp(n,4,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'dxp', error)
      allocate (eci(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'eci', error)
      allocate (eciev(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'eciev', error)
      allocate (fixpr(n,4), stat=error); if (error .gt. 0) call allocation_error (p, 'fixpr', error)
      allocate (hourp(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'hourp', error)
      allocate (lat_cf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lat_cf', error)
      allocate (lat_hdf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lat_hdf', error)
      allocate (latf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'latf', error)
      allocate (latp(n,0:it2), stat=error); if (error .gt. 0) call allocation_error (p, 'latp', error)
      allocate (latpgc(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'latpgc', error)
      allocate (lon_cf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lon_cf', error)
      allocate (lon_hdf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lon_hdf', error)
      allocate (lonf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lonf', error)
      allocate (lonp(n,0:it2), stat=error); if (error .gt. 0) call allocation_error (p, 'lonp', error)
      allocate (minp(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'minp', error)
      allocate (sdxhatc(n,4), stat=error); if (error .gt. 0) call allocation_error (p, 'sdxhatc', error)
      allocate (secp(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'secp', error)
      allocate (shatsqci(n), stat=error); if (error .gt. 0) call allocation_error (p, 'shatsqci', error)
      allocate (time_cf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'time_cf', error)
      allocate (time_hdf(n), stat=error); if (error .gt. 0) call allocation_error (p, 'time_hdf', error)
      allocate (timef(n), stat=error); if (error .gt. 0) call allocation_error (p, 'timef', error)
      allocate (vhatc(mtmax,mtmax), stat=error); if (error .gt. 0) call allocation_error (p, 'vhatc', error)
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)

      return
   
   end subroutine cluster_vectors_allocate


!*****************************************************************************************
end module declare_cluster_vectors

