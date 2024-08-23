!> Declaration and allocation of variables related to seismic stations

module declare_stations

   use declare_limits
   use mloclib_messages
   
   implicit none
   save

   integer, allocatable, dimension(:,:) :: kcode
   real, allocatable, dimension(:,:) :: ahgts
   real, allocatable, dimension(:,:) :: azes
   real, allocatable, dimension(:,:) :: azse
   real, allocatable, dimension(:,:) :: delt
   real, allocatable, dimension(:,:,:) :: dt
   real, allocatable, dimension(:,:) :: elcr
   real, allocatable, dimension(:,:) :: psd
   real, allocatable, dimension(:,:,:) :: s
   real, allocatable, dimension(:,:) :: sca0s
   real, allocatable, dimension(:,:) :: sca1s
   real, allocatable, dimension(:,:) :: sca2s
   real, allocatable, dimension(:,:) :: scb1s
   real, allocatable, dimension(:,:) :: scb2s
   real, allocatable, dimension(:,:) :: sd1s
   real, allocatable, dimension(:,:) :: sd2s
   real, allocatable, dimension(:,:) :: sdstcs
   real, allocatable, dimension(:,:) :: stacs
   real, allocatable, dimension(:,:) :: stladg
   real, allocatable, dimension(:,:) :: stlndg
   real, allocatable, dimension(:,:) :: ttcomp
   real, allocatable, dimension(:,:) :: tto
   
   integer, allocatable, dimension(:) :: jdate_off
   integer, allocatable, dimension(:) :: jdate_on
   integer, allocatable, dimension(:) :: kstat
   integer :: n_failed_date_range
   integer :: n_radf
   integer :: n_supp_stn
   integer :: n_supp_stn_file
   integer :: nkstat
   integer :: nstat1
   integer :: nstn_used
   real, allocatable, dimension(:) :: ahgtr
   real, allocatable, dimension(:) :: sca0
   real, allocatable, dimension(:) :: sd1
   real, allocatable, dimension(:) :: stalat
   real, allocatable, dimension(:) :: stalon
   logical :: suppstn
   logical :: read_ad
   character(len=8), allocatable, dimension(:) :: sta_author
   character(len=8), allocatable, dimension(:) :: sta_deployment
   character(len=5), allocatable, dimension(:) :: nstr1
   character(len=5), allocatable, dimension(:) :: sta_agency
   character(len=5), allocatable, dimension(:) :: radf_stn
   character(len=13), allocatable, dimension(:) :: stn_dcal_used 
   character(len=100), allocatable, dimension(:) :: suppfilnam
   character(len=2), allocatable, dimension(:) :: sta_loc ! station location code from ADSLC
   character(len=3), allocatable, dimension(:) :: sta_cha ! station channel code from ADSLC

   ! List of stations suspected of reporting bogus depth phases
   integer :: n_bdp
   character(len=100) :: bdp_filnam
   character(len=108), dimension(n_bdp_max) :: bdp_station
   logical :: bdp_list
   
   ! Duplicate entries in station lists
   integer :: n_dupe
   integer :: n_dupe_conflict
   integer :: n_dupe_minor
   integer :: n_dupe_pure
   integer :: n_dupe_significant
   character(len=40), dimension(nmax1) :: duplication

   ! Missing station codes and other station statistics
   integer, dimension(n_miss_sta_max) :: n_miss_sta
   integer :: n_miss_sta_total
   integer, allocatable, dimension(:,:) :: ndat
   integer, allocatable, dimension(:,:) :: ndatc
   integer, allocatable, dimension(:,:) :: ndatdl
   integer, allocatable, dimension(:) :: ndatfl
   integer, allocatable, dimension(:,:) :: ndatpr
   integer, allocatable, dimension(:) :: nmiss ! number of missing stations for each event
   character(len=20), allocatable, dimension(:,:) :: missta
   character(len=20), dimension(n_miss_sta_max) :: n_miss_sta_list

   ! NEIC station metadata
    logical :: nsmd

   ! ISF-FDSN station metadata
   logical :: ifsm

contains

!*****************************************************************************************
   subroutine stations_allocate ()
   
   ! Allocate variable arrays related to stations
   
      integer :: error
      integer :: it
      integer :: it1
      integer :: l
      integer :: m
      integer :: n
      character(len=32) :: p
      
      l = n_miss_sta_max
      n = nevmax ! number of events in the cluster
      m = npmax ! maximum number of phase records for any event in the cluster
      it = itmax ! maximum number of iterations
      it1 = itmax1
      p = 'stations_allocate' ! procedure

!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (ahgtr(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'ahgtr', error)
      allocate (ahgts(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ahgts', error)
      allocate (azes(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'azes', error)
      allocate (azse(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'azse', error)
      allocate (delt(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'delt', error)
      allocate (dt(n,m,0:it), stat=error); if (error .gt. 0) call allocation_error (p, 'dt', error)
      allocate (elcr(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'elcr', error)
      allocate (jdate_off(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'jdate_off', error)
      allocate (jdate_on(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'jdate_on', error)
      allocate (kcode(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'kcode', error)
      allocate (kstat(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'kstat', error)
      allocate (missta(n,l), stat=error); if (error .gt. 0) call allocation_error (p, 'missta', error)
      allocate (ndat(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'ndat', error)
      allocate (ndatc(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'ndatc', error)
      allocate (ndatdl(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'ndatdl', error)
      allocate (ndatfl(n), stat=error); if (error .gt. 0) call allocation_error (p, 'ndatfl', error)
      allocate (ndatpr(n,0:it1), stat=error); if (error .gt. 0) call allocation_error (p, 'ndatpr', error)
      allocate (nmiss(n), stat=error); if (error .gt. 0) call allocation_error (p, 'nmiss', error)
      allocate (nstr1(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'nstr1', error)
      allocate (psd(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'psd', error)
      allocate (radf_stn(n_radf_max), stat=error); if (error .gt. 0) call allocation_error (p, 'radf_stn', error)
      allocate (s(n,m,0:it), stat=error); if (error .gt. 0) call allocation_error (p, 's', error)
      allocate (sca0(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sca0', error)
      allocate (sca0s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sca0s', error)
      allocate (sca1s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sca1s', error)
      allocate (sca2s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sca2s', error)
      allocate (scb1s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'scb1s', error)
      allocate (scb2s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'scb2s', error)
      allocate (sd1(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sd1', error)
      allocate (sd1s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sd1s', error)
      allocate (sd2s(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sd2s', error)
      allocate (sdstcs(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sdstcs', error)
      allocate (sta_agency(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sta_agency', error)
      allocate (sta_author(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sta_author', error)
      allocate (sta_cha(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sta_cha', error)
      allocate (sta_deployment(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sta_deployment', error)
      allocate (sta_loc(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'sta_loc', error)
      allocate (stacs(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'stacs', error)
      allocate (stalat(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'stalat', error)
      allocate (stalon(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'stalon', error)
      allocate (stladg(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'stladg', error)
      allocate (stlndg(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'stlndg', error)
      allocate (stn_dcal_used(nmax1), stat=error); if (error .gt. 0) call allocation_error (p, 'stn_dcal_used', error)
      allocate (suppfilnam(n_supp_stn_file_max), stat=error); if (error .gt. 0) call allocation_error (p, 'suppfilnam', error)
      allocate (ttcomp(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ttcomp', error)
      allocate (tto(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'tto', error)
      
      return
   
   end subroutine stations_allocate


!*****************************************************************************************
end module declare_stations

