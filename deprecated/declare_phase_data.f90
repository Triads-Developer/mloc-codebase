!> Declaration and allocation of global variables related to the phase arrival time data
module declare_phase_data

   use declare_limits
   use mloclib_messages
   
   implicit none
   save
   
   integer, allocatable, dimension(:,:) :: ipad
   integer, allocatable, dimension(:,:) :: ipah
   integer, allocatable, dimension(:,:) :: ipam
   integer, allocatable, dimension(:,:) :: ipamo
   integer, allocatable, dimension(:,:) :: ipay
   integer, allocatable, dimension(:,:) :: iptim
   integer, allocatable, dimension(:,:) :: mnf_line
   integer, allocatable, dimension(:) :: nst
   
   real, allocatable, dimension(:,:) :: pas
   real, allocatable, dimension(:,:) :: resisc
   real, allocatable, dimension(:,:) :: sdread
   real, allocatable, dimension(:,:) :: ttoff
   real, allocatable, dimension(:,:) :: ttsprd
   real, allocatable, dimension(:,:) :: weight
   
   character(len=52), allocatable, dimension(:,:) :: adslcaa
   character(len=5), allocatable, dimension(:,:) :: agency
   character(len=3), allocatable, dimension(:,:) :: channel
   character(len=8), allocatable, dimension(:,:) :: deployment
   character(len=10), allocatable, dimension(:,:) :: dist_az
   character(len=1), allocatable, dimension(:,:) :: fcode
   character(len=2), allocatable, dimension(:,:) :: location
   character(len=8), allocatable, dimension(:,:) :: phase
   character(len=8), allocatable, dimension(:,:) :: phase0
   character(len=8), allocatable, dimension(:,:) :: readsrc
   character(len=20), allocatable, dimension(:,:) :: sad
   character(len=5), allocatable, dimension(:,:) :: stname
   
   logical, allocatable, dimension(:,:) :: connected
   logical :: dflag ! if true, use data flags (x, p, d, k, l)
   logical, allocatable, dimension(:,:) :: fltrcdelt
   logical, allocatable, dimension(:,:) :: fltrcres
   logical, allocatable, dimension(:,:) :: fltrcflag
   logical, allocatable, dimension(:,:) :: fltrc
   logical, allocatable, dimension(:,:) :: fltrhdelt
   logical, allocatable, dimension(:,:) :: fltrhres
   logical, allocatable, dimension(:,:) :: fltrhflag
   logical, allocatable, dimension(:,:) :: fltrh
   logical :: no_xflags ! Ignore 'x' flags
   logical, allocatable, dimension(:,:) :: rel_depth_phase ! Relative depth phase flag
   logical, allocatable, dimension(:,:) :: rel_phase ! Relative phase flag

   ! Timing error corrections
   integer :: n_terr
   real, dimension(n_terr_max) :: terr_corr
   character(len=5), dimension(n_terr_max) :: terr_stn
   logical :: timing_error_correction
   
   ! Phase-specific default reading error values
   integer :: n_psdre
   real, dimension(n_psdre_max) :: psdre
   character(len=8), dimension(n_psdre_max) :: psdre_phase
   
   ! Relative depth phases
   integer, allocatable, dimension(:,:) :: h_rdp
   integer :: h_rdp1
   integer :: h_rdp2
   integer :: n_rdpp
   integer, allocatable, dimension(:) :: n_samples 
   real, allocatable, dimension(:,:) :: hrmsi
   real, allocatable, dimension(:,:) :: z_hrdp
   real, allocatable, dimension(:) :: z_test
   character(len=30), allocatable, dimension(:) :: rdpp_evt
   logical :: rdpp
   logical :: rdpp_all
   
   ! Differential time data
   integer, allocatable, dimension(:,:) :: diff_line
   integer, allocatable, dimension(:,:) :: idiff
   integer, allocatable, dimension(:) :: idiff0
   integer :: ndiff
   character(len=80) :: diffdatfilnam
   logical :: diffdat
   
   ! Reading errors
   integer, allocatable, dimension(:) :: indexq
   integer, allocatable, dimension(:,:) :: iqiev
   integer :: nqname
   real, allocatable, dimension(:) :: epa
   real, allocatable, dimension(:) :: ere
   real, allocatable, dimension(:) :: qelev
   real, allocatable, dimension(:) :: qlat
   real, allocatable, dimension(:) :: qlon
   real, allocatable, dimension(:,:) :: qres
   real :: rderr_loc_delt
   real :: rderr_loc_p
   real :: rderr_loc_s
   real :: rderr_min
   real :: rderr_min_depth
   real :: rderr_min_loc
   real, allocatable, dimension(:) :: rderr0
   real, allocatable, dimension(:) :: rdsigma
   character(len=21), allocatable, dimension(:) :: qname
   character(len=21), allocatable, dimension(:) :: qname1
   character(len=100) :: rderrfname
   logical :: read_rderr
   logical :: rels_set
   
contains

!*****************************************************************************************
   subroutine phase_data_allocate ()
   
   ! Allocate variable arrays related to phase data
   
      integer :: j
      integer :: k
      integer :: l
      integer :: m
      integer :: n
      integer :: error
      character(len=32) :: p
      
      j = n_qres_max
      k = nqmax
      l = n_hrdp_max
      m = npmax ! maximum number of phase records for any event in the cluster
      n = nevmax  ! number of events in the cluster
      p = 'phase_data_allocate'
      
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (adslcaa(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'adslcaa', error)
      allocate (agency(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'agency', error)
      allocate (channel(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'channel', error)
      allocate (connected(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'connected', error)
      allocate (deployment(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'deployment', error)
      allocate (diff_line(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'diff_line', error)
      allocate (dist_az(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'dist_az', error)
      allocate (epa(k), stat=error); if (error .gt. 0) call allocation_error (p, 'epa', error)
      allocate (ere(k), stat=error); if (error .gt. 0) call allocation_error (p, 'ere', error)
      allocate (fcode(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fcode', error)
      allocate (fltrc(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrc', error)
      allocate (fltrcdelt(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrcdelt', error)
      allocate (fltrcflag(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrcflag', error)
      allocate (fltrcres(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrcres', error)
      allocate (fltrh(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrh', error)
      allocate (fltrhdelt(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrhdelt', error)
      allocate (fltrhflag(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrhflag', error)
      allocate (fltrhres(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'fltrhres', error)
      allocate (h_rdp(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'h_rdp', error)
      allocate (hrmsi(n,l), stat=error); if (error .gt. 0) call allocation_error (p, 'hrmsi', error)
      allocate (idiff(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'idiff', error)
      allocate (idiff0(k), stat=error); if (error .gt. 0) call allocation_error (p, 'idiff0', error)
      allocate (indexq(k), stat=error); if (error .gt. 0) call allocation_error (p, 'indexq', error)
      allocate (ipad(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ipad', error)
      allocate (ipah(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ipah', error)
      allocate (ipam(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ipam', error)
      allocate (ipamo(n,0:m), stat=error); if (error .gt. 0) call allocation_error (p, 'ipamo', error)
      allocate (ipay(n,0:m), stat=error); if (error .gt. 0) call allocation_error (p, 'ipay', error)
      allocate (iptim(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'iptim', error)
      allocate (iqiev(k,j), stat=error); if (error .gt. 0) call allocation_error (p, 'iqiev', error)
      allocate (location(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'location', error)
      allocate (mnf_line(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'mnf_line', error)
      allocate (n_samples(n), stat=error); if (error .gt. 0) call allocation_error (p, 'n_samples', error)
      allocate (nst(n), stat=error); if (error .gt. 0) call allocation_error (p, 'nst', error)
      allocate (pas(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'pas', error)
      allocate (phase(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'phase', error)
      allocate (phase0(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'phase0', error)
      allocate (qelev(k), stat=error); if (error .gt. 0) call allocation_error (p, 'qelev', error)
      allocate (qlat(k), stat=error); if (error .gt. 0) call allocation_error (p, 'qlat', error)
      allocate (qlon(k), stat=error); if (error .gt. 0) call allocation_error (p, 'qlon', error)
      allocate (qname(k), stat=error); if (error .gt. 0) call allocation_error (p, 'qname', error)
      allocate (qname1(k), stat=error); if (error .gt. 0) call allocation_error (p, 'qname1', error)
      allocate (qres(k,j), stat=error); if (error .gt. 0) call allocation_error (p, 'qres', error)
      allocate (rderr0(k), stat=error); if (error .gt. 0) call allocation_error (p, 'rderr0', error)
      allocate (rdpp_evt(n), stat=error); if (error .gt. 0) call allocation_error (p, 'rdpp_evt', error)
      allocate (rdsigma(k), stat=error); if (error .gt. 0) call allocation_error (p, 'rdsigma', error)
      allocate (readsrc(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'readsrc', error)
      allocate (rel_depth_phase(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'rel_depth_phase', error)
      allocate (rel_phase(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'rel_phase', error)
      allocate (resisc(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'resisc', error)
      allocate (sad(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sad', error)
      allocate (sdread(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sdread', error)
      allocate (stname(n,0:m), stat=error); if (error .gt. 0) call allocation_error (p, 'stname', error)
      allocate (ttoff(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ttoff', error)
      allocate (ttsprd(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'ttsprd', error)
      allocate (weight(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'weight', error)
      allocate (z_hrdp(n,l), stat=error); if (error .gt. 0) call allocation_error (p, 'z_hrdp', error)
      allocate (z_test(n), stat=error); if (error .gt. 0) call allocation_error (p, 'z_test', error)
     
      return
   
   end subroutine phase_data_allocate


!*****************************************************************************************
end module declare_phase_data

