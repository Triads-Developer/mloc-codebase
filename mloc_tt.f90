module mloc_tt

   use mloc_declare
   use mloc_lib
   use mloc_taup
   use mloc_hyposat

   implicit none
   save
   
   
   contains

!*****************************************************************************************
   subroutine mlocsteq (iev, istat, rlatdg1, rlondg1, stladg1, stlndg1, hp, sth, phase1, t,&
    delt1, azes1, azse1, dtddx, dtdhx, elcr1, hgtcr, eqlat, bptc_print)
             
   !  Version of mlocsteq for use with tau-p software
   !  9/19/94 by EAB.
   
   !  for a given station-earthquake pair, locsteq returns the travel time,
   !  epicentral distance, azimuth and back-azimuth, ray parameter, and
   !  ellipticity and station elevation corrections to travel time.
   
   !  input:
   !     rlatdg1  earthquake latitude (geocentric degrees)
   !     rlondg1  earthquake longitude (geocentric degrees)
   !     stladg1  station latitude (geocentric degrees)
   !     stlndg1  station longitude (geocentric degrees)
   !     hp      earthquake depth (km down +)
   !     sth     station elevation (km up +)
   !     phase1   phase name
   !     eqlat epicenter latitude in geographic degress
   !     bptc_print logical to print output for bounce point corrections
   
   !  output:
   !     t       travel time
   !     delt1    epicentral distance (degrees)
   !     azes1    azimuth (degrees clockwize from north)
   !     azse1    back azimuth (degrees clockwize from north)
   !     dtddx   travel time derivative w.r.t. delta: ray parameter (sec/rad)
   !     dtdhx   travel time derivative w.r.t. depth.
   !     elcr1    ellipticity correction (dziewonski and gilbert)
   !     hgtcr   station elevation correction
   
   !  subroutines called:
   !    delaz
   !    trtm
   !    tt_mixed_model
   !    ellip

      ! Travel time variables
      integer :: nphase
      real, dimension(max_mpl) :: dtdd
      real, dimension(max_mpl) :: dtdh
      real, dimension(max_mpl) :: dddp
      real, dimension(max_mpl) :: tt ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd

      integer :: i
      integer :: iev
      integer :: ii
      integer :: istat
      real, intent(out) :: azes1
      real, intent(out) :: azse1
      real :: bp_crust_corr
      real :: bp_water_corr
      real :: bplat
      real :: bplon
      real :: crvel
      real, intent(out) :: delt1
      real, intent(out) :: dtddx
      real, intent(out) :: dtdhx
      real :: dum1
      real :: dum2
      real :: dum3
      real :: dum4
      real, intent(out) :: elcr1
!      real :: ellip
      real, intent(in) :: eqlat
      real, intent(out) :: hgtcr
      real, intent(in) :: hp
      real, intent(in) :: rlatdg1
      real, intent(in) :: rlondg1
      real, intent(in) :: sth
      real, intent(in) :: stladg1
      real, intent(in) :: stlndg1
      real, intent(out) :: t
      real :: t1
      real :: t2
      real :: t3
      real :: t4
      real :: topo
      character(len=1) :: bdp_flag
      character(len=132) :: msg
      character(len=8), intent(in) :: phase1
      character(len=8) :: phase2
      logical, intent(in) :: bptc_print
      logical :: error
      
      topo = 0.
      bp_crust_corr = 0.
      bp_water_corr = 0.
      bdp_flag = ' '
   
      !  Epicentral distance, azimuth, and back-azimuth
      t1 = rlatdg1*rpd  ! Convert to geocentric radians
      t2 = rlondg1*rpd  ! Convert to geocentric radians
      t3 = stladg1*rpd  ! Convert to geocentric radians
      t4 = stlndg1*rpd  ! Convert to geocentric radians
      call delaz (t1, t2, t3, t4, dum1, delt1, dum2, dum3, azes1, dum4, azse1, 1)
      if (delt1 .lt. 0. .or. delt1 .gt. 180.) then
         write (msg,'(a,f10.3,a,i3,2a)') 'mlocsteq: illegal value for delta: ', delt1, 'event ',&
          iev, ' Station ', stname(iev,istat)//' '//phase(iev,istat)
         call oops (trim(msg))
      end if
           
      !  Theoretical travel-time
      if (.not.locmod) then
         call trtm (delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
      else
         call tt_mixed_model (hp, delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
      end if
      if (nphase .eq. 0) then
         write (msg,'(a,f10.3)') 'mlocsteq: no phases returned for delt1 = ', delt1
         call warnings (trim(msg))
         if (log_tt) write (io_tt_log,'(a)') trim(msg)
         nphase = 1
         phcd(1) = 'CRAP    '
      end if
      
      ! PKP precursors: calculate travel time and derivatives for PKPdf. This does not change the name.
      ! pwP: calculate travel time and derivatives for pP. This does not change the name. The correction
      ! for propagation in the water column is done later.
      phase2 = phase1
      If (phase1 .eq. 'PKPdfpre') phase2 = 'PKPdf   '      
      if (phase1 .eq. 'pP-P    ') phase2 = 'pP      '
      if (phase1 .eq. 'sP-P    ') phase2 = 'sP      '
      If (phase1 .eq. 'pwP     ' .or. phase1 .eq. 'pwP-P   ') phase2 = 'pP      '
      
      ! Find the phase in the list of returned phases
      ii = 0
      do i = 1,nphase
         if (phase2 .eq. phcd(i)) then
            ii = i
            dtddx = dtdd(ii)
            dtdhx = dtdh(ii)
            t = tt(ii)
            exit ! This is necessary to avoid picking up duplicate phases at later times
         end if
      end do
      
      if (ii .gt. 0) then
      
         ! Correction to depth phase TT for bounce point topography.
         ! No correction for zero focal depth, even though there would be a theoretical pP and sP phase from non-zero topography.
         if (phase2 .eq. 'pP      ' .or. phase2 .eq. 'sP      ') then ! This catches pwP also.
            if (bptc .and. (hp .gt. 0.)) then
               call dpbp (phase2, rlatdg1, rlondg1, azes1, dtddx, bplat, bplon, error)
               call find_topo (bplat, bplon, topo)
               call topo_corr (phase2, topo, dtddx, bp_crust_corr, bp_water_corr)
               t = t + bp_crust_corr
               if (phase1 .eq. 'pwP     ' .or. phase1 .eq. 'pwP-P   ') t = t + bp_water_corr
            else
               if (phase1 .eq. 'pwP     ' .or. phase1 .eq. 'pwP-P   ') then 
                  write (msg,'(a)') 'mlocsteq: there are pwP phases but bouncepoint/water leg correction has not been made'
                  call warnings (trim(msg))
               end if      
            end if
         end if
         
         ! Reduction to relative depth phases
         if (rel_depth_phase(iev,istat)) then
            t = t - tt(1)
            if (bptc .and. bptc_print) then
               if (bdp(iev,istat)) then ! Check for stations known to report bad depth phases
                  bdp_flag = '*'
               else
                  bdp_flag = ' '
               end if
               write (io_bptc_log,'(1x,2a,1x,a,2f10.4,3f6.2,f7.2)',advance='no') bdp_flag,&
                stname(iev,istat), phase1, bplat, bplon, topo, bp_crust_corr, bp_water_corr, t
            end if
         end if
                  
      else ! The phase was not found in the standard phase list
      
         if (phase2(1:7) .ne. 'UNKNOWN' .and.&
             phase2(1:1) .ne. 'X' .and.&
             phase2(1:1) .ne. ' ' .and.&
             phase2(1:2) .ne. 'P_' .and.&
             phase2(1:2) .ne. 'Lg' .and.&
             phase2(1:2) .ne. 'Rg' .and.&
             phase2(1:2) .ne. 'T ' .and.&
             phase2(1:3) .ne. 'S-P' .and.&
             .not.skipp(phase2)) then
            if (log_tt) write (io_tt_log,'(a,i3,1x,a,1x,2a,f8.1,a)') 'mlocsteq: ', iev,&
             stname(iev,istat), phase1, ' not found at ', delt1, ' degrees'
         end if
         dtddx = 0.
         dtdhx = 0.
         t = 0.
         
         ! P_ travel time
         if (phase2(1:2) .eq. 'P_') then
            dtddx = p__b
            dtdhx = 0.
            t = p__a + delt1*p__b
         end if
         
         ! Lg travel time
         if (phase2(1:2) .eq. 'Lg') then
            dtddx = lg_b
            dtdhx = 0.
            t = lg_a + delt1*lg_b
         end if
         
         ! Rg travel time
         if (phase2(1:2) .eq. 'Rg') then
            dtddx = rg_b
            dtdhx = 0.
            t = rg_a + delt1*rg_b
         end if
         
         ! T-phase travel time
         if (phase2(1:2) .eq. 'T ') then
            dtddx = tphase_b
            dtdhx = 0.
            t = tphase_a + delt1*tphase_b
         end if
         
         ! S-P relative phase
         if (phase2(1:3) .eq. 'S-P') then
            ! Find the first S phase
            do i = 1,nphase
               if (stype(phcd(i))) exit
            end do
            t = tt(i) - tt(1)
            dtddx = dtdd(i) - dtdd(1)
            dtdhx = dtdh(i) - dtdh(1)
            if (debug) write (io_log,'(a8,f10.2,4x,a8,f10.2,3f10.1)') phcd(1), tt(1), phcd(i),&
             tt(i), dtdd(1), dtdd(i), dtddx
         end if
         
      end if
      
      if (debug)  write (io_log,*) phase2, t, dtddx, dtdhx
      
      ! convert to sec/radian
      dtddx = dtddx*57.29578
      
      !  ellipticity correction (sec).
      if (phase1 .eq. 'S-P') then
         elcr1 = 0.
      else if (phase1 .eq. 'pP-P    ') then
         elcr1 = 0.
      else if (phase1 .eq. 'sP-P    ') then
         elcr1 = 0.
      else if (phase1 .eq. 'pwP-P   ') then
         elcr1 = 0.
      else
         elcr1 = ellip (phase2, delt1, hp, eqlat, azes1, .false.)
      end if
      
      ! Station elevation correction (sec). Use S velocity crust if 'S' is last leg of phase name.
      if (phase1 .eq. 'S-P') then
         hgtcr = 0.
      else if (phase1 .eq. 'pP-P    ') then
         hgtcr = 0.
      else if (phase1 .eq. 'sP-P    ') then
         hgtcr = 0.
      else if (phase1 .eq. 'pwP-P   ') then
         hgtcr = 0.
      else if (phase1 .eq. 'P_      ') then
         hgtcr = 0.
      else if (phase1 .eq. 'Lg      ') then
         hgtcr = 0.
      else if (phase1 .eq. 'Rg      ') then
         hgtcr = 0.
      else if (phase1 .eq. 'T       ') then
         hgtcr = 0.
      else if (phase1(1:3) .eq. 'UNK') then ! UNKNOWN phase
         hgtcr = 0.
      else
         if (ptype(phase2)) then
            crvel = pcrvel
         else if (stype(phase2)) then
            crvel = scrvel
         else
            if (log_tt) write (io_tt_log,'(3a)') 'mlocsteq: phase ', trim(phase2),&
             ' not classified for elevation correction'
            hgtcr = 0.
            return
         end if
         call get_elev_corr (crvel, dtddx, sth, hgtcr)
      end if
      
      return
      
   end subroutine mlocsteq
   
   
!*****************************************************************************************
   subroutine get_elev_corr (vel, dtdd, elev, corr)
   
   ! Bondar's algorithm for station elevation correction in the ISC locator, converted from C code
   ! Input:
   ! vel: surface velocity in km/sec
   ! dtdd: ray parameter in sec/rad
   ! elev: station elevation in km
   ! Return:
   ! corr: elevation correction in sec
   
      real, intent(out) :: corr
      real, intent(in) :: dtdd
      real, intent(in) :: elev
      real, intent(in) :: vel
   
      corr = vel*(dtdd/6371.)
      corr = corr**2
      if (corr .gt. 1.) corr = 1./corr
      corr = sqrt(1. - corr)
      corr = corr*(elev/vel)
      
      return
   
   end subroutine get_elev_corr
   
   
!*****************************************************************************************
   subroutine tt_mixed_model (focal_depth, delta, nphase_out, tt_out, dtdd_out, dtdh_out, dddp_out, phcd_out)
         
   ! For a specified epicentral distance and focal depth, this subroutine returns a list of phases for
   ! which theoretical travel times (and derivatives) are calculated, for the case when there is a local
   ! crustal model available. In this case, the local model is used for all crustal phases out to the distance limit
   ! specified in the model file itself. If crustal arrivals (e.g., Pn) exist beyond this distance in the data set,
   ! they will be calculated with the global model. The global model is used for all teleseismic phases.
   ! Some "teleseismic" phases may be observed within the distance range associated with the local crustal model.
         
      ! Travel time variables
      integer :: nphase
      real, dimension(max_mpl) :: dtdd
      real, dimension(max_mpl) :: dtdh
      real, dimension(max_mpl) :: dddp
      real, dimension(max_mpl) :: tt ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd

      integer :: nphase_loc
      real, dimension(max_mpl) :: dtdd_loc
      real, dimension(max_mpl) :: dtdh_loc
      real, dimension(max_mpl) :: dddp_loc
      real, dimension(max_mpl) :: tt_loc ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd_loc
      
      integer :: i
      integer :: ierr
      integer, intent(out) :: nphase_out
      real, dimension(max_mpl), intent(out) :: dddp_out
      real, intent(in) :: delta
      real, dimension(max_mpl), intent(out) :: dtdd_out
      real, dimension(max_mpl), intent(out) :: dtdh_out
      real, intent(in) :: focal_depth
      real, dimension(max_mpl), intent(out) :: tt_out
      character(len=132) :: msg
      character(len=8), dimension(max_mpl), intent(out) :: phcd_out
      
      nphase_out = 0
      nphase_loc = 0
      if (log_tt) write (io_log,'(a,4f10.2)') 'tt_mixed_model: ', delta, dlimlocmod, focal_depth, zlimlocmod
   
      ! Local model
      if (delta .le. dlimlocmod .and. focal_depth .le. zlimlocmod) then
         call ttloc2 (focal_depth, delta, 'D', nphase_loc, tt_loc, dtdd_loc, dtdh_loc, dddp_loc, phcd_loc, ierr, 30000)
         if (log_tt) write (io_tt_log,'(a,f5.1,a,i3)') 'tt_mixed_model: Local model phases, delta = ', delta,&
          '; nphase_loc = ', nphase_loc
         do i = 1,nphase_loc
            if (crustal_phase(phcd_loc(i))) then
               nphase_out = nphase_out + 1
               phcd_out(nphase_out) = phcd_loc(i)
               tt_out(nphase_out) = tt_loc(i)
               dtdd_out(nphase_out) = dtdd_loc(i)
               dtdh_out(nphase_out) = dtdh_loc(i)
               dddp_out(nphase_out) = dddp_loc(i)
               if (log_tt) write (io_tt_log,'(i3,2x,a,f10.3,3e10.3)') nphase_out, phcd_out(nphase_out),&
                tt_out(nphase_out), dtdd_out(nphase_out), dtdh_out(nphase_out), dddp_out(nphase_out)
            end if
         end do
      end if
      if (log_tt) write (io_tt_log,'(a,i3)') 'nphase_out from local model = ', nphase_out
      
      
      ! Global model
      call trtm (delta, nphase, tt, dtdd, dtdh, dddp, phcd)
      if (log_tt) write (io_tt_log,'(a,f5.1,a,i3)') 'tt_mixed_model: Global model phases, delta = ', delta,&
          '; nphase = ', nphase
      do i = 1,nphase
         if (.not.crustal_phase(phcd(i)) .or. (crustal_phase(phcd(i)) .and. delta .gt. dlimlocmod)) then
            nphase_out = nphase_out + 1
            if (nphase_out .gt. max_mpl) then
               write (msg,'(a,i3,a)') 'tt_mixed_model: maximum number of phases (', max_mpl, ') exceeded'
               call oops (trim(msg))
            end if
            phcd_out(nphase_out) = phcd(i)
            tt_out(nphase_out) = tt(i)
            dtdd_out(nphase_out) = dtdd(i)
            dtdh_out(nphase_out) = dtdh(i)
            dddp_out(nphase_out) = dddp(i)
            if (log_tt) write (io_tt_log,'(i3,2x,a,f10.3,3e10.3)') nphase_out, phcd_out(nphase_out),&
             tt_out(nphase_out), dtdd_out(nphase_out), dtdh_out(nphase_out), dddp_out(nphase_out)
         else
            if (log_tt) write (io_tt_log,'(i3,1x,2a)') i, phcd(i), ' was skipped'
         end if
      end do
      if (log_tt) write (io_tt_log,'(a,i3)') 'nphase_out total = ', nphase_out
      
      return
      
   end subroutine tt_mixed_model
   
   
!*****************************************************************************************
   real function ellip (phasein, edist, edepth_in, slat, bazim, init)
   
   ! Ellipticity correction for any given phase using Dziewonski & Gilbert
   ! representation. The ellipticity corrections are found by linear
   ! interpolation in terms of values calculated for the ak135 model for a
   ! wide range of phases to match the output of the iasp software 
   
   ! Input Parameters:
   
   ! character  
   !  phase  : A string specifying the PHASE, e.g., P, ScP etc.  
   !            'phase' should have at least 8 characters length.
                                                        
   ! real 
   !  edist  : Epicentral distance to station (in degrees)     
   !  edepth_in : Depth of event (km)        
   !  slat   : Epicentral latitude of source (in degrees) 
   !  bazim  : Azimuth from source to station (in degrees)
                                
   ! Output:
   
   !  real
   !   ellip : Time correction for path to allow for ellipticity, added to
   !            AK135 travel times. Returns ellip=0. for unknown phases. 
   
   ! Usage:
   
   !  One call: tcor=ellip(phase,edist,edepth_in,slat,bazim) to initialize.
   !  Every next call computes ellip from input parameters.
   
   ! numph: number of phases supported by the tau-tables.
   ! numph should exactly match the number of phases in the tau-tables.
   ! To add a phase to the software proceed as follows:
   !   . increase numph by 1
   !   . add the phase code at the end of the phcod data statement
   !   . add the proper tau-table entries at the bottom of file 'tau.tables'
   ! The only syntax check on the tau.table file is that 
   ! the order in which the phases appear in the file should be 
   ! the same as the order of the phases in the phcod data statement
   
   ! B.L.N. Kennett RSES, ANU, May 1995. Based on earlier routine by D.J. Brown.  
   ! Modified for processing of large data sets by W. Spakman, Earth Sciences, Utrecht University, June 1996.
   ! Cleaned up and slightly modified for use with "mloc" by Eric Bergman, June 1999.
   ! Further clean-up, update to f90 by eab, April 21, 2013.
   
   ! Modified to take variable "init" as an argument 2024/5/7 by EAB.
   
      integer, parameter :: nd = 6 ! Number of depths in table
      integer, parameter :: numph = 57
      integer, parameter :: mdel = 50
      real, dimension(nd), parameter :: dpth = (/0.,100.,200.,300.,500.,700./)
   
      ! Formerly in ellip.inc to preserve values between calls
      integer, save, dimension(numph) :: phnch
      integer, save, dimension(numph) :: np
      real, save, dimension(numph,mdel) :: delta
      real, save, dimension(numph) :: di1
      real, save, dimension(numph) :: di2
      real, save, dimension(numph,mdel,6) :: t0
      real, save, dimension(numph,mdel,6) :: t1
      real, save, dimension(numph,mdel,6) :: t2
      character(len=8), save, dimension(numph) :: phcod

      integer :: i
      integer :: i1
      integer :: i2
      integer :: idist
      integer :: ios
      integer :: ip
      integer :: j
      integer :: jdepth
      integer :: k
      integer :: l
!      integer :: lunit
      integer :: lut
      integer :: m
      integer :: nc
!       integer, dimension(numph) :: np
!       integer, dimension(numph) :: phnch
      real :: a0
      real :: a1
      real :: a2
      real :: azim
      real :: b0
      real :: b1
      real :: b2
      real, intent(in) :: bazim
      real :: caz
      real :: cbz
      real :: d0
      real :: d1
      real :: d2
      real, parameter :: degrad = 0.01745329
      real, parameter :: deldst = 5.0
!       real, dimension(numph,mdel) :: delta
!       real, dimension(numph) :: di1
!       real, dimension(numph) :: di2
!       real, dimension(nd), parameter :: dpth = (/0.,100.,200.,300.,500.,700./)
      real :: e0
      real :: e1
      real :: e2
      real :: ecolat
      real, intent(in) :: edepth_in
      real :: edepth
      real, intent(in) :: edist
      real :: f0
      real :: f1
      real :: f2
      real :: g0
      real :: g1
      real :: g2
      real :: h0
      real :: h1
      real :: h2
      real :: sc0
      real :: sc1
      real :: sc2
      real, intent(in) :: slat
!       real, dimension(numph,mdel,nd) :: t0
!       real, dimension(numph,mdel,nd) :: t1
!       real, dimension(numph,mdel,nd) :: t2
      real :: tau0
      real :: tau1
      real :: tau2
      real :: tcor
      character(len=80) :: fname
      character(len=132) :: msg
      character(len=*), intent(in) :: phasein
!       character(len=8), dimension(numph) :: phcod
      character(len=8) :: phdum
      logical, intent(in) :: init

      data phcod/&
      "Pup     ","P       ","Pdiff   ","PKPab   ","PKPbc   ","PKPdf   ",&
      "PKiKP   ","pP      ","pPKPab  ","pPKPbc  ","pPKPdf  ","pPKiKP  ",&
      "sP      ","sPKPab  ","sPKPbc  ","sPKPdf  ","sPKiKP  ","PcP     ",&
      "ScP     ","SKPab   ","SKPbc   ","SKPdf   ","SKiKP   ","PKKPab  ",&
      "PKKPbc  ","PKKPdf  ","SKKPab  ","SKKPbc  ","SKKPdf  ","PP      ",&
      "P'P'    ","Sup     ","S       ","Sdiff   ","SKSac   ","SKSdf   ",&
      "pS      ","pSKSac  ","pSKSdf  ","sS      ","sSKSac  ","sSKSdf  ",&
      "ScS     ","PcS     ","PKSab   ","PKSbc   ","PKSdf   ","PKKSab  ",&
      "PKKSbc  ","PKKSdf  ","SKKSac  ","SKKSdf  ","SS      ","S'S'    ",&
      "SP      ","PS      ","PnS     "/
   
      ellip = 0.
   
      if (init) then ! Initialize at first call and return.
   
         ! Check on the length of phase
         l = len(phasein)
         if (l .lt. 8) then
            msg = 'ellip: phase name must be at least 8 characters'
            call oops (trim(msg))
         end if
   
         ! Initialize arrays
         do i = 1,numph
            phnch(i) = len_trim(phcod(i))
            np(i) = 0
            di1(i) = 0.
            di2(i) = 0.
            do j = 1,mdel
               delta(i,j) = 0.
               do k = 1,6
                  t0(i,j,k) = 0.
                  t1(i,j,k) = 0.
                  t2(i,j,k) = 0.
               end do
            end do
         end do
   
         ! Open tau.table
         if (verbose_screen) call fyi ('ellip: open and read ellipticity tau-table')
         lut = lunit() ! Find an unused logical unit number to read tau.table.
         fname = trim(ellip_path)//dirsym//'tau.table'
         call open_file (lut, fname, 'old')
         
         ! Read tau.table
         ip = 0
         do
            ip = ip + 1
            if (ip .gt. numph) exit ! stop reading numph: limit is reached'
      
            ! phase_code, number_of_epicentral_distance_entries, min_max_dist.
            read (lut,'(a8,i2,2f10.1)',iostat=ios) phdum, np(ip), di1(ip), di2(ip)
            if (ios .lt. 0) then ! EOF
               exit
            else if (ios .gt. 0) then
              write (msg,'(a,i6,a)') 'ellip: read error on phase record err= ', ios, ' on tau-table'
              call oops (trim(msg))
            end if
            if (debug) write (io_log,'(a,i3,1x,a8,1x,i2,2f10.1)') 'ellip: reading ', ip, phdum, np(ip), di1(ip), di2(ip)
            if (phdum(1:8) .ne. phcod(ip)) then
               msg = 'ellip: syntax of tau-table does not conform to syntax of phase data statement'
               call oops (trim(msg))
            end if
            do i = 1,np(ip)
               read (lut,*,iostat=ios) delta(ip,i)
               if (ios.gt. 0) then
                  write (msg,'(a,i6,a)') 'ellip: read error ',ios,' on tau-table'
                  call oops (trim(msg))
               end if
               read (lut,*,iostat=ios) (t0(ip,i,m),m=1,6)
               if (ios.gt. 0) then
                  write (msg,'(a,i6,a)') 'ellip: read error ',ios,' on tau-table'
                  call oops (trim(msg))
               end if
               read (lut,*,iostat=ios) (t1(ip,i,m),m=1,6)
               if (ios.gt. 0) then
                  write (msg,'(a,i6,a)') 'ellip: read error ',ios,' on tau-table'
                  call oops (trim(msg))
               end if
               read (lut,*,iostat=ios) (t2(ip,i,m),m=1,6)
               if (ios.gt. 0) then
                  write (msg,'(a,i6,a)') 'ellip: read error ',ios,' on tau-table'
                  call oops (trim(msg))
               end if
            end do
         end do
   
         if (log_tt) then
            write (io_tt_log,'(a,i5)') 'ellip: Number of phases: ', ip-1
            write (io_tt_log,'(a)') 'ellip: Phase codes     : '
            i = int(float(numph)/6.)
            do j = 1,i
               i1 = 1 + (j-1)*6
               i2 = i1 - 1 + 6
               write (io_tt_log,'(7a9)') (phcod(k),k=i1,i2)
            end do
            i = i*6 + 1
            j = numph - i
            if (j .gt. 0) then
               write (io_tt_log,'(7a9)') (phcod(k),k=i,numph)
            end if
            write (io_tt_log,'(/a/)') 'See also the phase aliases in routine phase_alias'
         end if
       
         close (lut)
         return
           
      end if
   
      ! Set up source dependent constants
   
      azim = bazim*degrad
      ecolat = (90.0-slat)*degrad
      call ellref (ecolat, sc0, sc1, sc2)
   
      if (debug) then
         write (io_log,'(a)') 'ellip:'
         write (io_log,'(a)') 'phase, edist, edepth, ecolat, azim'
         write (io_log,'(a8,1x,f8.4,1x,f6.2,1x,f6.4,1x,f6.4)') phasein(1:8), edist, edepth_in, ecolat, azim
         write (io_log,'(a,3(1x,e10.4))') 'sc0, sc1, sc2: ', sc0, sc1, sc2
      end if
   
      ! Select phase
      ! In addition to the phase names listed above a number of phase aliases are available in the routine phase_alias,
      ! e.g. Pn --> P etc. The input phase code is first checked against the phcod array and next against the phase aliases.
      ip = -1
      nc = min(len_trim(phasein),8)
      do i = 1,numph
         if (nc .ne. phnch(i)) cycle
         if (phasein(1:nc) .eq. phcod(i)(1:nc)) then
            ip = i
            exit
         end if
      end do
      if (ip .eq. -1) call phase_alias (phasein, edist, ip)
   
      if (ip .gt. 0) then
         if (debug) write (io_log,'(a,i3,a,a8,a,2f8.1)') 'ellip: ip: ',ip,' Selected phase: ',&
          phcod(ip),' Table distance range: ', di1(ip), di2(ip)
      else             
         if (debug) write (io_log,'(a,a8,a)') 'ellip: Selected phase: ', phasein(1:8), ' is not available'
         ellip = 0.    
         return
      end if
   
      ! distance index
      idist = 1 + int( (edist-di1(ip))/ deldst )
      if (edist .lt. di1(ip)) idist = 1
      if (edist .gt. di2(ip)) idist = np(ip) - 1
   
      ! depth index
      ! Check for a negative depth; this can happen when a shift in depth is applied for
      ! indirect calibration that takes a very shallow event to a negative value.
      if (edepth_in .lt. 0.) then
         write (msg,'(a,f6.2,a)') 'ellip: negative depth = ', edepth_in, ' adjusted to 1.0 km'
         call warnings (trim(msg))
         edepth = 1.0
      else
         edepth = edepth_in
      end if
      do j = 1,nd-1
         if ((dpth(j) .le. edepth) .and. (dpth(j+1) .ge. edepth)) then
            jdepth = j
            exit
         end if
      end do
   
      if (debug) write (io_log,*) 'idist, jdepth;', idist, jdepth
   
      ! Compute tau-values and ellip.
      ! Need to allow for zero entries (where phase description strongly depth dependent)
   
      ! tau0
      a0 = t0(ip,idist,jdepth)
      b0 = t0(ip,idist,jdepth+1)
      h0 = t0(ip,idist+1,jdepth+1)
      d0 = t0(ip,idist+1,jdepth)
      e0 = a0 + (d0-a0)*(edist-delta(ip,idist))/(delta(ip,idist+1)-delta(ip,idist))
      f0 = b0 + (h0-b0)*(edist-delta(ip,idist))/(delta(ip,idist+1)-delta(ip,idist))
      g0 = e0 + (f0-e0)*(edepth-dpth(jdepth))/ (dpth(jdepth+1)-dpth(jdepth))
      tau0 = g0
   
      ! tau1
      a1 = t1(ip,idist,jdepth)
      b1 = t1(ip,idist,jdepth+1)
      h1 = t1(ip,idist+1,jdepth+1)
      d1 = t1(ip,idist+1,jdepth)
      e1 = a1 + (d1-a1)*(edist-delta(ip,idist))/ (delta(ip,idist+1)-delta(ip,idist))
      f1 = b1 + (h1-b1)*(edist-delta(ip,idist))/ (delta(ip,idist+1)-delta(ip,idist))
      g1 = e1 + (f1-e1)*(edepth-dpth(jdepth))/ (dpth(jdepth+1)-dpth(jdepth))
      tau1 = g1
   
      ! tau2
      a2 = t2(ip,idist,jdepth)
      b2 = t2(ip,idist,jdepth+1)
      h2 = t2(ip,idist+1,jdepth+1)
      d2 = t2(ip,idist+1,jdepth)
      e2 = a2 + (d2-a2)*(edist-delta(ip,idist))/ (delta(ip,idist+1)-delta(ip,idist))
      f2 = b2 + (h2-b2)*(edist-delta(ip,idist))/ (delta(ip,idist+1)-delta(ip,idist))
      g2 = e2 + (f2-e2)*(edepth-dpth(jdepth))/ (dpth(jdepth+1)-dpth(jdepth))
      tau2 = g2
   
      caz = cos(azim)
      cbz = cos(2.0*azim)
   
      if (debug) then
         write (io_log,*) 'tau0, tau1, tau2:', tau0, tau1, tau2
         write (io_log,*) 'azim, caz, cbz', azim, caz, cbz    
      end if
   
      tcor = sc0*tau0 + sc1*cos(azim)*tau1 + sc2*cos(2.0*azim)*tau2
      ellip = tcor
   
      return
         
   end function ellip
         
   
!*****************************************************************************************
   subroutine ellref (ecolat, sc0, sc1, sc2)
      
      !Arguments
      real, intent(in) :: ecolat
      real, intent(out) :: sc0
      real, intent(out) :: sc1
      real, intent(out) :: sc2
      
      real, parameter :: s3 = 0.8660254 ! s3 = sqrt(3.0)/2.0
      
      sc0 = 0.25*(1.0+3.0*cos(2.0*ecolat))
      sc1 = s3*sin(2.0*ecolat)
      sc2 = s3*sin(ecolat)*sin(ecolat)
      
      return
      
   end subroutine ellref
         
   
!*****************************************************************************************
   subroutine redtab ()
   
   ! Read tau-p tables.
         
      character(len=8), dimension(10) :: phlst ! Phase list for tau-p
      logical, dimension(3) :: prflg = (/.false.,.false.,.true./) ! Information printed from tau-p calculations
      
      ! The keywords in array phlst do the following:
      !      P      gives P-up, P, Pdiff, PKP, and PKiKP
      !      P+     gives P-up, P, Pdiff, PKP, PKiKP, PcP, pP, pPdiff, pPKP,
      !             pPKiKP, sP, sPdiff, sPKP, and sPKiKP
      !      S+     gives S-up, S, Sdiff, SKS, sS, sSdiff, sSKS, pS, pSdiff, and pSKS
      !      basic  gives P+ and S+ as well as ScP, SKP, PKKP, SKKP, PP, and P'P'
      ! Note that generic S gives S-up, Sdiff, and SKS already and so doesn't require a keyword.
      
      data phlst(1)/'all'/
      
      if (verbose_screen) call fyi ('redtab: calling tabin...')
      call tabin (io_taup, taup_model)
      if (verbose_screen) call fyi ('redtab: returned from tabin; calling brnset...')
      call brnset (1, phlst, prflg)
      if (verbose_screen) call fyi ('redtab: returned from brnset...')
      
      return
      
   end subroutine redtab
         
         
!*****************************************************************************************
   subroutine ttsig (phase, delta, spread, offset)
   
   ! For a given phase and epicentral distance, returns a measure of the spread of readings. This is not reading error,
   ! but rather a measure of the heterogeneity of the real Earth. This spread is only used in the estimation of the hypocentroid,
   ! not the cluster vectors. No dependence on delta, so far. Offset is not yet calculated.
   ! Basic values are mostly taken from Kennett, Engdahl and Buland (1995), Fig 3.
   ! added 5/20/02 by eab.
         
      real, intent(in) :: delta
      real, intent(out) :: spread
      real, intent(out) :: offset
      character(len=8), intent(in) :: phase
      
      offset = 0.
      
      if (phase .eq. 'P       ') then
         if (delta .ge. 28.) then
            spread = 1.2
         else
            spread = 1.5
         end if
      else if (phase .eq. 'pP      ') then
         if (delta .ge. 28.) then
            spread = 1.5
         else
            spread = 2.0
         end if
      else if (phase .eq. 'sP      ') then
         if (delta .ge. 28.) then
            spread = 2.0
         else
            spread = 2.5
         end if
      else if (phase .eq. 'Pg      ') then
         spread = amax1(0.15, 0.8*delta)
      else if (phase .eq. 'Pb      ') then
         spread = amax1(0.15, 0.8*delta)
      else if (phase .eq. 'Pn      ') then
         spread = 2.5
      else if (phase .eq. 'PP      ') then
         spread = 2.0      
      else if (phase .eq. 'PcP     ') then
         spread = 1.7      
      else if (phase .eq. 'S       ') then
         spread = 2.0      
      else if (phase .eq. 'Sg      ') then
         spread = amax1(0.3, 1.6*delta)
      else if (phase .eq. 'Sb      ') then
         spread = amax1(0.3, 1.6*delta)
      else if (phase .eq. 'Sn      ') then
         spread = 3.0      
      else if (phase .eq. 'SS      ') then
         spread = 2.2      
      else if (phase .eq. 'ScS     ') then
         spread = 2.0      
      else if (phase .eq. 'SP      ') then
         spread = 2.1      
      else if (phase .eq. 'ScP     ') then
         spread = 1.9      
      else if (phase .eq. 'PKiKP   ') then
         spread = 1.8      
      else if (phase .eq. 'PKPdf   ') then
         spread = 1.8      
      else if (phase .eq. 'PKPbc   ') then
         spread = 1.8      
      else if (phase .eq. 'PKPab   ') then
         spread = 1.8      
      else if (phase .eq. 'PKKPbc  ') then
         spread = 1.8      
      else if (phase .eq. 'PKKPab  ') then
         spread = 1.9      
      else if (phase .eq. 'PKKPdf  ') then
         spread = 1.9      
      else if (phase .eq. 'SKSac   ') then
         spread = 2.0      
      else if (phase .eq. 'SKKSac  ') then
         spread = 2.0      
      else if (phase .eq. 'SKPdf   ') then
         spread = 2.0      
      else if (phase .eq. 'SKPbc   ') then
         spread = 1.9
      else if (phase .eq. 'P_      ') then
         spread = 5.0
      else if (phase .eq. 'Lg      ') then
         spread = 5.0
      else if (phase .eq. 'Rg      ') then
         spread = 5.0
      else if (phase .eq. 'T       ') then
         spread = 5.0
      else if (phase .eq. 'S-P     ') then
         spread = 0.50
      else
         spread = 2.0
      end if
      
      return
      
   end subroutine ttsig
   
   
!*****************************************************************************************
   subroutine ttsig2 (phase2, delta, spread, offset)
   
   ! Spread of observations for different phases, read from a .ttsprd file output from a previous run.
   ! No dependence on delta so far, except for Lg and Rg. See ttsig.
         
      integer :: i
      real, intent(inout) :: delta
      real, intent(out) :: offset
      real, intent(out) :: spread
      character(len=8), intent(inout) :: phase2
      
      do i = 1,nsprd
         if (phase2 .eq. sprdph(i)) then
            if (phase2(1:2) .eq. 'P_') then
               spread = amax1(delta*2.0,1.5)
               offset = ttoffset(i)
            else if (phase2(1:2) .eq. 'Lg') then
               spread = amax1(delta*2.0,1.5)
               offset = ttoffset(i)
            else if (phase2(1:2) .eq. 'Rg') then
               spread = amax1(delta*2.0,1.5)
               offset = ttoffset(i)
            else
               spread = sprd(i)
               offset = ttoffset(i)
            end if
            return
         end if
      end do
      
      ! Use defaults if phase was not found in the TTSPRD file
      call ttsig (phase2, delta, spread, offset)
            
      return
      
   end subroutine ttsig2
         
         
!*****************************************************************************************
   subroutine rdttsprd ()
   
   ! Get info on travel time spread
         
      integer :: i
      integer :: ios
      real :: offsetmax = 4.00
      real :: sprdmax = 3.50
      real :: sprdmin = 0.70
      character(len=132) :: msg
      logical :: op
      
      if (read_ttsprd) then
         
         call open_file (io_ttsprd, ttsprdfname, 'old')
         
         inquire (unit=io_ttsprd,opened=op,name=ttsprdfname)
         if (op) then
            nsprd = 0
            do i = 1,n_sprd_max
               read (io_ttsprd,'((a8,8x,2f10.3))',iostat=ios) sprdph(i), sprd(i), ttoffset(i)
               if (ios .lt. 0) exit
               sprd(i) = amax1(sprd(i),sprdmin) ! Minimum allowed TT spread
               sprd(i) = amin1(sprd(i),sprdmax) ! Maximum allowed TT spread
               ttoffset(i) = amin1(ttoffset(i),offsetmax) ! Maximum allowed TT offset
               nsprd = nsprd + 1
            end do
            close (io_ttsprd)
            if (verbose_screen) then
               write (msg,'(a,i3,2a)') 'rdttsprd: ', nsprd, ' phase spread entries read from ', trim(ttsprdfname)
               call fyi (trim(msg))
            end if
         else
            msg = 'rdttsprd: failed to open TTSPRD file '//trim(ttsprdfname)//'...using default values'
            call warnings (trim(msg))
            ttsprdfname = 'default'
            read_ttsprd = .false.
         end if
      else
         ttsprdfname = 'default'
         read_ttsprd = .false.
      end if
      
      return
      
   end subroutine rdttsprd
   
   
!*****************************************************************************************
   subroutine dpbp (phase_in, elat, elon, az, drp2, bplat, bplon, error)
   
   ! Bounce point for pP or sP
   ! Event latitiude and longitude (elat, elon) in geocentric coordinates
   
   ! Adapted from Bob Engdahl's code 'bounce' by EAB.
   
      integer :: ierr
!      integer :: iupcor
      real, intent(in) :: az
      real :: bp2
      real :: bpaz
      real :: bpdel
      real, intent(out) :: bplat
      real, intent(out) :: bplon
      real :: bptim
      real, intent(in) :: drp2
      real, intent(in) :: elat
      real, intent(in) :: elon
      character(len=8), intent(in) :: phase_in
      character(len=132) :: msg
      logical :: error
      
      error = .false.
      bplat = elat
      bplon = elon
      bp2 = abs(drp2)
      bpaz = az ! For depth phases, the azimuth of the bouncepoint is the same as for P
   
      if (phase_in .eq. 'pP      ') then
         ierr = iupcor ('P', bp2, bpdel, bptim)
      else if (phase_in .eq. 'sP      ') then
         ierr = iupcor ('S', bp2, bpdel, bptim)
      else
         write (msg,'(a)') 'dpbp: unsupported phase ('//trim(phase_in)//')'
         call warnings (trim(msg))
         error = .true.
         return
      end if
      
      if (ierr .lt. 0) then
         write (msg,'(a)') 'dpbp: iupcor failed'
         call warnings (trim(msg))
         if (log_tt) write (io_tt_log,'(a)') trim(msg)
         error = .true.
         return
      end if
      
      call givloc (elat, elon, bpdel, bpaz, bplat, bplon)
      
      return
         
   end subroutine dpbp
   
   
!*****************************************************************************************
   subroutine find_topo (xlat, xlon, top)
   
   ! Returns elevation in kilometers for given location in geographical coordinates.
   ! Elevation above sea level is taken positive, below sea level negative.
   
   ! Input:
   !  xlon,xlat = geographic longitude, latitude in degrees
   ! Output:
   !  top = topography in kilometers
   
   ! Adapted from Bob Engdahl's code 'findtopo' by EAB.
   
      integer :: i
      integer :: j
      integer :: k1
      integer :: k2
      integer :: m
      real :: a1
      real :: a2
      real :: res
      real, intent(out) :: top
      real :: topo1
      real :: topo2
      real :: topo3
      real :: topo4
      real, intent(in) :: xlat
      real :: xlat1
      real :: xlat2
      real, intent(in) :: xlon
      real :: xlon1
      real :: xlon2
      
      ! Look up indices and coordinates of grid points around the target location
      res = 60./5. ! Resolution 5 arc-sec
      i = int(res*(xlon + 180.))
      j = int(res*(90. - xlat))
      xlon1 = float(i)/res - 180.
      xlat1 = 90. - float(j)/res
      xlon2 = float(i + 1)/res - 180.
      xlat2 = 90. - float(j + 1)/res
      a1 = (xlon2 - xlon)/(xlon2 - xlon1)
      a2 = (xlat2 - xlat)/(xlat2 - xlat1)
    
      ! Take care of points outside the grid 
      k1 = i + 1
      k2 = i + 2
      m = j + 1
      if (i .le. 0 .or. i .ge. (n_bp_lon - 1)) then
          k1 = n_bp_lon
          k2 = 1
      end if
      if (j .le. 0) then
          m = 1
          a2 = 0.
      end if
      if (j .ge. (n_bp_lat - 1)) then
          m = n_bp_lat - 1
          a2 = 1.
      end if
    
      ! Bilinear interpolation 
      topo1 = float(bp_topo(k1,m))
      topo2 = float(bp_topo(k1,m+1))
      topo3 = float(bp_topo(k2,m))
      topo4 = float(bp_topo(k2,m+1))
      
      top = (1. - a1)*(1. - a2)*topo1 + a1*(1. - a2)*topo3 + (1. - a1)*a2*topo2 + a1*a2*topo4
      top = 1.0e-3*top ! convert to km
      
      return
      
   end subroutine find_topo
   
   
!*****************************************************************************************
   subroutine givloc (elat, elon, del, az, t1, p1)
   
   ! Returns the coordinates of the point which is a given distance and azimuth from a reference point
   
   ! input:
   !   elat,elon   = source latitude and longitude (geocentric coordinates, deg).
   !   del         = epicentral distance to bounce point (deg).
   !   az          = azimuth of bounce point (deg).         |
   ! output:
   !   t1 = bounce point latitude  ( + = N, - = S) in geographic coordinates (deg)
   !   p1 = bounce point longitude ( + = E, - = W) in geographic coordinates (deg)
   
   ! Adapted from Bob Engdahl's code 'givloc' by EAB.
   
      real, intent(in) :: az
      real :: azr
!      real :: bgeocen
      real :: colat
      real :: colon
      real :: cphi
      real :: ctheta
      real, intent(in) :: del
      real :: delr
      real, intent(in) :: elat
      real, intent(in) :: elon
!      real :: geogrf
      real, intent(out) :: p1
      real :: sphi
      real :: t0
      real, intent(out) :: t1
   
      ! conversion to epicentre colatitude and colongitude:
      colat = 90.0 - elat
      colon = elon
      if (elon .lt. 0.0) colon = colon + 360.0
   
      delr = del*rpd
      azr = az*rpd
      ! conversion to geocentric latitude:
      t0 = bgeocen(colat*rpd)
      t0 = elat*rpd
      ctheta = sin(delr)*sin(t0)*cos(azr) + cos(t0)*cos(delr)
      t1 = acos(ctheta)
   !    if (t0 .eq. 0.0) then ! original
      if (t0 .lt. 1.e-3) then ! eliminating the comparison of a real with zero
        p1 = az
      else if (t1 .lt. 1.e-3) then
        p1 = 0.0
      else
        sphi = sin(delr)*sin(azr)/sin(t1)
        cphi = (cos(delr) - cos(t0)*ctheta)/(sin(t0)*sin(t1))
        p1 = colon + atan2(sphi,cphi)*dpr
      end if
      ! convert colatitude to geographic latitude:
      !   assume p1 never > 720  
      t1 = 90.0 - geogrf(t1)*dpr 
      ! convert colongitude to longitude:
      if (p1 .gt. 360.0) p1 = p1 - 360.0   
      if (p1 .gt. 180.0) p1 = p1 - 360.0 
       
      return
         
   end subroutine givloc
   
   
!*****************************************************************************************
   real function bgeocen (arg)
   
   ! input:
   !   arg    = geographic colatitude (radians)
   ! output:
   !   bgeocen = geocentric colatitude (radians)
   ! (n.b. fac=(1-f)**2)
   
   ! Code by Bob Engdahl, cleaned up by EAB.
   
      real, intent(in) :: arg
      real :: fac = 0.993305621334896
      real :: pi2 = 1.570796326794895
   
      bgeocen = pi2 - atan(fac*cos(arg)/amax1(1.e-30,sin(arg)))
      
      return
         
   end function bgeocen
   
   
!*****************************************************************************************
   real function geogrf (arg)
   
   ! input:
   !   arg    = geocentric colatitude (radians)
   ! output:
   !   geogrf = geographic colatitude (radians
   ! (n.b. fac=(1-f)**2)
   
   ! Code by Bob Engdahl, cleaned up by EAB.
   
      real, intent(in) :: arg
      real :: fac = 0.993305621334896
      real :: pi2 = 1.570796326794895
   
      geogrf = pi2 - atan(cos(arg)/(fac*amax1(1.e-30,sin(arg))))
      
      return
      
   end function geogrf
   
   
!*****************************************************************************************
   subroutine topo_corr (phase, topo, dtddx, deltc, deltw)
   
   ! TOPography CORrection - returns the correction, in seconds, to be added
   ! to predicted pP, sP travel times.  
   ! The correction is for the bathymetry/elevation at the bounce point.  
   ! The ray parameter is used to calculate the incident angle at the surface, 
   ! so an accurate estimate of the ray length in the topography can be made.
   
   !  input:
   !    phase  = pP or sP
   !    topo = topography, positive above sealevel, in km
   !    rp2  = ray parameter [sec/km]
   !  output:
   !    deltc = crust travel time correction, s
   !    deltw = water travel time correction, s
   
   !  constants:
   !    vp   = p velocity at surface (km/s)
   !    vs   = s velocity at surface (km/s)
   !    vw   = water velocity (km/s)
   
   ! Adapted from Bob Engdahl's code 'topcor' by EAB.
   
      real :: bp2
      real, intent(out) :: deltc
      real, intent(out) :: deltw
      real, intent(in) :: dtddx
      real :: rp2
      real :: term
      real :: term1
      real :: term2
      real, intent(in) :: topo
      real, parameter :: vp = 5.80
      real, parameter :: vs = 3.46
      real, parameter :: vw = 1.5
      character(len=8), intent(in) :: phase
      character(len=132) :: msg
         
      rp2 = dtddx*dpr/radius
      bp2 = abs(rp2)
   
      if (abs(topo) .lt. 1.e-3) then ! No corrections needed
         deltc = 0.0
         deltw = 0.0
         return
      end if
            
      if (phase .eq. 'pP      ') then
         term = (vp*bp2)*(vp*bp2)
         if (term .gt. 1.0) term = 1.0
         deltc = 2.0*(topo/vp)*((1. - term)**0.5)
         if (topo .gt. 0.0) then ! Positive elevation, no water correction
            deltw = 0.0
         else
            term = (vw*bp2)*(vw*bp2)
            if (term .gt. 1.0) term = 1.0
            deltw = 2.0*(topo/vw)*((1. - term)**0.5)
            deltw = -deltw ! Water column correction is positive, added to pP time to get pwP time
         end if
      else if (phase .eq. 'sP      ') then
         term1 = (vp*bp2)*(vp*bp2)
         if (term1 .gt. 1.0) term1 = 1.0
         term2 = (vs*bp2)*(vs*bp2)
         if (term2 .gt. 1.0) term2 = 1.0
         deltc = (topo/vp)*((1. - term1)**0.5) + (topo/vs)*((1. - term2)**0.5)
         deltw = 0. ! No swP (rarely observed)
      else
         write (msg,'(a)') 'topo_corr: unsupported phase ('//trim(phase)//')'
         call warnings (trim(msg))
      end if
      
      ! Water layer correction is not done if water depth < 1.5 km
   !   if (topo .gt. -1.5) deltw = 0.
      
      return
      
   end subroutine topo_corr
   
   
!*****************************************************************************************
   subroutine rdp_depth_test (it)
   
   ! Trial over a range of depths to find the best fit of relative depth
   ! phases for each event that has relative depth phase data. At each depth,
   ! and for each relative depth phase time difference, the theoretical TT
   ! difference is calculated for all supported relative depth phase
   ! identifications: pP-P, sP-P and pwP-P. The minimum residual is selected,
   ! whether or not the phase names agree. The total RMS error from all
   ! selected residuals for each event is calculated. The preferred depth is
   ! the depth with  mnimum cumulative error. The asymmetric uncertainty in
   ! depth (at ~90% confidence level) is calculated using the zM test (e.g.,
   ! Langley, 'Practical Statistics, Simply Explained', Dover Publications,
   ! 1970, p 152.), to find the range of depths around the preferred depth in
   ! which the z statistic, relative to the value of z at the preferred
   ! depth, is less than the value at which the null hypothesis (no
   ! difference) would be rejected at the 10% level of probability.
   
      integer :: h
      integer :: h_minus
      integer :: h_plus
      integer :: iev
      integer :: ird
      integer, intent(in) :: it
      integer :: it1
      real :: azes_
      real :: azse_
      real :: cxlat
      real :: cxlon
      real :: delt_
      real :: depth_compare
      real :: depth1
      real :: depth2
      real :: depth3
      real :: dtpdh_
      real :: elcr_
      real :: hgtcr_
      real :: hrms_min
      real :: psr_
      real, dimension(n_event,n_hrdp_max) :: res_mean
      real :: res_min
      real :: res_min2
      real :: res1
      real :: res2
      real :: res3
      real :: sd
      real :: sxlat
      real :: sxlon
      real, dimension(3) :: ttcomp_
      real, dimension(2) :: usrc
      real :: xlat
      real :: xlon
      real :: z_90
      real :: z_min
      character(len=1) :: bdp_flag
      character(len=28) :: depth_pr
      character(len=132) :: msg
      character(len=8) :: phase_min
!      logical :: bdp
!      logical :: depth_constraint
      logical, dimension(n_event) :: has_rdp ! Event has relative depth phases
!      logical :: pass
      
      it1 = it + 1
      has_rdp = .false.
      hrmsi = 0.
      res_mean = 0.
      sd = 1. ! standard deviation of parent population
      n_samples = 0
      z_90 = 1.64 ! 10% probability of rejecting the null hypothesis
      bdp_flag = ' '
   
         
      if (verbose_log) write (io_depth_phase,'(a)') 'rdp_depth_test: Focal depth from the fit to relative depth phases:'
      do h = h_rdp1,h_rdp2 ! Depth loop
         do iev = 1,n_event ! Loop over events
            n_samples(iev) = 0
            do ird = 1,nst(iev) ! Loop over readings
               if (rel_depth_phase(iev,ird) .and. pass(fcode(iev,ird))) then
                  n_samples(iev) = n_samples(iev) + 1
                  has_rdp(iev) = .true.
                  if (indirect_cal) then
                    call geocen (latp_cal(iev), lonp_cal(iev), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
                    call depset (real(h), usrc)
                    call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                     real(h), ahgts(iev,ird), 'pP-P    ', ttcomp_(1), delt_,&
                     azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp_cal(iev), .false.)
                    call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                     real(h), ahgts(iev,ird), 'sP-P    ', ttcomp_(2), delt_,&
                     azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp_cal(iev), .false.)
                    call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                     real(h), ahgts(iev,ird), 'pwP-P   ', ttcomp_(3), delt_,&
                     azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp_cal(iev), .false.)
                  else
                    call geocen (latp(iev,it1), lonp(iev,it1), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
                    call depset (real(h), usrc)
                    call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                     real(h), ahgts(iev,ird), 'pP-P    ', ttcomp_(1), delt_,&
                     azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp(iev,it1), .false.)
                    call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                     real(h), ahgts(iev,ird), 'sP-P    ', ttcomp_(2), delt_,&
                     azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp(iev,it1), .false.)
                    call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                     real(h), ahgts(iev,ird), 'pwP-P   ', ttcomp_(3), delt_,&
                     azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp(iev,it1), .false.)
                  end if
                  res1 = tto(iev,ird) - ttcomp_(1)
                  res2 = tto(iev,ird) - ttcomp_(2)
                  res3 = tto(iev,ird) - ttcomp_(3)
                  res_min = min(abs(res1), abs(res2), abs(res3))
                  res_min2 = res_min*res_min
                  hrmsi(iev,h) = hrmsi(iev,h) + res_min2
                  res_mean(iev,h) = res_mean(iev,h) + res_min
                  if (verbose_log) then
                     if (bdp(iev,ird)) then ! Check for stations known to report bad depth phases
                        bdp_flag = '*'
                     else
                        bdp_flag = ' '
                     end if
                     write (io_depth_phase,'(i3,1x,2a,1x,a,1x,i5,1x,i3,4f10.2)') iev, bdp_flag,&
                      stname(iev,ird), phase(iev,ird), mnf_line(iev,ird), h, res1, res2, res3, res_min2
                  end if
               end if
            end do ! End loop over readings
            ! z-statistic for uncertainty of depth
            if (n_samples(iev) .gt. 0) then
               hrmsi(iev,h) = sqrt(hrmsi(iev,h)/real(n_samples(iev)))
               res_mean(iev,h) = res_mean(iev,h)/real(n_samples(iev))
               z_hrdp(iev,h) = sqrt(real(n_samples(iev)))*res_mean(iev,h)/sd
            end if
         end do ! End loop over events
      end do ! End loop over depths
      
      ! Find depth of minimum residual for each event
      do iev = 1,n_event
         if (.not.has_rdp(iev)) cycle
         h_rdp(iev,2) = -999
         hrms_min = huge(hrms_min)
         do h = h_rdp1,h_rdp2
            if (hrmsi(iev,h) .lt. hrms_min) then
               hrms_min = hrmsi(iev,h)
               h_rdp(iev,2) = h
            end if
            if (verbose_log) write (io_depth_phase,'(t6,i3,1x,i3,2e12.3,i3)') iev, h, hrmsi(iev,h), hrms_min, h_rdp(iev,2)
         end do
         if (verbose_log) write (io_depth_phase,'(a,i3,1x,a,1x,i4,a,i3,a)') 'Event ', iev, evtnam(iev),&
          h_rdp(iev,2), ' km on ', n_samples(iev), ' samples'
      end do
      
      ! Find range of depth uncertainty
      depth_compare = 1.0e-2
      do iev = 1,n_event
         if (.not.has_rdp(iev)) cycle
         z_min = z_hrdp(iev,h_rdp(iev,2))
         z_test(iev) = z_min + z_90
         h_rdp(iev,1) = h_rdp(iev,2)
         h_rdp(iev,3) = h_rdp(iev,2)
         do h = h_rdp(iev,2),h_rdp1,-1
            if (z_hrdp(iev,h) .lt. z_test(iev)) then
               h_rdp(iev,1) = h
            end if
         end do
         if (h_rdp(iev,1) .gt. h_rdp1) h_rdp(iev,1) = h_rdp(iev,1) - 1 ! Err on the conservative side
         do h = h_rdp(iev,2),h_rdp2
            if (z_hrdp(iev,h) .lt. z_test(iev)) then
               h_rdp(iev,3) = h
            end if
         end do
         if (h_rdp(iev,3) .lt. h_rdp2) h_rdp(iev,3) = h_rdp(iev,3) + 1 ! Err on the conservative side
         h_plus = h_rdp(iev,3) - h_rdp(iev,2)
         h_minus = h_rdp(iev,2) - h_rdp(iev,1)
         if (depth_cf(iev,1) .ge. 0.) then ! Depth settings taken from command file
            depth_pr = 'command file '//depset_pr(iev)
            depth1 = depth_cf(iev,1)
            depth2 = depth_cf(iev,2)
            depth3 = depth_cf(iev,3)
         else if (depth_hdf(iev,1) .ge. 0. .and. depth_constraint(depth_hdf_c(iev))) then ! Depth settings
          ! taken from a .hdf file if they are constrained
            depth_pr = 'hdf file '//depset_pr(iev)
            depth1 = depth_hdf(iev,1)
            depth2 = depth_hdf(iev,2)
            depth3 = depth_hdf(iev,3)
         else if (depth_inp(iev,1) .ge. 0. .and. depth_constraint(depth_inp_c(iev))) then ! Depth settings
          ! taken from the input data file if depth is constrained
            depth_pr = 'input file, constrained '//depset_pr(iev)
            depth1 = depth_inp(iev,1)
            depth2 = depth_inp(iev,2)
            depth3 = depth_inp(iev,3)
         else if (depth_default .ge. 0.) then ! Depth settings taken from the cluster default values
            depth_pr = 'default '//depset_pr(iev)
            depth1 = depth_default
            depth2 = 99.9
            depth3 = 99.9      
         else if (depth_inp(iev,1) .ge. 0.) then ! Depth settings taken from the input data file
            depth_pr = 'input file '//depset_pr(iev)
            depth1 = depth_inp(iev,1)
            depth2 = depth_inp(iev,2)
            depth3 = depth_inp(iev,3)           
         else ! This should not happen
            write (msg,'(a,i3,a)') 'rdp_depth_test: no depth has been set for event', iev, evtnam(iev)
            call oops (trim(msg))
         end if
         if (abs(h_rdp(iev,2) - depth1) .lt. depth_compare .and.&
             abs(h_plus - depth2) .lt. depth_compare .and.&
             abs(h_minus - depth3) .lt. depth_compare) then
            write (io_depth_phase,'(i3,1x,2a,i3,a,i3,a,i3,a,i3,a,3i3,a,i3,a)') iev, evtnam(iev),&
             ' preferred depth = ', h_rdp(iev,2), ' km (', h_rdp(iev,1), ' to ', h_rdp(iev,3),&
             ') on ', n_samples(iev), ' samples (depd ', h_rdp(iev,2), h_plus, h_minus, ' ! on ',&
             n_samples(iev), ' samples); current ('//trim(depth_pr)//') = same'      
         else
            write (io_depth_phase,'(i3,1x,2a,i3,a,i3,a,i3,a,i3,a,3i3,a,i3,a,3f5.1)') iev, evtnam(iev),&
             ' preferred depth = ', h_rdp(iev,2), ' km (', h_rdp(iev,1), ' to ', h_rdp(iev,3),&
             ') on ', n_samples(iev), ' samples (depd ', h_rdp(iev,2), h_plus, h_minus, ' ! on ',&
             n_samples(iev), ' samples); current ('//trim(depth_pr)//') = ', depth1, depth2, depth3
         end if
      end do
      
      ! Residuals against each theoretical phase for the preferred depth
      write (io_depth_phase,'(/a)') 'Residuals against each theoretical phase for the depth with minimum misfit'
      write (io_depth_phase,'(t25,a,t34,a,t44,a,t53,a,t63,a)') 'MNF', 'pP-P', 'sP-P', 'pwP-P', 'Err**2'
      do iev = 1,n_event ! Loop over events
         if (.not.has_rdp(iev)) cycle
         h = h_rdp(iev,2)
         write (io_depth_phase,'(a,i3,1x,a,1x,i4,a)') 'Event ', iev, evtnam(iev), h_rdp(iev,2), ' km'
         do ird = 1,nst(iev) ! Loop over readings, including ones that have been flagged (except duplicates)
            if (rel_depth_phase(iev,ird)) then
               if (fcode(iev,ird) .eq. 'd') cycle
               if (indirect_cal) then
                 call geocen (latp_cal(iev), lonp_cal(iev), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
                 call depset (real(h), usrc)
                 call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                  real(h), ahgts(iev,ird), 'pP-P    ', ttcomp_(1), delt_,&
                  azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp_cal(iev), .false.)
                 call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                  real(h), ahgts(iev,ird), 'sP-P    ', ttcomp_(2), delt_,&
                  azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp_cal(iev), .false.)
                 call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                  real(h), ahgts(iev,ird), 'pwP-P   ', ttcomp_(3), delt_,&
                  azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp_cal(iev), .false.)
               else
                 call geocen (latp(iev,it1), lonp(iev,it1), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
                 call depset (real(h), usrc)
                 call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                  real(h), ahgts(iev,ird), 'pP-P    ', ttcomp_(1), delt_,&
                  azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp(iev,it1), .false.)
                 call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                  real(h), ahgts(iev,ird), 'sP-P    ', ttcomp_(2), delt_,&
                  azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp(iev,it1), .false.)
                 call mlocsteq (iev, ird, xlat, xlon, stladg(iev,ird), stlndg(iev,ird),&
                  real(h), ahgts(iev,ird), 'pwP-P   ', ttcomp_(3), delt_,&
                  azes_, azse_, psr_, dtpdh_, elcr_, hgtcr_, latp(iev,it1), .false.)
               end if
               res1 = tto(iev,ird) - ttcomp_(1)
               res2 = tto(iev,ird) - ttcomp_(2)
               res3 = tto(iev,ird) - ttcomp_(3)
               res_min = min(abs(res1), abs(res2), abs(res3))
               if (abs(res_min - abs(res1)) .lt. 1.0e-3) then
                  phase_min = 'pP-P    '
               else if (abs(res_min - abs(res2)) .lt. 1.0e-3) then
                  phase_min = 'sP-P    '
               else if (abs(res_min - abs(res3)) .lt. 1.0e-3) then
                  phase_min = 'pwP-P   '
               else
                  phase_min = 'UNKNOWN '
               end if
               res_min2 = res_min*res_min
               if (bdp(iev,ird)) then ! Check for stations known to report bad depth phases
                  bdp_flag = '*'
               else
                  bdp_flag = ' '
               end if
               if (phase(iev,ird) .eq. phase_min) then
                  write (io_depth_phase,'(1x,2a,f6.2,1x,a,i5,4f10.2,1x,a)') bdp_flag, stname(iev,ird), delt_,&
                   phase(iev,ird), mnf_line(iev,ird), res1, res2, res3, res_min2, fcode(iev,ird)
               else
                  write (io_depth_phase,'(1x,2a,f6.2,1x,a,i5,4f10.2,2(1x,a))') bdp_flag, stname(iev,ird), delt_,&
                   phase(iev,ird), mnf_line(iev,ird), res1, res2, res3, res_min2, fcode(iev,ird), phase_min
                  ! Output file for depth phase name changes
                  write (io_dpnc,'(i3,1x,a5,1x,a8,1x,a8,1x,i5,1x,a8,1x,a)') iev, stname(iev,ird),&
                   readsrc(iev,ird), phase0(iev,ird), mnf_line(iev,ird), phase_min, infile20(iev)
               end if
            end if
         end do ! End loop over readings
      end do ! End loop over events
      
      return
      
   end subroutine rdp_depth_test


!*****************************************************************************************
   subroutine phase_alias (phase, delta, ip)
         
   ! Check for alternative phase names
   ! Input: phase, delta
   ! Output: ip (index of phcod)
   
      integer, intent(out) :: ip
      real, intent(in) :: delta
      character(len=*), intent(in) :: phase
   
      if (phase(1:3).eq.'Pn ') then ! phase='P       '
        ip=2
      else if (phase(1:3).eq.'Sn ') then ! phase='S       '
        ip=33
      else if (phase(1:4).eq.'pPn ') then ! phase='pP      '
        ip=8
      else if (phase(1:4).eq.'pwP ') then ! phase='pP      '
        ip=8
      else if (phase(1:5).eq.'pwPn ') then ! phase='pP      '
        ip=8
      else if (phase(1:4).eq.'sPn ') then ! phase='sP      '
        ip=13
      else if (phase(1:4).eq.'pSn ') then ! phase='pS      '
        ip=37
      else if (phase(1:4).eq.'sSn ') then ! phase='sS      '
        ip=40
      else if (phase(1:4).eq.'SPn ') then ! phase='SP      '
        ip=55
      else if (phase(1:4).eq.'SnP ') then ! phase='SP      '
        ip=55
      else if (phase(1:4).eq.'PSn ') then ! phase='PS      '
        ip=56
      else if (phase(1:5).eq.'PnPn ') then ! phase='PP      '
        ip=30
      else if (phase(1:5).eq.'SnSn ') then ! phase='SS      '
        ip=53
      else if (phase(1:2).eq.'p ') then ! phase='Pup     '
        ip=1  
      else if (phase(1:2).eq.'s ') then ! phase='Sup     '
        ip=32 
      else if (phase(1:6).eq."P'P'ab") then ! phase="P'P'    '
        ip=31
      else if (phase(1:6).eq."P'P'bc") then ! phase="P'P'    '
        ip=31
      else if (phase(1:6).eq."P'P'df") then ! phase="P'P'    '
        ip=31
      else if (phase(1:6).eq."S'S'ac") then ! phase="S'S'    '
        ip=54
      else if (phase(1:6).eq."S'S'df") then ! phase="S'S'    '
        ip=54
      else if (delta.le.100.0.and.phase.eq.'pPdiff  ') then ! phase='pP      '
        ip=8
      else if (delta.le.100.0.and.phase.eq.'pwPdiff ') then ! phase='pP      '
        ip=8
      else if (delta.le.100.0.and.phase.eq.'sPdiff  ') then ! phase='sP      '
        ip=13
      else if (delta.le.100.0.and.phase.eq.'pSdiff  ') then ! phase='pS      '
        ip=37
      else if (delta.le.100.0.and.phase.eq.'sSdiff  ') then ! phase='sS      '
        ip=40
      else if (delta.le.165.0.and.phase.eq.'PKPdiff ') then ! phase='PKPbc   '
        ip=5
      else if (delta.le.165.0.and.phase.eq.'pPKPdiff') then ! phase='pPKPbc  '
        ip=10
      else if (delta.le.165.0.and.phase.eq.'sPKPdiff') then ! phase='sPKPbc  '
        ip=15
      else
        ip=-1
      end if
      
      return
      
   end subroutine phase_alias
   
   
!*****************************************************************************************
end module mloc_tt

