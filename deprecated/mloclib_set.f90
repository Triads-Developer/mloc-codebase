!> Procedures related to mlocset, the main driver for the relocation

module mloclib_set

   use declare_limits
   use declare_calibration
   use declare_cluster_vectors
   use declare_configuration
   use declare_events
   use declare_hdc
   use declare_hypocentroid
   use declare_lun
   use declare_output
   use declare_phase_data
   use declare_phases
   use declare_stations
   use mloclib_geog
   use mloclib_messages
   use mloclib_statistics

   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine model_init (mt)
         
   ! Initialize location parameters and text strings for output.
   ! Convert coordinates to geocentric.
   ! Initial estimate of hypocentroid.
   ! Number of free parameters for each event and total.
   ! Initialize bias correction term.
   ! Initialize ellipticity function subroutine.
         
      integer :: i
      integer :: iev
      integer, intent(out) :: mt ! Total number of free parameters
      real :: cxlat
      real :: cxlon
      real :: ellip
      real, external :: hms2s
      real :: openaz
      real :: openazlim
      real :: sxlat
      real :: sxlon
      real :: tcor
      real :: x1sum
      real :: x2sum
      real :: x3sum
      real :: x4sum
      real :: xlat
      real :: xlon
      character(len=60) :: depth_set1
      character(len=60) :: depth_set2
     
      openazlim = 300.
      
      if (verbose_screen) write (*,'(/a/)') 'Starting model and free parameters:'
      
      do iev = 1,nev
      
         call azgap (iev, openaz)
         if (verbose_screen) write (*,'(i3,2(2x,a),a,i3)') iev, trim(evtnam(iev)), trim(infile(iev)),&
          ' / Open azimuth = ', nint(openaz)
         
         latpr(iev) = '   '
         lonpr(iev) = '   '
         depthpr(iev) = '     '
         timepr(iev) = '           '
         
         ! Latitude free or fixed
         if (latf(iev) .and. (openaz .lt. openazlim)) then ! Free
            mindx(iev,1) = 1
            latpr(iev) = 'lat'
         else ! Fixed
            if (verbose_screen) write (*,'(t5,a,f7.3)') 'latitude fixed at ', latp(iev,0)
            mindx(iev,1) = 0
         end if
         
         ! Longitude free or fixed
         if (lonf(iev) .and. (openaz .lt. openazlim)) then ! Free
            mindx(iev,2) = 1
            lonpr(iev) = 'lon'
         else ! Fixed
            if (verbose_screen) write (*,'(t5,a,f8.3)') 'longitude fixed at ', lonp(iev,0)
            mindx(iev,2) = 0
         end if
                  
   !  Make sure all longitudes are given in the same range.
   !      if (iev .gt. 1) then
   !         diflon = lonp(iev,0) - lonp(iev-1,0)
   !         if (diflon .gt. 100.) then
   !            write (msg,'(a,f5.1,a,i3,1x,a,2f7.1)') 'model_init: longitude range violation: difference = ', diflon,&
   !             '° for event ', iev, evtnam(iev), lonp(iev,0), lonp(iev-1,0)
   !            call oops (trim(msg))
   !         end if
   !      end if
         
         ! Origin time free or fixed
         if (timef(iev)) then ! Free
            mindx(iev,4) = 1
            timepr(iev) = 'origin time'
         else ! Fixed
            mindx(iev,4) = 0
         end if
         otsp(iev,0) = hms2s(hourp(iev,0), minp(iev,0), secp(iev,0))
         
         ! Depth
         if (depset_pr(iev) .eq. 'c') then
            depth_set1 = ', from cluster default depth'
         else if (depset_pr(iev) .eq. 'd') then
            depth_set1 = ', from depth phases'
         else if (depset_pr(iev) .eq. 'e') then
            depth_set1 = ', from engineering information'
         else if (depset_pr(iev) .eq. 'f') then
            depth_set1 = ', from fault modeling (e.g., InSAR, GPS)'
         else if (depset_pr(iev) .eq. 'i') then
            depth_set1 = ', from input data file preferred hypocenter'
         else if (depset_pr(iev) .eq. 'l') then
            depth_set1 = ', from local-distance readings'
         else if (depset_pr(iev) .eq. 'm') then
            depth_set1 = ', from mloc solution with free depth'
         else if (depset_pr(iev) .eq. 'n') then
            depth_set1 = ', from near-source readings'
         else if (depset_pr(iev) .eq. 'r') then
            depth_set1 = ', from relocation with free depth outside mloc'
         else if (depset_pr(iev) .eq. 'u') then
            depth_set1 = ', unconstrained'
         else if (depset_pr(iev) .eq. 'w') then
            depth_set1 = ', from waveform analysis'
         else
            depth_set1 = ', from unknown source with code '//depset_pr(iev)
         end if
         if (hdepthshift .lt. 0.01) then
            depth_set2 = ' '
         else
            depth_set2 = ', perturbed in command file'
         end if
         ! Free or fixed
         if (depthf(iev)) then ! Free
            mindx(iev,3) = 1
            depthpr(iev) = 'depth'
         else ! Fixed
            if (verbose_screen) write (*,'(t5,a,f5.1,2a)') 'depth fixed at ', depthp(iev,0), trim(depth_set1), trim(depth_set2)
            mindx(iev,3) = 0
         end if
         
      end do      
   
      ! Convert event coordinates to geocentric.
      do iev = 1,nev
         call geocen (latp(iev,0), lonp(iev,0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
         latpgc(iev,0) = xlat
         lonp(iev,0) = xlon
      end do
   
      ! Initial estimate of hypocentroid
      x1sum = 0.
      x2sum = 0.
      x3sum = 0.
      x4sum = 0.
      do iev = 1,nev
         x1sum = x1sum + latp(iev,0)
         x2sum = x2sum + lonp(iev,0)
         x3sum = x3sum + depthp(iev,0)
         x4sum = x4sum + otsp(iev,0)
      end do
      lath(0) = x1sum/real(nev)
      lonh(0) = x2sum/real(nev)
      call geocen (lath(0), lonh(0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
      lathgc(0) = xlat
      depthh(0) = x3sum/real(nev)
      otsh(0) = x4sum/real(nev)
      call timecr (otsh(0), hourh(0), minh(0), sech(0)) 
   
      ! Number of free parameters for each event and total.
      mt = 0
      do iev = 1,nev
         mtiev(iev) = 0
         do i = 1,4
            if (mindx(iev,i) .ne. 0) then
               mtiev(iev) = mtiev(iev) + 1
               mindx(iev,i) = mtiev(iev)
            end if
         end do
         mt = mt + mtiev(iev)
      end do
   
      ! Initialize bias correction term
      bcorr = 0.
      
      ! Initialize ellipticity function subroutine
      if (verbose_screen) write (*,'(/a/)') 'Initialize ellip'
      tcor = ellip ('P       ', 30., 30., 30., 30., .true.)
      
      return
      
   end subroutine model_init
      
   
!*****************************************************************************************
   subroutine stflt (iev, it)
   
   ! Completely rewritten 3/3/06 by eab.
   ! Logical variables are defined for each reading to determine if it passes criteria on
   ! duplicate readings, size of residual, and flags for both the cluster vector and the hypocentroid.
   
   ! .TRUE. for one of the logical variables means the reading has been filtered for that reason,
   ! i.e., it should not be used.
   
   ! Epicentral distance tests are based on the distance ranges defined by commands CLIM and HLIM.
   ! Residual size tests are based on the window defined by command WIND. In the range used for
   ! direct calibration (defined by windloclim) the residual limits for filtering are set explicitly
   ! and they are rather large because local-distance readings can have large residulas if the
   ! starting location is poor.
   ! If a reading is flagged, it is not counted for epicentral distance or residual problems.
   
      integer :: i
      integer, intent(in) :: iev ! event number
      integer :: ird
      integer, intent(in) :: it ! iteration index
      real :: dtmp
      real :: dts
      real :: rtmp
      character(len=160) :: msg
   
      ndat(iev,it) = 0
      ndatdl(iev,it) = 0
      ndatpr(iev,it) = 0
      ndatfl(iev) = 0
      do ird = 1,nst(iev)
         dtmp = delt(iev,ird)
         rtmp = wind2*ttsprd(iev,ird)
         if (direct_cal .and. dtmp .le. windloclim) rtmp = rtmp * 2.
         
         ! Filter for flags
         fltrcflag(iev,ird) = .false.
         fltrhflag(iev,ird) = .false.
         if (fcode(iev,ird) .eq. 'x') then ! Outlier (large residual), either cluster or absolute.
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         if (fcode(iev,ird) .eq. 'p') then ! The phase is problematic.
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         if (fcode(iev,ird) .eq. 'd') then ! Duplicate reading.
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         if (fcode(iev,ird) .eq. 's') then ! Flagged on the basis of phase name, station and/or author by the SKIP command.
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         if (fcode(iev,ird) .eq. 'm') then ! Missing station coordinates
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         if (fcode(iev,ird) .eq. 't') then ! Timing is suspect
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         if (fcode(iev,ird) .eq. '*') then ! Engdahl's flag.
            fltrcflag(iev,ird) = .true.
            fltrhflag(iev,ird) = .true.
         end if
         
         ! Filter for epicentral distance
   !     if (fltrcflag(iev,ird)) then
   !        fltrcdelt(iev,ird) = .false.
   !     else
   !        fltrcdelt(iev,ird) = .true.
   !     end if
         fltrcdelt(iev,ird) = .true.
         do i = 1,nclim
            if (dtmp .ge. clim(i,1) .and. dtmp .le. clim(i,2)) fltrcdelt(iev,ird) = .false.
         end do
   !     if (fltrhflag(iev,ird)) then
   !        fltrhdelt(iev,ird) = .false.
   !     else
   !        fltrhdelt(iev,ird) = .true.
   !     end if
         fltrhdelt(iev,ird) = .true.
         do i = 1,nhlim
            if (dtmp .ge. hlim(i,1) .and. dtmp .le. hlim(i,2)) fltrhdelt(iev,ird) = .false.
         end do
         
         ! Filter for residual.
         dts = dt(iev,ird,it) - s(iev,ird,it) - ttoff(iev,ird)
         fltrcres(iev,ird) = .true.
         fltrhres(iev,ird) = .true.
         if (abs(dts) .le. rtmp) then
            fltrcres(iev,ird) = .false.
            fltrhres(iev,ird) = .false.
         end if
         
         
         fltrc(iev,ird) = (fltrcdelt(iev,ird) .or. fltrcres(iev,ird) .or. fltrcflag(iev,ird))               
         fltrh(iev,ird) = (fltrhdelt(iev,ird) .or. fltrhres(iev,ird) .or. fltrhflag(iev,ird))
         
         ! Extra filtering when using only P arrivals for hypocentroid.
         if (ponly .and. .not.fltrh(iev,ird)) then
            if (phase(iev,ird) .ne. 'P       ' .and.&
                phase(iev,ird) .ne. 'Pg      ' .and.&
                phase(iev,ird) .ne. 'Pb      ' .and.&
                phase(iev,ird) .ne. 'Pn      ' .and.&
                phase(iev,ird) .ne. 'S-P     ' .and.&
                .not.(depthfh .and. phase(iev,ird) .eq. 'pP      ') .and.&
                .not.(depthfh .and. phase(iev,ird) .eq. 'sP      ')) fltrh(iev,ird) = .true.
         end if
         
         ! Filter for differential time data (never used for hypocentroid).
         if (idiff(iev,ird) .gt. 0) fltrh(iev,ird) = .true.
         
         ! Update counters
         if (.not.fltrh(iev,ird)) ndat(iev,it) = ndat(iev,it) + 1
         if (fltrhdelt(iev,ird)) ndatdl(iev,it) = ndatdl(iev,it) + 1
         if (fltrhres(iev,ird)) ndatpr(iev,it) = ndatpr(iev,it) + 1
         if (fltrhflag(iev,ird) .and. fcode(iev,ird) .ne. 'd') ndatfl(iev) = ndatfl(iev) + 1
         
         if (debug) write (io_log,'(a,i3,i5,1x,a6,1x,a8,1x,f7.2,1x,f7.2,1x,a1,1x,4l1,1x,4l1)') 'stflt: ',&
             iev, ird, stname(iev,ird), phase(iev,ird), dtmp, rtmp, fcode(iev,ird),&
             fltrc(iev,ird), fltrcdelt(iev,ird), fltrcres(iev,ird), fltrcflag(iev,ird),&
             fltrh(iev,ird), fltrhdelt(iev,ird), fltrhres(iev,ird), fltrhflag(iev,ird)
         
      end do
      
      return
      
   end subroutine stflt
   
   
!*****************************************************************************************
   subroutine hdf_read ()
   
   ! To set the starting locations for a run, the final locations from a previous run can be
   ! read from an HDF file. This routine reads the HDF file. It will not be compatible with
   ! very old HDF files because of changes to the format.
      
      integer :: hdf_hour
      integer :: hdf_min
      integer :: ios
      integer :: nhdf
      real :: hdf_sec
      real, external :: hms2s
      character(len=146) :: linein
      character(len=132) :: msg
      logical :: loop
      
      call open_file (io_rhdf, rhdf_filnam, 'old')
      write (io_log,'(/a/a)') 'Starting locations read from:', trim(rhdf_filnam)
      nhdf = 0
      loop = .true.
      do while (loop)
         read (io_rhdf,'(a)',iostat=ios) linein
         if (ios .lt. 0) exit
         write (io_log,'(a)') linein(1:50)
         nhdf = nhdf + 1
         if (nhdf .gt. nhdfmax) then
            call warnings ('hdf_read: nhdfmax exceeded')
            nhdf = nhdf - 1
            exit
         end if
         ! From HDF format v9.8.2 on
         read (linein(1:4),'(i4)') hdf_yr4(nhdf)
         read (linein(6:7),'(i2)') hdf_mon(nhdf)
         read (linein(9:10),'(i2)') hdf_day(nhdf)
         read (linein(12:13),'(i2)') hdf_hour
         read (linein(15:16),'(i2)') hdf_min
         read (linein(18:22),'(f5.2)') hdf_sec
         read (linein(24:32),'(f9.5)') hdf_lat(nhdf)
         read (linein(34:43),'(f10.5)') hdf_lon(nhdf)
         read (linein(45:50),'(f6.2)') hdf_dep(nhdf)
         hdf_dep_code(nhdf) = linein(52:53)
         hdf_evid(nhdf) = linein(67:76)
         read (linein(106:109),'(f4.1)') hdf_dep_plus(nhdf)
         read (linein(111:114),'(f4.1)') hdf_dep_minus(nhdf)
         if (debug) write (io_log,'(a,i3,f6.1,1x,a2,2f6.1)') 'hdf_read: ', nhdf, hdf_dep(nhdf),&
          hdf_dep_code(nhdf), hdf_dep_plus(nhdf), hdf_dep_minus(nhdf)
         hdf_time(nhdf) = hms2s(hdf_hour, hdf_min, hdf_sec)
      end do
      nrhdf = nhdf
      if (verbose_screen) write (*,'(/i3,a/a)') nrhdf, ' starting locations read from:', trim(rhdf_filnam)
      
      close (io_rhdf)
      
      return
      
   end subroutine hdf_read
   
   
!*****************************************************************************************
   subroutine read_bdp ()
   
   ! Read an optional file with a listing of stations suspected of reporting bogus depth phases.
   ! This information will be reflected in the type 6 ('_tt6') plot of depth phases, in which these
   ! suspected stations will be plotted with a smaller symbol. Such readings will also be flagged
   ! with "*" in output files dealing with depth phases.
   
      integer :: ios
      integer :: isstn
      character(len=108) :: linein
      character(len=132) :: msg
      
      write (io_log,'(/2a)') 'Bogus Depth Phase Stations read from: ', trim(bdp_filnam)
      call open_file (io_bdp, bdp_filnam, 'old')
      
      n_bdp = 0
      
      ! As of release date 1/15/2019 (v10.4.7) the generic station file format (type 3) must be used.
      read (io_bdp,'(i1)') isstn
      if (isstn .ne. 3) then
         call warnings ('read_bdp: incorrect file format (must be type 3, generic station file format)')
         return
      end if
         
      do
         read (io_bdp,'(a)',iostat=ios) linein
         if (ios .lt. 0) exit
         if (n_bdp .lt. n_bdp_max) then
            n_bdp = n_bdp + 1
            bdp_station(n_bdp) = linein
            write (io_log,'(a)') linein
         else
            write (msg,'(a,i3)') 'read_bdp: maximum number of entries read: ', n_bdp_max
            call warnings (trim(msg))
            exit
         end if
      end do
      if (verbose_screen) write (*,'(i3,2a)') n_bdp, ' entries read from ', trim(bdp_filnam)
      
      close (io_bdp)
      
      return
      
   end subroutine read_bdp
   
   
!*****************************************************************************************
   subroutine depth_spread (it)
   
   ! Robust estimate of spread of constrained focal depths
   
      real, dimension(nevmax) :: cfd
      integer :: iev
      integer, intent(in) :: it ! iteration index
      character(len=132) :: msg
      
      do iev = 1,nev
         if (depset_pr(iev) .ne. 'c' .and. depset_pr(iev) .ne. 'u' .and. depset_pr(iev) .ne. 'i' .and.&
             depset_pr(iev) .ne. ' ') then
            n_cfd = n_cfd + 1
            cfd(n_cfd) = depthp(iev,it+1)
         end if
      end do
      
      if (n_cfd .ge. 3) then
         call croux (cfd, n_cfd, rescfd)
         write (msg,'(a,f4.1,a,i3,a)') 'depth_spread: spread of constrained focal depths = ', rescfd,&
          ' on ', n_cfd, ' samples'
      else
         msg = 'depth_spread: too few samples'
      end if
      
      call fyi (trim(msg))
      write (io_log,'(/a)') trim(msg)
      
      return
      
   end subroutine depth_spread
   

!*****************************************************************************************
end module mloclib_set

