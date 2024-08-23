!> The main routine that runs mloc

module mloc_set

   use mloc_calibration
   use mloc_declare
   use mloc_inv
   use mloc_lib
   use mloc_math
   use mloc_mnf
   use mloc_out
   use mloc_phases
   use mloc_plots
   use mloc_stations
   use mloc_taup
   use mloc_tt
   
   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine mlocset ()
      
   ! Set up location problem: establish free parameters and compute partial derivatives and residuals.
   ! Call inversion, update parameters, and iterate until convergence. Make calibration shift if
   ! necessary. Produce all requested output files and plots.
   
   ! January 26, 1989 by eric bergman, numerous changes since...
      
      integer :: atsp_hr
      integer :: atsp_min
      integer :: i
      integer :: ic
      integer :: id0
      integer :: ielev_m
      integer :: iev
      integer :: ip
      integer :: ird
      integer :: it
      integer :: itlim_flag
      integer :: j
      integer :: k
      integer :: k0
      integer :: l
      integer :: mt
      integer :: nphreid
      integer :: nr
      real :: arrtime
      real :: atsp
      real :: atsp_sec
      real :: cxlat
      real :: cxlon
      real :: ddep
      real :: depth_max
      real :: depth_min
      real :: dlat
      real :: dlon
      real :: dtim
      real :: dtpdh
      real :: dts
      real :: glat
      real :: glon
      real :: hgtcr
      real :: minimum_weight
      real :: p1
      real :: p2
      real :: psr
      real :: res_test
      real :: rlatdg1
      real :: rlondg1
      real :: sxlat
      real :: sxlon
      real :: theta
      real :: ttobs
      real, dimension(2) :: usrc
      real :: x1sum
      real :: x2sum
      real :: x3sum
      real :: x4sum
      real :: xlat
      real :: xlon
      character(len=160) :: command_line
      character(len=30) :: blank30
      character(len=5) :: blank5
      character(len=132) :: dum1
      character(len=132) :: dum2
      character(len=132) :: file_folder
      character(len=132) :: msg
      character(len=100) :: outfil
      character(len=100) :: outfildat0 
      character(len=8) :: phtest
      logical :: confidence_ellipse_plot
      logical :: convrg
      logical :: ex
      
      data itlim_flag /2/ ! Iteration limit for changing filter flags
      data minimum_weight /0.01/
      blank5 = ' '
      blank30 = ' '
            
      ! Output file for .dat0 file (event data as read in), only in cluster mode
      outfildat0 = trim(outfile)//'_dat0.mnf'
      call open_file (io_dat0, outfildat0, 'new')
      write (io_dat0,'(a,3x,a,t125,a)') 'B', 'original disposition of '//trim(outfile), ' '
      write (io_dat0,'(a,t125,a)') 'F   MNF v1.3.4', ' '
      
      ! Read station file.
      write (*,'(/a)') 'Station coordinates...'
      call redsta 
      
      ! Read travel time and ellipticity correction datasets
      write (*,'(/a)') 'Travel time tables...'
      call redtab
      
      ! Travel time spread for individual phases
      write (*,'(/a)') 'Travel times spreads...'
      call rdttsprd ()
      
      ! Reading errors from previous run or use default values.
      write (*,'(/a)') 'Reading errors...'
      call getrderr ()
      
      ! Get starting locations from an HDF file
      if (read_hdf) then
         write (*,'(/a)') 'Starting locations from an HDF file...'
         call hdf_read ()
      end if
      
      ! List of stations suspected of reporting bogus depth phases
      if (bdp_list) then
         write (*,'(/a)') 'Stations suspected of reporting bogus depth phases...'
         call read_bdp ()
      end if
            
      ! For each event, read epicenter information and phase data.
      write (*,'(/a)') 'Current MNF format: 1.3.4'
      write (*,'(/a)') 'Phase data for each event...'
      write (io_log,'(/a)') 'Phase data for each event...'
      if (star_auto) n_star = 0 ! Counter for auto-plotting stars based on magnitude
      n_miss_sta_total = 0
      do iev = 1,n_event
         write (msg,'(a,i3,2(2x,a))') 'mlocset: event ', iev, trim(evtnam(iev)), trim(infile(iev))
         write (io_log,'(a)') trim(msg)
         write (*,'(a)') trim(msg)
         call open_file (io_in, infile(iev), 'old')
         call read_mnf (iev, io_in)
         close (io_in) 
      end do
      write (io_dat0,'(a,t110,a)') 'EOF', ' '
      close (io_dat0)
      
      if (nsmd) then
         write (msg,'(a,i8,2a)') 'mlocset: ', n_nsmd_written, ' NEIC station metadata entries written to ',&
          trim(nsmd_stn_file)
         call fyi (trim(msg))
         call warnings ('mlocset: A supplemental station file must be referenced to use these entries')
         call nsmd_deallocate () ! Deallocate arrays, we don't need them anymore
         close (io_nsmd)
      end if
      
      if (ifsm) then
         write (msg,'(a,i8,2a)') 'mlocset: ', n_ifsm_written, ' ISC-FDSN station metadata entries written to ',&
          trim(ifsm_stn_file)
         call fyi (trim(msg))
         call warnings ('mlocset: A supplemental station file must be referenced to use these entries')
         call ifsm_deallocate () ! Deallocate arrays, we don't need them anymore
         close (io_ifsm)
      end if

      ! Differential time data
      if (diffdat) then
            call open_file (io_diffdat, diffdatfilnam, 'old')
            call read_mnf (0, io_diffdat)
            close (io_diffdat)
      end if
      
      ! For indirect calibration if calibration data have been given in both the input file and the command file.
      ! Command file values take precedence.
      if (indirect_cal) then
         do iev = 1,n_event
            if (cal_event(iev,1)) then ! Values from command file
               cal_lat(iev,3) = cal_lat(iev,1)
               cal_lon(iev,3) = cal_lon(iev,1)
               cal_dep(iev,3) = cal_dep(iev,1)
               cal_hr(iev,3) = cal_hr(iev,1)
               cal_min(iev,3) = cal_min(iev,1)
               cal_sec(iev,3) = cal_sec(iev,1)
               rcv(iev,3,1,1) = rcv(iev,1,1,1)
               rcv(iev,3,1,2) = rcv(iev,1,1,2)
               rcv(iev,3,2,1) = rcv(iev,1,2,1)
               rcv(iev,3,2,2) = rcv(iev,1,2,2)
               rcv(iev,3,3,3) = rcv(iev,1,3,3)
               rcv(iev,3,4,4) = rcv(iev,1,4,4)
               cal_event(iev,3) = .true.
               ncal(3) = ncal(3) + 1
            else if (cal_event(iev,2)) then ! Values from input file
               cal_lat(iev,3) = cal_lat(iev,2)
               cal_lon(iev,3) = cal_lon(iev,2)
               cal_dep(iev,3) = cal_dep(iev,2)
               cal_hr(iev,3) = cal_hr(iev,2)
               cal_min(iev,3) = cal_min(iev,2)
               cal_sec(iev,3) = cal_sec(iev,2)
               rcv(iev,3,1,1) = rcv(iev,2,1,1)
               rcv(iev,3,1,2) = rcv(iev,2,1,2)
               rcv(iev,3,2,1) = rcv(iev,2,2,1)
               rcv(iev,3,2,2) = rcv(iev,2,2,2)
               rcv(iev,3,3,3) = rcv(iev,2,3,3)
               rcv(iev,3,4,4) = rcv(iev,2,4,4)            
               cal_event(iev,3) = .true.
               ncal(3) = ncal(3) + 1
            end if
         end do
      end if
      
      ! Readings that failed the station file date range
      if (n_failed_date_range .gt. 0) then
         write (msg,'(a,i4,a)') 'mlocset: ', n_failed_date_range, ' readings failed the station file date range'
         call fyi (trim(msg))
         write (io_stn_log,'(a)') trim(msg)
      end if
      
      ! Summary of missing stations
      if (n_miss_sta_total .gt. 0) then
         write (msg,'(a,i4,a)') 'mlocset: ', n_miss_sta_total, ' codes missing from the station files'
         call fyi (trim(msg))
         write (io_stn_log,'(/a)') 'Stations in the dataset that are missing from the station files:'
         do i = 1,n_miss_sta_total
            write (io_stn_log,'(a,1x,i3,a)') n_miss_sta_list(i), n_miss_sta(i), ' instances'
         end do
      end if
      
      ! Stations from supplemental station files that are present in the dataset
      call supp_stations_used ()
      
      ! Stations from supplemental station files that are not present in the dataset
      call supp_stations_not_used ()
      
      ! Stations used
      write (io_stn_log,'(/i6,a)') nkstat, ' stations in the arrival time dataset for which coordinates were found:'
      do i = 1,nkstat
         call geogra (stalat(kstat(i)), glat) ! Geographic latitude
         glon = stalon(kstat(i))
         call set_longitude_range (glon, 0)
         ielev_m = int(ahgtr(kstat(i))*1.0e3) ! Elevation in m
         write (msg,'(2a)') sta_author(kstat(i)), trim(duplication(kstat(i)))
         call supp_station_file_write (io_stn_log, nstr1(kstat(i)), sta_agency(kstat(i)), sta_deployment(kstat(i)),&
          glat, glon, ielev_m, 0, jdate_on(kstat(i)), jdate_off(kstat(i)), trim(msg))
      end do
      
      ! Initialize
      call model_init (mt)
      
      ! Output file for depth phase information
      if (verbose_screen) write (*,'(/a)') 'Depth phases log file'
      outfil = trim(outfile)//'.depth_phases'
      call open_file (io_depth_phase, outfil, 'new')
      
      ! Output file for depth phase name changes
      if (verbose_screen) write (*,'(/a)') 'Depth phase name change log file'
      outfil = trim(outfile)//'.dpnc'
      call open_file (io_dpnc, outfil, 'new')
      
      if (bptc) write (io_bptc_log,'(t3,a,t9,a,t20,a,t30,a,t39,a,t44,a,t50,a,t59,a,t68,a,t78,a,t87,a)')&
       'STA', 'PHASE', 'BP_LAT', 'BP_LON', 'TOPO', 'CRUST', 'WATER', 'T', 'Tobs', 'RES', 'MNF'
       
!***********************************************************************
      !  location iteration loop starts here                                  
      it = 0
      convrg = .false.
      
      do while (it .le. nstep .and. .not.convrg) ! Iteration loop
         write (*,'(/a,i1,a)') ' Beginning iteration ', it, '...'
         write (io_log,'(/a)')       '**************************************************'
         write (io_log,'(a,i1,a)') ' Beginning iteration ', it, '...'
         write (io_log,'(/a)') 'Hypocentroid readings from each event'
         write (io_log,'(a)') 'iev ndat    ndatdl  ndatpr'
         if (bptc) write (io_bptc_log,'(/a,i1/)') 'Iteration ', it
         depth_min = 999.
         depth_max = 0.
         do iev = 1,n_event ! Loop over events
            call depset (depthp(iev,it), usrc)
            rlatdg1 = latpgc(iev,it)
            rlondg1 = lonp(iev,it)
            call geogra (rlatdg1, glat)
            if (bptc) write (io_bptc_log,'(a,i3,1x,a,2f10.4,1x,f5.1)') 'Event ', iev, evtnam(iev), glat, rlondg1, depthp(iev,it)
            do i = 1,nst(iev) ! Loop over readings
            
               ! Theoretical travel time, azimuth, back-azimuth, delta, ray parameter,
               ! elipticity correction, station height correction. All relative to current
               ! event position.
               call mlocsteq (iev, i, rlatdg1, rlondg1, stladg(iev,i), stlndg(iev,i),&
                depthp(iev,it), ahgts(iev,i), phase(iev,i), ttcomp(iev,i), delt(iev,i),&
                azes(iev,i), azse(iev,i), psr, dtpdh, elcr(iev,i), hgtcr, latp(iev,it), .true.)
               psd(iev,i) = psr*rpd ! Convert ray parameter to sec/degree
               if (plogout) write (io_plog,'(a,2i4,1x,a6,1x,4f10.3,f6.1,f6.3,1x,a,f8.2,f10.2)') 'mlocsteq (ttcomp): ',&
                iev, i, stname(iev,i), rlatdg1, rlondg1, stladg(iev,i), stlndg(iev,i), depthp(iev,it), ahgts(iev,i),&
                phase(iev,i), delt(iev,i), ttcomp(iev,i)
                              
               ! If the named phase doesn't exist at this distance (ttcomp = 0.), and we're
               ! allowing phases to be renamed, reset the phase name
               if (ttcomp(iev,i) .lt. 0.1 .and. phidird(iev,i)) then
                  if (ptype(phase(iev,i))) then
                     if (phase(iev,i)(1:3) .ne. 'PPP') then
                        phase(iev,i) = 'UNKNOWNP' ! Unknown P-type
                        fcode(iev,i) = 'p'
                     end if
                  else if (stype(phase(iev,i))) then
                     if (phase(iev,i)(1:3) .ne. 'SSS') then
                        phase(iev,i) = 'UNKNOWNS' ! Unknown S-type
                        fcode(iev,i) = 'p'
                     end if
                  end if
               end if
               
               ! Spread of the travel time model for a given phase and epicentral distance
               if (read_ttsprd) then ! Use estimates from a previous run
                  call ttsig2 (phase(iev,i), delt(iev,i), ttsprd(iev,i), ttoff(iev,i))
               else ! default values of spread, all zero-offset
                  call ttsig  (phase(iev,i), delt(iev,i), ttsprd(iev,i), ttoff(iev,i))
               end if
               
               ! Partial derivatives
               a(iev,i,1) = psr*cos(azes(iev,i)*rpd)/radius ! geocentric
               a(iev,i,2) = -psr*sin(azes(iev,i)*rpd)/radius
               a(iev,i,3) = dtpdh
               if (rel_phase(iev,i) .or. rel_depth_phase(iev,i)) then ! Can't use relative phases to determine OT
                  a(iev,1,4) = 0.
               else
                  a(iev,i,4) = 1.
               end if
               
               !write (*,'(i3,1x,2a,3f10.3)') iev, stname(iev,i), phase (iev,i), a(iev,i,1), a(iev,i,2), a(iev,i,3)
               
               ! Station corrections 
               select case (tt_corr)
                  case (0) ! None
                     stacs(iev,i) = 0.0
                     sdstcs(iev,i) = 1.0
                  case (1) ! Station elevation corrections
                     stacs(iev,i) = hgtcr
                     sdstcs(iev,i) = 1.0
               end select
               
               ! Observed travel time in seconds, relative to current origin time for each event
               if (rel_phase(iev,i) .or. rel_depth_phase(iev,i)) then
                  arrtime = hms2s(ipah(iev,i), ipam(iev,i), pas(iev,i))
                  ttobs = arrtime
               else
                  arrtime = hms2s(ipah(iev,i), ipam(iev,i), pas(iev,i))
                  if (ipah(iev,i) .lt. hourp(iev,it)) arrtime = arrtime + 86400.
                  ttobs = arrtime - otsp(iev,it)
               end if
               tto(iev,i) = ttobs
               
               s(iev,i,it) = elcr(iev,i) + stacs(iev,i) ! Path anomalies: ellipticity and station corrections
               dt(iev,i,it) = ttobs - ttcomp(iev,i) ! Travel time residuals (uncorrected for path anomalies)
               
               ! List bounce point corrections for relative depth phases
               if (bptc .and. rel_depth_phase(iev,i)) write (io_bptc_log,'(2f10.2,1x,a,1x,i5)')&
                tto(iev,i), dt(iev,i,it), fcode(iev,i), mnf_line(iev,i)
                     
               ! Check for one-minute errors, first iteration only, only for unflagged readings
               if (it .eq. 0 .and. fcode(iev,i) .eq. ' ') then
                  res_test = abs(dt(iev,i,it))
                  if (res_test .ge. 55. .and. res_test .le. 65.) then
                     write (msg,'(a, i3,1x,a,1x,a,1x,i5, 1x,f6.2)') 'mlocset: possible one-minute error for ',&
                      iev, stname(iev,i), phase(iev,i), mnf_line(iev,i), dt(iev,i,it)
                     call warnings (trim(msg))
                  end if
               end if

               if (plogout) then
                  atsp = otsp(iev,it) + ttcomp(iev,i) + s(iev,i,it)
                  call timecr (atsp, atsp_hr, atsp_min, atsp_sec)
                  write (io_plog,'(a,2i3,f6.2)') 'mlocsteq (theoretical arrival time):', atsp_hr, atsp_min, atsp_sec
               end if
               
               ! Weights on the basis of absolute residual, i.e., relative to the window for that phase.
               ! For each phase the window is defined by a multiple (command WIND) of the spread calculated for that phase
               ! from a previous run (.ttsprd file). There is a roll-off at the edge if inverse weighting has been
               ! selected (command WEIG); otherwise it is a strict cut-off. The window is also offset by the
               ! mean residual from the .ttsprd file. 
               ! The window is expanded at local distance (delta < windloclim) for direct calibration,
               ! to help avoid losing important readings because of a poor starting location.
               weight(iev,i) = 1.
               p1 = wind1*ttsprd(iev,i)
               p2 = wind2*ttsprd(iev,i)
               if (direct_cal .and. delt(iev,i) .lt. windloclim) then
                  p1 = p1*2.
                  p2 = p2*2.
               end if
               dts = dt(iev,i,it) - s(iev,i,it) - ttoff(iev,i)
               if (debug) write (io_log,'(a,3i5,f10.2)') 'iev,i,it,dts: ', iev, i, it, dts
               if (data_weight) then
                  if (abs(dts) .le. p1) then
                     weight(iev,i) = 1.0
                  else if (abs(dts) .gt. p2) then
                     weight(iev,i) = 0.
                  else if (p2 .gt. p1) then
                     theta = pi*0.5*(abs(dts)-p1)/(p2-p1)
                     weight(iev,i) = 1.0 - sin(theta)
                     ! Make sure the weight is not too small. This can lead to a divide by zero
                     ! problem in mlocinv from round-off error.
                     if (weight(iev,i) .lt. minimum_weight) weight(iev,i) = minimum_weight
                  end if
               else
                  if (abs(dts) .le. p2) then
                     weight(iev,i) = 1.
                  else
                     weight(iev,i) = 0.
                  end if
               end if
               
            end do ! End of loop over readings
            
            ! Set flags (fltrh-, fltrc-) to determine if a station will be used to estimate an improved hypocentroid.
            ! To avoid convergence problems because of stations slipping in and out of the data set on alternate
            ! iterations, the flag is not reset after the iteration defined by "itlim_flag".
            if (it .le. itlim_flag) then
               call stflt (iev, it)
            else
               ndat(iev,it) = ndat(iev,it-1)
               ndatdl(iev,it) = ndatdl(iev,it-1)
               ndatpr(iev,it) = ndatpr(iev,it-1)
            end if
            write (io_log,'(i3,1x,3(i4,4x))') iev, ndat(iev,it), ndatdl(iev,it), ndatpr(iev,it)
            
            ! Minimum and maximum event depths, needed for depth testing with relative depth phases
            depth_min = min(depthp(iev,it),depth_min)
            depth_max = amax1(depthp(iev,it),depth_max)
            
         end do ! End of loop over events
         
         if (verbose_screen) then
            write (msg,'(a,f6.1,a,f6.1,a)') 'mlocset: event depth range is ', depth_min, ' to ', depth_max, ' km'
            call fyi (trim(msg))
         end if

         ! Adjust the depth range for testing with relative depth phase data
         h_rdp1 = int(depth_min - 40.)
         if (h_rdp1 .le. 0) h_rdp1 = 1
         h_rdp2 = int(depth_max + 40.)
         if (h_rdp2 .ge. n_hrdp_max) h_rdp2 = n_hrdp_max - 1
         if (verbose_screen) then
            write (msg,'(a,i3,a,i3,a)') 'mlocset: relative depth phase test range is ', h_rdp1, ' to ', h_rdp2, ' km'
            call fyi (trim(msg))
         end if
                  
         ! Determine the number of other events in the cluster with which each station-phase is associated,
         ! counting only those cases in which fltrh(iev,i) .eq. .false. It only matters that this number
         ! is at least 1, so we cut the search off when ic=1 for any station-phase. fltrc is reset whenever
         ! fltrh is (see above).
         if (it .le. itlim_flag) then
            if (debug) write (io_log,'(a)') ' iev    j stat  phase    h c  kcode   idf connected'
            do iev = 1,n_event ! Outer loop over events
               do j = 1,nst(iev) ! Outer loop over readings
                  phtest = phase(iev,j)
                  ic = 0
                  connected(iev,j) = .false.
                  if (.not.fltrc(iev,j)) then
                     k0 = kcode(iev,j)
                     id0 = idiff(iev,j)
                     do l = 1,n_event ! Inner loop over events
                        if (l .eq. iev) cycle ! Don't consider the same event, to avoid problems
                         ! with multiple readings for the same station-phase.
                        do k = 1,nst(l) ! Inner loop over readings
                           if (k0 .eq. kcode(l,k) .and. phtest .eq. phase(l,k) .and. id0 .eq. idiff(l,k)) then
                              if (.not.fltrc(l,k)) then
                                 ic = ic + 1
                                 if (ic .eq. 1) then ! That's all we need to know
                                    connected(iev,j) = .true.
                                    go to 100
                                 end if
                              end if
                           end if
                        end do ! End inner loop over phases
                     end do ! End inner loop over events
                  end if
  100             continue
                  if (debug) write (io_log,'(i4,i5,1x,a,1x,a,1x,l1,1x,l1,1x,i6,1x,i4,1x,l1)')&
                   iev, j, stname(iev,j), phase(iev,j), fltrh(iev,j), fltrc(iev,j),&
                   kcode(iev,j), idiff(iev,j), connected(iev,j)
               end do ! End outer loop over phases
            end do ! End outer loop over events
         end if
         
         if (nstep .eq. 0 .or. mt .eq. 0) then
            call write_summary (0) ! Summary data for the forward problem
            call write_phase_data (0) ! Phase data for the forward problem
            if (tt5) call tt_local (0, blank30, blank5, dum1, dum2) ! Local data
            if (tt5e) call single_event_tt5 (0) ! Single-event tt5 plots
            if (tt5s) call single_station_tt5 (0) ! Single-station tt5 plots
            if (bptc) close (io_bptc_log)
            return
         else
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            call mlocinv (it) ! Least squares inversion  (solve a * dx = dt)
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         end if
         
         ! Convert standard errors in epicenter from km to degrees
         sdxhath(1) = sdxhath(1)*dgkmla
         sdxhath(2) = sdxhath(2)*dgkmlogc(lathgc(it))
         do iev = 1,n_event
            sdxhatc(iev,1) = sdxhatc(iev,1)*dgkmla
            sdxhatc(iev,2) = sdxhatc(iev,2)*dgkmlogc(latpgc(iev,it))
         end do
         
         ! Update hypocentroid
         lathgc(it+1) = lathgc(it) + delx0(1,it)*dgkmla
         lonh(it+1) = lonh(it) + delx0(2,it)*dgkmlogc(lathgc(it))
         depthh(it+1) = depthh(it) + delx0(3,it)
         otsh(it+1) = otsh(it) + delx0(4,it)
         call timecr (otsh(it+1), hourh(it+1), minh(it+1), sech(it+1))
         
         ! Update cluster event hypocenters
         ! If a cluster vector parameter (e.g., depth) is fixed but the corresponding parameter for the hypocentroid is free,
         ! the cluster vector is updated with the change in the hypocentroid parameter.
         do iev = 1,n_event ! Loop over events
         
            if (mindx(iev,1) .ne. 0 .or. latfh) then ! Latitude
               latpgc(iev,it+1) = latpgc(iev,it) + (delx0(1,it) + dxp(iev,1,it))*dgkmla
               call geogra (latpgc(iev,it+1), latp(iev,it+1)) ! Update latitude in geographic coordinates.
            else
               latpgc(iev,it+1) = latpgc(iev,it)
            end if
            
            if (mindx(iev,2) .ne. 0 .or. lonfh) then ! Longitude
               lonp(iev,it+1) = lonp(iev,it) + (delx0(2,it) + dxp(iev,2,it))*dgkmlogc(latpgc(iev,it))
            else
               lonp(iev,it+1) = lonp(iev,it)
            end if
            
            if (mindx(iev,3) .ne. 0 .or. depthfh) then ! Depth
               depthp(iev,it+1) = depthp(iev,it) + delx0(3,it) + dxp(iev,3,it)
               if (depthp(iev,it+1) .le. 0.) then
                  write (msg,'(a,i3,a)') 'mlocset: depth for event ', iev,' went negative - held at previous value'
                  call warnings (trim(msg))
                  depthp(iev,it+1) = depthp(iev,it)
               end if
               ddep = depthp(iev,it+1) - depthp(iev,it)
               if (abs(ddep) .gt. 50.) then
                  depthp(iev,it+1) = depthp(iev,it+1) + ddep/2.
                  write (msg,'(a,i3,a)') 'mlocset: depth change for event ', iev,' was damped'
                  call warnings (trim(msg))
               end if
            else
               depthp(iev,it+1) = depthp(iev,it)
            end if
            
            if (mindx(iev,4) .ne. 0 .or. timefh) then ! OT
               otsp(iev,it+1) = otsp(iev,it) + delx0(4,it) + dxp(iev,4,it)
            else
               otsp(iev,it+1) = otsp(iev,it)
            end if
             
         end do ! End of loop over events
         
         ! Convert OT in seconds to hours-minutes-seconds
         do iev = 1,n_event
            call timecr (otsp(iev,it+1), hourp(iev,it+1), minp(iev,it+1), secp(iev,it+1))
         end do
         
         ! Print changes in cluster vectors
         write (io_log,'(/a)') 'Cluster vector changes (geocentric latitude)'
         write (io_log,'(a)') ' iev   dtim(s)   dlat(¡)   dlon(¡)  ddep(km)'
         do iev = 1,n_event
            dtim = dxp(iev,4,it)
            dlat = dxp(iev,1,it)*dgkmla
            dlon = dxp(iev,2,it)*dgkmlogc(latpgc(iev,it))
            ddep = dxp(iev,3,it)
            write (io_log,'(i4,4f10.4)') iev, dtim, dlat, dlon, ddep
         end do
         
         ! Recalculate the hypocentroid. This is necessary to keep the hypocentroid synchronized with the cluster vectors
         ! in the case where some parameter (e.g., depth) is fixed for some events and not others.
         x1sum = 0.
         x2sum = 0.
         x3sum = 0.
         x4sum = 0.
         do iev = 1,n_event
            x1sum = x1sum + latpgc(iev,it+1)
            x2sum = x2sum + lonp(iev,it+1)
            x3sum = x3sum + depthp(iev,it+1)
            x4sum = x4sum + otsp(iev,it+1)
         end do
         lathgc(it+1) = x1sum/real(n_event)
         call geogra (lathgc(it+1),lath(it+1)) ! Update in geographic coordinates
         lonh(it+1) = x2sum/real(n_event)
         depthh(it+1) = x3sum/real(n_event)
         otsh(it+1) = x4sum/real(n_event)
         call timecr (otsh(it+1), hourh(it+1), minh(it+1), sech(it+1))
         
         ! Hypocentroid for calculation of empirical path anomalies
         lath_epa = lathgc(it+1)
         lonh_epa = lonh(it+1)
         depthh_epa = depthh(it+1)
         
         ! Print changes in hypocentroid
         dtim = otsh(it+1) - otsh(it)
         dlat = lathgc(it+1) - lathgc(it)
         dlon = lonh(it+1) - lonh(it)
         ddep = depthh(it+1) - depthh(it)
         write (io_log,'(/a)') 'Hypocentroid changes (geocentric latitude)'
         write (io_log,'(t5,a)') '   dtim(s)   dlat(¡)   dlon(¡)  ddep(km)'
         write (io_log,'(t5,4f10.4)') dtim, dlat, dlon, ddep
         
         ! Check for convergence
         convrg = .true.
         select case (convergence_test_index)
         
            case (0) ! Convergence tests are done separately on the cluster vectors and hypocentroid
               write (io_log,'(/a)') 'Convergence tests are done separately on the cluster vectors and hypocentroid'
               write (io_log,'(a,f10.3,a,f10.3,a)') 'Epicenter criterion:   ', cl_epi_h, 'on hypocentroid; ',&
                cl_epi_c, ' on cluster vectors (km)'
               write (io_log,'(a,f10.3,a,f10.3,a)') 'Depth criterion:       ', cl_dep_h, 'on hypocentroid; ',&
                cl_dep_c, ' on cluster vectors (km)'
               write (io_log,'(a,f10.3,a,f10.3,a)') 'Origin time criterion: ', cl_ot_h, 'on hypocentroid; ',&
                cl_ot_c, ' on cluster vectors (s)'
               if (abs(lathgc(it+1) - lathgc(it)) .gt. cl_epi_h) then
                  convrg = .false.
                  write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lat', abs(lathgc(it+1) - lathgc(it)), ' > ', cl_epi_h
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Hypocentroid latitude change: ', lathgc(it+1) - lathgc(it), ' deg'
               end if
               if (abs(lonh(it+1) - lonh(it)) .gt. cl_epi_h) then
                  convrg = .false.
                  write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lon', abs(lonh(it+1) - lonh(it)), ' > ', cl_epi_h
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Hypocentroid longitude change: ', lonh(it+1) - lonh(it), ' deg'
               end if
               if (abs(depthh(it+1) - depthh(it)) .gt. cl_dep_h) then
                  convrg = .false.
                  write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' dep', abs(depthh(it+1) - depthh(it)), ' > ', cl_dep_h
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Hypocentroid depth change: ', depthh(it+1) - depthh(it), ' km'
               end if
               if (abs(otsh(it+1) - otsh(it)) .gt. cl_ot_h) then
                  convrg = .false.
                  write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' tim', abs(otsh(it+1) - otsh(it)), ' > ', cl_ot_h
                  if (it .eq. nstep) write (*,'(a,f8.4,a)') ' Origin time change: ', otsh(it+1) - otsh(it), ' sec'
               end if
               do iev = 1,n_event
                  if (abs(dxp(iev,1,it)) .gt. cl_epi_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lat', abs(dxp(iev,1,it)), ' > ', cl_epi_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' latitude change: ', dxp(iev,1,it), ' km'
                  end if
                  if (abs(dxp(iev,2,it)) .gt. cl_epi_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lon', abs(dxp(iev,2,it)), ' > ', cl_epi_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' longitude change: ', dxp(iev,2,it), ' km'
                  end if
                  if (abs(dxp(iev,3,it)) .gt. cl_dep_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' dep', abs(dxp(iev,3,it)), ' > ', cl_dep_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' depth change: ', dxp(iev,3,it), ' km'
                  end if
                  if (abs(dxp(iev,4,it)) .gt. cl_ot_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' tim', abs(dxp(iev,4,it)), ' > ', cl_ot_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' origin time change: ', dxp(iev,4,it), ' sec'
                  end if
               end do
               
            case (1) ! Convergence tests are done on final event hypocenters
               write (io_log,'(/a)') 'Convergence tests are done on final event hypocenters'
               write (io_log,'(a,f10.3,a)') 'Epicenter criterion:   ', cl_epi_c, ' (km)'
               write (io_log,'(a,f10.3,a)') 'Depth criterion:       ', cl_dep_c, ' (km)'
               write (io_log,'(a,f10.3,a)') 'Origin time criterion: ', cl_ot_c, ' (s)'
               do iev = 1, n_event
                  dlat = (latpgc(iev,it+1) - latpgc(iev,it))/dgkmla
                  if (abs(dlat) .gt. cl_epi_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lat', abs(dlat), ' > ', cl_epi_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' latitude change: ', dlat, ' km'
                  else
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lat', abs(dlat), ' < ', cl_epi_c
                  end if
                  dlon = (lonp(iev,it+1) - lonp(iev,it))/dgkmlogc(latpgc(iev,it))
                  if (abs(dlon) .gt. cl_epi_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lon', abs(dlon), ' > ', cl_epi_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' longitude change: ', dlon, ' km'
                  else
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' lon', abs(dlon), ' < ', cl_epi_c
                  end if
                  ddep = depthp(iev,it+1) - depthp(iev,it)
                  if (abs(ddep) .gt. cl_dep_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' dep', abs(ddep), ' > ', cl_dep_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' depth change: ', ddep, ' km'
                  else
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' dep', abs(ddep), ' < ', cl_dep_c
                  end if
                  dtim = otsp(iev,it+1) - otsp(iev,it)
                  if (abs(dtim) .gt. cl_ot_c) then
                     convrg = .false.
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' tim', abs(dtim), ' > ', cl_ot_c
                     if (it .eq. nstep) write (*,'(a,i3,a,f8.4,a)') ' Event ', iev,' origin time change: ', dtim, ' sec'
                  else
                     write (io_log,'(i4,a,f10.3,a,f10.3)') iev, ' tim', abs(dtim), ' < ', cl_ot_c
                  end if
               end do
            
         end select
         
         ! Even if convergence criteria are met we require at least one iteration in order
         ! to generate values for eci (cluster vectors).
         if (convrg .and. it .eq. 0) convrg = .false.

         if (convrg .or. it .eq. nstep) then
            if (convrg) then
               write (*,'(/a,i1,a/)') 'Converged after ', it, ' iterations'
               write (io_log,'(/a)') '*****************************************************'
               write (io_log,'(a,i1,a/)') 'Converged after ', it, ' iterations'
               write (io_log,'(a)') '*****************************************************'
            else
               write (*,'(/a,i1,a/)') 'Not converged after ', it, ' iterations'
               write (io_log,'(/a)') '*****************************************************'
               write (io_log,'(a,i1,a/)') 'Not converged after ', it, ' iterations'
               write (io_log,'(a)') '*****************************************************'
            end if
            
            ! Update geographic latitude for hypocentroid and cluster events
            call geogra (lathgc(it+1), lath(it+1))
            do iev = 1,n_event
               call geogra (latpgc(iev,it+1), latp(iev,it+1))
            end do
            
            if (.not. bias_corr) then
               lathgc(it+1) = lathgc(it+1) + bcorr(1)*dgkmla
               call geogra (lathgc(it+1), lath(it+1)) ! Convert to geographic coordinates
               lonh(it+1) = lonh(it+1) + bcorr(2)*dgkmlogc(lathgc(it+1))
               depthh(it+1) = depthh(it+1) + bcorr(3)
               otsh(it+1) = otsh(it+1) + bcorr(4)
               call timecr (otsh(it+1), hourh(it+1), minh(it+1), sech(it+1)) 
               do iev = 1,n_event
                  latpgc(iev,it+1) = latpgc(iev,it+1) + bcorr(1)*dgkmla
                  call geogra (latpgc(iev,it+1), latp(iev,it+1)) ! Convert to geographic coordinates
                  lonp(iev,it+1) = lonp(iev,it+1) + bcorr(2)*dgkmlogc(latpgc(iev,it+1))
                  depthp(iev,it+1) = depthp(iev,it+1) + bcorr(3)
                  otsp(iev,it+1) = otsp(iev,it+1) + bcorr(4)
               end do
               do iev = 1,n_event
                  call timecr (otsp(iev,it+1), hourp(iev,it+1),minp(iev,it+1), secp(iev,it+1))
               end do
            end if
            
            ! Calibration events
            if (indirect_cal) then
            
               ! Calculate shift vector to match calibration locations
               call cal_shift (1, it)
               do iev = 1,n_event ! Add shift vector to all cluster events
                  latp_cal(iev) = latp(iev,it+1) + del_cal_lat
                  lonp_cal(iev) = lonp(iev,it+1) + del_cal_lon
                  depthp_cal(iev) = depthp(iev,it+1) + del_cal_dep
                  ! Don't let focal depth be shifted to a negative value
                  if (depthp_cal(iev) .lt. 0.) then
                     write (msg,'(a,i3,a,f6.2,a)') 'mlocset: calibration shift wants to make depth negative for event',&
                      iev, ' (', depthp_cal(iev), '); reset to 0 km'
                     call warnings (trim(msg))
                     depthp_cal(iev) = 0.0
                  end if
                  otsp_cal(iev) = otsp(iev,it+1) + del_cal_tim
               end do
               ! Update calibrated hypocentroid for calculation of empirical path anomalies.
               lath_epa = lath_epa + del_cal_lat ! Geocentric
               lonh_epa = lonh_epa + del_cal_lon
               depthh_epa = depthh_epa + del_cal_dep
               
               ! Residuals relative to indirect calibrated (shifted) locations
               if (bptc) write (io_bptc_log,'(/a/)') 'Indirect calibration shifted locations'
               do iev = 1,n_event ! Loop over events
                  call geocen (latp_cal(iev), lonp_cal(iev), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
                  call depset (depthp_cal(iev), usrc)
                  if (bptc) write (io_bptc_log,'(a,i3,1x,a,2f10.4,1x,f5.1)') 'Event ',&
                   iev, evtnam(iev), latp_cal(iev), lonp_cal(iev), depthp_cal(iev)
                  do i = 1,nst(iev) ! Loop over readings
                     call mlocsteq (iev, i, xlat, xlon, stladg(iev,i), stlndg(iev,i),&
                      depthp_cal(iev), ahgts(iev,i), phase(iev,i), ttcomp(iev,i), delt_cal(iev,i),&
                      azes_cal(iev,i), azse_cal(iev,i), psr, dtpdh, elcr_cal(iev,i), hgtcr,&
                      latp_cal(iev), .true.)
                     psd(iev,i) = psr*rpd ! Convert ray parameter to sec/degree
                     
                     ! TT corrections
                     select case (tt_corr)
                     case (0) ! None
                        stacs_cal(iev,i) = 0.0
                        sdstcs_cal(iev,i) = 1.0
                     case (1) ! Station elevation corrections
                        stacs_cal(iev,i) = hgtcr
                        sdstcs_cal(iev,i) = 1.0
                     case (2) ! Patch corrections + station elevation corrections
                        if (delt_cal(iev,i) .ge. 25.) then ! Patch corrections valid for teleseismic only
                           stacs_cal(iev,i) = sca0s(iev,i) + hgtcr
                           sdstcs_cal(iev,i) = sd1s(iev,i)
                        else
                           stacs_cal(iev,i) = hgtcr
                           sdstcs(iev,i) = 1.0
                        end if
                     end select
                     
                     arrtime = hms2s(ipah(iev,i), ipam(iev,i), pas(iev,i))
                     if (ipah(iev,i) .lt. hourp(iev,it)) arrtime = arrtime + 86400.
                     
                     if (rel_phase(iev,i) .or. rel_depth_phase(iev,i)) then
                        ttobs = real(ipam(iev,i)*60) + pas(iev,i)
                     else
                        ttobs = arrtime - otsp_cal(iev)
                     end if
                     tto(iev,i) = ttobs
                     s_cal(iev,i) = elcr_cal(iev,i) + stacs_cal(iev,i) ! Path anomalies: ellipticity and station corrections
                     dt_cal(iev,i) = ttobs - ttcomp(iev,i) ! Travel time residuals (uncorrected for path anomalies)
                     
                     ! List bounce point corrections for relative depth phases
                     if (bptc .and. rel_depth_phase(iev,i)) write (io_bptc_log,'(2f10.2,1x,a,1x,i5)')&
                      tto(iev,i), dt_cal(iev,i), fcode(iev,i), mnf_line(iev,i)
                     
                  end do ! End of loop over readings
               end do ! End of loop over events
                              
            end if
            
            ! Depth analysis using relative depth phases
            call fyi ('mlocset: depth determination based on relative depth phases')
            call rdp_depth_test (it)
            close (io_dpnc)
            
            ! Spread of focal depths
            call depth_spread (it)
               
            exit ! Break out of iteration loop
            
         else ! Not converged
         
            it = it + 1 ! Increment iteration counter
            
            ! Re-identify phases based on new location and OT
            if (phid .and. it .le. 1) then
               nphreid = 0
               do iev = 1,n_event
                  do ird = 1,nst(iev)
                     if (phidird(iev,ird)) then
                        call phreid3 (it, iev, ird, nr)
                        nphreid = nphreid + nr
                     end if
                  end do
               end do
               write (msg,'(a,i6,a,i2)') 'mlocset: ', nphreid, ' phases re-identified for iteration ', it
               if (verbose_screen) call fyi (trim(msg))
               write (io_log,'(/a)') trim(msg)
            end if
                        
         end if
         
      end do ! End of iteration loop
      
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      if (bptc) close (io_bptc_log)

      ! Printed output
      call fyi ('mlocset: Summary file')
      call write_summary (it) ! Summary data (.summary), covariance (.cv)
      if (direct_cal) then
         call fyi ('mlocset: Direct Calibration file')
         call dircal (it) ! Direct calibration (must be called before mlocout_hdf!)
      end if
      call fyi ('mlocset: HDF file(s)')
      call write_hdf (it) ! .hdf and .hdf_cal if there is calibration data
      call fyi ('mlocset: Phase Data file')
      call write_phase_data (it) ! Phase data (.phase_data), large residuals (.xdat and .lres)
      call fyi ('mlocset: RDERR file')
      call write_rderr () ! Empirical reading errors (.rderr)
      if (diffdat) call write_rderr_diff () ! Empirical reading errors for differential time
       ! data (.rderr_diff)
      call fyi ('mlocset: TTSPRD file')
      call write_ttsprd (it) ! TT spread for different phases (.ttsprd)
            
      ! GMT scripts
      call fyi ('mlocset: Base plots...')
      if (gmtp) then
         confidence_ellipse_plot = (n_event .gt. 1)
         call base_map (it, 0, .true., confidence_ellipse_plot, .true., '_base') ! Base plot
         if (eplt) call base_map (it, 0, .false., .true., .false., '_ell') ! Confidence ellipses,
          ! no event numbers or vectors
         if (splt) call base_map (it, 0, .false., .false., .true., '_seis') ! Seismicity plot,
          ! no event numbers or confidence ellipses
         if (direct_cal) call dcal_map (it) ! direct calibration raypath plot
         do ip = 1,n_selected_max
            do iev = 1,n_event
               if (plot(iev,ip)) then
                  call base_map (it, ip, .true., .true., .true., '_sel') ! Selected events, base plot
                  exit
               end if
            end do
         end do
         if (n_xsec .gt. 0) then
            do i = 1,n_xsec
               call xsec_gmt (it, i)
            end do
         end if
      else
         call warnings ('mlocset: Plotting has been disabled')
      end if
      
      ! If it doesn't already exist, make a folder for the symbols used in the .kml file
      inquire (file=trim(datadir)//dirsym//'_kml/a_blue.png',exist=ex)
      if (.not. ex) then
         call fyi ('mlocset: Creating _kml directory')
         command_line = 'cp -R tables/kml '//trim(datadir)//'/_kml'
         call execute_command_line (trim(command_line))
      end if      
      call fyi ('mlocset: KML file')
      call write_kml (it) ! KML file for display in Google Earth
      
      if (gccelout) then
         call fyi ('mlocset: GCCEL files')
         call write_gccel (it) ! GCCEL format (~_gccel.dat)
         call gccel_hypocentroid (it) ! Entry for doc.kml
         call gccel_epicenters () ! Entry for epicenters.kml
      end if

      if (blocout) then
         call fyi ('mlocset: BayesLoc file')      
         call write_bloc (it) ! BAYESLOC format (.bloc)
      end if
                                          
      if (nitomo .gt. 0) then ! Tomography output files (.TOMO_PHASE.ITOMO.tomo)
         ! Make a folder for the files
         file_folder = trim(datadir)//dirsym//trim(basename)//'_tomo'
         command_line = 'mkdir '//trim(file_folder)
         call execute_command_line (trim(command_line))
         call fyi ('mlocset: Tomography file(s)')      
         do i = 1,nitomo
            call tomo (i, it, tomo_phase(i), itomo(i), file_folder)
         end do
      end if
      
      if (ttou) then
         if (.not.indirect_cal .and. .not.direct_cal)&
          call warnings ('mlocset: Empirical TTs are from an uncalibrated cluster')
         call fyi ('mlocset: Empirical TT file(s)')      
         call write_tt (it) ! Output empirical TT data for specific phases
      end if
      
      ! .datf file (event data files with flags and phase IDs as actually used)
      if (datfout) then
         call fyi ('mlocset: DATF file')
         outfil = trim(outfile)//'_datf.mnf'
         call open_file (io_datf,  outfil, 'new')
         write (io_datf,'(a,3x,a,t125,a)') 'B', 'final disposition of '//trim(outfile), ' '
         write (io_datf,'(a,t125,a)') 'F   MNF v1.3.4', ' '
         do iev = 1,n_event
            call write_mnf_134 (io_datf, iev, it)
         end do
         write (io_datf,'(a,t125,a)') 'EOF', ' '
         close (io_datf)
      end if
      
      ! Travel time and other plots
      call fyi ('mlocset: Travel-time plots')
      if (gmtp) then
         if (tt1) call tt_summary (it) ! Summary travel time plot
         if (tt2) call tt_teleseismic_p (it) ! Teleseismic P
         if (tt3) call tt_pkp_caustic2 (it) ! PKP caustic
         if (tt4) call tt_near_source (it) ! Near source readings
         if (tt5) call tt_local (it, blank30, blank5, dum1, dum2) ! Local distance
         if (tt5e) call single_event_tt5 (it) ! Single-event tt5 plots
         if (tt5s) call single_station_tt5 (it) ! Single-station tt5 plots
         If (tt6) call tt_local_regional (it) ! Local-regional distance
         if (tt7) call tt_local_regional_s (it) ! Local-regional shear phases
         if (tt8) call tt_rdp_summary (it) ! Relative depth phase summary plot
         if (rdpp) call rdp_single (it) ! Relative depth phase plots for single events
         if (tt9) call tt_s_minus_p () ! S-P times
         if (epa_plot) call epa_plot_driver (it) ! Empirical path anomaly plots
         if (fdhp) call focal_depth_histogram (it) ! Histogram of focal depths
         if (phrp) call phrp_plot_driver (it) ! Residuals for a specified phase
         if (evrp) call event_res_plot (it)

         ! For GCCEL output combine all plots into a single PDF, then delete the postscript files
         if (gccelout) then
            command_line = 'gmt psconvert -TF -Fmloc_plots -A1 '//trim(gcat_folder)//'/*.ps'
            call execute_command_line (trim(command_line))
            command_line = 'mv mloc_plots.pdf '//trim(gcat_folder)//dirsym//trim(basename)//'_plots.pdf'
            call execute_command_line (trim(command_line))
            command_line = 'rm '//trim(gcat_folder)//'/*.ps'
            call execute_command_line (trim(command_line))
         end if
      else
         call warnings ('mlocset: Plotting has been disabled')
      end if    
      
      return
      
      end subroutine mlocset
      

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
!      real :: ellip
!      real :: hms2s
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
      character(len=132) :: msg
     
      openazlim = 300.
      
      if (verbose_screen) write (*,'(/a/)') 'Starting model and free parameters:'
      
      do iev = 1,n_event
      
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
   !             '¡ for event ', iev, evtnam(iev), lonp(iev,0), lonp(iev-1,0)
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
      do iev = 1,n_event
         call geocen (latp(iev,0), lonp(iev,0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
         latpgc(iev,0) = xlat
         lonp(iev,0) = xlon
      end do
   
      ! Initial estimate of hypocentroid
      x1sum = 0.
      x2sum = 0.
      x3sum = 0.
      x4sum = 0.
      do iev = 1,n_event
         x1sum = x1sum + latp(iev,0)
         x2sum = x2sum + lonp(iev,0)
         x3sum = x3sum + depthp(iev,0)
         x4sum = x4sum + otsp(iev,0)
      end do
      lath(0) = x1sum/real(n_event)
      lonh(0) = x2sum/real(n_event)
      call geocen (lath(0), lonh(0), xlat, sxlat, cxlat, xlon, sxlon, cxlon)
      lathgc(0) = xlat
      depthh(0) = x3sum/real(n_event)
      otsh(0) = x4sum/real(n_event)
      call timecr (otsh(0), hourh(0), minh(0), sech(0)) 
   
      ! Number of free parameters for each event and total.
      mt = 0
      do iev = 1,n_event
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
      tcor = ellip ('P       ', 30., 30., 30., 30., .true.)
      if (verbose_screen) then
         write (msg,'(a,f5.2)') 'model_init: Initialize ellipticity, tcor = ', tcor
         call fyi (trim(msg))
      end if
      
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
         if (nhdf .gt. n_hdf_max) then
            call warnings ('hdf_read: n_hdf_max exceeded')
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
      if (verbose_screen) then
         write (msg,'(a,i3,2a)') 'hdf_read: ', nrhdf, ' starting locations read from:',&
          trim(rhdf_filnam)
         call fyi (trim(msg))
      end if
      
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
      if (verbose_screen) then
         write (msg,'(a,i3,2a)') 'read_bdp: ', n_bdp, ' entries read from ', trim(bdp_filnam)
         call fyi (trim(msg))
      end if
      
      close (io_bdp)
      
      return
      
   end subroutine read_bdp
   
   
!*****************************************************************************************
   subroutine depth_spread (it)
   
   ! Robust estimate of spread of constrained focal depths
   
      real, dimension(n_event) :: cfd
      integer :: iev
      integer, intent(in) :: it ! iteration index
      character(len=132) :: msg
      
      do iev = 1,n_event
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
end module mloc_set

