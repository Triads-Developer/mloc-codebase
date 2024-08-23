!> Procedures related to tomography output files

module mlocout_tomo

   use declare_limits
   use declare_calibration
   use declare_cluster_vectors
   use declare_configuration
   use declare_constants
   use declare_environment
   use declare_events
   use declare_lun
   use declare_output
   use declare_phase_data
   use declare_plotting
   use declare_stations
   use declare_tt
   use mloclib_geog
   use mloclib_statistics
   
   implicit none
   save

   integer :: hourpr
   integer :: minpr
   real :: secpr
   real :: latout
   real :: lonout
   real :: depout
   
contains

!*****************************************************************************************
   subroutine tomo (ntomo, it, tomo_phase_in, itomo_in, folder)
         
   ! Output file for use in tomography. "tomo_phase" is a phase name for which data are extracted and "itomo" is
   ! a flag for the type of readings included:
   !  itomo = 1 Extract all readings of the specified phase
   !  itomo = 2 Extract only readings which were used for the cluster vectors
   !  itomo = 3 Extract empirical path anomalies
   
      integer :: i
      integer, dimension(npmax) :: idum
      integer :: iev
      integer :: ii
      integer, dimension(npmax) :: indx
      integer, dimension(20) :: indx2
      integer :: it
      integer :: it1
      integer :: itomo_in
      integer :: j
      integer :: jj
      integer :: kk
      integer :: ntomo
      real, dimension(npmax) :: deltiev
      real, dimension(20) :: sortime
      character(len=4) :: cell_lat_inc_chr
      character(len=4) :: cell_lon_inc_chr
      character(len=132) :: folder
      character(len=1) :: itomo_char
      character(len=132) :: msg
      character(len=100) :: outfil
      character(len=5) :: stnam
      character(len=8) :: tomo_phase_in
      logical :: loop
            
      write (itomo_char, '(i1)') itomo_in
      outfil = trim(folder)//dirsym//trim(basename)//'_'//trim(tomo_phase_in)//'_'//itomo_char//'.tomo'
      if (tomo_sub(ntomo)) then
         write (cell_lat_inc_chr, '(f4.2)') cell_lat_inc(ntomo)
         write (cell_lon_inc_chr, '(f4.2)') cell_lon_inc(ntomo)
         do i = 1,4
            if (cell_lat_inc_chr(i:i) .eq. ' ') cell_lat_inc_chr(i:i) = '_'
            if (cell_lon_inc_chr(i:i) .eq. ' ') cell_lon_inc_chr(i:i) = '_'
         end do
         outfil = trim(outfil)//'g_'//cell_lat_inc_chr//'_'//cell_lon_inc_chr
      end if
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_tomo: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      call open_file (io_out, outfil, 'new')
      
      it1 = it + 1
      
      if (itomo_in .eq. 1 .or. itomo_in .eq. 2) then
         do iev = 1,nev ! Loop over events
         
            ! Hypocenter info
            ! Uncalibrated cluster or direct calibration
            hourpr = hourp(iev,it1)
            minpr = minp(iev,it1)
            secpr = secp(iev,it1)
            latout = latp(iev,it1)
            lonout = lonp(iev,it1)
            depout = depthp(iev,it1)
            if (indirect_cal) then ! Indirect calibration trumps direct calibration
               call timecr (otsp_cal(iev), hourpr, minpr, secpr)
               latout = latp_cal(iev)
               lonout = lonp_cal(iev)
               depout = depthp_cal(iev)
            end if
            call set_longitude_range (lonout, longitude_range)
            
            ! Index table on delta for sorting phase reading lines
            do i = 1,nst(iev)
               deltiev(i) = delt(iev,i)
            end do
            call indexx (nst(iev), deltiev, indx)
            
            ! Check that the sorting has not disrupted the time order of multiple readings from the same station
            do ii = 1,nst(iev) - 1
               i = indx(ii)
               stnam = stname(iev,i)
               loop = .true.
               j = 1
               do while (loop)
                  if (stname(iev,indx(ii+j)) .eq. stnam) then
                     if ((j+ii) .lt. nst(iev)) then
                        loop = .true.
                        j = j + 1
                     else if ((j + ii) .eq. nst(iev)) then ! End of the list
                        loop = .false.
                     else
                        loop = .false.
                     end if
                  else ! New station name encountered
                     loop = .false.
                  end if
               end do
               if (j .gt. 1) then
                  do jj = 1,j
                     kk = indx(ii + jj - 1)
                     sortime(jj) = ipah(iev,kk)*3600. + ipam(iev,kk)*60. + pas(iev,kk)
                  end do
                  call indexx (j, sortime, indx2)
                  do jj = 1,nst(iev)
                     idum(jj) = indx(jj)
                  end do
                  do jj = 1,j
                     indx(ii + jj - 1) = idum(ii + indx2(jj) - 1)
                  end do
               end if
            end do
            
            ! Write phase readings
            
            do ii = 1,nst(iev)
               i = indx(ii)
               if (phase(iev,i) .ne. tomo_phase_in) cycle
               if (itomo_in .eq. 1) then
   !            if (fcode (iev,i) .eq. ' ') call tomo2 (it, iev, i)
                  if (.not.fltrhres(iev,i) .and. .not.fltrhflag(iev,i) .and. .not.fltrcres(iev,i) .and. &
                   .not.fltrcflag(iev,i)) call tomo2 (it, iev, i)
               else if (itomo_in .eq. 2) then
                  if (connected(iev,i)) call tomo2 (it, iev, i)
               end if
            end do
            
         end do
      else if (itomo_in .eq. 3) then
         if (tomo_sub(ntomo)) then
            call tomo_grid (it, ntomo)
            call tomo3grid (it, tomo_phase_in, ntomo)
         else
            call tomo3 (it, tomo_phase_in)
         end if
      end if
      
      close (io_out)
      
      return
      
   end subroutine tomo
         
         
!*****************************************************************************************
   subroutine tomo2 (it, iev, ii)
   
   ! Write a tomo format phase reading line.
         
      integer :: i
      integer :: iev
      integer :: igt
      integer :: ii
      integer :: it
      integer :: nobsere
      real :: deltout
      real :: dts
      real :: stladg_geog
      character(len=21) :: qtest
      
      deltout = delt(iev,ii)
      if (indirect_cal)  deltout = delt_cal(iev,ii)
      
      ! Travel time residual
      dts = dt(iev,ii,it) - s(iev,ii, it) ! Direct calibration or uncalibrated
      if (indirect_cal) dts = dt_cal(iev,ii) - s_cal(iev,ii) ! Indirect calibration trumps others
      
      call geogra (stladg(iev,ii), stladg_geog) ! Convert from geocentric to geographic latitude
      
      ! Empirical reading error and number of observations
      qtest=stname(iev,ii)//deployment(iev,ii)//phase(iev,ii)
      nobsere = 1
      do i = 1,nqc
         if (qtest .eq. qname1(i)) then
            nobsere = indexq(i)
            exit
         end if
      end do
      
      ! Calibration level (Indirect trumps direct)
      if (indirect_cal) then
         igt = nint(xl2cg(iev))
      else if (direct_cal) then
         igt = nint(xl2dc(iev))
      else
         igt = 99 ! No calibration
      end if
      
      write (io_out,'(f7.3,f7.2,f7.2,1x,i3,1x,i4,2i2.2,1x,2i2.2,f5.2,2x,i4,2i2.2,1x,2i2.2,f5.2,f8.2,'//&
      '2f8.3,f6.1,1x,2f8.3,f6.1,1x,a,1x,i3,1x,i3,1x,a)')&
       deltout,& ! Epicentral distance (deg)
       dts,& ! TT residual (s)
       sdread(iev,ii),& ! Reading error (s)
       nobsere,& ! Number of observations of this station-phase
       iyre(iev),& ! Reading year
       mone(iev),& ! Reading month
       idye(iev),& ! Reading day
       ipah(iev,ii),& ! Reading hour
       ipam(iev,ii),& ! Reading minute
       pas(iev,ii),& ! Reading seconds
       iyre(iev),& ! Origin year
       mone(iev),& ! Origin month
       idye(iev),& ! Origin day
       hourpr,& ! Origin hour
       minpr,& ! Origin minute
       secpr,& ! Origin seconds
       tto(iev,ii),& ! Observed travel time (s)
       latout,& ! Event latitude
       lonout,& ! Event longitude
       depout,& ! Event depth (km)
       stladg_geog,& ! Station latitude
       stlndg(iev,ii),& ! Station longitude
       ahgts(iev,ii),& ! Station elevation
       qtest,& ! Station-phase name
       iev,& ! Event number
       igt,& ! Calibration level (99 for no uncalibrated)
       trim(basename) ! Title
       
      return
      
   end subroutine tomo2
   
   
!*****************************************************************************************
   subroutine tomo3 (it, tomo_phase_in)
   
   ! Tomography output file of empirical path anomalies. Based on subroutine mlocout_rderr.
   ! The full travel time associated with an empirical path anomaly is calculated by adding the
   ! empirical path anomaly to the theoretical travel time between the hypocentroid and the
   ! appropriate station. In this case, the hypocentroid is not of the entire cluster, but the
   ! hypocentroid of those events which contribute to this station-phase. The entries for origin
   ! time and reading time are ignored.
   
      integer, parameter :: max = 60
      
      integer :: i
      integer :: iev
      integer :: igt
      integer :: it
      integer :: j
      real :: delt1
      real :: depth_epa
      real :: dum1
      real :: dum2
      real :: dum3
      real :: dum4
      real :: dum5
      real :: dum6
      real :: lat_epa
      real :: lat_epa_gc
      real :: lon_epa
      real :: qlatdg
      real :: qlondg
      real :: t1
      real :: t2
      real :: t3
      real :: t4
      real :: ttepa
      real, dimension(2) :: usrc
      real :: x1sum
      real :: x2sum
      real :: x3sum
      character(len=8) :: tomo_phase_in
      
      write (io_log,'(/a)') 'tomo3:'
      
      do i = 1,nqc ! Loop over station-phases with data
         if (indexq(i) .ge. 2 .and. idiff0(i) .eq. 0) then
            if (qname1(i)(14:21) .ne. tomo_phase_in) cycle
            
            ! Calculate hypocentroid of events contributing to this empirical path anomally
            x1sum = 0.
            x2sum = 0.
            x3sum = 0.
            do j = 1,indexq(i)
               iev = iqiev(i,j)
               x1sum = x1sum + latpgc(iev,it+1)
               x2sum = x2sum + lonp(iev,it+1)
               x3sum = x3sum + depthp(iev,it+1)
            end do
            lat_epa_gc = x1sum/real(indexq(i))
            call geogra (lat_epa_gc,lat_epa) ! Update in geographic coordinates
            lon_epa = x2sum/real(indexq(i))
            depth_epa = x3sum/real(indexq(i))
            call set_longitude_range (lon_epa, longitude_range)
            write (io_log,'(2i4,1x,a,1x,2f10.3,f6.1)') i, indexq(i), qname1(i), lat_epa, lon_epa, depth_epa
            write (io_log,'(4x,10i4)') (iqiev(i,j),j=1,indexq(i))
            
            ! Epicentral distance
            t1 = lat_epa_gc*rpd  ! Convert to geocentric radians
            t2 = lon_epa*rpd  ! Convert to geocentric radians
            t3 = qlat(i)*rpd  ! Convert to geocentric radians
            t4 = qlon(i)*rpd  ! Convert to geocentric radians
            call delaz (t1, t2, t3, t4, dum1, delt1, dum2, dum3, dum4, dum5, dum6, 1)
            call geogra (qlat(i), qlatdg)
            qlondg = qlon(i)
            call set_longitude_range (qlondg, 0)
            
            !  Theoretical travel-time
            call depset (depth_epa, usrc)
            if (.not.locmod) then
               call trtm (delt1)
            else
               if (verbose_log) write (io_log,'(a)') 'Calling tt_mixed_model from tomo3'
               call tt_mixed_model (depth_epa, delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
   !               if (delt1 .le. dlimlocmod .and. depth_epa .le. zlimlocmod) then
   !                  call ttloc2 (depth_epa, delt1, 'D', nphase, tt, dtdd, dtdh, dddp, phcd, ierr, 30000)
   !               else
   !                  call trtm (delt1)
   !               end if
            end if
            if (nphase .eq. 0) then
               if (verbose_log) write (io_log,'(a,f10.3)') ' tomo3: no phases returned for delt1 = ', delt1
               nphase = 1
               phcd(1) = 'CRAP    '
            end if
            ttepa = 999.99
            do j = 1,nphase
               if (tomo_phase_in .eq. phcd(j)) then
                  ttepa = tt(j) + epa(i)
                  exit
               end if
            end do
            
            ! Calibration level (Indirect trumps direct)
            if (indirect_cal) then
               igt = nint(s2gt)
            else if (direct_cal) then
               igt = nint(xl2h)
            else
               igt = 99 ! No calibration
            end if
            
            write (io_out,'(f7.3,f7.2,f7.2,1x,i3,1x,i4,2i2.2,1x,2i2.2,f5.2,2x,i4,2i2.2,1x,2i2.2,f5.2,f8.2,'//&
            '2f8.3,f6.1,1x,2f8.3,f6.1,1x,a,1x,i3,1x,i3,1x,a)')&
             delt1,& ! Epicentral distance (deg)
             epa(i),& ! Empirical path anomaly (s)
             ere(i),& ! Empirical reading error (s)
             indexq(i),& ! Number of observations of this station-phase
             9999,& ! Reading year (N/A for empirical path anomalies)
             99,& ! Reading month (N/A for empirical path anomalies)
             99,& ! Reading day (N/A for empirical path anomalies)
             99,& ! Reading hour (N/A for empirical path anomalies)
             99,& ! Reading minute (N/A for empirical path anomalies)
             99.99,& ! Reading seconds (N/A for empirical path anomalies)
             9999,& ! Origin year (N/A for empirical path anomalies)
             99,& ! Origin month (N/A for empirical path anomalies)
             99,& ! Origin day (N/A for empirical path anomalies)
             99,& ! Origin hour (N/A for empirical path anomalies)
             99,& ! Origin minute (N/A for empirical path anomalies)
             99.99,& ! Origin seconds (N/A for empirical path anomalies)
             ttepa,& ! Empirical path anomaly travel time (s)
             lat_epa,& ! Hypocentroid latitude
             lon_epa,& ! Hypocentroid longitude
             depth_epa,& ! Hypocentroid depth (km)
             qlatdg,& ! Station latitude
             qlondg,& ! Station longitude
             qelev(i),& ! Station elevation
             qname1(i),& ! Station-deployment-phase name
             0,& ! Event number = 0 for empirical path anomalies
             igt,& ! Hypocentroid calibration level (99 = no calibration)
             trim(basename) ! Title
            
         end if
      end do ! End of loop over station-phases with data
      
      return
      
   end subroutine tomo3
   
   
!*****************************************************************************************
   subroutine tomo3grid (it, tomo_phase_in, ntomo)
   
   ! Tomography output file of empirical path anomalies. Based on subroutine mlocout_rderr.
   ! The full travel time associated with an empirical path anomaly is calculated by adding the
   ! empirical path anomaly to the theoretical travel time between the hypocentroid and the
   ! appropriate station. 
   
   ! This version of tomo3 breaks up the output into cell-based subclusters. For this output
   ! I don't enforce a minimum number of samples. Some cells will have only a single sample.
   ! This value is carried in each data line so it's easy for the user to decide on a threshhold.
   ! The hypocentroid is recalculated for the events in each cell that reported each station-phase.
   
      integer, parameter :: max = 60
      
      integer, dimension(nevmax) :: cell_events ! list of event numbers for a cell
      integer :: column ! grid index for the current event
      integer :: i
      integer :: icell ! grid index for loop over entire grid
      integer :: iev ! index for event number
      integer :: igt ! calibration level
      integer :: ird ! index for phase readings in an event
      integer :: isp ! index for loop over station-phases with data
      integer, intent(in) :: it ! iteration number at convergence
      integer :: it1
      integer :: j
      integer :: jcell ! grid index for loop over entire grid
      integer :: n_cell_events ! number of events in a given cell
      integer :: n_obs ! number of observations of the desired station-phase from a cell
      integer, intent(in) :: ntomo ! index of call to command tomo, defines the grid
      integer, dimension(nevmax,2) :: readings ! event and reading numbers for observations of the station-phase in a cell
      integer :: row ! grid index for the current event
      real :: cell_hypocentroid_dep
      real :: cell_hypocentroid_lat
      real :: cell_hypocentroid_lon
      real :: delt1 ! epicentral distance
      real :: dep_test
      real :: dts
      real :: dum1 ! dummy variable for subroutine delaz
      real :: dum2 ! dummy variable for subroutine delaz
      real :: dum3 ! dummy variable for subroutine delaz
      real :: dum4 ! dummy variable for subroutine delaz
      real :: dum5 ! dummy variable for subroutine delaz
      real :: dum6 ! dummy variable for subroutine delaz
      real :: epa_cell
      real :: lat_test
      real :: lon_test
      real :: qlatdg
      real :: qlondg
      real, dimension(nevmax):: residual
      real :: spread
      real :: sum_data
      real :: t1 ! input for subroutine delaz
      real :: t2 ! input for subroutine delaz
      real :: t3 ! input for subroutine delaz
      real :: t4 ! input for subroutine delaz
      real :: ttepa ! sum of theoretical travel time and empirical path anomaly
      real, dimension(2) :: usrc
      real :: x1sum ! sum of latitudes for hypocentroid
      real :: x2sum ! sum of longitudes for hypocentroid
      real :: x3sum ! sum of depths for hypocentroid
      real :: xlat ! cell hypocentroid latitude in geocentric coordinates
      real :: xlon ! cell hypocentroid longitude in geocentric coordinates
      character(len=8), intent(in) :: tomo_phase_in ! phase name desired for tomographic output
      
      write (io_log,'(/a)') 'tomo3grid:'
      
      it1 = it + 1
      
      n_cell_events = 0
      do icell = 1, n_cell_lat(ntomo) ! loop over rows
         do jcell = 1, n_cell_lon(ntomo) ! loop over columns
            write (io_log,'(t3,a,2i4)') 'cell coordinates = ', icell, jcell
            do iev = 1, nev
               call cell_indices (it, iev, ntomo, row, column) ! get cell indices for this event
               if (row .eq. icell .and. column .eq. jcell) then
                  n_cell_events = n_cell_events + 1
                  cell_events(n_cell_events) = iev
                  write (io_log,'(t6,i4,1x,a)') iev, evtnam(iev)
               end if
            end do
            write (io_log,'(t3,i3,a,2i4)') n_cell_events, ' events in cell ', icell, jcell
            if (n_cell_events .eq. 0) cycle
            
            write (io_log,'(t3,a)') 'Loop over station-phases with data'
            do isp = 1, nqc ! loop over station-phases with data
               if (idiff0(isp) .ne. 0) cycle ! No differential time phases
               if (qname1(isp)(14:21) .ne. tomo_phase_in) cycle
               write (io_log,'(t3,i5,1x,a)') isp, qname1(isp)
               n_obs = 0
               do i = 1, n_cell_events ! loop over events in the cell
                  iev = cell_events(i)
                  do ird = 1, nst(iev) ! loop over all readings in this event
                     if (fcode(iev,ird) .ne. ' ') cycle ! don't use flagged readings
                     if (stname(iev,ird) .eq. qname1(isp)(1:5) .and. phase(iev,ird) .eq. tomo_phase_in) then
                        if (n_obs .lt. nevmax) then
                           n_obs = n_obs + 1
                           readings(n_obs,1) = iev
                           readings(n_obs,2) = ird
                           if (.not.fltrhres(readings(n_obs,1),readings(n_obs,2)) .and.&
                               .not.fltrcres(readings(n_obs,1),readings(n_obs,2))) then ! Not filtered for large residual ("PRES")
                              if (indirect_cal) then
                                 dts = dt_cal(readings(n_obs,1),readings(n_obs,2)) - s_cal(readings(n_obs,1),readings(n_obs,2))
                              else
                                 dts = dt(readings(n_obs,1),readings(n_obs,2),it) - s(readings(n_obs,1),readings(n_obs,2),it)
                              end if
                              ! Check for unflagged gross outliers caused by calibration shift in indirect calibration
                              if (abs(dts) .lt. 20.) then
                                 write (io_log,'(t6,i3,i4,1x,a,i6,1x,a,1x,a)') n_obs,&
                                  iev, evtnam(iev), ird, stname(iev,ird), phase(iev,ird)
                              else
                                 n_obs = n_obs - 1
                              end if
                           else
                              n_obs = n_obs - 1
                           end if
                        else if (n_obs .eq. nevmax) then
                           call warnings ('tomo3grid: maximum number of observations (nevmax) reached for '//qname1(isp))
                        end if
                     end if
                  end do ! end loop over all readings in this event
               end do ! end loop over events in the cell
               if (n_obs .eq. 0) cycle
               
               ! calculate hypocentroid of events that have this station-phase in this cell
               write (io_log,'(t3,3a,2i4)') 'Calculate hypocentroid for ', qname1(isp), ' in cell ', icell, jcell
               x1sum = 0.
               x2sum = 0.
               x3sum = 0.
               do j = 1, n_obs
                  iev = readings(j,1)
                  if (indirect_cal) then ! indirect calibration
                     lat_test = latp_cal(iev) ! Geographic latitude
                     lon_test = lonp_cal(iev) ! Geographic longitude
                     call set_longitude_range (lon_test, longitude_range)
                     dep_test = depthp_cal(iev)
                  else if (direct_cal) then ! direct calibration
                     lat_test = latp(iev,it1) ! Geographic latitude
                     lon_test = lonp(iev,it1) ! Geographic longitude
                     call set_longitude_range (lon_test, longitude_range)
                     dep_test = depthp(iev,it+1)
                  else ! uncalibrated
                     lat_test = latp(iev,it1) ! Geographic latitude
                     lon_test = lonp(iev,it1) ! Geographic longitude
                     call set_longitude_range (lon_test, longitude_range)
                     dep_test = depthp(iev,it+1)
                  end if
                  x1sum = x1sum + lat_test
                  x2sum = x2sum + lon_test
                  x3sum = x3sum + dep_test
               end do
               cell_hypocentroid_lat = x1sum/real(n_obs)
               cell_hypocentroid_lon = x2sum/real(n_obs)
               call set_longitude_range (cell_hypocentroid_lon, longitude_range)
               cell_hypocentroid_dep = x3sum/real(n_obs)
               write (io_log,'(t6,i3,3f10.3)') n_obs, cell_hypocentroid_lat, cell_hypocentroid_lon, cell_hypocentroid_dep
               
               ! Convert cell hypocentroid to geocentric coordinates
               call geocen (cell_hypocentroid_lat, cell_hypocentroid_lon, xlat, dum1, dum2, xlon, dum3, dum4)
               
               t1 = xlat*rpd  ! Convert hypocentroid latitude to geocentric radians
               t2 = xlon*rpd  ! Convert hypocentroid longitude to geocentric radians
               t3 = qlat(isp)*rpd  ! Convert station latitude (already geocentric) to geocentric radians
               t4 = qlon(isp)*rpd  ! Convert station longitude (already geocentric) to geocentric radians
               call delaz (t1, t2, t3, t4, dum1, delt1, dum2, dum3, dum4, dum5, dum6, 1)
               call geogra (qlat(isp), qlatdg)
               qlondg = qlon(isp)
               call set_longitude_range (qlondg, 0)
               
               ! Empirical path anomaly and spread of the group of observations
               if (n_obs .eq. 1) then
                  if (indirect_cal) then
                     dts = dt_cal(readings(1,1),readings(1,2)) - s_cal(readings(1,1),readings(1,2))
                  else
                     dts = dt(readings(1,1),readings(1,2),it) - s(readings(1,1),readings(1,2),it)
                  end if
                  epa_cell = dts ! empirical path anomaly is just the residual
                  spread = ere(isp) ! empirical reading error of all observations of this station-phase
               else
                  sum_data = 0.
                  do j = 1, n_obs
                     if (indirect_cal) then
                        dts = dt_cal(readings(j,1),readings(j,2)) - s_cal(readings(j,1),readings(j,2))
                     else
                        dts = dt(readings(j,1),readings(j,2),it) - s(readings(j,1),readings(j,2),it)
                     end if
                     residual(j) = dts
                     sum_data = sum_data + residual(j)
                  end do
                  epa_cell = sum_data/real(n_obs)
                  call croux (residual, min(n_obs,1000), spread) ! Robust scale estimator Sn
               end if
               
               !  Theoretical travel-time
               call depset (cell_hypocentroid_dep, usrc)
               if (.not.locmod) then
                  call trtm (delt1)
               else
                  call tt_mixed_model (cell_hypocentroid_dep, delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
               end if
               if (nphase .eq. 0) then
                  nphase = 1
                  phcd(1) = 'CRAP    '
               end if
               ttepa = 999.99
               do j = 1,nphase
                  if (tomo_phase_in .eq. phcd(j)) then
                     ttepa = tt(j) + epa_cell
                     exit
                  end if
               end do
               
               ! Calibration level (Indirect trumps direct)
               if (indirect_cal) then
                  igt = nint(s2gt)
               else if (direct_cal) then
                  igt = nint(xl2h)
               else
                  igt = 99 ! No calibration
               end if
               
               write (io_out,'(f7.3,f7.2,f7.2,1x,i3,1x,i4,2i2.2,1x,2i2.2,f5.2,2x,i4,2i2.2,1x,2i2.2,f5.2,f8.2,'//&
               '2f8.3,f6.1,1x,2f8.3,f6.1,1x,a,1x,i3,1x,i3,1x,a)')&
                delt1,& ! Epicentral distance (deg)
                epa_cell,& ! Empirical path anomaly (s)
                spread,& ! Empirical reading error (s)
                n_obs,& ! Number of observations of this station-phase
                9999,& ! Reading year (N/A for empirical path anomalies)
                99,& ! Reading month (N/A for empirical path anomalies)
                99,& ! Reading day (N/A for empirical path anomalies)
                99,& ! Reading hour (N/A for empirical path anomalies)
                99,& ! Reading minute (N/A for empirical path anomalies)
                99.99,& ! Reading seconds (N/A for empirical path anomalies)
                9999,& ! Origin year (N/A for empirical path anomalies)
                99,& ! Origin month (N/A for empirical path anomalies)
                99,& ! Origin day (N/A for empirical path anomalies)
                99,& ! Origin hour (N/A for empirical path anomalies)
                99,& ! Origin minute (N/A for empirical path anomalies)
                99.99,& ! Origin seconds (N/A for empirical path anomalies)
                ttepa,& ! Empirical path travel time (s)
                cell_hypocentroid_lat,& ! Hypocentroid latitude
                cell_hypocentroid_lon,& ! Hypocentroid longitude
                cell_hypocentroid_dep,& ! Hypocentroid depth (km)
                qlatdg,& ! Station latitude
                qlondg,& ! Station longitude
                qelev(isp),& ! Station elevation
                qname1(isp),& ! Station-deployment-phase name
                0,& ! Event number = 0 for empirical path anomalies
                igt,& ! Hypocentroid calibration level (99 = no calibration)
                trim(basename) ! Title
                
            end do ! end loop over station-phases with data
            
         end do ! end loop over columns
      end do ! end loop over rows
      
      return
      
   end subroutine tomo3grid
   
   
!*****************************************************************************************
   subroutine tomo_grid (it, ntomo)
   
   ! Defines a grid of cells for tomographic output if the option to subdivide a large cluster
   ! has been chosen in a particular call to command "tomo". The grid can be different for
   ! each call if the requested size of cells is different. The origin of the grid is always
   ! the same, at the lower left corner, based on the minimum latitude of any event and the
   ! minimum longitude of any event in the cluster. 
   
      integer :: iev
      integer, intent(in) :: it
      integer :: it1
      integer :: n_cells
      integer, intent(in) :: ntomo ! index of call to command tomo
      real :: latmax
      real :: latmin
      real :: lon_test
      real :: lonmax
      real :: lonmin
      character(len=132) :: msg
      
      it1 = it + 1
      
      ! Get geographic limits of the cluster
      latmin = 90.
      latmax = -90.
      lonmin = 360
      lonmax = -180
      do iev = 1,nev
         if (indirect_cal) then ! indirect calibration
            lon_test = lonp_cal(iev)
            call set_longitude_range (lon_test, longitude_range)
            if (latp_cal(iev) .lt. latmin) latmin = latp_cal(iev) ! Geographic latitude.
            if (lon_test .lt. lonmin) lonmin = lon_test ! Geographic longitude.
            if (latp_cal(iev) .gt. latmax) latmax = latp_cal(iev) ! Geographic latitude.
            if (lon_test .gt. lonmax) lonmax = lon_test ! Geographic longitude.
         else if (direct_cal) then ! direct calibration
            lon_test = lonp(iev,it1)
            call set_longitude_range (lon_test, longitude_range)
            if (latp(iev,it1) .lt. latmin) latmin = latp(iev,it1) ! Geographic latitude.
            if (lon_test .lt. lonmin) lonmin = lon_test ! Geographic longitude.
            if (latp(iev,it1) .gt. latmax) latmax = latp(iev,it1) ! Geographic latitude.
            if (lon_test .gt. lonmax) lonmax = lon_test ! Geographic longitude.
         else ! uncalibrated
            lon_test = lonp(iev,it1)
            call set_longitude_range (lon_test, longitude_range)
            if (latp(iev,it1) .lt. latmin) latmin = latp(iev,it1) ! Geographic latitude.
            if (lon_test .lt. lonmin) lonmin = lon_test ! Geographic longitude.
            if (latp(iev,it1) .gt. latmax) latmax = latp(iev,it1) ! Geographic latitude.
            if (lon_test .gt. lonmax) lonmax = lon_test ! Geographic longitude.
         end if
      end do
      grid_origin_lat(ntomo) = latmin - 0.1
      grid_origin_lon(ntomo) = lonmin - 0.1
      
      n_cell_lat(ntomo) = ceiling((latmax - latmin) / cell_lat_inc(ntomo))
      n_cell_lon(ntomo) = ceiling((lonmax - lonmin) / cell_lon_inc(ntomo))
      n_cells = n_cell_lat(ntomo)*n_cell_lon(ntomo)
      
      if (n_cells .gt. n_cell_max) then
         write (msg,'(4(a,i4),a)') 'tomo_grid: ', n_cell_lat(ntomo), ' x ', n_cell_lon(ntomo),&
          ' = ', n_cells, '; exceeds the maximum number (', n_cell_max, ')'
         call oops (trim(msg))
      else
         write (io_log,'(3(a,i4))') 'tomo_grid: ', n_cell_lat(ntomo), ' x ', n_cell_lon(ntomo),&
          ' = ', n_cells
      end if
      
      write (io_log,'(a,2f10.3)') 'tomo_grid: grid origin = ', grid_origin_lat(ntomo), grid_origin_lon(ntomo)
      write (io_log,'(a,i2,a,f5.2)') 'tomo_grid: ', n_cell_lat(ntomo), ' cells in latitude; increment = ', cell_lat_inc(ntomo)
      write (io_log,'(a,i2,a,f5.2)') 'tomo_grid: ', n_cell_lon(ntomo), ' cells in longitude; increment = ', cell_lon_inc(ntomo)
      
      return
      
   end subroutine tomo_grid
   
   
!*****************************************************************************************
   subroutine cell_indices (it, iev, ntomo, row, column)
   
   ! Determines the indices (row and column) of the cell occupied by an event. The grid, and
   ! therefore the cell indices, can be different for different calls to command tomo.
   
      integer, intent(out) :: column ! of grid(ntomo)
      integer, intent(in) :: iev ! event index
      integer, intent(in) :: it ! iteration at convergence
      integer :: it1
      integer, intent(in) :: ntomo ! Index for grid associated with a call to command tomo
      integer, intent(out) :: row ! of grid(ntomo)
      real :: lat_test
      real :: lon_test
      
      it1 = it + 1
      if (indirect_cal) then ! indirect calibration
         lon_test = lonp_cal(iev)
         call set_longitude_range (lon_test, longitude_range)
         lat_test = latp_cal(iev) ! Geographic latitude.
      else if (direct_cal) then ! direct calibration
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         lat_test = latp(iev,it1) ! Geographic latitude.
      else ! uncalibrated
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         lat_test = latp(iev,it1) ! Geographic latitude.
      end if
      
      row = ceiling((lat_test - grid_origin_lat(ntomo)) / cell_lat_inc(ntomo))
      column = ceiling((lon_test - grid_origin_lon(ntomo)) / cell_lon_inc(ntomo))
      
      write (io_log,'(a,i3,a,2i4)') 'cell_indices for event ', iev, ' = ', row, column
      
      return
      
   end subroutine cell_indices
   
   
!*****************************************************************************************
end module mlocout_tomo

