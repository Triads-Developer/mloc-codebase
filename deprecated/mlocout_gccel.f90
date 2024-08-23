!> GCCEL output files

module mlocout_gccel

   use declare_calibration
   use declare_environment
   use declare_lun
   use declare_output
   use declare_tt
   use mloclib
   use mloclib_messages
   
   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine gccel (it)
   
   ! Write ~_gccel.dat file for distribution in GCCEL. Similar to a _datf.mnf file, with
   ! addition of information about crustal model, empirical reading errors, station
   ! coordinates, and commentary.
   
   ! 11/16/2017: update to MNF v1.4.2, adding deployment code to phase records
   
      integer :: ecl
      integer :: eclmax
      integer :: eclmin
      integer :: i
      integer :: ielev_m
      integer :: iev
      integer :: io_commentary
      integer :: ios
      integer, intent(in) :: it
      integer :: it1
      integer, external :: lunit
      integer, dimension(8) :: values
      real :: depmax
      real :: depmin
      real :: glat
      real :: glon
      real :: latmax
      real :: latmin
      real :: lon_test
      real :: lonmax
      real :: lonmin
      real :: magmax
      real :: magmin
      character(len=4) :: cal_code
      character(len=100) :: commentary_path
      character(len=8) :: date
      character(len=8) :: datefirst
      character(len=8) :: datelast
      character(len=132) :: line132
      character(len=34) :: line34
      character(len=132) :: msg
      character(len=100) :: outfil
      character(len=10) :: time
      character(len=5) :: zone
   
      call fyi ('mlocout_gccel: writing ~_gccel.dat file')
      outfil = trim(gcat_folder)//dirsym//trim(basename)//'_gccel.dat'
      call open_file (io_gccel, outfil, 'new')
      write (io_gccel,'(a,3x,a,t121,a)') 'B', 'Calibrated Earthquake Cluster: '//trim(basename), ' '
      write (io_gccel,'(a,t121,a)') 'F   1.4.2', ' '
   
      ! Attributions
      write (io_gccel,'(2a)') '# Program: ', trim(mloc_version)
      write (io_gccel,'(2a)') '# Run: ', trim(basename)
      write (io_gccel,'(2a)') '# Author: ', trim(mloc_author)
      call date_and_time (date, time, zone, values)
      write (io_gccel,'(4a)') '# Date: ', date, ' ', time
      write (io_gccel,'(a,t121,a)') '# ', ' '
   
      ! Add a commentary file, if it exists. The filename (assumed to be in the data directory)
      ! is given as an argument to the gcat command. The file should be hard-wrapped
      ! at a reasonable line length (say, 72 characters), with a maximum of 132 characters.
      commentary_path = trim(datadir)//dirsym//trim(commentary_fname)
      io_commentary = lunit()
      call open_file (io_commentary, commentary_path, 'old')
      do
         read (io_commentary,'(a)',iostat=ios) line132
         if (ios .lt. 0) exit
         write (io_gccel,'(2a)') '# ', trim(line132)
      end do
      close (io_commentary)
      write (io_gccel,'(a,t121,a)') '# ', ' '
   
      write (io_gccel,'(a,i3)') '# Number of events: ', nev
      write (commentary_buffer(1),'(a,i3)') 'Number of events: ', nev
      
      ! Calibration type and calibration statistics
      eclmin = 99
      eclmax = 0
      if (indirect_cal) then ! Indirect calibration
         write (io_gccel,'(a,i3,a,f4.1,a)') '# Calibration type: indirect calibration on ', ncal(3),&
          ' calibration events; hypocentroid calibration level = ', cal_level, ' km'
         write (commentary_buffer(2),'(a,i3,a)') 'Calibration type: indirect calibration on ', ncal(3),&
          ' calibration events'
         write (commentary_buffer(3),'(a,f4.1,a)') 'Hypocentroid calibration level = ', cal_level, ' km'
         do iev = 1,nev
            call calibration_code (iev, 2, cal_code)
            read (cal_code(3:4),'(i2)') ecl
            if (ecl .lt. eclmin) eclmin = ecl
            if (ecl .gt. eclmax) eclmax = ecl
         end do
         write (io_gccel,'(a,i2,a,i2,a)') '# Epicentral calibration range: ', eclmin, ' - ', eclmax, ' km'
         write (commentary_buffer(4),'(a,i2,a,i2,a)') 'Epicentral calibration range: ', eclmin, ' - ', eclmax, ' km'
      else if (direct_cal) then ! Direct calibration
         write (io_gccel,'(a,f4.1,a,f4.1,a)') '# Calibration type: direct calibration using data to ',&
          hlim(1,2), ' degrees; hypocentroid calibration level = ', xl2h, ' km'     
         write (commentary_buffer(2),'(a,f4.1,a,f4.1,a)') 'Calibration type: direct calibration using data to ',&
          hlim(1,2), ' degrees'     
         write (commentary_buffer(3),'(a,f4.1,a)') 'Hypocentroid calibration level = ', xl2h, ' km'     
         do iev = 1,nev
            call calibration_code (iev, 1, cal_code)
            read (cal_code(3:4),'(i2)') ecl
            if (ecl .lt. eclmin) eclmin = ecl
            if (ecl .gt. eclmax) eclmax = ecl
         end do
         write (io_gccel,'(a,i2,a,i2,a)') '# Epicentral calibration range: ', eclmin, ' - ', eclmax, ' km'
         write (commentary_buffer(4),'(a,i2,a,i2,a)') 'Epicentral calibration range: ', eclmin, ' - ', eclmax, ' km'
      else
         write (io_gccel,'(a)') '# Calibration type: uncalibrated'
         write (commentary_buffer(2),'(a)') 'Calibration type: uncalibrated'
         write (commentary_buffer(3),'(a)') 'Hypocentroid calibration level = N/A'
         write (commentary_buffer(4),'(a)') 'Epicentral calibration range: N/A'
         msg = 'mlocout_gccel: this cluster is uncalibrated. GCCEL output is normally only invoked for calibrated clusters.'
         call warnings (trim(msg))
      end if
      
      ! Date range (assumes events are in chronological order)
      write (datefirst,'(i4,2i2.2)') iyre(1), mone(1), idye(1)
      write (datelast,'(i4,2i2.2)') iyre(nev), mone(nev), idye(nev)
      write (io_gccel,'(4a)') '# Date range: ', datefirst, ' - ', datelast
      write (commentary_buffer(5),'(4a)') 'Date range: ', datefirst, ' - ', datelast
      
      ! Latitude, longitude, depth and magnitude range
      latmin = 90.
      latmax = -90.
      lonmin = 360
      lonmax = -180
      depmin = 999.
      depmax = 0.
      magmin = 10.
      magmax = 0.
      it1 = it + 1
      do iev = 1,nev
         if (indirect_cal) then ! indirect calibration
            lon_test = lonp_cal(iev)
            call set_longitude_range (lon_test, longitude_range)
            if (latp_cal(iev) .lt. latmin) latmin = latp_cal(iev) ! Geographic latitude.
            if (lon_test .lt. lonmin) lonmin = lon_test ! Geographic longitude.
            if (latp_cal(iev) .gt. latmax) latmax = latp_cal(iev) ! Geographic latitude.
            if (lon_test .gt. lonmax) lonmax = lon_test ! Geographic longitude.
            if (depthp_cal(iev) .lt. depmin) depmin = depthp_cal(iev)
            if (depthp_cal(iev) .gt. depmax) depmax = depthp_cal(iev)
         else if (direct_cal) then ! direct calibration
            lon_test = lonp(iev,it1)
            call set_longitude_range (lon_test, longitude_range)
            if (latp(iev,it1) .lt. latmin) latmin = latp(iev,it1) ! Geographic latitude.
            if (lon_test .lt. lonmin) lonmin = lon_test ! Geographic longitude.
            if (latp(iev,it1) .gt. latmax) latmax = latp(iev,it1) ! Geographic latitude.
            if (lon_test .gt. lonmax) lonmax = lon_test ! Geographic longitude.
            if (depthp(iev,it1) .lt. depmin) depmin = depthp(iev,it1)
            if (depthp(iev,it1) .gt. depmax) depmax = depthp(iev,it1)
         else ! uncalibrated
            lon_test = lonp(iev,it1)
            call set_longitude_range (lon_test, longitude_range)
            if (latp(iev,it1) .lt. latmin) latmin = latp(iev,it1) ! Geographic latitude.
            if (lon_test .lt. lonmin) lonmin = lon_test ! Geographic longitude.
            if (latp(iev,it1) .gt. latmax) latmax = latp(iev,it1) ! Geographic latitude.
            if (lon_test .gt. lonmax) lonmax = lon_test ! Geographic longitude.
            if (depthp(iev,it1) .lt. depmin) depmin = depthp(iev,it1)
            if (depthp(iev,it1) .gt. depmax) depmax = depthp(iev,it1)
         end if
         if (rmag(iev) .gt. 0. .and. rmag(iev) .lt. magmin) magmin = rmag(iev)
         if (rmag(iev) .gt. magmax) magmax = rmag(iev)
      end do
      write (io_gccel,'(a,f7.3,a,f7.3)') '# Latitude range: ', latmin, ' - ', latmax
      write (commentary_buffer(6),'(a,f7.3,a,f7.3)') 'Latitude range: ', latmin, ' - ', latmax
      write (io_gccel,'(a,f8.3,a,f8.3)') '# Longitude range: ', lonmin, ' - ', lonmax
      write (commentary_buffer(7),'(a,f8.3,a,f8.3)') 'Longitude range: ', lonmin, ' - ', lonmax
      write (io_gccel,'(a,f5.1,a,f5.1)') '# Depth range: ', depmin, ' - ', depmax
      write (commentary_buffer(8),'(a,f5.1,a,f5.1)') 'Depth range: ', depmin, ' - ', depmax
      if (magmax .gt. 0.) then
         write (io_gccel,'(a,f3.1,a,f3.1)') '# Magnitude range: ', magmin, ' - ', magmax
         write (commentary_buffer(9),'(a,f3.1,a,f3.1)') 'Magnitude range: ', magmin, ' - ', magmax
      else
         write (io_gccel,'(a)') '# Magnitude range: not available'
         write (commentary_buffer(9),'(a)') 'Magnitude range: not available'
      end if
   
      ! Local velocity model
      if (locmod) then
         call open_file (io_locmod, trim(locmodfname), 'old')
         ! Epicentral range line
         read (io_locmod,'(a)',iostat=ios) line34
         ! Layer lines
         do
            read (io_locmod,'(a)',iostat=ios) line34
            if (ios .lt. 0) exit
            write (io_gccel,'(2a)') 'L   ', line34(1:30)
         end do
         close (io_locmod)
      end if
      
      ! Station coordinates
      write (io_gccel,'(a)') '# Stations used.'
      do i = 1,nkstat
         call geogra (stalat(kstat(i)), glat) ! Geographic latitude
         glon = stalon(kstat(i))
         call set_longitude_range (glon, 0)
         ielev_m = int(ahgtr(kstat(i))*1.0e3) ! Elevation in m
         write (io_gccel,'(2a,1x,a5,1x,a8,1x,2f10.4,i8)') 'C ', nstr1(kstat(i)), sta_agency(kstat(i)),&
          sta_deployment(kstat(i)), glat, glon, ielev_m
      end do
      
      ! Event data
      do iev = 1,nev
         call write_mnf_14 (io_gccel, iev, it)
      end do
      write (io_gccel,'(a,t121,a)') 'EOF', ' '
      close (io_gccel)
   
      return
      
   end subroutine gccel
   
   
!*****************************************************************************************
end module mlocout_gccel

