!> Procedures related to seismograph stations.

module mloc_stations

   use mloc_declare
   use mloc_math
   use mloc_lib
   
   implicit none
   save   

   ! ISC-FDSN station metadata file
   integer :: n_ifsm_read = 0 ! Entries read from the main ISC-FDSN station metadata file
   integer :: n_ifsm_written = 0 ! Matched entries written to the run-specific output file
   integer, allocatable, dimension(:) :: ifsm_off ! Operational epoch off, Year + Julian date
   integer, dimension(n_ifsm_max) :: ifsm_off_matched ! Matched entries, operational epoch off, Year + Julian date
   integer, allocatable, dimension(:) :: ifsm_on ! Operational epoch on, Year + Julian date
   integer, dimension(n_ifsm_max) :: ifsm_on_matched ! Matched entries, operational epoch on, Year + Julian date
   character(len=7), allocatable, dimension(:) :: ds ! Deployment, Station Code, main file
   character(len=7), dimension(n_ifsm_max) :: ds_matched ! Deployment, Station Code, matched entries
   character(len=100) :: ifsm_stn_file = 'ifsm_stn_file.dat' ! Output file for a supplemental station file
   character(len=94), allocatable, dimension(:) :: ifsm_line
   
   ! NEIC station metadata file 
   integer :: n_nsmd_read = 0 ! Entries read from the main NEIC station metadata file
   integer :: n_nsmd_written = 0 ! Matched entries written to the run-specific output file
   integer, allocatable, dimension(:) :: nsmd_off ! Operational epoch off, Year + Julian date
   integer, dimension(n_nsmd_max) :: nsmd_off_matched ! Matched entries, operational epoch off, Year + Julian date
   integer, allocatable, dimension(:) :: nsmd_on ! Operational epoch on, Year + Julian date
   integer, dimension(n_nsmd_max) :: nsmd_on_matched ! Matched entries, operational epoch on, Year + Julian date
   character(len=12), allocatable, dimension(:) :: dslc ! Deployment, Station Code, Location, Channel, main file
   character(len=12), dimension(n_nsmd_max) :: dslc_matched ! Deployment, Station Code, Location, Channel, matched entries
   character(len=100) :: nsmd_stn_file = 'nsmd_stn_file.dat' ! Output file for a supplemental station file
   character(len=102), allocatable, dimension(:) :: nsmd_line

contains

!*****************************************************************************************
   subroutine stafind (iev, ird)
   
   ! Find coordinates for a given station code. The algorithm now makes use of date codes
   ! for the start and end of operations of a given station at specified coordinates, allowing
   ! for changes in station location through time. If the deployment field for a reading is
   ! non-blank, then the matching algorithm honors it and requires a match for the deployment,
   ! as well as the station code, in the station list. Agency is not considered.
   
      integer :: i
      integer, intent(in) :: iev ! event number
      integer, intent(in) :: ird ! reading number
      integer :: juldat_test
      integer :: julday
      integer :: k
      integer :: k_match
      character(len=132) :: msg
!      logical :: date_range
      logical :: failed
      logical :: match
   
      failed = .false.
   
      ! Julian date for current event
      call juldat (iyre(iev), mone(iev), idye(iev), julday, 0)
      juldat_test = iyre(iev)*1000 + julday
   
      ! Search station data list
      do k = 1,nstat1
         if (nstr1(k) .eq. sad(iev,ird)(1:5) .and. sta_deployment(k) .eq. sad(iev,ird)(13:20)) then
            k_match = k
            if (date_range(juldat_test,jdate_on(k),jdate_off(k))) then
               stladg(iev,ird) = stalat(k) ! Geocentric coordinates
               stlndg(iev,ird) = stalon(k) ! Geocentric coordinates
               ahgts(iev,ird) = ahgtr(k)
               sca0s(iev,ird) = sca0(k)
               sd1s(iev,ird) = sd1(k)
               kcode(iev,ird) = k
               call stations_used (k)
               return
            else
               failed = .true.
            end if
         end if
      end do
      
      if (failed) then
         n_failed_date_range = n_failed_date_range + 1
         write (msg,'(a,2i4,1x,a,1x,3i8)') 'stafind: failed date_range: ',&
          iev, ird, stname(iev,ird), juldat_test, jdate_on(k_match), jdate_off(k_match)
         if (verbose_screen) call warnings (trim(msg))
         write (io_stn_log,'(3x,a)') msg
      end if
   
      ! station stname(iev,ird) is not in the station set.
      if (verbose_screen) call warnings ('stafind: missing station code: '//sad(iev,ird))
      kcode(iev,ird) = -1
      fcode(iev,ird) = 'm'
      stladg(iev,ird) = 0.
      stlndg(iev,ird) = 0.
      ahgts(iev,ird) = 0.
      sca0s(iev,ird) = 0.
      sd1s(iev,ird) = 0.
      if (nmiss(iev) .lt. n_miss_sta_max) then
         nmiss(iev) = nmiss(iev) + 1
      else
         write (msg,'(a,i4)') 'stafind: maximum number of missing stations reached: ', n_miss_sta_max
         call warnings (trim(msg))
      end if
      
      missta(iev,nmiss(iev)) = sad(iev,ird)
      if (n_miss_sta_total .eq. 0) then
         n_miss_sta(1) = 1
         n_miss_sta_total = 1
         n_miss_sta_list(1) = missta(iev,nmiss(iev))
      else
         match = .false.
         do i = 1,n_miss_sta_total
            if (sad(iev,ird) .eq. n_miss_sta_list(i)) then
               match = .true.
               n_miss_sta(i) = n_miss_sta(i) + 1
            end if
         end do
         if (.not.match) then
            if (n_miss_sta_total .lt. n_miss_sta_max) then
               n_miss_sta_total = n_miss_sta_total + 1
               n_miss_sta(n_miss_sta_total) = 1
               n_miss_sta_list(n_miss_sta_total) = sad(iev,ird)
            else
               write (msg,'(a,i4)') 'stafind: maximum number of missing stations reached: ', n_miss_sta_max
               call warnings (trim(msg))
            end if
         end if
      end if
   
      return
      
   end subroutine stafind
   
   
!*****************************************************************************************
   subroutine stations_used (kindex)
   
   ! Keeps a record of all stations in the dataset for which station coordinates were found
   
      integer :: i
      integer, intent(in) :: kindex ! index from the station list
      character(len=132) :: msg
   
      if (nkstat .gt. 0) then
         if (nkstat .gt. n_sta_max) then
            msg = 'stations_used: n_sta_max exceeded'
            call oops (trim(msg))
         end if
         do i = 1,nkstat
            if (kindex .eq. kstat(i)) return
         end do
         nkstat = nkstat + 1
         kstat(nkstat) = kindex
      else if (nkstat .eq. 0) then
         nkstat = 1
         kstat(1) = kindex
      else
         write (msg,'(a,i8)') 'stations_used: illegal value for kindex: ', kindex
         call oops (trim(msg))
      end if
   
      return
      
   end subroutine stations_used
   
   
!*****************************************************************************************
   subroutine supp_stations_used ()
   
   ! List stations from supplemental station files that are used in the dataset. The list is
   ! written in the generic supplemental station file format (isstn = 3), including the header
   ! line, so this section can be copied and pasted into a new supplemental station file that
   ! is specific to this cluster.
   
      integer :: i
      integer :: ielev_m
      integer :: ioff
      integer :: ion
      integer :: j
      real :: glat
      real :: glon
      character(len=5) :: supp_agency
      character(len=8) :: supp_deployment
      character(len=5) :: supp_sta
      logical :: used
      
      if (.not.suppstn) return
      
      write (io_stn_log,'(/a)') 'Stations from supplemental station files that are present in the dataset, generic format:'
      write (io_stn_log,'(a,a)') '3 Supplemental station file for ', trim(basename)
      
      do i = 1,n_supp_stn
         used = .false.
         do j = 1,nkstat
            if (i .eq. kstat(j)) then
               used = .true.
               supp_sta = nstr1(kstat(j))
               supp_agency = sta_agency(kstat(j))
               supp_deployment = sta_deployment(kstat(j))
               call geogra (stalat(kstat(j)), glat) ! Geographic latitude
               glon = stalon(kstat(j))
               call set_longitude_range (glon, 0)
               ielev_m = int(ahgtr(kstat(j))*1.0e3) ! Elevation in m, corrected for depth of burial
               ion = jdate_on(kstat(j))
               ioff = jdate_off(kstat(j))
               exit
            end if
         end do
         if (used) call supp_station_file_write (io_stn_log, supp_sta, supp_agency, supp_deployment,&
          glat, glon, ielev_m, 0, ion, ioff, sta_author(i))
      end do
   
      return
      
   end subroutine supp_stations_used
   
   
!*****************************************************************************************
   subroutine supp_stations_not_used ()
   
   ! List stations from supplemental station files that were not used in the dataset. The list is
   ! written in the generic supplemental station file format (isstn = 3), but without the header
   ! line, since there is no good reason to use this as a supplemental station file.
   
      integer :: i
      integer :: ielev_m
      integer :: ioff
      integer :: ion
      integer :: j
      real :: glat
      real :: glon
      character(len=5) :: supp_agency
      character(len=8) :: supp_deployment
      character(len=5) :: supp_sta
      logical :: used
      
      if (.not.suppstn) return
      
      write (io_stn_log,'(/a)') 'Stations from supplemental station files that are not present in the dataset:'
      
      do i = 1,n_supp_stn
         used = .false.
         do j = 1,nkstat
            if (i .eq. kstat(j)) then
               used = .true.
               exit
            end if
         end do
         if (.not.used) then
            supp_sta = nstr1(i)
            supp_agency = sta_agency(i)
            supp_deployment = sta_deployment(i)
            call geogra (stalat(i), glat) ! Geographic latitude
            glon = stalon(i)
            call set_longitude_range (glon, 0)
            ielev_m = int(ahgtr(i)*1.0e3) ! Elevation in m, corrected for depth of burial
            ion = jdate_on(i)
            ioff = jdate_off(i)
            call supp_station_file_write (io_stn_log, supp_sta, supp_agency, supp_deployment,&
             glat, glon, ielev_m, 0, ion, ioff, sta_author(i))
         end if
      end do
   
      return
      
   end subroutine supp_stations_not_used
   
   
!*****************************************************************************************
   subroutine redsta ()
   
   ! Read station locations, elevations, and patch corrections from the main station data file and any
   ! supplemental station data files that have been declared. Supplemental station data files are read
   ! first so that they take precedence when searching the list of stations. Several formats are
   ! defined for supplemental station data files.
   
   ! All latitudes and longitudes are returned in geocentric coordinates.
   
      integer :: day_oe
      integer :: dob
      integer :: i
      integer :: increment
      integer :: ios
      integer :: isstn
      integer :: k
      integer :: latdeg
      integer :: latmin
      integer :: latsec
      integer :: londeg
      integer :: lonmin
      integer :: lonsec
      integer :: lun
!      integer :: lunit
      integer :: month_oe
      integer :: n
      integer :: n_master
      integer :: station_elevation
      integer :: year_oe
      real :: cxlat
      real :: cxlon
      real :: glat
      real :: glon
      real :: sxlat
      real :: sxlon
      real :: xlat
      real :: xlatmin
      real :: xlatsec
      real :: xlon
      real :: xlonmin
      real :: xlonsec
      character(len=1) :: char1
      character(len=132) :: line
      character(len=132) :: filename
      character(len=132) :: msg
      character(len=8) :: sta_author_default
!      logical :: radf_stname
        
      n = 0
      n_supp_stn = 0
      isstn = 99
      nstr1 = ' '
      sta_author = ' '
      sta_agency = ' '
      sta_deployment = ' '
      sta_loc = ' '
      sta_cha = ' '
   
      ! Read one or more supplemental station data files, for local (non-registered) stations.
      ! If there are any conflicts, the first instance encountered takes precedence (see stafind).
      ! The first column of the first line contains an index that specifies one of several standard formats.
      ! The rest of the first line can contain a comment about the source of the station data. All station
      ! codes and coordinates defined in supplemental station data files are written to the .stn file.
      if (suppstn) then
         do i = 1,n_supp_stn_file
            write (char1,'(i1)') i
            sta_author_default = 'SSTN-'//char1//'  ' ! Author based on index of the supplemental station file
            lun = lunit()
            call open_file (lun, suppfilnam(i), 'old')
            if (verbose_screen) write (*,'(/a)') trim(suppfilnam(i))
            write (io_stn_log,'(/a)') trim(suppfilnam(i))
            read (lun,'(a)') line
            if (verbose_screen) write (*,'(a)') trim(line)
            write (io_stn_log,'(a)') trim(line)
            read (line(1:1),'(i1)') isstn
            do k = 1,n_sta_max
               read (lun,'(a)',iostat=ios) line
               if (ios .lt. 0) exit
               if (line(1:1) .eq. '#') cycle ! Comment line
               n = n + 1
               if (n .gt. n_sta_max) call oops ('redsta: too many stations')
               select case (isstn)
            
                  case (0) ! New master file format from April 28, 2014
                     nstr1(n) = line(1:5)
                     read (line(7:15),'(f9.5)') glat
                     read (line(17:26),'(f10.5)') glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     read (line(28:32),'(i5)') station_elevation
                     read (line(34:37),'(i4)') dob
                     ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                     sta_author(n) = line(39:46)
                     if (read_ad .and. radf_stname(nstr1(n))) then
                        sta_agency(n) = line(48:52)
                        sta_deployment(n) = line(54:61)
                     end if
                     if (sta_author(n)(1:1) .eq. ' ') sta_author(n) = sta_author_default
                     read (line(66:72),'(i7)') jdate_on(n)
                     read (line(74:80),'(i7)') jdate_off(n)
                     sca0(n) = 0. ! No patch correction
                     sd1(n) = 1.0 ! Don't use patch correction variances.
                     if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,2i6,f8.3,1x,a8)') nstr1(n), glat, stalat(n),&
                      stalon(n), station_elevation, dob, ahgtr(n), sta_author(n)
   
                  case (1) ! ISC FFB format (geographic coordinates in deg-min-sec*10, elevation in m)
                     read (line(15:19),'(a5)') nstr1(n)
                     read (line(62:63),'(i2)') latdeg
                     read (line(64:65),'(i2)') latmin
                     read (line(66:68),'(i3)') latsec
                     glat = real(latdeg) + (real(latmin)/60.) + (real(latsec)/36000.)
                     if (line(69:69) .eq. 'S') glat = -glat
                     read (line(70:72),'(i3)') londeg
                     read (line(73:74),'(i2)') lonmin
                     read (line(75:77),'(i3)') lonsec
                     glon = real(londeg) + (real(lonmin)/60.) + (real(lonsec)/36000.)
                     if (line(78:78) .eq. 'W') glon = -glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     read (line(79:82),'(i4)') station_elevation
                     ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     sca0(n) = 0. ! No patch corrections.
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  
                  case (2) ! SEISAN station list, with 4-character station code starting in column 3
                  ! The standard format has been modified by adding optional fields for on- and off-date.
                  ! If column 2 is non-blank, reads columns 2:6 for support of 5-character codes.
                     if (line(2:2) .eq. ' ') then
                        read (line(3:6),'(a4)') nstr1(n)(1:4)
                        nstr1(n)(5:5) = ' '
                     else
                        read (line(2:6),'(a5)') nstr1(n)
                     end if
                     read (line(7:8),'(i2)') latdeg
                     read (line(9:13),'(f5.2)') xlatmin
                     glat = real(latdeg) + xlatmin/60.
                     if (line(14:14) .eq. 'S') glat = -glat
                     read (line(15:17),'(i3)') londeg
                     read (line(18:22),'(f5.2)') xlonmin
                     glon = real(londeg) + xlonmin/60.
                     if (line(23:23) .eq. 'W') glon = -glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     read (line(24:27),'(i4)') station_elevation
                     ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     if (len_trim(line) .eq. 40) read (line(34:40),'(i7)') jdate_on(n) ! On-date only
                     if (len_trim(line) .eq. 48) then
                        read (line(34:40),'(i7)') jdate_on(n)
                        read (line(42:48),'(i7)') jdate_off(n)
                     end if
                     sca0(n) = 0. ! No patch corrections.
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4,2i8)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n),&
                      jdate_on(n), jdate_off(n)
                  
                  case (3) ! Generic format, decimal degrees, geographic coordinates
                  ! Format changed on November 6, 2018 with v10.4.6
                  ! Required fields are:
                  ! Station code in columns 1:5
                  ! Decimal geographic latitude in columns 22:29 in f8.4 format
                  ! Decimal geographic longitude in columns 31:39 in f9.4 format
                  ! Integer elevation in meters in columns 41:45 in i5 format
                  ! The format also can carry agency, deployment, depth of burial, and operational dates as
                  ! optional fields.
                     if (line(11:11) .eq. '.') then
                        write (msg,'(a)') 'redsta: format error, generic format (type 3) changed on November 6, 2018'
                        call oops (trim(msg))
                     end if
                     nstr1(n) = line(1:5)
                     if (read_ad .and. radf_stname(nstr1(n))) then
                        sta_agency(n) = line(7:11)
                        sta_deployment(n) = line(13:20)
                     end if
                     read (line(22:29),'(f8.4)') glat
                     read (line(31:39),'(f9.4)') glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     read (line(41:45),'(i5)') station_elevation ! In meters
                     dob = 0
                     if (len_trim(line) .ge. 51) read (line(47:51),'(i5)') dob ! Depth of burial in meters (positive number)
                     ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     if (len_trim(line) .ge. 59) read (line(53:59),'(i7)') jdate_on(n)
                     if (len_trim(line) .ge. 67) read (line(61:67),'(i7)') jdate_off(n)
                     sca0(n) = 0. ! No patch corrections
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,2i6,f8.3,2i8)') nstr1(n), glat, stalat(n), stalon(n),&
                      station_elevation, dob, ahgtr(n), jdate_on(n), jdate_off(n)
                  
                  case (4) ! China Seismic Bureau format
                  ! Lower case seismic code in columns 1:3
                  ! Elevation in meters in columns 5:8 (i3)
                  ! Geographic latitude in degrees-minute-secs in columns 10:21
                  ! Longitude in degrees-minute-secs in columns 25:37
                     nstr1(n) = line(1:3)//'  '
                     read (line(5:8),'(i4)') station_elevation
                     ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                     read (line(10:11),'(i2)') latdeg
                     read (line(14:15),'(i2)') latmin
                     read (line(18:21),'(f4.1)') xlatsec
                     glat = real(latdeg) + (real(latmin)/60.) + xlatsec/3600.
                     read (line(25:27),'(i3)') londeg
                     read (line(30:31),'(i2)') lonmin
                     read (line(34:37),'(f4.1)') xlonsec
                     glon = real(londeg) + (real(lonmin)/60.) + xlonsec/3600.
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     sca0(n) = 0. ! No patch corrections
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  
                  case (5) ! NEIC format
                     nstr1(n) = line(4:8)
                     sta_loc(n) = line(10:11)
                     sta_cha(n) = line(13:15)
                     read (line(40:47),'(f8.4)') glat
                     read (line(48:57),'(f10.4)') glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     read (line(58:62),'(i5)') station_elevation ! In meters
                     dob = 0 ! Depth of burial in meters (positive number)
                     ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                     ! Operational epoch
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     if (line(83:84) .eq. 'to') then
                        read (line(66:75),'(i4,1x,i2,1x,i2)') year_oe, month_oe, day_oe
                        call juldat (year_oe, month_oe, day_oe, jdate_on(n), 0)
                        jdate_on(n) = year_oe*1000 + jdate_on(n)
                        read (line(86:95),'(i4,1x,i2,1x,i2)') year_oe, month_oe, day_oe
                        call juldat (year_oe, month_oe, day_oe, jdate_off(n), 0)
                        jdate_off(n) = year_oe*1000 + jdate_off(n)
                     end if
                     sca0(n) = 0. ! No patch corrections
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (read_ad .and. radf_stname(nstr1(n))) then
                        sta_agency(n) = 'FDSN '
                        sta_deployment(n) = line(1:2)//'      '
                     end if
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4,2i8)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n),&
                      jdate_on(n), jdate_off(n)
                  
                  case (6) ! MSU format
                  ! Seismic code in columns 1:5
                  ! Elevation in meters in columns 30:33 (i4)
                  ! Geographic latitude in degrees-minute-secs in columns 6:15
                  ! Longitude in degrees-minute-secs in columns 17:27
                     nstr1(n) = line(1:5)
                     read (line(30:33),'(i4)') station_elevation
                     ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                     read (line(6:7),'(i2)') latdeg
                     read (line(9:10),'(i2)') latmin
                     read (line(12:15),'(f4.1)') xlatsec
                     glat = real(latdeg) + (real(latmin)/60.) + xlatsec/36000.
                     if (line(16:16) .eq. 'S') glat = -glat
                     read (line(17:19),'(i3)') londeg
                     read (line(21:22),'(i2)') lonmin
                     read (line(24:27),'(f4.1)') xlonsec
                     glon = real(londeg) + (real(lonmin)/60.) + xlonsec/36000.
                     if (line(28:28) .eq. 'W') glon = -glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     sca0(n) = 0. ! No patch corrections
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                  
                  case (7) ! ISC-FDSN format
                     nstr1(n) = line(5:9)
                     sta_loc(n) = line(41:42)
                     sta_cha(n) = line(37:39)
                     read (line(12:19),'(f8.4)') glat
                     read (line(21:29),'(f9.4)') glon
                     call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                     stalat(n) = xlat
                     stalon(n) = xlon
                     read (line(31:35),'(i5)') station_elevation ! In meters
                     dob = 0 ! Depth of burial in meters (positive number)
                     ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
                     ! Operational epoch
                     if (line(68:68) .eq. 'Z') then
                        read (line(49:58),'(i4,1x,i2,1x,i2)') year_oe, month_oe, day_oe
                        call juldat (year_oe, month_oe, day_oe, jdate_on(n), 0)
                        jdate_on(n) = year_oe*1000 + jdate_on(n)
                     else
                        jdate_on(n) = 0
                     end if
                     if (line(94:94) .eq. 'Z') then
                        read (line(75:84),'(i4,1x,i2,1x,i2)') year_oe, month_oe, day_oe
                        call juldat (year_oe, month_oe, day_oe, jdate_off(n), 0)
                        jdate_off(n) = year_oe*1000 + jdate_off(n)
                     else
                        jdate_off(n) = 0
                     end if
                     sca0(n) = 0. ! No patch corrections
                     sd1(n) = 1.0 ! No patch correction variances
                     sta_author(n) = sta_author_default
                     if (read_ad .and. radf_stname(nstr1(n))) then
                        if (line(2:3) .eq. 'IR') then
                           sta_agency(n) = 'IR   '
                           sta_deployment(n) = line(2:3)//'      '
                        else ! FDSN
                           sta_agency(n) = 'FDSN '
                           sta_deployment(n) = line(2:3)//'      '
                        end if
                     end if
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4,2i8)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n),&
                      jdate_on(n), jdate_off(n)
                  
                  case (9) ! Same format as the former (up to v10.0.0) master station file
                     nstr1(n) = line(2:6)
                     read (line(8:16),'(f9.3)') stalat(n)
                     read (line(17:26),'(f10.3)') stalon(n)
                     read (line(27:36),'(f10.3)') ahgtr(n)
                     jdate_on(n) = 0
                     jdate_off(n) = 0
                     read (line(37:46),'(f10.3)') sca0(n)
                     sd1(n) = 1.0 ! Don't use patch correction variances.
                     call geogra (stalat(n), glat)
                     sta_author(n) = sta_author_default ! This format should only show up as a supplemental station file now
                     if (verbose_log) write (io_stn_log,'(a,1x,4f10.4)') nstr1(n), glat, stalat(n), stalon(n), ahgtr(n)
                     
                  case (99) ! No index number was supplied
                     call warnings ('redsta: no value supplied for isstn')
                     exit
                  
                  case default
                     write (msg,'(a,i3)') 'redsta: illegal value for isstn: ', isstn
                     call warnings (trim(msg))
                     exit
                  
               end select
            end do
            
            increment = n - n_supp_stn
            write (io_stn_log,'(i5,a)') increment, ' stations read'
            n_supp_stn = n ! Updated only after all entries are read from a supplemental file
            close (lun)
         end do
      end if
   
      ! Read the main station data file.
      ! From April 28, 2014, the file has a completely new format
      ! The former format is still supported as isstn = 9
      filename = trim(station_path)//dirsym//trim(station_master)
      if (verbose_screen) write (*,'(/a)') trim(filename)
      write (io_stn_log,'(/a)') trim(filename)
      lun = lunit()
      call open_file (lun, filename, 'old')
      read (lun,'(a)') line
      if (verbose_screen) write (*,'(a)') trim(line)
      write (io_stn_log,'(a)') trim(line)
      read (line(1:1),'(i1)') isstn
      select case (isstn)
      
         case (0) ! New master file format from April 28, 2014
            do k = 1,n_sta_max
               read (lun,'(a)',iostat=ios) line
               if (ios .lt. 0) exit
               if (line(1:1) .eq. '#') cycle ! Comment line
               n = n + 1
               if (n .gt. n_sta_max) call oops ('redsta: too many stations')
               nstr1(n) = line(1:5)
               read (line(7:15),'(f9.5)') glat
               read (line(17:26),'(f10.5)') glon
               call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
               stalat(n) = xlat
               stalon(n) = xlon
               read (line(28:32),'(i5)') station_elevation
               read (line(34:37),'(i4)') dob
               ahgtr(n) = real(station_elevation - dob)*1.e-3 ! Convert to km
               sta_author(n) = line(39:46)
               if (read_ad .and. radf_stname(nstr1(n))) then
                  sta_agency(n) = line(48:52)
                  sta_deployment(n) = line(54:61)
               end if
               read (line(66:72),'(i7)') jdate_on(n)
               read (line(74:80),'(i7)') jdate_off(n)
               sca0(n) = 0. ! No patch correction
               sd1(n) = 1.0 ! Don't use patch correction variances.
               if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,2i6,f8.3,1x,a8)') nstr1(n), glat, stalat(n), stalon(n),&
                station_elevation, dob, ahgtr(n), sta_author(n)
            end do
            
         case (9) ! Former master station file format, decimal geographic coordinates, elevation in m, operational epochs
            do k = 1,n_sta_max
               read (lun,'(a)',iostat=ios) line
               if (ios .lt. 0) exit
               if (line(2:2) .ne. ' ') then
                  n = n + 1
                  if (n .gt. n_sta_max) call oops ('redsta: too many stations')
                  nstr1(n) = line(2:6)
                  read (line(60:67),'(f8.4)') glat
                  read (line(69:77),'(f9.4)') glon
                  call geocen (glat, glon, xlat, sxlat, cxlat, xlon, sxlon, cxlon) ! Convert to geocentric coordinates
                  stalat(n) = xlat
                  stalon(n) = xlon
                  read (line(80:84),'(i5)') station_elevation
                  ahgtr(n) = real(station_elevation)*1.e-3 ! Convert to km
                  read (line(90:96),'(i7)') jdate_on(n)
                  read (line(99:105),'(i7)') jdate_off(n)
                  sca0(n) = 0. ! No patch correction
                  sd1(n) = 1.0 ! Don't use patch correction variances.
                  sta_author(n) = sta_author_default
                  if (verbose_log) write (io_stn_log,'(a,1x,3f10.4,i6,f8.3)') nstr1(n), glat, stalat(n), stalon(n),&
                   station_elevation, ahgtr(n)
               end if
            end do
            
         case default
            write (msg,'(a,i1)') 'redsta: illegal value for isstn for master station file: ', isstn
            call oops (trim(msg))
            
      end select
      
      nstat1 = n ! Number of stations read
      
      n_master = nstat1 - n_supp_stn
      write (io_stn_log,'(i5,a)') n_master, ' stations read'
      if (suppstn) then
         write (io_stn_log,'(/i5,a)') nstat1, ' total stations read'
         if (verbose_screen) then
            write (msg,'(a,i6)') 'redsta: total number of stations        ', nstat1
            call fyi (trim(msg))
         end if
      end if
            
      ! Check for duplicate station codes and conflicts
      call station_duplicate ()
         
      if (suppstn) then
         if (verbose_screen) then
            write (msg,'(a,i6)') 'redsta: duplicate stations              ', n_dupe
            call fyi (trim(msg))
            if (n_dupe .gt. 0) then
               write (msg,'(a,i6)') 'redsta: ...pure duplicates              ', n_dupe_pure
               call fyi (trim(msg))
               write (msg,'(a,i6)') 'redsta: ...with minor differences       ', n_dupe_minor
               call fyi (trim(msg))
               write (msg,'(a,i6)') 'redsta: ...with significant differences ', n_dupe_significant
               call fyi (trim(msg))
               write (msg,'(a,i6)') 'redsta: ...station conflicts            ', n_dupe_conflict
               call fyi (trim(msg))
            end if
         end if
         write (io_stn_log,'(/a)') 'Duplicate stations:'
         if (n_dupe .gt. 0) then
            write (io_stn_log,'(i5,a)') n_dupe_pure, ' pure duplicates'
            write (io_stn_log,'(i5,a)') n_dupe_minor, ' with minor differences'
            write (io_stn_log,'(i5,a)') n_dupe_significant, ' with significant differences'
            write (io_stn_log,'(i5,a)') n_dupe_conflict, ' station conflicts'
         end if
         write (io_stn_log,'(i5,a)') n_dupe, ' total'
      end if
      
      close (lun)
         
      return
      
   end subroutine redsta
         
         
!*****************************************************************************************
   subroutine station_duplicate ()
   
   ! Check for duplicate station codes, looping through all supplemental
   ! stations and checking against the entire list of stations. When an
   ! NEIC-formatted supplemental station file is used, typically following
   ! use of the "NSMD" command, or a supplemental station file based on a
   ! search of the ISC-FDSN station file (command IFSM) many duplicates may
   ! be noted. This happens
   ! because the selection of entries for the supplemental file is based on
   ! matching deployment, station code, location and channel fields, as well
   ! as operational date as eash phase reading is read from the MNF data
   ! files. When station coordinates are assigned for a phase reading,
   ! however, location and channel fields are not considered, so supplemental
   ! station entries that differ only by those fields are considered
   ! duplicates. These are noted by an FYI message during runtime. These
   ! warnings can usually be ignored, because differences in location and
   ! channel codes rarely if ever lead to different station coordinates.
   ! The ~.stn output file contains detailed information that can be used to
   ! investigate, and if necessary, correct instances of duplication or
   ! incorrectly assigned station coordinates.
   
      integer :: j
      integer :: k
      integer :: n_dupe_j
      real :: glatj
      real :: glatk
      real :: diff_lat
      real :: diff_lon
      real :: diff_hgt
      real :: lonj
      real :: lonk
      character(len=72) :: form1
      character(len=72) :: form2
      character(len=132) :: msg
      character(len=40) :: result
         
      form1 = '(a,1x,a,1x,a,1x,f10.4,10x,f10.4,f10.3,1x,a)'
      form2 = '(21x,f10.4,10x,f10.4,f10.3,a/)'
      result = ' '
      n_dupe = 0 ! Entries with the same station code
      n_dupe_pure = 0 ! Identical coordinates
      n_dupe_minor = 0 ! Minor difference in coordinates, but probably the same station
      n_dupe_significant = 0 ! Significant difference in coordinates, but probably the same station
      n_dupe_conflict = 0 ! Different stations with the same code
      
      if (.not.suppstn) return ! No check unless there are supplemental station files
      
      do j = 1,n_supp_stn ! loop through stations from all supplemental station files
         call geogra (stalat(j), glatj)
         lonj = stalon(j)
         call set_longitude_range (lonj, 0)
         n_dupe_j = 0 ! number of duplicates for this station code
         do k = j+1,nstat1
            if (nstr1(k) .eq. nstr1(j)) then
               n_dupe = n_dupe + 1
               n_dupe_j = n_dupe_j + 1
               call geogra (stalat(k), glatk)
               diff_lat = glatj - glatk
               lonk = stalon(k)
               call set_longitude_range (lonk, 0)
               diff_lon = lonj - lonk
               diff_hgt = ahgtr(j) - ahgtr(k)
               ! Region of equivalnce is about 100 m
               if (abs(diff_lat) .gt. 1.1 .or. abs(diff_lon) .gt. 1.1) then
                  result = ' conflicting station codes'
                  n_dupe_conflict = n_dupe_conflict + 1
               else if (abs(diff_lat) .gt. 1.0e-3 .or. abs(diff_lon) .gt. 1.0e-3 .or. abs(diff_hgt) .gt. 1.0e-1) then
                  result = ' significant difference in coordinates'
                  n_dupe_significant = n_dupe_significant + 1
               else if (abs(diff_lat) .gt. 1.0e-4 .or. abs(diff_lon) .gt. 1.0e-4 .or. abs(diff_hgt) .gt. 1.0e-2) then
                  result = ' minor difference in coordinates'
                  n_dupe_minor = n_dupe_minor + 1
               else ! Same coordinates
                  result = ' pure duplicate'
                  n_dupe_pure = n_dupe_pure + 1
               end if
               write (io_stn_log,form1) nstr1(j), sta_agency(j), sta_deployment(j), glatj,&
                lonj, ahgtr(j), sta_author(j) ! First instance
               write (io_stn_log,form1) nstr1(k), sta_agency(k), sta_deployment(k), glatk,&
                lonk, ahgtr(k), sta_author(k) ! Second instance
               write (io_stn_log,form2) diff_lat, diff_lon, diff_hgt, result ! Difference
               if (n_dupe_j .eq. 1) then
                  duplication(j) = result
               else
                  write (msg,'(5a)') 'station_duplicate: multiple duplicates for ', nstr1(j), '(', trim(result), ')'
                  call warnings (trim(msg))
               end if
            end if
         end do
      end do
   
      return
      
   end subroutine station_duplicate
   
   
!*****************************************************************************************
   subroutine station_check (iev, ird)
   
   ! Check for a few known problems with certain stations, mainly derived from Bob Engdahl's
   ! experiences with the EHB and GEM processing. There are two primary situations:
   ! 1) Certain stations are known to have timing problems during certain periods, and readings
   !    in these periods are flagged (fcode = 't').
   ! 2) Certain stations were moved significant distances without changing the station code 
   !    and a new code for the new location has been assigned retroactively.
   
      integer, intent(in) :: iev ! event number
      integer, intent(in) :: ird ! reading index
      character(len=132) :: msg
      
      ! These checks are only made for stations with deployment code IR or blank.
      if (deployment(iev,ird) .ne. '        ' .and. deployment(iev,ird) .ne. 'IR      ') return 
         
      if (stname(iev,ird) .eq. 'DAVO ') then
         if (ipay(iev,ird) .lt. 1999) then
            stname(iev,ird) = 'DAVOS'
         else
            stname(iev,ird) = 'DAVON'
         end if
         sad(iev,ird)(1:5) = stname(iev,ird)
         call fyi ('station_check: station code DAVO changed to '//stname(iev,ird))
      else if (stname(iev,ird) .eq. 'FCC  ' .and. (ipay(iev,ird) .ge. 1992 .and. ipay(iev,ird) .le. 1993)) then
         fcode(iev,ird) = 't'
         write (msg,'(a,i4)') 'station_check: station FCC had timing problems in ', ipay(iev,ird)
         call fyi (trim(msg))
      else if (stname(iev,ird) .eq. 'FITZ*') then
         stname(iev,ird) = 'FITZ '
         sad(iev,ird)(1:5) = stname(iev,ird)
         call fyi ('station_check: station code FITZ* changed to FITZ')
      else if (stname(iev,ird) .eq. 'MAG  ' .and. ipay(iev,ird) .gt. 1970) then
         stname(iev,ird) = 'MGD  '
         sad(iev,ird)(1:5) = stname(iev,ird)
         call fyi ('station_check: station MAG moved without changing codes; subsequent readings changed to MGD')
      else if (stname(iev,ird) .eq. 'NRI  ' .and. ipay(iev,ird) .gt. 1993) then
         stname(iev,ird) = 'NRIS '
         sad(iev,ird)(1:5) = stname(iev,ird)
         call fyi ('station_check: station NRI moved without changing codes; subsequent readings changed to NRIS')
      else if (stname(iev,ird) .eq. 'TUC  ' .and. (ipay(iev,ird) .ge. 1992 .and. ipay(iev,ird) .le. 1995)) then
         fcode(iev,ird) = 't'
         write (msg,'(a,i4)') 'station_check: station TUC had timing problems in ', ipay(iev,ird)
         call fyi (trim(msg))
      else if (stname(iev,ird) .eq. 'WIT  ' .and. (ipay(iev,ird) .ge. 1964 .and. ipay(iev,ird) .le. 1979)) then
         fcode(iev,ird) = 't'
         write (msg,'(a,i4)') 'station_check: station WIT had timing problems in ', ipay(iev,ird)
         call fyi (trim(msg))
      else if (stname(iev,ird) .eq. 'YKA  ' .and. (ipay(iev,ird) .ge. 1991 .and. ipay(iev,ird) .le. 1995)) then
         fcode(iev,ird) = 't'
         write (msg,'(a,i4)') 'station_check: station YKA had timing problems in ', ipay(iev,ird)
         call fyi (trim(msg))
      end if
   
      return
      
   end subroutine station_check
   
   
!*****************************************************************************************
   subroutine supp_station_file_write (io, station, agency, deployment, lat, lon, elev, dob, ion, ioff, comment)
   
   ! Write station information in the format of the generic supplemental station file (isstn = 3)
   
   ! io = logical unit
   ! station = station code (5 characters max)
   ! agency = ADSLC agency code
   ! deployment = ADSLC deployment code
   ! lat = geographic latitude
   ! lon = geographic longitude
   ! elev = station elevation in meters
   ! dob = depth of burial in meters (positive down)
   ! ion = date on (year in columns 1:4 + day of year in columns 5:7)
   ! ioff = date off (year in columns 1:4 + day of year in columns 5:7)
   ! comment = optional text written at the end of the line
   
   ! dob, ion and ioff are printed as blanks if they are zero
   
      integer, intent(in) :: dob
      integer, intent(in) :: elev
      integer, intent(in) :: io
      integer, intent(in) :: ioff
      integer, intent(in) :: ion
      real, intent(in) :: lat
      real, intent(in) :: lon
      character(len=5), intent(in) :: agency
      character(len=*), intent(in) :: comment
      character(len=8), intent(in) :: deployment
      character(len=5) :: dob_pr
      character(len=80) :: form
      character(len=7) :: ion_pr
      character(len=7) :: ioff_pr
      character(len=5), intent(in) :: station
      
      if (dob .eq. 0) then
         dob_pr = ' '
      else if (dob .lt. 0) then
         call warnings ('supp_station_file_write: depth of burial must be positive')
         write (dob_pr,'(i7)') abs(dob)
      else
         write (dob_pr,'(i7)') dob
      end if
      
      if (ion .eq. 0) then
         ion_pr = ' '
      else
         write (ion_pr,'(i7)') ion
      end if
      
      if (ioff .eq. 0) then
         ioff_pr = ' '
      else
         write (ioff_pr,'(i7)') ioff
      end if
         
      form = '(a5,1x,a5,1x,a8,1x,f8.4,1x,f9.4,1x,i5,1x,a5,1x,a7,1x,a7,1x,a)'
      write (io,form) station, agency, deployment, lat, lon, elev, dob_pr, ion_pr, ioff_pr, trim(comment)
      
      return
      
   end subroutine supp_station_file_write
   
   
!*****************************************************************************************
   logical function radf_stname (stname_in)
   
      integer :: i
      character(len=5), intent(in) :: stname_in
      
      radf_stname = .false.
      do i = 1,n_radf
         if (stname_in .eq. radf_stn(i)) then
            radf_stname = .true.
            exit
         end if
      end do
      
      return
   
   end function radf_stname
   

!*****************************************************************************************
   subroutine read_ifsm (io_ifsm, filename, outfile)
      
   ! Read the main ISC-FDSN station metadata file

      integer :: day_oe1 ! operational epoch day on
      integer :: day_oe2 ! operational epoch day off
      integer :: i
      integer, intent(in) :: io_ifsm ! logical unit number
      integer :: ios ! I/O status
      integer :: jdate ! Julian date
      integer :: month_oe1 ! operational epoch month on
      integer :: month_oe2 ! operational epoch month off
      integer :: n_entries ! number of entries in the main ISC-FDSN station metadata file
      integer :: year_oe1 ! operational epoch year on
      integer :: year_oe2 ! operational epoch year off
      character(len=*), intent(in) :: filename ! main ISC-FDSN station metadata filename
      character(len=94) :: line94
      character(len=132) :: msg
      character(len=*) :: outfile ! run-specific output filename

      !Initialize
      n_ifsm_written = 0
      
      call open_file (io_ifsm, filename, 'old')
      
      ! Count the entries
      ! The first line in the file is a creation date but it's not a problem to allocate
      ! an extra entry.
      n_entries = 0
      do
         read (io_ifsm,'(a)',iostat=ios) line94
         if (ios .lt. 0) exit
         n_entries = n_entries + 1
      end do
      rewind (io_ifsm)
      
      ! Allocate variable arrays
      call ifsm_allocate (n_entries)
      
      read (io_ifsm,'(a)') line94 ! First line is the file creation date
      do i = 1, n_entries-1
         read (io_ifsm,'(a)',iostat=ios) line94
         ifsm_line(i) = line94
         ds(i) = line94(2:3)//line94(5:9) ! Deployment, Station Code
         if (line94(49:49) .ne. ' ') then
            read (line94(49:58),'(i4,1x,i2,1x,i2)') year_oe1, month_oe1, day_oe1
         else
            year_oe1 = 1900
            month_oe1 = 1
            day_oe1 = 1
         end if
         if (line94(75:75) .ne. ' ') then
            read (line94(75:84),'(i4,1x,i2,1x,i2)') year_oe2, month_oe2, day_oe2
         else
            year_oe2 = 2099
            month_oe2 = 12
            day_oe2 = 31
         end if
         call juldat (year_oe1, month_oe1, day_oe1, jdate, 0)
         ifsm_on(i) = year_oe1*1000 + jdate
         call juldat (year_oe2, month_oe2, day_oe2, jdate, 0)
         ifsm_off(i) = year_oe2*1000 + jdate
      end do
      n_ifsm_read = i
      write (msg,'(a,i8,a)') 'read_ifsm: ', n_ifsm_read, ' ISC-FDSN station metadata entries read: '
      call fyi (trim(msg))
      close (io_ifsm)
      
      ! Open a new file for the selected station entries
      ifsm_stn_file = trim(outfile)//'_ifsm_stn.dat'
      call open_file (io_ifsm, ifsm_stn_file, 'new')
      write (io_ifsm,'(a)') '3 Supplemental station file from ISC-FDSN Station Metadata'
      
      return

   end subroutine read_ifsm
   

!*****************************************************************************************
   subroutine write_ifsm (io_ifsm, ds_test, jdate)

   ! Generate a supplemental station file entry from the ISC-FDSN station metadata file
   ! for stations missing from the current station list (supplemental lists + master list).
   ! The output file is in the "generic" (type 3) format.
   ! A match, which causes a new entry in the supplemental station file, is based on
   ! "deployment/station code" from the MNF data file, and the date of the event must
   ! be within the operational date range of the entry in the ISC-FDSN station metadata
   ! file. For entries that carry no operational date, a default (1900-2099) is used.

      ! Arguments
      integer, intent(in) :: io_ifsm
      integer, intent(in) :: jdate
      character(len=7), intent(in) :: ds_test

      integer :: i
      integer :: elev_ifsm
      real :: lat_ifsm, lon_ifsm
      character(len=5) :: agency_ifsm
      character(len=132) :: msg
      logical :: matched

      matched = .false.
      
      ! First check against entries that have already been matched in earlier events
      if (n_ifsm_written .ge. 1) then
         do i = 1,n_ifsm_written
            if (ds_test .ne. ds_matched(i)) cycle ! Network/deployment code
            if (date_range(jdate,ifsm_on_matched(i),ifsm_off_matched(i))) then
               matched = .true.
               exit
            end if
         end do
      end if
      
      ! Test against the full ISC-FDSN data file
      if (.not.matched) then
         if (n_ifsm_written .lt. n_ifsm_max) then
            do i = 1, n_ifsm_read
               if (ds_test .eq. ds(i)) then
                  if (date_range(jdate,ifsm_on(i),ifsm_off(i))) then
                     n_ifsm_written = n_ifsm_written + 1
                     read (ifsm_line(i)(12:19),'(f8.4)') lat_ifsm
                     read (ifsm_line(i)(21:29),'(f9.4)') lon_ifsm
                     read (ifsm_line(i)(31:35),'(i5)') elev_ifsm
                     if (ds(i)(1:2) .eq. 'IR') then
                        agency_ifsm = 'ISC  '
                     else
                        agency_ifsm = 'FDSN '
                     end if
                     call supp_station_file_write (io_ifsm, ds(i)(3:7), agency_ifsm,&
                      ds(i)(1:2)//'      ', lat_ifsm, lon_ifsm, elev_ifsm, 0,&
                      ifsm_on(i), ifsm_off(i), ' ') 
                     ds_matched(n_ifsm_written) = ds(i)
                     ifsm_on_matched(n_ifsm_written) = ifsm_on(i)
                     ifsm_off_matched(n_ifsm_written) = ifsm_off(i)
                     exit
                  end if
               end if
            end do
         else if (n_ifsm_written .eq. n_ifsm_max) then
            write (msg,'(a,i4,a)') 'read_mnf_134: maximum number of matched ISC-FDSN metadata entries (',&
             n_ifsm_max, ') written'
            call warnings (trim(msg))
         end if
      end if
      
      return
            
   end subroutine write_ifsm
   
   
!*****************************************************************************************
!> Read the main NEIC station metadata file
   
   subroutine read_nsmd (io_nsmd, filename, outfile)
   
      integer :: day_oe1 ! operational epoch day on
      integer :: day_oe2 ! operational epoch day off
      integer :: i
      integer, intent(in) :: io_nsmd ! logical unit number
      integer :: ios ! I/O status
      integer :: jdate ! Julian date
      integer :: month_oe1 ! operational epoch month on
      integer :: month_oe2 ! operational epoch month off
      integer :: n_entries ! number of entries in the main NEIC station metadata file
      integer :: year_oe1 ! operational epoch year on
      integer :: year_oe2 ! operational epoch year off
      character(len=*), intent(in) :: filename ! main NEIC station metadata filename
      character(len=102) :: line102
      character(len=132) :: msg
      character(len=*) :: outfile ! run-specific output filename
      
      !Initialize
      n_nsmd_written = 0
      
      call open_file (io_nsmd, filename, 'old')
      
      ! Count the entries
      n_entries = 0
      do
         read (io_nsmd,'(a)',iostat=ios) line102
         if (ios .lt. 0) exit
         n_entries = n_entries + 1
      end do
      rewind (io_nsmd)
      
      ! Allocate arrays
      call nsmd_allocate (n_entries)
      
      do i = 1,n_entries
         read (io_nsmd,'(a)') line102
         nsmd_line(i) = line102
         dslc(i) = line102(1:2)//line102(4:8)//line102(10:11)//line102(13:15) ! Deployment, Station Code, Location, Channel
         if (line102(63:63) .eq. '.') then
            read (line102(66:75),'(i4,1x,i2,1x,i2)') year_oe1, month_oe1, day_oe1
            read (line102(86:95),'(i4,1x,i2,1x,i2)') year_oe2, month_oe2, day_oe2
         else if (line102(64:64) .eq. '.') then
            read (line102(67:76),'(i4,1x,i2,1x,i2)') year_oe1, month_oe1, day_oe1
            read (line102(87:96),'(i4,1x,i2,1x,i2)') year_oe2, month_oe2, day_oe2
         else
            write (msg,'(3a,i8)') 'read_nsmd: Error reading ', trim(filename), ' line ', i
            call oops (trim(msg))
         end if
         call juldat (year_oe1, month_oe1, day_oe1, jdate, 0)
         nsmd_on(i) = year_oe1*1000 + jdate
         call juldat (year_oe2, month_oe2, day_oe2, jdate, 0)
         nsmd_off(i) = year_oe2*1000 + jdate
      end do
      n_nsmd_read = n_entries
      write (msg,'(a,i8,a)') 'read_nsmd: ', n_nsmd_read, ' NEIC station metadata entries read: '
      call fyi (trim(msg))
      close (io_nsmd)
      
      ! Open a new run-specific file for the matched station entries and write the header line.
      ! Matched entries will not be written until the event data files are read.
      
      nsmd_stn_file = trim(outfile)//'_nsmd_stn.dat'
      call open_file (io_nsmd, nsmd_stn_file, 'new')
      write (io_nsmd,'(a)') '5 Supplemental station file from NEIC Station Metadata'
      
      return

   end subroutine read_nsmd
   
   
!*****************************************************************************************
!> Check against NEIC station metadata file as phase lines as they're read, write matches
!  to a supplemental station file.

   subroutine write_nsmd (io_nsmd, dslc_test, jdate)
   
      integer :: i
      integer, intent(in) :: io_nsmd
      integer, intent(in) :: jdate
      character(len=12), intent(in) :: dslc_test
      character(len=132) :: msg
      logical :: matched

      matched = .false.
      
      ! First check against entries that have already been matched
      if (n_nsmd_written .ge. 1) then
         do i = 1,n_nsmd_written
            if (dslc_test(1:7) .ne. dslc_matched(i)(1:7)) cycle ! Network/deployment code
            if (date_range(jdate,nsmd_on_matched(i),nsmd_off_matched(i))) then
               matched = .true.
               exit
            end if
         end do
      end if
      
      ! Test against the full data file
      if (.not.matched) then
         if (n_nsmd_written .lt. n_nsmd_max) then
            do i = 1, n_nsmd_read
               if (dslc_test(1:7) .eq. dslc(i)(1:7)) then
                  if (date_range(jdate,nsmd_on(i),nsmd_off(i))) then
                     n_nsmd_written = n_nsmd_written + 1
                     write (io_nsmd,'(a)') nsmd_line(i)
                     dslc_matched(n_nsmd_written)(1:7) = dslc(i)(1:7)
                     nsmd_on_matched(n_nsmd_written) = nsmd_on(i)
                     nsmd_off_matched(n_nsmd_written) = nsmd_off(i)
                     exit
                  end if
               end if
            end do
         else if (n_nsmd_written .eq. n_nsmd_max) then
            write (msg,'(a,i4,a)') 'read_mnf_134: maximum number of matched NEIC metadata entries (',&
             n_nsmd_max, ') written'
            call warnings (trim(msg))
         end if
      end if
      
      return
   
   end subroutine write_nsmd
   

!*****************************************************************************************
   subroutine ifsm_allocate (n)
      
      integer :: error
      integer, intent(in) :: n
      character(len=40) :: p
      
      p = 'ifsm_allocate'
      
      allocate (ifsm_off(n), stat=error); if (error .gt. 0) call allocation_error (p, 'ifsm_off', error)
      allocate (ifsm_on(n), stat=error); if (error .gt. 0) call allocation_error (p, 'ifsm_on', error)
      allocate (ds(n), stat=error); if (error .gt. 0) call allocation_error (p, 'ds', error)
      allocate (ifsm_line(n), stat=error); if (error .gt. 0) call allocation_error (p, 'ifsm_line', error)
      
      return
      
   end subroutine ifsm_allocate
   
      
!*****************************************************************************************
   subroutine ifsm_deallocate ()
      
      integer :: error
      character(len=40) :: p
      
      p = 'ifsm_deallocate'
      
      deallocate (ifsm_off, stat=error); if (error .gt. 0) call deallocation_error (p, 'ifsm_off', error)
      deallocate (ifsm_on, stat=error); if (error .gt. 0) call deallocation_error (p, 'ifsm_on', error)
      deallocate (ds, stat=error); if (error .gt. 0) call deallocation_error (p, 'ds', error)
      deallocate (ifsm_line, stat=error); if (error .gt. 0) call deallocation_error (p, 'ifsm_line', error)
      
      return
      
   end subroutine ifsm_deallocate
   
      
!*****************************************************************************************
   subroutine nsmd_allocate (n)
      
      integer :: error
      integer :: n
      character(len=40) :: p
      
      p = 'nsmd_allocate'
      
      allocate (nsmd_off(n), stat=error); if (error .gt. 0) call allocation_error (p, 'nsmd_off', error)
      allocate (nsmd_on(n), stat=error); if (error .gt. 0) call allocation_error (p, 'nsmd_on', error)
      allocate (dslc(n), stat=error); if (error .gt. 0) call allocation_error (p, 'dslc', error)
      allocate (nsmd_line(n), stat=error); if (error .gt. 0) call allocation_error (p, 'nsmd_line', error)
      
      return
      
   end subroutine nsmd_allocate
   
      
!*****************************************************************************************
   subroutine nsmd_deallocate ()
      
      integer :: error
      character(len=40) :: p
      
      p = 'nsmd_deallocate'
      
      deallocate (nsmd_off, stat=error); if (error .gt. 0) call deallocation_error (p, 'nsmd_off', error)
      deallocate (nsmd_on, stat=error); if (error .gt. 0) call deallocation_error (p, 'nsmd_on', error)
      deallocate (dslc, stat=error); if (error .gt. 0) call deallocation_error (p, 'dslc', error)
      deallocate (nsmd_line, stat=error); if (error .gt. 0) call deallocation_error (p, 'nsmd_line', error)
      
      return
      
   end subroutine nsmd_deallocate
   
      
!*****************************************************************************************
end module mloc_stations

