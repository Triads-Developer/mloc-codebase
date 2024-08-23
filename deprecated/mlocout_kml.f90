!> Procedures realted to constructing KML files

module mlocout_kml

   use declare_calibration
   use declare_cluster_vectors
   use declare_configuration
   use declare_environment
   use declare_events
   use declare_hypocentroid
   use declare_lun
   use declare_output
   use declare_phases
   use mloclib_geog

   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine write_kml (it)
   
   ! Creates a kml file of the epicenters that can be displayed in Google Earth
   ! Expects to find icons in a folder named "tables/kml" under the executable
   ! mloc file. Icons are color-coded by depth for events with depth constraint:
   !   0 -  9 km - red
   !  10 - 19 km - green
   !  20 - 29 km - skyblue
   !  30+     km - blue
   ! A yellow icon is used for events set to the cluster default depth
   
      integer :: iev
      integer :: irange
      integer :: it
      integer :: it1
      real :: azeqst
      real :: azesdg
      real :: azsedg
      real :: azsteq
      real :: delta
      real :: deltdg
      real :: deltkm
      real :: xdep
      real :: xlat
      real :: xlath
      real :: xlon
      real :: xlonh
      real :: xlonmax
      real :: xlonmin
      real :: xmagms
      real :: xmagmw
      real :: ylatmax
      real :: ylatmin
      character(len=24) :: cdate
      character(len=20) :: color
      character(len=30) :: dep_con
      character(len=132) :: msg
      character(len=100) :: outfil
      
      it1 = it + 1
      color = ' '
            
      outfil = trim(outfile)//'.kml'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_kml: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      call open_file (io_out, outfil, 'new')
      
      ! Focus point for KML is the iterated hypocentroid. Don't worry about correction in the case of indirect calibration
      xlath = lath(it1)
      xlonh = lonh(it1)
      call set_longitude_range (xlonh, 0) ! Google Earth requires longitude -180 < xlon < 180
      
      ! Calculate range for KML "LookAt" command from the diagonal length of the cluster bounding box
      call map_boundaries (it1, 0, xlonmin, xlonmax, ylatmin, ylatmax)
      call delaz (ylatmin, xlonmin, ylatmax, xlonmax, delta, deltdg, deltkm, azeqst, azesdg, azsteq, azsedg, 0)
      irange = nint(deltkm*1.e3)
            
      call kml_prelude (io_out, xlath, xlonh, irange)
      
      do iev = 1,nev ! Loop over events
   
         write (cdate,'(i4,a,i2,a,i2,1x,i2,a,i2,a,f6.3)') iyre(iev), '/', mone(iev), '/', idye(iev),&
          hourp(iev,it1), ':', minp(iev,it1), ':', secp(iev,it1)
   
         if (indirect_cal) then
            xlat = latp_cal(iev)
            xlon = lonp_cal(iev)
            xdep = depthp_cal(iev)
         else
            xlat = latp(iev,it1)
            xlon = lonp(iev,it1)
            xdep = depthp(iev,it1)
         end if
         
         call set_longitude_range (xlon, 0) ! Google Earth requires longitude -180 < xlon < 180
         
         if (depset_pr(iev) .eq. 'c') then
            color = 'yellow'
            dep_con = ' cluster default'
         else
            dep_con = ' constrained ('//depset_pr(iev)//')'
            if (xdep .lt. 10.0) then
               color = 'red'
            else if (xdep .lt. 20.0) then
               color = 'green'
            else if (xdep .lt. 30.0) then
               color = 'skyblue'
            else
               color = 'blue'
            end if
         end if
    
         xmagms = 0.
         xmagmw = 0.
   
         call kml_point (io_out, iev, cdate, xlat, xlon, xdep, rmag(iev), xmagms, xmagmw, color, dep_con)
      
      end do
   
      call kml_coda (io_out)
      
      close (io_out)
   
      return
   
   end subroutine write_kml
 

!*****************************************************************************************
   subroutine kml_point (iout, iev, cdate, xlat, xlon, xdep, xmagmb, xmagms, xmagmw, color, dep_con)
    
      integer :: iev
      integer :: iout
      real :: scale
      real :: scale_factor
      real :: scale0
      real :: xdep
      real :: xlat
      real :: xlon
      real :: xmag
      real :: xmagmb
      real :: xmagms
      real :: xmagmw
      character(len=30) :: dep_con
      character(len=24) :: cdate
      character(len=21) :: icon
      character(len=20) :: color
      character(len=8) :: cxlat
      character(len=8) :: cxlon
      character(len=5) :: cxdep
      character(len=6) :: cscale
      character(len=3) :: cxiev
      character(len=3) :: cxmagmb
      character(len=3) :: cxmagms
      character(len=3) :: cxmagmw
      
      scale0 = 0.3
      scale_factor = 1.75
   
      icon = '#'//color
      
      xmag = amax1(xmagmb, xmagms, xmagmw)
      if (xmag .gt. 2.) then
         scale = scale0 + ((xmag - 2.0)/5.0)*scale_factor
      else
         scale = scale0
      end if
      write (cscale,'(f6.3)') scale
      
      write (cxiev,'(i3)') iev
      write (cxlat,'(f8.3)') xlat
      write (cxlon,'(f8.3)') xlon
      write (cxdep,'(f5.1)') xdep
      
      write (iout,'(a)')  '   <Placemark>'
      write (iout,'(3a)') '      <name>', trim(adjustl(cxiev)), '</name>'
      write (iout,'(3a)') '      <description><![CDATA[Event : <b>', trim(adjustl(cxiev)), '</b>'
      write (iout,'(3a)')   '      <br>Time : <b>', cdate, '</b>'
      if (xmagmb .lt. 0.1 .and. xmagms .lt. 0.1 .and. xmagmw .lt. 0.1) then
         write (iout,'(a)') '      <br>Magnitude : <b>Unknown</b>'
      else
         if (xmagmb .gt. 0.) then
            write (cxmagmb,'(f3.1)') xmagmb
            write (iout,'(4a)') '      <br>Magnitude : <b>', trim(cxmagmb), ' mb', '</b>'
         end if
         if (xmagms .gt. 0.) then
            write (cxmagms,'(f3.1)') xmagms
            write (iout,'(4a)') '      <br>Magnitude : <b>', trim(cxmagmb), ' Ms', '</b>'
         end if
         if (xmagmw .gt. 0.) then
            write (cxmagmw,'(f3.1)') xmagmw
            write (iout,'(4a)') '      <br>Magnitude : <b>', trim(cxmagmb), ' Mw', '</b>'
         end if
      end if
      write (iout,'(3a)') '      <br>Latitude : <b>', trim(adjustl(cxlat)), '</b>'     
      write (iout,'(3a)') '      <br>Longitude : <b>', trim(adjustl(cxlon)), '</b>'     
      write (iout,'(4a)') '      <br>Depth : <b>', trim(adjustl(cxdep)), trim(dep_con), '</b>]]></description>'     
      write (iout,'(3a)') '      <styleUrl>', trim(icon), '</styleUrl>'
      write (iout,'(a)')  '      <Style>'
      write (iout,'(a)')  '         <IconStyle>'
      write (iout,'(3a)') '            <scale>', trim(adjustl(cscale)), '</scale>'
      write (iout,'(a)')  '            <hotSpot x="0.5" y="0.5" xunits="fraction" yunits="fraction" />'
      write (iout,'(a)')  '         </IconStyle>'
      write (iout,'(a)')  '      </Style>'
      write (iout,'(a)')  '      <Point>'
      write (iout,'(5a)') '         <coordinates>', trim(adjustl(cxlon)), ',', trim(adjustl(cxlat)), ',0</coordinates>'
      write (iout,'(a)')  '      </Point>'
      write (iout,'(a)')  '   </Placemark>'
   
      return
      
   end subroutine kml_point
 

!*****************************************************************************************
   subroutine kml_prelude (iout, xlath, xlonh, irange)
   
      integer :: iout
      integer :: irange
      real :: xlath
      real :: xlonh
      character(len=8) :: cirange
      character(len=8) :: cxlath
      character(len=8) :: cxlonh
      
      write (cxlath,'(f8.3)') xlath
      write (cxlonh,'(f8.3)') xlonh
      write (cirange,'(i8)') irange
      
      write (iout,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
      write (iout,'(a)') '<kml xmlns="http://www.opengis.net/kml/2.2">'
      write (iout,'(a)') '<Document>'
      
      write (iout,'(a)') '   <Style id="i_red">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/i_red.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_red">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/a_red.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="red">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_red</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_red</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '  </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_yellow">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/i_yellow.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_yellow">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/a_yellow.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="yellow">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_yellow</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_yellow</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_green">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/i_green.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_green">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/a_green.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="green">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_green</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_green</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_blue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/i_blue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_blue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/a_blue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="blue">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_blue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_blue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_skyblue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/i_skyblue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_skyblue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>_kml/a_skyblue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="skyblue">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_skyblue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_skyblue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <LookAt>'
      write (iout,'(3a)') '      <longitude>',trim(adjustl(cxlonh)),'</longitude>'
      write (iout,'(3a)') '      <latitude>',trim(adjustl(cxlath)),'</latitude>'
      write (iout,'(3a)') '      <range>',trim(adjustl(cirange)),'</range>'
      write (iout,'(a)') '      <tilt>0</tilt>'
      write (iout,'(a)') '      <heading>0</heading>'
      write (iout,'(a)') '   </LookAt>'
  
      return
      
   end subroutine kml_prelude
    

!*****************************************************************************************
   subroutine kml_coda (iout)
   
      integer :: iout
  
      write (iout,'(a)') '</Document>'
      write (iout,'(a)') '</kml>'
      
      return
      
   end subroutine kml_coda


!*****************************************************************************************
   subroutine gccel_hypocentroid (it)
   
   ! Writes an entry for the "doc.kml" file used in GCCEL distributions to display hypocentroids
   ! of calibrated clusters. By itself, this entry is not a legal kml file.
   
      integer :: i
      integer :: io_commentary
      integer :: ios
      integer :: it
      integer :: it1
      integer, external :: lenb
      integer, external :: lunit
      real :: lon_test
      character(len=100) :: cluster_name
      character(len=100) :: commentary_path
      character(len=8) :: cxlat
      character(len=9) :: cxlon
      character(len=132) :: line132
      character(len=132) :: msg
      character(len=100) :: outfil
      logical :: op
      
      cluster_name = ' '
      it1 = it + 1
      
      call fyi ('gccel_hypocentroid: writing doc.kml entry')
      outfil = trim(gcat_folder)//dirsym//trim(basename)//'_doc.kml'
      call open_file (io_gccel, outfil, 'new')
      inquire (unit=io_gccel,opened=op)
      if (.not.op) then
         msg = 'gccel_hypocentroid: file '//trim(outfil)//' was not opened'
         call oops (trim(msg))
      end if
      
      ! Cluster name
      i = lenb(datadir) - 1 ! Get the basic cluster name by removing the trailing digit
      cluster_name = datadir(1:i)
      cluster_name(1:1) = achar(iachar(cluster_name(1:1))-32) ! Capitalize first letter of the name
      
      ! Hypocntroid
      write (cxlat,'(f8.4)') lath(it1)
      lon_test = lonh(it1)
      call set_longitude_range (lon_test, longitude_range)
      write (cxlon,'(f9.4)') lon_test
      
      write (io_gccel,'(3a,5a)') (achar(9),i=1,3), '<Folder id="',trim(cluster_name),'"><name>',&
       trim(cluster_name),'</name>'
      write (io_gccel,'(4a,a)') (achar(9),i=1,4), '<LookAt>'
      write (io_gccel,'(5a,3a)') (achar(9),i=1,5), '<longitude>',trim(adjustl(cxlon)),'</longitude>'
      write (io_gccel,'(5a,3a)') (achar(9),i=1,5), '<latitude>',trim(adjustl(cxlat)),'</latitude>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '<range>300000</range>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '<tilt>0</tilt>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '<heading>0</heading>'
      write (io_gccel,'(4a,a)') (achar(9),i=1,4), '</LookAt>'
      write (io_gccel,'(4a,a)') (achar(9),i=1,4), '<Placemark>'
      write (io_gccel,'(5a,3a)') (achar(9),i=1,5), '<name>',trim(cluster_name),'</name>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '<description>'
      write (io_gccel,'(5a,2a)') (achar(9),i=1,5), '<![CDATA[', trim(basename)
      
      ! Add a commentary file, if it exists. The filename (assumed to be in the data directory)
      ! is given as an argument to the gcat command. The file should be hard-wrapped
      ! at a reasonable line length (say, 72 characters), with a maximum of 132 characters.
      commentary_path = trim(datadir)//dirsym//trim(commentary_fname)
      io_commentary = lunit()
      call open_file (io_commentary, commentary_path, 'old')
      do
         read (io_commentary,'(a)',iostat=ios) line132
         if (ios .lt. 0) exit
         write (io_gccel,'(a)') trim(line132)
      end do
      close (io_commentary)
      do i=1,9
         write (io_gccel,'(a)') trim(commentary_buffer(i))
      end do
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), ']]>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '</description>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '<styleUrl>#hypocenter</styleUrl>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '<Point>'
      write (io_gccel,'(6a,5a)') (achar(9),i=1,6), '<coordinates>',trim(adjustl(cxlon)),',',&
       trim(adjustl(cxlat)),'</coordinates>'
      write (io_gccel,'(5a,a)') (achar(9),i=1,5), '</Point>'
      write (io_gccel,'(4a,a)') (achar(9),i=1,4), '</Placemark>'
      write (io_gccel,'(3a,a)') (achar(9),i=1,3), '</Folder>'
      
      close (io_gccel)
      
      return
      
   end subroutine gccel_hypocentroid
   
   
!*****************************************************************************************
   subroutine gccel_epicenters (it)
   
   ! Writes an entry for the "epicenters.kml" file used in GCCEL distributions to display
   ! epicenters of calibrated clusters.
   
      integer :: iev
      integer :: it
      integer :: it1
      character(len=116) :: ev_name
      character(len=100) :: outfil
      
      it1 = it + 1
      
      call fyi ('gccel_epicenters: writing epicenters.kml entry')
      outfil = trim(gcat_folder)//dirsym//trim(basename)//'_epicenters.kml'
      call open_file (io_gccel, outfil, 'new')
      
      do iev = 1,nev
         ev_name = ' '
         if (iev .lt. 10) then
            write (ev_name,'(3a,i1)') 'cec_', trim(basename), '_', iev
         else if (iev .lt. 100) then
            write (ev_name,'(3a,i2)') 'cec_', trim(basename), '_', iev
         else
            write (ev_name,'(3a,i3)') 'cec_', trim(basename), '_', iev
         end if
         
         call write_epicenter_kml (io_gccel, ev_name, hypocenter_list(iev))
      end do
      
      close (io_gccel)
      return
      
   end subroutine gccel_epicenters
   
   
!*****************************************************************************************
   subroutine write_epicenter_kml (unit, ev_name, hypocenter)
   
      integer, intent(in) :: unit
      real :: scale
      real :: scale0
      real :: scale_factor
      real :: xmag
      character(len=4) :: calibration
      character(len=6) :: cscale
      character(len=8) :: cxlat
      character(len=9) :: cxlon
      character(len=5) :: depth
      character(len=116), intent(in) :: ev_name
      character(len=80), intent(in) :: hypocenter
      character(len=3) :: mag
      character(len=22) :: ot
      
      ot = hypocenter(5:26)
      cxlat = hypocenter(28:35)
      cxlon = hypocenter(37:45)
      depth = hypocenter(47:51)
      calibration = hypocenter(53:56)
      mag = hypocenter(58:60)
      read (mag,'(f3.1)') xmag
      
      scale0 = 0.3
      scale_factor = 1.0
      if (xmag .gt. 2.) then
         scale = scale0 + ((xmag - 2.0)/5.0)*scale_factor
      else
         scale = scale0
      end if
      write (cscale,'(f6.3)') scale
   
      write (unit,'(a)') '<Placemark>'
      write (unit,'(2a)') '<description><![CDATA[Event : ',trim(ev_name)
      write (unit,'(2a)') '<br>Date-OT : ', ot
      write (unit,'(4a)') '<br>Lat/Lon : ', cxlat, ', ', cxlon
      write (unit,'(2a)') '<br>Depth : ', depth
      write (unit,'(2a)') '<br>Calibration : ', calibration
      write (unit,'(2a)') '<br>Magnitude : ', mag
      write (unit,'(a)') ']]></description>'
      write (unit,'(a)') '<styleUrl>#event</styleUrl>'
      write (unit,'(3a)') '<Style><IconStyle><scale>',trim(adjustl(cscale)),'</scale></IconStyle></Style>'
      write (unit,'(5a)') '<Point><coordinates>',trim(adjustl(cxlon)),',',trim(adjustl(cxlat)),',0</coordinates></Point>'
      write (unit,'(a)') '</Placemark>'
      
      return
      
   end subroutine write_epicenter_kml


!*****************************************************************************************
end module mlocout_kml

