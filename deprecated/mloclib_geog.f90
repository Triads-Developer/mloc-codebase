!> Procedures related to geographical calculations

module mloclib_geog

   use declare_limits
   use declare_configuration
   use declare_constants
   use declare_lun
   use declare_phase_data
   use declare_stations
   
   implicit none
   save
   
contains

!*****************************************************************************************
real function dgkmlogc (atmp)

   real :: atmp

   dgkmlogc = dgkmla/sin(atmp*rpd) ! geocentric coordinates

   return
   
end function dgkmlogc

      
!*****************************************************************************************
real function dgkmlo (atmp)

   real :: atmp

   dgkmlo = dgkmla/cos(atmp*rpd) ! geographic coordinates

   return
   
end function dgkmlo


!*****************************************************************************************
subroutine azgap (iev, openaz)

! Largest open angle from a set of azimuths. Only uses azimuth
! values for which the corresponding value of fltr is negative. In single
! mode, the open azimuth is based on fltrh. For a cluster, it is based
! on fltrc.

   integer :: i
   integer :: iev
   integer, dimension(npmax) :: indx
   integer :: j
   integer :: k
   integer :: naz
   real :: azdif
   real, dimension(npmax) :: aztest
   real :: gap
   real :: openaz
   character(len=132) :: msg
   logical :: fltrtest
      
   naz = 0
   do i = 1, nst(iev) ! Form array of azimuths to test
      fltrtest = fltrc(iev,i)
      if (.not.fltrtest) then
         naz = naz + 1
         aztest(naz) = azes(iev,i)
         if (aztest(naz) .lt. 0. .or. aztest(naz) .ge. 360.) then
            if (aztest(naz) .lt. 0.) then
               do while (aztest(naz) .lt. 0.)
                  aztest(naz) = aztest(naz) + 360.
               end do
            else if (aztest(naz) .ge. 360.) then
               do while (aztest(naz) .ge. 360.)
                  aztest(naz) = aztest(naz) - 360.
               end do
            end if
         end if
      end if
   end do
   if (naz .le. 1) then
      openaz = 360.
      write (msg,'(a,i3,a)') 'azgap: naz = ', naz, '; openaz = 360'
      call warnings (trim(msg))
      return
   end if

   if (naz .ge. 1) then
      call indexx (naz, aztest, indx)
   else
      write (msg,'(a,i4,a)') 'azgap: illegal value for naz (', naz, ')'
      call oops (trim(msg))
   end if

   openaz = 0.
   do i = 1,naz-1
      j = indx(i)
      k = indx(i+1)
      azdif = aztest(k) - aztest(j)
      openaz = amax1(openaz,azdif)
   end do

   ! Check gap between largest and smallest azimuths
   gap = aztest(indx(1)) + 360. - aztest(indx(naz))
   openaz = amax1(openaz,gap)

   return
   
end subroutine azgap


!*****************************************************************************************
subroutine set_longitude_range (lon, longitude_range)

! Keep longitude in the desired range, defined by the value of input variable
! 'longitude_range', which gives the center longitude of the desired range.

   integer, intent(in) :: longitude_range
   real, intent(inout) :: lon
   real :: temp
   
   if (lon .lt. real(longitude_range-180)) then
      temp = lon
      lon = lon + 360.
      if (verbose_log) write (io_log,'(a,f8.3,a,f8.3)') 'set_longitude_range: ', temp,&
       ' changed to ', lon
   end if
   
   if (lon .ge. real(longitude_range+180)) then
      temp = lon
      lon = lon - 360.
      if (verbose_log) write (io_log,'(a,f8.3,a,f8.3)') 'set_longitude_range: ', temp,&
       ' changed to ', lon
   end if
         
   return
   
end subroutine set_longitude_range
 

!*****************************************************************************************
subroutine geocen (glat, glon, elat, elats, elatc, elon, elons, elonc)

!  Converts geographic to geocentric coordinates. Latitude is measured 
!  positive to the south from the north pole, longitude is positive to the
!  east from the prime meridian.
     
!  Input:
!    glat   latitude in geographic coordinates
!    glon   longitude in geographic coordinates

!  Output:
!    elat   geocentric latitude
!    elats  sin(geocentric latitude)
!    elatc  cos(geocentric latitude)
!    elon   geocentric longitude
!    elons  sin(geocentric longitude)
!    elonc  cos(geocentric longitude)

   real, intent(out) :: elat
   real, intent(out) :: elatc
   real, intent(out) :: elats
   real, intent(out) :: elon
   real, intent(out) :: elonc
   real, intent(out) :: elons
   real, intent(in) :: glat
   real, intent(in) :: glon
         
   elat = 90.0 - dpr*atan(0.9932773*(sin(glat*rpd)/cos(glat*rpd)))
   elats = sin(elat*rpd)
   elatc = cos(elat*rpd)
   elon = glon
   if (glon .lt. 0.0) elon = 360.0 + glon
   elons = sin(elon*rpd)
   elonc = cos(elon*rpd)
   
   return
   
end subroutine geocen
      
      
!*****************************************************************************************
subroutine geogra (elat, glat)
      
! Input: elat    geocentric latitude
! Output glat    geographic latitude
! Engdahl's code

   real, intent(in) :: elat
   real :: elatc
   real :: elats
   real, intent(out) :: glat

   elats = sin(elat*rpd)
   elatc = cos(elat*rpd)
   glat = dpr*atan(1.0067682*(elatc/elats))
   
   return
   
end subroutine geogra


!*****************************************************************************************
subroutine delaz (eqlt, eqln, stlt, stln, delta, deltdg, deltkm, azeqst, azesdg, azsteq, azsedg, i)

!  compute distance and azimuth from earthquake (eq) to station (st).

!  delta  = distance between (eq) and (st) in radians.
!  azeqst = azimuth from (eq) to (st) clockwise from north in radians
!  azsteq = azimuth from (st) to (eq) clockwise from north in radians
!  cos(delta)  = aa
!  sin(delta)  = bb
!  sin(azeqst) = cc
!  sin(azsteq) = dd

!  i=0 if input coordinates are geographic degrees.
!  i=1 if input coordinates are geocentric radians.

!  subroutines called:
!    coortr (appended)

   integer :: deltm ! epicentral distance in meters
   integer, intent(in) :: i
   real :: angle
   real, intent(out) :: azeqst
   real, intent(out) :: azesdg
   real, intent(out) :: azsedg
   real, intent(out) :: azsteq
   real :: cc
   real :: dd
   real, intent(out) :: delta
   real, intent(out) :: deltdg
   real, intent(out) :: deltkm
   real, external :: dgkmlo
   real :: dx
   real :: dy
   real, intent(inout) :: eqln
   real :: eqlnrd
   real, intent(inout) :: eqlt
   real :: eqltrd
   real :: eqstln
   real, intent(inout) :: stln
   real :: stlnrd
   real, intent(inout) :: stlt
   real :: stltrd
   double precision :: aa
   double precision :: bb
   double precision :: dce
   double precision :: dces
   double precision :: dcs
   double precision :: deqlt
   double precision :: desln
   double precision :: dse
   double precision :: dses
   double precision :: dss
   double precision :: dstlt
   character(len=132) :: msg
   
   ! Convert to spherical polar coordinates in radians.
   if (i .eq. 0) then
      call coortr (eqltrd, eqlnrd, eqlt, eqln, i)
      call coortr (stltrd, stlnrd, stlt, stln, i)
      eqltrd = pi/2. - eqltrd
      stltrd = pi/2. - stltrd
      dy = (stlt - eqlt)/dgkmla ! latitude difference in km
      dx = (stln - eqln)/dgkmlo(eqlt) ! longitude difference in km
   else if (i .eq. 1) then
!         eqltrd=pi/2.-eqlt
      eqltrd = eqlt
      eqlnrd = eqln
!         stltrd=pi/2.-stlt
      stltrd = stlt
      stlnrd = stln
      dy = (stltrd - eqltrd)*radius ! latitude difference in km
      dx = (stlnrd - eqlnrd)*radius ! longitude difference in km
   else
      write (msg,'(a,i3)') 'delaz: illegal index: ', i
      call oops (trim(msg))
   end if
   deltkm = sqrt(dx*dx + dy*dy)
   
   if (deltkm .lt. 1.e-3) then ! Bail out if the separation is less than 1 meter
   
      azeqst = 0.
      azsteq = 0.
      azesdg = 0.
      azsedg = 0.
      delta = 0.
      deltdg = 0.
      deltkm = 0.
      if (verbose_screen) call fyi ('delaz: identical eq and st; delta and azimuth set to zero')
      return
      
   else if (deltkm .lt. 1.) then ! Do short distances (< 1 km) in plane geometry
   
      angle = atan2(dy,dx)
      azeqst = (pi*0.5) - angle ! Convert to azimuth, clockwise from north, still in radians
      if (azeqst .lt. 0.) azeqst = pi*2. + azeqst
      azsteq = azeqst - pi
      if (azsteq .lt. 0.) azsteq = pi*2. + azsteq
      azesdg = azeqst*dpr
      azsedg = azsteq*dpr
      deltdg = deltkm*dgkmla
      delta = deltdg*rpd
      deltm = nint(deltkm*1.0e3)
      if (verbose_screen) then
         write (msg,'(a,i4,a,2(f5.3,a))') 'delaz: separation is only ', deltm, ' meters (',&
          dx, ',', dy, ')'
         call fyi (trim(msg))
         write (*,'(t5,4f10.4)') eqlt, eqln, stlt, stln
      end if

   else ! Traditional calculation in spherical coordinates

      eqstln = stlnrd - eqlnrd
      desln = dble(eqstln)
      deqlt = dble(eqltrd)
      dstlt = dble(stltrd)
      dse = dsin(deqlt)
      dce = dcos(deqlt)
      dss = dsin(dstlt)
      dcs = dcos(dstlt)
      dses = dsin(desln)
      dces = dcos(desln)
      aa = dce*dcs + dse*dss*dces
      bb = dsqrt(1.0d0-aa*aa)
      cc = sngl(dss*dses/bb)
      dd = sngl(-dse*dses/bb)
      delta = sngl(datan2(bb,aa))
      azeqst = asin(cc)
      azsteq = asin(dd)
      if ((dse*dcs - dce*dss*dces) .lt. 0.) azeqst = pi - azeqst
      if (azeqst .lt. 0.) azeqst = pi*2. + azeqst
      if ((dce*dss - dse*dcs*dces) .le. 0.) azsteq = pi - azsteq
      if (azsteq .lt. 0.) azsteq = pi*2. + azsteq
      deltdg = delta*dpr
      azesdg = azeqst*dpr
      azsedg = azsteq*dpr
      deltkm = delta*radius
   
   end if
   
   return
   
end subroutine delaz


!*****************************************************************************************
subroutine coortr (alatrd, alonrd, alatdg, alondg, i)
      
!  alatrd (geocentric latitude in radians)
!  alonrd (geocentric longitude in radians)
!  alatdg (geographical latitude in degrees)
!  alogdg (geographical latitude in degrees)

!  if i=1, transformation from radians to degrees
!  if i=0, transformation from degrees to radians

!  no subroutines called

   integer, intent(in) :: i
   real :: aaa
   real :: alat2
   real, intent(inout) :: alatdg
   real, intent(inout) :: alatrd
   real, intent(inout) :: alondg
   real, intent(inout) :: alonrd
   real :: bbb
   character(len=132) :: msg

   if (i .eq. 0) then ! Degrees to radians
      alatrd = alatdg*rpd
      alonrd = alondg*rpd
      bbb = abs(alatdg)
      if (bbb .lt. 89.9) then
         aaa = 0.9933056*tan(alatrd)
         alatrd = atan(aaa)
      end if
   else if (i .eq. 1) then ! Radians to degrees
      bbb = abs(alatrd)
      if (bbb .lt. 1.57) then
         aaa = tan(alatrd)/0.9933056
         alat2 = atan(aaa)
      else 
         alat2 = alatrd
      end if
      alatdg = alat2*dpr
      alondg = alonrd*dpr
   else
      write (msg,'(a,i9)') 'coortr: illegal value for i: ', i
      call oops (trim(msg))
   end if

   return
   
end subroutine coortr


!*****************************************************************************************
subroutine surface_focus (phasein, depth, delta, dtddel, sfd)

! Input: phasein (phase name)
!        depth (focal depth in km)
!        delta (epicentral distance in degrees)
!        dtdel (ray parameter, in seconds per degree)

! Output: sfd (epicentral distance with distance to first bouncepoint removed)

   integer :: ierr
   integer, external :: iupcor
   real :: bpdel
   real :: bptim
   real, intent(in) :: delta
   real, intent(in) :: depth
   real, intent(in) :: dtddel
   real, intent(out) :: sfd
   character(len=8), intent(in) :: phasein

   sfd = delta
   if (depth .gt. 0.) then
      ierr = iupcor(phasein(1:1), abs(dtddel), bpdel, bptim)
      if (ierr .lt. 0) call warnings ('surface_focus: iupcor failed')
      if (dtddel .lt. 0.0) bpdel = -bpdel
      if (phasein(1:1) .eq. 'p' .or. phasein(1:1) .eq. 's') bpdel = -bpdel
      sfd = delta + bpdel
   end if

   !print *,'surface_focus: phasein depth delta sfd = ', phasein, depth, delta, sfd

   return
   
end subroutine surface_focus
      

!*****************************************************************************************
real function dlttrn (depth)

!  for a ray with a turning depth of 'depth', 'dlttrn' is the epicentral
!  distance from the turning point to where the ray reaches the earths
!  surface (degrees). it is calculated with respect to the 1968 herrin
!  tables for depths less than 800 km.

   integer :: i
   integer :: j
   real, dimension(17) :: d = (/0.0,3.0,6.5,7.6,8.3,8.7,9.2,9.7,10.1,10.4,10.7,11.2,11.8,12.2,13.2,14.2,16.5/)
   real, intent(in) :: depth
   real :: x
   
   x = depth/50. + 1.
   i = int(x)
   j = i + 1
   x = x - i
   dlttrn = d(i) + x*(d(j)-d(i))

   return
   
end function dlttrn


!*****************************************************************************************
end module mloclib_geog

