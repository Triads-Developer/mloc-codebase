!> Procedures to write output text files.

module mloc_out

   use mloc_declare
   use mloc_lib
   use mloc_math
   use mloc_taup
   use mloc_plots
   use mloc_tt
   use mloc_hyposat
   use mloc_mnf

   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine write_summary (it)
   
   !   Printed summary output (.summary)
   
   ! January 26, 1989 by eric bergman
   ! Modified 9/20/94 by eab
   ! Individual event data importances added, 11/14/2001 by eab

      integer :: hourhr
      integer :: i
      integer :: iacvc
      integer :: iacvf
      integer :: iacvi
      integer :: iazc
      integer :: iazes
      integer :: iazh
      integer :: iazp
      integer :: iazt
      integer :: ibin
      integer :: iev
      integer :: it
      integer :: it1
      integer :: j
      integer :: minhr
      integer :: ndatd
      integer, dimension(0:n_it_max1) :: ndatdlh
      integer :: ndatflh
      integer, dimension(0:n_it_max1) :: ndath
      integer :: ndatm
      integer :: ndatp
      integer, dimension(0:n_it_max1) :: ndatprh
      integer :: ndats
      integer :: ndatt
      integer :: ndatx
      integer :: ndtotal
      integer, dimension(8) :: values
      real :: azeqst
      real :: azesdg
      real :: azsedg
      real :: azsteq
      real :: cal_input_dep
      real :: cal_input_lat
      real :: cal_input_lon
      real :: cal_input_s
      real :: cv1f
      real :: cv1i
      real :: cv2f
      real :: cv2i
      real :: dcal_dep
      real :: dcal_lat
      real :: dcal_lon
      real :: dcal_ot
      real :: dcvc
      real :: dcvf
      real :: dcvi
      real :: delta
      real :: deltdg
      real :: deltkm
!      real :: dgkmlo
      real :: dlatkm
      real :: dlatpkm
      real :: dlonkm
      real :: dlonpkm
      real :: hdc
      real :: hdh
      real :: hdp
      real :: hdt
!      real :: hms2s
      real :: lon_test
      real :: lon_test1
      real :: lon_test2
      real :: nsvsprd
      real :: sechr
      real :: sol_input_dep
      real :: sol_input_lat
      real :: sol_input_lon
      real :: sol_input_s
      real, dimension(4) :: sumdx
      real, dimension(4) :: sumdx0
      real, dimension(4) :: sumdxp
      real :: xsigma
      character(len=3), dimension(12) :: binlbl
      character(len=160) :: command_line
      character(len=8) :: date
      character(len=5) :: depthprh
      character(len=12) :: dimpcode
      character(len=3) :: latprh
      character(len=3) :: lonprh
      character(len=100) :: outfil
      character(len=10) :: time
      character(len=11) :: timeprh
      character(len=5) :: zone
      
      data binlbl/'NNE',' NE','ENE','ESE',' SE','SSE','SSW',' SW', 'WSW','WNW',' NW','NNW'/
      data dimpcode/'            '/
      
      it1 = it + 1
      
      outfil = trim(outfile)//'.summary'
      call open_file (io_out, outfil, 'new')
      
      ! Various data sets and procedures used in the analysis.
      if (debug) then
         write (io_out,'(4a)') ' Program: ', trim(mloc_version), '; status: '//trim(mloc_status), '; running in test mode'
      else
         write (io_out,'(3a)') ' Program: ', trim(mloc_version), '; status: '//trim(mloc_status)
      end if   
      write (io_out,'(2a)') ' Run: ', trim(basename)
      write (io_out,'(2a)') ' Author: ', trim(mloc_author)
      call date_and_time (date, time, zone, values)
      write (io_out,'(4a)') ' Date: ', date, ' ', time
      if (locmod) then
         write (io_out,'(5a)') ' Travel times: ', trim(taup_model), ' (', trim(locmodfname), ')'
      else
         write (io_out,'(2a)') ' Travel times: ', trim(taup_model)
      end if
      write (io_out,'(a,f6.2,a,f6.2,a,f6.2,a)') ' P_ travel times: ', p__a, ' + ', p__b, ' * delta from ', p__min, ' degrees'
      write (io_out,'(a,f6.2,a,f6.2,a,f6.2,a)') ' Lg travel times: ', lg_a, ' + ', lg_b, ' * delta from ', lg_min, ' degrees'
      write (io_out,'(a,f6.2,a,f6.2,a,f6.2,a)') ' Rg travel times: ', rg_a, ' + ', rg_b, ' * delta from ', rg_min, ' degrees'
      write (io_out,'(a,f6.2,a,f6.2,a)') ' T-phase travel times: ', tphase_a, ' + ', tphase_b, ' * delta'
      select case (tt_corr)
         case (0)
            write (io_out,'(2a)') ' Station elevation correction: off',&
             ' (focal depth referenced to surface in direct calibration)'
         case (1)
            write (io_out,'(2a)') ' Station elevation correction: on',&
             ' (focal depth referenced to geoid in direct calibration)'
      end select
      if (suppstn) then
         do i = 1,n_supp_stn_file
            write (io_out,'(2a)') ' Supplemental station file: ', trim(suppfilnam(i))
         end do
      end if
      if (n_skip .gt. 0) then
         do i = 1,n_skip
            write (io_out,'(a,1x,a8,1x,a8,1x,a8)') ' Skipping station/phase/author :', (skip_params(i,j),j=1,3)
         end do
      end if
      write (io_out,'(2a)') ' Data flags used: ', dflag_pr
      write (io_out,'(2a)') ' Residuals weighted by reading errors: ', data_weight_pr
      write (io_out,'(2a)') ' Assume perfect theoretical TTs for hypocentroid (ttsprd = 0): ', pttt_pr
      if (data_weight) write (io_out,'(2a)') ' Reading errors: ', trim(rderrfname)
      write (io_out,'(a,3(f5.2,a))') ' Minimum allowed reading errors: ', rderr_min_loc,&
       ' (local) /', rderr_min, ' (general) /', rderr_min_depth, ' (depth phases)'
      if (rels_set) then
         write (io_out,'(a,2(f5.2,a))') ' Reading errors for local phases: ', rderr_loc_p,&
          ' (P) /', rderr_loc_s, ' (S)'
      else
         write (io_out,'(a)') ' Reading errors for local phases: Not set'      
      end if
      write (io_out,'(a,f6.2)') ' Cluster vector fudge factor (non-gaussian): ', radius_cvff
      write (io_out,'(2a)') ' Travel time spread: ', trim(ttsprdfname)
      write (io_out,'(a,3f5.1)') ' Windowing: ', wind1, wind2, windloclim
      if (damping) then
         write (io_out,'(a)') ' Cluster vector damping: on'
      else
         write (io_out,'(a)') ' Cluster vector damping: off'
      end if
      if (tikhonov_factor .gt. 0.) write (io_out,'(a,f5.2)') ' Tikhonov regularization factor: ', tikhonov_factor
      if (bias_corr) then
         write (io_out,'(a)') ' Hypocentroid bias correction: on'
      else
         write (io_out,'(a)') ' Hypocentroid bias correction: off'
      end if
      if (read_hdf) write (io_out,'(2a)') ' Starting locations: ', trim(rhdf_filnam)
      
      if (direct_cal) then
         write (io_out,'(a)') ' Direct calibration: on'
      else
         write (io_out,'(a)') ' Direct calibration: off'
      end if
      
      ! Limits used to select arrival time data.
      write (io_out,'(a)') ' Distance ranges:'
      if (ponly) then
         if (nhlim .eq. 1) write (io_out,'(4x,a,1(2f6.1,3x),a)')&
          'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim), ' using P only'
         if (nhlim .eq. 2) write (io_out,'(4x,a,2(2f6.1,3x),a)')&
          'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim), ' using P only'
         if (nhlim .eq. 3) write (io_out,'(4x,a,3(2f6.1,3x),a)')&
          'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim), ' using P only'
      else
         write (io_out,'(4x,a,3(2f6.1,3x))') 'Hypocentroid   : ', (hlim(i,1),hlim(i,2),i=1,nhlim)
      end if
      write (io_out,'(4x,a,3(2f6.1,3x))') 'Cluster vectors: ', (clim(i,1),clim(i,2),i=1,nclim)
      
! Hypocentroid data

      ! Location parameters.
      if (latfh) then
         latprh = 'lat'
      else
         latprh = '   '
      End if
      if (lonfh) then
         lonprh = 'lon'
      else
         lonprh = '   '
      end if
      if (depthfh) then
         depthprh = 'depth'
      else
         depthprh = '     '
      end if
      if (timefh) then
         timeprh = 'origin time'
      else
         timeprh = '           '
      end if
      write (io_out,'(a,2(a3,a),a5,a,a11,a)') ' Hypocentroid location parameters: (', latprh, ', ', lonprh, ', ', depthprh,&
       ', ', timeprh, ')'
      
      ! Starting point for relocation.
      write (io_out,'(/25x,a,8x,a,5x,a,7x,a,10x,a,3x,a)') 'TIME', 'LATITUDE', 'LONGITUDE', 'DEPTH', 'E**2/   NSTA', 'SHAT'
      lon_test = lonh(0)
      call set_longitude_range (lon_test, longitude_range)
      write (io_out,'(a,16x,2(i2,a),f5.2,f12.3,f14.3,f12.1/)') ' START', hourh(0), ':', minh(0), ':', sech(0), lath(0), lon_test,&
       depthh(0)
      
      ! Corrections at each iteration.
      do i = 0,it
         ndath(i) = 0
         ndatdlh(i) = 0
         ndatprh(i) = 0
         ndatflh = 0
         do iev = 1,n_event
            ndath(i) = ndath(i) + ndat(iev,i)
            ndatdlh(i) = ndatdlh(i) + ndatdl(iev,i)
            ndatprh(i) = ndatprh(i) + ndatpr(iev,i)
            ndatflh = ndatflh + ndatfl(iev)
         end do
         write (io_out,'(a,i2,20x,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,2x,f8.1,a,i7,2x,f5.2)') ' ITER', i, delx0(4,i),&
          ' S', -delx0(1,i), ' KM', delx0(2,i), ' KM', delx0(3,i), ' KM', ehatsqh(i), '/', ndath(i), shath(i)
      end do
      
      ! Bias correction term
      if (.not. bias_corr) then
         write (io_out,'(a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a)') ' BIAS', bcorr(4), ' S', bcorr(1), ' KM', bcorr(2), ' KM',&
          bcorr(3), ' KM'
      end if
      
      ! Cumulative difference: solution - starting point
      do i = 1,4
         sumdx0(i) = 0.
      end do
      do i = 0,it
         sumdx0(1) = sumdx0(1) - delx0(1,i) ! Sign change for geocentric latitude
         sumdx0(2) = sumdx0(2) + delx0(2,i)
         sumdx0(3) = sumdx0(3) + delx0(3,i)
         sumdx0(4) = sumdx0(4) + delx0(4,i)
      end do
      if (.not. bias_corr) then
         sumdx0(1) = sumdx0(1) - bcorr(1) ! Sign change for geocentric latitude
         sumdx0(2) = sumdx0(2) + bcorr(2)
         sumdx0(3) = sumdx0(3) + bcorr(3)
         sumdx0(4) = sumdx0(4) + bcorr(4)
      end if
      if ((sumdx0(1)**2 + sumdx0(2)**2) .ge. 0.01) then ! No vector calculated unless length will be at least 100 m
         hdh = sqrt(sumdx0(1)**2 + sumdx0(2)**2)
         iazh = nint(atan2(sumdx0(2),sumdx0(1))*dpr)
      else
         hdh = 0.
         iazh = 0
      end if  
      if (iazh .lt. 0) iazh = iazh + 360
      write (io_out,'(/a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,a,i4,a,4x,a,f4.1,a)') ' CUMULATIVE', sumdx0(4), ' S',&
       sumdx0(1), ' KM', sumdx0(2), ' KM', sumdx0(3), ' KM', 'AZIM = ', iazh, ' DEG', 'DIST = ', hdh, ' KM'
      
      ! Solution
      lon_test = lonh(it1)
      call set_longitude_range (lon_test, longitude_range)
      write (io_out,'(/a,13x,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' SOLUTION', hourh(it1), ':', minh(it1), ':', sech(it1), lath(it1),&
       lon_test, depthh(it1)
      
      ! Standard errors and confidence ellipse
      write (io_out,'(a,t28,f6.2,a,4x,f6.3,a,4x,f6.3,a,3x,f6.2,a,3(f6.1,a))') ' STANDARD ERRORS', sdxhath(4), ' S', sdxhath(1),&
       ' DEG', sdxhath(2), ' DEG', sdxhath(3), ' KM    ELLIPSE:', alphah, ' DEG', xl1h, ' AND ', xl2h, ' KM'
      
      ! Calibration shift
      if (indirect_cal) then
         dlatkm = del_cal_lat/dgkmla
         dlonkm = del_cal_lon/dgkmlo(lath(it1))
         if ((dlatkm**2 + dlonkm**2) .ge. 0.01) then ! No vector calculated unless length will be at least 100 m
            hdh = sqrt(dlatkm**2 + dlonkm**2)
            iazh = nint(atan2(dlonkm,dlatkm)*dpr)
         else
            hdh = 0.
            iazh = 0
         end if  
         if (iazh .lt. 0) iazh = iazh + 360
         write (io_out,'(/a)')&
          ' ***************************************************************************************************************'
         write (io_out,'(a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f4.1,a)') ' SHIFT', del_cal_tim, ' S',&
          del_cal_lat, ' DEG', del_cal_lon, ' DEG', del_cal_dep, ' KM', 'AZIM = ', iazh, ' DEG', 'DIST = ', hdh, ' KM'
          
         call timecr (otsh(it+1)+del_cal_tim, hourhr, minhr, sechr) 
         lon_test = lonh(it1) + del_cal_lon
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(a,t23,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' CALIBRATED', hourhr, ':', minhr, ':', sechr,&
          lath(it1)+del_cal_lat, lon_test, depthh(it1)+del_cal_dep
          
         ! Calibration event correction vectors
         write (io_out,'(/a/)') ' Calibration event correction vectors (weights):'
         do iev = 1,n_event
            if (cal_event(iev,3)) then ! Calculate difference between calibration and inversion solution
               dcal_lat = cal_lat(iev,3) - latp(iev,it1)
               lon_test = lonp(iev,it1)
               call set_longitude_range (lon_test, longitude_range)
               dcal_lon = cal_lon(iev,3) - lon_test
               dcal_dep = cal_dep(iev,3) - depthp(iev,it1)
               dcal_ot = hms2s(cal_hr(iev,3), cal_min(iev,3), cal_sec(iev,3)) - otsp(iev,it1)
               dlatpkm = dcal_lat/dgkmla
               dlonpkm = dcal_lon/dgkmlo(cal_lat(iev,3))
               if ((dlatpkm**2 + dlonpkm**2) .ge. 0.01) then ! No vector calculated unless length will be at least 100 m
                  hdp = sqrt(dlatpkm**2 + dlonpkm**2)
                  iazp = nint(atan2(dlonpkm,dlatpkm)*dpr)
               else
                  hdp = 0.
                  iazp = 0
               end if  
               if (iazp .lt. 0) iazp = iazp + 360
               write (io_out,&
                '(a7,i3,2x,a16,t30,f6.2,a,f4.2,a,t44,f6.3,a,f4.2,a,t58,f6.3,a,f4.2,a,t70,f7.2,a,f4.2,a,t84,a,i4,a,t103,a,f4.1,a)')&
                ' EVENT ', iev, evtnam(iev)(1:16), dcal_ot, '(', w4(iev), ')', dcal_lat, '(', w12(iev), ')',& 
                dcal_lon, '(', w12(iev), ')', dcal_dep, '(', w3(iev), ')', 'AZIM = ', iazp, ' DEG', 'DIST = ', hdp, ' KM'
            end if
         end do
         
      end if
     
! Cluster vectors

      if (n_event .eq. 1) go to 999
      
      write (io_out,'(/a)')&
       ' *************************************************************************************************************'
      write (io_out,'(a26,t30,a4,t42,a13,t57,a5)') ' cluster vector statistics:', 'ITER', 'E**2/NSTA-PHA', 'SHATC'
      do i = 0,it
         ndtotal = 0
         do iev = 1,n_event
            ndtotal = ndtotal + ndatc(iev,i)
         end do
         write (io_out,'(t32,i1,t38,f8.1,a1,i8,t57,f5.2)') i, ehatsqc(i), '/', ndtotal, shatc(i)
      end do
      call croux (shatsqci, n_event, nsvsprd)
      write (io_out,'(a,f6.3)') ' Spread of event normalized sample variances: ', nsvsprd
      
      ! Initial and final cluster vectors
      write (io_out,'(/t40,a,t69,a,t98,a)') 'INITIAL', 'CHANGE', 'FINAL'
      do iev = 1,n_event ! Loop over events
      
         ! Initial cluster vector
         cv1i = (latp(iev,0) - lath(0))/dgkmla
         cv2i = (lonp(iev,0) - lonh(0))/dgkmlo(lath(0))
         if ((cv1i**2 + cv2i**2) .ge. 0.1) then ! No vector calculated unless length is 10 m or more
            dcvi = sqrt(cv1i**2 + cv2i**2)
            iacvi = nint(atan2(cv2i,cv1i)*dpr)
         else
            dcvi = 0.
            iacvi = 0
         end if  
         if (iacvi .lt. 0) iacvi = iacvi + 360
         
         ! Change in cluster vector
         do i = 1,4
            sumdxp(i) = 0.
         end do
         do i = 0,it 
            sumdxp(1) = sumdxp(1) - dxp(iev,1,i) ! Sign change for geocentric latitude
            sumdxp(2) = sumdxp(2) + dxp(iev,2,i)
            sumdxp(3) = sumdxp(3) + dxp(iev,3,i)
            sumdxp(4) = sumdxp(4) + dxp(iev,4,i)
         end do
         if ((sumdxp(1)**2 + sumdxp(2)**2) .ge. 0.1) then ! No vector calculated unless length is 10 m or more
            dcvc = sqrt(sumdxp(1)**2 + sumdxp(2)**2)
            iacvc = nint(atan2(sumdxp(2),sumdxp(1))*dpr)
         else
            dcvc = 0.
            iacvc = 0
         end if 
         if (iacvc .lt. 0) iacvc = iacvc + 360
         
         ! Final cluster vector
         cv1f = (latp(iev,it1) - lath(it1))/dgkmla
         cv2f = (lonp(iev,it1) - lonh(it1))/dgkmlo(lath(it1))
         if ((cv1f**2 + cv2f**2) .ge. 0.1) then ! No vector calculated unless length is 10 m or more
            dcvf = sqrt(cv1f**2 + cv2f**2)
            iacvf = nint(atan2(cv2f,cv1f)*dpr)
         else
            dcvf = 0.
            iacvf = 0
         end if 
         if (iacvf .lt. 0) iacvf = iacvf + 360
         
         write (io_out,'(a7,i3,2x,a16,3(f8.1,a9,i3,a4,5x))') ' EVENT ', iev, evtnam(iev)(1:16), dcvi, ' KM  AT  ', iacvi, ' DEG',&
          dcvc, ' KM  AT  ', iacvc, ' DEG', dcvf, ' KM  AT  ', iacvf, ' DEG'
      end do
      
      ! Data importance distribution, by 30-degree wedges
      write (io_out,'(/a)')&
       ' *******************************************************************************************************************'
      write (io_out,'(a,t30,12(3x,a3))') ' DATA IMPORTANCE:', (binlbl(i),i=1,12)
      write (io_out,'(2x,a,t30,12(f6.3))') 'HYPOCENTROID:', (dimph(i),i=1,12)
      do iev = 1,n_event
         do ibin = 1,12
            if (dimpciev(iev,ibin) .lt. 0.010) then
               dimpcode(ibin:ibin) = '0'
            else
               dimpcode(ibin:ibin) = '-'
            end if
         end do
         write (io_out,'(2x,a,i4,2x,a16,t30,12(f6.3),3x,a)') 'EVENT', iev, evtnam(iev)(1:16), (dimpciev(iev,i),i=1,12), dimpcode
      end do
      write (io_out,'(2x,a,t30,12(f6.3))') 'CLUSTER VECTORS:', (dimpc(i),i=1,12)
      
      ! Cluster vector info about each event.
      shatsqc=shatc(it)*shatc(it)
      do iev=1,n_event ! Loop over events
         write (io_out,'(/2a)') ' ***********************************************************************************************',&
         '*******************************'
         write (io_out,'(a,i3,t20,a,t60,a)') ' CLUSTER EVENT ', iev, trim(evtnam(iev)), trim(infile(iev))
         if (rmag(iev) .gt. 0.) then
            write (io_out,'(a,f3.1)') ' M', rmag(iev)
         else
            write (io_out,'(a)') ' No magnitude'
         end if
         write (io_out,'(a,2(a3,a),a5,a,a11,a)') ' Cluster vector location parameters: (', latpr(iev), ', ', lonpr(iev),&
          ', ', depthpr(iev), ', ', timepr(iev), ')'
         
         ! Input hypocenter from data file and starting point for relocation
         write (io_out,'(t13,a,9x,a,8x,a,5x,a,7x,a,9x,a,4x,a)') 'DATE', 'TIME', 'LATITUDE', 'LONGITUDE',&
          'DEPTH', 'EC**2/NSTA', 'EH**2/NSTA'
         lon_test = lon_inp(iev)
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(a,t10,2(i2,a),i4,3x,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' INPUT', mone(iev), '/', idye(iev), '/', iyre(iev),&
          hour_inp(iev), ':', min_inp(iev), ':', sec_inp(iev), lat_inp(iev), lon_test, depth_inp(iev,1)
         lon_test = lonp(iev,0)
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(a,t10,2(i2,a),i4,3x,2(i2,a),f5.2,f12.3,f14.3,f12.1/)') ' START', mone(iev), '/', idye(iev), '/',&
          iyre(iev), hourp(iev,0), ':', minp(iev,0), ':', secp(iev,0), latp(iev,0), lon_test, depthp(iev,0)
         
         ! Corrections at each iteration
         do i = 0,it
            write (io_out,'(a,i2,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,f6.1,a,i4,3x,f6.1,a,i4)') ' ITER', i,&
             dxp(iev,4,i), ' S', -dxp(iev,1,i), ' KM', dxp(iev,2,i), ' KM', dxp(iev,3,i), ' KM', eciev(iev,i), '/',&
             ndatc(iev,i), ehiev(iev,i), '/', ndat(iev,i)
         end do
         
         ! Cumulative difference: solution - starting point
         do i = 1,4
            sumdxp(i) = 0.
         end do
         do i = 0,it 
            sumdxp(1) = sumdxp(1) - dxp(iev,1,i) ! Sign change for geocentric latitude
            sumdxp(2) = sumdxp(2) + dxp(iev,2,i)
            sumdxp(3) = sumdxp(3) + dxp(iev,3,i)
            sumdxp(4) = sumdxp(4) + dxp(iev,4,i)
         end do
         if ((sumdxp(1)**2 + sumdxp(2)**2) .ge. 0.01) then ! No vector calculated unless length is 100 m or more
            hdc = sqrt(sumdxp(1)**2 + sumdxp(2)**2)
            iazc = nint(atan2(sumdxp(2),sumdxp(1))*dpr)
         else
            hdc = 0.
            iazc = 0
         end if
         if (iazc .lt. 0) iazc = iazc + 360
         
         write (io_out,'(/a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,a,i4,a,4x,a,f4.1,a)') '  CUMULATIVE', sumdxp(4), ' S',&
          sumdxp(1), ' KM', sumdxp(2), ' KM', sumdxp(3), ' KM', 'AZIM = ', iazc, ' DEG', 'DIST = ', hdc, ' KM'
         
         ! Hypocentroid change
         write (io_out,'(a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a)') ' +HYPOCENTROID', sumdx0(4), ' S', sumdx0(1), ' KM',&
          sumdx0(2), ' KM', sumdx0(3), ' KM'
         
         ! Total change (cluster + hypocentroid)
         do i = 1,4
            sumdx(i) = sumdx0(i) + sumdxp(i)
         end do
         if ((sumdx(1)**2 + sumdx(2)**2) .ge. 0.01) then ! No vector calculated unless length is 100 m or more
            hdt = sqrt(sumdx(1)**2 + sumdx(2)**2)
            iazt = nint(atan2(sumdx(2),sumdx(1))*dpr)
         else
            hdt = 0.
            iazt = 0
         end if  
         if (iazt .lt. 0) iazt = iazt + 360
         
         write (io_out,'(a,t28,f6.2,a,3x,f6.2,a,5x,f6.2,a,4x,f7.2,a,4x,a,i4,a,4x,a,f4.1,a)') '  TOTAL', sumdx(4), ' S', sumdx(1),&
          ' KM', sumdx(2), ' KM', sumdx(3), ' KM', 'AZIM = ', iazt, ' DEG', 'DIST = ', hdt, ' KM'
         
         ! Solution
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         write (io_out,'(/a,2(i2,a),i4,2x,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' SOLUTION ', mone(iev), '/', idye(iev), '/',&
          iyre(iev), hourp(iev,it1), ':', minp(iev,it1), ':', secp(iev,it1), latp(iev,it1), lon_test, depthp(iev,it1)
         
         ! Standard errors and confidence ellipse
         write (io_out,'(a,t28,f6.2,a,4x,f6.3,a,4x,f6.3,a,3x,f6.2,a,3(f6.1,a))') ' STANDARD ERRORS', sdxhatc(iev,4), ' S',&
          sdxhatc(iev,1), ' DEG', sdxhatc(iev,2), ' DEG', sdxhatc(iev,3), ' KM    ELLIPSE:', alphac(iev), ' DEG', xl1c(iev),&
          ' AND ', xl2c(iev), ' KM'
         
         ! Change from input hypocenter (data file)
         sol_input_s = (real(minp(iev,it1)*60) + secp(iev,it1)) - (real(min_inp(iev)*60) + sec_inp(iev))
         sol_input_lat = latp(iev,it1) - lat_inp(iev)
         lon_test1 = lonp(iev,it1)
         call set_longitude_range (lon_test1, longitude_range)
         lon_test2 = lon_inp(iev)
         call set_longitude_range (lon_test2, longitude_range)
         sol_input_lon = lon_test1 - lon_test2
         sol_input_dep = depthp(iev,it1) - depth_inp(iev,1)
         if (verbose_log) then
            write (io_log,'(a,i4,4f10.4)') 'write_summary: ', iev, lat_inp(iev), lon_test2, latp(iev,it1), lon_test1
            write (io_log,'(a,i4, 2f10.4)') 'write_summary: ', iev, sol_input_lat, sol_input_lon
         end if
         call delaz (lat_inp(iev), lon_test2, latp(iev,it1), lon_test1, delta, deltdg, deltkm, azeqst, azesdg, azsteq,&
          azsedg, 0)
         if (azesdg .lt. 0) azesdg = azesdg + 360.
         iazes = nint(azesdg)
         if (verbose_log) write (io_log,'(a,i4,f10.4,i4)') 'write_summary: ', iev, deltkm, iazes
         write (io_out,'(/a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f5.1,a)') ' SOLUTION-INPUT', sol_input_s,&
          ' S', sol_input_lat, ' DEG', sol_input_lon, ' DEG', sol_input_dep, ' KM', 'AZIM = ', iazes, ' DEG', 'DIST =',&
         deltkm, ' KM'
         
         ! Calibration
         if (indirect_cal) then
            write (io_out,'(/a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f4.1,a)') ' SHIFT', del_cal_tim, ' S',&
             del_cal_lat, ' DEG', del_cal_lon, ' DEG', del_cal_dep, ' KM', 'AZIM = ', iazh, ' DEG', 'DIST = ', hdh, ' KM'
            call timecr (otsp_cal(iev), hourhr, minhr, sechr) 
            lon_test = lonp_cal(iev)
            call set_longitude_range (lon_test, longitude_range)
            write (io_out,'(a,t23,2(i2,a),f5.2,f12.3,f14.3,f12.1)') ' CALIBRATED', hourhr, ':', minhr, ':', sechr,&
             latp_cal(iev), lon_test, depthp_cal(iev)
            ! Change from input hypocenter (data file)
!            cal_input_s = (real(minp(iev,it1)*60) + secp(iev,it1)) - (real(min_inp(iev)*60) + sec_inp(iev))
            cal_input_s = (real(minhr*60) + sechr) - (real(min_inp(iev)*60) + sec_inp(iev))
            cal_input_lat = latp_cal(iev) - lat_inp(iev)
            lon_test1 = lonp_cal(iev)
            call set_longitude_range (lon_test1, longitude_range)
            lon_test2 = lon_inp(iev)
            call set_longitude_range (lon_test2, longitude_range)
            cal_input_lon = lon_test1 - lon_test2
            cal_input_dep = depthp_cal(iev) - depth_inp(iev,1)
            call delaz (lat_inp(iev), lon_inp(iev), latp_cal(iev), lonp_cal(iev), delta, deltdg, deltkm, azeqst, azesdg,&
             azsteq, azsedg, 0)
            if (azesdg .lt. 0) azesdg = azesdg + 360.
            iazes = nint(azesdg)
            write (io_out,'(/a,t28,f6.2,a,t40,f6.3,a,t54,f6.3,a,t66,f7.2,a,t80,a,i4,a,4x,a,f5.1,a)') ' CALIBRATED-INPUT',&
             cal_input_s, ' S', cal_input_lat, ' DEG', cal_input_lon, ' DEG', cal_input_dep, ' KM', 'AZIM = ', iazes, ' DEG',&
             'DIST =', deltkm, ' KM'
         end if
         
         ! Normalized sample error for cluster vector and how many sigma out from expected value
         xsigma = (shatsqci(iev) - 1.)/nsvsprd
         write (io_out,'(/a,f6.3,a,f7.2,a)') ' Normalized sample variance for cluster vector: ', shatsqci(iev), ' (', xsigma, ')'
         
         ! Flag statistics
         ndatx = 0
         ndatp = 0
         ndats = 0
         ndatm = 0
         ndatt = 0
         ndatd = 0
         do i = 1,nst(iev)
            if (fcode(iev,i) .eq. 'x') then
               ndatx = ndatx + 1
            else if (fcode(iev,i) .eq. 'p') then
               ndatp = ndatp + 1
            else if (fcode(iev,i) .eq. 's') then
               ndats = ndats + 1
            else if (fcode(iev,i) .eq. 'm') then
               ndatm = ndatm + 1
            else if (fcode(iev,i) .eq. 't') then
               ndatt = ndatt + 1
            else if (fcode(iev,i) .eq. 'd') then
               ndatd = ndatd + 1
            end if
         end do
         write (io_out,'(/a)') ' Flagged readings:'
         if (ndatx .gt. 0) write (io_out,'(2x,a,i6,a)') 'x', ndatx, ' for large residual (absolute or cluster)'
         if (ndatp .gt. 0) write (io_out,'(2x,a,i6,a)') 'p', ndatp, ' for problematic phase'
         if (ndats .gt. 0) write (io_out,'(2x,a,i6,a)') 's', ndats, ' by a SKIP command'
         if (ndatm .gt. 0) write (io_out,'(2x,a,i6,a)') 'm', ndatm, ' for missing station coordinates'
         if (ndatt .gt. 0) write (io_out,'(2x,a,i6,a)') 't', ndatt, ' for a timing problem'
         if (ndatd .gt. 0) write (io_out,'(2x,a,i6,a)') 'd', ndatd, ' as a duplicate reading'

      end do ! End loop over events
            
999   continue

      close (io_out)

    ! Put a copy of the ~.summary file in the ~_gccel/ subdirectory
      if (gccelout) then
         command_line = 'cp '//trim(outfil)//' '//trim(gcat_folder)//dirsym
         call execute_command_line (trim(command_line))
      end if
      
      return
      
   end subroutine write_summary
      

!*****************************************************************************************
   subroutine write_hdf (it)
   
   ! Writes one or more files summarizing final locations and uncertainties. There are several scenarios:
   
   ! (1) No calibration: writes a .hdf file, uncertainties for relative location only (cluster vectors).
   ! (2) Direct calibration: writes a .hdf_dcal file, uncertainties for absolute location (cluster vector plus hypocentroid).
   ! (3) Indirect calibration: writes both a .hdf file (as above) and a .hdf_cal file, with full uncertainties.
   ! (4) Direct and indirect calibration: Writes both a .hdf_dcal and a .hdf_cal file. Indirect calibration takes precedence.
   
   ! This subroutine also writes output to the .log file for the command 'subc', creating the main body of a command file
   ! for a subcluster that is especially well-suited for direct calibration.
   
      integer :: calibration_type
      integer :: i
      integer :: iev
      integer :: ii
      integer :: ilast
      integer, dimension(n_arrtim_max) :: indx
      integer, intent(in) :: it
      integer :: it1
      integer :: j
      integer :: jj
      integer :: n_subc
      integer :: n_subc_delt
      integer :: ndatx
      real :: alpha
      real :: avh
      real :: azdif
      real, dimension(n_arrtim_max) :: azesiev
      real :: cstadelc
      real :: ddep
      real :: dot
      real :: lon_test
      real :: min_depth_uncertainty = 0.1 ! To avoid rounding to zero on output
      real :: openaz_limit = 200.
      real :: openazc
      real :: rstadelc
      real :: saaz1
      real :: saaz2
      real :: xl1
      real :: xl2
      character(len=4) :: cal_code
      character(len=1) :: cal2
      character(len=160) :: command_line
      character(len=4) :: fdm
      character(len=4) :: fdp
      character(len=132) :: fmt
      character(len=1) :: free_depth
      character(len=132) :: msg
      character(len=100) :: outfil
      character(len=3) :: rmag_pr
      
      hypocenter_list = ' '
      
      if (direct_cal) then
         outfil = trim(outfile)//'.hdf_dcal' ! Hypocentroid + cluster vector uncertainties
      else
         outfil = trim(outfile)//'.hdf' ! Cluster vector uncertainties only
      end if
      call open_file (io_out, outfil, 'new')
      
      it1 = it + 1
      fmt = '(i4,4i3,f6.2,f10.5,f11.5,f7.2,1x,a1,a1,f6.2,1x,a3,a2,1x,a10,3i5,f7.2,f6.2,2(1x,a4),3f6.1,2(i4,f6.2),f7.1,1x,a4)'
      
      ! Main body of the command file for a subcluster for direct calibration
      if (subc_set) then
         write (io_log,'(/a)') 'Subcluster for direct calibration (command "subc")'
         write (io_log,'(a,f6.2)') ' Epicentral distance limit = ', subc_delt
         write (io_log,'(a,i3)') ' Minimum # of readings within that distance = ', subc_nmin
         write (io_log,'(a,i3)') ' Minimum connectivity with the rest of the cluster = ', subc_nconnect
         n_subc = 0
         do iev = 1,n_event
            n_subc_delt = 0
            do i = 1,nst(iev)
               if (fcode(iev,i) .eq. ' ' .and. delt(iev,i) .le. subc_delt) n_subc_delt = n_subc_delt + 1
            end do
            if (n_subc_delt .ge. subc_nmin .and. ndatc(iev,it) .ge. subc_nconnect) then
               n_subc = n_subc + 1
               write (io_log,'(a)') 'memb'
               write (io_log,'(2a)') 'even ', trim(evtnam(iev))
               write (io_log,'(2a)') 'inpu ', trim(infile20(iev))
            end if
         end do   
         write (msg,'(a,i3,a)') 'mlocout_hdf: ', n_subc, ' events selected by command "subc" for subcluster'
         call fyi (trim(msg))
      end if
      
      if (direct_cal) write (io_log,'(/a)') 'cal_ commands from direct calibration:' ! Logging cal_ commands
      
      do iev = 1,n_event ! Loop over events
      
         if (direct_cal) then
            alpha = alphadc(iev)
            xl1 = xl1dc(iev)
            xl2 = xl2dc(iev)
            ddep = ddepdc(iev)
            dot = dotdc(iev)
            calibration_type = 1
         else
            alpha = alphac(iev)
            xl1 = xl1c(iev)
            xl2 = xl2c(iev)
            ddep = sdxhatc(iev,3)
            dot = sdxhatc(iev,4)
            calibration_type = 0
         end if
         
         ! Calibration code
         call calibration_code (iev, calibration_type, cal_code)
         cal2 = cal_code(2:2)
         call uctolc (cal2,-1)
         
         ! Uncertainty of focal depth
         fdp = ' '
         fdm = ' '
         if (mindx(iev,3) .gt. 0) then ! Free depth solution
            free_depth = 'f'
            if (ddep .le. 99.) then
               write (fdp,'(f4.1)') max(min_depth_uncertainty,ddep)
               write (fdm,'(f4.1)') max(min_depth_uncertainty,ddep)
            end if
         else ! Take from default or assigned uncertainties
            free_depth = ' '
            if (depthp_plus(iev) .ge. min_depth_uncertainty .and. depthp_plus(iev) .le. 99.)&
             write (fdp,'(f4.1)') depthp_plus(iev)
            if (depthp_minus(iev) .ge. min_depth_uncertainty .and. depthp_minus(iev) .le. 99.)&
             write (fdm,'(f4.1)') depthp_minus(iev)
            ddep = max(depthp_minus(iev),depthp_plus(iev))
         end if
         
         ! Magnitude
         if (rmag(iev) .gt. 0.) then
            write (rmag_pr,'(f3.1)') rmag(iev)
         else
            rmag_pr = '   '
         end if
         
         avh = xl1*xl2*3.1415 ! Geometric area of 90% confidence ellipse
         
         ! Make sure the semi-axis azimuths are between 0 and 360 degrees
         saaz1 = alpha
         if (saaz1 .lt. 0.) saaz1 = saaz1 + 360.
         saaz2 = saaz1 + 90.
         if (saaz2 .gt. 360.) saaz2 = saaz2 - 360.
         
         ! Open azimuth based on defining phases for cluster vector
         do i = 1,nst(iev) ! Index table on azimuth for sorting
            azesiev(i) = azes(iev,i)
         end do
         call indexx (nst(iev), azesiev, indx)
         openazc = 0.
         do ii = 1,nst(iev)-1
            i = indx(ii)
            ilast = i
            if (connected(iev,i)) then
               do jj = ii+1,nst(iev)
                  j = indx(jj)
                  if (connected(iev,j)) then
                     azdif = azesiev(j) - azesiev(i)
                     openazc = max(openazc,azdif)
                     ilast = j
                     exit
                  end if
               end do
            end if
         end do
         ! Gap between the largest and smallest azimuth
         do jj = 1,nst(iev)-1
            j = indx(jj)
            if (connected(iev,j)) then
               azdif = azesiev(j) + 360. - azesiev(ilast)
               openazc = max(openazc,azdif)
               exit
            end if
         end do
         if (openazc .ge. openaz_limit) then
            write (msg,'(a,i3,1x,a,1x,f5.1)') 'mlocout_hdf: large open azimuth, event ', iev, trim(evtnam(iev)), openazc
            call fyi (trim(msg))
         end if
         
         ! Closest and furthest station, based on defining stations for cluster vector
         rstadelc = 180.
         cstadelc = 0.
         do i = 1,nst(iev)
            if (connected(iev,i)) then
               rstadelc = min(rstadelc,delt(iev,i))
               cstadelc = max(cstadelc,delt(iev,i))
            end if
         end do
         
         ! Number of readings flagged as outliers (fcode = 'x')
         ndatx = 0
         do i = 1,nst(iev)
            if (fcode(iev,i) .eq. 'x') ndatx = ndatx + 1
         end do
         
         lon_test = lonp(iev,it1)
         call set_longitude_range (lon_test, longitude_range)
         
         ! Log the parameters for the cal_ command in case we want to use these events later as calibration events
         if (calibration_type .gt. 0) then
            write (io_log,'(i3,1x,i4,2i3,1x,a4,2i3,f6.2,f10.5,f11.5,f7.2,i4,2f6.2,f6.1,f6.2)')&
             iev,& ! Event number
             iyre(iev),& ! Origin year.
             mone(iev),& ! Origin month.
             idye(iev),& ! Origin day.
             'cal'//cal2,& ! cal_ command
             hourp(iev,it1),& ! Origin hour.
             minp(iev,it1),& ! Origin minute.
             secp(iev,it1),& ! Origin seconds.
             latp(iev,it1),& ! Geographic latitude.
             lon_test,& ! Geographic longitude.
             depthp(iev,it1),& ! Final depth (km).
             nint(saaz1),& ! Semi-axis azimuth.
             xl1,& ! Semi-axis length, km.
             xl2,& ! Semi-axis length, km.
             ddep,& ! Uncertainty in depth (km)
             dot ! Uncertainty in origin time (sec).
         end if
         
         ! Simplified hypocentral data written to the log file for easy import into a document.
         ! To avoid interference with the 'cal_' command output, it is written to a character array
         ! and printed later
         write (hypocenter_list(iev),'(i3,1x,i4,4i3,f6.2,f9.4,f10.4,f6.1,1x,a4,1x,a3,a2)')&
          iev,& ! Event number
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          hourp(iev,it1),& ! Origin hour.
          minp(iev,it1),& ! Origin minute.
          secp(iev,it1),& ! Origin seconds.
          latp(iev,it1),& ! Geographic latitude.
          lon_test,& ! Geographic longitude.
          depthp(iev,it1),& ! Final depth (km).
          cal_code,& ! Calibration code (GTCNU)
          rmag_pr,& ! magnitude.
          mmag(iev) ! magnitude scale
   
         ! HDF file for uncalibrated or direct calibration      
         write (hdfline(iev),fmt)&
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          hourp(iev,it1),& ! Origin hour.
          minp(iev,it1),& ! Origin minute.
          secp(iev,it1),& ! Origin seconds.
          latp(iev,it1),& ! Geographic latitude.
          lon_test,& ! Geographic longitude.
          depthp(iev,it1),& ! Final depth (km).
          depset_pr(iev),& ! How depth was set
          free_depth,& ! Free depth flag
          depth_inp(iev,1),& ! Depth from input file.
          rmag_pr,& ! magnitude.
          mmag(iev),& ! magnitude scale
          evid(iev),& ! Event ID, right-most 10 characters
          ndat(iev,it),& ! Number of observations contributed to hypocentroid estimation.
          ndatc(iev,it),& ! Number of observations used for cluster vector.
          ndatx,& ! Number of observations flagged as outliers (fcode = 'x')
          shatsqci(iev),& ! Normalized sample variance for the cluster vector 
          dot,& ! Uncertainty in origin time (sec). Formerly SE of position.
          fdp,& ! + uncertainty in depth (deeper), in km.
          fdm,& ! - uncertainty in depth (shallower), in km.
          rstadelc,& ! Epicentral distance of nearest station for cluster vector
          cstadelc,& ! Epicentral distance of farthest station for cluster vector
          openazc,& ! Largest open azimuth for cluster vector.
          nint(saaz1),& ! Semi-axis azimuth.
          xl1,& ! Semi-axis length, km.
          nint(saaz2),& ! Semi-axis azimuth.
          xl2,& ! Semi-axis length, km.
          avh,& ! Area of confidence ellipse, km**2.
          cal_code ! Calibration code (GTCNU)
         write (io_out,'(2a)') trim(hdfline(iev)), ' '//annotation(iev)
      end do
      
      write (io_log,'(/a)') 'Basic hypocenter list, uncalibrated or direct calibration:'
      do iev = 1,n_event
         write (io_log,'(a)') hypocenter_list(iev)
      end do
      
      close (io_out)
      
      ! Put a copy of the ~.hdf_dcal file in the ~_gccel/ subdirectory, but not if indirect
      ! calibration is also being done
      if (gccelout .and. .not.indirect_cal) then
         command_line = 'cp '//trim(outfil)//' '//trim(gcat_folder)//dirsym
         call execute_command_line (trim(command_line))
      end if
   
      ! HDF file for indirect calibrated locations
      if (indirect_cal) then
         hypocenter_list = ' '
         outfil=trim(outfile)//'.hdf_cal'
         call open_file (io_out, outfil, 'new')
         write (io_log,'(/a)') 'From Indirect calibration:' ! Logging cal_ commands
   
         do iev = 1,n_event ! Loop over events
         
            ! Calibration code
            call calibration_code (iev, 2, cal_code)
            cal2 = cal_code(2:2)
            call uctolc (cal2,-1)
            
            ! Uncertainty of focal depth
            fdp = ' '
            fdm = ' '
            if (mindx(iev,3) .gt. 0) then ! Free depth solution
               free_depth = 'f'
               if (ddep .le. 99.) then
                  write (fdp,'(f4.1)') max(min_depth_uncertainty,sqrt(accv(iev,3,3))) ! Uncertainty of depth, calibration-shifted
                  write (fdm,'(f4.1)') max(min_depth_uncertainty,sqrt(accv(iev,3,3))) ! Uncertainty of depth, calibration-shifted
               end if
            else ! Take from default or assigned uncertainties
               free_depth = ' '
               if (depthp_plus(iev) .le. 99.) write (fdp,'(f4.1)') max(min_depth_uncertainty,depthp_plus(iev))
               if (depthp_minus(iev) .le. 99.) write (fdm,'(f4.1)') max(min_depth_uncertainty,depthp_minus(iev))
            end if
            
            ! Magnitude
            if (rmag(iev) .gt. 0.) then
               write (rmag_pr,'(f3.1)') rmag(iev)
            else
               rmag_pr = '   '
            end if
            
            avh = xl1cg(iev)*xl2cg(iev)*pi ! Geometric area of 90% confidence ellipse (including calibration shift uncertainty)
            call timecr (otsp_cal(iev), hourpr, minpr, secpr)
            
            ! Make sure the semi-axis azimuths are between 0 and 360 degrees. 
            saaz1 = alphacg(iev)
            if (saaz1 .lt. 0.) saaz1 = saaz1 + 360.
            saaz2 = saaz1 + 90.
            if (saaz2 .gt. 360.) saaz2 = saaz2 - 360.
               
            ! Open azimuth based on defining phases for cluster vector
            do i = 1,nst(iev) ! Index table on azimuth for sorting
               azesiev(i) = azes(iev,i)
            end do
            call indexx (nst(iev), azesiev, indx)
            openazc = 0.
            do ii = 1,nst(iev)-1
               i = indx(ii)
               ilast = i
               if (connected(iev,i)) then
                  do jj = ii+1,nst(iev)
                     j = indx(jj)
                     if (connected(iev,j)) then
                        azdif = azesiev(j) - azesiev(i)
                        openazc = max(openazc,azdif)
                        ilast = j
                        exit
                     end if
                  end do
               end if
            end do
            ! Gap between the largest and smallest azimuth
            do jj = 1,nst(iev)-1
               j = indx(jj)
               if (connected(iev,j)) then
                  azdif = azesiev(j) + 360. - azesiev(ilast)
                  openazc = max(openazc,azdif)
                  exit
               end if
            end do
   !         if (openazc .ge. openaz_limit) then
   !            write (msg,'(a,i3,1x,a,1x,f5.1)') 'mlocout_hdf: large open azimuth, event ', iev, trim(evtnam(iev)), openazc
   !            call fyi (trim(msg))
   !         end if
            
            ! Closest and furthest station, based on defining stations for cluster vector
            rstadelc = 180.
            cstadelc = 0.
            do i = 1,nst(iev)
               if (connected(iev,i)) then
                  rstadelc = min(rstadelc,delt(iev,i))
                  cstadelc = max(cstadelc,delt(iev,i))
               end if
            end do
            
            ! Number of readings flagged as outliers (fcode = 'x')
            ndatx = 0
            do i = 1,nst(iev)
               if (fcode(iev,i) .eq. 'x') ndatx = ndatx + 1
            end do
            
            lon_test = lonp_cal(iev)
            call set_longitude_range (lon_test, longitude_range)
            
            ! Log the parameters for the cal_ command in case we want to use these events later as calibration events
            write (io_log,'(i3,1x,i4,2i3,1x,a4,2i3,f6.2,f10.5,f11.5,f7.2,i4,2f6.2,f6.1,f6.2)')&
             iev,& ! Event number
             iyre(iev),& ! Origin year.
             mone(iev),& ! Origin month.
             idye(iev),& ! Origin day.
             'cal'//cal2,& ! cal_ command
             hourpr,& ! Origin hour.
             minpr,& ! Origin minute.
             secpr,& ! Origin seconds.
             latp_cal(iev),& ! Geographic latitude.
             lon_test,& ! Geographic longitude.
             depthp_cal(iev),& ! Final depth (km).
             nint(saaz1),& ! Semi-axis azimuth.
             xl1cg(iev),& ! Semi-axis length, km.
             xl2cg(iev),& ! Semi-axis length, km.
             ddep,& ! Uncertainty in depth (km)
             dot ! Uncertainty in origin time (sec).
   
            ! Simplified hypocentral data written to the log file for easy import into a document.
            ! To avoid interference with the 'cal_' command output, it is written to a character array
            ! and printed later
            write (hypocenter_list(iev),'(i3,1x,i4,4i3,f6.2,f9.4,f10.4,f6.1,1x,a4,1x,a3,a2)')&
             iev,& ! Event number
             iyre(iev),& ! Origin year.
             mone(iev),& ! Origin month.
             idye(iev),& ! Origin day.
             hourp(iev,it1),& ! Origin hour.
             minp(iev,it1),& ! Origin minute.
             secp(iev,it1),& ! Origin seconds.
             latp(iev,it1),& ! Geographic latitude.
             lon_test,& ! Geographic longitude.
             depthp(iev,it1),& ! Final depth (km).
             cal_code,& ! Calibration code (GTCNU)
             rmag_pr,& ! magnitude.
             mmag(iev) ! magnitude scale
   
            write (hdfline(iev),fmt)&
             iyre(iev),& ! Origin year.
             mone(iev),& ! Origin month.
             idye(iev),& ! Origin day.
             hourpr,& ! Origin hour.
             minpr,& ! Origin minute.
             secpr,& ! Origin seconds.
             latp_cal(iev),& ! Geographic latitude.
             lon_test,& ! Geographic longitude.
             depthp_cal(iev),& ! Final depth (km).
             depset_pr(iev),& ! How depth was set
             free_depth,& ! Free depth flag
             depth_inp(iev,1),& ! Depth from input file.
             rmag_pr,& ! magnitude.
             mmag(iev),& ! magnitude scale
             evid(iev),& ! Event ID, right-most 10 characters
             ndat(iev,it),& ! Number of observations contributed to hypocentroid estimation.
             ndatc(iev,it),& ! Number of observations used for cluster vector.
             ndatx,& ! Number of observations flagged as large outliers (fcode = 'x')
             shatsqci(iev),& ! Normalized sample variance for the cluster vector ! Adopted here on 1/19/2012
             sqrt(accv(iev,4,4)),& ! Uncertainty of OT, calibration shifted.
             fdp,& ! + uncertainty in depth (deeper), in km.
             fdm,& ! - uncertainty in depth (shallower), in km.
             rstadelc,& ! Epicentral distance of nearest station for cluster vector
             cstadelc,& ! Epicentral distance of farthest station for cluster vector
             openazc,& ! Largest open azimuth for cluster vector.
             nint(saaz1),& ! Semi-axis azimuth.
             xl1cg(iev),& ! Semi-axis length, km.
             nint(saaz2),& ! Semi-axis azimuth.
             xl2cg(iev),& ! Semi-axis length, km.
             avh,& ! Area of confidence ellipse, km**2.
             cal_code ! Calibration code (GTCNU)
            write (io_out,'(2a)') trim(hdfline(iev)), ' '//annotation(iev)
         end do ! End of loop over events
            
         write (io_log,'(/a)') 'Basic hypocenter list from indirect calibration:'
         do iev = 1,n_event
            write (io_log,'(a)') hypocenter_list(iev)
         end do
      
         close (io_out)
   
         ! Put a copy of the ~.hdf_cal file in the ~_gccel/ subdirectory
         if (gccelout) then
            command_line = 'cp '//trim(outfil)//' '//trim(gcat_folder)//dirsym
            call execute_command_line (trim(command_line))
         end if
   
      end if
         
      return
      
   end subroutine write_hdf
      

!*****************************************************************************************
   subroutine write_phase_data (it)
   
   ! Standard printed output files:
   !   Phase data (~.phase_data)
   !   Large absolute residuals (~.xdat)
   !   Large cluster residuals (~.lres)
   !   Depth phases (~.depth_phases)
   ! Optional:
   !   direct calibration hypocentroid data (~.phase_data)
   !   indirect calibration hypocentroid data after direct calibration (~.cal_phase_data)
   !   all shifted phase data after indirect calibration (~.phase_data_cal)
   !   Limited distance range (.oldr)
   
   ! Also fills the arrays to calculate the reading error for output to .rderr
   
   ! January 26, 1989 by eric bergman
   ! Modified 9/20/94 by eab.
   ! 11/17/2017: format change for v10.4.0, adding deployment code
   ! 3/28/2024: Major re-write, breaking off many repetitive parts into subroutines
   
      integer :: i
      integer :: iev
      integer :: ii
      integer, dimension(n_arrtim_max) :: indx
      integer :: io_hypo_cal
      integer :: io_hypo_dcal
      integer :: io_out_cal
      integer, intent(in) :: it ! iteration #
      integer :: it1
      integer :: j
!      integer :: lunit
      integer :: n_res
      real, dimension(n_arrtim_max) :: deltiev
      real :: deltkm
      real, dimension(0:n_it_max1) :: dts
      real :: dtsit
      real :: dtslim
      real :: h_mult
      real :: mean_res1
      real :: mean_res2
      real, parameter :: s_parent1 = 1.0 ! Standard deviation of parent population for zM test, based on raw residuals
      real, parameter :: s_parent2 = 1.0 ! Standard deviation of parent population for zM test, based on eci (cluster residuals)
      real :: sum_res1
      real :: sum_res2
      real :: z_stat1
      real :: z_stat2      
      character(len=1) :: bdp_flag
      character(len=120) :: bfmt
      character(len=60) :: depth_set
      character(len=120) :: gfmt
      character(len=21) :: qtest
      character(len=4) :: reason
      character(len=4), parameter :: rpres = 'PRES'
!      logical :: bdp
      
      it1 = it + 1
      if (nstep .eq. 0) it1 = 0 ! Forward modelling
      
      ! Initialize reading errors
      indexq = 0
      qres = 0.
   
      ! Open output files (~.lres is opened in the main program)
      call open_file (io_out, trim(outfile)//'.phase_data', 'new') ! Main phase data file
      call open_file (io_xdat, trim(outfile)//'.xdat', 'new') ! List of unflagged readings that fail PRES
      if (oldr_out) call open_file (io_oldr, trim(outfile)//'.oldr', 'new') ! Limited distance range
      if (direct_cal) then ! Direct calibration data for hypocentroid
         io_hypo_dcal = lunit()
         call open_file (io_hypo_dcal, trim(outfile)//'.dcal_phase_data', 'new')
         if (indirect_cal) then ! If direct calibration is followed by indirect calibration
            io_hypo_cal = lunit()
            call open_file (io_hypo_cal, trim(outfile)//'.cal_phase_data', 'new')
         end if
      end if
      if (indirect_cal) then
         io_out_cal = lunit()
         call open_file (io_out_cal, trim(outfile)//'.phase_data_cal', 'new') ! shifted phase data after indirect calibration
      end if
      
      ! Format assignments for individual phase data, good and bad
      call format_assignment (1, it, gfmt)
      call format_assignment (0, it, bfmt)
      
      ! Some information about event depths
      call depth_statistics (it1, io_depth_phase)
      
      ! Column labels
      if (direct_cal) then
         call write_column_labels (io_hypo_dcal, 1)
         if (indirect_cal) call write_column_labels (io_hypo_cal, 1)
      end if
      if (oldr_out) call write_column_labels (io_oldr, 1)
      call write_column_labels (io_depth_phase, 1)
      
      if (direct_cal) then
         write (io_log,'(/a)') 'Mean residual and mean ECI of readings used for direct calibration'
         write (io_log,'(a,t26,a,1x,a,1x,a,3x,a,6x,a,1x,a,3x,a)')&
          'Event', 'N', 'Mean Res', 'Sres', 'zres', 'Mean ECI', 'Seci', 'zeci'
         write (*,'(/a)')&
          'Possible need for depth adjustment based on zM tests (of raw residuals and ECI values):'
         write (*,'(a,t26,a,1x,a,1x,a,3x,a,6x,a,1x,a,3x,a)')&
          'Event', 'N', 'Mean Res', 'Sres', 'zres', 'Mean ECI', 'Seci', 'zeci'
      end if
   
      do iev = 1,n_event ! Loop over events for good data
            
         ! How depth was set
         call depth_set_text (iev, depth_set)
      
         ! Event separator
         write (io_out,'(1x,164a1)') ('*',i=1,164)
         if (indirect_cal) write (io_out_cal,'(1x,164a1)') ('*',i=1,164)
         
         ! Heading for each event
         write (io_out,'(a,i3,t30,a30,5x,a9,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev), 'GOOD DATA',' '
         if (indirect_cal) write (io_out_cal,'(a,i3,t30,a30,5x,a9,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev), 'GOOD DATA',' '
         call cluster_event_header (io_depth_phase, iev, depthp(iev,it1), depth_set, 1)
         if (direct_cal) then
            call cluster_event_header (io_hypo_dcal, iev, depthp(iev,it1), depth_set, 1)
            if (indirect_cal) call cluster_event_header (io_hypo_cal, iev, depthp_cal(iev), depth_set, 1)
         end if
         if (oldr_out) call cluster_event_header (io_oldr, iev, depthp(iev,it1), depth_set, 1)
   
         write (io_out,'(1x,a,t165,a)') trim(infile(iev)), ' ' ! Data file
         if (indirect_cal) write (io_out_cal,'(1x,a,t165,a)') trim(infile(iev)), ' '
         
         ! Starting and final hypocenters
         call start_final_location (it1, iev, io_out)
         ! Header line
         write (io_out,'(t165,a)') ' '
         call write_column_labels (io_out, 1)
         write (io_out,'(t165,a)') ' '
         if (indirect_cal) then
            call start_final_location (it1, iev, io_out_cal)
            write (io_out_cal,'(t165,a)') ' '
            call write_column_labels (io_out_cal, 1)
            write (io_out_cal,'(t165,a)') ' '
         end if
         
         ! Index table on delta for sorting
         do i=1,nst(iev)
            deltiev(i) = delt(iev,i)
         end do
         call indexx (nst(iev), deltiev, indx)
         
         h_mult = depthp(iev,it1)*3. ! Multiple of focal depth, used to select near-source readings for depth control
         
         sum_res1 = 0. ! Sum of residuals for readings used in direct calibration of the hypocentroid
         sum_res2 = 0. ! Sum of eci values for readings used in direct calibration of the hypocentroid
         n_res = 0
         
         do ii = 1,nst(iev) ! Loop over readings
         
            i = indx(ii)
            deltkm = delt(iev,i)*111. ! Used to select near-source readings for depth control 
            
            if (.not.fltrhres(iev,i) .and.&
                .not.fltrhflag(iev,i) .and.&
                .not.fltrcres(iev,i) .and.&
                .not.fltrcflag(iev,i)) then ! Not filtered
               
               do j = 0,it
                  dts(j) = dt(iev,i,j) - s(iev,i,j)
               end do
               
               call write_good_phase_line (io_out, it, iev, i, gfmt, ' ', .false., .true.)
               if (indirect_cal) call write_good_phase_line (io_out_cal, it, iev, i, gfmt, ' ', .true., .true.)
               
               ! Near-source readings written to the depth_phases file
               if (deltkm .le. max(h_mult,100.)) call write_good_phase_line (io_depth_phase, it, iev, i, gfmt, ' ', .false., .true.)
                  
               ! Depth phases, including relative depth phases, written to the depth_phases file                   
               if (phase(iev,i) .eq. 'pP      ' .or. phase(iev,i) .eq. 'sP      ' .or. phase(iev,i) .eq. 'pwP     ' .or.&
                   phase(iev,i) .eq. 'pP-P    ' .or. phase(iev,i) .eq. 'sP-P    ' .or. phase(iev,i) .eq. 'pwP-P   ') then
                  if (bdp(iev,i)) then ! Check for stations known to report bad depth phases
                     bdp_flag = '*'
                  else
                     bdp_flag = ' '
                  end if
                  call write_good_phase_line (io_depth_phase, it, iev, i, gfmt, bdp_flag, .false., .false.)
               end if
                
               ! Direct calibration readings written to dcal_phase_data file
               ! Increment sum of residuals.
               if (direct_cal) then
                  if (.not.fltrh(iev,i)) then
                     call write_good_phase_line (io_hypo_dcal, it, iev, i, gfmt, ' ', .false., .true.) 
                     sum_res1 = sum_res1 + dts(it)/sdread(iev,i)
                     sum_res2 = sum_res2 + eci(iev,i)
                     n_res = n_res + 1
                     if (indirect_cal) call write_good_phase_line (io_hypo_cal, it, iev, i, gfmt, ' ', .true., .true.)
                  end if
               end if
               
               ! Output for limited distance range file
               if (oldr_out .and. delt(iev,i) .ge. dist1 .and. delt(iev,i) .le. dist2) then
                  call write_good_phase_line (io_oldr, it, iev, i, gfmt, ' ', .false., .true.)
               end if
                  
               ! Output of readings with large cluster residuals
               if (lresout .and. it .gt. 0) then
                  if (abs(eci(iev,i)) .gt. lres) then
                     write (io_lres,'(i3,1x,a20,1x,a8,1x,a8,i5,1x,f8.2,1x,a,2x,a,1x,i6)') iev, sad(iev,i),&
                      readsrc(iev,i), phase0 (iev,i), mnf_line(iev,i), eci(iev,i), infile20(iev), phase(iev,i),&
                      diff_line(iev,i)
                  end if
               end if
            end if
               
            ! Fill in array of residuals for calculation of reading errors
            if (connected(iev,i)) then
               qtest = stname(iev,i)//deployment(iev,i)//phase(iev,i)
               do j = 1,nqc
                  if (qtest .eq. qname1(j) .and. idiff(iev,i) .eq. 0 .and. indexq(j) .lt. n_qres_max) then
                     indexq(j) = indexq(j) + 1
                     if (indirect_cal) then ! Indirect calibration
                        qres(j,indexq(j)) = dt_cal(iev,i) - s_cal(iev,i)
   !                        print *, iev, i, qtest, j, indexq(j), qres(j,indexq(j)), dt_cal(iev,i), s_cal(iev,i)
                     else ! Direct calibration or uncalibrated
                        qres(j,indexq(j)) = dt(iev,i,it) - s(iev,i,it)
                     end if
                     iqiev(j,indexq(j)) = iev ! array of event numbers association with residuals, needed for tomo3
                     exit
                  end if
               end do
            end if
        
         end do ! End of loop over readings
         write (io_out,'(t165,a)') ' '
         write (io_out,'(t165,a)') ' '
         if (indirect_cal) then
            write (io_out_cal,'(t165,a)') ' '
            write (io_out_cal,'(t165,a)') ' '
         end if      
         ! zM test statistic for depth resolution
         if (n_res .gt. 0) then
            mean_res1 = sum_res1 / real(n_res)
            mean_res2 = sum_res2 / real(n_res)
            z_stat1 = sqrt(real(n_res)) * abs(mean_res1) / s_parent1
            z_stat2 = sqrt(real(n_res)) * abs(mean_res2) / s_parent2
            write (io_log,'(2x,i3,1x,a16,1x,i3,1x,f8.4,1x,f4.2,1x,f6.2,6x,f8.4,1x,f4.2,1x,f6.2)')&
             iev, evtnam(iev)(1:16), n_res, mean_res1, s_parent1, z_stat1, mean_res2, s_parent2, z_stat2
            if (z_stat1 .gt. 2.0 .or. z_stat2 .gt. 2.0) then
               if (z_stat1 .gt. 2.0 .and. z_stat2 .gt. 2.0) then
                  write (*,'(2x,i3,1x,a16,1x,i3,1x,f8.4,1x,f4.2,1x,f6.2,6x,f8.4,1x,f4.2,1x,f6.2,a)')&
                   iev, evtnam(iev)(1:16), n_res, mean_res1, s_parent1, z_stat1, mean_res2, s_parent2, z_stat2, ' *'
               else
                  write (*,'(2x,i3,1x,a16,1x,i3,1x,f8.4,1x,f4.2,1x,f6.2,6x,f8.4,1x,f4.2,1x,f6.2)')&
                   iev, evtnam(iev)(1:16), n_res, mean_res1, s_parent1, z_stat1, mean_res2, s_parent2, z_stat2
               end if
            end if
         end if
         
      end do ! End of loop over events for good data
      
      !*********************************
      ! Bad phase data for each event
      
      if (debug) print *, 'Bad data'
      do iev = 1,n_event ! Loop over events
      
         ! How depth was set
         if (debug) write (io_log,'(i4,a)') iev, ' How depth was set'
         call depth_set_text (iev, depth_set)
      
         ! Event separator
         write (io_out,'(1x,164a1)') ('*',i=1,164)
         if (indirect_cal) write (io_out_cal,'(1x,164a1)') ('*',i=1,164)
         
         ! Heading for each event
         if (debug) write (io_log,'(i4,a)') iev, ' Heading for each event'
         write (io_out,'(a,i3,t30,a30,5x,a8,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev), 'BAD DATA',' '
         if (indirect_cal) write (io_out_cal,'(a,i3,t30,a30,5x,a8,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev), 'BAD DATA',' '
         call cluster_event_header (io_depth_phase, iev, depthp(iev,it1), depth_set, 0)
         if (direct_cal) then
            call cluster_event_header (io_hypo_dcal, iev, depthp(iev,it1), depth_set, 0)
            if (indirect_cal) call cluster_event_header (io_hypo_cal, iev, depthp_cal(iev), depth_set, 0)
         end if
         if (oldr_out) call cluster_event_header (io_oldr, iev, depthp(iev,it1), depth_set, 0)
   
         write (io_out,'(1x,a,t165,a)') trim(infile(iev)), ' ' ! Data file
         if (indirect_cal) write (io_out_cal,'(1x,a,t165,a)') trim(infile(iev)), ' '
         ! Starting and final hypocenters
         call start_final_location (it1, iev, io_out)
         if (indirect_cal) call start_final_location (it1, iev, io_out_cal)
   
         ! Missing stations
         if (debug) write (io_log,'(2i4,a)') iev, nmiss(iev), ' Missing stations'
         if (nmiss(iev) .gt. 0) then
            write (io_out,'(t165,a)') ' '
            write (io_out,'(a,t165,a)') ' Missing stations:', ' '
            do i = 1,nmiss(iev) 
               write (io_out,'(5x,a,t165,a)') missta(iev,i), ' '
            end do
            if (indirect_cal) then
               write (io_out_cal,'(t165,a)') ' '
               write (io_out_cal,'(a,t165,a)') ' Missing stations:', ' '
               do i = 1,nmiss(iev) 
                  write (io_out_cal,'(5x,a,t165,a)') missta(iev,i), ' '
               end do
            end if
         end if
         
         ! Header line 
         write (io_out,'(t165,a)') ' '
         call write_column_labels (io_out, 0)
         write (io_out,'(t165,a)') ' '
         if (indirect_cal) then
            write (io_out_cal,'(t165,a)') ' '
            call write_column_labels (io_out_cal, 0)
            write (io_out_cal,'(t165,a)') ' '
         end if
         
         ! Index table on delta for sorting
         if (debug) write (io_log,'(a)') 'Index table on delta for sorting'
         do i = 1,nst(iev)
            deltiev(i) = delt(iev,i)
         end do
         call indexx (nst(iev), deltiev, indx)
         
         do ii=1,nst(iev) ! Loop over readings
            i = indx(ii)
            
            if (fcode(iev,i) .eq. 'm') cycle ! Don't print readings from missing stations
            if (fcode(iev,i) .eq. 'd') cycle ! Don't print duplicate readings
            
            ! Reasons for badness
            ! Only consider if a reading has been flagged or has a large residual
            if (fltrhres(iev,i) .or. fltrhflag(iev,i) .or. fltrcres(iev,i) .or. fltrcflag(iev,i)) then
               do j = 0,it
                  dts(j) = dt(iev,i,j) - s(iev,i,j)
               end do
               
               if (fltrhres(iev,i) .or. fltrcres(iev,i)) then ! Large residual
                  reason = rpres
                  if (.not.fltrhflag(iev,i) .and. .not.fltrcflag(iev,i)) then ! Not already flagged
                     write (io_xdat,'(i3,1x,a6,1x,a8,1x,a8,1x,i5,2(1x,a))') iev, stname(iev,i),&
                      readsrc(iev,i), phase0(iev,i), mnf_line(iev,i), 'x', infile20(iev)
                  end if
               end if
               
               if (fltrhflag(iev,i) .or. fltrcflag(iev,i)) then ! Flagged for some reason
                  dtslim = 2.0*ttsprd(iev,i)
                  dtsit = dt(iev,i,it) - s(iev,i,it)
                  if (abs(dtsit-ttoff(iev,i)) .le. dtslim .and. (fcode(iev,i) .eq. 'x' .or. fcode(iev,i) .eq. 's')) then
                     reason = '?   ' ! Maybe investigate why small residual is flagged
                  else
                     reason = '    ' ! Leave blank for flagged phases
                  end if
               end if
               
               call write_bad_phase_line (io_out, it, iev, i, bfmt, ' ', .false., .true., reason)
               if (indirect_cal) call write_bad_phase_line (io_out_cal, it, iev, i, bfmt, ' ', .true., .true., reason)
                
               ! Near-source readings written to the depth_phases file
               if (deltkm .le. max(h_mult,100.)) then
                  call write_bad_phase_line (io_depth_phase, it, iev, i, bfmt, ' ', .false., .true., reason)
               end if
                
               ! Depth phases written to the depth_phases file
               if (phase(iev,i) .eq. 'pP      ' .or. phase(iev,i) .eq. 'sP      ' .or. phase(iev,i) .eq. 'pwP     ' .or.&
                   phase(iev,i) .eq. 'pP-P    ' .or. phase(iev,i) .eq. 'sP-P    ' .or. phase(iev,i) .eq. 'pwP-P   ') then
                  if (bdp(iev,i)) then ! Check for stations known to report bad depth phases
                     bdp_flag = '*'
                  else
                     bdp_flag = ' '
                  end if
                  call write_bad_phase_line (io_depth_phase, it, iev, i, bfmt, bdp_flag, .false., .false., reason)
               end if
                
               ! Direct calibration readings written to dcal_phase_data file
               if (direct_cal) then
                  if (.not.fltrhdelt(iev,i) .and. (fltrhres(iev,i) .or. fltrhflag(iev,i)))&
                   call write_bad_phase_line (io_hypo_dcal, it, iev, i, bfmt, ' ', .false., .true., reason)
                  if (indirect_cal) then
                     if (.not.fltrhdelt(iev,i) .and. (fltrhres(iev,i) .or. fltrhflag(iev,i)))&
                      call write_bad_phase_line (io_hypo_cal, it, iev, i, bfmt, ' ', .true., .true., reason)
                  end if
               end if
               
               ! Output for limited distance range file
               if (oldr_out .and. delt(iev,i) .ge. dist1 .and. delt(iev,i) .le. dist2) then
                  call write_bad_phase_line (io_oldr, it, iev, i, bfmt, ' ', .false., .true., reason)
               end if
               
            end if
            
         end do ! End of loop over readings
         
         write (io_out,'(t165,a)') ' '
         write (io_out,'(t165,a)') ' '
         if (indirect_cal) then
            write (io_out_cal,'(t165,a)') ' '
            write (io_out_cal,'(t165,a)') ' '
         end if
         
      end do ! End of loop over events
      
      close (io_out)
      close (io_xdat)
      close (io_depth_phase)
      if (indirect_cal) close (io_out_cal)
      if (lresout) close (io_lres)
      if (direct_cal) then
         close (io_hypo_dcal)
         if (indirect_cal) close (io_hypo_cal)
      end if
      if (oldr_out) close (io_oldr)
      
      return
      
   end subroutine write_phase_data
   
   
!*****************************************************************************************
   subroutine format_assignment (ifmt, it, fmt)
   
   ! Assigns a format for good or bad data, depending on iteration number
   ! ifmt = 1 for good data
   ! ifmt = 0 for bad data
   
      integer, intent(in) :: ifmt
      integer, intent(in) :: it
      character(len=120), intent(out) :: fmt
      character(len=132) :: msg
      character(len=80) :: part0g = '                ,t134,a8,1x,a3,1x,a8,1x,2i5)'
      character(len=80) :: part1 = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,'
      character(len=80) :: part2b = ',t112,a4        ,t134,a8,1x,a3,1x,a8,1x,2i5)'
      character(len=80) :: part2g = ',t112,2f7.4,f7.2,t134,a8,1x,a3,1x,a8,1x,2i5)'
      
      if (ifmt .eq. 1) then ! good data
         select case (it)
            case (0)
               fmt = trim(part1)//'1f7.2'//trim(part0g)
            case (1)
               fmt = trim(part1)//'2f7.2'//trim(part2g)
            case (2)
               fmt = trim(part1)//'3f7.2'//trim(part2g)
            case (3)
               fmt = trim(part1)//'4f7.2'//trim(part2g)
            case (4)
               fmt = trim(part1)//'5f7.2'//trim(part2g)
            case default
               write (msg,'(a,i3)') 'format_assignment: illegal value for "it": ', it
               call oops (trim(msg))
         end select
      end if 
   
      if (ifmt .eq. 0) then ! bad data
         select case (it)
            case (0)
               fmt = trim(part1)//'1f7.2'//trim(part2b)
            case (1)
               fmt = trim(part1)//'2f7.2'//trim(part2b)
            case (2)
               fmt = trim(part1)//'3f7.2'//trim(part2b)
            case (3)
               fmt = trim(part1)//'4f7.2'//trim(part2b)
            case (4)
               fmt = trim(part1)//'5f7.2'//trim(part2b)
            case default
               write (msg,'(a,i3)') 'format_assignment: illegal value for "it": ', it
               call oops (trim(msg))
         end select
      end if 
   
      return
      
   end subroutine format_assignment
   
   
!*****************************************************************************************
   subroutine write_column_labels (io_unit, itype)
   
   ! Write the column labels for a "phase_data" file to a specified logical unit
   ! Argument "itype" controls the heading for good data or bad data
   
      integer, intent(in) :: io_unit
      integer, intent(in) :: itype
      character(len=120) :: line1
      character(len=120) :: line1a
      character(len=120) :: line2
      character(len=120) :: line2a
      character(len=120) :: line2b
      character(len=120) :: line3a
      character(len=120) :: line3b
      character(len=132) :: msg
      logical :: op
      
      data line1  /' STA   NETWORK    PHASE    READ  DELTA AZIM  RAY     WGT   STA  P RESIDUALS FOR ITERATION #'/
      data line1a /'                    WHY'/
      data line2  /' CODE                       ERR             PARAM         CORR  *INP*   *0*    *1*    *2*    *3*    *4*'/
      data line2a /'          DTMPH  DTMPC    ECI'/
      data line2b /'        BAD'/
      data line3a /' AUTHOR   CHA PHASE0     MNF IDIF'/
      data line3b /'                   AUTHOR   CHA PHASE0     MNF IDIF'/
      
      inquire (unit=io_unit,opened=op)
      if (.not.op) then
         write (msg,'(a,i3,a)') 'write_column_labels: unit ', io_unit, ' is not open'
         call oops (trim(msg))
      end if
      
      if (itype .eq. 1) then
         write (io_unit,'(a,t165,a)') trim(line1), ' '
         write (io_unit,'(3a)') trim(line2),trim(line2a),trim(line3a)
      else
         write (io_unit,'(2a,t165,a)') trim(line1), trim(line1a), ' '
         write (io_unit,'(3a)') trim(line2),trim(line2b),trim(line3b)
      end if
      
      return
      
   end subroutine write_column_labels
   
   
!*****************************************************************************************
   subroutine depth_set_text (iev, depth_set)
   
   ! Information on how depths were set
   
      ! Arguments
      integer, intent(in) :: iev
      character(len=60), intent(out) :: depth_set
      
      character(len=8) :: depth_free_fixed
      character(len=60) :: depth_set1
      character(len=60) :: depth_set2
   
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
      if (mindx(iev,3) .gt. 0) then
         depth_free_fixed = ', free'
         depth_set = trim(depth_free_fixed)
      else
         depth_free_fixed = ', fixed'
         depth_set = trim(depth_free_fixed)//trim(depth_set1)//trim(depth_set2)
      end if
   
      return
      
   end subroutine depth_set_text
   
   
!*****************************************************************************************
   subroutine cluster_event_header (io_unit, iev, depth_value, depth_set, itype)
   
   ! Header for each event's data written to a specified logical unit
   
      ! Arguments
      integer, intent(in) :: io_unit ! logical unit
      integer, intent(in) :: iev ! event #
      real, intent(in) :: depth_value ! focal depth
      character(len=60), intent(in) :: depth_set ! how focal depth was set
      integer, intent(in) :: itype ! 1 = good data, otherwise bad data
      
      character(len=48) :: fmt
      
      data fmt /'(/a,i3,t30,a30,t65,a,t80,a,f5.1,a,t165,a)'/
      
      if (itype .eq. 1) then
         write (io_unit,fmt) ' CLUSTER EVENT ', iev, evtnam(iev), 'GOOD DATA','Depth = ', depth_value, depth_set, ' '
      else
         write (io_unit,fmt) ' CLUSTER EVENT ', iev, evtnam(iev), 'BAD DATA','Depth = ', depth_value, depth_set, ' '
      end if
         
      return
      
   end subroutine cluster_event_header
      
   
!*****************************************************************************************
   subroutine write_good_phase_line (io_unit, it, iev, i, fmt, c1, shifted, path_corrected)
   
   ! Write a phase line of good data to the specified logical unit. Argument 'i' is an index
   ! to the reading in the event data file, sorted by epicentral distance. 
   
      integer, intent(in) :: i
      integer, intent(in) :: iev
      integer, intent(in) :: io_unit
      integer, intent(in) :: it
      integer :: j
      real :: azes_pr
      real :: delt_pr
      real, dimension(0:n_it_max1) :: dtj
      real :: stacs_pr
      character(len=1), intent(in) :: c1
      character(len=120), intent(in) :: fmt
      logical, intent(in) :: path_corrected
      logical, intent(in) :: shifted
      
      do j = 0,it
         if (shifted) then
            dtj(j) = dt_cal(iev,i) - s_cal(iev,i)
         else
            if (path_corrected) then
               dtj(j) = dt(iev,i,j) - s(iev,i,j)
            else
               dtj(j) = dt(iev,i,j)
            end if
         end if
      end do
      
      if (shifted) then
         delt_pr = delt_cal(iev,i)
         azes_pr = azes_cal(iev,i)
         stacs_pr = stacs_cal(iev,i)
      else
         delt_pr = delt(iev,i)
         azes_pr = azes(iev,i)
         stacs_pr = stacs(iev,i)
      end if
      
      if (it .eq. 0) then ! Forward problem
         write (io_unit,fmt) c1, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
          phase(iev,i), sdread(iev,i), delt_pr, nint(azes_pr), psd(iev,i), weight(iev,i),&
          stacs_pr, resisc(iev,i), dtj(0), readsrc(iev,i), channel(iev,i), phase0(iev,i),&
          mnf_line(iev,i), diff_line(iev,i)
      else ! Iterated results
         write (io_unit,fmt) c1, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
          phase(iev,i), sdread(iev,i), delt_pr, nint(azes_pr), psd(iev,i), weight(iev,i),&
          stacs_pr, resisc(iev,i), (dtj(j), j=0,it), dtmph(iev,i), dtmpc(iev,i), eci(iev,i),&
          readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
      end if
      
      return
      
   end subroutine write_good_phase_line
   
   
!*****************************************************************************************
   subroutine write_bad_phase_line (io_unit, it, iev, i, fmt, c1, shifted, path_corrected, reason)
   
   ! Write a phase line of bad data to the specified logical unit. Argument 'i' is an index
   ! to the reading in the event data file, sorted by epicentral distance. 
   
      integer, intent(in) :: i
      integer, intent(in) :: iev
      integer, intent(in) :: io_unit
      integer, intent(in) :: it
      integer :: j
      real :: azes_pr
      real :: delt_pr
      real :: dtj(0:n_it_max1)
      real :: stacs_pr
      character(len=1), intent(in) :: c1
      character(len=120), intent(in) :: fmt
      character(len=4), intent(in) :: reason
      logical, intent(in) :: path_corrected
      logical, intent(in) :: shifted
      
      do j = 0,it
         if (shifted) then
            dtj(j) = dt_cal(iev,i) - s_cal(iev,i)
         else
            if (path_corrected) then
               dtj(j) = dt(iev,i,j) - s(iev,i,j)
            else
               dtj(j) = dt(iev,i,j)
            end if
         end if
      end do
      
      if (shifted) then
         delt_pr = delt_cal(iev,i)
         azes_pr = azes_cal(iev,i)
         stacs_pr = stacs_cal(iev,i)
      else
         delt_pr = delt(iev,i)
         azes_pr = azes(iev,i)
         stacs_pr = stacs(iev,i)
      end if
   
      write (io_unit,fmt) c1, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
       phase(iev,i), sdread(iev,i), delt_pr, nint(azes_pr), psd(iev,i), weight(iev,i),&
       stacs_pr, resisc(iev,i), (dtj(j), j=0,it), reason,&
       readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
      
      return
      
   end subroutine write_bad_phase_line
   
   
!*****************************************************************************************
   subroutine start_final_location (it, iev, io_unit)
   
      ! Arguments
      integer, intent(in) :: iev ! event #
      integer, intent(in) :: io_unit ! logical unit
      integer, intent(in) :: it ! iteration #
      
      real :: lon_out
      character(len=60) :: fmt
      
      data fmt /'(1x,a6,i4,4(i3),f5.1,f8.3,f9.3,f6.1,f4.1,t165,a)'/
      
      ! Start location
      lon_out = lon_inp(iev)
      call set_longitude_range (lon_out, longitude_range)
      write (io_unit,fmt) 'Input ', iyre (iev), mone(iev), idye(iev), hour_inp(iev), min_inp(iev),&
       sec_inp(iev), lat_inp(iev), lon_out, depth_inp(iev,1), rmag(iev), ' '
      
      ! Final location, unshifted
      lon_out = lonp(iev,it)
      call set_longitude_range (lon_out, longitude_range)
      write (io_unit,fmt) 'Final ', iyre (iev), mone(iev), idye(iev), hourp(iev,it), minp(iev,it),&
       secp(iev,it), latp(iev,it), lon_out, depthp(iev,it), rmag(iev), ' '
         
   
      return
      
   end subroutine start_final_location
   
   
!*****************************************************************************************
   subroutine depth_statistics (it, io_unit)
   
   ! Some information about event depths
   
      ! Arguments
      integer, intent(in) :: io_unit ! logical unit
      integer, intent(in) :: it ! iteration number
      
      integer :: iev
      integer :: n_constrained
!       integer :: n_input
      real :: cluster_default_depth
      real, dimension(n_event) :: x_dep_constrained
      real, dimension(n_event) :: x_dep_input
      real :: xmed
      
      n_constrained = 0
!       n_input = 0
      cluster_default_depth = -99.
      do iev = 1, n_event
         if (depset_pr(iev) .ne. 'c' .and.&
             depset_pr(iev) .ne. 'u' .and.&
             depset_pr(iev) .ne. 'i' .and.&
             depset_pr(iev) .ne. ' ') then
            n_constrained = n_constrained + 1
            x_dep_constrained(n_constrained) = depthp(iev,it)
         end if
         if (depset_pr(iev) .eq. 'c') then
            cluster_default_depth = depthp(iev,it)
         end if
         x_dep_input(iev) = depth_inp(iev,1)
      end do
      if (cluster_default_depth .gt. 0.) write (io_unit,'(a,f6.1)') 'Cluster default depth: ', cluster_default_depth
      if (n_constrained .ge. 3) then
         call mdian1 (x_dep_constrained, n_constrained, xmed)
         write (io_unit,'(a,f6.1)') 'Median of constrained depths: ', xmed
         median_constrained_depths = xmed
      end if
      if (n_event .ge. 3) then
         call mdian1 (x_dep_input, n_event, xmed)
         write (io_unit,'(a,f6.1)') 'Median of input file depths: ', xmed
      end if
      
      return
      
   end subroutine depth_statistics


!*****************************************************************************************
subroutine write_rderr ()

! Summary file of updated reading errors.
! The measure of spread is 'sn', a robust, efficient scale estimator (Croux & Rousseeuw, 1992)
! The average is also estimated, for empirical path anomalies (but these are uncalibrated).
! The last column is the reading error actually used in this run.
! This is only done for cluster mode runs.
! 11/17/2017: format changed to add deployment code
      
   integer :: i
   integer :: j
   integer :: k
   real :: azeqst
   real :: delt1
   real :: dum1
   real :: dum2
   real :: dum3
   real :: dum5
   real :: dum6
   real :: lath_epa_g
   real :: lonh_epa_g
   real :: qlatdg
   real :: qlondg
   real, dimension(n_qres_max) :: res_data
   real :: sumdata
   real :: t1
   real :: t2
   real :: t3
   real :: t4
   character(len=132) :: msg
   character(len=100) :: outfil
   
   outfil=trim(outfile)//'.rderr'
   call open_file (io_out, outfil, 'new')
   
   ! First line carries the hypocentroid (calibrated if available) for plotting empirical path anomalies
   call geogra (lath_epa,lath_epa_g) ! hypocentroid in geographic coordinates
   lonh_epa_g = lonh_epa
   call set_longitude_range (lonh_epa_g, 0)
   write (io_out,'(a,3f10.3)') 'Hypocentroid: ', lath_epa_g, lonh_epa_g, depthh_epa

   do i = 1,nqc ! Loop over station-phases with data
      if (indexq(i) .ge. 2 .and. idiff0(i) .eq. 0) then
         sumdata = 0.
         k = 0
         do j = 1,indexq(i)
            ! On rare occasions, the calibration shift during indirect calibration causes the epicentral distance
            ! of a reading to fall off the phase branch and the residual becomes very large, upsetting the calculation of 
            ! mean residual (epa) and, if the sample size is small, the empirical reading error (ere). The test on
            ! absolute value of qres is meant to remove such readings from the calculation.
            if (abs(qres(i,j)) .lt. 20.) then
               k = k + 1
               res_data(k) = qres(i,j)
               sumdata = sumdata + res_data(k)
               if (debug) write (io_log,'(a,3x,i3,2f10.3)') qname1(i), j, qres(i,j), sumdata
            end if
         end do
         if (k .gt. 1) then
            epa(i) = sumdata/real(k)
            call croux (res_data, min(k,1000), ere(i)) ! Robust scale estimator Sn
         else
            epa(i) = 0.0
            ere(i) = 1.0
            write (msg,'(a,i1,2a)') 'mlocout_rderr: k = ', k, ' for ', qname1(i)
            call warnings (trim(msg))
         end if
         if (debug) write (io_log,'(a,3x,2f10.3)') qname1(i), epa(i), ere(i)
         t1 = lath_epa*rpd  ! Convert to geocentric radians
         t2 = lonh_epa*rpd  ! Convert to geocentric radians
         t3 = qlat(i)*rpd  ! Convert to geocentric radians
         t4 = qlon(i)*rpd  ! Convert to geocentric radians
         call delaz (t1, t2, t3, t4, dum1, delt1, dum2, dum3, azeqst, dum5, dum6, 1)
         call geogra (qlat(i), qlatdg)
         qlondg = qlon(i)
         call set_longitude_range (qlondg, 0)
         write (io_out,'(a,1x,i3,7f10.3)') qname1(i), k, epa(i), ere(i), rderr0(i), qlatdg, qlondg, delt1, azeqst
         if (abs(epa(i)) .gt. 10.0 .or. ere(i) .gt. 3.0) then
            write (msg,'(2a,2f10.3)') 'mlocout_rderr: possible outliers ', qname1(i), epa(i), ere(i)
            call fyi (trim(msg))
         end if
      end if
   end do ! End of loop over station-phases with data
   close (io_out)
   
   return
   
end subroutine write_rderr


!***********************************************************************************************************************************
subroutine write_rderr_diff ()

! Summary file of updated reading errors for differential time data. The estimate is made from
! all differential times for the same station and phase. Empirical path anomalies are not calculated.
! The measure of spread is 'sn', a robust, efficient scale estimator (Croux & Rousseeuw, 1992)
      
   integer :: i
   integer :: j
   integer :: k
   real, dimension(n_qres_max) :: res_data
   character(len=132) :: msg
   character(len=100) :: outfil
   
   outfil=trim(outfile)//'.rderr_diff'
   call open_file (io_out, outfil, 'new')
   
   do i = 1,nqc ! Loop over station-phases with data
      if (indexq(i) .ge. 2 .and. idiff0(i) .ne. 0) then
         k = 0
         do j = 1,indexq(i)
            ! On rare occasions, the calibration shift during indirect calibration causes the epicentral distance
            ! of a reading to fall off the phase branch and the residual becomes very large, upsetting the calculation of 
            ! the empirical reading error (ere). The test on
            ! absolute value of qres is meant to remove such readings from the calculation.
            if (abs(qres(i,j)) .lt. 20.) then
               k = k + 1
               res_data(k) = qres(i,j)
               if (debug) write (io_log,'(a,3x,i3,f10.3)') qname1(i), j, qres(i,j)
            end if
         end do
         if (k .gt. 1) then
            call croux (res_data, min(k,1000), ere(i)) ! Robust scale estimator Sn
         else
            ere(i) = 1.0
            write (msg,'(a,i1,2a)') 'mlocout_rderr_diff: k = ', k, ' for ', qname1(i)
            call warnings (trim(msg))
         end if
         if (debug) write (io_log,'(a,3x,f10.3)') qname1(i), ere(i)
         write (io_out,'(a,1x,i3,2f10.3)') qname1(i), k, ere(i), rderr0(i)
         if (ere(i) .gt. 0.5) then ! Tolerance for outliers is smaller for differential time data
            write (msg,'(2a,f10.3)') 'mlocout_rderr_diff: possible outliers ', qname1(i), ere(i)
            call fyi (trim(msg))
         end if
      end if
   end do ! End of loop over station-phases with data
   close (io_out)
   
   return
   
end subroutine write_rderr_diff


!*****************************************************************************************
   subroutine write_ttsprd (it)
   
   ! Runs through all possible phases in the tau-p software, collects the readings,
   ! calculates the spread (as robust scale estimator Sn), and writes an output file.
   ! Since the calculation of spread by subroutine croux uses only the first 1000
   ! readings of a given phase, that is the limit set in nphdat_max.
         
      integer, parameter :: nphases_max = 100 ! Maximum number of different phases that can be processed
      integer, parameter :: nphdat_max = 1000 ! Maximum number of readings of a given phase that will be used for the calculation
      
      integer :: i
      integer :: iev
      integer :: ird
      integer :: it
      integer :: j
      integer :: jout
!      integer :: lunit
      integer, dimension(nphases_max) :: nph
      integer :: nphases
      real :: average
      real, dimension(nphdat_max) :: d
      real :: dts
      real, dimension(nphases_max,nphdat_max) :: phdat
      real :: sn
      real :: total
      character(len=132) :: msg
      character(len=100) :: outfil
      character(len=8), dimension(nphases_max) :: ph
      logical :: newphase
      
      ph(1) = 'P       '
      nphases = 1
      nph(1) = 0
      
      do iev = 1,n_event
         do ird = 1,nst(iev)
            if (idiff(iev,ird) .gt. 0) cycle ! Skip differential time data for this purpose, since absolute times are dummy values
            if (.not.fltrh(iev,ird) .or. .not.fltrc(iev,ird)) then ! Only consider data used for the cluster vectors or hypocentroid
               dts = dt(iev,ird,it) - s(iev,ird,it)
               do j = 1,nphases
                  newphase = .true.
                  if (phase(iev,ird)(1:7) .eq. ph(j)(1:7)) then
                     newphase = .false.
                     if (nph(j) .lt. nphdat_max) then
                        nph(j) = nph(j) + 1
                        phdat(j,nph(j)) = dts
                     end if
                     exit
                  end if
               end do
               if (newphase) then
                  if (nphases .lt. nphases_max) then
                     nphases = nphases + 1
                     ph(nphases) = phase(iev,ird)
                     nph(nphases) = 1
                     phdat(nphases,nph(nphases)) = dts
                  else
                     write (msg,'(a,i3,2a)') 'mlocout_ttsprd: maximum number of phases (',&
                      nphases_max,') reached: ', phase(iev,ird)
                     call warnings (trim(msg))
                     if (verbose_log) write (io_log,'(a)') msg
                  end if
               end if
            end if
         end do
      end do
   
      ! Output file
      outfil = trim(outfile)//'.ttsprd'
      jout = lunit()
      call open_file (jout, outfil, 'new')
      do i = 1,nphases
         total = 0.
         if (nph(i) .ge. 5) then
            do j = 1,nph(i)
               d(j) = phdat(i,j)
               total = total + d(j)
            end do
            average = total/real(nph(i))
            call croux (d, min(nph(i),1000), sn)
            write (jout,'(a8,i8,2f10.3)') ph(i), nph(i), sn, average
         end if
      end do
      close (jout)
      
      return
      
   end subroutine write_ttsprd
         
         
!*****************************************************************************************
   subroutine write_tt (it)
   
   ! Outout of empirical TT data for specific phases (command "ttou")
         
      integer :: i
      integer :: iev
      integer :: ii
      integer, dimension(n_arrtim_max) :: indx
      integer :: it
      integer :: j
      integer :: jj
      real :: azes_pr
      real :: delt_pr
      real, dimension(n_arrtim_max) :: deltiev
      real :: dts
      character(len=160) :: command_line
      character(len=132) :: file_folder
      character(len=132) :: msg
      character(len=132) :: tt_file
      character(len=160) :: tt_file_path
      
      ! Make a folder for the files
      file_folder = trim(datadir)//dirsym//trim(basename)//'_tt'
      command_line = 'mkdir '//trim(file_folder)
      call execute_command_line (trim(command_line))
         
      ! Loop over requested phases
      do j = 1, n_ttou
         if (j .gt. 1) then
            do jj = 1, j-1
               if (ttou_phase(jj) .eq. ttou_phase(j)) then
                  msg = 'mlocout_tt: repeated phase (,'//trim(ttou_phase(j))//')'
                  call warnings (trim(msg))
                  cycle
               end if
            end do
         end if
         tt_file = trim(basename)//'_'//trim(ttou_phase(j))//'.ttou'
         tt_file_path = trim(file_folder)//dirsym//trim(tt_file)
         call open_file (io_ttou, tt_file_path, 'new')
         
         if (ttou_phase(j) .eq. 'S-P     ') then ! S-P is handled separately
            call mlocout_sp (it)
         else
            ! Loop over events
            do iev = 1, n_event
               ! Index table on delta for sorting
               do i = 1,nst(iev)
                  if (indirect_cal) then
                     deltiev(i) = delt_cal(iev,i)
                  else
                     deltiev(i) = delt(iev,i)
                  end if
               end do
               call indexx (nst(iev), deltiev, indx)
               ! Loop over readings, sorted by delta
               do ii = 1,nst(iev)
                  i = indx(ii)
                  if (indirect_cal) then
                     dts = dt_cal(iev,i) - s_cal(iev,i)
                     delt_pr = delt_cal(iev,i)
                     azes_pr = azes_cal(iev,i)
                  else
                     dts = dt(iev,i,it) - s(iev,i,it)
                     delt_pr = delt(iev,i)
                     azes_pr = azes(iev,i)
                  end if
                  if (phase(iev,i) .eq. ttou_phase(j))&
                   write (io_ttou,'(i3,1x,a16,1x,f5.1,1x,i5,1x,a20,1x,f8.3,1x,i3,1x,3f10.2,3x,a)')&
                   iev, evtnam(iev)(1:16), depthp(iev,it), mnf_line(iev,i), sad(iev,i), delt_pr,&
                   nint(azes_pr), tto(iev,i), dts, sdread(iev,i), fcode(iev,i)
               end do
            end do
         end if
         
         close (io_ttou)
      end do
      
      return
   
   end subroutine write_tt
   
   
!*****************************************************************************************
   subroutine mlocout_sp (it)
   
   ! Output file with S-P data
   
      integer :: iev
      integer :: ird
      integer :: it
      integer :: lastiev
      integer :: nevsp
      integer :: nsp
      integer :: nspx
      character(len=132) :: msg
         
      nevsp = 0
      nsp = 0
      nspx = 0
      lastiev = 0
         
      do iev = 1,n_event ! Loop over events
         do ird = 1,nst(iev) ! Loop over readings
            if (phase(iev,ird) .eq. 'S-P     ') then
               call sp_write (iev, ird, it)
               if (iev .gt. lastiev) then
                  nevsp = nevsp + 1
                  lastiev = iev
               end if
               nsp = nsp + 1
               if (fcode (iev,ird) .ne. ' ') nspx = nspx + 1
            end if
         end do
      end do
         
      write (msg,'(a,i3,a)') 'mlocout_sp: ', nevsp, ' events with S-P readings'
      call fyi (trim(msg))
      write (io_log,'(/a)') trim(msg)
      write (msg,'(a,i3,a,i3,a)') 'mlocout_sp: ', nsp, ' S-P readings (', nspx, ' flagged)'
      call fyi (trim(msg))
      write (io_log,'(a)') trim(msg)
      
      return
      
   end subroutine mlocout_sp
   
   
!*****************************************************************************************
   subroutine sp_write (iev, ird, it)
   
   ! Write an output file for S-P data
   
      ! Arguments
      integer, intent(in) :: iev
      integer, intent(in) :: ird
      integer, intent(in) :: it
      
      integer :: it1
      real :: azesout
      real :: azseout
      real :: ddep
      real :: deltout
      real :: dot
      real :: dts
      real :: fdm
      real :: fdp
      real :: scale_length
      real :: stladg_geog
      real :: stlndg180
!       real :: xl1
!       real :: xl2
      character(len=164) :: fmt
      
      it1 = it + 1
      
      ! Uncalibrated cluster, uses cluster vector confidence ellipse
!       xl1 = xl1c(iev)
!       xl2 = xl2c(iev)
      ddep = sdxhatc(iev,3)
      dot = sdxhatc(iev,4)
      hourpr = hourp(iev,it1)
      minpr = minp(iev,it1)
      secpr = secp(iev,it1)
      latout = latp(iev,it1)
      lonout = lonp(iev,it1)
      depout = depthp(iev,it1)
      if (direct_cal) then ! Direct calibration
!          xl1 = xl1dc(iev)
!          xl2 = xl2dc(iev)
         ddep = ddepdc(iev)
         dot = dotdc(iev)
         hourpr = hourp(iev,it1)
         minpr = minp(iev,it1)
         secpr = secp(iev,it1)
         latout = latp(iev,it1)
         lonout = lonp(iev,it1)
         depout = depthp(iev,it1)
      end if
      if (indirect_cal) then ! Indirect calibration trumps direct calibration
         ddep = sqrt(accv(iev,3,3))
         dot = sqrt(accv(iev,4,4))
         call timecr (otsp_cal(iev), hourpr, minpr, secpr)
         latout = latp_cal(iev)
         lonout = lonp_cal(iev)
         depout = depthp_cal(iev)
      end if
      call set_longitude_range (lonout, longitude_range)
         
      ! Epicentral distance and azimuth from event to station
      deltout = delt(iev,ird)
      azesout = azes(iev,ird)
      azseout = azse(iev,ird)
      if (indirect_cal) then
         deltout = delt_cal(iev,ird)
         azesout = azes_cal(iev,ird)
         azseout = azse_cal(iev,ird)
      end if   
   
      ! Scale length ! Semi-major axis of the confidence ellipse
      if (direct_cal) then ! Direct calibration
         scale_length = xl2dc(iev)
      else if (indirect_cal) then ! Indirect calibration
         scale_length = xl2cg(iev)
      else
         scale_length = xl2c(iev) ! Uncalibrated
      end if
      scale_length = scale_length / 111.2 ! Convert km to degrees
   
      ! Travel time residual
      dts = dt(iev,ird,it) - s(iev,ird,it) ! Direct calibration or uncalibrated
      if (indirect_cal) dts = dt_cal(iev,ird) - s_cal(iev,ird) ! Indirect calibration trumps others
   
      ! Convert from geocentric to geographic latitude
      call geogra (stladg(iev,ird), stladg_geog) 
      
      ! Station longitude runs from -180 to 180
      stlndg180 = stlndg(iev,ird)
      call set_longitude_range (stlndg180, 0)
   
      ! Uncertainty of focal depth
      if (mindx(iev,3) .gt. 0) then ! Free depth solution
         fdp = ddep
         fdm = ddep
      else ! Take from default or assigned uncertainties
         fdp = depthp_plus(iev)
         fdm = depthp_minus(iev)
      end if
      
      fmt = '(i3,1x,a5,1x,f6.2,1x,f5.2,1x,f6.2,f6.3,2f6.1,1x,f5.3,1x,i4,2i2.2,1x,2(i2,1x),f4.1,1x,&
       &f4.2,f9.4,f10.4,f6.1,2f5.1,f9.4,f10.4,f6.3,1x,a1)'
      write (io_ttou,fmt)&
       iev,& ! Event number 
       stname(iev,ird),& ! Station code
       tto(iev,ird),& ! S-P time, s
       sdread(iev,ird),& ! Reading error (s)
       dts,& ! TT residual (s)
       deltout,& ! Epicentral distance (deg)
       azseout,& ! Back-azimuth, station to event
       azesout,& ! Azimuth, event to station
       scale_length,& ! Uncertainty in epicentral distance, degrees
       iyre(iev),& ! Origin year.
       mone(iev),& ! Origin month.
       idye(iev),& ! Origin day.
       hourpr,& ! Origin hour.
       minpr,& ! Origin minute.
       secpr,& ! Origin seconds.
       dot,& ! Standard error in origin time (s).
       latout,& ! Event latitude
       lonout,& ! Event longitude
       depout,& ! Depth
       fdp,& ! Error in depth, positive (km).
       fdm,& ! Error in depth, minus (km).
       stladg_geog,& ! Station latitude
       stlndg180,& ! Station longitude
       ahgts(iev,ird),& ! Station elevation, km
       fcode(iev,ird) ! flag
   
      return
      
   end subroutine sp_write
   
   
!*****************************************************************************************
   subroutine write_gccel (it)
   
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
!      integer :: lunit
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
   
      write (io_gccel,'(a,i3)') '# Number of events: ', n_event
      write (commentary_buffer(1),'(a,i3)') 'Number of events: ', n_event
      
      ! Calibration type and calibration statistics
      eclmin = 99
      eclmax = 0
      if (indirect_cal) then ! Indirect calibration
         write (io_gccel,'(a,i3,a,f4.1,a)') '# Calibration type: indirect calibration on ', ncal(3),&
          ' calibration events; hypocentroid calibration level = ', cal_level, ' km'
         write (commentary_buffer(2),'(a,i3,a)') 'Calibration type: indirect calibration on ', ncal(3),&
          ' calibration events'
         write (commentary_buffer(3),'(a,f4.1,a)') 'Hypocentroid calibration level = ', cal_level, ' km'
         do iev = 1,n_event
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
         do iev = 1,n_event
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
      write (datelast,'(i4,2i2.2)') iyre(n_event), mone(n_event), idye(n_event)
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
      do iev = 1,n_event
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
      do iev = 1,n_event
         call write_mnf_14 (io_gccel, iev, it)
      end do
      write (io_gccel,'(a,t121,a)') 'EOF', ' '
      close (io_gccel)
   
      return
      
   end subroutine write_gccel
   
   
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
      
      do iev = 1,n_event ! Loop over events
   
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
!      integer :: lenb
!      integer :: lunit
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
   subroutine gccel_epicenters ()
   
   ! Writes an entry for the "epicenters.kml" file used in GCCEL distributions to display
   ! epicenters of calibrated clusters.

      integer :: iev
!       integer :: it1
      character(len=116) :: ev_name
      character(len=100) :: outfil
      
!       it1 = it + 1
      
      call fyi ('gccel_epicenters: writing epicenters.kml entry')
      outfil = trim(gcat_folder)//dirsym//trim(basename)//'_epicenters.kml'
      call open_file (io_gccel, outfil, 'new')
      
      do iev = 1,n_event
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
      character(len=*), intent(in) :: hypocenter
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
   subroutine tomo (ntomo, it, tomo_phase_in, itomo_in, folder)
         
   ! Output file for use in tomography. "tomo_phase" is a phase name for which data are extracted and "itomo" is
   ! a flag for the type of readings included:
   !  itomo = 1 Extract all readings of the specified phase
   !  itomo = 2 Extract only readings which were used for the cluster vectors
   !  itomo = 3 Extract empirical path anomalies
   
      integer :: i
      integer, dimension(n_arrtim_max) :: idum
      integer :: iev
      integer :: ii
      integer, dimension(n_arrtim_max) :: indx
      integer, dimension(100) :: indx2
      integer :: it
      integer :: it1
      integer :: itomo_in
      integer :: j
      integer :: jj
      integer :: kk
      integer :: ntomo
      real, dimension(n_arrtim_max) :: deltiev
      real, dimension(100) :: sortime
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
         do iev = 1,n_event ! Loop over events
         
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
   
      ! Travel time variables
      integer :: nphase
      real, dimension(max_mpl) :: dtdd
      real, dimension(max_mpl) :: dtdh
      real, dimension(max_mpl) :: dddp
      real, dimension(max_mpl) :: tt ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd

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
               call trtm (delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
            else
               if (verbose_log) write (io_log,'(a)') 'Calling tt_mixed_model from tomo3'
               call tt_mixed_model (depth_epa, delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
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
   
      ! Travel time variables
      integer :: nphase
      real, dimension(max_mpl) :: dtdd
      real, dimension(max_mpl) :: dtdh
      real, dimension(max_mpl) :: dddp
      real, dimension(max_mpl) :: tt ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd

      integer, dimension(n_event) :: cell_events ! list of event numbers for a cell
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
      integer, dimension(n_event,2) :: readings ! event and reading numbers for observations of the station-phase in a cell
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
      real, dimension(n_event):: residual
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
            do iev = 1, n_event
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
                        if (n_obs .lt. n_event) then
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
                        else if (n_obs .eq. n_event) then
                           call warnings ('tomo3grid: maximum number of observations (n_event) reached for '//qname1(isp))
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
                  call trtm (delt1, nphase, tt, dtdd, dtdh, dddp, phcd)
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
      do iev = 1,n_event
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
   subroutine write_bloc (it)
         
   ! Output file in BAYESLOC format.
   ! Used to import mloc calibrated locations as priors in BayesLoc
   
      ! Argument
      integer, intent(in) :: it

      integer, parameter :: real_kind = selected_real_kind(16,30)
      real(kind=real_kind) :: xot, xot_idc

      integer :: calibration_type
      integer :: iev
      integer :: it1
!       real :: alpha
!       real :: alpha_idc
      real :: ddep
      real :: ddep_idc
      real :: depout_idc
      real :: dot
      real :: dot_idc
      real :: latout_idc
      real :: lonout_idc
      real :: secpr_idc
!       real :: xl1
!       real :: xl1_idc
!       real :: xl2
!       real :: xl2_idc
      character(len=4) :: cal_code
      character(len=1) :: cal2
      character(len=4) :: fdu_pr
      character(len=132) :: fmt
      character(len=102) :: header
      character(len=132) :: msg
      character(len=100) :: outfil
   
      outfil = trim(outfile)//'.bloc'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_bloc: opening ', trim(outfil), ' on unit ', io_bayes
         call fyi (trim(msg))
      end if
      call open_file (io_bayes, outfil, 'new')
   
      it1 = it + 1
      fmt = '(i4,2i2.2,a1,2i2.2,a1,i2.2,f10.5,f11.5,3x,f10.3,f7.1,3x,a4,x,f20.3,f6.2)' 
      header = 'ev_id lat_mean lon_mean dist_sd depth_mean depth_sd time_mean time_sd'
      
      if (indirect_cal) then ! Indirect calibration
         write (io_bayes,'(a)') 'Indirect calibration'
         write (io_bayes,'(a)') header
         do iev = 1,n_event
         
            ! Hypocentral uncertainties
!             xl1_idc = xl1cg(iev)
!             xl2_idc = xl2cg(iev)
!             alpha_idc = alphacg(iev)
            ddep_idc = sqrt(accv(iev,3,3))
            dot_idc = sqrt(accv(iev,4,4))
            call timecr (otsp_cal(iev), hourpr, minpr, secpr_idc)
            latout_idc = latp_cal(iev)
            lonout_idc = lonp_cal(iev)
            depout_idc = depthp_cal(iev)
            calibration_type = 2
   
   
            ! Calibration code
            call calibration_code (iev, calibration_type, cal_code)
            cal2 = cal_code(2:2)
            call uctolc (cal2,-1)  
            
            call focal_depth_uncertainty (iev, mindx(iev,3), ddep_idc, fdu_pr)    
               
            ! Unix time
            call unix_time (iyre(iev), mone(iev), idye(iev), hourp(iev,it1), minp(iev,it1), secpr_idc, xot_idc)
             
            write (io_bayes,fmt)&
             iyre(iev),& ! Origin year.
             mone(iev),& ! Origin month.
             idye(iev),& ! Origin day.
             '.',&
             hourpr,& ! Origin hour.
             minpr,& ! Origin minute.
             '.',&
             int(secpr_idc),& ! Origin seconds.
             latout_idc, & ! Geographic latitude.
             lonout_idc,& ! Geographic longitude.
             eqr(iev),& !Uncertanity in lat & long (km)
             depout_idc,& ! Final depth (km).
             fdu_pr,& ! Uncertainty in depth (km)
             xot_idc,& !Origin time in epoch (sec)
             dot_idc  ! Uncertainty in origin time (sec).
         end do
         close (io_bayes)
         return
      end if
      
      ! Output for direct calibration or uncalibrated
      if (direct_cal) then
         write (io_bayes,'(a)') 'Direct calibration'
      else
         write (io_bayes,'(a)') 'Uncalibrated'
      end if
      write (io_bayes,'(a)') header
   
      do iev = 1,n_event
         if (direct_cal) then ! Direct calibration
!             xl1 = xl1dc(iev)
!             xl2 = xl2dc(iev)
!             alpha = alphadc(iev)
            ddep = ddepdc(iev)
            dot = dotdc(iev)
            hourpr = hourp(iev,it1)
            minpr = minp(iev,it1)
            secpr = secp(iev,it1)
            latout = latp(iev,it1)
            lonout = lonp(iev,it1)
            depout = depthp(iev,it1)
            calibration_type = 1
         else ! Uncalibrated
!             xl1 = xl1c(iev)
!             xl2 = xl2c(iev)
!             alpha = alphac(iev)
            ddep = sdxhatc(iev,3)
            dot = sdxhatc(iev,4)
            hourpr = hourp(iev,it1)
            minpr = minp(iev,it1)
            secpr = secp(iev,it1)
            latout = latp(iev,it1)
            lonout = lonp(iev,it1)
            depout = depthp(iev,it1)
            calibration_type = 0
         end if
         
   
         ! Calibration code
         call calibration_code (iev, calibration_type, cal_code)
         cal2 = cal_code(2:2)
         call uctolc (cal2,-1)
            
         call focal_depth_uncertainty (iev, mindx(iev,3), ddep, fdu_pr)
               
         ! Unix time
         call unix_time (iyre(iev), mone(iev), idye(iev), hourp(iev,it1), minp(iev,it1), secp(iev,it1), xot)
          
         write (io_bayes,fmt)&
          iyre(iev),& ! Origin year.
          mone(iev),& ! Origin month.
          idye(iev),& ! Origin day.
          '.',&
          hourpr,& ! Origin hour.
          minpr,& ! Origin minute.
          '.',&
          int(secpr),& ! Origin seconds.
          latout,& ! Geographic latitude.
          lonout,& ! Geographic longitude.
          eqr(iev),& !Uncertanity in lat & long (km)
          depout,& ! Final depth (km).
          fdu_pr,& ! Uncertainty in depth (km)
          xot,& !Origin time in epoch (sec)
          dot ! Uncertainty in origin time (sec).
   
      end do
   
      close (io_bayes)      
      return
   
   end subroutine write_bloc
   
   
!*****************************************************************************************
   subroutine focal_depth_uncertainty (iev, mindx3, ddep, fdu_pr)
   
      integer, intent(in) :: iev
      integer, intent(in) :: mindx3
      real, intent(in) :: ddep
      real, parameter :: default_depth_uncertainty = 5. ! For events set to cluster default depth
      real :: fdu
      real, parameter :: min_depth_uncertainty = 0.1 ! To avoid rounding to zero on output
      character(len=4), intent(out) :: fdu_pr
   
      if (mindx3 .gt. 0) then ! Free depth solution
            fdu = ddep
      else ! Take from default or assigned uncertainties
         if (depthp_plus(iev) .le. 99.) then
            fdu = max(depthp_plus(iev),depthp_minus(iev))
         else
            fdu = default_depth_uncertainty ! default value for events set to cluster default depth
         end if
      end if
      write (fdu_pr,'(f4.1)') max(fdu,min_depth_uncertainty)
      
      return
      
   end subroutine focal_depth_uncertainty
   
   


!*****************************************************************************************
end module mloc_out

