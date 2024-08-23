!> Procedures to allocate sets of variable arrays related to specific parts of mloc.

module mloc_allocate

   use mloc_declare

   implicit none
   save
   
   integer :: error ! allocation status error code
   character(len=40) :: p ! procedure name
   
contains

!*****************************************************************************************
   subroutine events_allocate ()
   
   ! Allocate variable arrays related to events, some are initialized too
   
      integer :: i
      integer :: it
      integer :: jt
      integer :: n
!      integer :: m
      
      it = n_it_max1
      jt = n_it_max2
      n = n_event
!      m = n_hdf_max
      p = 'events_allocate' ! procedure
      
      allocate (alic(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'alic', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': alic allocated (', n, ')'
      end if
      
      allocate (alphac(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'alphac', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': alphac allocated (', n, ')'
      end if
      
      allocate (annotation(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'annotation', error)
      else
         annotation = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': annotation allocated (', n, ')'
      end if
      
      allocate (blic(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'blic', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': blic allocated (', n, ')'
      end if
      
      allocate (cal_code_input(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_code_input', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_code_input allocated (', n, ')'
      end if
      
      allocate (cal_event(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_event', error)
      else
         cal_event = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_event allocated (', n, ',3)'
      end if
      
      allocate (cal_par(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_par', error)
      else
         cal_par = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_par allocated (', n, ')'
      end if
      
      allocate (ccv(n,4,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ccv', error)
      else
         ccv = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': ccv allocated (', n, ',4,4)'
      end if
      
      allocate (depset_pr(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depset_pr', error)
      else
         depset_pr = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depset_pr allocated (', n, ')'
      end if
      
      allocate (depth_cf(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depth_cf', error)
      else
         do i = 1,n_event
            depth_cf(i,1) = -999.
            depth_cf(i,2) = 99.9
            depth_cf(i,3) = 99.9
         end do
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depth_cf allocated (', n, ',3)'
      end if
      
      allocate (depth_cf_c(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depth_cf_c', error)
      else
         depth_cf_c = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depth_cf_c allocated (', n, ')'
      end if
      
      allocate (depth_hdf_c(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depth_hdf_c', error)
      else
         depth_hdf_c = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depth_hdf_c allocated (', n, ')'
      end if
      
      allocate (depth_hdf(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depth_hdf', error)
      else
         do i = 1,n_event
            depth_hdf(i,1) = -999.
            depth_hdf(i,2) = 99.9
            depth_hdf(i,3) = 99.9
         end do
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depth_hdf allocated (', n, ',3)'
      end if
      
      allocate (depth_inp_c(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depth_inp_c', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depth_inp_c allocated (', n, ')'
      end if
      
      allocate (depth_inp(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depth_inp', error)
      else
         do i = 1,n_event
            depth_inp(i,1) = -999.
            depth_inp(i,2) = 99.9
            depth_inp(i,3) = 99.9
         end do
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depth_inp allocated (', n, ',3)'
      end if
      
      allocate (depthf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depthf', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depthf allocated (', n, ')'
      end if
      
      allocate (depthp(n,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depthp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i1,a)') trim(p)//': depthp allocated (0:', it, ')'
      end if
      
      allocate (depthp_minus(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depthp_minus', error)
      else
         depthp_minus = 99.9
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depthp_minus allocated (', n, ')'
      end if
      
      allocate (depthp_plus(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depthp_plus', error)
      else
         depthp_plus = 99.9
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depthp_plus allocated (', n, ')'
      end if
      
      allocate (depthpr(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depthpr', error)
      else
         depthpr = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depthpr allocated (', n, ')'
      end if
      
      allocate (event_comment(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'event_comment', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': event_comment allocated (', n, ')'
      end if
      
      allocate (evid_source(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'evid_source', error)
      else
         evid_source = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': evid_source allocated (', n, ')'
      end if
      
      allocate (evid(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'evid', error)
      else
         evid = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': evid allocated (', n, ')'
      end if
      
      allocate (evtnam(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'evtnam', error)
      else
         evtnam = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': evtnam allocated (', n, ')'
      end if
      
      allocate (fixpr(n,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fixpr', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': fixpr allocated (', n, ',4)'
      end if
      
      allocate (hour_inp(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hour_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hour_inp allocated (', n, ')'
      end if
      
      allocate (hourp(n,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hour_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hour_inp allocated (', n, ')'
      end if
      
      allocate (hypo_author(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hypo_author', error)
      else
         hypo_author = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hypo_author allocated (', n, ')'
      end if
      
      allocate (hypocenter_list(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hypocenter_list', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hypocenter_list allocated (', n, ')'
      end if
      
      allocate (idye(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'idye', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': idye allocated (', n, ')'
      end if
      
      allocate (infile(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'infile', error)
      else
         infile = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': infile allocated (', n, ')'
      end if
      
      allocate (infile20(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'infile20', error)
      else
         infile = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': infile20 allocated (', n, ')'
      end if
      
      allocate (iyre(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'iyre', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': iyre allocated (', n, ')'
      end if
      
      allocate (kcritc(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'kcritc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': kcritc allocated (', n, ')'
      end if
      
      allocate (lat_cf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lat_cf', error)
      else
         lat_cf = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lat_cf allocated (', n, ')'
      end if
      
      allocate (lat_hdf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lat_hdf', error)
      else
         lat_hdf = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lat_hdf allocated (', n, ')'
      end if
      
      allocate (lat_inp(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lat_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lat_inp allocated (', n, ')'
      end if
      
      allocate (latf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'latf', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': latf allocated (', n, ')'
      end if
      
      allocate (latp(n,0:jt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'latp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': latp allocated (', n, ',0:', jt, ')'
      end if
      
      allocate (latpgc(n,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'latpgc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': latpgc allocated (', n, ',0:', it, ')'
      end if
      
      allocate (latpr(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'latpr', error)
      else
         latpr = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': latpr allocated (', n, ')'
      end if
      
      allocate (lon_hdf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lon_hdf', error)
      else
         lon_hdf = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lon_hdf allocated (', n, ')'
      end if
      
      allocate (lon_inp(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lon_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lon_inp allocated (', n, ')'
      end if
      
      allocate (lon_cf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lon_cf', error)
      else
         lon_cf = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lon_cf allocated (', n, ')'
      end if
      
      allocate (lonf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lonf', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lonf allocated (', n, ')'
      end if
      
      allocate (lonp(n,0:jt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lonp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': lonp allocated (', n, ',0:', jt, ')'
      end if
      
      allocate (lonpr(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lonpr', error)
      else
         lonpr = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lonpr allocated (', n, ')'
      end if
      
      allocate (magnitude_author(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'magnitude_author', error)
      else
         magnitude_author = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': magnitude_author allocated (', n, ')'
      end if
      
      allocate (min_inp(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'min_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': min_inp allocated (', n, ')'
      end if
      
      allocate (mindx(n,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'mindx', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': mindx allocated (', n, ',4)'
      end if
      
      allocate (minp(n,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'minp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': minp allocated (', n, ',0:', it, ')'
      end if
      
      allocate (mmag(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'mmag', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': mmag allocated (', n, ')'
      end if
      
      allocate (mone(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'mone', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': mone allocated (', n, ')'
      end if
      
      allocate (mtiev(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'mtiev', error)
      else
         mtiev = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': mtiev allocated (', n, ')'
      end if
      
      allocate (nevt(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'nevt', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': nevt allocated (', n, ')'
      end if
      
      allocate (otsp(n,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'otsp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': otsp allocated (', n, ',0:', it, ')'
      end if
      
      allocate (rmag(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rmag', error)
      else
         rmag = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': rmag allocated (', n, ')'
      end if
      
      allocate (secp(n,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'secp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': secp allocated (', n, ',0:', it, ')'
      end if
      
      allocate (sec_inp(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sec_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': sec_inp allocated (', n, ')'
      end if
      
      allocate (time_cf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'time_cf', error)
      else
         time_cf = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': time_cf allocated (', n, ')'
      end if
      
      allocate (time_hdf(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'time_hdf', error)
      else
         time_hdf = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': time_hdf allocated (', n, ')'
      end if
      
      allocate (time_inp(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'time_inp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': time_inp allocated (', n, ')'
      end if
      
      allocate (timef(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'timef', error)
      else
         timef = .true.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': timef allocated (', n, ')'
      end if
      
      allocate (timepr(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'timepr', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': timepr allocated (', n, ')'
      end if
            
      allocate (xl1c(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'xl1c', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': xl1c allocated (', n, ')'
      end if
      
      allocate (xl2c(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'xl2c', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': xl2c allocated (', n, ')'
      end if
      
      return
   
   end subroutine events_allocate


!*****************************************************************************************
   subroutine hdf_allocate ()
   
   ! Allocation and initialization of arrays related to HDF files
   
      integer :: m
      
      p = 'hdf_allocate'
      
      m = n_hdf_max
   
      allocate (hdf_day(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_day', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_day allocated (', m, ')'
      end if
      
      allocate (hdf_dep(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_dep', error)
      else
         hdf_dep = -999.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_dep allocated (', m, ')'
      end if
      
      allocate (hdf_mon(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_mon', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_mon allocated (', m, ')'
      end if
      
      allocate (hdf_yr4(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_yr4', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_yr4 allocated (', m, ')'
      end if
      
      allocate (hdf_dep_code(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_dep_code', error)
      else
         hdf_dep_code = 'u '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_dep_code allocated (', m, ')'
      end if
      
      allocate (hdf_dep_minus(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_dep_minus', error)
      else
         hdf_dep_minus = 99.9
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_dep_minus allocated (', m, ')'
      end if
      
      allocate (hdf_dep_plus(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_dep_plus', error)
      else
         hdf_dep_plus = 99.9
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_dep_plus allocated (', m, ')'
      end if
      
      allocate (hdf_evid(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_evid', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_evid allocated (', m, ')'
      end if
      
      allocate (hdf_lat(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_lat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_lat allocated (', m, ')'
      end if
      
      allocate (hdf_lon(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_lon', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_lon allocated (', m, ')'
      end if
      
      allocate (hdf_time(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdf_time', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdf_time allocated (', m, ')'
      end if
      
      allocate (hdfline(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hdfline', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': hdfline allocated (', m, ')'
      end if
      
      return
   
   end subroutine hdf_allocate


!*****************************************************************************************
   subroutine direct_cal_allocate ()
   
   ! Allocation and initialization of arrays related to direct calibration that are used in several
   ! modules. See procedures in module 'mloc_calibration' for allocation and deallocation
   ! of arrays used only within that module.
      
!      integer :: m
      integer :: n
!      integer :: n2
      
      n = n_event
!      n2 = n_event*2
!      m = n_arrtim_max
      p = 'direct_cal_allocate' ! procedure name

      allocate (alphadc(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'alphadc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': alphadc allocated (', n, ')'
      end if
      
      allocate (ddepdc(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ddepdc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': ddepdc allocated (', n, ')'
      end if
      
      allocate (dotdc(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dotdc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': dotdc allocated (', n, ')'
      end if
      
      ! eqr can be allocated for direct or indirect calibration
      if (.not.allocated(eqr)) then
         allocate (eqr(n), stat=error)
         if (error .gt. 0) then
            call allocation_error (p, 'eqr', error)
         else
            if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': eqr allocated (', n, ')'
         end if
      end if
      
      allocate (xl1dc(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'xl1dc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': xl1dc allocated (', n, ')'
      end if
      
      allocate (xl2dc(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'xl2dc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': xl2dc allocated (', n, ')'
      end if
      
      return
      
   end subroutine direct_cal_allocate


!*****************************************************************************************
   subroutine indirect_cal_allocate ()
   
   ! Allocation and initialization of arrays related to indirect calibration that are used in several
   ! modules. See procedures in module 'mloc_calibration' for allocation and deallocation
   ! of arrays used only within that module.
      
      integer :: m
      integer :: n
!      integer :: n2
      
      n = n_event
!      n2 = n_event*2
      m = n_arrtim_max
      p = 'indirect_cal_allocate' ! procedure name
   
      allocate (accv(n,4,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'accv', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': accv allocated (', n, ',4,4)'
      end if
      
      allocate (alphacg(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'alphacg', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': alphacg allocated (', n, ')'
      end if
      
      allocate (azes_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'azes_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': azes_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (azse_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'azse_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': azse_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (cal_lat(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_lat', error)
      else
         cal_lat = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_lat allocated (', n, ',3)'
      end if
      
      allocate (cal_lon(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_lon', error)
      else
         cal_lon = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_lon allocated (', n, ',3)'
      end if
      
      allocate (cal_min(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_min', error)
      else
         cal_min = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_min allocated (', n, ',3)'
      end if
      
      allocate (cal_sec(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_sec', error)
      else
         cal_sec = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_sec allocated (', n, ',3)'
      end if
      
      allocate (cal_dep(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_dep', error)
      else
         cal_dep = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_dep allocated (', n, ',3)'
      end if
      
      allocate (cal_hr(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'cal_hr', error)
      else
         cal_hr = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': cal_hr allocated (', n, ',3)'
      end if
      
      allocate (delt_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'delt_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': delt_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (depthp_cal(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'depthp_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': depthp_cal allocated (', n, ')'
      end if
      
      allocate (dt_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dt_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': dt_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (elcr_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'elcr_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': elcr_cal allocated (', n, ',', m, ')'
      end if
      
      ! eqr can be allocated for direct or indirect calibration
      if (.not.allocated(eqr)) then
         allocate (eqr(n), stat=error)
         if (error .gt. 0) then
            call allocation_error (p, 'eqr', error)
         else
            if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': eqr allocated (', n, ')'
         end if
      end if
      
      allocate (latp_cal(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'latp_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': latp_cal allocated (', n, ')'
      end if
      
      allocate (lonp_cal(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'lonp_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': lonp_cal allocated (', n, ')'
      end if
      
      allocate (otsp_cal(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'otsp_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': otsp_cal allocated (', n, ')'
      end if
      
      allocate (rcv(n,3,4,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rcv', error)
      else
         rcv = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': rcv allocated (', n, ',3,4,4)'
      end if
      
      allocate (s_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 's_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': s_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (sdstcs_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sdstcs_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sdstcs_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (stacs_cal(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stacs_cal', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': stacs_cal allocated (', n, ',', m, ')'
      end if
      
      allocate (w12(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'w12', error)
      else
         w12 = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': w12 allocated (', n, ')'
      end if
      
      allocate (w3(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'w3', error)
      else
         w3 = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': w3 allocated (', n, ')'
      end if
      
      allocate (w4(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'w4', error)
      else
         w4 = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': w4 allocated (', n, ')'
      end if
      
      allocate (xl1cg(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'xl1cg', error)
      else
         w4 = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': xl1cg allocated (', n, ')'
      end if
      
      allocate (xl2cg(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'xl2cg', error)
      else
         w4 = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': xl2cg allocated (', n, ')'
      end if
      
      return
      
   end subroutine indirect_cal_allocate
   
   
!*****************************************************************************************
   subroutine phase_data_allocate ()
   
   ! Allocate variable arrays related to phase data
   
      integer :: j
      integer :: k
      integer :: l
      integer :: m
      integer :: n
      
      j = n_qres_max
      k = n_q_max
      l = n_hrdp_max
      m = n_arrtim_max ! maximum number of phase records for any event in the cluster
      n = n_event  ! number of events in the cluster
      p = 'phase_data_allocate'
      
      allocate (adslcaa(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'adslcaa', error)
      else
         adslcaa = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': adslcaa allocated (', n, ',', m, ')'
      end if
      
      allocate (agency(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'agency', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': agency allocated (', n, ',', m, ')'
      end if
      
      allocate (channel(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'channel', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': channel allocated (', n, ',', m, ')'
      end if
      
      allocate (connected(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'connected', error)
      else
         connected = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': connected allocated (', n, ',', m, ')'
      end if
      
      allocate (deployment(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'deployment', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': deployment allocated (', n, ',', m, ')'
      end if
      
      allocate (diff_line(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'diff_line', error)
      else
         diff_line = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': diff_line allocated (', n, ',', m, ')'
      end if
      
      allocate (dist_az(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dist_az', error)
      else
         dist_az = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': dist_az allocated (', n, ',', m, ')'
      end if
      
      allocate (epa(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'epa', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': epa allocated (', k, ')'
      end if
      
      allocate (ere(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ere', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': ere allocated (', k, ')'
      end if
      
      allocate (fcode(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fcode', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fcode allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrc(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrc', error)
      else
         fltrc = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrc allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrcdelt(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrcdelt', error)
      else
         fltrcdelt = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrcdelt allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrcflag(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrcflag', error)
      else
         fltrcflag = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrcflag allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrcres(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrcres', error)
      else
         fltrcres = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrcres allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrh(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrh', error)
      else
         fltrh = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrh allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrhdelt(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrhdelt', error)
      else
         fltrhdelt = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrhdelt allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrhflag(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrhflag', error)
      else
         fltrhflag = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrhflag allocated (', n, ',', m, ')'
      end if
      
      allocate (fltrhres(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'fltrhres', error)
      else
         fltrhres = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': fltrhres allocated (', n, ',', m, ')'
      end if
      
      allocate (h_rdp(n,3), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'h_rdp', error)
      else
         fltrhres = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': h_rdp allocated (', n, ',3)'
      end if
      
      allocate (hrmsi(n,l), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'hrmsi', error)
      else
         fltrhres = .false.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': hrmsi allocated (', n, ',', l, ')'
      end if
      
      allocate (idiff(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'idiff', error)
      else
         idiff = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': idiff allocated (', n, ',', m, ')'
      end if
      
      allocate (idiff0(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'idiff0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': idiff0 allocated (', k, ')'
      end if
      
      allocate (indexq(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'indexq', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': indexq allocated (', k, ')'
      end if
      
      allocate (ipad(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ipad', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ipad allocated (', n, ',', m, ')'
      end if
      
      allocate (ipah(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ipah', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ipah allocated (', n, ',', m, ')'
      end if
      
      allocate (ipam(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ipam', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ipam allocated (', n, ',', m, ')'
      end if
      
      allocate (ipamo(n,0:m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ipamo', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ipamo allocated (', n, ',0:', m, ')'
      end if
      
      allocate (ipay(n,0:m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ipay', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ipay allocated (', n, ',0:', m, ')'
      end if
      
      allocate (iptim(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'iptim', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': iptim allocated (', n, ',', m, ')'
      end if
      
      allocate (iqiev(k,j), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'iqiev', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i3,a)') trim(p)//': iqiev allocated (', k, ',', j, ')'
      end if
      
      allocate (location(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'location', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': location allocated (', n, ',', m, ')'
      end if
      
      allocate (mnf_line(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'mnf_line', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': mnf_line allocated (', n, ',', m, ')'
      end if
      
      allocate (n_samples(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'n_samples', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': n_samples allocated (', n, ')'
      end if
      
      allocate (no_phid(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'no_phid', error)
      else
         no_phid = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': no_phid allocated (', n, ',', m, ')'
      end if
      
      allocate (nst(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'nst', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': nst allocated (', n, ')'
      end if
      
      allocate (pas(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'pas', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': pas allocated (', n, ',', m, ')'
      end if
      
      allocate (phase(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'phase', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': phase allocated (', n, ',', m, ')'
      end if
      
      allocate (phase0(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'phase0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': phase0 allocated (', n, ',', m, ')'
      end if
      
      allocate (phidird(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'phidird', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': phidird allocated (', n, ',', m, ')'
      end if
      
      allocate (qelev(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qelev', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qelev allocated (', k, ')'
      end if
      
      allocate (qlat(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qlat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qlat allocated (', k, ')'
      end if
      
      allocate (qlon(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qlon', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qlon allocated (', k, ')'
      end if
      
      allocate (qname(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qname', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qname allocated (', k, ')'
      end if
      
      allocate (qname1(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qname1', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qname1 allocated (', k, ')'
      end if
      
      allocate (qres(k,j), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qres', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i3,a)') trim(p)//': qres allocated (', k, ',', j, ')'
      end if
      
      allocate (rderr0(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rderr0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': rderr0 allocated (', k, ')'
      end if
      
      allocate (rdsigma(k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rdsigma', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': rdsigma allocated (', k, ')'
      end if
      
      allocate (readsrc(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'readsrc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': readsrc allocated (', n, ',', m, ')'
      end if
      
      allocate (rel_depth_phase(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rel_depth_phase', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': rel_depth_phase allocated (', n, ',', m, ')'
      end if
      
      allocate (rel_phase(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rel_phase', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': rel_phase allocated (', n, ',', m, ')'
      end if
      
      allocate (resisc(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'resisc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': resisc allocated (', n, ',', m, ')'
      end if
      
      allocate (sad(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sad', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sad allocated (', n, ',', m, ')'
      end if
      
      allocate (sdread(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sdread', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sdread allocated (', n, ',', m, ')'
      end if
      
      allocate (stname(n,0:m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stname', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': stname allocated (', n, ',0:', m, ')'
      end if
      
      allocate (ttoff(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ttoff', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ttoff allocated (', n, ',', m, ')'
      end if
      
      allocate (ttsprd(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ttsprd', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ttsprd allocated (', n, ',', m, ')'
      end if
      
      allocate (weight(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'weight', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': weight allocated (', n, ',', m, ')'
      end if
      
      allocate (z_hrdp(n,l), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'z_hrdp', error)
      else
         z_hrdp = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': z_hrdp allocated (', n, ',', l, ')'
      end if
      
      allocate (z_test(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'z_test', error)
      else
         z_test = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': z_test allocated (', n, ')'
      end if
     
      return
   
   end subroutine phase_data_allocate


!*****************************************************************************************
   subroutine stations_allocate ()
   
   ! Allocate variable arrays related to stations
   
      integer :: it
      integer :: it1
      integer :: l
      integer :: m
      integer :: n
      
      l = n_miss_sta_max
      n = n_event ! number of events in the cluster
      m = n_arrtim_max ! maximum number of phase records for any event in the cluster
      it = n_it_max ! maximum number of iterations
      it1 = n_it_max1
      p = 'stations_allocate' ! procedure

      allocate (ahgtr(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ahgtr', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': ahgtr allocated (', n_sta_max, ')'
      end if
      
      allocate (ahgts(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ahgts', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ahgts allocated (', n, ',', m, ')'
      end if
      
      allocate (azes(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'azes', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': azes allocated (', n, ',', m, ')'
      end if
      
      allocate (azse(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'azse', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': azse allocated (', n, ',', m, ')'
      end if
      
      allocate (delt(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'delt', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': delt allocated (', n, ',', m, ')'
      end if
      
      allocate (dt(n,m,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dt', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a,i1,a)') trim(p)//': dt allocated (', n, ',', m, ',0:', it, ')'
      end if
      
      allocate (elcr(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'elcr', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': elcr allocated (', n, ',', m, ')'
      end if
      
      allocate (jdate_off(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'jdate_off', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': jdate_off allocated (', n_sta_max, ')'
      end if
      
      allocate (jdate_on(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'jdate_on', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': jdate_on allocated (', n_sta_max, ')'
      end if
      
      allocate (kcode(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'kcode', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': kcode allocated (', n, ',', m, ')'
      end if
      
      allocate (kstat(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'kstat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': kstat allocated (', n_sta_max, ')'
      end if
      
      allocate (missta(n,l), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'missta', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i3,a)') trim(p)//': missta allocated (', n, ',', l, ')'
      end if
      
      allocate (ndat(n,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ndat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': ndat allocated (', n, ',0:', it1, ')'
      end if
      
      allocate (ndatc(n,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ndatc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': ndatc allocated (', n, ',0:', it1, ')'
      end if
      
      allocate (ndatdl(n,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ndatdl', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': ndatdl allocated (', n, ',0:', it1, ')'
      end if
      
      allocate (ndatfl(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ndatfl', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': ndatfl allocated (', n, ')'
      end if
      
      allocate (ndatpr(n,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ndatpr', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': ndatpr allocated (', n, ',0:', it1, ')'
      end if
      
      allocate (nmiss(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'nmiss', error)
      else
         nmiss = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': nmiss allocated (', n, ')'
      end if
      
      allocate (nstr1(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'nstr1', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': nstr1 allocated (', n_sta_max, ')'
      end if
      
      allocate (psd(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'psd', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': psd allocated (', n, ',', m, ')'
      end if
      
      allocate (s(n,m,0:it), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 's', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a,i1,a)') trim(p)//': s allocated (', n, ',', m, ',0:', it, ')'
      end if
      
      allocate (sca0(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sca0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sca0 allocated (', n_sta_max, ')'
      end if
      
      allocate (sca0s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sca0s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sca0s allocated (', n, ',', m, ')'
      end if
      
      allocate (sca1s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sca1s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sca1s allocated (', n, ',', m, ')'
      end if
      
      allocate (sca2s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sca2s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sca2s allocated (', n, ',', m, ')'
      end if
      
      allocate (scb1s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'scb1s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': scb1s allocated (', n, ',', m, ')'
      end if
      
      allocate (scb2s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'scb2s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': scb2s allocated (', n, ',', m, ')'
      end if
      
      allocate (sd1(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sd1', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sd1 allocated (', n_sta_max, ')'
      end if
      
      allocate (sd1s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sd1s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sd1s allocated (', n, ',', m, ')'
      end if
      
      allocate (sd2s(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sd2s', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sd2s allocated (', n, ',', m, ')'
      end if
      
      allocate (sdstcs(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sdstcs', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': sdstcs allocated (', n, ',', m, ')'
      end if
      
      allocate (sta_agency(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sta_agency', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sta_agency allocated (', n_sta_max, ')'
      end if
      
      allocate (sta_author(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sta_author', error)
      else
         sta_author = ' '
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sta_author allocated (', n_sta_max, ')'
      end if
      
      allocate (sta_cha(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sta_cha', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sta_cha allocated (', n_sta_max, ')'
      end if
      
      allocate (sta_deployment(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sta_deployment', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sta_deployment allocated (', n_sta_max, ')'
      end if
      
      allocate (sta_loc(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sta_loc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': sta_loc allocated (', n_sta_max, ')'
      end if
      
      allocate (stacs(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stacs', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': stacs allocated (', n, ',', m, ')'
      end if
      
      allocate (stalat(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stalat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': stalat allocated (', n_sta_max, ')'
      end if
      
      allocate (stalon(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stalon', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': stalon allocated (', n_sta_max, ')'
      end if
      
      allocate (stladg(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stladg', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': stladg allocated (', n, ',', m, ')'
      end if
      
      allocate (stlndg(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stlndg', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': stlndg allocated (', n, ',', m, ')'
      end if
      
      allocate (stn_dcal_used(n_sta_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'stn_dcal_used', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': stn_dcal_used allocated (', n_sta_max, ')'
      end if
      
      allocate (ttcomp(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ttcomp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': ttcomp allocated (', n, ',', m, ')'
      end if
      
      allocate (tto(n,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tto', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': tto allocated (', n, ',', m, ')'
      end if
      
      return
   
   end subroutine stations_allocate


!*****************************************************************************************
   subroutine plot_allocate ()
   
   ! Allocate variable arrays related to plotting
   
   integer :: i
   integer :: m
   integer :: n
   
   n = n_event
   m = n_star_max
   p = 'plot_allocate'
   
      allocate (evrp_delta_max(n_evrp), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'evrp_delta_max', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': evrp_delta_max allocated (', n_evrp, ')'
      end if
      
      allocate (evrp_delta_min(n_evrp), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'evrp_delta_min', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': evrp_delta_min allocated (', n_evrp, ')'
      end if
      
      allocate (evrp_event(n_evrp), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'evrp_event', error)
      else
         evrp_event = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': evrp_event allocated (', n_evrp, ')'
      end if
      
      allocate (iev_star(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'iev_star', error)
      else
         iev_star = 0
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': iev_star allocated (', m, ')'
      end if
      
      allocate (pslp_color(n_pslp), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'pslp_color', error)
      else
         pslp_color = ' '
         if (verbose_log) write (io_alloc_log,'(a,i2,a)') trim(p)//': pslp_color allocated (', n_pslp, ')'
      end if
      
      allocate (pslp_filnam(n_pslp), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'pslp_filnam', error)
      else
         pslp_filnam = ' '
         if (verbose_log) write (io_alloc_log,'(a,i2,a)') trim(p)//': pslp_filnam allocated (', n_pslp, ')'
      end if
      
      allocate (plot(n,0:n_selected_max), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'plot', error)
      else
         plot = .false.
         do i = 1,n_event
            plot(i,0) = .true.
         end do
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': plot allocated (', n, ',0:', n_selected_max, ')'
      end if
      
      allocate (rdpp_evt(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'rdpp_evt', error)
      else
         rdpp_evt = ' '
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': rdpp_evt allocated (', n, ')'
      end if
      
      allocate (star_size(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'star_size', error)
      else
         star_size = 0.
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': star_size allocated (', m, ')'
      end if
      
      allocate (tt5e_evt(n_tt5e), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tt5e_evt', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': tt5e_evt allocated (', n_tt5e, ')'
      end if
      
      return
   
   end subroutine plot_allocate
   
   
!*****************************************************************************************   
   subroutine hdc_allocate ()
   
   integer :: it1
   integer :: j
   integer :: k
   integer :: l
   integer :: m
   integer :: mt
   integer :: n
   
   it1 = n_it_max1
   j = n_event*n_q_max
   k = n_arrtim_max
   l = n_qi_max
   m = n_q_max
   mt = n_fp_max
   n = n_event
   p = 'hdc_allocate'

      allocate (a(n,k,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'a', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': a allocated (', n, ',', k, ',4)'
      end if
      
      allocate (dimpciev(n,12), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dimpciev', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': dimpciev allocated (', n, ',12)'
      end if
      
      allocate (dtmpc(n,k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dtmpc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': dtmpc allocated (', n, ',', k, ')'
      end if
      
      allocate (dtmph(n,k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dtmph', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': dtmph allocated (', n, ',', k, ')'
      end if
      
      allocate (dxp(n,4,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dxp', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': dxp allocated (', n, ',4,0:', it1, ')'
      end if
      
      allocate (eci(n,k), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'eci', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': eci allocated (', n, ',', k, ')'
      end if
      
      allocate (eciev(n,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'eciev', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': eciev allocated (', n, ',0:', it1, ')'
      end if
      
      allocate (e_dsvd(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'e_dsvd', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': e_dsvd allocated (', mt, ')'
      end if
      
      allocate (ehiev(n,0:it1), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ehiev', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i1,a)') trim(p)//': ehiev allocated (', n, ',0:', it1, ')'
      end if
      
      allocate (jb(j), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'jb', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i8,a)') trim(p)//': jb allocated (', j, ')'
      end if
      
      allocate (ntq(l,m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ntq', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i3,a,i4,a)') trim(p)//': ntq allocated (', l, ',', m, ')'
      end if
      
      allocate (ntqi(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ntqi', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': ntqi allocated (', m, ')'
      end if
      
      allocate (sdxhatc(n,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sdxhatc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': sdxhatc allocated (', n, ',4)'
      end if
      
      allocate (shatsqci(n), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'shatsqci', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': shatsqci allocated (', n, ')'
      end if
      
      allocate (sighatj(j), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sighatj', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i8,a)') trim(p)//': sighatj allocated (', j, ')'
      end if
      
      allocate (vhatc(mt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'vhatc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': vhatc allocated (', mt, ',', mt, ')'
      end if
      
      allocate (vnhat(j), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'vnhat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i8,a)') trim(p)//': vnhat allocated (', j, ')'
      end if
      
      allocate (wq(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'wq', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': wq allocated (', m, ')'
      end if
      
      allocate (wq2(m), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'wq2', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': wq2 allocated (', m, ')'
      end if
      
      return
      
   end subroutine hdc_allocate

!*****************************************************************************************
end module mloc_allocate

