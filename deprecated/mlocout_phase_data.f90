!> Procedures for writing output files related to phase data

module mlocout_phase_data

   use declare_limits
   use declare_calibration
   use declare_cluster_vectors
   use declare_configuration
   use declare_environment
   use declare_events
   use declare_hypocentroid
   use declare_lun
   use declare_output
   use declare_phase_data
   use declare_phases
   use declare_stations
   use mloclib_stations
      
   implicit none
   save
   
contains

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
      integer, dimension(npmax) :: indx
      integer :: io_hypo_cal
      integer :: io_hypo_dcal
      integer :: io_out_cal
      integer, intent(in) :: it ! iteration #
      integer :: it1
      integer :: j
      integer, external :: lunit
      integer :: n_res
      real, dimension(npmax) :: deltiev
      real :: deltkm
      real, dimension(0:itmax1) :: dts
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
      logical :: bdp ! logical function in mloclib_stations.f90
      
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
   
      do iev = 1,nev ! Loop over events for good data
            
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
                     if (indirect_cal) call write_good_phase_line (io_hypo_cal, 0, iev, i, gfmt, ' ', .true., .true.)
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
      
      do iev = 1,nev ! Loop over events
      
         ! How depth was set
         call depth_set_text (iev, depth_set)
      
         ! Event separator
         write (io_out,'(1x,164a1)') ('*',i=1,164)
         if (indirect_cal) write (io_out_cal,'(1x,164a1)') ('*',i=1,164)
         
         ! Heading for each event
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
      real, dimension(0:itmax1) :: dtj
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
      real :: dtj(0:itmax1)
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
      integer :: n_input
      real :: cluster_default_depth
      real, dimension(nevmax) :: x_dep_constrained
      real, dimension(nevmax) :: x_dep_input
      real :: xmed
      
      n_constrained = 0
      n_input = 0
      cluster_default_depth = -99.
      do iev = 1, nev
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
      if (nev .ge. 3) then
         call mdian1 (x_dep_input, nev, xmed)
         write (io_unit,'(a,f6.1)') 'Median of input file depths: ', xmed
      end if
      
      return
      
   end subroutine depth_statistics


!*****************************************************************************************
end module mlocout_phase_data

