!*****************************************************************************************
program mloc

! Main program and command processor to perform multiple event relocation on a cluster of
! earthquakes, using the Hypocentroidal Decomposition (HD) algorithm of Jordan
! and Sverdrup (BSSA, 71, 1105-1130, 1981). The code has been heavily modified
! to support research in calibrated locations, but the core HD algorithm is unchanged.

! See the file 'mloc version history' for documentation on major changes.
! There are also extensive comments in the code.

! Version number and release date are carried in "declare_environment.f90".

! Dr. Eric A. Bergman
! Global Seismological Services
! 1900 19th St., Golden, CO 80401, USA
! +1 (720) 400-7835
! bergman@seismo.com
! http://www.seismo.com

   use mloc_allocate
   use mloc_calibration
   use mloc_commands
   use mloc_declare
   use mloc_init
   use mloc_lib
   use mloc_set
   use mloc_stations

   implicit none
   save

   integer :: i
   integer :: iev = 0
   integer :: io_bptc
   integer :: io_temp
   integer :: ios
   integer :: j
   integer :: jev
   integer :: luopen
   integer :: n_hdf_lines
   integer :: npiev
   real :: cpu_finish
   real :: cpu_start
   character(len=100) :: alloc_log_file ! log file for array allocation
   character(len=32) :: bp_file = 'etopo5_bed_g_i2.txt' ! topography data file for bounce point time correction
   character(len=32) :: bp_path = 'tables/tau-p' ! path to bp_file
   character(len=100) :: bptc_log_file ! bounce point time correction log file
   character(len=100) :: command_line
   character(len=100) :: filename
   character(len=32) :: ifsm_file = 'ISC-FDSN_stn.dat' ! ISC-FDSN station metadata file
   character(len=32) :: junk_phases_file = 'junk_phases.dat' ! list of phase names to ignore
   character(len=80) :: line = ' '
   character(len=80) :: line_temp = ' '
   character(len=12) :: line12
   character(len=100) :: line132
   character(len=34) :: line34
   character(len=8) :: line8
   character(len=80) :: logfile
   character(len=132) :: msg
   character(len=32) :: nsmd_file = 'neic_stn.dat' ! NEIC station metadata file
   character(len=100) :: open_file_name
   character(len=100) :: outfil
   character(len=32) :: phases_path = 'tables/phases' ! path to phase-related files
   character(len=100) :: plog_file ! phase re-identification log file
   character(len=32) :: psdre_file = 'psdre.dat' ! Phase-specific default reading errors
   character(len=32) :: psdre_path = 'tables/spread' ! path to psdre.dat
   character(len=4) :: record_type
   character(len=100) :: stn_log_file ! station data log file
   character(len=100) :: tlog_file ! log file for travel-time calculations
   logical :: command_help
   logical :: ex
   logical :: op

   call cpu_time (cpu_start)

   ! Initialization
   call initialize ()
   
   write (*,'(/2a/)') trim(mloc_version), '; status: '//trim(mloc_status)
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
   ! Configuration file
   call open_file (1, 'mloc.conf', 'old')
   do
      line132 = ' '
      line = ' '
      read (1,'(a)',iostat=ios) line132
      if (ios .lt. 0) exit
      j = len_trim(line132)
      if (line132(1:13) .eq. 'WORKING_DIR: ') then
         mloc_path = line132(14:j) ! Absolute pathname to the directory containing the mloc executable
      else if (line132(1:16) .eq. 'STATION_MASTER: ') then
         station_master = line132(17:j) ! master station list
      else if (line132(1:9) .eq. 'GMT_VER: ') then
         read (line132(10:10),'(i1)') gmt_version
         if (gmt_version .eq. 6) then
            write (*,'(a/11x,a)') 'mloc_main: GMT v6.5 is incompatible with this version of mloc',&
             'GMT v6.4 or lower is recommended'
            cycle
         else if (gmt_version .eq. 5) then
            call fyi ('mloc_main: GMT v5 is still supported but it would be advisable to upgrade to v6 soon.')
         else
            write (msg,'(a,i1,a)') 'mloc_main: illegal value (', gmt_version, ') for GMT version (GMT v5 or v6 is required)'
            call oops (trim(msg))
         end if
      else if (line132(1:8) .eq. 'AUTHOR: ') then
         line = line132(9:j) ! ID for the person running mloc (maximum 8 characters)
         if (len_trim(line) .le. 8) then
            mloc_author = trim(line)
         else
            mloc_author = line(1:8)
            write (msg,'(a,a8)') 'mloc_main: "mloc_author" truncated to 8 characters: ', mloc_author
            call warnings (trim(msg))
         end if
      else if (line132(1:7) .eq. 'SHELL: ') then
         line = line132(8:j)
         if (line(1:4) .eq. 'bash') then
            ishell = 1
         else if (line(1:3) .eq. 'csh') then
            ishell = 2
         else if (line(1:3) .eq. 'zsh') then
            ishell = 3
         else
            write (msg,'(2a)') 'mloc_main: unknown argument for the SHELL parameter in mloc.conf: ', trim(line132)
            call oops (trim(msg))
         end if
      else if (line132(1:8) .eq. 'SAMPLE: ') then
         call oops ('mloc_main: mloc.conf must be edited prior to running mloc!') 
      else
         write (msg,'(2a)') 'mloc_main: unknown keyword in mloc.conf: ', trim(line132)
         call oops (trim(msg))
      end if
   end do
   close (1)
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! Interactive start
   write (*,'(/a)',advance='no') 'Enter a basename for this run: '
   read (*,'(a)') basename
   write (*,'(a)',advance='no') 'Enter the name of the data directory: '
   read (*,'(a)') datadir
   outfile = trim(datadir)//dirsym//trim(basename)
   
   ! Check if this basename already exists, open log file
   logfile = trim(outfile)//'.log'
   call open_file (io_log, trim(logfile), 'new')

   ! Control of mloc is based on the use of 'commands' which can be entered interactively or
   ! through a command_file. Commands have three or four letter codes, issued in lowercase.
   ! Some commands require additional parameters; some can take additional parameters as
   ! options; some take no additional parameters.
   
   ! A command_file is required to be read before any interactive commands are issued but
   ! it is only required to contain the basic event definition blocks:
   !   memb
   !   even {event name}
   !   inpu {event data filename}
   ! These commands (and KILL) can only be issued from the command_file. The same is true
   ! for some other commands as well (e.g., RHDF, CAL_, DCAL) because associated arrays
   ! must be allocated before the command file and interactive commands are fully processed. 
   ! It is recommended to put most of the intended commands in the command_file to prevent
   ! errors.
   
   ! Get the command_file name
   write (*,'(/a)', advance='no') 'Enter a command_file name (or "help"): '
   read (*,'(a)') line
   if (line(1:4) .eq. 'cfil') then ! catch old-style cfil command
      line_temp = line(6:lenb(line))
      line = line_temp
      call warnings ('mloc_main: command "cfil" is no longer needed, just the filename')
   end if
   if (trim(line) .eq. 'help') then
      call proc_cfil (trim(line), .true.)
   else
      call proc_cfil (trim(line), .false.)
   end if
   
   ! Process the command file for number of events, calibration commands, number of phase
   ! records in each event, number of entries in an hdf file to be read, verbose_log, 
   ! verbose_screen and various limits nneded for allocation.
   n_star = 0
   kill_count = 0
   do
      read (io_cfil,'(a)',iostat=ios) line
      if (ios .lt. 0) then
         rewind (io_cfil)
         exit
      else
         call parse (line)
         
         ! Skip all commands except 'kill' when doing a block kill
         if (kill_all .and. comd .ne. 'kill') cycle
      
         ! Skip all commands except 'memb' and 'kill' when killing a single event with the 'memb' command
         if (kill_one .and. comd .ne. 'memb' .and. comd .ne. 'kill') cycle
   
         if (comd(1:3) .eq. 'cal') then
            indirect_cal = .true.
         else if (comd .eq. 'dcal' .and. params(1:2) .eq. 'on') then
            direct_cal = .true.
         else if (comd .eq. 'evrp') then
            evrp = .true.
            if (params(1:3) .eq. 'all') then ! "all" option is handled after the command file has been read
               evrp_all = .true.
            else
               n_evrp = n_evrp + 1
            end if
         else if (comd .eq. 'inpu') then
            if (verbose_screen) write (*,'(2a)') 'main: ', trim(line)
            call open_file (io_in, trim(datadir)//dirsym//trim(params), 'old')
            npiev = 0
            do
               read (io_in,'(a4)',iostat=ios) record_type
               if (ios .lt. 0) exit
               if (record_type(1:1) .eq. 'P') npiev = npiev + 1 ! Count flagged readings too
            end do
            close (io_in)
            n_arrtim_max = max(npiev, n_arrtim_max)
         else if (comd .eq. 'kill') then
            call proc_kill (.false.)
         else if (comd .eq. 'memb') then
            call proc_memb (iev, .false.)
         else if (comd .eq. 'pslp') then
            n_pslp = n_pslp + 1
         else if (comd .eq. 'rhdf') then
            call open_file (io_rhdf, trim(datadir)//dirsym//trim(params), 'old')
            n_hdf_lines = 0
            do
               read (io_rhdf,'(a)',iostat=ios) line
               if (ios .lt. 0) exit
               n_hdf_lines = n_hdf_lines + 1
            end do
            close (io_rhdf)
         else if (comd .eq. 'star') then
            n_star = n_star + 1
            star_plot = .true.
            if (params(1:4) .eq. 'auto') star_auto = .true.
         else if (comd .eq. 'tt5e') then
            tt5e = .true.
            if (params(1:3) .eq. 'all') then ! "all" option is handled after the command file has been read
               tt5e_all = .true.
            else
               n_tt5e = n_tt5e + 1
            end if
         else if (comd .eq. 'vlog' .and. params(1:2) .eq. 'on') then
            verbose_log = .true.
         else if (comd .eq. 'vscr' .and. params(1:2) .eq. 'on') then
            verbose_screen = .true.
         end if
      end if
   end do
   
   if (kill_one) kill_one = .false.
   if (kill_all) kill_all = .false.
   
   n_event = iev

   ! Update limits derived from number of events
   n_fp_max = n_event*4 ! max number of free parameters
   n_hdf_max = max(n_hdf_lines,n_event) ! in case no hdf is read
   if (star_plot) then
      if (star_auto) then
         n_star_max = n_event
      else
         n_star_max = n_star
      end if
   end if
   n_qi_max = n_event*4 ! max number of instances of a particular station-phase
   ! This is a guess. It may need to be increased for some datasets. It is only used for
   ! one array (ntq) so it's not worth the trouble to estimate accurately.
   if (evrp_all) n_evrp = n_event
   if (tt5e_all) n_tt5e = n_event
      
   write (io_log,'(a)') 'mloc_main: calculated limits:'
   write (io_log,'(2x,a,i4)') 'mloc_main: number of events = ', n_event
   write (io_log,'(2x,a,i5)') 'mloc_main: max # arrivals for any event = ', n_arrtim_max
   write (io_log,'(2x,a,i4)') 'mloc_main: # hdf file entries to read = ', n_hdf_max
   
   write (io_log,'(/a)') 'mloc_main: derived limits:'
   write (io_log,'(2x,a,i4)') 'mloc_main: max # free parameters = ', n_fp_max
   write (io_log,'(2x,a,i4)') 'mloc_main: max # distinct station-phases (default) = ', n_q_max
   write (io_log,'(2x,a,i4)') 'mloc_main: max # instances of a particular station-phase = ', n_qi_max
   if (star_plot) write (io_log,'(2x,a,i4/)') 'mloc_main: max # stars to plot = ', n_star_max
   if (evrp) write (io_log,'(2x,a,i4)') 'mloc_main: # event residual plots = ', n_evrp
   if (tt5e) write (io_log,'(2x,a,i4)') 'mloc_main: # tt5e plots = ', n_tt5e
   
   ! Allocation
   ! Logging of allocation statements is done when verbose_log = .true.
   if (verbose_log) then
      alloc_log_file = trim(outfile)//'.alloc'
      call open_file (io_alloc_log, alloc_log_file, 'new')
   end if
   
   call events_allocate ()
   call hdf_allocate ()
   call plot_allocate ()
   if (direct_cal) then
      call direct_cal_allocate ()
      call direct_calibration_allocate ()
   end if
   if (indirect_cal) then
      call indirect_cal_allocate ()
      call indirect_calibration_allocate ()
   end if
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Re-read the command_file, processing all commands and take interactive commands to set
!  up the relocation
   write (*,'(/a,21(/t3,a))') 'The commands are:',&
    'ANNO',&
    'bdps bias bloc bptc',&
    'CAL  clim comm corr cptf ctyp cvff cvtt',&
    'damp datf dbug DCAL dem1 dem2 DEP  diff',&
    'ellp epap eplt EVEN evrp',&
    'fdhp flag fmap frec freh',&
    'gcat gmtp',&
    'help hlim',&
    'ifsm INPU',&
    'KILL',&
    'LAT  lgtt lmod LONG lonr lres',&
    'mare MEMB',&
    'noxf nsmd',&
    'oldr',&
    'p_tt pert phid phrp phyp plog PLOT pltt ppri pslp pttt',&
    'radf rdpp rels revi rfil rgtt RHDF rpth run ',&
    'shcl skip splt sstn STAR stat step stop subc',&
    'taup terr test tfil tikh TIME tlog tomo tptt tt5e tt5s ttou',&
    'vect VLOG vscr',&
    'weig wind',&
    'xsec'
   
   write (*,'(/a)') 'Commands in all caps can only be issued in a command file.'

   write (*,'(/a/)') 'Processing the command_file'
   
   iev = 0
   n_pslp = 0
   n_star = 0
   kill_count = 0
   n_evrp = 0
   n_tt5e = 0
   do ! Loop for commands, both command_file and interactive

      ! Read a line from a command file or interactive input
      if (cmndfil) then ! Reading from a command file
         read (io_cfil,'(a)',iostat=ios) line
         if (ios .ge. 0) then
            if (line(1:1) .eq. '!') cycle ! alternative to the COMM command
            if (line(1:1) .eq. '#') cycle ! alternative to the COMM command
            call parse (line)
            if (verbose_log) write (io_log,'(a,a,1x,a)') 'main: ', trim(comd), trim(params)
         else
            close (io_cfil)
            cmndfil = .false. ! Stop reading from the command_file and take interactive commands
         
            write (msg,'(a,i3,a)') 'mloc_main: ', n_event, ' events will be relocated'
            call fyi (trim(msg))
            
            ! Summary of killed events
            if (kill_count .gt. 0) then
               write (msg,'(a,i3,a)') 'mloc_main: ', kill_count, ' events from the command file were killed'
               call fyi (trim(msg))
            end if
         
            ! If the last event in the command file was killed by memb/kill, turn killing off now
            if (kill_one) kill_one = .false.
         
            write (*,'(/a)') 'End of command_file reached, continue with interactive commands'
            write (*,'(a)') 'For more information, follow the "help" command with a command name'
            write (*,'(/a/)') 'Enter a command:'
         
            ! Start interactive input
            write (*,'(a)',advance='no') ': '
            read (*,'(a)') line
            call parse (line)
         end if
      else ! Interactive input
         write (*,'(a)',advance='no') ': '
         read (*,'(a)') line
         call parse (line)
      end if
      
      ! Log the "inpu" argument for killed events
      if ((kill_one .or. kill_all) .and. comd .eq. 'inpu') then ! Keep track of killed events
         kill_count = kill_count + 1
         if (kill_one) then
            write (io_log,'(2a,1x,a)') 'killed ', trim(params), trim(kill_reason)
         else if (kill_all) then
            write (io_log,'(2a,1x,a)') 'killed ', trim(params), 'block kill'
         end if
      end if
      
      ! Skip all commands except 'kill' when doing a block kill
      if (kill_all .and. comd .ne. 'kill') cycle
   
      ! Skip all commands except 'memb' and 'kill' when killing a single event with the 'memb' command
      if (kill_one .and. comd .ne. 'memb' .and. comd .ne. 'kill') cycle
   
      ! Help system for specific commands
      command_help = .false.
      if (comd .eq. 'help' .and. params .ne. ' ') then
         comd = params(1:4)
         command_help = .true.
      end if
      
      ! Process commands
      if (comd .eq. 'anno') then
         call proc_anno (iev, command_help)
      else if (comd .eq. 'bdps') then
         call proc_bdps (command_help)
      else if (comd .eq. 'bias') then
         call proc_bias (command_help)
      else if (comd .eq. 'bloc') then
         call proc_bloc (command_help)
      else if (comd .eq. 'bptc') then
         call proc_bptc (command_help)
      else if (comd(1:3) .eq. 'cal') then
         if (.not.allocated(cal_event)) call indirect_cal_allocate ()
         call proc_cal (iev, command_help)
      else if (comd .eq. 'ccat') then
         call warnings ('mloc_main: "ccat" has been replaced by "gcat"')
      else if (comd .eq. 'clim') then
         call proc_clim (command_help)
      else if (comd .eq. 'comm') then
         call proc_comm (command_help)
      else if (comd .eq. 'corr') then
         call proc_corr (command_help)
      else if (comd .eq. 'cptf') then
         call proc_cptf (command_help)
      else if (comd .eq. 'ctyp') then
         call proc_ctyp (command_help)
      else if (comd .eq. 'cvff') then
         call proc_cvff (command_help)
      else if (comd .eq. 'cvtt') then
         call proc_cvtt (command_help)
      else if (comd .eq. 'damp') then
         call proc_damp (command_help)
      else if (comd .eq. 'datf') then
         call proc_datf (command_help)
      else if (comd .eq. 'dbug') then
         call proc_dbug (command_help)
      else if (comd .eq. 'dcal') then
         call proc_dcal (command_help)
      else if (comd .eq. 'dem1') then
         call proc_dem1 (command_help)
      else if (comd .eq. 'dem2') then
         call proc_dem2 (command_help)
      else if (comd(1:3) .eq. 'dep') then
         call proc_dep (iev, command_help)
      else if (comd .eq. 'diff') then
         call proc_diff (command_help)
      else if (comd .eq. 'ellp') then
         call proc_ellp (command_help)
      else if (comd .eq. 'epap') then
         call proc_epap (command_help)
      else if (comd .eq. 'eplt') then
         call proc_eplt (command_help)
      else if (comd .eq. 'even') then
         call proc_even (iev, command_help)
      else if (comd .eq. 'evrp') then
         call proc_evrp (command_help)
      else if (comd .eq. 'fdhp') then
         call proc_fdhp (command_help)
      else if (comd .eq. 'flag') then
         call proc_flag (command_help)
      else if (comd .eq. 'fmap') then
         call proc_fmap (command_help)
      else if (comd .eq. 'frec') then
         call proc_frec (iev, command_help)
      else if (comd .eq. 'freh') then
         call proc_freh (command_help)
      else if (comd .eq. 'gcat') then
         call proc_gcat (command_help)
      else if (comd .eq. 'gmtp') then
         call proc_gmtp (command_help)
      else if (comd .eq. 'help') then
         call proc_help (command_help)
      else if (comd .eq. 'hlim') then
         call proc_hlim (command_help)
      else if (comd .eq. 'ifsm') then
         call proc_ifsm (command_help)
      else if (comd .eq. 'inpu') then
         call proc_inpu (iev, command_help)
      else if (comd .eq. 'kill') then
         call proc_kill (command_help)
      else if (comd(1:3) .eq. 'lat') then
         call proc_lat (iev, command_help)
      else if (comd .eq. 'lgou') then
         call warnings ('mloc_main: "lgou" has been replaced by "ttou"')
      else if (comd .eq. 'lgtt') then
         call proc_lgtt (command_help)
      else if (comd .eq. 'lmod') then
         call proc_lmod (command_help)
      else if (comd .eq. 'long') then
         call proc_long (iev, command_help)
      else if (comd .eq. 'lonr') then
         call proc_lonr (command_help)
      else if (comd .eq. 'lres') then
         call proc_lres (command_help)
      else if (comd .eq. 'mare') then
         call proc_mare (command_help)
      else if (comd .eq. 'memb') then
         call proc_memb (iev, command_help)
      else if (comd .eq. 'noxf') then
         call proc_noxf (command_help)
      else if (comd .eq. 'nsmd') then
         call proc_nsmd (command_help)
      else if (comd .eq. 'oldr') then
         call proc_oldr (command_help)
      else if (comd .eq. 'outp') then
         call warnings ('mloc main: command "outp" is no longer needed')
      else if (comd .eq. 'p_tt') then
         call proc_p_tt (command_help)
      else if (comd .eq. 'pert') then
         call proc_pert (command_help)
      else if (comd .eq. 'phid') then
         call proc_phid (command_help)
      else if (comd .eq. 'phrp') then
         call proc_phrp (command_help)
      else if (comd .eq. 'phyp') then
         call proc_phyp (command_help)
      else if (comd .eq. 'plog') then
         call proc_plog (command_help)
      else if (comd .eq. 'plot') then
         call proc_plot (iev, command_help)
      else if (comd .eq. 'pltt') then
         call proc_pltt (command_help)
      else if (comd .eq. 'ppri') then
         call proc_ppri (command_help)
      else if (comd .eq. 'pslp') then
         call proc_pslp (command_help)
      else if (comd .eq. 'pttt') then
         call proc_pttt (command_help)
      else if (comd .eq. 'radf') then
         call proc_radf (command_help)
      else if (comd .eq. 'rdpp') then
         call proc_rdpp (command_help)
      else if (comd .eq. 'rels') then
         call proc_rels (command_help)
      else if (comd .eq. 'resp') then
         call warnings ('mloc main: "resp" has been replaced by "evrp"')
         call proc_evrp (command_help)
      else if (comd .eq. 'revi') then
         call proc_revi (command_help)
      else if (comd .eq. 'rfil') then
         call proc_rfil (command_help)
      else if (comd .eq. 'rgtt') then
         call proc_rgtt (command_help)
      else if (comd .eq. 'rhdf') then
         call proc_rhdf (command_help)
      else if (comd .eq. 'rpth') then
         call proc_rpth (command_help)
      else if (comd(1:3) .eq. 'run') then
         call proc_run (command_help)
         if (.not.command_help) exit ! Leave the interactive command loop and start the run
      else if (comd .eq. 'shcl') then
         call proc_shcl (command_help)
      else if (comd .eq. 'skip') then
         call proc_skip (command_help)
      else if (comd .eq. 'splt') then
         call proc_splt (command_help)
      else if (comd .eq. 'sstn') then
         call proc_sstn (command_help)
      else if (comd .eq. 'star') then
         call proc_star (iev, command_help)
      else if (comd .eq. 'stat') then
         call proc_stat (command_help)
      else if (comd .eq. 'step') then
         call proc_step (command_help)
      else if (comd .eq. 'subc') then
         call proc_subc (command_help)
      else if (comd .eq. 'stop') then
         call proc_stop (command_help)
      else if (comd .eq. 'taup') then
         call proc_taup (command_help)
      else if (comd .eq. 'terr') then
         call proc_terr (command_help)
      else if (comd .eq. 'test') then
         call proc_test (command_help)
      else if (comd .eq. 'tfil') then
         call proc_tfil (command_help)
      else if (comd .eq. 'tikh') then
         call proc_tikh (command_help)
      else if (comd .eq. 'time') then
         call proc_time (iev, command_help)
      else if (comd .eq. 'tlog') then
         call proc_tlog (command_help)
      else if (comd .eq. 'tomo') then
         call proc_tomo (command_help)
      else if (comd .eq. 'topo') then
         call warnings ('mloc_main: "topo" has been replaced by "dem1"')
         call proc_dem1 (command_help)
      else if (comd .eq. 'tpou') then
         call warnings ('mloc_main: "tpou" has been replaced by "ttou"')
      else if (comd .eq. 'tptt') then
         call proc_tptt (command_help)
      else if (comd .eq. 'tt5e') then
         call proc_tt5e (command_help)
      else if (comd .eq. 'tt5s') then
         call proc_tt5s (command_help)
      else if (comd .eq. 'ttou') then
         call proc_ttou (command_help)
      else if (comd .eq. 'vect') then
         call proc_vect (command_help)
      else if (comd .eq. 'vlog') then
         call proc_vlog (command_help)
      else if (comd .eq. 'vscr') then
         call proc_vscr (command_help)
      else if (comd .eq. 'weig') then
         call proc_weig (command_help)
      else if (comd .eq. 'wind') then
         call proc_wind (command_help)
      else if (comd .eq. 'xsec') then
         call proc_xsec (command_help)
      else
         call warnings ('mloc_main: '//comd//' not found')
      end if
   
   end do
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! Preparations for the run
   
   ! Allocate arrays
   call phase_data_allocate ()
   call stations_allocate ()
   call hdc_allocate ()
   
   ! Turn off ~.lres output file for a single event
   if (n_event .eq. 1) lresout = .false.

   ! Make a directory for the GMT scripts
   gmt_script_dir = trim(datadir)//dirsym//trim(basename)//'_gmt_scripts'
   command_line = 'mkdir '//trim(gmt_script_dir)
   call execute_command_line (trim(command_line))

   ! Check for necessary file specifications
   do jev = 1,n_event
      if (infile(jev)(1:1) .eq. ' ') then
         write (msg,'(a,i3)') 'mloc_main: no input file specified for event ', jev
         call oops (trim(msg))
      end if
   end do
   
   ! Open the LRES file
   if (lresout) then
      outfil = trim(outfile)//'.lres'
      call open_file (io_lres, outfil, 'new')
      write (io_lres,'(f6.2)') lres
   end if
   
   ! Log convergence limits
   write (io_log,'(/a)') 'Convergence limits'
   write (io_log,'(a,f6.3,a,2(f4.1,a))') 'Hypocentroid: ', cl_epi_h, ' (deg), ', cl_dep_h, '(km), ', cl_ot_h, ' (sec)'
   write (io_log,'(a,3(f4.1,a))') 'Cluster vectors: ', cl_epi_c, ' (km),  ', cl_dep_c, '(km), ', cl_ot_c, ' (sec)'
   
   ! Log the local velocity model
   if (locmod) then 
      call open_file (io_locmod, trim(locmodfname), 'old')
      write (io_log,'(/2a)') 'Local velocity model: ', trim(locmodfname)
      ! Epicentral range line
      read (io_locmod,'(a)',iostat=ios) line34
      write (io_log,'(a)') line34
      ! Layer lines
      do
         read (io_locmod,'(a)',iostat=ios) line34
         if (ios .lt. 0) exit
         write (io_log,'(a)') line34
      end do
      close (io_locmod)
   end if

   ! Log file for phase reidentification
   if (plogout) then
      plog_file = trim(outfile)//'.plog'
      call open_file (io_plog, plog_file, 'new')
   end if

   ! Log file for station data
   stn_log_file = trim(outfile)//'.stn'
   call open_file (io_stn_log, stn_log_file, 'new')
   
   ! Log file for travel-time calculations
   if (log_tt) then
      tlog_file = trim(outfile)//'.tlog'
      call open_file (io_tt_log, tlog_file, 'new')
   end if

   ! Skipped stations
   if (n_skip .gt. 0) then
      do i = 1,n_skip
         if (skip_params(i,1) .ne. '*') write (io_stn_log,'(a,1x,a,1x,a,1x,a)')&
          'skipping station ', (skip_params(i,j),j=1,3)
      end do
   end if
   
   ! Supplemental station files
   if (n_supp_stn_file .gt. 0) then
      write (io_stn_log,'(a)') 'Supplemental station files: '
      do i = 1,n_supp_stn_file
         write (io_stn_log,'(a,i1,1x,a)') 'SSTN-', i, trim(suppfilnam(i))
      end do
   end if
   
   ! NEIC station metadata
   If (nsmd) then
      filename = trim(station_path)//dirsym//trim(nsmd_file)
      call read_nsmd (io_nsmd, filename, outfile)
   end if
         
   ! ISC-FDSN station metadata
   If (ifsm) then
      filename = trim(station_path)//dirsym//trim(ifsm_file)
      call read_ifsm (io_ifsm, filename, outfile)
   end if

   ! Phase-specific default reading errors
   filename = trim(psdre_path)//dirsym//trim(psdre_file)
   inquire (file=filename,exist=ex)
   if (ex) then
      io_temp = lunit()
      call open_file (io_temp, filename, 'old')
      i = 0
      do
         read (io_temp,'(a)',iostat=ios) line12
         if (ios .ge. 0) then
            i = i + 1
            if (i .le. n_psdre_max) then
               psdre_phase(i) = line12(1:8) 
               read (line12(10:12),'(f3.1)') psdre(i)
            else
               call warnings ('mloc_main: maximum number of phase-specific default reading errors reached')
               n_psdre = i - 1
               exit
            end if
         else ! EOF
            n_psdre = i
            exit
         end if
      end do
      close (io_temp)
      if (verbose_screen) then
         write (msg,'(a,i2,a)') 'mloc_main: ', n_psdre, ' phase-specific default reading errors read'
         call fyi (trim(msg))
      end if
   else
      call oops ('mloc_main: '//trim(filename)//' not found')
   end if
   
   ! Topography file used for bounce point corrections (pP, sP, pwP)
   filename = trim(bp_path)//dirsym//trim(bp_file)
   io_bptc = lunit()
   call open_file (io_bptc, filename, 'old')
   do j = 1,n_bp_lat
      read (io_bptc,'(4321i7)') (bp_topo(i,j),i=1,n_bp_lon)
   end do
   close(io_bptc)
   if (verbose_screen) call fyi ('mloc_main: bounce point topography data read from '//trim(filename))
   
   ! Bounce point correction log file
   if (bptc) then
      bptc_log_file = trim(outfile)//'.bptc'
      call open_file (io_bptc_log, bptc_log_file, 'new')
   end if
   
   ! Relative depth phase plots when the "all" argument is given
   if (rdpp_all) then
      n_rdpp = n_event
      do iev = 1,n_event
         rdpp_evt(iev) = evtnam(iev)
      end do
   end if
   
   ! tt5e plots when the "all" argument is given
   if (tt5e_all) then
      n_tt5e = n_event
      do iev = 1,n_event
         tt5e_evt(iev) = evtnam(iev)
      end do
   end if
   
   ! Residual plots when the "all" argument is given
   if (evrp_all) then
      n_evrp = n_event
      do iev = 1,n_event
         evrp_event(iev) = evtnam(iev)
         evrp_delta_min(iev) = evrp_delta_min_all
         evrp_delta_max(iev) = evrp_delta_max_all
      end do
   end if   
   
   ! Junk phase list
   filename = trim(phases_path)//dirsym//trim(junk_phases_file)
   call open_file (io_junk, filename, 'old')
   do
     read (io_junk,'(a)',iostat=ios) line8
     if (ios .lt. 0) exit
     n_junk = n_junk + 1
     if (n_junk .gt. n_junk_max) then
        call warnings ('mloc_main: too many junk phases read; increase n_junk_max')
        n_junk = n_junk -1
        exit
     end if
     junk_phase(n_junk) = trim(line8)
   end do
   close (io_junk)
   if (verbose_screen) then
      write (msg,'(a,i3,a)') 'mloc_main: ', n_junk, ' junk phases read from '//trim(filename)
      call fyi (trim(msg))
   end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
   call mlocset ! Set up problem, do linearized inversion, write output files and plots
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
   
   ! Close files
   close (io_log)
   close (io_alloc_log)
   close (io_stn_log)
   close (io_taup)
   if (plogout) close (io_plog)
   if (log_tt) close (io_tt_log) 
   ! Catch any files still open
   do luopen = 10,99
      inquire (luopen,opened=op)
      if (op) then
         inquire (luopen,name=open_file_name)
         write (msg,'(a,i2,1x,a)') 'mloc_main: file still open on logical unit ', luopen, trim(open_file_name)
         call fyi (trim(msg))
         close (luopen)
      end if
   end do
   
   call cpu_time (cpu_finish)
   write (msg,'(a,f12.3,a)') 'mloc_main: Total CPU usage =  ', cpu_finish-cpu_start, ' seconds'
   call fyi (trim(msg))
      
   write (*,'(/a/)') '*** Run completed ***'

   stop
   
end program mloc

