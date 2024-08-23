!> Miscellaneous procedures, plus several procedures that would normally be in a different
! module because of their functionality (especially mloc_phases) but need to be compiled
! earlier to avoid circularity between modules.

module mloc_lib

   use mloc_declare
   use mloc_math
   
   implicit none
   save
   
contains

!*****************************************************************************************
subroutine open_file (io_unit, file_name, file_type)

! Opens a file (old or new) with the specified name on the requested logical unit.
! Files are assumed to be sequential access and formatted.

   integer, intent(in) :: io_unit
   character(len=*), intent(in) :: file_name
   character(len=*), intent(in) :: file_type ! old or new
   character(len=132) :: msg
   logical :: op
   logical :: ex
   
   if (file_type .eq. 'new') then
      inquire (file=file_name,exist=ex)
      if (ex) then
         write (msg,'(3a)') 'open_file: ', file_name, ' already exists'
         call oops(trim(msg))
      end if
   else if (file_type .ne. 'old') then
      write (msg,'(3a)') 'open_file: unrecognized file_type "', trim(file_type), '"'
      call oops(trim(msg))
   end if
   inquire (unit=io_unit,opened=op)
   if (op) then
      write (msg,'(a,i2,a)') 'open_file: ', io_unit, ' is already attached to a file'
      call oops(trim(msg))
   end if
   open (io_unit,file=file_name,status=file_type)
   inquire (file=file_name,opened=op)
   if (op) then
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'open_file: ', trim(file_name), ' opened on unit ', io_unit
         call fyi (trim(msg))
      end if
   else
      write (msg,'(3a,i3)') 'open_file: ', trim(file_name), ' failed to open on unit ', io_unit
      call oops (trim(msg))
   end if
   
   return
   
end subroutine open_file


!*****************************************************************************************
   subroutine calibration_code (iev, calibration_type, cal_code)
   
   ! Constructs the GTCU code for calibration level
   !  calibration_type = 0 for uncalibrated
   !                     1 for direct calibration
   !                     2 for indirect calibration
   
      integer, intent(in) :: calibration_type
      integer, intent(in) :: iev
      integer :: igt
      character(len=4), intent(out) :: cal_code
      character(len=2) :: scale_length_pr
      logical :: epicenter_free
      logical :: depth_calibrated
      
      ! Calibrated events must have both latitude and longitude as free parameters
      if (mindx(iev,1) .gt. 0 .and. mindx(iev,2) .gt. 0) then
         epicenter_free = .true.
      else
         epicenter_free = .false.
      end if
      
      ! Depth calibration.
      ! Uncalibrated events can still have calibrated depth (e.g., from depth phases).
      ! If depth is a free parameter it is assumed to be calibrated, but certain values of
      ! depset_pr are also taken to mean a calibrated depth. Depth is also considered to be
      ! calibrated if an event is a calibration event with cal_par = 'h'.
      depth_calibrated = .false.
      if (cal_event(iev,3) .and. (cal_par(iev) .eq. 'f' .or. cal_par(iev) .eq. 'h')) depth_calibrated = .true.
      if (mindx(iev,3) .gt. 0) then ! Free depth
         depth_calibrated = .true.
      else ! Depth fixed but calibrated in some way
         if (depset_pr(iev) .eq. 'd') then ! Depth phases
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'e') then ! engineering information
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'f') then ! from fault modeling (InSAR, GPS, etc.)
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'l') then ! from local-distance readings
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'm') then ! from an mloc run with free depth
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'n') then ! from near-source readings
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'r') then ! from relocation outside of mloc with free depth
            depth_calibrated = .true.
         else if (depset_pr(iev) .eq. 'w') then ! from waveform analysis
            depth_calibrated = .true.
         end if
      end if
      
      ! Scale length ! Nearest integer to the semi-major axis length, capped at 99.
      ! Always written as a 2-digit integer.
      select case (calibration_type)
      
         case (0) ! Uncalibrated
            igt = nint(xl2c(iev))
            
         case (1) ! Direct Calibration
            igt = nint(xl2dc(iev))
            
         case (2) ! Indirect Calibration
            igt = nint(xl2cg(iev))
            
      end select
      
      if (igt .le. 99) then
         write (scale_length_pr,'(i2.2)') igt
      else
         scale_length_pr = '99'
      end if
      
      select case (calibration_type)
      
         case (0) ! Uncalibrated
            if (depth_calibrated) then
               cal_code = 'UF'//scale_length_pr
            else
               cal_code = 'UE'//scale_length_pr
            end if
         
         case (1) ! Direct calibration
            if (epicenter_free) then
               if (depth_calibrated) then
                  cal_code = 'CH'//scale_length_pr
               else
                  cal_code = 'CT'//scale_length_pr
               end if
            else
               if (depth_calibrated) then
                  cal_code = 'UF'//scale_length_pr
               else
                  cal_code = 'UE'//scale_length_pr
               end if
            end if
         
         case (2) ! Indirect calibration
            if (epicenter_free) then
               if (depth_calibrated) then
                  if (ot_cal) then
                     cal_code = 'CH'//scale_length_pr
                  else
                     cal_code = 'CF'//scale_length_pr
                  end if
               else
                  if (ot_cal .or. cal_par(iev) .eq. 't') then
                     cal_code = 'CT'//scale_length_pr
                  else
                     cal_code = 'CE'//scale_length_pr
                  end if
               end if
            else
               if (depth_calibrated) then
                  cal_code = 'UF'//scale_length_pr
               else
                  cal_code = 'UE'//scale_length_pr
               end if
            end if
         
      end select
         
      return
      
   end subroutine calibration_code
   
   
!*****************************************************************************************
   subroutine readerr (iev, ird)
   
   ! Sets default reading errors for phase arrivals. There is a hierarchy of methods:
   ! 1) A set of phase-specific default reading errors is automatically read during start-up.
   ! 2) If a reading is not found in that phase list, the old hard-wired algorithm is used.
   ! 3) Whatever value results from these first two methods is over-ridden by the value of
   !    empirical reading error from a previous run (.rderr file), if it has been specified.
   !    A minimum allowed value is still enforced.
   ! 4) Reading errors for local-distance stations may be further reset by the 'rels' command,
   !    with no minimum value enforced.
   ! 5) Finally, a check is made to ensure that the reading error is not less than the standard
   !    deviation of the uniform rectangular continuous distribution over the range defined
   !    by the precision of the reading.
         
      integer :: i
      integer, intent(in) :: iev
      integer, intent(in) :: ird
      real :: readerror
      character(len=132) :: msg
      character(len=21) :: qtest
!       logical :: ptype
!       logical :: stype
      
      ! Reading errors for differential time data are not set here.
      if (idiff(iev,ird) .gt. 0) return
   
      if (.not.data_weight) then ! No weighting by reading precision.
         sdread(iev,ird) = 1.0
         return
      end if
      
      ! Initialize
      sdread(iev,ird) = 0.
      
      ! See if a value can be found in the phase-specific default reading error list.
      do i = 1,n_psdre
         if (phase(iev,ird) .eq. psdre_phase(i)) then
            sdread(iev,ird) = psdre(i)
         else ! Use the old hard-wired values
      
            ! First cut
            if (phase(iev,ird) .eq. 'P       ' .or. phase(iev,ird)(1:3) .eq. 'PKP') then
               sdread(iev,ird) = 0.5
            else
               sdread(iev,ird) = 0.8
            end if
      
            ! Depth phases
            if (phase(iev,ird)(1:1) .eq. 'p') then
               sdread(iev,ird) = 1.0
            else if (phase(iev,ird)(1:1) .eq. 's') then
               sdread(iev,ird) = 1.5
            end if
      
            ! S-P
            if (phase(iev,ird)(1:3) .eq. 'S-P') sdread(iev,ird) = 0.20
        
            ! Increase reading error for phases with S-type.
            if (stype(phase(iev,ird))) sdread(iev,ird) = sdread(iev,ird)*2.0
            
            ! P_
            if (phase(iev,ird)(1:2) .eq. 'P_') sdread(iev,ird) = 4.0
      
            ! Lg
            if (phase(iev,ird)(1:2) .eq. 'Lg') sdread(iev,ird) = 4.0
      
            ! Rg
            if (phase(iev,ird)(1:2) .eq. 'Rg') sdread(iev,ird) = 4.0
    
             ! T-phase
            if (phase(iev,ird)(1:2) .eq. 'T ') sdread(iev,ird) = 4.0
         
         end if
      end do
      
      ! Set reading errors according to observed distribution of
      ! residuals from a previous run. This is done by searching in a summary file ("rderr" file)
      ! of station/phase residuals. A minimum allowed reading error is still enforced.
      if (read_rderr) then
         qtest = stname(iev,ird)//deployment(iev,ird)//phase(iev,ird)
         readerror = -1.
         do i = 1,nqname
            if (qtest .eq. qname(i)) then
               readerror = rdsigma(i)
               if (phase(iev,ird)(2:2) .eq. 'g' .or. phase(iev,ird)(2:2) .eq. 'b') then
                  sdread(iev,ird) = amax1(readerror,rderr_min_loc)
               else if (phase(iev,ird)(1:1) .eq. 'p' .or. phase(iev,ird)(1:1) .eq. 's') then
                  sdread(iev,ird) = amax1(readerror,rderr_min_depth)
               else
                  sdread(iev,ird) = amax1(readerror,rderr_min)
               end if
               exit
            end if
         end do
      end if
   
      ! Set reading errors for crustal phases to fixed values.
      if (rels_set) then
         if ((phase(iev,ird)(2:2) .eq. 'g' .or. phase(iev,ird)(2:2) .eq. 'b') .and. delt(iev,ird) .le. rderr_loc_delt) then
            if (ptype(phase(iev,ird))) sdread(iev,ird) = rderr_loc_p
            if (stype(phase(iev,ird))) sdread(iev,ird) = rderr_loc_s
         end if
      end if
      
      ! Reading errors cannot be smaller than the standard deviation of the uniform rectangular
      ! continuous distribution over the range defined by the precision of the reading. This is
      ! irrelevant for iptim=-1 or iptim=-2, since the theoretical standard deviation is smaller
      ! than any plausible minimum reading error.
      if (iptim(iev,ird) .eq. 0) then ! Precision to nearest second.
         sdread(iev,ird) = max(sdread(iev,ird),0.29)
      else if (iptim(iev,ird) .eq. 1) then ! Precision to nearest 10 seconds.
         sdread(iev,ird) = max(sdread(iev,ird),2.9)
      else if (iptim(iev,ird) .eq. 2) then ! Precision to nearest minute.
         sdread(iev,ird) = max(sdread(iev,ird),17.0)
      else if (iptim(iev,ird) .eq. 3) then ! Precision to nearest tenth of a minute.
         sdread(iev,ird) = max(sdread(iev,ird),1.7)
      end if
      
      ! Catch anything that fell through the cracks with zero reading error
      if (sdread(iev,ird) .lt. 0.01) then
         write (msg,'(a,i3,a,1x,a,a,e12.4,a)') 'readerr: ', iev, stname(iev,ird), phase(iev,ird),&
          ': reading error of ', sdread(iev,ird), ' not allowed, set to 1.0'
         call warnings (trim(msg))
         sdread(iev,ird) = 1.0
      end if
            
      return
      
   end subroutine readerr
         
   
!*****************************************************************************************
   subroutine getrderr ()
   
   ! Get empirical reading errors from a previous run (.rderr file), or use default values.
   ! As of v8.3, the first line of a .rderr file contains the hypocentroid used for plotting empirical path anomalies.
         
      integer :: i
      integer :: ios
      character(len=12) :: line
      character(len=132) :: msg
      logical :: op
      
      if (read_rderr) then
         
         call open_file (io_rderr,  rderrfname,  'old')
         
         inquire (unit=io_rderr,opened=op,name=rderrfname)
         if (op) then
            ! First line is the hypocentroid used for plotting empirical path anomalies. Skip it.
            read (io_rderr,'(a)') line
            if (line(1:12) .ne. 'Hypocentroid') then
               call warnings ('getrderr: attempt to read an invalid .rderr file! Using default values.')
               close (io_rderr)
               read_rderr = .false.
               rderrfname = 'default'
               return
            end if
            do i = 1,n_q_max
               read (io_rderr,'(a,t37,f10.3)',iostat=ios) qname(i), rdsigma(i)
               if (ios .lt. 0) exit
            end do
            close (io_rderr)
            nqname = i - 1
            if (verbose_screen) then
               write (msg,'(a,i4,2a)') 'getrderr: ', nqname, ' station-phase reading errors read from:',&
                trim(rderrfname)
               call fyi (trim(msg))
            end if
         else
            msg = 'getrderr: failed to open file '//trim(rderrfname)//'...using default values'
            call warnings (trim(msg))
            read_rderr = .false.
            rderrfname = 'default'
         end if
      else
         read_rderr = .false.
         rderrfname = 'default'
      end if
      
      return
      
   end subroutine getrderr
   
   
!*****************************************************************************************
   integer function lenb (string)
   
   !  finds index of last nonblank character in string
   
      integer :: i
      integer :: n
      character(len=*), intent(in) :: string
   
      n = len(string)
      do i = n,1,-1
         if (string(i:i) .ne. ' ') then
            lenb = i
            return
         end if
      end do
      lenb = 0
   
      return
      
   end function lenb
   
   
!*****************************************************************************************
   integer function lennb (string)
   
   !  finds index of first nonblank character in string
   
      integer :: i
      character(len=*) :: string
   
      do i = 1,len(string)
         if (string(i:i) .ne. ' ') then
            lennb = i
            return
         end if
      end do
      lennb = 1
   
      return
      
   end function lennb
   
   
!*****************************************************************************************
   subroutine dataflags (dflag, fcode)
   
   ! Set data flags
         
      character(len=1), parameter :: blank = ' '
      character(len=1) :: fcode
      logical :: dflag
      
      ! My flags
      if (.not. dflag) then
         if (fcode .eq. 'x') fcode = blank ! Don't use because of large residual
         if (fcode .eq. 'p') fcode = blank ! Don't use because of the phase
         if (fcode .eq. 'd') fcode = blank ! Don't use because of duplicate reading
         if (fcode .eq. 's') fcode = blank ! Don't use because this station is to skipped
         if (fcode .eq. 'a') fcode = blank ! Don't use because readings by this author are to be skipped
         if (fcode .eq. 'm') fcode = blank ! Don't use because this station is missing in the station list
         if (fcode .eq. 't') fcode = blank ! Don't use because timing is suspect
      end if
   
      return
      
   end subroutine dataflags
   
   
!*****************************************************************************************
   logical function bdp (iev, ird)
   
   ! Returns "true" if the station is on a list of stations that are suspected of reporting
   ! depth phases based on theoretical arrival times from a preliminary location such as the PDE.
   ! The test also considers the date of the event. It could consider the deployment as well,
   ! but for now it does not.
   
      integer :: i
      integer, intent(in) :: iev ! event number
      integer, intent(in) :: ird ! reading index
      integer :: jdate_begin
      integer :: jdate_end
      integer :: juldat_test
      integer :: julday
!      logical :: date_range
      
      bdp = .false.
      if (.not.bdp_list) return
      
      ! Julian date for current event
      call juldat (iyre(iev), mone(iev), idye(iev), julday, 0)
      juldat_test = iyre(iev)*1000 + julday
   
      do i = 1,n_bdp
         if (stname(iev,ird) .eq. bdp_station(i)(1:5)) then
            read (bdp_station(i)(53:59),'(i7)') jdate_begin
            read (bdp_station(i)(61:67),'(i7)') jdate_end
            if (date_range(juldat_test,jdate_begin,jdate_end)) bdp = .true.
            exit
         end if
      end do
      
      return
   
   end function bdp
   
   
!*****************************************************************************************
   logical function ptype (phasein)
   
   ! returns .true. if phasein is a P phase
   ! The first leg of depth phases is ignored.
   
!      integer :: lenb
      integer :: k
      character(len=2) :: c2
      character(len=3) :: c3
      character(len=4) :: c4
      character(len=5) :: c5
      character(len=8), intent(in) :: phasein
      logical :: dphase
   
      if (phasein .eq. '        ') then
         ptype = .false.
         return
      end if
         
      if (phasein .eq. 'UNKNOWNP') then ! Unknown P-type phase
         ptype = .true.
         return
      end if
   
      if (phasein(1:2) .eq. 'T ') then ! T is considered an P phase.
         ptype = .true.
         return
      end if
   
      k = lenb(phasein)
      c2 = '  '
      c3 = '   '
      c4 = '    '
      c5 = '     '
   
      dphase = .false.
      if (phasein(1:1) .eq. 'p' .or. phasein(1:1) .eq. 's') dphase = .true. ! Depth phase
      if (dphase) then
         if (k .ge. 3) c2 = phasein(2:3)
         if (k .ge. 4) c3 = phasein(2:4)
         if (k .ge. 5) c4 = phasein(2:5)
         if (k .ge. 6) c5 = phasein(2:6)
      else
         if (k .ge. 2) c2 = phasein(1:2)
         if (k .ge. 3) c3 = phasein(1:3)
         if (k .ge. 4) c4 = phasein(1:4)
         if (k .ge. 5) c5 = phasein(1:5)
      end if
   
      if (phasein(k:k) .eq. 'P') then ! Check last leg
         ptype = .true.
      else if (c2 .eq. 'Pn' .or.&
               c2 .eq. 'Pb' .or.&
               c2 .eq. 'P*' .or.&
               c2 .eq. 'Pg') then
         ptype = .true.
      else if (c3 .eq. 'PKP' .or.&
               c3 .eq. 'SKP') then
         ptype = .true.
      else if (c4 .eq. "P'P'" .or.&
               c4 .eq. 'PKKP' .or.&
               c4 .eq. 'SKKP') then
         ptype = .true.
      else if (c5 .eq. 'Pdiff') then
         ptype = .true.
      else
         ptype = .false.
      end if
   
      return
      
   end function ptype
   
               
!*****************************************************************************************
   logical function stype (phasein)
   
   ! Returns .true. if phasein is an S phase.
   ! The first leg of depth phases is ignored.
   
!      integer :: lenb
      integer :: k
      character(len=2) :: c2
      character(len=3) :: c3
      character(len=4) :: c4
      character(len=5) :: c5
      character(len=8), intent(in) :: phasein
      logical :: dphase
   
      if (phasein .eq. '        ') then
         stype = .false.
         return
      end if
         
      if (phasein .eq. 'UNKNOWNS') then ! Unknown S-type phase
         stype = .true.
         return
      end if
   
      if (phasein(1:2) .eq. 'P_') then ! P_ is considered a P phase.
         stype = .false.
         return
      end if
   
      if (phasein(1:2) .eq. 'Lg') then ! Lg is considered an S phase.
         stype = .true.
         return
      end if
   
      if (phasein(1:2) .eq. 'Rg') then ! Rg is considered a P phase.
         stype = .false.
         return
      end if
   
      if (phasein(1:2) .eq. 'T ') then ! T is considered an P phase.
         stype = .false.
         return
      end if
   
      k = lenb(phasein)
      c2 = '  '
      c3 = '   '
      c4 = '    '
      c5 = '     '
   
      dphase = .false.
      if (phasein(1:1) .eq. 'p' .or.&
          phasein(1:1) .eq. 's') dphase = .true. ! Depth phase
      if (dphase) then
         if (k .ge. 3) c2 = phasein(2:3)
         if (k .ge. 4) c3 = phasein(2:4)
         if (k .ge. 5) c4 = phasein(2:5)
         if (k .ge. 6) c5 = phasein(2:6)
      else
         if (k .ge. 2) c2 = phasein(1:2)
         if (k .ge. 3) c3 = phasein(1:3)
         if (k .ge. 4) c4 = phasein(1:4)
         if (k .ge. 5) c5 = phasein(1:5)
      end if
   
      if (phasein(k:k) .eq. 'S') then ! Check last leg
         stype = .true.
      else if (c2 .eq. 'Sn' .or.&
               c2 .eq. 'Sb' .or.&
               c2 .eq. 'S*' .or.&
               c2 .eq. 'Sg') then
         stype = .true.
      else if (c3 .eq. 'SKS' .or.&
               c3 .eq. 'PKS') then
         stype = .true.
      else if (c4 .eq. "S'S'" .or.&
               c4 .eq. 'SKKS' .or.&
               c4 .eq. 'PKKS') then
         stype = .true.
      else if (c5 .eq. 'Sdiff') then
         stype = .true.
      else
         stype = .false.
      end if
   
      return
      
   end function stype
   
         
!*****************************************************************************************
   logical function crustal_phase (phasein)
         
   ! Returns "true" if the candidate phase matches an entry in the list of crustal phases.
   ! Basic list of crustal phases from Storchak et al., The IASPEI Standard Seismic Phase List,
   ! Seismological Research Letters, v. 74, No. 6, 2003. Crustal surface waves (Lg, Rg) 
   ! are not considered, nor are multiple Moho reflections, e.g., PmPN. Some additional
   ! phases have been added (PbPb, SbSb, depth phases).
         
      integer, parameter :: n_crustal_phases = 28
      
      integer :: i
      character(len=8), intent(in) :: phasein
      character(len=8), dimension(n_crustal_phases) :: crustal_phases
      
      data crustal_phases /'Pg      ','Pb      ','Pn      ','Sg      ','Sb      ','Sn      ',&
						   'pPg     ','pPb     ','pPn     ','pSg     ','pSb     ','pSn     ',&
						   'sPg     ','sPb     ','sPn     ','sSg     ','sSb     ','sSn     ',&
						   'PgPg    ','PbPb    ','PnPn    ','SgSg    ','SbSb    ','SnSn    ',&
						   'PmP     ','PmS     ','SmS     ','SmP     '/
      
      crustal_phase = .false.
      do i = 1,n_crustal_phases
         if (phasein .eq. crustal_phases(i)) then
            crustal_phase = .true.
            exit
         end if
      end do
      
      return
      
   end function crustal_phase
   
   
!*****************************************************************************************
   logical function skipp (phasein)
   
   ! Returns a logical value controlling whether a phase should be kept
   ! in the dataset or not. Mostly used to get rid of various forms of
   ! surface waves, and amplitude measurements.
   
   ! Some exotic phases are skipped because they are not in the ak135 phase set.
   
      character(len=1) :: c
      character(len=8), intent(in) :: phasein
   
      skipp = .false.
      if (phasein(1:2) .eq. 'P_') return
      if (phasein(1:2) .eq. 'Lg') return
      if (phasein(1:2) .eq. 'Rg') return
      if (phasein(1:2) .eq. 'T ') return
   
      if (phasein(1:3) .eq. 'amp') then
         skipp = .true.
         return
      end if
   
      if (phasein(1:3) .eq. 'PPP' .or. phasein(1:3) .eq. 'SSS') then ! Not in ak135 phase set
         skipp = .true.
         return
      end if
   
      if (phasein(1:4) .eq. 'P3KP') then ! Not in ak135 phase set
         skipp = .true.
         return
      end if
   
      c = phasein(1:1)
      if (c .eq. 'L') then ! Catches, e.g., L, LR. LQ
         skipp = .true.
      else if (c .eq. 'M') then ! Catches, e.g.,  M, MLR, MAXIMUM
         skipp = .true.
      else if (c .eq. 'R') then ! Catches, e.g., R, RM
         skipp = .true.
      else if (c .eq. 'Q') then ! Catches, e.g., Q, QM
         skipp = .true.
      else if (c .eq. 'A') then ! Catches AMP
         skipp = .true.
      else if (c .eq. 'F') then ! Catches FINAL
         skipp = .true.
      end if
   
      if (trim(phasein) .eq. 'Pmax') skipp = .true.
      if (trim(phasein) .eq. 'Smax') skipp = .true.
      if (trim(phasein) .eq. 'max') skipp = .true.
   
      return
      
   end function skipp
   

!*****************************************************************************************
logical function depth_constraint (c)

! Based on the depth code, returns .true. if the associated focal depth is considered constrained.
! Returns .false. if not.

   character(len=1) :: c
      
   if (c .eq. 'c') then ! cluster default depth 
      depth_constraint = .false.
   else if (c .eq. 'd') then ! depth phases
      depth_constraint = .true.
   else if (c .eq. 'e') then ! engineered (man-made explosion)
      depth_constraint = .true.
   else if (c .eq. 'f') then ! fault model (InSAR, GPS, etc.)
      depth_constraint = .true.
   else if (c .eq. 'i') then ! input data file
      depth_constraint = .false.
   else if (c .eq. 'l') then ! local distance readings (more than 2-3 focal depths)
      depth_constraint = .true.
   else if (c .eq. 'm') then ! mloc solution (with free depth)
      depth_constraint = .true.
   else if (c .eq. 'n') then ! near-source station readings
      depth_constraint = .true.
   else if (c .eq. 'r') then ! relocation (outside mloc) with free depth
      depth_constraint = .true.
   else if (c .eq. 'u') then ! unknown
      depth_constraint = .false.
   else if (c .eq. 'w') then ! waveform analysis
      depth_constraint = .true.
   else if (c .eq. ' ') then ! blank
      depth_constraint = .false.
   else
      call warnings ('depth_constraint: unknown depth code "'//trim(c)//'"')
      depth_constraint = .false.
   end if
                      
   return
   
end function depth_constraint

   
!*****************************************************************************************
   logical function pass (fcode)
      
      character(len=1) :: fcode
      
      pass = (fcode .ne. 'x' .and. &
              fcode .ne. 'd' .and. &
              fcode .ne. 'p' .and. &
              fcode .ne. 's' .and. &
              fcode .ne. 'a' .and. &
              fcode .ne. 'm' .and. &
              fcode .ne. 't')
      
      return
      
   end function pass


!*****************************************************************************************
end module mloc_lib

