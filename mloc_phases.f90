!> Procedures related to seismic phases

module mloc_phases

   use mloc_declare
   use mloc_taup
   use mloc_lib
   use mloc_math
   use mloc_tt

   implicit none
   save
   
   ! No-duplicates phase list
   integer :: nphasend ! Number of entries in the no-duplicates phase list
   real, dimension(max_mpl) :: ttnd ! travel times in the no-duplicate phase 
   character(len=8), dimension(max_mpl) :: phcdnd ! phase codes in the no-duplicates list
   
contains

!*****************************************************************************************
   subroutine phreid3 (it, iev, ird, nr)
         
   ! Phase identification for single readings. No concept of groups of readings.
   ! This version uses a traditional "best fit" approach, but modified to take advantage of
   ! information on the probability distribution functions (PDFs) of different phases.
   ! The target arrival time is tested against all phases in the theoretical TT model for the
   ! corresponding focal depth and epicentral distance, regardless arrival time order. Probability
   ! is calculated for each possible phase ID, based on the candidate phases's PDF
   ! and the choice is made on the basis of highest probability. There are a couple of rules:
   !  1) Depth phases are handled separately
   !  2) Phase type (P or S) is honored if a phase name has been provided
   ! A phase name is not changed unless the probability of the new phase ID is 0.05 or greater.
   ! The necessary data of coefficients to describe the PDFs of all the needed
   ! phases is not yet available. Therefore, at this time all PDS are the same and the choice boils
   ! down to the classical "smallest residual" criterion.
   
      ! Travel time variables
      integer :: nphase
      real, dimension(max_mpl) :: dtdd
      real, dimension(max_mpl) :: dtdh
      real, dimension(max_mpl) :: dddp
      real, dimension(max_mpl) :: tt ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd

      integer :: i
      integer, intent(in) :: iev ! event number
      integer, intent(in) :: ird ! reading index
      integer, intent(in) :: it ! iteration number
      integer, intent(out) :: nr ! =1 if re-identified, 0 otherwise
      real :: arrtime
!      real :: hms2s
      real :: otspievit
!      real :: prob
      real :: probi
      real :: probmax
      real :: ttobs
      real, dimension(2) :: usrc
      character(len=8) :: phasein
      character(len=8) :: phaseout
!       logical :: ptype
!       logical :: stype
      
      nr = 0
            
      ! Theoretical travel-time
      call depset (depthp(iev,it), usrc)
      if (.not.locmod) then
         call trtm (delt(iev,ird), nphase, tt, dtdd, dtdh, dddp, phcd)
      else
         call tt_mixed_model (depthp(iev,it), delt(iev,ird), nphase, tt, dtdd, dtdh, dddp, phcd)
      end if
      
!       call nodupe (nphase, tt, phcd, nphasend, ttnd, phcdnd) ! Get the equivalent phase list without duplicates
!       if (nphasend .eq. 0) then
!          if (plogout) write (io_plog,'(a,f10.3)') 'phreid3: no phases returned delta = ', delt(iev,ird)
!          nphasend = 1
!          phcdnd(1) = 'CRAP    '
!       end if
      
      ! Some cases in which phase identification is not done. 
      if (phase(iev,ird) .eq. phcd(1)) return
      if (phase(iev,ird)(1:2) .eq. 'T ') return
      if (phase(iev,ird)(1:3) .eq. 'pwP') return
      
      phasein = phase(iev,ird)
      if (phase(iev,ird) .eq. 'P_      ' .and. delt(iev,ird) .le. p__min) phasein = 'Pg      '
      if (phase(iev,ird) .eq. 'Lg      ' .and. delt(iev,ird) .le. lg_min) phasein = 'Sg      '
      if (phase(iev,ird) .eq. 'Rg      ' .and. delt(iev,ird) .le. rg_min) phasein = 'Pg      '
      
      ! Observed travel time in seconds, relative to current origin time
      otspievit = hms2s(hourp(iev,it),minp(iev,it),secp(iev,it))
      arrtime = hms2s(ipah(iev,ird),ipam(iev,ird),pas(iev,ird))
      if (ipah(iev,ird) .lt. hourp(iev,it)) arrtime = arrtime + 86400.
      ttobs = arrtime - otspievit
      
      if (nphase .eq. 0) then ! Phase ID not changed.
         phaseout = phasein
      else if (nphase .ge. 1) then
         ! Probability of candidate phases
         ! Depth phases can only be identified as depth phases.
         probmax = 0.
         if (phasein(1:1) .eq. 'p' .or. phasein(1:1) .eq. 's') then
            do i = 1,nphase
               if (phcd(i)(1:1) .ne. 'p' .and. phcd(i)(1:1) .ne. 's') cycle
               probi = prob(phcd(i), 0, ttobs, delt(iev,ird), depthp(iev,it), phcd, tt, nphase)
               if (probi .gt. probmax) then
                  probmax = probi
                  phaseout = phcd(i)
               end if
            end do
         else
            do i = 1,nphase
               if (phcd(i)(1:1) .eq. 'p' .or. phcd(i)(1:1) .eq. 's') cycle
               if (ptype(phasein) .and. .not.ptype(phcd(i))) cycle
               if (stype(phasein) .and. .not.stype(phcd(i))) cycle
               probi = prob(phcd(i), 0, ttobs, delt(iev,ird), depthp(iev,it), phcd, tt, nphase)
               if (probi .gt. probmax) then
                  probmax = probi
                  phaseout = phcd(i)
               end if
            end do
            
            ! P_ is not in the standard phase list so check it separately.
            ! Also check if a Pg or Pn might be re-identified as P_ 
            if (phasein .eq. 'P_      ' .or. phasein .eq. 'Pg      ' .or. phasein .eq. 'Pn      ') then
               if (delt(iev,ird) .gt. p__min) then
                  probi = prob('P_      ', 0, ttobs, delt(iev,ird), depthp(iev,it), phcd, tt, nphase)
                  if (probi .gt. probmax) then
                     probmax = probi
                     phaseout = 'P_      '
                  end if
               end if
            end if
            
            ! Lg is not in the standard phase list so check it separately.
            ! Also check if an Sg or Sn might be re-identified as Lg
            if (phasein .eq. 'Lg      ' .or. phasein .eq. 'Sg      ' .or. phasein .eq. 'Sn      ') then
               if (delt(iev,ird) .gt. lg_min) then
                  probi = prob('Lg      ', 0, ttobs, delt(iev,ird), depthp(iev,it), phcd, tt, nphase)
                  if (probi .gt. probmax) then
                     probmax = probi
                     phaseout = 'Lg      '
                  end if
               end if
            end if
            
            ! Rg is not in the standard phase list so check it separately.
            ! Also check if a Pg might be re-identified as Rg
            if (phasein .eq. 'Rg      ' .or. phasein .eq. 'Pg      ') then
               if (delt(iev,ird) .gt. rg_min) then
                  probi = prob('Rg      ', 0, ttobs, delt(iev,ird), depthp(iev,it), phcd, tt, nphase)
                  if (probi .gt. probmax) then
                     probmax = probi
                     phaseout = 'Rg      '
                  end if
               end if
            end if
            
         end if
         if (probmax .lt. 0.05) phaseout = phasein ! Leave it unchanged
         if (plogout .and. verbose_log) write (io_plog,'(2a,2f8.2,2f8.3,1x,a)') 'phreid3: ',&
          phasein, delt(iev,ird), ttobs, probi, probmax, phaseout
      end if
         
      phase(iev,ird) = phaseout
      
      if (phaseout .ne. phasein) then
         call readerr (iev,ird) ! Assign the correct reading error for the new phase ID
         nr = nr + 1
         if (plogout) write (io_plog,'(a,i3,1x,a,2(1x,a),f7.2,1x,a)') 'phreid3: ', iev, stname(iev,ird),&
          phasein, phaseout, delt(iev,ird), infile20(iev)
      end if
                               
      return
      
   end subroutine phreid3
   
   
!*****************************************************************************************
   real function prob (phasein, itype, t, d, h, phcd, tt, nphase)
         
   ! For a given phase arrival, return the probability of being observed.
   ! In general, the calculation is done on the basis of:
   !   1) relative observability of that phase
   !   2) median (absolute or relative to AK135) as a function of delta
   !   3) spread as a function of delta
   !   4) a composit PDF which defines the shape of the distribution 
   ! The details of the algorithm can be different for each phase
   
   ! Input: itype - 0 for travel time, 1 for residual
   !        t - time, either travel time or a residual
   !        d - epicentral distance
   !        h - focal depth
   
      ! Travel time variables
      integer :: nphase
      real, dimension(max_mpl) :: tt ! travel time, s
      character(len=8), dimension(max_mpl) :: phcd

      integer :: i
      integer, intent(in) :: itype
      integer, parameter :: mb = 8
      real, parameter, dimension(mb) :: b = (/-0.07,0.88,-0.40,0.58,0.13,2.00,1.10,0.10/)
      real, intent(in) :: d
      real :: f
      real, intent(in) :: h
      real :: psprd
      real :: r
      real :: robs
      real, intent(in) :: t
      real :: ttlg
      real :: ttp_
      real :: ttrg
      real :: tttp
      real :: x
      real :: xmed
      character(len=1) :: c
      character(len=8), intent(in) :: phasein
!       logical :: ptype
!       logical :: stype
            
      if (itype .eq. 0) then ! Convert travel time to residual
         do i = 1,nphase
            if (phcd(i) .eq. phasein) then
               r = t - tt(i)
               go to 10 ! to avoid later arrivals of the same phase
            end if
         end do
         if (phasein(1:2) .eq. 'P_') then
            ttp_ = p__a + d*p__b
            r = t - ttp_
         end if
         if (phasein(1:2) .eq. 'Lg') then
            ttlg = lg_a + d*lg_b
            r = t - ttlg
         end if
         if (phasein(1:2) .eq. 'Rg') then
            ttrg = rg_a + d*rg_b
            r = t - ttrg
         end if
         if (phasein(1:2) .eq. 'T ') then
            tttp = tphase_a + d*tphase_b
            r = t - tttp
         end if
   10    continue
      else
         r = t
      end if
      
      ! Default values
      xmed = 0.
      if (ptype(phasein)) then
         psprd = 2.5
      else if (stype(phasein)) then
         psprd = 4.
      else
         psprd = 3.
      end if
      c = phasein(1:1)
      robs = 1.
      if (c .eq. 'p' .or. c .eq. 's') then
         if (h .le. 15) then
            robs = robs*0.4
         else
            robs = robs*0.6
         end if
      end if
      
      ! Use offset from a .ttsprd file if it exists. This helps remove bias due to baseline offset
      ! for different phases.
      ! It is problematic to use the spread until I have distance-dependent parameterization. Use
      ! of a single value for a phase such as Pn, which covers a large epicentral distance range
      ! and has greatly increasing spread with distance, creates problems at both ends. Even
      ! the offset creates problems in certain situations (direct calibration) if it is not
      ! distance-dependent
      
   !      if (read_ttsprd) then
   !         do i = 1,nsprd
   !            if (phasein .eq. sprdph(i)) then
   !               xmed = ttoffset(i)
   !               psprd = sprd(i)
   !               exit
   !            end if
   !         end do
   !      end if
         
      x = (r - xmed)/psprd
      call gswt (x, mb, b, f)
      prob = f*robs
      
      if (plogout) then
         if (prob .gt. 1.e-4) write (io_plog,'(2a,1x,f4.1,2f10.2,3f8.3)') 'prob: ', phasein, psprd, r, x, f, robs, prob
      end if
      
      return
      
   end function prob
   
         
!*****************************************************************************************
   subroutine gswt (x, m, b, f0)
   
   ! Derived from GSWT2, by removing calculation for derivatives, and entrypoint.
   ! gswt evaluates a function at a given location. The function is composed
   ! of one or more PDFs (a Gaussian plus one or two Cauchy PDFs).
   
   ! Input:      
   ! x  = the location at which the function is evaluated.
   ! m  = number of coefficients needed to describe the desired combination of PDFs which is evaluated:
   !    = 2 for fitting only a Gaussian
   !    = 5 for fitting a Gaussian plus a Cauchy
   !    = 8 for fitting a Gaussian plus 2 Cauchy distributions
   ! b  = coefficients for the PDFs
   
   ! Output:
   ! f0 = the predicted value at x for the desired combination of PDFs
   
      integer, intent(in) :: m
      real :: arg0
      real :: arg1
      real :: arg2
      real, dimension(8), intent(in) :: b
      real, intent(out) :: f0
      real :: g0
      real :: g1
      real :: g2
      real, parameter :: pii = 0.3183098 ! 1/pi
      real, parameter :: rpi2 = 0.3989422 ! 1/sqrt(2pi)
      real :: t0
      real :: t1
      real :: t2
      real, intent(in) :: x
      
      ! Guassian PDF (see Evans et al., p. 145)
      ! b(1) = location parameter, the mean
      ! b(2) = scale parameter, >0, the standard deviation
   
      arg0 = (x - b(1))/b(2)
      if (abs(arg0) .le. 4.) then
         t0 = rpi2*exp(-.5*arg0*arg0)/b(2) ! Gaussian PDF
      else ! To avoid FPU underflow
         t0 = 0.
      end if
      if (m .gt. 2) go to 1
      g0 = t0
      f0 = g0
      return
         
      ! Cauchy PDF (see Evans et al., p. 48)
      ! b(3) = location parameter a, the median
      ! b(4) = scale parameter b, b>0
      ! b(5) is the relative percentage of this Cauchy
         
    1 arg1 = (x - b(3))/b(4)
      t1 = pii/(b(4)*(1. + arg1*arg1)) ! Cauchy PDF
      g1 = b(5)*t1
      if (m .gt. 5) go to 2
      g0 = (1. - b(5))*t0 ! Scale the Gaussian back
      f0 = g0+g1
      return
         
      ! A second Cauchy PDF
      ! b(6) = location parameter a, the median
      ! b(7) = scale parameter b, b>0
      ! b(8) is the relative percentage of this Cauchy
         
    2 arg2 = (x - b(6))/b(7)
      t2 = pii/(b(7)*(1. + arg2*arg2)) ! Cauchy PDF
      g2 = b(8)*t2
      g0 = (1. - b(5) - b(8))*t0 ! Scale the Gaussian back
      f0 = g0 + g1 + g2
               
      return
      
   end subroutine gswt
   
!*****************************************************************************************
   subroutine phase_utility (iev, ird, psta, p_arrtime)
   
   ! Several steps applied on the basis of phase name:
   !   Setting a flag to determine if phase re-identification can be done
   !   Processing of S-P relative phases
   !   Clean-up of phase names
   !   Processing for relative depth phases (pP-P, pwP-P and sP-P)
   !   Processing PKPdf precursors
   !   Flagging certain phases that are always skipped in mloc
   
      integer :: i
      integer, intent(in) :: iev
      integer, intent(in) :: ird
      real :: delay_time
      real :: dp_arrtime
!      real :: hms2s
      real, intent(in) :: p_arrtime
      character(len=5), intent(in) :: psta
!      logical :: skipp
         
      ! Phase re-identification flag
      phidird(iev,ird) = phid
      if (no_phid(iev,ird) .eq. '!') phidird(iev,ird) = .false. ! This phase ID cannot be changed
      if (fcode(iev,ird) .eq. 'm') phidird(iev,ird) = .false. ! Don't re-identify because we don't know where the station is
      if (fcode(iev,ird) .eq. 't') phidird(iev,ird) = .false. ! Don't re-identify because we don't trust the timing
      if (plogout) write (io_plog,'(a,1x,l1)') 'phase re-identification ', phidird(iev,ird)
      
      ! S-P relative phase
      if (phase0(iev,ird)(1:7) .eq. 'S-P    ') then
         rel_phase(iev,ird) = .true.
         phidird(iev,ird) = .false. ! Don't re-identify this phase
      else
         rel_phase(iev,ird) = .false.
      end if
      if (plogout) write (io_plog,'(a,1x,l1)') 'relative phase ', rel_phase(iev,ird)
      if (plogout) write (io_plog,'(a,1x,l1)') 'phase re-identification ', phidird(iev,ird)
      
      ! Clean up phase name (original saved in phase0)
      call pnclean (phase0(iev,ird), phase(iev,ird))
      if (plogout) then
         if (phase(iev,ird) .ne. phase0(iev,ird)) then
             write (io_plog,'(a,i3,1x,a)') 'pnclean changed '//phase0(iev,ird)//' to '//phase(iev,ird)//' for event ',&
             iev, stname(iev,ird)
         end if
      end if
      
      ! Relative depth phases
      ! pP, pwP and sP phases are converted to pP-P, pwP-P and sP-P times, referenced to the
      ! most recent P arrival from the same station. 
      if (phase0(iev,ird) .eq. 'pP      ' .or. phase0(iev,ird) .eq. 'sP      ' .or. phase0(iev,ird) .eq. 'pwP     ') then
         if (stname(iev,ird) .eq. psta) then
            dp_arrtime = hms2s(ipah(iev,ird), ipam(iev,ird), pas(iev,ird))
            if (ipah(iev,ird) .lt. hourp(iev,0)) dp_arrtime = dp_arrtime + 86400.
            delay_time = dp_arrtime - p_arrtime
            call timecr (delay_time, ipah(iev,ird), ipam(iev,ird), pas(iev,ird)) ! New "arrival time" for the relative phases
            if (phase0(iev,ird) .eq. 'pP      ') then
               phase(iev,ird) = 'pP-P    '
            else if (phase0(iev,ird) .eq. 'sP      ') then
               phase(iev,ird) = 'sP-P    '
            else if (phase0(iev,ird) .eq. 'pwP     ') then
               phase(iev,ird) = 'pwP-P   '
            end if
            rel_depth_phase(iev,ird) = .true. ! This is a relative depth phase
            phidird(iev,ird) = .false. ! Don't re-identify this phase
            if (plogout) write (io_plog,'(a,1x,a,1x,a,1x,f10.3,2i3,f7.3,1x,l1,1x,l1)') stname(iev,ird), phase0(iev,ird),&
             phase(iev,ird), dp_arrtime, ipah(iev,ird), ipam(iev,ird), pas(iev,ird), rel_depth_phase(iev,ird), phidird(iev,ird)
         end if
      end if
      
      ! List of phases that cannot be renamed (command PPRI)
      if (n_no_phreid .gt. 0) then
         do i = 1,n_no_phreid
            if (trim(no_phreid(i)) .eq. trim(phase0(iev,ird))) phidird(iev,ird) = .false. ! This phase ID cannot be changed.
         end do
      end if
      
      ! PKPdf precursors are not used, but their residual will be calculated relative to the PKPdf phase.
      ! They cannot be re-identified.
      if (phase0(iev,ird) .eq. 'PKPdfpre' .or. phase0(iev,ird) .eq. 'PKPpre  ') then
         fcode(iev,ird) = 'p'
         phidird(iev,ird) = .false.
      end if
      
      ! Some phases are always skipped
      if (skipp(phase(iev,ird))) then
         if (plogout) write (io_plog,'(3a)') 'phase_utility: phase always skipped - ', stname(iev,ird), phase(iev,ird)
         fcode(iev,ird) = 'p'
         phidird(iev,ird) = .false.
      end if
         
      return
      
   end subroutine phase_utility
   
         
!*****************************************************************************************
!    subroutine nodupe ()
!    
!    ! Removes duplicate phases from the tau-p phase list (as given in the
!    ! common block /taup/, and returns
!    ! the equivalent phase list in the common block /taupnd/.
!    
!       integer :: i
!       integer :: j
!       integer :: nd
!       logical :: dupe
!    
!       ! Initialize
!       do i = 1,max_mpl
!          phcdnd(i) = '        '
!          nphasend = 0
!          ttnd(i) = 0.
!       end do
!    
!       if (nphase .eq. 1) then
!          phcdnd(1) = phcd(1)
!          nphasend = 1
!          ttnd(1) = tt(1)
!       else if (nphase .gt. 1) then
!          phcdnd(1) = phcd(1)
!          ttnd(1) = tt(1)
!          nd = 1
!          do i = 2,nphase
!             dupe = .false.
!             do j = 1,nd
!                if (phcdnd(j) .eq. phcd(i)) dupe = .true.
!             end do
!             if (.not.dupe) then
!                nd = nd + 1
!                phcdnd(nd) = phcd(i)
!                ttnd(nd) = tt(i)
!             end if
!          end do
!          nphasend = nd
!       end if
!    
!       return
!       
!    end subroutine nodupe
         
         
!*****************************************************************************************
   character(len=8) function iscphase (n)
         
   !  Given a phase code number, returns the ISC phase identification.
   
      integer, intent(in) :: n
      character(len=8) :: phase(0:100)
      
      data phase/'P       ','PP      ','PPP     ','PCP     ','PKP     ',&
                 'PKP2    ','PKPPKP  ','PCPPKP  ','PS      ','PPS     ',&
                 'PCS     ','PKS     ','PKKS    ','PCSPKP  ','PKPPKS  ',&
                 'PKPSKS  ','PKKP    ','3PKP    ','PKIKP   ','PP2     ',&
                 'PPP2    ','PKS2    ','PSS     ','PSS2    ','SSP2    ',&
                 'PCPPKP2 ','PCSPKP2 ','SS2     ','PKKP2   ','PKKS2   ',&
                 'SCSPKP3 ','SCSPKP2 ','SCSP2   ','SKSP2   ','SSS2    ',&
                 'S       ','SS      ','SSS     ','SCS     ','SKS     ',&
                 'SKKS    ','SKKKS   ','SCSPKP  ','SKSSKS  ','SCSP    ',&
                 'SKSP    ','SCP     ','SP      ','SKP     ','SKKP    ',&
                 'SKPPKP  ','SSP     ','SKP2    ','SKS2    ','SKKS2   ',&
                 'SKKS3   ','SKKKS2  ','*SPKP2  ','*PPCP   ','*PPKP   ',&
                 '*PP     ','*PPP    ','*SP     ','*SPKP   ','*SS     ',&
                 '*SSS    ','*SPP    ','*SPCP   ','*SSCS   ','*PPKP2  ',&
                 'P*      ','S*      ','PG      ','SG      ','PN      ',&
                 'SN      ','PGPG    ','SGSG    ','LR      ','LQ      ',&
                 'L       ','PKKP3   ','PKKS3   ','SPP     ','PHASE84 ',&
                 'P DIFF  ','QM      ','RM      ','T       ','T(max)  ',&
                 'NORTH   ','SOUTH   ','EAST    ','WEST    ','UP      ',&
                 'DOWN    ','E       ','I       ','MAXIMUM ','FINAL   ',&
                 '        '/
     
      iscphase = phase(n)
      
      return
      
   end function iscphase
        
   
!*****************************************************************************************
   subroutine pnclean (pin, pout)
   
   ! Clean up phase names.
   ! Remove initial 'X', 'e', 'i', and 'q' and trailing 'c', ,'d' or 'r' from phase name.
   ! Remove the !-flag in character position 8 that fixes phase ID.
   ! Also removes parentheses
   
!      integer :: lenb
      integer :: j
      integer :: k
      integer :: i
      character(len=8), parameter :: b8 = ' ' ! blank phase name
      character(len=1) :: c
      character(len=8), intent(in) :: pin
      character(len=8), intent(out) :: pout
      logical :: trim
         
      pout = pin
      if (pout .eq. b8 .or. pout .eq. 'UNK     ') then
         pout = 'UNKNOWN '
         return
      end if
      if (pout(8:8) .eq. '!') pout(8:8) = ' '
      if (pout .eq. 'PKPpre  ') pout = 'PKPdfpre'
      if (pout .eq. 'PKPdfpre') return
      if (pout .eq. 'S-P     ') return
   
      j = lenb(pout)
      if (j .ge. 1) then
   
         ! Trim junk off the beginning
         trim = .true.
         do while (trim)
            if (pout(1:1) .eq. 'X' .or.& 
                pout(1:1) .eq. ' ' .or.& 
                pout(1:1) .eq. 'e' .or.& 
                pout(1:1) .eq. 'E' .or.& 
                pout(1:1) .eq. 'i' .or.& 
                pout(1:1) .eq. 'I' .or.& 
                pout(1:1) .eq. 'q' .or.& 
                pout(1:1) .eq. '(') then
               if (j .ge. 2) then
                  do i = 2,j
                     k = i - 1
                     pout(k:k) = pout(i:i)
                  end do
                  pout(j:j) = ' '
                  j = j - 1
               else
                  pout = 'UNKNOWN '
                  return
               end if
            else
               trim = .false.
            end if
         end do
            
         ! Trim junk off the end
         if (j .ge. 2) then
            k = j - 1
            if (pout(k:j) .ne. 'bc' .and. pout(k:j) .ne. 'ac') then ! Watch out for ac and bc branches
               c = pout(j:j)
               if (c .eq. 'c' .or. c .eq. 'd' .or. c .eq. 'r' .or. c .eq. ')') pout(j:j) = ' '
            end if
         else if (j .eq. 1) then
            c = pout(1:1)
            if (c .ne. 'P') return
            if (c .ne. 'S') return
            if (c .ne. 'L') return
            if (c .ne. 'M') return
            if (c .ne. 'Q') return
            if (c .ne. 'R') return
            pout = 'UNKNOWN '
            return
         end if
      
         j = lenb(pout)
      
         ! Fix second character of some crustal phases
         if (j .eq. 2) then
            if (pout(2:2) .eq. 'N') pout(2:2) = 'n'
            if (pout(2:2) .eq. 'G') pout(2:2) = 'g'
            if (pout(2:2) .eq. 'B') pout(2:2) = 'b'
            if (pout(2:2) .eq. '*') pout(2:2) = 'g'
         end if
      
         ! Outer core reflections
         if (j .eq. 3 .and. pout(2:2) .eq. 'C') pout(2:2) = 'c'
   
         ! Inner-core reflections
         if (j .eq. 5 .and. pout(3:3) .eq. 'I') pout(3:3) = 'i'
      
         ! Fix any remaining blanks inside the phase name (trailing blanks are OK)
         do i = 1,j
            if (pout(i:i) .eq. ' ') pout(i:i) = '_'
         end do
      
      end if
   
      return
      
   end subroutine pnclean
   
   
!*****************************************************************************************
   logical function phase_ok (phasein)
   
   ! Catches "junk phases" that have nothing to do with the location problem
   
      integer :: i
      character(len=8), intent(in) :: phasein
      
      phase_ok = .true.
      if (n_junk .ge. 1) then
         do i=1,n_junk
            if (phasein .eq. junk_phase(i)) then
               phase_ok = .false.
               exit
            end if
         end do
      end if
      
      return
      
   end function phase_ok
   
   
!*****************************************************************************************
end module mloc_phases

