!> Declaration of variables related to logical unit numbers

module declare_lun

   implicit none
   save

   integer :: io_bayes = 11 ! BAYESLOC output
   integer :: io_bdp = 12 ! Bad depth phase station list
   integer :: io_bptc_log = 13 ! Bounce point correction log file
   integer :: io_cal = 14 ! .cal file
   integer :: io_cfil = 15 ! command file
   integer :: io_dat0 = 16 ! .dat0 file
   integer :: io_datf = 17 ! .datf file
   integer :: io_depth_phase = 18 ! depth phases output
   integer :: io_diffdat = 19 ! Differential time data
   integer :: io_dpnc = 20 ! Output file for depth phase name changes
   integer :: io_gccel = 21 ! GCCEL output
   integer :: io_gmt = 22 ! GMT scripts
   integer :: io_ifsm = 23 ! ISC-FDSN station metadata file
   integer :: io_in = 24 ! Event data files
   integer :: io_junk = 25 ! Input file of junk phase names that will not be read
   integer :: io_locmod = 26 ! Custom crustal velocity model
   integer :: io_log = 27 ! log file
   integer :: io_lres = 28 ! .lres file (large cluster residuals)
   integer :: io_nsmd = 29 ! NEIC station metadata file
   integer :: io_oldr = 30 ! Output for limited distance range
   integer :: io_out = 31 ! Used for several output files
   integer :: io_pdf = 32 ! Probability density function file
   integer :: io_plog = 33 ! phase re-identification log
   integer :: io_rderr = 34 ! .rderr (reading errors file)
   integer :: io_rhdf = 35 ! HDF file to set starting locations
   integer :: io_shifted = 36 ! Shifted phase_data after indirect calibration
   integer :: io_stn_log = 37 ! station data log
   integer :: io_taup = 38 ! .hed and .tbl files (ak135)
   integer :: io_taup_log = 39 ! Log file for tau-p debugging
   integer :: io_tt = 40 ! Empirical TT data for specific phases
   integer :: io_ttsprd = 41 ! .ttsprd file
   integer :: io_xdat = 42 ! .xdat file
   
contains

!*****************************************************************************************
   integer function lunit ()
   
   ! find an unopened fortran unit number between 10 and 100.
         
      logical :: lopen
      
      do lunit = 43,99
         inquire (unit=lunit,opened=lopen)
         if (.not.lopen) return
      end do
      lunit = 0
      write (*,'(t3,a)') 'Oops! Fatal error in lunit: failed to find an unopened unit between 43 and 99'
      stop
      
   end function lunit

      
!*****************************************************************************************
end module declare_lun

