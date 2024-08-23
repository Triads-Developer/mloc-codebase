!> Declaration and allocation of variables and code related to the NEIC Station Metadata
!>  file and the run-specific list of matched entries.

module mloclib_nsmd

   use declare_limits
   use declare_lun
   use mloclib_date_time
   
   implicit none
   save
   
   integer :: n_nsmd_read ! Entries read from the main NEIC station metadata file
   integer :: n_nsmd_written ! Matched entries written to the run-specific output file
   integer, allocatable, dimension(:) :: nsmd_off ! Operational epoch off, Year + Julian date
   integer, dimension(n_nsmd_matched_max) :: nsmd_off_matched ! Matched entries, operational epoch off, Year + Julian date
   integer, allocatable, dimension(:) :: nsmd_on ! Operational epoch on, Year + Julian date
   integer, dimension(n_nsmd_matched_max) :: nsmd_on_matched ! Matched entries, operational epoch on, Year + Julian date
   character(len=12), allocatable, dimension(:) :: dslc ! Deployment, Station Code, Location, Channel, main file
   character(len=12), dimension(n_nsmd_matched_max) :: dslc_matched ! Deployment, Station Code, Location, Channel, matched entries
   character(len=100) :: nsmd_stn_file
   character(len=102), allocatable, dimension(:) :: nsmd_line
    
contains

!*****************************************************************************************
!> Read the main NEIC station metadata file
   
   subroutine read_nsmd (io_nsmd, filename, outfile)
   
      integer :: day_oe1 ! operational epoch day on
      integer :: day_oe2 ! operational epoch day off
      integer :: allocate_error ! allocation status variable
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
      
      ! Allocate variable arrays
      allocate (nsmd_off(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_nsmd: space for ', n_entries,&
          ' for array "nsmd_off" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      allocate (nsmd_on(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_nsmd: space for ', n_entries,&
          ' for array "nsmd_on" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      allocate (dslc(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_nsmd: space for ', n_entries,&
          ' for array "dslc" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      allocate (nsmd_line(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_nsmd: space for ', n_entries,&
          ' for array "nsmd_line" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
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
      logical :: matched, date_range

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
         if (n_nsmd_written .lt. n_nsmd_matched_max) then
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
         else if (n_nsmd_written .eq. n_nsmd_matched_max) then
            write (msg,'(a,i4,a)') 'read_mnf_134: maximum number of matched NEIC metadata entries (',&
             n_nsmd_matched_max, ') written'
            call warnings (trim(msg))
         end if
      end if
      
      return
   
   end subroutine write_nsmd
   

!*****************************************************************************************
end module mloclib_nsmd

