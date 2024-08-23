!> Declaration and allocation of variables and code related to the ISC-FDSN Station
!>  Metadata file and the run-specific list of matched entries.

module mloclib_ifsm

   use declare_limits
   use declare_lun
   use mloclib_date_time
   
   implicit none
   save
   
   integer :: n_ifsm_read ! Entries read from the main ISC-FDSN station metadata file
   integer :: n_ifsm_written ! Matched entries written to the run-specific output file
   integer, allocatable, dimension(:) :: ifsm_off ! Operational epoch off, Year + Julian date
   integer, dimension(n_ifsm_matched_max) :: ifsm_off_matched ! Matched entries, operational epoch off, Year + Julian date
   integer, allocatable, dimension(:) :: ifsm_on ! Operational epoch on, Year + Julian date
   integer, dimension(n_ifsm_matched_max) :: ifsm_on_matched ! Matched entries, operational epoch on, Year + Julian date
   character(len=7), allocatable, dimension(:) :: ds ! Deployment, Station Code, main file
   character(len=7), dimension(n_ifsm_matched_max) :: ds_matched ! Deployment, Station Code, matched entries
   character(len=100) :: ifsm_stn_file
   character(len=94), allocatable, dimension(:) :: ifsm_line
    
contains

!*****************************************************************************************
!> Read the main ISC-FDSN station metadata file

   subroutine read_ifsm (io_ifsm, filename, outfile)
      
      integer :: day_oe1 ! operational epoch day on
      integer :: day_oe2 ! operational epoch day off
      integer :: allocate_error ! allocation status variable
      integer :: i
      integer, intent(in) :: io_ifsm ! logical unit number
      integer :: ios ! I/O status
      integer :: jdate ! Julian date
      integer :: month_oe1 ! operational epoch month on
      integer :: month_oe2 ! operational epoch month off
      integer :: n_entries ! number of entries in the main ISC-FDSN station metadata file
      integer :: year_oe1 ! operational epoch year on
      integer :: year_oe2 ! operational epoch year off
      character(len=*), intent(in) :: filename ! main ISC-FDSN station metadata filename
      character(len=94) :: line94
      character(len=132) :: msg
      character(len=*) :: outfile ! run-specific output filename

      !Initialize
      n_ifsm_written = 0
      
      call open_file (io_ifsm, filename, 'old')
      
      ! Count the entries
      n_entries = 0
      do
         read (io_ifsm,'(a)',iostat=ios) line94
         if (ios .lt. 0) exit
         n_entries = n_entries + 1
      end do
      rewind (io_ifsm)
      
      ! Allocate variable arrays
      allocate (ifsm_off(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_ifsm: space for ', n_entries,&
          ' for array "ifsm_off" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      allocate (ifsm_on(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_ifsm: space for ', n_entries,&
          ' for array "ifsm_on" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      allocate (ds(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_ifsm: space for ', n_entries,&
          ' for array "ds" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      allocate (ifsm_line(n_entries), stat=allocate_error)
      if (allocate_error .gt. 0) then
         write (msg,'(a,i8,a,i2)') 'read_ifsm: space for ', n_entries,&
          ' for array "ifsm_line" could not be allocated; stat = ', allocate_error
         call oops (trim(msg))
      end if
      
      read (io_ifsm,'(a)') line94 ! First line is the file creation date    
      do i = 1, n_entries-1
         read (io_ifsm,'(a)',iostat=ios) line94
         ifsm_line(i) = line94
         ds(i) = line94(2:3)//line94(5:9) ! Deployment, Station Code
         if (line94(49:49) .ne. ' ') then
            read (line94(49:58),'(i4,1x,i2,1x,i2)') year_oe1, month_oe1, day_oe1
         else
            year_oe1 = 1900
            month_oe1 = 1
            day_oe1 = 1
         end if
         if (line94(75:75) .ne. ' ') then
            read (line94(75:84),'(i4,1x,i2,1x,i2)') year_oe2, month_oe2, day_oe2
         else
            year_oe2 = 2099
            month_oe2 = 12
            day_oe2 = 31
         end if
         call juldat (year_oe1, month_oe1, day_oe1, jdate, 0)
         ifsm_on(i) = year_oe1*1000 + jdate
         call juldat (year_oe2, month_oe2, day_oe2, jdate, 0)
         ifsm_off(i) = year_oe2*1000 + jdate
      end do
      n_ifsm_read = i
      write (msg,'(a,i8,a)') 'read_ifsm: ', n_ifsm_read, ' ISC-FDSN station metadata entries read: '
      call fyi (trim(msg))
      close (io_ifsm)
      
      ! Open a new file for the selected station entries
      ifsm_stn_file = trim(outfile)//'_ifsm_stn.dat'
      call open_file (io_ifsm, ifsm_stn_file, 'new')
      write (io_ifsm,'(a)') '3 Supplemental station file from ISC-FDSN Station Metadata'
      
      return

   end subroutine read_ifsm
   

!*****************************************************************************************
   subroutine write_ifsm (io_ifsm, ds_test, jdate)

   ! Generate a supplemental station file entry from the ISC-FDSN station metadata file
   ! for stations missing from the current station list (supplemental lists + master list).
   ! The output file is in the "generic" (type 3) format.
   ! A match, which causes a new entry in the supplemental station file, is based on
   ! "deployment/station code" from the MNF data file, and the date of the event must
   ! be within the operational date range of the entry in the ISC-FDSN station metadata
   ! file. At this time no operational dates are carried in the file so the default
   ! (1900-2099) is used, meaning there is no check for operational epoch.

      integer :: i
      integer :: elev_ifsm
      integer, intent(in) :: io_ifsm
      integer, intent(in) :: jdate
      real :: lat_ifsm, lon_ifsm
      character(len=5) :: agency_ifsm
      character(len=7), intent(in) :: ds_test
      character(len=132) :: msg
      logical :: matched, date_range

      ! First check against entries that have already been matched in earlier events
      if (n_ifsm_written .ge. 1) then
         do i = 1,n_ifsm_written
            if (ds_test .ne. ds_matched(i)) cycle ! Network/deployment code
            if (date_range(jdate,ifsm_on_matched(i),ifsm_off_matched(i))) then
               matched = .true.
               exit
            end if
         end do
      end if
      
      ! Test against the full ISC-FDSN data file
      if (.not.matched) then
         if (n_ifsm_written .lt. n_ifsm_matched_max) then
            do i = 1, n_ifsm_read
               if (ds_test .eq. ds(i)) then
                  if (date_range(jdate,ifsm_on(i),ifsm_off(i))) then
                     n_ifsm_written = n_ifsm_written + 1
                     read (ifsm_line(i)(12:19),'(f8.4)') lat_ifsm
                     read (ifsm_line(i)(21:29),'(f9.4)') lon_ifsm
                     read (ifsm_line(i)(31:35),'(i5)') elev_ifsm
                     if (ds(i)(1:2) .eq. 'IR') then
                        agency_ifsm = 'ISC  '
                     else
                        agency_ifsm = 'FDSN '
                     end if
                     call supp_station_file_write (io_ifsm, ds(i)(3:7), agency_ifsm,&
                      ds(i)(1:2)//'      ', lat_ifsm, lon_ifsm, elev_ifsm, 0,&
                      ifsm_on(i), ifsm_off(i), ' ') 
                     ds_matched(n_ifsm_written) = ds(i)
                     ifsm_on_matched(n_ifsm_written) = ifsm_on(i)
                     ifsm_off_matched(n_ifsm_written) = ifsm_off(i)
                     exit
                  end if
               end if
            end do
         else if (n_ifsm_written .eq. n_ifsm_matched_max) then
            write (msg,'(a,i4,a)') 'read_mnf_134: maximum number of matched ISC-FDSN metadata entries (',&
             n_ifsm_matched_max, ') written'
            call warnings (trim(msg))
         end if
      end if
      
      return
            
   end subroutine write_ifsm
   
   
!*****************************************************************************************
end module mloclib_ifsm

