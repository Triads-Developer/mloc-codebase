!> Procedure for ~.ttsprd output file

module mlocout_ttsprd

   use declare_configuration
   use declare_environment
   use declare_lun
   use declare_phase_data
   use declare_phases
   use declare_stations
   use mloclib_statistics
   
   implicit none
   save
   
contains

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
      integer, external :: lunit
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
      
      do iev = 1,nev
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
end module mlocout_ttsprd

