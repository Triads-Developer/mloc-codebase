!> Direct calibration

module dcal

   use declare_limits
   use declare_lun, only: io_cal
   use declare_calibration
   
   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine dircal (it)
      
      use declare_limits
      
      implicit none
      
      integer :: i
      integer :: iev
      integer :: igt
      integer :: it
      integer :: j
      integer :: kbayes
      integer :: mf
      integer :: ncum
      integer :: nf
      integer, dimension(0:11) :: ngt
      integer :: ngt5
      real :: al
      real :: alpha
      real :: area
      real :: avc
      real :: bl
      real, dimension(2,2) :: cv22
      real, allocatable, dimension(:,:,:) :: dccv
      real :: ddep
      real :: dot
      real :: f
      real :: kcrit1
      real :: kcrit2
      real :: pc
      real, allocatable, dimension(:,:,:) :: sccv
      real, dimension(4,4) :: shcv
      real :: t11
      real :: t12
      real :: t22
      real :: xl
      real :: yl
      character(len=100) :: outfil
      
      data pc/0.90/
      
      outfil = trim(outfile)//'.dcal'
      call open_file (io_cal, outfil, 'new')
      
      sccv = 0.
      shcv = 0.
      dccv = 0.
      ngt = 0
      ncum = 0
      ngt5 = 0
      avc = 0.
      
      do iev = 1,nev
      
         ! Scaled covariance matrices
         if (mindx(iev,1) .ne. 0 .and. mindx(iev,2) .ne. 0) then
            mf = 2
            kbayes = 3
            nf = kbayes + ndatc(iev,it) - mf
            call fstat1 (mf, nf, pc, f)
            kcrit2 = mf*((real(kbayes)+eciev(iev,it))/real(nf))*f 
            cv22(1,1) = ccv(iev,1,1)
            cv22(1,2) = ccv(iev,1,2)
            cv22(2,1) = ccv(iev,1,2)
            cv22(2,2) = ccv(iev,2,2)
            call elips (cv22, alpha, al, bl)
            xl = sqrt(kcrit2/al)
            yl = sqrt(kcrit2/bl)
            avc = xl*yl*pi ! Geometric area of 90% confidence ellipse for cluster vector
            
            ! Convert 90% confidence ellipse back to matrix form under the asumption of kcrit = 1.
            call ell2cv (1.0, alpha, xl, yl, t11, t12, t22)
            sccv(iev,1,1) = t11
            sccv(iev,1,2) = t12
            sccv(iev,2,1) = t12
            sccv(iev,2,2) = t22
         end if
         
         ! Statistics for univariate parameters
         mf = 1
         kbayes = 3
         nf=kbayes+ndatc(iev,it)-mf
         call fstat1 (mf, nf, pc, f)
         kcrit1 = mf*((real(kbayes)+eciev(iev,it))/real(nf))*f 
!         write (io_cal,'(i3,3x,a,f10.3)') iev, 'kcrit1 = ', kcrit1
         
         ! Depth uncertainty
         if (mindx(iev,3) .ne. 0) then
            ddep = sqrt(kcrit1*ccv(iev,3,3))
            sccv(iev,3,3) = ddep*ddep
         end if
         
         ! Origin time uncertainty
         if (mindx(iev,4) .ne. 0) then
            dot = sqrt(kcrit1*ccv(iev,4,4))
            sccv(iev,4,4) = dot*dot
         end if
         
      end do ! End of loop over events
      
      ! Scaled hypocentroid covariance matrix
      call ell2cv (1.0, alphah, xl1h, xl2h, t11, t12, t22)
      shcv(1,1) = t11
      shcv(1,2) = t12
      shcv(2,1) = t12
      shcv(2,2) = t22
      shcv(3,3) = sdxhath(3)*sdxhath(3)
      shcv(4,4) = sdxhath(4)*sdxhath(4)
      
      write (io_cal,'(/a)') 'Cumulative uncertainty of direct-calibrated events'
      write (io_cal,'(a)') 'iev     alpha        xl        yl      area       eqr      ddep       dot'
      
      do iev = 1,nev
         if (mindx(iev,1) .ne. 0 .and. mindx(iev,2) .ne. 0) then
            do i = 1,4
               do j = 1,4
                  dccv(iev,i,j) = shcv(i,j) + sccv(iev,i,j)
               end do
            end do
   !         write (io_cal,'(i3,4f10.3)') iev, dccv(iev,1,1), dccv(iev,1,2), dccv(iev,2,1), dccv(iev,2,2)
            cv22(1,1) = dccv(iev,1,1)
            cv22(1,2) = dccv(iev,1,2)
            cv22(2,1) = dccv(iev,2,1)
            cv22(2,2) = dccv(iev,2,2)
            call elips (cv22, alpha, al, bl)
            xl = sqrt(1./al)
            yl = sqrt(1./bl)
            igt = nint(yl)
            if (igt .ge. 0 .and. igt .le. 10) then
               ngt(igt) = ngt(igt) + 1
            else
               ngt(11) = ngt(11) + 1 ! This holds the number of events worse than CE10
            end if
            area = xl*yl*pi ! Geometric area
            eqr(iev) = sqrt(area/pi)
            ddep = sqrt(dccv(iev,3,3))
            dot = sqrt(dccv(iev,4,4))
            write (io_cal,'(i3,7f10.3)') iev, alpha, xl, yl, area, eqr(iev), ddep, dot
            alphadc(iev) = alpha
            xl1dc(iev) = xl
            xl2dc(iev) = yl
            ddepdc(iev) = ddep
            dotdc(iev) = dot
         else if (mindx(iev,3) .ne. 0 .and. mindx(iev,4) .ne. 0) then ! Depth and OT are free parameters
            ddep = sqrt(dccv(iev,3,3))
            dot = sqrt(dccv(iev,4,4))
            write (io_cal,'(i3,50x,2f10.3)') iev, ddep, dot
            alphadc(iev) = 0.
            xl1dc(iev) = 99.
            xl2dc(iev) = 99.
            ddepdc(iev) = ddep
            dotdc(iev) = dot
         else if (mindx(iev,4) .ne. 0) then ! Only OT is free
            dot = sqrt(dccv(iev,4,4))
            write (io_cal,'(i3,60x,f10.3)') iev, dot
            alphadc(iev) = 0.
            xl1dc(iev) = 99.
            xl2dc(iev) = 99.
            ddepdc(iev) = 99.
            dotdc(iev) = dot
         end if
      end do
      
      ! Summary of calibration levels
      if (mindx(1,1) .ne. 0 .and. mindx(1,2) .ne. 0) then
         write (*,'(/a)') 'CE    N   Cumulative'
         write (io_cal,'(/a)') 'CE    N   Cumulative'
         do i=0,11
            ncum = ncum + ngt(i)
            if (i .eq. 5) ngt5 = ncum
            write (*,'(i2,2x,i3,7x,i4)') i, ngt(i), ncum
            write (io_cal,'(i2,2x,i3,7x,i4)') i, ngt(i), ncum
         end do
         write (*,'(/a,i3,a,i3,a/)') 'Direct calibration: ', ngt5, ' events out of ', nev, ' are CE05 or better'
         write (io_cal,'(/a,i3,a,i3,a/)') 'Direct calibration: ', ngt5, ' events out of ', nev, ' are CE05 or better'
      end if
      
      close (io_cal)
      
      return
      
   end subroutine dircal
   
   
!*****************************************************************************************
   subroutine dcal_allocate ()
   
   ! Allocate variable arrays related to direct calibration
   
      integer :: error
      integer :: n
      character(len=32) :: p
      
      n = nevmax
      p = 'dcal_allocate' ! procedure
      
      allocate (dccv(n,4,4), stat=error); if (error .gt. 0) call allocation_error (p, 'dccv', error)
      allocate (sccv(n,4,4), stat=error); if (error .gt. 0) call allocation_error (p, 'sccv', error)
      
      return
   
   end subroutine dcal_allocate   
   

!*****************************************************************************************
end module dcal

