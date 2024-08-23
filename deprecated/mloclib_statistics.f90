!> Statistical procedures

module mloclib_statistics

   use declare_lun
   use mloclib_messages

   implicit none
   save
   
contains

!*****************************************************************************************
   subroutine sort (n, ra)
   
   ! Heapsort, from Numerical Recipes p. 231
   
      integer :: i
      integer :: ir
      integer :: j
      integer :: l
      integer, intent(in) :: n
      real, intent(inout), dimension(n) :: ra
      real :: rra
   
      l = n/2 + 1
      ir = n
   10 continue
      if (l .gt. 1) then
         l = l - 1
         rra = ra(l)
      else
         rra = ra(ir)
         ra(ir) = ra(1)
         ir = ir -1
         if (ir .eq. 1) then
            ra(1) = rra
            return
         end if
      end if
      i = l
      j = l + l
   20 if (j .le. ir) then
         if (j .lt. ir) then
            if (ra(j) .lt. ra(j+1)) j = j + 1
         end if
         if (rra .lt. ra(j)) then
            ra(i) = ra(j)
            i = j
            j = j + j
         else
            j = ir + 1
         end if
         go to 20
      end if
      ra(i) = rra
      go to 10
   
   end subroutine sort
   
   
!*****************************************************************************************
   subroutine mdian1 (x, n, xmed)
   
   ! Median, from Numerical Recipes, p. 460
   
      integer, intent(in) :: n
      integer :: n2
      real, dimension(n), intent(inout) :: x
      real, intent(out) :: xmed
   
      call sort (n,x)
      n2 = n/2
      if (2*n2 .eq. n) then
         xmed = 0.5*(x(n2) + x(n2 + 1))
      else
         xmed = x(n2 + 1)
      end if
   
      return
      
   end subroutine mdian1
   
   
!*****************************************************************************************
   subroutine fstat1 (nu1, nu2, p0, f)
   
   !  Returns F(nu1,nu2) for confidence level p0.  This value of F will be
   !  exceeded with a probability of 1-p0.
   !  This version of fstat is the original one used since the initial coding of mloc
   
      integer ::  i
      integer :: j
      integer, intent(in) :: nu1
      integer, intent(in) :: nu2
      real :: a1
      real :: a2
      real :: betai
      real :: diff
      real, intent(out) :: f
      real :: factor
      real :: ftest
      real, dimension(101) :: p
      real, intent(in) :: p0
      real :: x
      real :: xnu1
      real :: xnu2
      character(len=132) :: msg
   
      a1 = real(nu1)/2.
      a2 = real(nu2)/2.
      xnu1 = real(nu1)
      xnu2 = real(nu2)
      f = 1.
      factor = 1.
      diff = 1.
      do while (diff .gt. .001)
         do i = 0,100
            ftest = f + (i*factor)
            x = xnu2/(xnu2+(xnu1*ftest))
            if (x .lt. 0. .or. x .gt. 1.) then
               write (msg,'(4(a,f10.3,2x))') 'fstat1: x = ', x, 'xnu2 = ', xnu2, 'xnu1 = ',&
                xnu1, 'ftest = ', ftest
               call warnings (msg)
            end if
            p(i+1) = 1. - betai(a2, a1, x)
         end do
         call locate (p, 101, p0, j)
         if (j .lt. 101) then
            diff = abs(p0-p(j))
            f = f + ((j-1)*factor)
            factor = factor*0.01
         else
            factor = factor*100.
         end if
      end do
   
      return
      
   end subroutine fstat1
   
   
!*****************************************************************************************
   subroutine locate (xx, n, x, j)
         
      integer, intent(out) :: j
      integer :: jl
      integer :: jm
      integer :: ju
      integer, intent(in) :: n
      real, intent(in) :: x
      real, dimension(n), intent(in) :: xx
   
      jl = 0
      ju = n + 1
   10 if (ju-jl .gt. 1) then
         jm = (ju+jl)/2
         if ((xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
         go to 10   
      end if
      j = jl
   
      return
      
   end subroutine locate
   
   
!*****************************************************************************************
   real function gammln (xx)
   
      integer :: j
      real, intent(in) :: xx
      double precision, parameter, dimension(6) :: cof = &
       (/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5/)
      double precision, parameter :: fpf = 5.5d0
      double precision, parameter :: half = 0.5d0
      double precision, parameter :: one = 1.0d0
      double precision :: ser
      double precision, parameter :: stp = 2.50662827465d0
      double precision :: tmp
      double precision :: x
   
      x = xx - one
      tmp = x + fpf
      tmp = (x+half)*log(tmp) - tmp
      ser = one
      do j = 1,6
         x = x + one
         ser = ser + cof(j)/x
      end do
      gammln = sngl(tmp + log(stp*ser))
   
      return
      
   end function gammln
   
   
!*****************************************************************************************
   real function betai (a, b, x)
   
   ! Incomplete beta function, based on BETAI (Section 6.3, Numerical Recipes). Modified to
   ! avoid FPU exceptions, equality tests for real numbers, etc.
         
      real, intent(in) :: a
      real :: arg
      real, intent(in) :: b
      real :: betacf
      real :: bt
      real :: eps
      real :: gammln
      real, intent(in) :: x
      character(len=132) :: msg
      
      eps = 2.*epsilon(x)
   
      ! a,b > 0
      if (a .lt. eps .or. b .lt. eps .or. x .lt. 0. .or. x .gt. 1.) then
         write (msg,'(a,3e12.4)') 'betai: bad argument: a, b, x = ', a, b, x
         call oops (trim(msg))
      end if
      if (x .lt. eps .or. (1.-x) .lt. eps) then
         bt = 0.
      else
         arg = gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.-x)
         if (arg .lt. -10.) then
            bt = 0.
         else
            bt = exp(arg)
         end if
      end if
      if (x .lt. (a + 1.)/(a + b + 2.)) then
         betai = bt*betacf(a,b,x)/a
      else
         betai = 1. - bt*betacf(b,a,1.-x)/b
      end if
      
      return
   
   end function betai
   
   
!*****************************************************************************************
   real function betacf (a, b, x)
   
   ! Continued fraction for incomplete beta function, based on BETACF (section 6.3 of Numerical
   ! Recipes). Modified from original code to deal better with pathological situations, i.e.,
   ! failure to converge.
   
      integer, parameter :: itmax = 200
      integer :: m
      real, intent(in) :: a
      real :: am
      real :: aold
      real :: ap
      real :: app
      real :: az
      real, intent(in) :: b
      real :: bm
      real :: bp
      real :: bpp
      real :: bz
      real :: d
      real :: em
      real :: eps
      real :: qab
      real :: qam
      real :: qap
      real :: tem
      real, intent(in) :: x
      character(len=132) :: msg
      
      eps = 2.*epsilon(x)
   
      am = 1.
      bm = 1.
      az = 1.
      qab = a + b
      qap = a + 1.
      qam = a - 1.
      bz = 1. - qab*x/qap
      do m = 1,itmax
         em = real(m)
         tem = em + em
         d = em*(b-em)*x/((qam+tem)*(a+tem))
         ap = az + d*am
         bp = bz + d*bm
         d = -(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
         app = ap + d*az
         bpp = bp + d*bz
         aold = az
         am = ap/bpp
         bm = bp/bpp
         az = app/bpp
         bz = 1.
         betacf = az
         if (abs(az-aold) .lt. eps*abs(az)) return
      end do
      
      write (io_log,'(a,3f10.3,2e10.3)') 'betacf: ', a, b, x, betacf, aold
      call warnings ('betacf: a or b too big, or itmax too small; check the log file for details')
   
      return
      
   end function betacf
   
   
!*****************************************************************************************
   subroutine moment2 (data_in, n, adev, sdev)
   
   ! Basic statistical parameters for an input array. Adapted from
   ! Numerical Recipes subroutine MOMENT
         
      integer :: j 
      integer, intent(in) :: n
      real, intent(out) :: adev
      real :: ave
      real, dimension(n), intent(in) :: data_in
      real :: p
      real :: s
      real, intent(out) :: sdev
      real :: var
      character(len=132) :: msg
      
      if (n .le. 1) then
         write (msg,'(a,i6)') 'moment2: illegal value for n: ', n
         call warnings (trim(msg))
         adev = 1.0
         sdev = 1.0
         return
      end if 
      
      s = 0.
      do j = 1,n
         s = s + data_in(j)
      end do
      
      ave=s/n
      adev=0.
      var=0.
      do j = 1,n
         s = data_in(j) - ave
         adev = adev + abs(s)
         p = s*s
         var = var + p
      end do
      adev = adev/n
      var = var/(n-1)
      sdev = sqrt(var)
      
      return
      
   end subroutine moment2
   
   
!*****************************************************************************************
   subroutine croux (x, nin, sn)
   
   ! Implementation of the "naive" algorithm for Sn in:
   ! "Time Efficient algorithms for two highly robust estimators of scale" by Croux & Rousseuw
   ! Computational Statistics, V1 (1992), Dodge and Whittaker, ed., Physica-Verlag, Heidleberg, pp. 411-428.
   ! Sn is a robust estimator for the spread of a sample distribution that does not need an estimate
   ! of central location. It is well-behaved even with small sample size.
   
   ! The original formulation behaves badly with n=3 and two values close to each other. 
   ! the trivial difference (i = j) always yeilds a zero value and if one of the other
   ! differences is also small, the estimate of Sn implodes. After consulting
   ! Christophe Croux, I altered the algorithm for n=3 so that the inner loop takes the
   ! average instead of the lomed of the three differences.
   
   ! It follows that the constant cn(3) should be recalculated. I did the same experiment as
   ! reported in Croux & Rousseuw with the new algorithm and found cn(3) = 1.172.
         
      integer, parameter :: nmax = 1000 ! Maximum size of input array "x"
      integer :: i
      integer :: ihimed
      integer :: ilomed
      integer, dimension(nmax) :: indx
      integer :: j
      integer :: n
      integer, intent(in) :: nin
      real, dimension(nmax) :: a1
      real, dimension(nmax) :: a2
      real, intent(out) :: sn
      real, external :: cn
      real, dimension(nin), intent(in) :: x
      character(len=132) :: msg
      
      if (nin .le. nmax .and. nin .ge. 2) then
         n = nin
      else if (nin .lt. 2) then
         write (msg,'(a,i6)') 'croux: illegal value for n: ', nin
         call warnings (trim(msg))
         sn = 1.0
         return   
      else if (nin .gt. nmax) then
         write (msg,'(a,i6)') 'croux: nin exceeds maximum value, set to ', nmax
         call fyi (trim(msg))
         n = nmax
      end if
      
      ! Equation 1
      do i = 1,n
         do j = 1,n
            a1(j) = abs(x(i) - x(j))
         end do
         call indexx(n, a1, indx)
         ihimed = indx((n/2)+1)
         a2(i) = a1(ihimed)
         if (n .eq. 3) a2(i) = (a1(1)+a1(2)+a1(3))/3.
      end do
      
      call indexx(n, a2, indx)
      ilomed = indx((n+1)/2)
      sn = cn(n)*1.1926*a2(ilomed)
      
      return
      
   end subroutine croux
         
         
!*****************************************************************************************
   real function cn (n)
   
   ! Small sample correction terms for subroutine croux
         
      integer, intent(in) :: n
      
      if (n .eq. 2) then
         cn = 0.743
      else if (n .eq. 3) then
         ! cn = 1.851   
         cn = 1.172 ! Special correction for modified algorithm using average of differences
      else if (n .eq. 4) then
         cn = 0.954
      else if (n .eq. 5) then
         cn = 1.351
      else if (n .eq. 6) then
         cn = 0.993
      else if (n .eq. 7) then
         cn = 1.198
      else if (n .eq. 8) then
         cn = 1.005
      else if (n .eq. 9) then
         cn = 1.131
      else if (n .ge. 10) then
         if (mod(n,2) .eq. 0) then ! n even
            cn = 1.
         else ! n odd
            cn = real(n)/(real(n)-0.9)
         end if
      else
         cn = 1.
      end if
      
      return
      
   end function cn
   
   
!*****************************************************************************************
end module mloclib_statistics
