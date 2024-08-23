!> Procedures related to array manipulation and calculation of confidence ellipses, mainly
! used in the hypocentroidal decomposition

module mloclib_inv

   use declare_constants
   use declare_hdc
   
   implicit none
   save

   double precision, allocatable, dimension(:) :: e
   
contains


!*****************************************************************************************      
   subroutine ell2cv (kcrit, alpha, sminor, smajor, t11, t12, t22)
      
      ! Given a confidence ellipse and a critical value, returns the 
      ! equivalent 2x2 covariance matrix
      
      ! Input
      ! kcrit   critical value
      ! alpha   strike of semi-minor axis (positive clockwise from north)
      ! sminor  semi-minor axis length
      ! smajor  semi-major axis length
      
      ! Output
      ! t11     CV(1,1)
      ! t12     CV(1,2) and CV(2,1)
      ! t22     CV(2,2)
      
      real, intent(in) :: alpha
      real :: cf2
      real :: csf
      real :: fi
      real, intent(in) :: kcrit
      real :: lambda1
      real :: lambda2
      real :: sf2
      real, intent(in) :: smajor
      real, intent(in) :: sminor
      real, intent(out) :: t11
      real, intent(out) :: t12
      real, intent(out) :: t22
      
      lambda1 = (sminor*sminor)/kcrit
      lambda2 = (smajor*smajor)/kcrit
      
      fi = rpd*alpha
      
      cf2 = cos(fi)*cos(fi)
      sf2 = sin(fi)*sin(fi)
      csf = cos(fi)*sin(fi)
      
      t11 = lambda1*cf2 + lambda2*sf2
      t22 = lambda1*sf2 + lambda2*cf2
      t12 = (lambda1-lambda2)*csf

      return
      
   end subroutine ell2cv


!*****************************************************************************************
   subroutine elips (t22a, alphad, al, bl)
   
   ! Confidence ellipse
   
      real, intent(out) :: al
      real :: alpha
      real, intent(out) :: alphad
      real :: arg1
      real :: arg2
      real, intent(out) :: bl
      real :: c2a
      real :: s2a
      real :: sca
      real, dimension(2,2), intent(in) :: t22a
      real, dimension(2,2) :: t22b
   
      call smatinv2 (t22a, t22b, 2, 2, 2, 2)
      arg1 = 2.*t22b(1,2)
      arg2 = t22b(1,1) - t22b(2,2)
      alpha = 0.5*atan2(arg1,arg2)
      alphad = alpha/rpd
      c2a = cos(alpha)*cos(alpha)
      s2a = sin(alpha)*sin(alpha)
      sca = sin(alpha)*cos(alpha)
      al = t22b(1,1)*c2a + 2.*t22b(1,2)*sca + t22b(2,2)*s2a
      bl = t22b(1,1)*s2a - 2.*t22b(1,2)*sca + t22b(2,2)*c2a 
   
      return
      
   end subroutine elips

      
!*****************************************************************************************
   subroutine delips (t22a, alphad, al, bl)
   
   ! Confidence ellipse, using double precision
   
      real, intent(out) :: al
      real :: alpha
      real, intent(out) :: alphad
      real :: arg1
      real :: arg2
      real, intent(out) :: bl
      double precision :: c2a
      double precision :: s2a
      double precision :: sca
      double precision, dimension(2,2), intent(in) :: t22a
      double precision, dimension(2,2) :: t22b
   
      call dmatinv2 (t22a, t22b, 2, 2, 2, 2)
      arg1 = 2.*sngl(t22b(1,2))
      arg2 = sngl(t22b(1,1) - t22b(2,2))
      alpha = 0.5*atan2(arg1,arg2)
      alphad = alpha/rpd
      c2a = dble(cos(alpha)*cos(alpha))
      s2a = dble(sin(alpha)*sin(alpha))
      sca = dble(sin(alpha)*cos(alpha))
      al = sngl(t22b(1,1)*c2a + 2.0d0*t22b(1,2)*sca + t22b(2,2)*s2a)
      bl = sngl(t22b(1,1)*s2a - 2.0d0*t22b(1,2)*sca + t22b(2,2)*c2a) 
   
      return
      
   end subroutine delips


!*****************************************************************************************
   subroutine postq (ain, na, ma, n, m, bout)
   
   ! post-multiplies qbhat by matrix ain, without forming qbhat.
   ! qbhat * ain = bout
   ! dimension of ain and bout in calling program is (na,ma).
   ! dimensions used in multiplication are (n,m)
   ! Uses double precision.
   
      integer :: i
      integer :: ii
      integer :: j
      integer :: jbi
      integer :: k
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: n
      integer, intent(in) :: na
      double precision, dimension(na,ma), intent(in) :: ain
      double precision, dimension(na,ma), intent(out) :: bout
      double precision :: boutij 
      double precision :: qbhatik
      double precision :: tempij
   
      do j = 1,m
         do i = 1,n
            jbi = jb(i)
            tempij = dble(wq2(jbi)*sighatj(i))
            boutij = 0.0d0
            do k = 1,ntqi(jbi)
               ii = ntq(k,jbi)
               qbhatik = -1.0d0/(tempij*dble(sighatj(ii)))
               if (ii .eq. i) qbhatik = qbhatik + 1.0d0
               boutij = boutij + (qbhatik*ain(ii,j))
            end do
            bout(i,j) = boutij
         end do
      end do
   
      return
      
   end subroutine postq


!*****************************************************************************************
   subroutine postqv (ain, na, n, bout)
   
   ! post-multiplies qbhat by vector ain, without forming qbhat.
   ! qbhat * ain = bout
   ! dimension of ain and bout in calling program is na.
   ! vector length used in multiplication is n.
   ! Uses double precision.
   
      integer :: i
      integer :: ii
      integer :: jbi
      integer :: k
      integer, intent(in) :: n
      integer, intent(in) :: na
      double precision, dimension(na), intent(in) :: ain
      double precision, dimension(na), intent(out) :: bout
      double precision :: tempi
      double precision :: bouti
      double precision :: qbhatik
   
      do i = 1,n
         jbi = jb(i)
         tempi = dble(wq2(jbi)*sighatj(i))
         bouti = 0.0d0
         do k = 1,ntqi(jbi)
            ii = ntq(k,jbi)
            qbhatik = -1.0d0/(tempi*dble(sighatj(ii)))
            if (ii .eq. i) qbhatik = qbhatik + 1.0d0
            bouti = bouti + (qbhatik*ain(ii))
         end do
         bout(i) = bouti
      end do
   
      return
      
   end subroutine postqv


!*****************************************************************************************
   subroutine postpv (ain, na, n, bout)
   
   ! post-multiplies pbhat by vector ain, without forming qbhat.
   ! pbhat * ain = bout
   ! dimension of ain and bout in calling program is (na).
   ! vector length used in multiplication is n.
   ! Uses double precision.
   
      integer :: i
      integer :: ii
      integer :: jbi
      integer :: k
      integer, intent(in) :: n
      integer, intent(in) :: na
      double precision, dimension(na), intent(in) :: ain
      double precision, dimension(na), intent(out) :: bout
      double precision :: bouti
      double precision :: pbhatik
      double precision :: tempi
   
      do i = 1,n
         jbi = jb(i)
         tempi = wq2(jbi)*dble(sighatj(i))
         bouti = 0.0d0
         do k = 1,ntqi(jbi)
            ii = ntq(k,jbi)
            pbhatik = 1.0d0/(tempi*dble(sighatj(ii)))
            bouti = bouti + (pbhatik*ain(ii))
         end do
         bout(i) = bouti
      end do
   
      return
      
   end subroutine postpv

!*****************************************************************************************
   subroutine indexx (n, arrin, indx)
   
   ! indexes array 'arrin' of length 'n', i.e., outputs the array 'indx'
   ! such that arrin(indx(j)) is in ascending order for j=1,2,...,n.  the
   ! input quantities n and arrin are not changed.
   
      integer :: i
      integer, dimension(n), intent(out) :: indx
      integer :: indxt
      integer :: ir
      integer :: j
      integer :: l
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: arrin
      real :: q
      
      if (n .eq. 0) call oops ('indexx: n = 0')
      
      do j = 1,n
         indx(j) = j
      end do
      if (n .eq. 1) return
      l = n/2+1
      ir = n
   10 continue
      if (l .gt. 1) then
         l = l - 1
         indxt = indx(l)
         q = arrin(indxt)
      else
         indxt = indx(ir)
         q = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir-1
         if (ir .eq. 1) then
            indx(1) = indxt
            return
         end if
      end if
      i = l
      j = l + l
   20 if (j .le. ir) then
         if (j .lt. ir) then
            if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
         end if
         if (q .lt. arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
         else
            j = ir + 1
         end if
         go to 20 
      end if 
      indx(i) = indxt
      go to 10
   
   end subroutine indexx

         
!*****************************************************************************************
   subroutine dindexx (n, arrin, indx)
   
   ! indexes double precision array 'arrin' of length 'n', i.e., outputs 
   ! the array 'indx'
   ! such that arrin(indx(j)) is in ascending order for j=1,2,...,n.  the
   ! input quantities n and arrin are not changed.
   
      integer :: i
      integer, dimension(n), intent(out) :: indx
      integer :: indxt 
      integer :: ir
      integer :: j
      integer :: l
      integer, intent(in) :: n
      double precision, dimension(n), intent(in) :: arrin
      double precision :: q
   
      do j = 1,n
         indx(j) = j
      end do
      if (n .eq. 1) return
      l = n/2+1
      ir = n
   10 continue
      if (l .gt. 1) then
         l = l - 1
         indxt = indx(l)
         q = arrin(indxt)
      else
         indxt = indx(ir)
         q = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir .eq. 1) then
            indx(1) = indxt
            return
         end if
      end if
      i = l
      j = l + l
   20 if (j .le. ir) then
         if (j .lt. ir) then
            if (arrin(indx(j)) .lt. arrin(indx(j+1))) j = j + 1
         end if
         if (q .lt. arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
         else
            j = ir + 1
         end if
         go to 20 
      end if 
      indx(i) = indxt
      go to 10
   
   end subroutine dindexx
         
!*****************************************************************************************
   subroutine rank (n, indx, irank)
   
   ! given 'indx' of length 'n' as output from subroutine 'indexx', this
   ! subroutine returns an array 'irank', the corresponding table of ranks.
   
      integer, dimension(n), intent(in) :: indx
      integer, dimension(n), intent(out) :: irank
      integer :: j
      integer, intent(in) :: n
   
      do j = 1,n
         irank(indx(j)) = j
      end do
      
      return
      
   end subroutine rank

!*****************************************************************************************
   subroutine dot1 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)
   
   ! dot product of two matrices, single precision.
   ! a * b = c
   ! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
   ! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
   ! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).
   
   ! no subroutines called.
   
      integer :: i
      integer :: j
      integer :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: mb
      integer, intent(in) :: mc
      integer, intent(in) :: n
      integer, intent(in) :: na
      integer, intent(in) :: nb
      integer, intent(in) :: nc
      real, dimension(ma,na), intent(in) :: a
      real, dimension(mb,nb), intent(in) :: b
      real :: bjk
      real, dimension(mc,nc), intent(out) :: c
   
      do k = 1,n
         do i = 1,l
            c(i,k) = 0.
         end do
      end do 
   
      do k = 1,n
         do j =1,m
            bjk = b(j,k)
            do i = 1,l
               c(i,k) = c(i,k) + a(i,j)*bjk
            end do
         end do
      end do
   
      return
      
   end subroutine dot1


!*****************************************************************************************
   subroutine dot2 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)
   
   ! dot product of two matrices, single precision.
   ! a**tr * b = c
   ! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
   ! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
   ! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).
   
   ! no subroutines called.
   
      integer :: i
      integer :: j
      integer :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: mb
      integer, intent(in) :: mc
      integer, intent(in) :: n
      integer, intent(in) :: na
      integer, intent(in) :: nb
      integer, intent(in) :: nc
      real, dimension(ma,na), intent(in) :: a
      real, dimension(mb,nb), intent(in) :: b
      real, dimension(mc,nc), intent(out) :: c
      real :: cik
   
      do k = 1,n
         do i = 1,l
            c(i,k) = 0.
         end do
      end do 
   
      do k = 1,n
         do i = 1,l
            cik = c(i,k)
            do j = 1,m
               cik = cik + a(j,i)*b(j,k)
            end do
            c(i,k) = cik
         end do
      end do
   
      return
      
   end subroutine dot2


!*****************************************************************************************
   subroutine dot3 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)
   
   ! dot product of two matrices, single precision.
   ! a * b**tr = c
   ! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
   ! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
   ! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).
   
   ! no subroutines called.
   
      integer :: i
      integer :: j
      integer :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: mb
      integer, intent(in) :: mc
      integer, intent(in) :: n
      integer, intent(in) :: na
      integer, intent(in) :: nb
      integer, intent(in) :: nc
      real, dimension(ma,na), intent(in) :: a
      real, dimension(mb,nb), intent(in) :: b
      real :: bkj
      real, dimension(mc,nc), intent(out) :: c
   
      do k = 1,n
         do i = 1,l
            c(i,k) = 0.
         end do
      end do 
   
      do j = 1,m
         do k = 1,n
            bkj = b(k,j)
            do i = 1,l
               c(i,k) = c(i,k) + a(i,j)*bkj
            end do
         end do
      end do
   
      return
      
   end subroutine dot3
      
      
!*****************************************************************************************
   subroutine ddot1 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)
   
   ! dot product of two matrices, double precision.
   ! a * b = c
   ! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
   ! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
   ! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).
   
   ! no subroutines called.
   
      integer :: i
      integer :: j
      integer :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: mb
      integer, intent(in) :: mc
      integer, intent(in) :: n
      integer, intent(in) :: na
      integer, intent(in) :: nb
      integer, intent(in) :: nc
      double precision, dimension(ma,na), intent(in) :: a
      double precision, dimension(mb,nb), intent(in) :: b
      double precision :: bjk
      double precision, dimension(mc,nc), intent(out) :: c
   
      do k = 1,n
         do i = 1,l
            c(i,k) = 0.0d0
         end do
      end do 
   
      do k = 1,n
         do j = 1,m
            bjk = b(j,k)
            do i = 1,l
               c(i,k) = c(i,k) + a(i,j)*bjk
            end do
         end do
      end do
   
      return
      
   end subroutine ddot1


!***********************************************************************************************************************************
   subroutine ddot2 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)
   
   ! dot product of two matrices, double precision.
   ! a**tr * b = c
   ! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
   ! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
   ! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).
   
   ! no subroutines called.
   
      integer :: i
      integer :: j
      integer :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: mb
      integer, intent(in) :: mc
      integer, intent(in) :: n
      integer, intent(in) :: na
      integer, intent(in) :: nb
      integer, intent(in) :: nc
      double precision, dimension(ma,na), intent(in) :: a
      double precision, dimension(mb,nb), intent(in) :: b
      double precision, dimension(mc,nc), intent(out) :: c
      double precision :: cik
   
      do k = 1,n
         do i = 1,l
            c(i,k) = 0.0d0
         end do
      end do 
   
      do k = 1,n
         do i = 1,l
            cik = c(i,k)
            do j = 1,m
               cik = cik + a(j,i)*b(j,k)
            end do
            c(i,k) = cik
         end do
      end do
   
      return
      
   end subroutine ddot2


!***********************************************************************************************************************************
   subroutine ddot3 (a, b, c, ma, na, mb, nb, mc, nc, l, m, n)
   
   ! dot product of two matrices, double precision.
   ! a * b**tr = c
   ! dimensions of arrays in calling program are: a(ma,na), b(mb,nb), and
   ! c(mc,nc). for example, a vector b is dimensioned b(mb,1). the sizes of
   ! the arrays to be multiplied are: a(l,m)*b(m,n)=c(l,n).
   
   ! no subroutines called.
   
      integer :: i
      integer :: j
      integer :: k
      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(in) :: ma
      integer, intent(in) :: mb
      integer, intent(in) :: mc
      integer, intent(in) :: n
      integer, intent(in) :: na
      integer, intent(in) :: nb
      integer, intent(in) :: nc
      double precision, dimension(ma,na), intent(in) :: a
      double precision, dimension(mb,nb), intent(in) :: b
      double precision :: bkj
      double precision, dimension(mc,nc), intent(out) :: c
   
      do k = 1,n
         do i = 1,l
            c(i,k) = 0.0d0
         end do
      end do 
   
      do j = 1,m
         do k = 1,n
            bkj = b(k,j)
            do i = 1,l
               c(i,k) = c(i,k) + a(i,j)*bkj
            end do
         end do
      end do
   
      return
      
   end subroutine ddot3


!*****************************************************************************************
   subroutine smatinv2 (a, b, idim, jdim, ia, ja)
   
   ! inverse of square matrix a returned in b.
   ! single precision.  a is unchanged.
   
      integer :: i
      integer, intent(in) :: ia
      integer, intent(in) :: idim
      integer, dimension(500) :: ik
      integer :: j
      integer, intent(in) :: ja
      integer, intent(in) :: jdim
      integer, dimension(500) :: jk
      integer :: k
      integer :: l
      real, dimension(idim,jdim), intent(in) :: a
      real, dimension(idim,jdim), intent(out) :: b
      real :: amax
      real :: save_me
      real :: eps
      character(len=132) :: msg
      
      eps = 2.*epsilon(amax)
   
      if (ia .ne. ja) then
        call warnings ('smatinv2: non-square matrix dimensions')
        return
      end if
   
      do j = 1,ja
         do i = 1,ia
            b(i,j) = a(i,j)
         end do
      end do
   
      do k = 1,ia
   
         ! find largest element b(i,j) in rest of matrix.
         amax = 0.0
         do j = k,ia
            do i = k,ia
               if (abs(amax) .lt. abs(b(i,j))) then
                  amax = b(i,j)
                  ik(k) = i
                  jk(k) = j
               end if
            end do
         end do
   
         ! interchange rows and columns of matrix to put amax in b(k,k)
         ! (ie. do partial pivoting to improve accuracy).
         if (abs(amax) .lt. eps) then
           write (msg,'(a,10e12.3)') 'smatinv2: matrix is singular; ', a
           call oops (trim(msg))
         end if
   
         i = ik(k)
         if (i .ne. k) then
            do j = 1,ia
               save_me = b(k,j)
               b(k,j) = b(i,j)
               b(i,j) = -save_me
            end do
         end if
   
         j = jk(k)
         if (j .ne. k) then
            do i = 1,ia
               save_me = b(i,k)
               b(i,k) = b(i,j)
               b(i,j) = -save_me
            end do
         end if
   
         ! accumulate elements of the inverse matrix.
         save_me = b(k,k)
         do i = 1,ia
            b(i,k) = -b(i,k)/amax
         end do
         b(k,k) = save_me
   
         do j = 1,ia
            if (k .ne. j) then
               save_me = b(k,j)
               do i = 1,ia
                  b(i,j) = b(i,j) + b(i,k)*save_me
               end do
               b(k,j) = save_me
            end if
         end do
         do j = 1,ia
            b(k,j) = b(k,j)/amax
         end do
         b(k,k) = 1.0/amax
      end do 
   
      ! restore ordering of the matrix.
      do l = 1,ia
         k = ia - l + 1
         j = ik(k)
         if (j .gt. k) then
            do i = 1,ia
               save_me = b(i,k)
               b(i,k) = -b(i,j)
               b(i,j) = save_me
            end do
         end if
         i = jk(k)
         if (i .gt. k) then
            do j = 1,ia
               save_me = b(k,j)
               b(k,j) = -b(i,j)
               b(i,j) = save_me
            end do
         end if
      end do
   
      return
      
   end subroutine smatinv2
      
      
!*****************************************************************************************
   subroutine dmatinv2 (a, b, idim, jdim, ia, ja)
   
   ! inverse of square matrix a returned in b.
   ! double precision.  a is unchanged.
   
      integer :: i
      integer, intent(in) :: ia
      integer, intent(in) :: idim
      integer, dimension(500) :: ik
      integer :: j
      integer, intent(in) :: ja
      integer, intent(in) :: jdim
      integer, dimension(500) :: jk
      integer :: k
      integer :: l
      double precision, dimension(idim,jdim), intent(in) :: a
      double precision, dimension(idim,jdim), intent(out) :: b
      double precision :: amax
      double precision :: save_this
      double precision :: eps

      eps = 2.*epsilon(amax)
   
      if (ia .ne. ja) then
        call warnings ('dmatinv2: non-square matrix dimensions')
        return
      end if
   
      do j = 1,ja
         do i = 1,ia
            b(i,j) = a(i,j)
         end do
      end do
   
      do k = 1,ia
   
         ! find largest element b(i,j) in rest of matrix.
         amax = 0.0d0
         do j = k,ia
            do i = k,ia
               if (abs(amax) .lt. abs(b(i,j))) then
                  amax = b(i,j)
                  ik(k) = i
                  jk(k) = j
               end if
            end do
         end do
   
         ! interchange rows and columns of matrix to put amax in b(k,k)
         ! (ie. do partial pivoting to improve accuracy).
         if (abs(amax) .lt. eps) then
           call warnings ('dmatinv2: matrix is singular')
           return
         end if
   
         i = ik(k)
         if (i .ne. k) then
            do j = 1,ia
               save_this = b(k,j)
               b(k,j) = b(i,j)
               b(i,j) = -save_this
            end do
         end if
   
         j = jk(k)
         if (j .ne. k) then
            do i = 1,ia
               save_this = b(i,k)
               b(i,k) = b(i,j)
               b(i,j) = -save_this
            end do
         end if
   
         ! accumulate elements of the inverse matrix.
         save_this = b(k,k)
         do i = 1,ia
            b(i,k) = -b(i,k)/amax
         end do
         b(k,k) = save_this
   
         do j = 1,ia
            if (k .ne. j) then
               save_this = b(k,j)
               do i = 1,ia
                  b(i,j) = b(i,j) + b(i,k)*save_this
               end do
               b(k,j) = save_this
            end if
         end do
         do j = 1,ia
            b(k,j) = b(k,j)/amax
         end do
         b(k,k) = 1.0/amax
      end do 
   
      ! restore ordering of the matrix.
      do l = 1,ia
         k = ia - l + 1
         j = ik(k)
         if (j .gt. k) then
            do i = 1,ia
               save_this = b(i,k)
               b(i,k) = -b(i,j)
               b(i,j) = save_this
            end do
         end if
         i = jk(k)
         if (i .gt. k) then
            do j = 1,ia
               save_this = b(k,j)
               b(k,j) = -b(i,j)
               b(i,j) = save_this
            end do
         end if
      end do
   
      return
      
   end subroutine dmatinv2


!*****************************************************************************************
   subroutine dsvd (ain, m, n, mp, np, u, v, q)

   ! singular value decomposition.
   
   ! further attempts made to optimize 7/8/88 eab
   
   ! double precision
   
   ! for algol program see wilkinson+reinsch: handbook for automatic
   ! computation vol 2 - linear algebra, pp. 140-144. translated from
   ! algol by r.l.parker.
   
   ! cleaned up by e. bergman - 
   
   ! the matrix ain(m,n) is decomposed. singular values in q, pre-matrix in u,
   ! post-matrix in v. 
   
   ! program altered by p. silver to handle unpacked arrays.
   ! mp and np are dimensions in main routine. m and n are actual
   ! dimensions to be used in the subroutine.
   
   ! no subroutines called.
   
   ! eps = machine accuracy: smallest double precision value which,
   !       when added to 1.0, yields a result different from 1.0.
   ! tol = smallest double precision floating point number considered to be
   !       greater than 0.
   
   ! the include file is only used to bring in the dimension of the
   ! error vector 'e' from the main program.

      integer :: i
      integer :: i1d
      integer :: iback
      integer :: j
      integer :: k
      integer :: kback
      integer :: l
      integer :: l1
      integer :: lback
      integer :: lm1
      integer :: lplus
      integer, intent(in) :: m
      integer, intent(in) :: mp
      integer, intent(in) :: n
      integer, intent(in) :: np
      double precision, intent(in) :: ain(mp,np)
      double precision :: c 
      double precision :: eps
      double precision :: f
      double precision :: g
      double precision :: ginv
      double precision :: h
      double precision :: hinv
      double precision, intent(out) :: q(np)
      double precision :: scale
      double precision :: tol
      double precision, intent(out) :: u(mp,np)
      double precision, intent(out) :: v(np,np)
      double precision :: x
      double precision :: y
      double precision :: z

      eps=1.0d-15
      tol=1.0d-100

      do j=1,n
         do i=1,m
            u(i,j)=ain(i,j)
         end do
      end do

      ! householder reduction to bi-diagonal form

      g=0.0d0
      x=0.0d0
      do i=1,n

         e(i)=g
         scale=0.0d0
         l=i+1
         do j=i,m
            scale=u(j,i)*u(j,i)+scale
         end do

         if (scale .ge. tol) then
            f=u(i,i)
            g=-dsign(dsqrt(scale),f)
            h=f*g-scale
            u(i,i)=f-g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  scale=0.0d0
                  do k=i,m
                     scale=u(k,i)*u(k,j)+scale
                  end do
                  f=scale*hinv
                  do k=i,m
                     u(k,j)=u(k,j)+f*u(k,i)
                  end do 
               end do
            end if
         else
            g=0.0d0
         end if

         q(i)=g
         scale=0.0d0
         if (l .le. n) then
            do j=l,n
               scale=u(i,j)*u(i,j)+scale
            end do
         end if

         if (scale .ge. tol) then
            i1d = i+1
            f=u(i,i1d)
            g=-dsign(dsqrt(scale),f)
            h=f*g-scale
            u(i,i+1)=f-g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  e(j)=u(i,j)*hinv
               end do
            end if
            if (l .le. m) then
               do j=l,m
                  scale=0.0d0
                  if (l .le. n) then
                     do k=l,n
                        scale=u(j,k)*u(i,k)+scale
                     end do
                     do k=l,n
                        u(j,k)=u(j,k)+scale*e(k)
                     end do
                  end if
               end do
            end if
         else
            g=0.0d0
         end if

         y=dabs(q(i))+dabs(e(i))
         if (y .gt. x) x=y

      end do

      ! accumulation of right-hand transforms (v)

      do iback=1,n
         i=n+1-iback
         if (abs(g) .gt. eps) then
            h=u(i,i+1)*g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  v(j,i)=u(i,j)*hinv
               end do
               do j=l,n
                  scale=0.0d0
                  do k=l,n
                     scale=u(i,k)*v(k,j)+scale
                  end do
                  do k=l,n
                     v(k,j)=v(k,j)+scale*v(k,i)
                  end do
               end do
            end if
         end if
         if (l .le. n) then
            do j=l,n
               v(j,i)=0.0d0
               v(i,j)=0.0d0
            end do
         end if
         v(i,i)=1.0d0
         g=e(i)
         l=i
      end do

      ! accumulation of left-hand transforms

      do iback=1,n
         i=n+1-iback
         l=i+1
         g=q(i)
         if (l .le. n) then
            do j=l,n
               u(i,j)=0.0d0
            end do
         end if
         if (abs(g) .gt. eps) then
            h=u(i,i)*g
            if (l .le. n) then
               hinv=1.0d0/h
               do j=l,n
                  scale=0.0d0
                  do k=l,m
                     scale=u(k,i)*u(k,j)+scale
                  end do
                  f=scale*hinv
                  do k=i,m
                     u(k,j)=u(k,j)+f*u(k,i)
                  end do
               end do
            end if
            ginv=1.0d0/g
            do j=i,m
               u(j,i)=u(j,i)*ginv
            end do
         else
            do j=i,m
               u(j,i)=0.0d0
            end do
         end if
         u(i,i)=u(i,i)+1.0d0
      end do

      ! diagonalization of bi-diagonal form

      eps=eps*x
      do kback=1,n
         k=n+1-kback

         ! test f-splitting

 5000    continue
         do lback=1,k
            l=k+1-lback
            if (dabs(e(l)) .le. eps) go to 6500
            lm1 = l-1
            if (dabs(q(lm1)) .le. eps) go to 6000
         end do

         ! cancellation of e(l) if l .gt. 1

 6000    c=0.0d0
         scale=1.0d0
         l1=l-1
         do i=l,k
            f=scale*e(i)
            e(i)=c*e(i)
            if (dabs(f) .le. eps) go to 6500
            g=q(i)
            q(i)=dsqrt(f*f+g*g)
            h=q(i)
            c=g/h
            scale=-f/h
            do j=1,m
               y=u(j,l1)
               z=u(j,i)
               u(j,l1)=y*c+z*scale
               u(j,i)=-y*scale+z*c
            end do
         end do

         ! test f-convergence

 6500    z=q(k)
         if (l .ne. k) then
            x=q(l)
            y=q(k-1)
            g=e(k-1)
            h=e(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
            g=sqrt(f*f+1.0d0)
            f=((x-z)*(x+z)+h*(y/(f+dsign(g,f))-h))/x

            ! next q-r transformation

            c=1.0d0
            scale=1.0d0
            lplus=l + 1
            do i=lplus,k
               g=e(i)
               y=q(i)
               h=scale*g
               g=c*g
               z=dsqrt(f*f+h*h)
               e(i-1)=z
               c=f/z
               scale=h/z
               f=x*c+g*scale
               g=-x*scale+g*c
               h=y*scale
               y=y*c
               do j=1,n
                  x=v(j,i-1)
                  z=v(j,i)
                  v(j,i-1)=x*c+z*scale
                  v(j,i)=-x*scale+z*c
               end do
               z=dsqrt(f*f+h*h)
               q(i-1)=z
               c=f/z
               scale=h/z
               f=c*g+scale*y
               x=-scale*g+c*y
               do j=1,m
                  y=u(j,i-1)
                  z=u(j,i)
                  u(j,i-1)=y*c+z*scale
                  u(j,i)=-y*scale+z*c
               end do
            end do
            e(l)=0.0d0
            e(k)=f
            q(k)=x
            go to 5000
         end if

         ! convergence

         if (z .lt. 0.0d0) then
            q(k)=-z
            do j=1,n
               v(j,k)=-v(j,k)
            end do
         end if

      end do

      return
      
   end subroutine dsvd
      

!*****************************************************************************************
   subroutine dsvd_allocate ()
   
   ! Allocate variable arrays related to direct calibration
   
      integer :: error
      integer :: m
      character(len=32) :: p
      
      m = mtmax
      p = 'dsvd_allocate' ! procedure
      
      allocate (e(m), stat=error); if (error .gt. 0) call allocation_error (p, 'e', error)
      
      return
   
   end subroutine dsvd_allocate   
   

!*****************************************************************************************
end module mloclib_inv

