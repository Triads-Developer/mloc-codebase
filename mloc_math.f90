!> Math, statistical procedures, date and time, and geographical calculations

module mloc_math

   use mloc_declare

   implicit none
   save
      
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
   
      integer, intent(in) :: n
      integer :: i
      integer, dimension(n), intent(out) :: indx
      integer :: indxt
      integer :: ir
      integer :: j
      integer :: l
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
   
      integer, intent(in) :: n
      integer :: i
      integer, dimension(n), intent(out) :: indx
      integer :: indxt 
      integer :: ir
      integer :: j
      integer :: l
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
   
      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: indx
      integer, dimension(n), intent(out) :: irank
      integer :: j
   
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

         e_dsvd(i)=g
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
                  e_dsvd(j)=u(i,j)*hinv
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
                        u(j,k)=u(j,k)+scale*e_dsvd(k)
                     end do
                  end if
               end do
            end if
         else
            g=0.0d0
         end if

         y=dabs(q(i))+dabs(e_dsvd(i))
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
         g=e_dsvd(i)
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
            if (dabs(e_dsvd(l)) .le. eps) go to 6500
            lm1 = l-1
            if (dabs(q(lm1)) .le. eps) go to 6000
         end do

         ! cancellation of e_dsvd(l) if l .gt. 1

 6000    c=0.0d0
         scale=1.0d0
         l1=l-1
         do i=l,k
            f=scale*e_dsvd(i)
            e_dsvd(i)=c*e_dsvd(i)
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
            g=e_dsvd(k-1)
            h=e_dsvd(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
            g=sqrt(f*f+1.0d0)
            f=((x-z)*(x+z)+h*(y/(f+dsign(g,f))-h))/x

            ! next q-r transformation

            c=1.0d0
            scale=1.0d0
            lplus=l + 1
            do i=lplus,k
               g=e_dsvd(i)
               y=q(i)
               h=scale*g
               g=c*g
               z=dsqrt(f*f+h*h)
               e_dsvd(i-1)=z
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
            e_dsvd(l)=0.0d0
            e_dsvd(k)=f
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
!      real :: betai
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
!       real :: betacf
      real :: bt
      real :: eps
!       real :: gammln
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
   
      integer, parameter :: n_it_max = 200
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
      
      eps = 2.*epsilon(x)
   
      am = 1.
      bm = 1.
      az = 1.
      qab = a + b
      qap = a + 1.
      qam = a - 1.
      bz = 1. - qab*x/qap
      do m = 1,n_it_max
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
      call warnings ('betacf: a or b too big, or n_it_max too small; check the log file for details')
   
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
         
      integer, parameter :: n_croux = 1000 ! Maximum size of input array "x"
      integer :: i
      integer :: ihimed
      integer :: ilomed
      integer, dimension(n_croux) :: indx
      integer :: j
      integer :: n
      integer, intent(in) :: nin
      real, dimension(n_croux) :: a1
      real, dimension(n_croux) :: a2
      real, intent(out) :: sn
!      real :: cn
      real, dimension(nin), intent(in) :: x
      character(len=132) :: msg
      
      if (nin .le. n_croux .and. nin .ge. 2) then
         n = nin
      else if (nin .lt. 2) then
         write (msg,'(a,i6)') 'croux: illegal value for n: ', nin
         call warnings (trim(msg))
         sn = 1.0
         return   
      else if (nin .gt. n_croux) then
         write (msg,'(a,i6)') 'croux: nin exceeds maximum value, set to ', n_croux
         call fyi (trim(msg))
         n = n_croux
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
real function hms2s (i, j, atmp)

! Convert hours-minutes-seconds to seconds

   integer :: i
   integer :: j
   real :: atmp

   hms2s = real(3600*i + 60*j) + atmp

   return
   
end function hms2s
      

!*****************************************************************************************
subroutine timecr (ots0, hour, min, sec)

!  convert number of seconds since beginning of the day (ots0) to
!  hour-minute-seconds time. if ots0 < 0, hour < 0. if ots0 > 86400,
!  hour > 23.

   integer :: hour
   integer :: min
   real :: dsec
   real :: ots
   real :: ots0
   real :: sec
   
   ots = ots0

   dsec = ots
   hour = 0
   min = 0
   if (dsec .ge. 60.) then
      hour = 0
      min = 0
      do while (dsec .ge. 60.)
         dsec = dsec - 60.
         min = min + 1
         if (min .eq. 60) then
            min = 0
            hour = hour + 1
         end if
      end do
   else if (dsec .lt. 0.) then
      hour = -1
      min = 59
      do while (dsec .lt. -60.)
         dsec = dsec + 60.
         min = min - 1
         if (min .eq. -1) then
            min = 59
            hour = hour - 1
         end if
      end do
      dsec = 60. + dsec
   end if
 
   sec = dsec

   return
   
end subroutine timecr

!*****************************************************************************************
subroutine juldat (year, month, day, julday, iflg)

! Conversion between Gregorian and day-of-year (Julian day). The conversion direction
! is controlled by input parameter "iflg". Based on code by R. Buland, re-structured
! by eab.

! iflg = 0:
!   Given the usual year, month, and day (as integer numbers) subroutine
!   juldat returns integer Julian day JULDAY.  It properly accounts for
!   the leap years through 2099.

! iflg = 1:
!   Given the integer year and julian day, subroutine gredat returns the
!   correct integer month and day of the month.  Correctly treats
!   leap years through 2099.

   integer :: day
   integer, dimension(12) :: dpm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer :: feb = 28
   integer :: i
   integer, intent(in) :: iflg
   integer :: im
   integer :: julday
   integer :: month
   integer :: year
   character(len=132) :: message
   
   select case (iflg)
      case (0)
         dpm(2) = feb
         if ((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or. (mod(year,400) .eq. 0)) dpm(2) = feb + 1
         if (month .ge. 1 .and. month .le. 12 .and. day .ge. 1 .and. day .le. dpm(month)) then
            julday = day
            if (month .ge. 2) then
               im = month - 1
               do i = 1,im
                  julday = julday + dpm(i)
               end do
            end if
         else
            write (message,'(a,2i8)') 'juldat: illegal value for MONTH or DAY: ', month, day
            call oops (message)
         end if
      
      case (1)
         if (julday .ge. 1 .and. julday .le. (366 - min0(mod(year,4),1))) then
            dpm(2) = feb
            if (mod(year,4) .eq. 0) dpm(2) = feb + 1
            day = julday
            month = 1
            do i = 1,12
               if (day .le. dpm(i)) exit
               day = day - dpm(i)
               month=month + 1
            end do
         else
            write (message,'(a,i8)') 'juldat: illegal value for JULDAY: ', julday
            call oops (message)
         end if
      
      case default
         write (message,'(a,i8)') 'juldat: illegal value for IFLG: ', iflg
         call oops (message)
      
   end select

   return
   
end subroutine juldat
     
!*****************************************************************************************
logical function date_range (julday_test, julday_start, julday_end)

! Test if a Julian date is within a given date range.
! Date format is YYYYDDD

   integer :: julday_end
   integer :: julday_start
   integer :: julday_test
   character(len=132) :: message

   date_range = .false.

   if (julday_start .eq. 0 .and. julday_end .eq. 0) then ! No info on date range
      date_range = .true.
   else if (julday_start .eq. 0 .and. julday_end .gt. 0) then ! No info on start of date range
      if (julday_test .le. julday_end) date_range = .true.
   else if (julday_start .gt. 0 .and. julday_end .eq. 0) then ! No info on end of date range
      if (julday_test .ge. julday_start) date_range = .true.
   else if (julday_start .gt. 0 .and. julday_end .gt. 0) then ! Both ends of date range are defined
      if (julday_test .ge. julday_start .and. julday_test .le. julday_end) date_range = .true.
   else
      write (message,'(a,3i8)') 'date_range: illegal dates ', julday_test, julday_start, julday_end
      call oops (message)
   end if

   return
   
end function date_range
      
      
!*****************************************************************************************
subroutine unix_time (year, month, day, hr, minute, sec, eot)

! Unix time converted from 'calendar' or 'human' time
! Referenced to the Unix 'epoch': 00:00:00 UTC on 1 January 1970
! Times for dates before the epoch are negative.
! Fractional seconds are supported.
! This algorithm is valid for dates from 1/1/1902.
! Although the Unix epoch is January 1, 1970, the algorithm calculates number of days to
! the target date from January 1, 1900, and then subtracts the number of days between
! January 1, 1900 and January 1, 1970.

! Written by EAB (2016/03/10)

   integer, parameter :: real_kind = selected_real_kind(16,30)
   real(kind=real_kind) :: eot
   
   integer :: day
   integer :: day_days
   integer, dimension(12) :: dpm = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! days per month (non-leap-year)
   integer :: hr
   integer :: i
   integer :: minute
   integer :: month
   integer :: month_days
   integer, parameter :: offset_days = 25568 ! Number of days between 1/1/1900 and 1/1/1970
   integer, parameter :: spd = 86400 ! seconds per day
   integer :: whole_days
   integer :: whole_days_corrected
   integer :: year
   integer :: year_days
   real :: sec
   
   if (year .lt. 1902) then
      write (*,'(a)') 'unix_time: invalid date'
      eot = 0.0_real_kind
      return
   end if
   
   ! Number of days in whole years up to the year of the date of interest, starting from 1900.
   ! 1900 was not a leap-year.
   year_days = 0
   if (year .ge. 1901) then
      do i = 1900,year-1
         year_days = year_days + 365
         if (mod(i,4) .eq. 0 .and. year .ne. 1900) year_days = year_days + 1 ! leap-year
      end do
   end if
   
   ! Number of days in whole months of the target year, up to the date of interest
   month_days = 0
   if (month .ge. 2) then
      do i = 1,month-1
         if (mod(year,4) .eq. 0 .and. year .ne. 1900 .and. i .eq. 2) then ! February of a leap year
            month_days = month_days + dpm(i) + 1
         else
            month_days = month_days + dpm(i)
         end if
      end do
   end if
      
   day_days = day - 1 ! Number of days up to the current day of interest  
   whole_days = year_days + month_days + day_days ! Number of whole days since 1/1/1900
   whole_days_corrected = whole_days - offset_days ! Correct for the offset from 1/1/1900 to 1/1/1970
   
   ! Convert to seconds and add the time elapsed on the current day
   eot = real(whole_days_corrected*spd,real_kind) + real(hr*3600) + real(minute*60) + sec
   
   return
   
end subroutine unix_time


!*****************************************************************************************
real function dgkmlogc (atmp)

   real :: atmp

   dgkmlogc = dgkmla/sin(atmp*rpd) ! geocentric coordinates

   return
   
end function dgkmlogc

      
!*****************************************************************************************
real function dgkmlo (atmp)

   real :: atmp

   dgkmlo = dgkmla/cos(atmp*rpd) ! geographic coordinates

   return
   
end function dgkmlo


!*****************************************************************************************
subroutine azgap (iev, openaz)

! Largest open angle from a set of azimuths. Only uses azimuth
! values for which the corresponding value of fltr is negative. In single
! mode, the open azimuth is based on fltrh. For a cluster, it is based
! on fltrc.

   integer :: i
   integer :: iev
   integer, dimension(n_arrtim_max) :: indx
   integer :: j
   integer :: k
   integer :: naz
   real :: azdif
   real, dimension(n_arrtim_max) :: aztest
   real :: gap
   real :: openaz
   character(len=132) :: msg
   logical :: fltrtest
      
   naz = 0
   do i = 1, nst(iev) ! Form array of azimuths to test
      fltrtest = fltrc(iev,i)
      if (.not.fltrtest) then
         naz = naz + 1
         aztest(naz) = azes(iev,i)
         if (aztest(naz) .lt. 0. .or. aztest(naz) .ge. 360.) then
            if (aztest(naz) .lt. 0.) then
               do while (aztest(naz) .lt. 0.)
                  aztest(naz) = aztest(naz) + 360.
               end do
            else if (aztest(naz) .ge. 360.) then
               do while (aztest(naz) .ge. 360.)
                  aztest(naz) = aztest(naz) - 360.
               end do
            end if
         end if
      end if
   end do
   if (naz .le. 1) then
      openaz = 360.
      write (msg,'(a,i3,a)') 'azgap: naz = ', naz, '; openaz = 360'
      call warnings (trim(msg))
      return
   end if

   if (naz .ge. 1) then
      call indexx (naz, aztest, indx)
   else
      write (msg,'(a,i4,a)') 'azgap: illegal value for naz (', naz, ')'
      call oops (trim(msg))
   end if

   openaz = 0.
   do i = 1,naz-1
      j = indx(i)
      k = indx(i+1)
      azdif = aztest(k) - aztest(j)
      openaz = amax1(openaz,azdif)
   end do

   ! Check gap between largest and smallest azimuths
   gap = aztest(indx(1)) + 360. - aztest(indx(naz))
   openaz = amax1(openaz,gap)

   return
   
end subroutine azgap


!*****************************************************************************************
subroutine set_longitude_range (lon, longitude_range)

! Keep longitude in the desired range, defined by the value of input variable
! 'longitude_range', which gives the center longitude of the desired range.

   integer, intent(in) :: longitude_range
   real, intent(inout) :: lon
   real :: temp
   
   if (lon .lt. real(longitude_range-180)) then
      temp = lon
      lon = lon + 360.
      if (debug) write (io_log,'(a,f8.3,a,f8.3)') 'set_longitude_range: ', temp,&
       ' changed to ', lon
   end if
   
   if (lon .ge. real(longitude_range+180)) then
      temp = lon
      lon = lon - 360.
      if (debug) write (io_log,'(a,f8.3,a,f8.3)') 'set_longitude_range: ', temp,&
       ' changed to ', lon
   end if
         
   return
   
end subroutine set_longitude_range
 

!*****************************************************************************************
subroutine geocen (glat, glon, elat, elats, elatc, elon, elons, elonc)

!  Converts geographic to geocentric coordinates. Latitude is measured 
!  positive to the south from the north pole, longitude is positive to the
!  east from the prime meridian.
     
!  Input:
!    glat   latitude in geographic coordinates
!    glon   longitude in geographic coordinates

!  Output:
!    elat   geocentric latitude
!    elats  sin(geocentric latitude)
!    elatc  cos(geocentric latitude)
!    elon   geocentric longitude
!    elons  sin(geocentric longitude)
!    elonc  cos(geocentric longitude)

   real, intent(out) :: elat
   real, intent(out) :: elatc
   real, intent(out) :: elats
   real, intent(out) :: elon
   real, intent(out) :: elonc
   real, intent(out) :: elons
   real, intent(in) :: glat
   real, intent(in) :: glon
         
   elat = 90.0 - dpr*atan(0.9932773*(sin(glat*rpd)/cos(glat*rpd)))
   elats = sin(elat*rpd)
   elatc = cos(elat*rpd)
   elon = glon
   if (glon .lt. 0.0) elon = 360.0 + glon
   elons = sin(elon*rpd)
   elonc = cos(elon*rpd)
   
   return
   
end subroutine geocen
      
      
!*****************************************************************************************
subroutine geogra (elat, glat)
      
! Input: elat    geocentric latitude
! Output glat    geographic latitude
! Engdahl's code

   real, intent(in) :: elat
   real :: elatc
   real :: elats
   real, intent(out) :: glat

   elats = sin(elat*rpd)
   elatc = cos(elat*rpd)
   glat = dpr*atan(1.0067682*(elatc/elats))
   
   return
   
end subroutine geogra


!*****************************************************************************************
subroutine delaz (eqlt, eqln, stlt, stln, delta, deltdg, deltkm, azeqst, azesdg, azsteq, azsedg, i)

!  compute distance and azimuth from earthquake (eq) to station (st).

!  delta  = distance between (eq) and (st) in radians.
!  azeqst = azimuth from (eq) to (st) clockwise from north in radians
!  azsteq = azimuth from (st) to (eq) clockwise from north in radians
!  cos(delta)  = aa
!  sin(delta)  = bb
!  sin(azeqst) = cc
!  sin(azsteq) = dd

!  i=0 if input coordinates are geographic degrees.
!  i=1 if input coordinates are geocentric radians.

!  subroutines called:
!    coortr (appended)

   integer :: deltm ! epicentral distance in meters
   integer, intent(in) :: i
   real :: angle
   real, intent(out) :: azeqst
   real, intent(out) :: azesdg
   real, intent(out) :: azsedg
   real, intent(out) :: azsteq
   real :: cc
   real :: dd
   real, intent(out) :: delta
   real, intent(out) :: deltdg
   real, intent(out) :: deltkm
!   real :: dgkmlo
   real :: dx
   real :: dy
   real, intent(inout) :: eqln
   real :: eqlnrd
   real, intent(inout) :: eqlt
   real :: eqltrd
   real :: eqstln
   real, intent(inout) :: stln
   real :: stlnrd
   real, intent(inout) :: stlt
   real :: stltrd
   double precision :: aa
   double precision :: bb
   double precision :: dce
   double precision :: dces
   double precision :: dcs
   double precision :: deqlt
   double precision :: desln
   double precision :: dse
   double precision :: dses
   double precision :: dss
   double precision :: dstlt
   character(len=132) :: msg
   
   ! Convert to spherical polar coordinates in radians.
   if (i .eq. 0) then
      call coortr (eqltrd, eqlnrd, eqlt, eqln, i)
      call coortr (stltrd, stlnrd, stlt, stln, i)
      eqltrd = pi/2. - eqltrd
      stltrd = pi/2. - stltrd
      dy = (stlt - eqlt)/dgkmla ! latitude difference in km
      dx = (stln - eqln)/dgkmlo(eqlt) ! longitude difference in km
   else if (i .eq. 1) then
!         eqltrd=pi/2.-eqlt
      eqltrd = eqlt
      eqlnrd = eqln
!         stltrd=pi/2.-stlt
      stltrd = stlt
      stlnrd = stln
      dy = (stltrd - eqltrd)*radius ! latitude difference in km
      dx = (stlnrd - eqlnrd)*radius ! longitude difference in km
   else
      write (msg,'(a,i3)') 'delaz: illegal index: ', i
      call oops (trim(msg))
   end if
   deltkm = sqrt(dx*dx + dy*dy)
   
   if (deltkm .lt. 1.e-3) then ! Bail out if the separation is less than 1 meter
   
      azeqst = 0.
      azsteq = 0.
      azesdg = 0.
      azsedg = 0.
      delta = 0.
      deltdg = 0.
      deltkm = 0.
      if (verbose_screen) call fyi ('delaz: identical eq and st; delta and azimuth set to zero')
      return
      
   else if (deltkm .lt. 1.) then ! Do short distances (< 1 km) in plane geometry
   
      angle = atan2(dy,dx)
      azeqst = (pi*0.5) - angle ! Convert to azimuth, clockwise from north, still in radians
      if (azeqst .lt. 0.) azeqst = pi*2. + azeqst
      azsteq = azeqst - pi
      if (azsteq .lt. 0.) azsteq = pi*2. + azsteq
      azesdg = azeqst*dpr
      azsedg = azsteq*dpr
      deltdg = deltkm*dgkmla
      delta = deltdg*rpd
      deltm = nint(deltkm*1.0e3)
      if (verbose_screen) then
         write (msg,'(a,i4,a,2(f5.3,a))') 'delaz: separation is only ', deltm, ' meters (',&
          dx, ',', dy, ')'
         call fyi (trim(msg))
         write (*,'(t5,4f10.4)') eqlt, eqln, stlt, stln
      end if

   else ! Traditional calculation in spherical coordinates

      eqstln = stlnrd - eqlnrd
      desln = dble(eqstln)
      deqlt = dble(eqltrd)
      dstlt = dble(stltrd)
      dse = dsin(deqlt)
      dce = dcos(deqlt)
      dss = dsin(dstlt)
      dcs = dcos(dstlt)
      dses = dsin(desln)
      dces = dcos(desln)
      aa = dce*dcs + dse*dss*dces
      bb = dsqrt(1.0d0-aa*aa)
      cc = sngl(dss*dses/bb)
      dd = sngl(-dse*dses/bb)
      delta = sngl(datan2(bb,aa))
      azeqst = asin(cc)
      azsteq = asin(dd)
      if ((dse*dcs - dce*dss*dces) .lt. 0.) azeqst = pi - azeqst
      if (azeqst .lt. 0.) azeqst = pi*2. + azeqst
      if ((dce*dss - dse*dcs*dces) .le. 0.) azsteq = pi - azsteq
      if (azsteq .lt. 0.) azsteq = pi*2. + azsteq
      deltdg = delta*dpr
      azesdg = azeqst*dpr
      azsedg = azsteq*dpr
      deltkm = delta*radius
   
   end if
   
   return
   
end subroutine delaz


!*****************************************************************************************
subroutine coortr (alatrd, alonrd, alatdg, alondg, i)
      
!  alatrd (geocentric latitude in radians)
!  alonrd (geocentric longitude in radians)
!  alatdg (geographical latitude in degrees)
!  alogdg (geographical latitude in degrees)

!  if i=1, transformation from radians to degrees
!  if i=0, transformation from degrees to radians

!  no subroutines called

   integer, intent(in) :: i
   real :: aaa
   real :: alat2
   real, intent(inout) :: alatdg
   real, intent(inout) :: alatrd
   real, intent(inout) :: alondg
   real, intent(inout) :: alonrd
   real :: bbb
   character(len=132) :: msg

   if (i .eq. 0) then ! Degrees to radians
      alatrd = alatdg*rpd
      alonrd = alondg*rpd
      bbb = abs(alatdg)
      if (bbb .lt. 89.9) then
         aaa = 0.9933056*tan(alatrd)
         alatrd = atan(aaa)
      end if
   else if (i .eq. 1) then ! Radians to degrees
      bbb = abs(alatrd)
      if (bbb .lt. 1.57) then
         aaa = tan(alatrd)/0.9933056
         alat2 = atan(aaa)
      else 
         alat2 = alatrd
      end if
      alatdg = alat2*dpr
      alondg = alonrd*dpr
   else
      write (msg,'(a,i9)') 'coortr: illegal value for i: ', i
      call oops (trim(msg))
   end if

   return
   
end subroutine coortr


!*****************************************************************************************
real function dlttrn (depth)

!  for a ray with a turning depth of 'depth', 'dlttrn' is the epicentral
!  distance from the turning point to where the ray reaches the earths
!  surface (degrees). it is calculated with respect to the 1968 herrin
!  tables for depths less than 800 km.

   integer :: i
   integer :: j
   real, dimension(17) :: d = (/0.0,3.0,6.5,7.6,8.3,8.7,9.2,9.7,10.1,10.4,10.7,11.2,11.8,12.2,13.2,14.2,16.5/)
   real, intent(in) :: depth
   real :: x
   
   x = depth/50. + 1.
   i = int(x)
   j = i + 1
   x = x - i
   dlttrn = d(i) + x*(d(j)-d(i))

   return
   
end function dlttrn


!*****************************************************************************************
end module mloc_math

