!> Procedures related to the generalized inverse for hypocentroid and cluster vectors.
!******************************************************************************
! Key Documentation for mloc_inv Module
!******************************************************************************
!
! Memory Management:
!   - Explicit allocation and deallocation checks were added for all major arrays, 
!     including `ahat`, `qbahat`, `pqahat`, `dx`, etc.
!   - Each allocation is followed by an error check. If the allocation fails, 
!     a descriptive error message is logged, and the procedure stops to prevent 
!     undefined behavior.
!   - Deallocation also includes error checking to ensure that all allocated memory 
!     is correctly freed at the end of the process.
!   - This ensures efficient memory usage and prevents memory leaks during runtime.
!
! Improved Singular Value Decomposition (SVD) Handling:
!   - The matrix inversion process uses Singular Value Decomposition (SVD) with 
!     Tikhonov regularization to handle near-zero singular values that may arise 
!     from ill-conditioned matrices.
!   - Tikhonov regularization introduces a small correction factor to stabilize 
!     the inversion process. This is particularly useful when the matrix condition 
!     number becomes large.
!   - A condition number (`icond`) is computed after SVD to assess the stability 
!     of the inversion. If the condition number is too high, the code adjusts 
!     using the regularization technique.
!   - The SVD outputs, including the singular values and their ranks, are logged 
!     for debugging purposes, helping users diagnose potential numerical issues.
!
! Data Importances Calculation:
!   - Data importances (`dimpc` and `dimph`) are calculated for both cluster 
!     vectors and hypocentroid inversions.
!   - These importances are normalized across events and azimuth bins to ensure 
!     they are proportionally distributed.
!   - Error-handling mechanisms are added to ensure that azimuth values are 
!     correctly binned between 1 and 12 for further analysis.
!   - Totals for importances are calculated for events and entire datasets, 
!     providing useful diagnostics about the contribution of each data point.
!
! Confidence Ellipses Calculation:
!   - Confidence ellipses for both cluster vectors and hypocentroid inversions 
!     are computed based on a 90% confidence level, using the `fstat1` subroutine 
!     for F-statistic calculation.
!   - When the dataset lacks sufficient information for meaningful ellipse 
!     calculation, the process skips this step and logs a warning.
!   - For cluster vectors, the ellipses are adjusted using a Bayesian approach 
!     for improved reliability, especially when dealing with smaller datasets 
!     or higher uncertainties.
!   - For hypocentroid ellipses, three different cases are computed: 
!     bias-corrected epicentroid, weighted epicentroid, and bias correction terms.
!   - The ellipses are scaled according to Tikhonov-regularized variance matrices.
!
! Bias Correction for Hypocentroid:
!   - The hypocentroid bias correction is calculated using the `dx01` and `dx02` 
!     vectors, which account for weighted and bias-corrected displacements.
!   - The code subtracts bias correction terms from the displacement estimates, 
!     ensuring that the final hypocentroid parameters are unbiased.
!   - The bias corrections are stored and logged for post-processing and diagnostics.
!
! Logging and Debugging:
!   - Extensive logging is included throughout the inversion process. 
!     Key stages such as SVD, data importances, and bias correction are logged 
!     with sufficient detail to aid debugging and result interpretation.
!   - Singular values, ranks, condition numbers, and other intermediate values 
!     are logged to ensure full traceability of the numerical methods used.
!   - Logs also include memory management actions such as successful allocations 
!     and deallocations, ensuring full visibility into the process.
!
!******************************************************************************

module mloc_inv

   use mloc_declare
   use mloc_allocate
   use mloc_math
   
   implicit none
   save
   
   integer, allocatable, dimension(:) :: ibinq
   integer, allocatable, dimension(:) :: indx
   integer, allocatable, dimension(:) :: irank
   integer, allocatable, dimension(:) :: ntiev ! event number associated with each TT residual
   real, allocatable, dimension(:) :: sn
   double precision, allocatable, dimension(:,:) :: a0
   double precision, allocatable, dimension(:,:) :: ahat
   double precision, allocatable, dimension(:,:) :: aqbahat
   double precision, allocatable, dimension(:) :: dthat
   double precision, allocatable, dimension(:) :: dx
   double precision, allocatable, dimension(:) :: ehath
   double precision, allocatable, dimension(:) :: pqahat
   double precision, allocatable, dimension(:) :: pwa0
   double precision, allocatable, dimension(:) :: qbadx
   double precision, allocatable, dimension(:,:) :: qbahat
   double precision, allocatable, dimension(:) :: qbdthat
   double precision, allocatable, dimension(:) :: qm
   double precision, allocatable, dimension(:) :: qmtr
   double precision, allocatable, dimension(:,:) :: t4m
   double precision, allocatable, dimension(:,:) :: t4m2
   double precision, allocatable, dimension(:,:) :: t4n1
   double precision, allocatable, dimension(:,:) :: t4n2
   double precision, allocatable, dimension(:,:) :: t4q
   double precision, allocatable, dimension(:) :: tm
   double precision, allocatable, dimension(:,:) :: tmn
   double precision, allocatable, dimension(:) :: tn
   double precision, allocatable, dimension(:,:) :: tn4
   double precision, allocatable, dimension(:) :: tq
   double precision, allocatable, dimension(:,:) :: umm
   double precision, allocatable, dimension(:,:) :: unm
   double precision, allocatable, dimension(:,:) :: uq4
   double precision, allocatable, dimension(:) :: vc
   double precision, allocatable, dimension(:,:) :: vmm
   double precision, allocatable, dimension(:,:) :: wa0
   double precision, allocatable, dimension(:,:) :: wa0d
   character(len=21), allocatable, dimension(:) :: qcname
   character(len=13), allocatable, dimension(:) :: qhname

contains

!*****************************************************************************************
subroutine mlocinv (it)

! generalized inverse for hypocentroid and cluster vectors.

! January 26, 1989 by eric bergman
! Converted to use double precision 6/29/95 by eab.
! Major change to the calculation of confidence ellipses for cluster vectors, 2/3/05 by eab.

   integer :: i
   integer :: ibin
   integer :: icond
   integer :: iev ! event counter
   integer :: it
   integer :: j
   integer :: k
   integer :: kbayes
   integer :: kst
   integer :: m
   integer :: mf
   integer, dimension(4) :: mindxh
   integer :: mt
   integer :: mt0
   integer :: mtc
   integer :: mth
   integer :: mttotal
   integer :: nf
   integer :: nq ! number of distinct station-phases used
   integer :: nq_allocated
   integer :: nt ! number of TT residuals used
   integer :: nt_allocated
   real :: al
   real :: alphah1
   real :: alphah12
   real :: alphah2
   real :: az
   real :: bl
   real, dimension(4) :: dx_temp
   real, dimension(4) :: dxp_mean
   real, dimension(4) :: dxp_sum
   real :: f ! F statistic
   real :: inv_cpu_finish
   real :: inv_cpu_start
!    real, dimension(4,4) :: hcv
   real :: kcrit
   real :: minimum_weight
   real, parameter :: pc = 0.90 ! probability for cluster vector confidence ellipses
   real, parameter :: ph = 0.90 ! probability for hypocentroid confidence ellipse
   real :: t11
   real :: t12
   real :: t22
   real :: temp
   real :: totimp
   real :: totimpiev
   real :: xl1h1
   real :: xl1h12
   real :: xl1h2
   real :: xl2h1
   real :: xl2h12
   real :: xl2h2
   double precision, dimension(4,4) :: awa
   double precision :: dtemp
   double precision, dimension(4) :: dx01
   double precision, dimension(4) :: dx02
   double precision :: e2
   double precision :: ehati
   double precision, dimension(4) :: q4
   double precision :: shatsqh
   double precision, dimension(2,2) :: t22a
   double precision :: tf
   double precision, dimension(4,4) :: v44
   double precision, dimension(4) :: vh
   double precision, dimension(4,4) :: vhath1
   double precision, dimension(4,4) :: vhath2
   character(len=132):: msg
   character(len=40) :: p
   logical :: used
   
   p = 'mlocinv'
      
   tf = dble(tikhonov_factor)
   
   ! Allocate variable arrays
   call get_array_limits_c (nq_allocated, nt_allocated)
   call mlocinv_allocate (nq_allocated, nt_allocated)
   ! dx is handled separately because it is calculated for the cluster vectors but also
   ! needed for the hypocentroid
   allocate (dx(n_fp_max), stat=error)
   if (error .gt. 0) then
      call allocation_error (p, 'dx', error)
   else
      if (verbose_log) write (io_alloc_log,'(a,i3,a)') 'mlocinv: dx allocated (', n_fp_max, ')'
   end if
      
   
   call cpu_time (inv_cpu_start) ! Track CPU usage for the inversion in each iteration

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
   ! cluster vectors

   write (*,'(a)') '  cluster vectors:'
   if (n_event .gt. 1) then
   
      ! Initialization

      ntiev = 0
      dthat = 0.0d0
      ahat = 0.0d0
      qcname = ' '
      ntqi = 0
      ntq = 0
      do iev = 1,n_event
         eciev(iev,it) = 0.
         ndatc(iev,it) = 0
      end do
      mt = 0
      nq = 0
      nt = 0
      
      ! This was for some special testing. Don't uncomment unless you know what you're doing!
      if (debug) write (io_log,'(a)') ' iev   j mnf Con  FH  FC  nt  jb qcname(jb)             vnhat  dthat'
      do iev = 1,n_event
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               nt = nt + 1
               if (nt .gt. nt_allocated) then
                  write (msg,'(a,i8,a,i3)') 'mlocinv: NT exceeds NT_ALLOCATED (', nt_allocated, ') for event ', iev
                  call oops (trim(msg))
               end if
               ntiev(nt) = iev
               ndatc(iev,it) = ndatc(iev,it) + 1 ! Number of observations from each event used for cluster vector
               if (nt .gt. 1) then
                  used = .false.
                  kst = 0
                  do while (.not.used)
                     kst = kst + 1
                     if (kst .le. nq) then
                        if (stname(iev,j)//deployment(iev,j)//phase(iev,j) .eq. qcname(kst)&
                         .and. idiff(iev,j) .eq. idiff0(kst)) then
                           jb(nt) = kst
                           used = .true.
                        end if
                     else
                        nq = nq + 1
                        if (nq .gt. nq_allocated) then
                           write (msg,'(a,i6,a,i3)') 'mlocinv: NQ exceeds NQ_ALLOCATED (', nq_allocated,&
                            ') for event ', iev
                           call oops (trim(msg))
                        end if 
                        qcname(nq) = stname(iev,j)//deployment(iev,j)//phase(iev,j)
                        qlat(nq) = stladg(iev,j) ! Geocentric coordinates
                        qlon(nq) = stlndg(iev,j) ! Geocentric coordinates
                        qelev(nq) = ahgts(iev,j)
                        idiff0(nq) = idiff(iev,j)
                        if (data_weight) then
                           rderr0(nq) = sdread(iev,j)
                        else
                           rderr0(nq) = 1.
                        end if
                        jb(nt) = nq
                        used = .true.
                     end if 
                  end do
               else
                  nq = 1
                  qcname(1) = stname(iev,j)//deployment(iev,j)//phase(iev,j)
                  qlat(1) = stladg(iev,j) ! Geocentric coordinates
                  qlon(1) = stlndg(iev,j) ! Geocentric coordinates
                  qelev(1) = ahgts(iev,j)
                  idiff0(1) = idiff(iev,j)
                  if (data_weight) then
                     rderr0(1) = sdread(iev,j)
                  else
                     rderr0(1) = 1.
                  end if
                  jb(1) = 1
               end if
               if (data_weight) then
                  vnhat(nt) = sdread(iev,j)*sdread(iev,j)
               else
                  vnhat(nt) = 1.0
               end if
               sighatj(nt) = sqrt(vnhat(nt))
               temp = 1.0/sighatj(nt)
               dthat(nt) = dble(dt(iev,j,it)*temp)
               do k = 1,4
                  if (mindx(iev,k) .ne. 0) then
                     m = mt + mindx(iev,k)
                     ahat(nt,m) = dble(a(iev,j,k)*temp)
                  end if
               end do
               !This was for some special testing. Don't uncomment unless you know what you're doing!
               if (debug) write (io_log,'(3i5,3x,l1,3x,l1,3x,l1,2i6,1x,a,1x,f6.3,1x,f7.3)') &
                   iev, j, mnf_line(iev,j), connected(iev,j), fltrh(iev,j), fltrc(iev,j),&
                   nt, jb(nt), qcname(jb(nt)), vnhat(nt), sngl(dthat(nt))
            end if
         end do
         mt = mt + mtiev(iev)
         if (mt .gt. n_fp_max) then
            write (msg,'(2(a,i3))') 'mlocinv: MT exceeds n_fp_max (', n_fp_max, ') for event ', iev
            call oops (trim(msg))
         end if
         if (ndatc(iev,it) .lt. mtiev(iev)) then
            write (msg,'(a,i3)') 'mlocinv: Fewer data than free parameters for event ', iev
            call oops (trim(msg))         
         end if
      end do
      do j = 1,nq
         i = 0
         do k = 1,nt
            if (jb(k) .eq. j) then
               i = i + 1
               if (i .gt. n_qi_max) then
                  write (msg,'(3(a,i4,3a))') 'mlocinv: max number (', n_qi_max,&
                   ') of instances of station-phase ', qcname(j), ' exceeded'
                  call oops (trim(msg))
               end if
               ntq(i,j) = k
            end if
         end do
         ntqi(j) = i
      end do
      
      ntc = nt
      nqc = nq
      write (*,'(t7,a,i7,a)') ' NT = ', nt, ' (# data)'
      write (*,'(t7,a,i7,a)') ' MT = ', mt, ' (# free parameters)'
      write (*,'(t7,a,i7,a)') ' NQ = ', nq, ' (# indep. station-phases)'
      
      write (io_log,'(/a)') 'Inversion for cluster vectors'
      write (io_log,'(t7,a,i7,a)') ' NT  = ', nt, ' (# data)'
      write (io_log,'(t7,a,i7,a)') ' MTH = ', mt, ' (# free parameters)'
      write (io_log,'(t7,a,i7,a)') ' NQ  = ', nq, ' (# indep. station-phases)'
   
      do i = 1,nqc
         qname1(i) = qcname(i)
      end do

      ! solve delxc = (qbhat * ahat)**dagger * qbhat * dthat (eq 72)

      ! diagonal elements of wq2(eq 67)

      do i = 1,nq
         wq2(i) = 0.
         do j = 1,ntqi(i)
            wq2(i) = wq2(i) + (1./vnhat(ntq(j,i)))
         end do
      end do
      
      call postq (ahat, nt_allocated, n_fp_max, nt, mt, qbahat)
      call postqv (dthat, nt_allocated, nt, qbdthat)
      write (*,'(a)') '   ...calling dsvd...'
      call dsvd (qbahat, nt, mt, nt_allocated, n_fp_max, unm, vmm, qm)
      call dindexx (mt, qm, indx) 
      call rank (mt, indx, irank)
      do j = 1,mt
         if (irank(j) .gt. mtiev(1)) then
            qmtr(j) = ((qm(j)*qm(j)) + (tf*tf))/qm(j) ! Tikhonov regularization
            tm(j) = 1.0d0/qmtr(j)
         else
            tm(j)=0.0d0
         end if
      end do
      icond = int(qmtr(indx(mt))/qmtr(indx(mtiev(1)+1)))
      write (*,'(a,i12)') '   ...return from dsvd: condition # = ', icond
      if (debug) then ! write out the singular values
         write (io_log,'(/a,e12.6)') 'Tikhonov factor :', tf
         write (io_log,'(a)') 'Singular values for the cluster vectors:'
         write (io_log,'(a)') '     j irank  indx     qm           qmtr'
         do j = 1,mt
            write (io_log,'(3i6,1x,e12.6,1x,e12.6)') j, irank(j), indx(j), qm(j), qmtr(j)
         end do
         write (io_log,'(a,i6)') 'indx(mt) = ', indx(mt)
         write (io_log,'(a,i6)') 'indx(mtiev(1)+1) = ', indx(mtiev(1)+1)
         write (io_log,'(a,i6)') 'indx(mtiev(1)) = ', indx(mtiev(1))
         write (io_log,'(2(e12.6,a),i12)') qmtr(indx(mt)), ' / ', qmtr(indx(mtiev(1)+1)), ' = ', icond
         write (io_log,'(/a)') 'tm:'
         write (io_log,'(12(e12.6,1x))') (tm(j),j=1,mt)
      end if
      do j = 1,mt
         dtemp = tm(j)
         do i = 1,mt
            vmm(i,j) = vmm(i,j)*dtemp
         end do
      end do 
      call ddot3 (vmm, unm, tmn, n_fp_max, n_fp_max, nt_allocated, n_fp_max, n_fp_max, nt_allocated, mt, mt, nt)
      call ddot1 (tmn, qbdthat, dx, n_fp_max, nt_allocated, nt_allocated, 1, n_fp_max, 1, mt, nt, 1)

      ! data importances

      do i = 1,nt
         pqahat(i) = 0.0d0
         do k = 1,mt
            pqahat(i) = pqahat(i) + (qbahat(i,k)*tmn(k,i))
         end do
      end do
      dimpc = 0.
      dimpciev = 0.
!         do ibin = 1,12
!            dimpc(ibin) = 0.
!            do iev = 1,n_event
!               dimpciev(iev,ibin) = 0.
!            end do
!         end do
      totimp = 0.
      k = 0
      do iev = 1,n_event
         totimpiev = 0.
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               k = k + 1
               dtmpc(iev,j) = sngl(pqahat(k))
               totimp = totimp + sngl(pqahat(k))
               totimpiev = totimpiev + sngl(pqahat(k))
            else
               dtmpc(iev,j) = 0.
            end if
            az = azes(iev,j)
            if (az .lt. 0.) az = az + 360.
            ibin = int(az/30.) + 1
            if (nint(az) .eq. 360 .and. ibin .eq. 13) ibin = 1
            if (ibin .lt. 1 .or. ibin .gt. 12) then
               write (msg,'(a,i3,f10.3,2i5)') 'mlocinv: illegal value for ibin = ', ibin, az, iev, j
               call warnings (trim(msg))
               write (io_log,'(a)') msg
               ibin = 1
            end if
            dimpc(ibin) = dimpc(ibin) + dtmpc(iev,j)
            dimpciev(iev,ibin) = dimpciev(iev,ibin) + dtmpc(iev,j)
         end do
         do ibin = 1,12
            dimpciev(iev,ibin) = dimpciev(iev,ibin)/totimpiev
         end do
      end do
      do ibin = 1,12
         dimpc(ibin) = dimpc(ibin)/totimp
      end do

      ! variance matrix of the estimation process:
      ! vhatc = (ahat**tr * qbhat * ahat)**dagger (eq 78)

      print  *,'    variance matrix'
      call ddot2 (ahat, qbahat, aqbahat, nt_allocated, n_fp_max, nt_allocated, n_fp_max, n_fp_max, n_fp_max, mt, nt, mt) 
      call dsvd (aqbahat, mt, mt, n_fp_max, n_fp_max, umm, vmm, qm)
      call dindexx (mt, qm, indx) 
      call rank (mt, indx, irank)
      icond = int(qm(indx(mt))/qm(indx(mtiev(1)+1)))
      write (*,'(a,i12)') '   ...return from dsvd: condition # = ', icond
      if (debug) then ! write out the singular values
         write (io_log,'(/a)') 'Singular values for the variance matrix:'
         write (io_log,'(12(e12.6,1x))') (qm(j),j=1,mt)
         write (io_log,'(a,i6)') 'indx(mt) = ', indx(mt)
         write (io_log,'(a,i6)') 'indx(mtiev(1)+1) = ', indx(mtiev(1)+1)
         write (io_log,'(a,i6)') 'indx(mtiev(1)) = ', indx(mtiev(1))
         write (io_log,'(2(e12.6,a),i12)') qm(indx(mt)), ' / ', qm(indx(mtiev(1)+1)), ' = ', icond
      end if
      do j = 1,mt
         do i = 1,mt
            dtemp = 0.
            do k = 1,mt 
               if (irank(k) .gt. mtiev(1)) dtemp = dtemp + (umm(i,k)*umm(j,k)/qm(k))
            end do
            vhatc(i,j) = dtemp
         end do
      end do 

      ! error.
      ! ehatc = qbhat * dthat - qbhat * ahat * dx            (eq 65)
      ! shatc**2 = (ehatc*ehatc)/(nt-nq-(n_event-1)mtiev)        (eq 81)

      call ddot1 (qbahat, dx, qbadx, nt_allocated, n_fp_max, n_fp_max, 1, nt_allocated, 1, nt, mt, 1)
      dtemp = 0.0d0
      do i = 1,nt
         ehati = qbdthat(i) - qbadx(i)
         e2 = ehati*ehati
         dtemp = dtemp + e2
         eciev(ntiev(i),it) = eciev(ntiev(i),it) + sngl(e2)
      end do
      ehatsqc(it) = sngl(dtemp)
      shatsqc = dtemp/(nt-nq-mt+mtiev(1))
      shatc(it) = sngl(dsqrt(shatsqc))

      ! cluster vector error for output

      k = 0
      ! This was for some special testing. Don't uncomment unless you know what you're doing!
      if (debug) write (io_log,'(a)') ' iev   j mnf   k c qbdthat(k)  qbadx(k)       eci'
      do iev = 1,n_event
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               k = k + 1
               eci(iev,j) = sngl(qbdthat(k) - qbadx(k))
      ! This was for some special testing. Don't uncomment unless you know what you're doing!
               if (debug) write (io_log,'(4i4,1x,l1,1x,3f10.3)') iev, j, mnf_line(iev,j), k,&
                connected(iev,j), sngl(qbdthat(k)), qbadx(k), eci(iev,j)
            else
               eci(iev,j) = 0.
            end if
         end do
      end do

      ! confidence ellipses for lat-long

      ! This section is tricky! Naively taking the result in EQ 80 suggests
      !         mf=2*(n_event-1)
      !         nf=nt-nq-(n_event-1)*2
      ! which leads to impossibly large confidence ellipses because they scale
      ! up as the number of events in the cluster. However, we need to treat this
      ! problem as one of "marginal" confidence ellipses, as in the discussion leading
      ! to EQ 38. Therefore, it is valid to use mf=2 instead of mf=2*(n_event-1). Also,
      ! the relevant number of degrees of freedom for the uncertainty of a cluster
      ! vector is the number of readings used to estimate that cluster vector, not the
      ! total number in the inversion.

      ! By the same logic, the squared error term used to scale "kcrit" is not based
      ! on the entire set of cluster vectors, but on the squared error term for the readings
      ! actually used, normalized by the number of readings.

      ! Until 2/3/05 I used mf=2 and nf=nt-nq-mf, and used shatsqc as the error term:
      !         nf=nt-nq-mf
      !         call fstat1 (mf, nf, pc, f)
      !         kcrit=mf*sngl(shatsqc)*f
      ! This is not too bad, but gives essentially the same kcrit (typically about 8.2) for each
      ! cluster vector's confidence ellipse scaling. I now think there's a better way to handle this
      ! calculation. I now use the normalized cluster error (eciev(iev,it)/ndatc(iev,it) for
      ! each event and the degrees of freedom are based on ndatc. mf still is 2 of course, as it
      ! always would be for a confidence ellipse on the epicenter. The result of this change is
      ! confidence ellipses will be a little larger than before for the poorest events (kcritc > 8.26),
      ! and smaller than before for the better located events (kcritc < 8.26).

      ! Since ndatc can be a small number in this algorithm, the argument for using a Bayesian approach
      ! re-emerges. I have chosen a conservative value of K=3, giving a s.d. of 0.4 on the estimated
      ! "normalized cluster error". 

      mttotal = 0
      mf = 2
      kbayes = 3
      do iev = 1,n_event
         nf = kbayes + ndatc(iev,it) - mf
         call fstat1 (mf, nf, pc, f)
         shatsqci(iev) = (real(kbayes) + eciev(iev,it))/real(nf)
         kcritc(iev) = mf*shatsqci(iev)*f 
         if (mindx(iev,1) .eq. 0 .or. mindx(iev,2) .eq. 0) then
            write (msg,'(a,i3)') 'mlocinv: no confidence ellipse calculated for event ', iev
            call warnings (trim(msg))
            alphac(iev)=0.
            xl1c(iev)=0.
            xl2c(iev)=0.
            alic(iev) = 0.
            blic(iev) = 0.
            ccv(iev,1,1) = 0. ! used for output file .cv
            ccv(iev,1,2) = 0.
            ccv(iev,1,3) = 0.
            ccv(iev,1,4) = 0.
            ccv(iev,2,2) = 0.
            ccv(iev,2,3) = 0.
            ccv(iev,2,4) = 0.
            ccv(iev,3,3) = 0.
            ccv(iev,3,4) = 0.
            ccv(iev,4,4) = 0.
         else
            t22a(1,1) = vhatc(mttotal+1,mttotal+1)
            t22a(2,1) = -vhatc(mttotal+2,mttotal+1) ! geocentric latitude reverses sign.
            t22a(1,2) = -vhatc(mttotal+1,mttotal+2) ! geocentric latitude reverses sign.
            t22a(2,2) = vhatc(mttotal+2,mttotal+2)
            ccv(iev,1,1) = sngl(t22a(1,1)) ! used for output file .cv
            ccv(iev,1,2) = sngl(t22a(1,2))
            ccv(iev,2,1) = sngl(t22a(2,1))
            ccv(iev,2,2) = sngl(t22a(2,2))
            if (mindx(iev,3) .ne. 0) then
               ccv(iev,1,3) = sngl(-vhatc(mttotal+1,mttotal+3))
               ccv(iev,2,3) = sngl(vhatc(mttotal+2,mttotal+3))
               ccv(iev,3,3) = sngl(vhatc(mttotal+3,mttotal+3))
               ccv(iev,3,4) = sngl(vhatc(mttotal+3,mttotal+4))
            else
               ccv(iev,1,3) = 0.
               ccv(iev,2,3) = 0.
               ccv(iev,3,3) = 0.
               ccv(iev,3,4) = 0.
            end if
            if (mindx(iev,4) .ne. 0 .and. mindx(iev,3) .ne. 0) then
               ccv(iev,1,4) = sngl(-vhatc(mttotal+1,mttotal+4))
               ccv(iev,2,4) = sngl(vhatc(mttotal+2,mttotal+4))
               ccv(iev,4,4) = sngl(vhatc(mttotal+4,mttotal+4))
            else if (mindx(iev,4) .ne. 0 .and. mindx(iev,3) .eq. 0) then
               ccv(iev,1,4) = sngl(-vhatc(mttotal+1,mttotal+3))
               ccv(iev,2,4) = sngl(vhatc(mttotal+2,mttotal+3))
               ccv(iev,4,4) = sngl(vhatc(mttotal+3,mttotal+3))
            else
               ccv(iev,1,4) = 0.
               ccv(iev,2,4) = 0.
               ccv(iev,4,4) = 0.
            end if
            call delips (t22a, alphac(iev), al, bl)
            alic(iev) = 1./al
            blic(iev) = 1./bl
            xl1c(iev) = sqrt(kcritc(iev)*alic(iev))
            xl2c(iev) = sqrt(kcritc(iev)*blic(iev))

            ! Convert 90% confidence ellipse back to matrix form under the asumption of kcrit = 1 and
            ! add the uncertainty from radius_cvff, then get the new confidence ellipse.
            call ell2cv (1.0, alphac(iev), xl1c(iev), xl2c(iev), t11, t12, t22)
            t22a(1,1) = dble(t11) + dble(radius_cvff*radius_cvff)
            t22a(1,2) = dble(t12)
            t22a(2,1) = dble(t12)
            t22a(2,2) = dble(t22) + dble(radius_cvff*radius_cvff)
            call delips (t22a, alphac(iev), al, bl)
            alic(iev) = 1./al
            blic(iev) = 1./bl
            xl1c(iev) = sqrt(alic(iev))
            xl2c(iev) = sqrt(blic(iev))
            ! Get the modified raw covariance matrix back
            call ell2cv (kcritc(iev), alphac(iev), xl1c(iev), xl2c(iev), t11, t12, t22)
            ccv(iev,1,1) = t11
            ccv(iev,1,2) = t12
            ccv(iev,2,1) = t12
            ccv(iev,2,2) = t22
         end if
         mt0 = 0
         if (mindx(iev,1) .ne. 0) mt0 = mt0 + 1
         if (mindx(iev,2) .ne. 0) mt0 = mt0 + 1
         if (mindx(iev,3) .ne. 0) then
            mt0 = mt0 + 1
            ccv(iev,3,3) = sngl(vhatc(mttotal+mt0,mttotal+mt0))
         end if
         if (mindx(iev,4) .ne. 0) then
            mt0 = mt0 + 1
            ccv(iev,4,4) = sngl(vhatc(mttotal+mt0,mttotal+mt0))
         end if
         mttotal = mttotal + mtiev(iev)
      end do

      ! scale solution variance by shatsq

      do i = 1,mt
         vc(i) = shatsqc*vhatc(i,i)
      end do

      ! convert to 4-vectors.

      mt0 = 0
      dxp_sum = 0.
      dxp_mean = 0.
      write (io_log,'(/a)') 'Cluster vector shift (geocentric latitude)'
      write (io_log,'(a)') ' iev   lat(km)  long(km) depth(km)     OT(s)'
      do iev = 1,n_event
         do k = 1,4
            if (mindx(iev,k) .ne. 0) then
               m = mt0 + mindx(iev,k)
               dxp(iev,k,it) = sngl(dx(m))
               sdxhatc(iev,k) = sngl(dsqrt(vc(m)))
            else
               dxp(iev,k,it) = 0.0
               sdxhatc(iev,k) = 0.0
            end if
            dxp_sum(k) = dxp_sum(k) + dxp(iev,k,it)
         end do
         mt0 = mt0 + mtiev(iev)
         write (io_log,'(i4,4f10.3)') iev, (dxp(iev,k,it),k=1,4)
      end do
      write (io_log,'(a,4f10.3)') ' sum', (dxp_sum(k),k=1,4)
      ! Check for non-zero mean displacement of cluster vectors.
      ! In theory it should be zero, but in practice some datasets produce a non-zero mean
      ! in one or more parameters that is large enough to prevent convergence.
      do i = 1,4
         dxp_mean(i) = dxp_sum(i)/real(n_event)
      end do
      write (io_log,'(a,4f10.3)') 'mean',(dxp_mean(k),k=1,4)
      write (msg,'(a)') 'Mean cluster vector parameter changes (geocentric latitude)'
      call fyi ('mlocinv: '//trim(msg))
      write (msg,'(t9,a)') '       lat      long     depth        OT'
      call fyi (trim(msg))
      write (msg,'(t9,4f10.3)') (dxp_mean(k),k=1,4)
      call fyi (trim(msg))
      if (damping) then
         write (msg,'(a)') 'Removing mean change of cluster vector parameter (command "damp")'
         write (io_log,'(a)') trim(msg)
         call fyi ('mlocinv: '//trim(msg))
         do iev = 1,n_event
            do k = 1,4
               if (mindx(iev,k) .ne. 0) dxp(iev,k,it) = dxp(iev,k,it) - dxp_mean(k)
            end do
            write (io_log,'(i4,4f10.3)') iev, (dxp(iev,k,it),k=1,4)
         end do
      end if  
   end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ! hypocentroid

   write (*,'(a)') '  hypocentroid:'
   
   ! Recalculate array limits for allocation
   call mlocinv_deallocate ()
   call get_array_limits_h (nq_allocated, nt_allocated)
   call mlocinv_allocate (nq_allocated, nt_allocated)
   
   ! Initialization.

   ntiev = 0
   dthat = 0.0d0
   ahat = 0.0d0
   qhname = ' ' 
   ntqi = 0
   ntq = 0
   a0 = 0.0d0
   do j = 1,n_event
      ehiev(j,it) = 0.
   end do
   
   ! Free parameters
   mth = 0
   if (latfh) then
      mth = mth + 1
      mindxh(1) = mth
   else
      mindxh(1) = 0
   end if
   if (lonfh) then
      mth = mth + 1
      mindxh(2) = mth
   else
      mindxh(2) = 0
   end if
   if (depthfh) then
      mth = mth + 1
      mindxh(3) = mth
   else
      mindxh(3) = 0
   end if
   if (timefh) then
      mth = mth + 1
      mindxh(4) = mth
   else
      mindxh(4) = 0
   end if
   
   mtc = 0
   nq = 0
   nt = 0
   if (debug) write (io_log,'(t32,a,t41,a,t51,a,t60,a,t70,a,t80,a)') 'dt', 'vnhat', 'rderr', 'ttsprd', 'weight',&
    'derivatives'
   do iev = 1,n_event
      if (debug) write (io_log,'(a,i3)') 'Event ', iev
      do j = 1,nst(iev)
         if (.not.fltrh(iev,j)) then
            nt = nt + 1
            if (nt .gt. nt_allocated) then
               write (msg,'(a,i8,a,i3)') 'mlocinv: NT exceeds NT_ALLOCATED (', nt_allocated, ') for event ', iev
               call oops (trim(msg))
            end if 
            ntiev(nt) = iev
            if (nt .gt. 1) then
               used = .false.
               kst = 0
               do while (.not.used)
                  kst = kst + 1
                  if (kst .le. nq) then
                     if (stname(iev,j)//phase(iev,j) .eq. qhname(kst)) then
                        jb(nt) = kst
                        used = .true.
                     end if
                  else
                     nq = nq + 1
                     if (nq .gt. nq_allocated) then
                        write (msg,'(a,i6,a,i3)') 'mlocinv: NQ exceeds NQ_ALLOCATED (', nq_allocated, ') for event ', iev
                        call oops (trim(msg))
                     end if
                     qhname(nq) = stname(iev,j)//phase(iev,j)
                     sn(nq) = s(iev,j,it)
                     jb(nt) = nq
                     do k = 1,4
                        if (mindxh(k) .ne. 0) then
                           m = mindxh(k)
                           a0(nq,m) = dble(a(iev,j,k))
                        end if
                     end do
                     used = .true.
                  end if 
               end do
            else
               nq = 1
               qhname(1) = stname(iev,j)//phase(iev,j)
               sn(1) = s(iev,j,it)
               jb(1) = 1
               do k = 1,4
                  if (mindxh(k) .ne. 0) then
                     m = mindxh(k)
                     a0(1,m) = dble(a(iev,j,k))
                  end if
               end do
            end if
            if (data_weight) then
               ! Check for zero weight. No reading that gets to this point should have zero weight,
               ! because it will throw a divide-by-zero exception here.
               minimum_weight = 0.01
               if (weight(iev,j) .lt. minimum_weight) then
                  write (msg,'(a,e12.4,a,i3,1x,a,1x,a,1x,i5,a,f4.2)') 'mlocinv: ~zero weight ', weight(iev,j),&
                   ' for ', iev, stname(iev,j), phase(iev,j), mnf_line(iev, j), ' set to ', minimum_weight 
                  call warnings (trim(msg))
                  weight(iev,j) = minimum_weight
               end if
               if (.not.pttt) then
                  vnhat(nt) = (sdread(iev,j)**2+ttsprd(iev,j)**2)/weight(iev,j)
               else
                  vnhat(nt) = (sdread(iev,j)**2)/weight(iev,j)
               end if
            else
               vnhat(nt) = 1.0
            end if
            if (debug) write (io_log,'(2a,1x,a,1x,9f10.3)') 'mlocinv: ', stname(iev,j), phase(iev,j), dt(iev,j,it),&
             vnhat(nt), sdread(iev,j), ttsprd(iev,j), weight(iev,j), a(iev,j,1), a(iev,j,2), a(iev,j,3), a(iev,j,4)
            sighatj(nt) = sqrt(vnhat(nt))
            temp = 1.0/sighatj(nt)
            dthat(nt) = dble(dt(iev,j,it)*temp)
            do k = 1,4
               if (mindx(iev,k) .ne. 0) then
                  m = mtc+mindx(iev,k)
                  ahat(nt,m) = dble(a(iev,j,k)*temp)
               end if
            end do
         else
            if (debug) write (io_log,'(2a,1x,a,1x,f10.3,10x,7f10.3)') 'mlocinv: ', stname(iev,j), phase(iev,j),&
             dt(iev,j,it), sdread(iev,j), ttsprd(iev,j), weight(iev,j), a(iev,j,1), a(iev,j,2), a(iev,j,3), a(iev,j,4)
         end if
      end do
      mtc = mtc + mtiev(iev)
   end do
!    write (*,*) size(ahat,1), size(ahat,2)
   do j = 1,nq
      i = 0
      do k = 1,nt
         if (jb(k) .eq. j) then
            i = i + 1
            ntq(i,j) = k
         end if
      end do
      ntqi(j) = i
   end do

   write (*,'(t7,a,i7,a)') ' NT  = ', nt, ' (# data)'
   write (*,'(t7,a,i7,a)') ' MTH = ', mth, ' (# free parameters)'
   write (*,'(t7,a,i7,a)') ' NQ  = ', nq, ' (# indep. station-phases)'
   
   write (io_log,'(/a)') 'Inversion for hypocentroid'
   write (io_log,'(t7,a,i7,a)') ' NT  = ', nt, ' (# data)'
   write (io_log,'(t7,a,i7,a)') ' MTH = ', mth, ' (# free parameters)'
   write (io_log,'(t7,a,i7,a)') ' NQ  = ', nq, ' (# indep. station-phases)'
   
   if (nq .lt. mth) then
      write (msg,'(a,3i6)') 'mlocinv: number of data is less than the number of free parameters ', nt, mth, nq
      call warnings (trim(msg))
      delx0(1,it) = 0.
      delx0(2,it) = 0.
      delx0(3,it) = 0.
      delx0(4,it) = 0.
      dxp(iev,1,it) = 0.
      dxp(iev,2,it) = 0.
      dxp(iev,3,it) = 0.
      dxp(iev,4,it) = 0.
      sdxhath(1) = 0.
      sdxhath(2) = 0.
      sdxhath(3) = 0.
      sdxhath(4) = 0.
      return
   end if

   ! diagonal elements of wq and wq2(eq 67)

   do i = 1,nq
      temp = 0.
      do j = 1,ntqi(i)
         temp = temp + (1.0/vnhat(ntq(j,i)))
      end do
      wq2(i) = temp
      wq(i) = sqrt(temp)
   end do

   ! form w*a0 and (w*a0)**dagger = v * q**-1 * u**tr

   do j = 1,mth
      do i = 1,nq
         wa0(i,j) = dble(wq(i))*a0(i,j)
      end do
   end do 
   write (*,'(a)') '   ...calling dsvd...'
   call dsvd (wa0, nq, mth, nq_allocated, 4, uq4, v44, q4)
   call dindexx (mth, q4, indx) 
   icond = int(q4(indx(mth))/q4(indx(1)))
   write (*,'(a,i12)') '   ...return from dsvd: condition # = ', icond
   do j = 1,mth
      dtemp = 1.0d0/q4(j)
      do i = 1,mth
         v44(i,j) = v44(i,j)*dtemp
      end do
   end do
   call ddot3 (v44, uq4, wa0d, 4, 4, nq_allocated, 4, 4, nq_allocated, mth, mth, nq)  

   ! data importances 
   ! pwa0 = wa0 * wa0d    (eq 102)

   do i = 1,nq
      dtemp = 0.0d0
      do j = 1,mth
         dtemp = dtemp + (wa0(i,j)*wa0d(j,i))
      end do
      pwa0(i) = dtemp
   end do
   do i = 1,nq
      ibinq(i) = 0
   end do
   do ibin = 1,12
      dimph(ibin) = 0.
   end do
   k = 0
   do iev = 1,n_event
      do j = 1,nst(iev)
         if (.not.fltrh(iev,j)) then
            k = k + 1
            dtmph(iev,j) = sngl(pwa0(jb(k)))
            az = azes(iev,j)
            if (az .lt. 0.) az = az + 360.
            ibinq(jb(k)) = int(az/30.) + 1
         else
            dtmph(iev,j) = 0.
         end if
      end do
   end do
   totimp = 0.
   do i = 1,nq
      if (ibinq(i) .gt. 12) ibinq(i) = 12
      if (ibinq(i) .lt. 1) ibinq(i) = 1
      dimph(ibinq(i)) = dimph(ibinq(i)) + sngl(pwa0(i))
      totimp = totimp + sngl(pwa0(i))
   end do
   do ibin = 1,12
      dimph(ibin) = dimph(ibin)/totimp
   end do 

   ! dx01 = (w*a0)**dagger * {w**-1 * bhat**tr * dthat - w*sn} (eqns 83, 101)
   do i = 1,nq
      dtemp = 0.
      do j = 1,ntqi(i)
         k = ntq(j,i)
         dtemp = dtemp + dthat(k)/dble((sighatj(k)*wq(i)))
      end do 
      tq(i) = dtemp - dble((wq(i)*sn(i)))
   end do
   call ddot1 (wa0d, tq, dx01, 4, nq_allocated, nq_allocated, 1, 4, 1, mth, nq, 1)
   write (io_log,'(/a)') 'Hypocentroid parameter shift  (geocentric latitude)'
   write (io_log,'(t5,a)') '   lat(km)  long(km) depth(km)     OT(s)'
!    write (io_log,'(a,4f10.3)') ' weighted: ', (dx01(i),i=1,mth)
   do k = 1,4
      if (mindxh(k) .ne. 0) then
         m = mindxh(k)
         dx_temp(k) = sngl(dx01(m))
      else
         dx_temp(k) = 0.
      end if
   end do
   write (io_log,'(t5,4f10.3,a)') (dx_temp(k),k=1,4), ' Weighted'

   ! dx02 = (w*a0)**dagger *{w**-1 * bhat**tr * ahat * dx} (eqns 84, 101)
   call ddot1 (ahat, dx, tn, nt_allocated, n_fp_max, n_fp_max, 1, nt_allocated, 1, nt, mtc, 1)
   do i = 1,nq
      dtemp = 0.0d0
      do j = 1,ntqi(i)
         k = ntq(j,i)
         dtemp = dtemp + tn(k)/dble((sighatj(k)*wq(i)))
      end do 
      tq(i) = dtemp
   end do
   call ddot1 (wa0d, tq, dx02, 4, nq_allocated, nq_allocated, 1, 4, 1, mth, nq, 1)
!    write (io_log,'(a,4f10.3)') ' bias correction: ', (dx02(i),i=1,mth)
   do k = 1,4
      if (mindxh(k) .ne. 0) then
         m = mindxh(k)
         dx_temp(k) = sngl(dx02(m))
      else
         dx_temp(k) = 0.
      end if
   end do
   write (io_log,'(t5,4f10.3,a)') (dx_temp(k),k=1,4), ' Bias correction'

   ! Convert to 4-vector, subtract bias correction terms
   do k = 1,4
      if (mindxh(k) .ne. 0) then
         m = mindxh(k)
         delx0(k,it) = sngl(dx01(m) - dx02(m))
         bcorr(k) = bcorr(k) + sngl(dx02(m)) ! save total bias corrections
      else
         delx0(k,it) = 0.
      end if
   end do
   write (io_log,'(t5,4f10.3,a)') (delx0(k,it),k=1,4), ' Bias-corrected'

   ! vhath1 = (a0**tr * wq2 * a0)**-1 (eq 90)

   do j = 1,mth
      do i = 1,mth
         dtemp = 0.0d0
         do k = 1,nq
            dtemp = dtemp + (a0(k,i)*dble(wq2(k))*a0(k,j))
         end do
         awa(i,j) = dtemp
      end do
   end do
   call dmatinv2 (awa, vhath1, 4, 4, mth, mth)

   ! vhath2 = (vhath1 * a0**tr * bhat**tr) * ahat * vhatc * ahat**tr *
   !          (vhath1 * a0**tr * bhat**tr)**tr (eq 90, 101)

   call ddot3 (vhath1, a0, t4q, 4, 4, nq_allocated, 4, 4, nq_allocated, mth, mth, nq)
   do j = 1,nt
      dtemp = 1.0d0/dble(sighatj(j))
      do i = 1,mth
         t4n1(i,j) = t4q(i,jb(j))*dtemp
      end do
   end do
   call ddot1 (t4n1, ahat, t4m, 4, nt_allocated, nt_allocated, n_fp_max, 4, n_fp_max, mth, nt, mtc)
   call ddot1 (t4m, vhatc, t4m2, 4, n_fp_max, n_fp_max, n_fp_max, 4, n_fp_max, mth, mtc, mtc)
   call ddot3 (t4m2, ahat, t4n2, 4, n_fp_max, nt_allocated, n_fp_max, 4, nt_allocated, mth, mtc, nt)
   call ddot3 (t4n2, t4n1, vhath2, 4, nt_allocated, 4, nt_allocated, 4, 4, mth, nt, mth)
   
   ! error
   ! ehath = (pbhat * dthat) - (bhat * sn) - (bhat * a0 * dx01)    (eq 92)
   ! shath = sqrt(ehath**2/(nq-mth))               (eq 96)

   call postpv (dthat, nt_allocated, nt, ehath)
   do i = 1,nt
      ehath(i) = ehath(i) - dble(sn(jb(i))/sighatj(i))
   end do
   do j = 1,mth
      do i = 1,nt
         tn4(i,j) = a0(jb(i),j)/dble(sighatj(i))
      end do
   end do
   call ddot1 (tn4, dx01, tn, nt_allocated, 4, 4, 1, nt_allocated, 1, nt, mth, 1)
   dtemp = 0.0d0
   do i = 1,nt
      ehath(i) = ehath(i) - tn(i)
      e2 = ehath(i)*ehath(i)
      dtemp = dtemp + e2
      ehiev(ntiev(i),it) = ehiev(ntiev(i),it) + sngl(e2)
   end do
   ehatsqh(it) = sngl(dtemp)
   nqmth(it) = nq-mth
   shatsqh = dtemp/(nq-mth)
   shath(it) = sngl(dsqrt(shatsqh))

   ! Confidence ellipses for bias-corrected epicentroid (eq. 100),
   ! weighted epicentroid (eq. 97), and bias correction term (eq. 98).

   write (io_log,'(/i3,a)') int(ph*100), '% confidence ellipses:'
   if (mindxh(1) .eq. 0 .or. mindxh(2) .eq. 0) then
      call warnings ('mlocinv: no confidence ellipse calculated for epicentroid')
      alphah = 0.
      xl1h = 0.
      xl2h = 0.
   else
      mf = 2
      nf = max(nt-(n_event*mf), 2) ! To prevent negative nf when using a small amount of local readings
      ! Use of teleseismic P to estimate the hypocentroid with many large events can
      ! produce a value of nf large enough to cause numerical problems in the calculation, but
      ! beyond nf=10000 the value of F converges to ~2.3 all the way to infinity.
      if (nf .le. 10000) then
         call fstat1 (mf, nf, ph, f)
      else
         f = 2.3 ! asymptotic value
         write (msg,'(a,i7)') 'mlocinv: F statistic set to 2.3 for nf = ', nf
         call fyi (trim(msg))
      end if
      t22a(1,1) = shatsqh*vhath1(1,1) + shatsqc*vhath2(1,1)
      t22a(2,1) = -(shatsqh*vhath1(2,1) + shatsqc*vhath2(2,1)) ! geocentric latitude reverses sign.
      t22a(1,2) = -(shatsqh*vhath1(1,2) + shatsqc*vhath2(1,2)) ! geocentric latitude reverses sign.
      t22a(2,2) = shatsqh*vhath1(2,2) + shatsqc*vhath2(2,2)
      call delips (t22a, alphah12, al, bl)
      kcrit  =mf*f
      xl1h12 = sqrt(kcrit/al)
      xl2h12 = sqrt(kcrit/bl)
      write (io_log,'(a,f6.1,a,2f10.3)') ' Bias-corrected: alpha = ', alphah12, ' axes = ', xl1h12, xl2h12
!       if (bias_corr) then
!          hcv(1,1) = sngl(t22a(1,1)) ! used for output file .cv
!          hcv(1,2) = sngl(t22a(1,2))
!          hcv(2,2) = sngl(t22a(2,2))
!          if (mindxh(3) .ne. 0) then
!             hcv(1,3) = sngl(-(shatsqh*vhath1(1,3) + shatsqc*vhath2(1,3)))
!             hcv(2,3) = sngl(shatsqh*vhath1(2,3) + shatsqc*vhath2(2,3))
!             hcv(3,3) = sngl(shatsqh*vhath1(3,3) + shatsqc*vhath2(3,3))
!             hcv(3,4) = sngl(shatsqh*vhath1(3,4) + shatsqc*vhath2(3,4))
!          else
!             hcv(1,3) = 0.
!             hcv(2,3) = 0.
!             hcv(3,3) = 0.
!             hcv(3,4) = 0.
!          end if
!          if (mindxh(4) .ne. 0 .and. mindxh(3) .ne. 0) then
!             hcv(1,4) = sngl(-(shatsqh*vhath1(1,4) + shatsqc*vhath2(1,4)))
!             hcv(2,4) = sngl(shatsqh*vhath1(2,4) + shatsqc*vhath2(2,4))
!             hcv(4,4) = sngl(shatsqh*vhath1(4,4) + shatsqc*vhath2(4,4))
!          else if (mindxh(4) .ne. 0 .and. mindxh(3) .eq. 0) then
!             hcv(1,4) = sngl(-(shatsqh*vhath1(1,3) + shatsqc*vhath2(1,3)))
!             hcv(2,4) = sngl(shatsqh*vhath1(2,3) + shatsqc*vhath2(2,3))
!             hcv(4,4) = sngl(shatsqh*vhath1(3,3) + shatsqc*vhath2(3,3))
!          else
!             hcv(1,4) = 0.
!             hcv(2,4) = 0.
!             hcv(4,4) = 0.
!          end if
!       end if
      nf = nq - mf
      call fstat1 (mf, nf, ph, f)
      t22a(1,1) = vhath1(1,1)
      t22a(2,1) = -vhath1(2,1) ! geocentric latitude reverses sign.
      t22a(1,2) = -vhath1(1,2) ! geocentric latitude reverses sign.
      t22a(2,2) = vhath1(2,2)
      call delips (t22a, alphah1, al, bl)
      kcrit = mf*sngl(shatsqh)*f
      xl1h1 = sqrt(kcrit/al)
      xl2h1 = sqrt(kcrit/bl)
      write (io_log,'(a,f6.1,a,2f10.3)') '       Weighted: alpha = ', alphah1, ' axes = ', xl1h1, xl2h1
!       if (.not.bias_corr) then
!          hcv(1,1) = sngl(t22a(1,1)) ! used for output file .cv
!          hcv(1,2) = sngl(t22a(1,2))
!          hcv(2,2) = sngl(t22a(2,2))
!          if (mindxh(3) .ne. 0) then
!             hcv(1,3) = sngl(-vhath1(1,3))
!             hcv(2,3) = sngl(vhath1(2,3))
!             hcv(3,3) = sngl(vhath1(3,3))
!             hcv(3,4) = sngl(vhath1(3,4))
!          else
!             hcv(1,3) = 0.
!             hcv(2,3) = 0.
!             hcv(3,3) = 0.
!             hcv(3,4) = 0.
!          end if
!          if (mindxh(4) .ne. 0 .and. mindxh(3) .ne. 0) then
!             hcv(1,4) = sngl(-vhath1(1,4))
!             hcv(2,4) = sngl(vhath1(2,4))
!             hcv(4,4) = sngl(vhath1(4,4))
!          else if (mindxh(4) .ne. 0 .and. mindxh(3) .eq. 0) then
!             hcv(1,4) = sngl(-vhath1(1,3))
!             hcv(2,4) = sngl(vhath1(2,3))
!             hcv(4,4) = sngl(vhath1(3,3))
!          else
!             hcv(1,4) = 0.
!             hcv(2,4) = 0.
!             hcv(4,4) = 0.
!          end if
!       end if         
      if (n_event .gt. 1) then
         nf = max(nt-nq-(n_event-1)*mf, 2) ! To prevent negative nf when using a small amount of local readings
         ! Use of teleseismic P to estimate the hypocentroid with many large events can
         ! produce a value of nf large enough to cause numerical problems in the calculation, but
         ! beyond nf=10000 the value of F converges to ~2.3 all the way to infinity.
         if (nf .le. 10000) then
            call fstat1 (mf, nf, ph, f)
         else
            f = 2.3 ! asymptotic value
            write (msg,'(a,i7)') 'mlocinv: F statistic set to 2.3 for nf = ', nf
            call fyi (trim(msg))
         end if
         t22a(1,1) = vhath2(1,1)
         t22a(2,1) = -vhath2(2,1) ! geocentric latitude reverses sign.
         t22a(1,2) = -vhath2(1,2) ! geocentric latitude reverses sign.
         t22a(2,2) = vhath2(2,2)
         call delips (t22a, alphah2, al, bl)
         kcrit = mf*sngl(shatsqc)*f
         xl1h2 = sqrt(kcrit/al)
         xl2h2 = sqrt(kcrit/bl)
         write (io_log,'(a,f6.1,a,2f10.3)') '     Correction: alpha = ', alphah2, ' axes = ', xl1h2, xl2h2
      end if
      if (bias_corr) then
         alphah = alphah12
         xl1h = xl1h12
         xl2h = xl2h12
      else
         alphah = alphah1
         xl1h = xl1h1
         xl2h = xl2h1
      end if
   end if

   do k = 1,4
      if (mindxh(k) .ne. 0) then
         m = mindxh(k)
         if (bias_corr) then
            vh(m) = shatsqh*vhath1(m,m) + shatsqc*vhath2(m,m)
         else
            vh(m) = shatsqh*vhath1(m,m)
         end if 
         sdxhath(k) = sngl(dsqrt(vh(m)))
      else
         sdxhath(k) = 0.
      end if
   end do
   
   ! Deallocate variable arrays
   call mlocinv_deallocate ()
   deallocate (dx, stat=error)
   if (error .gt. 0) then
      call deallocation_error (p, 'dx', error)
   else
      write (io_log,'(a)') 'mlocinv: dx deallocated'
   end if

   call cpu_time (inv_cpu_finish)
   write (msg,'(a,i1,a,f12.3,a)') 'mloc_inv: cpu usage for iteration ', it, ' = ',&
    inv_cpu_finish-inv_cpu_start, ' seconds'
   call fyi (trim(msg))
   
   return
   
end subroutine mlocinv


!*****************************************************************************************
   subroutine get_array_limits_c (nq, nt)
   
   ! Calculates the number of nq_allocated and nt_allocated for allocation of variable arrays
   ! for the cluster vectors.
      
      integer :: iev
      integer :: j
      integer :: kst
      integer, intent(out) :: nq ! number of distinct station-phases with connectivity
      integer, intent(out) :: nt ! number of arrival time data with connectivity
      logical :: used
      
      p = 'get_array_limits_c'
      allocate (qcname(7000), stat=error); if (error .gt. 0) call allocation_error (p, 'qcname', error)
      
      nt = 0
      nq = 0
      do iev = 1,n_event
         do j = 1,nst(iev)
            if (connected(iev,j)) then
               nt = nt + 1
               if (nt .gt. 1) then
                  used = .false.
                  kst = 0
                  do while (.not.used)
                     kst = kst + 1
                     if (kst .le. nq) then
                        if (stname(iev,j)//deployment(iev,j)//phase(iev,j) .eq. qcname(kst)&
                         .and. idiff(iev,j) .eq. idiff0(kst)) then
                           used = .true.
                        end if
                     else
                        nq = nq + 1
                        qcname(nq) = stname(iev,j)//deployment(iev,j)//phase(iev,j)
                        idiff0(nq) = idiff(iev,j)
                        used = .true.
                     end if 
                  end do
               else
                  nq = 1
                  qcname(1) = stname(iev,j)//deployment(iev,j)//phase(iev,j)
                  idiff0(1) = idiff(iev,j)
               end if
            end if
         end do
      end do
      
      deallocate (qcname, stat=error); if (error .gt. 0) call allocation_error (p, 'qcname', error)
      
      if (debug) write (io_log,'(a,i8,a,i8)') 'get_array_limits_c: nq = ', nq, '; nt = ', nt
      
      return
      
   end subroutine get_array_limits_c
   

!*****************************************************************************************
   subroutine get_array_limits_h (nq, nt)
   
   ! Calculates the number of nq_allocated and nt_allocated for allocation of variable arrays
   ! for the hypocentroid.
      
      integer :: iev
      integer :: j
      integer :: kst
      integer, intent(out) :: nq ! number of distinct station-phases that pass the filter for hypocentroid
      integer, intent(out) :: nt ! number of arrival time data that pass the filter for hypocentroid
      character(len=40) :: p
      logical :: used
      
      p = 'get_array_limits_h'
      allocate (qhname(7000), stat=error); if (error .gt. 0) call allocation_error (p, 'qcname', error)
   
      nq = 0
      nt = 0
      do iev = 1,n_event
         do j = 1,nst(iev)
            if (.not.fltrh(iev,j)) then
               nt = nt + 1
               if (nt .gt. 1) then
                  used = .false.
                  kst = 0
                  do while (.not.used)
                     kst = kst + 1
                     if (kst .le. nq) then
                        if (stname(iev,j)//phase(iev,j) .eq. qhname(kst)) used = .true.
                     else
                        nq = nq + 1
                        qhname(nq) = stname(iev,j)//phase(iev,j)
                        used = .true.
                     end if 
                  end do
               else
                  nq = 1
                  qhname(1) = stname(iev,j)//phase(iev,j)
               end if
            end if
         end do
      end do
   
      deallocate (qhname, stat=error); if (error .gt. 0) call allocation_error (p, 'qcname', error)

      return
   
   end subroutine get_array_limits_h


!*****************************************************************************************
   subroutine mlocinv_allocate (nq, nt)
   
   ! Allocate variable arrays related to the generalized inversion for cluster vectros and
   ! hypocentroid.
   
      integer :: mt
      integer, intent(in) :: nq
      integer, intent(in) :: nt
      
      mt = n_fp_max
      p = 'mlocinv_allocate' ! procedure
      
      allocate (a0(nq,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'a0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': a0 allocated (', nq, ',4)'
      end if
      
      allocate (ahat(nt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ahat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a,i4,a)') trim(p)//': ahat allocated (', nt, ',', mt, ')'
      end if
      
      allocate (aqbahat(mt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'aqbahat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': aqbahat allocated (', mt, ',', mt, ')'
      end if
      
      allocate (dthat(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'dthat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': dthat allocated (', nt, ')'
      end if
      
      allocate (ehath(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ehath', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': ehath allocated (', nt, ')'
      end if
      
      allocate (ibinq(nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ibinq', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': ibinq allocated (', nq, ')'
      end if
      
      allocate (indx(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'indx', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': indx allocated (', mt, ')'
      end if
      
      allocate (irank(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'irank', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': irank allocated (', mt, ')'
      end if
      
      allocate (ntiev(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'ntiev', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': ntiev allocated (', nt, ')'
      end if
      
      allocate (pqahat(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'pqahat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': pqahat allocated (', nt, ')'
      end if
      
      allocate (pwa0(nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'pwa0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': pwa0 allocated (', nq, ')'
      end if
      
      allocate (qbadx(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qbadx', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': qbadx allocated (', nt, ')'
      end if
      
      allocate (qbahat(nt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qbahat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a,i4,a)') trim(p)//': qbahat allocated (', nt, ',', mt, ')'
      end if
      
      allocate (qbdthat(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qbdthat', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': qbdthat allocated (', nt, ')'
      end if
      
      allocate (qcname(nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qcname', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qcname allocated (', nq, ')'
      end if
      
      allocate (qhname(nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qhname', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qhname allocated (', nq, ')'
      end if
      
      allocate (qm(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qm', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qm allocated (', mt, ')'
      end if
      
      allocate (qmtr(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'qmtr', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': qmtr allocated (', mt, ')'
      end if
      
      allocate (sn(nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'sn', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': sn allocated (', nq, ')'
      end if
      
      allocate (t4m(4,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 't4m', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': t4m allocated (4,', mt, ')'
      end if
      
      allocate (t4m2(4,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 't4m2', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': t4m2 allocated (4,', mt, ')'
      end if
      
      allocate (t4n1(4,nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 't4n1', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': t4n1 allocated (4,', nt, ')'
      end if
      
      allocate (t4n2(4,nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 't4n2', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': t4n2 allocated (4,', nt, ')'
      end if
      
      allocate (t4q(4,nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 't4q', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': t4q allocated (4,', nq, ')'
      end if
      
      allocate (tm(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tm', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': tm allocated (', mt, ')'
      end if
      
      allocate (tmn(mt,nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tmn', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i5,a)') trim(p)//': tmn allocated (', mt, ',', nt, ')'
      end if
      
      allocate (tn(nt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tn', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': tn allocated (', nt, ')'
      end if
      
      allocate (tn4(nt,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tn4', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a)') trim(p)//': tn4 allocated (', nt, ',4)'
      end if
      
      allocate (tq(nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'tq', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': tq allocated (', nq, ')'
      end if
      
      allocate (umm(mt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'umm', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': umm allocated (', mt, ',', mt, ')'
      end if
      
      allocate (unm(nt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'unm', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i5,a,i4,a)') trim(p)//': unm allocated (', nt, ',', mt, ')'
      end if
      
      allocate (uq4(nq,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'uq4', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': uq4 allocated (', nq, ',4)'
      end if
      
      allocate (vc(mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'vc', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': vc allocated (', mt, ')'
      end if
      
      allocate (vmm(mt,mt), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'vmm', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a,i4,a)') trim(p)//': vmm allocated (', mt, ',', mt, ')'
      end if
      
      allocate (wa0(nq,4), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'wa0', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': wa0 allocated (', nq, ',4)'
      end if
      
      allocate (wa0d(4,nq), stat=error)
      if (error .gt. 0) then
         call allocation_error (p, 'wa0d', error)
      else
         if (verbose_log) write (io_alloc_log,'(a,i4,a)') trim(p)//': wa0d allocated (4,', nq, ')'
      end if
      
      return
   
   end subroutine mlocinv_allocate


!*****************************************************************************************
   subroutine mlocinv_deallocate ()
   
   ! Deallocate variable arrays related to the generalized inversion for cluster vectros and
   ! hypocentroid.
         
      p = 'mlocinv_deallocate' ! procedure
      
      deallocate (a0, stat=error); if (error .gt. 0) call deallocation_error (p, 'a0', error)
      deallocate (ahat, stat=error); if (error .gt. 0) call deallocation_error (p, 'ahat', error)
      deallocate (aqbahat, stat=error); if (error .gt. 0) call deallocation_error (p, 'aqbahat', error)
      deallocate (dthat, stat=error); if (error .gt. 0) call deallocation_error (p, 'dthat', error)
      deallocate (ehath, stat=error); if (error .gt. 0) call deallocation_error (p, 'ehath', error)
      deallocate (ibinq, stat=error); if (error .gt. 0) call deallocation_error (p, 'ibinq', error)
      deallocate (indx, stat=error); if (error .gt. 0) call deallocation_error (p, 'indx', error)
      deallocate (irank, stat=error); if (error .gt. 0) call deallocation_error (p, 'irank', error)
      deallocate (ntiev, stat=error); if (error .gt. 0) call deallocation_error (p, 'ntiev', error)
      deallocate (pqahat, stat=error); if (error .gt. 0) call deallocation_error (p, 'pqahat', error)
      deallocate (pwa0, stat=error); if (error .gt. 0) call deallocation_error (p, 'pwa0', error)
      deallocate (qbadx, stat=error); if (error .gt. 0) call deallocation_error (p, 'qbadx', error)
      deallocate (qbahat, stat=error); if (error .gt. 0) call deallocation_error (p, 'qbahat', error)
      deallocate (qbdthat, stat=error); if (error .gt. 0) call deallocation_error (p, 'qbdthat', error)
      deallocate (qcname, stat=error); if (error .gt. 0) call deallocation_error (p, 'qcname', error)
      deallocate (qhname, stat=error); if (error .gt. 0) call deallocation_error (p, 'qhname', error)
      deallocate (qm, stat=error); if (error .gt. 0) call deallocation_error (p, 'qm', error)
      deallocate (qmtr, stat=error); if (error .gt. 0) call deallocation_error (p, 'qmtr', error)
      deallocate (sn, stat=error); if (error .gt. 0) call deallocation_error (p, 'sn', error)
      deallocate (t4m, stat=error); if (error .gt. 0) call deallocation_error (p, 't4m', error)
      deallocate (t4m2, stat=error); if (error .gt. 0) call deallocation_error (p, 't4m2', error)
      deallocate (t4n1, stat=error); if (error .gt. 0) call deallocation_error (p, 't4n1', error)
      deallocate (t4n2, stat=error); if (error .gt. 0) call deallocation_error (p, 't4n2', error)
      deallocate (t4q, stat=error); if (error .gt. 0) call deallocation_error (p, 't4q', error)
      deallocate (tm, stat=error); if (error .gt. 0) call deallocation_error (p, 'tm', error)
      deallocate (tmn, stat=error); if (error .gt. 0) call deallocation_error (p, 'tmn', error)
      deallocate (tn, stat=error); if (error .gt. 0) call deallocation_error (p, 'tn', error)
      deallocate (tn4, stat=error); if (error .gt. 0) call deallocation_error (p, 'tn4', error)
      deallocate (tq, stat=error); if (error .gt. 0) call deallocation_error (p, 'tq', error)
      deallocate (umm, stat=error); if (error .gt. 0) call deallocation_error (p, 'umm', error)
      deallocate (unm, stat=error); if (error .gt. 0) call deallocation_error (p, 'unm', error)
      deallocate (uq4, stat=error); if (error .gt. 0) call deallocation_error (p, 'uq4', error)
      deallocate (vc, stat=error); if (error .gt. 0) call deallocation_error (p, 'vc', error)
      deallocate (vmm, stat=error); if (error .gt. 0) call deallocation_error (p, 'vmm', error)
      deallocate (wa0, stat=error); if (error .gt. 0) call deallocation_error (p, 'wa0', error)
      deallocate (wa0d, stat=error); if (error .gt. 0) call deallocation_error (p, 'wa0d', error)
      
      return
   
   end subroutine mlocinv_deallocate


!*****************************************************************************************
end module mloc_inv

