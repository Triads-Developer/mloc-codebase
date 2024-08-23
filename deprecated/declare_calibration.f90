!> Declaration and allocation of variables related to calibration

module declare_calibration

   use declare_limits
   use mloclib_messages
   
   implicit none
   save
   
   ! Calibration type
   integer :: icaltype
   logical :: direct_cal
   logical :: indirect_cal
   
   ! Direct calibration
   real, allocatable, dimension(:) :: alphadc
   real, allocatable, dimension(:) :: ddepdc
   real, allocatable, dimension(:) :: dotdc
   real, allocatable, dimension(:) :: xl1dc
   real, allocatable, dimension(:) :: xl2dc
   character(len=3) :: dcal_pr
   
   ! Indirect calibration
   integer, allocatable, dimension(:,:) :: cal_hr
   integer, allocatable, dimension(:,:) :: cal_min
   integer, dimension(3) :: ncal
   real, allocatable, dimension(:,:,:) :: accv
   real :: alphagt
   real, allocatable, dimension(:,:) :: azes_cal
   real, allocatable, dimension(:,:) :: azse_cal
   real, allocatable, dimension(:,:) :: cal_dep
   real, allocatable, dimension(:,:) :: cal_lat
   real :: cal_level
   real, allocatable, dimension(:,:) :: cal_lon
   real, allocatable, dimension(:,:) :: cal_sec
   real :: del_cal_dep
   real :: del_cal_lat
   real :: del_cal_lon
   real :: del_cal_tim
   real, allocatable, dimension(:,:) :: delt_cal
   real :: depgt
   real :: depthh_epa
   real, allocatable, dimension(:) :: depthp_cal
   real, allocatable, dimension(:,:) :: dt_cal
   real, allocatable, dimension(:,:) :: elcr_cal
   real :: gtlevel
   real :: lath_epa
   real, allocatable, dimension(:) :: latp_cal
   real :: lonh_epa
   real, allocatable, dimension(:) :: lonp_cal
   real :: otgt
   real, allocatable, dimension(:) :: otsp_cal
   real, allocatable, dimension(:) :: otsr
   real, allocatable, dimension(:,:,:,:) :: rcv
   real :: rdbt
   real, allocatable, dimension(:,:) :: s_cal
   real :: s1gt
   real :: s2gt
   real, allocatable, dimension(:,:,:) :: scv
   real, allocatable, dimension(:,:) :: sdstcs_cal
   real, allocatable, dimension(:,:,:) :: sgcv
   real, allocatable, dimension(:,:) :: stacs_cal
   real, allocatable, dimension(:) :: w12
   real, allocatable, dimension(:) :: w3
   real, allocatable, dimension(:) :: w4
   logical, allocatable, dimension(:,:) :: cal_event
   logical :: ot_cal
   character(len=4), allocatable, dimension(:) :: cal_code_input
   character(len=1), allocatable, dimension(:) :: cal_par
   character(len=12) :: icaltype_pr

contains

!*****************************************************************************************
   subroutine declare_calibration_allocate ()
      
      integer :: error
      integer :: m
      integer :: n
      character(len=32) :: p
      
      n = nevmax
      m = npmax
      p = 'calibration_allocate'  ! procedure name
   
      allocate (accv(n,4,4), stat=error); if (error .gt. 0) call allocation_error (p, 'accv', error)
      allocate (alphadc(n), stat=error); if (error .gt. 0) call allocation_error (p, 'alphadc', error)
      allocate (azes_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'azes_cal', error)
      allocate (azse_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'azse_cal', error)
      allocate (cal_code_input(n), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_code_input', error)
      allocate (cal_dep(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_dep', error)
      allocate (cal_event(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_event', error)
      allocate (cal_hr(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_hr', error)
      allocate (cal_lat(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_lat', error)
      allocate (cal_lon(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_lon', error)
      allocate (cal_min(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_min', error)
      allocate (cal_par(n), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_par', error)
      allocate (cal_sec(n,3), stat=error); if (error .gt. 0) call allocation_error (p, 'cal_sec', error)
      allocate (ddepdc(n), stat=error); if (error .gt. 0) call allocation_error (p, 'ddepdc', error)
      allocate (delt_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'delt_cal', error)
      allocate (depthp_cal(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depthp_cal', error)
      allocate (dotdc(n), stat=error); if (error .gt. 0) call allocation_error (p, 'dotdc', error)
      allocate (dt_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'dt_cal', error)
      allocate (elcr_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'elcr', error)
      allocate (latp_cal(n), stat=error); if (error .gt. 0) call allocation_error (p, 'latp_cal', error)
      allocate (lonp_cal(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lonp_cal', error)
      allocate (otsp_cal(n), stat=error); if (error .gt. 0) call allocation_error (p, 'otsp_cal', error)
      allocate (otsr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'otsr', error)
      allocate (rcv(n,3,4,4), stat=error); if (error .gt. 0) call allocation_error (p, 'rcv', error)
      allocate (s_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 's_cal', error)
      allocate (scv(n,4,4), stat=error); if (error .gt. 0) call allocation_error (p, 'scv', error)
      allocate (sdstcs_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'sdstcs_cal', error)
      allocate (sgcv(n,4,4), stat=error); if (error .gt. 0) call allocation_error (p, 'sgcv', error)
      allocate (stacs_cal(n,m), stat=error); if (error .gt. 0) call allocation_error (p, 'stacs_cal', error)
      allocate (w12(n), stat=error); if (error .gt. 0) call allocation_error (p, 'w12', error)
      allocate (w3(n), stat=error); if (error .gt. 0) call allocation_error (p, 'w3', error)
      allocate (w4(n), stat=error); if (error .gt. 0) call allocation_error (p, 'w4', error)
      allocate (xl1dc(n), stat=error); if (error .gt. 0) call allocation_error (p, 'xl1dc', error)
      allocate (xl2dc(n), stat=error); if (error .gt. 0) call allocation_error (p, 'xl2dc', error)
!     allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
   
      return
      
   end subroutine declare_calibration_allocate
   
   
!*****************************************************************************************
end module declare_calibration
