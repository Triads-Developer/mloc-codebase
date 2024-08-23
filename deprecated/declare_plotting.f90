!> Declaration of variables related to plotting

module declare_plotting

   use declare_limits
   use mloclib_messages
   
   implicit none
   save
   
   ! Plot stars
   integer :: n_star
   integer, dimension(n_star_nmax) :: iev_star
   real, dimension(n_star_nmax) :: star_size
   logical :: star_plot
   
   ! Plot stations
   integer :: n_stat
   real, dimension(n_stat_nmax,3) :: stat_data
   logical :: stat_plot
   
   ! Plot ellipses
   integer :: n_ellipse
   real, dimension(n_ellipse_nmax,5) :: ellipse_data
   logical :: ellipse_plot
   
   ! Digital fault map
   integer :: n_fault_map
   character(len=100), dimension(n_fault_map_nmax) :: fault_map_filename
   character(len=30) :: fault_plot_style
   logical :: fault_map
   
   ! Plot topography in GMT
   character(len=132) :: dem2_filename
   logical :: plot_dem1
   logical :: plot_dem2
   
   ! Cross sections
   integer :: n_xsec
   real, dimension(n_xsec_max) :: xsec_az
   
   ! Single-station tt5 plots
   integer :: n_tt5s
   character(len=5), dimension(n_tt5s_max) :: tt5s_sta
   logical :: tt5s

   ! Single-event tt5 plots
   integer :: n_tt5e
   character(len=30), allocatable, dimension(:) :: tt5e_evt
   logical :: tt5e
   logical :: tt5e_all
   
   ! Phase residual plots
   integer :: n_phrp 
   real, dimension(n_phrp_max,2) :: phrp_delt
   character(len=8), dimension(n_phrp_max) :: phrp_phase
   logical :: phrp

   ! Empirical path anomaly plots
   integer :: n_epa_plot
   real, dimension(n_epa_plot_max) :: epa_plot_distance
   character(len=8), dimension(n_epa_plot_max) :: epa_plot_phase
   logical :: epa_plot
   logical, dimension(n_epa_plot_max) :: epa_plot_raypath

   ! Residual plots
   integer :: n_res_plot
   real, allocatable, dimension(:) :: res_plot_delta_min
   real, allocatable, dimension(:) :: res_plot_delta_max
   real :: res_plot_delta_min_all
   real :: res_plot_delta_max_all
   character(len=30), allocatable, dimension(:) :: res_plot_event
   logical :: resp
   logical :: resp_all
    
   ! Plot types and options
   logical :: eplt
   logical :: fdhp
   logical :: gmtp
   logical, allocatable, dimension(:,:) :: plot
   logical :: reduced
   logical :: rpth
   logical :: splt
   logical :: tt1
   logical :: tt2
   logical :: tt3
   logical :: tt4
   logical :: tt5
   logical :: tt6
   logical :: tt7
   logical :: tt8
   logical :: tt9
   logical, dimension(4) :: vectors
   
contains

!*****************************************************************************************
   subroutine plotting_allocate ()
   
   ! Allocate variable arrays related to seismic phases
   
      integer :: error
      integer :: m, n
      character(len=32) :: p
      
      m = nplot_max
      n = nevmax
      p = 'plotting_allocate' ! procedure
      
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (plot(n,0:m), stat=error); if (error .gt. 0) call allocation_error (p, 'plot', error)
      allocate (res_plot_delta_max(n), stat=error); if (error .gt. 0) call allocation_error (p, 'res_plot_delta_max', error)
      allocate (res_plot_delta_min(n), stat=error); if (error .gt. 0) call allocation_error (p, 'res_plot_delta_min', error)
      allocate (res_plot_event(n), stat=error); if (error .gt. 0) call allocation_error (p, 'res_plot_event', error)
      allocate (tt5e_evt(n), stat=error); if (error .gt. 0) call allocation_error (p, 'tt5e_evt', error)
      
      return
   
   end subroutine plotting_allocate


!*****************************************************************************************
end module declare_plotting

