!> Declaration of global variables related to output files

module declare_output

   use declare_limits
   use mloclib_messages
   
   implicit none
   save

   ! Tomography output files
   integer, dimension(nitomomax) :: itomo ! index of type of output
   integer, dimension(nitomomax) :: n_cell_lat ! # of cells in latitude
   integer, dimension(nitomomax) :: n_cell_lon ! # of cells in longitude
   integer :: nitomo ! number of tomo output files
   real, dimension(nitomomax) :: cell_lat_inc ! cell dimension for latitude in degrees
   real, dimension(nitomomax) :: cell_lon_inc ! cell dimension for longitude in degrees
   real, dimension(nitomomax) :: grid_origin_lat ! latitude of grid origin
   real, dimension(nitomomax) :: grid_origin_lon ! longitude of grid_origin
   character(len=8), dimension(nitomomax) :: tomo_phase ! seismic phase name
   logical, dimension(nitomomax) :: tomo_sub ! if subdividing is requested   
   
   character(len=20), allocatable, dimension(:) :: annotation
   character(len=64), allocatable, dimension(:) :: hypocenter_list
   character(len=180), allocatable, dimension(:) :: hdfline

   ! Print variables
   character(len=3), dimension(0:1) :: biaspr
   character(len=3) :: data_weight_pr
   character(len=5), allocatable, dimension(:) :: depthpr
   character(len=3) :: dflag_pr
   character(len=3), allocatable, dimension(:) :: latpr
   character(len=3), allocatable, dimension(:) :: lonpr
   character(len=3) :: pttt_pr
   character(len=80), dimension(2) :: stacpr
   character(len=20), dimension(2) :: tablpr
   character(len=11), allocatable, dimension(:) :: timepr
   character(len=3), dimension(0:1) :: weigpr

   ! Output file for large residuals
   real :: lres ! cutoff value for cluster residual.
   logical :: lresout ! .true. if LRES file should be created.

   ! Output file to pass for input to BAYESLOC
   real, allocatable, dimension(:) :: eqr
   logical :: blocout ! .true. if BAYESLOC file should be created.
   
   ! .datf file
   logical :: datfout ! .true. if .datf file should be created.
   
   ! GCCEL files
   character(len=80), dimension(9) :: commentary_buffer
   character(len=80) :: commentary_fname
   character(len=132) :: gcat_folder
   logical :: gccelout
   
   ! Output file for limited distance range (usually used for Pg/Pn cross-over range)
   real :: dist1
   real :: dist2
   logical :: oldr_out
   
contains
   
!*****************************************************************************************
   subroutine output_allocate ()
   
   ! Allocate variable arrays related to output files
   
      integer :: error
      integer :: n
      character(len=32) :: p
      
      n = nevmax
      p = 'output_allocate' ! procedure
      
!       allocate (, stat=error); if (error .gt. 0) call allocation_error (p, '', error)
      allocate (annotation(n), stat=error); if (error .gt. 0) call allocation_error (p, 'annotation', error)
      allocate (depthpr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'depthpr', error)
      allocate (eqr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'eqr', error)
      allocate (hdfline(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hdfline', error)
      allocate (hypocenter_list(n), stat=error); if (error .gt. 0) call allocation_error (p, 'hypocenter_list', error)
      allocate (latpr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'latpr', error)
      allocate (lonpr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'lonpr', error)
      allocate (timepr(n), stat=error); if (error .gt. 0) call allocation_error (p, 'timepr', error)

      return
   
   end subroutine output_allocate
   

!*****************************************************************************************
end module declare_output

