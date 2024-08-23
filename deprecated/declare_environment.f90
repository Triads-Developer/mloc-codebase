!> Declaration of variables relating to the mloc run environment that are used in multiple
! procedures. Many file names and path names are declared and initialized in the main
! program because they are only used there.

module declare_environment

   implicit none
   save
   
   character(len=50) :: mloc_version = 'mloc v11.0.1, release date May 5, 2024'

   integer :: gmt_version = 6
   integer :: ishell = 2 ! shell (1 for csh, 2 for bash or 3 for zsh) used in GMT scripts

   character(len=32) :: basename = 'mloc' ! name for the current run
   character(len=80) :: datadir = ' ' ! directory containing the data for this run
   character(len=1) :: dirsym = '/' ! symbol used to separate directories
   character(len=8) :: mloc_author = 'default ' ! Author
   
   ! Folder paths
   character(len=32) :: cpt_path = 'tables/gmt/cpt'
   character(len=32) :: dem_path = 'tables/gmt/dem'
   character(len=32) :: ellip_path = 'tables/ellipticity'
   character(len=132) :: gmt_script_dir
   character(len=32) :: mloc_path = ' '
   character(len=32) :: station_path = 'tables/stn'
   character(len=32) :: taup_path = 'tables/tau-p'

   ! File names
   character(len=32) :: cpt_file = 'topo.cpt'
   character(len=100) :: outfile= 'default' ! Base-name for all output files
   character(len=32) :: station_master = 'master_stn.dat'
   
end module declare_environment

