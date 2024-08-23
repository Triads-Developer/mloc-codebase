!> Declare certain limits for the relocation

module declare_limits

   implicit none
   save
   
   ! Parameters
   integer, parameter :: bp_nlat = 2161 ! Number of latitude values of bouncepoint topography data file
   integer, parameter :: bp_nlon = 4321 ! Number of longitude values of bouncepoint topography data file
   integer, parameter :: itmax = 4 ! Maximum number of iterations
   integer, parameter :: n_bdp_max = 60 ! Maximum number of stations in the 'bad depth phase' list
   integer, parameter :: n_cell_max = 100 ! Maximum number of cells in the grid for subdivided type 3 tomographic output
   integer, parameter :: n_ellipse_nmax = 20 ! Maximum number of individual ellipses that can plotted
   integer, parameter :: n_epa_plot_max = 6 ! Maximum number of empirical path anomaly plots
   integer, parameter :: n_fault_map_nmax = 4 ! Maximum number of digital fault map files that can be read
   integer, parameter :: n_hrdp_max = 760 ! Maximum depth for testing and plotting depth from relative depth phases
   integer, parameter :: n_junk_max = 100 ! Maximum number of entries in the junk phase list
   integer, parameter :: n_miss_sta_max = 400 ! Maximum number of missing stations for any one event or all events together
   integer, parameter :: n_no_phreid_max = 6 ! Maximum number of phases that can be prevented from phase re-identification
   integer, parameter :: n_phrp_max = 10 ! Maximum number of phase residual plots
   integer, parameter :: n_psdre_max = 20 ! Maximum number of phase-specific default reading error values
   integer, parameter :: n_qres_max = 400 ! Maximum number of samples allowed for estimating empirical reading error
   integer, parameter :: n_radf_max = 40 ! Maximum number of stations for which agency and deployment fields can be read
   integer, parameter :: n_skip_max = 10 ! Maximum number of stations or phases that can be skipped
   integer, parameter :: n_star_nmax = 30 ! Maximum number of individual stars that can plotted
   integer, parameter :: n_stat_nmax = 30 ! Maximum number of station locations that can plotted
   integer, parameter :: n_supp_stn_file_max = 8 ! Maximum number of supplemental station files that can be read
   integer, parameter :: n_terr_max = 4 ! Maximum number of stations for which a time correction can be made
   integer, parameter :: n_tt5s_max = 20 ! Maximum number of single-station tt5 plots
   integer, parameter :: n_ttou_max = 6 ! Maximum number of output files with TT data for specific phases
   integer, parameter :: n_xsec_max = 4 ! Maximum number of cross sections that can plotted
   integer, parameter :: ndiffmax = 500 ! Maximum number of differential times
   integer, parameter :: nitomomax = 10 ! Maximum number of tomography output files
   integer, parameter :: nmax1 = 30000 ! Maximum number of stations in the station list (plus any auxilliary list)
   integer, parameter :: nplot_max = 9 ! Maximum number of selected-event plots
   integer, parameter :: nsprdmax = 60 ! Maximum number of entries in the ttsprd file
   integer, parameter :: n_ifsm_matched_max = 1000 ! Maximum number of matched entries
   integer, parameter :: n_nsmd_matched_max = 2000 ! Maximum number of matched entries
   
   ! Limits updated during runtime
   integer :: nevmax =200 ! number of events in the current cluster
   integer :: npmax = 7000 ! maximum number of unfiltered phase records in any event of the cluster (proxy for ntmax0)
   integer :: nqmax = 7000 ! max number of distinct station-phases used
   
   ! Derived limits, based on nevmax, updated during runtime
   integer :: mtmax ! max number of free parameters
   integer :: mtmax2
   integer :: nev ! for backward compatibility in the code
   integer :: nhdfmax ! Maximum number of events in an HDF file for input as starting locations.
   integer :: nmax
   integer :: nmax2
   integer :: ntmax1 ! max number TT residuals used in inversion
   integer :: ntqmax ! max number of instances of a particular station-phase
   
   ! Derived limits, based on itmax
   integer, parameter :: itmax1 = itmax + 1
   integer, parameter :: itmax2 = itmax + 2
   
end module declare_limits
