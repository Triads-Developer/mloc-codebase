module mloc_init

! Initialization of some global variables declared in mloc_declare.

   use mloc_declare
   
   implicit none
   save
   
   
   contains

!*****************************************************************************************
subroutine initialize ()

   ! Initialize many of the variables defined in mloc_declare
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!   Environment

   basename = 'mloc' ! name for the current run
   datadir = ' ' ! directory containing the data for this run
   dirsym = '/' ! symbol used to separate directories
   gmt_version = 6
   ishell = 1 ! shell ( 1 for bash, 2 for csh or 3 for zsh) used in GMT scripts
   mloc_author = 'default ' ! Author
   mloc_status = 'beta'
   mloc_version = 'mloc v11.0.3, release date Aug 17, 2024'
   
   ! File paths
   cpt_path = 'tables/gmt/cpt'
   dem_path = 'tables/gmt/dem'
   ellip_path = 'tables/ellipticity'
   gmt_script_dir = ' '
   mloc_path = ' '
   star_auto_path = 'tables/gmt'
   station_path = 'tables/stn'
   taup_path = 'tables/tau-p'

   ! File names
   cpt_file = 'topo.cpt'
   outfile = 'default' ! basename with path
   star_auto_file = 'star_mag_size.dat'
   station_master = 'master_stn.dat'

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Configuration

   ! Logging and verbosity
   debug = .false.
   verbose_log = .false.
   verbose_screen = .false.

   ! Longitude basis
   longitude_range = 0 ! Longitudes will be converted to the range -180 to 180°

   ! Windowing residuals
   wind1 = 3.
   wind2 = 4.
   windloclim = 1.2 ! Windows are expanded at epicentral distances less than this

   ! Epicentral distance ranges for hypocentroid and cluster vectors
   clim(1,1) = 1.0
   clim(1,2) = 180.0
   hlim(1,1) = 30.
   hlim(1,2) = 90.
   nclim = 1 ! Number of distance windows in use for cluster vectors
   nhlim = 1 ! Number of distance windows in use for hypocentroid

   ! Data weighting and bias correction (in the Jordan & Sverdrup sense)
   bias_corr = .true. ! Bias correction of hypocentroid from events with many more readings
   data_weight = .true. ! Weighting data by inverse of error
   pttt = .false. ! perfect theoretical travel times assumed for hypocentroid

   ! Test mode
   test_mode = .false. ! activation of alternative code segements
   
   ! Selection of a direct calibration subcluster
   subc_delt = 0.0 ! Default distance for selection of subclusters
   subc_nconnect = 999 ! minimum number of connected readings within default distance for selection of subclusters
   subc_nmin = 999 ! minimum number of readings within default distance for selection of subclusters
   subc_set = .false. ! Criteria for selection of subclusters
   
   ! Convergence limits
   cl_dep_c = 0.5 ! Cluster vector depth, km
   cl_dep_h = 0.5 ! Hypocentroid depth, km
   cl_epi_c = 0.5 ! Cluster vector lat-lon, km
   cl_epi_h = 0.005 ! Hypocentroid lat-lon, degrees
   cl_ot_c = 0.1 ! Cluster vector OT, sec
   cl_ot_h = 0.1 ! Hypocentroid OT, sec
   convergence_test_index = 1
   
   ! Kill events
   kill_all = .false.
   kill_one = .false.
   kill_reason = ' '
   
   ! Skip specified stations, phases and authors
   n_skip = 0 ! Number of SKIP commands
   skip = .false.
   skip_params = ' '
      
   ! Take starting locations from HDF file of a previous run
   nrhdf = 0
   read_hdf = .false.
   rhdf_filnam = ' '

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
!  Inversion

   damping = .false. ! Subtracts average from cluster vector changes
   no_xflags = .false. ! Ignore 'x' flags
   nstep = n_it_max ! Number of iterations to perform
   tikhonov_factor = 0. ! Tikhonov regularization

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Command arguments

   cmndfil = .true. ! reading from a command_file
   params = ' '

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Hypocentroid

   depthfh = .true. ! depth as a free parameter for the hypocentroid
   fixprh = 'free '
   latfh = .true. ! latitude as a free parameter for the hypocentroid
   lonfh = .true. ! longitude as a free parameter for the hypocentroid
   ponly = .true. ! Use P only or all phase types for hypocentroid
   timefh = .true. ! origin time as a free parameter for the hypocentroid
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Cluster Vectors
   radius_cvff = 0. ! cluster vector fudge factor, km

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Events

   ! Depth-related
   depth_default = -999.
   depth_default_c = 'c' ! Depth code used for default cluster depth
   depth_default_minus = 99.
   depth_default_plus = 99.
   median_constrained_depths = -99.
   n_cfd = 0
   rescfd = -99.

   ! Ellipse-related
   alphah = 0.
   xl1h = 0.
   xl2h = 0.
   
   ! Hypocenters, magnitudes
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Phase data

   dflag = .true. ! use data flags (x, p, d, k, l)   
   
   ! Timing error corrections
   n_terr = 0
   terr_corr = 0.
   terr_stn = ' '
   timing_error_correction = .false. ! Timing correction
   
   ! Phase-specific default reading error values
   n_psdre = 0
   
   ! Relative depth phases
   h_rdp1 = 1
   h_rdp2 = n_hrdp_max - 1
   n_rdpp = 0
   rdpp = .false. ! Relative depth phase plot for a single event
   rdpp_all = .false. ! Relative depth phase plot for all events with data
   
   ! Differential time data
   diffdat = .false.
   ndiff = 0

   ! Reading errors
   nqname = 0 ! number of station-phase reading errors read from ~.rderr file
   rderr_loc_delt = 0.6 ! Distance limit for applying comand RELS
   rderr_loc_p = 0.25 ! Reading error for P at local stations (RELS)
   rderr_loc_s = 0.40 ! Reading error for S at local stations (RELS)
   rderr_min = 0.15 ! Minimum value allowed for reading errors, for phases beyond local distance, except depth phases (MARE)
   rderr_min_depth = 1.0 ! Minimum value allowed for reading errors of teleseismic depth phases (MARE)
   rderr_min_loc = 0.10 ! Minimum value allowed for reading errors at local distance (MARE)
   rderrfname = 'default'
   read_rderr = .false. ! Use a file of empirical reading errors for individual station-phases
   rels_set = .false. ! Status of setting of reading errors for local stations (RELS)

   ! Junk Phases
   junk_phase = ' ' ! phases to ignore
   n_junk = 0

   ! Phase re-identification
   phid = .true. ! phase identification
   n_no_phreid = 0
   plogout = .false. ! ~.plog file output, phase identification debugging

   ! Spread and baseline offset of specific phases
   nsprd = 0 ! number of phase spread entries read from a ~.ttsprd file
   read_ttsprd = .false. ! Use a file of empirical travel-time spreads for different phases
   ttsprdfname = 'default'

   ! Crustal heterogeneity factor for differential-time data
   cscale = 10. ! Spatial scale over which crustal velocities may vary
   vmr = 0.05 ! Percentage variation in crustal velocities

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Stations

   n_failed_date_range = 0
   n_radf = 0
   n_supp_stn = 0 ! total number of stations read from supplementary station files
   n_supp_stn_file = 0 ! Number of supplementary station files
   nkstat = 0
   nstat1 = 0 ! total number of stations read (master + supplementary)
   nstn_used = 0 ! number of stations used in direct calibration
   read_ad = .false. ! Read agency and deployment fields to resolve station code conflicts
   suppfilnam = ' '
   suppstn = .false. ! Use upplementary station code files

   ! List of stations suspected of reporting bogus depth phases
   bdp_filnam = ' '
   bdp_list = .false. ! Use list of stations known to report bad depth phases
   bdp_station = ' '
   n_bdp = 0 ! Number of entries in the 'bad depth phases' list

   ! Duplicate entries in station lists
   duplication = ' '
   n_dupe = 0
   n_dupe_conflict = 0
   n_dupe_minor = 0
   n_dupe_pure = 0
   n_dupe_significant = 0

   ! Missing station codes and other station statistics
   n_miss_sta_list = ' '

   ! NEIC station metadata
   nsmd = .false. ! Search of NEIC station metadata for missing station codes

   ! ISF-FDSN station metadata
   ifsm = .false. ! Search of ISC-FDSN station metadata for missing station codes

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Travel time

   ! Local velocity model
   dlimlocmod = 0. ! Epicentral distance limit for local model (LMOD)
   zlimlocmod = 0. ! Focal depth limit for local model (LMOD)
   locmod = .false. ! Use a local velocity model (LMOD)
   locmodfname = ' ' ! Pathname to local crustal model file (LMOD)

   ! Phases not in ak135
   lg_a = 0. ! Intercept, in sec (LGTT)
   lg_b = 31.5 ! Slope, in sec/degree (LGTT)
   lg_min = 2.5 ! Minimum epicentral distance (degrees) at which Lg travel times will be calculated (LGTT)
   p__a = 0. ! Intercept, in sec (P_TT)
   p__b = 19.0 ! Slope, in sec/degree (P_TT)
   p__min = 2.5 ! Minimum epicentral distance (degrees) at which P_ travel times will be calculated (P_TT)
   rg_a = 0. ! Intercept, in sec (RGTT)
   rg_b = 36.5 ! Slope, in sec/degree (RGTT)
   rg_min = 2.5 ! Minimum epicentral distance (degrees) at which Rg travel times will be calculated (RGTT)
   tphase_a = 15. ! Intercept, in sec (TPTT)
   tphase_b = 75. ! Slope, in sec/degree (TPTT)

   ! Station elevation corrections
   pcrvel = 5.8 ! P velocity for station elevation correction
   scrvel = 3.46 ! S velocity for station elevation correction
   tt_corr = 1 ! Station elevation correction (CORR)

   ! tau-P
   log_tt = .false. ! Log file for intermediate tau-p results
   prnt = .true. ! common /prtflc/
   segmsk = .true. ! (common /prtflc/)
   taup_model = 'ak135' ! (TAUP)

   ! Main phase list (mpl), travel times and derivatives

   ! Bounce point topography
   bptc = .true. ! Bounce point topography correction, for pP, sP and pwP

   ! Empirical TT output files
   n_ttou = 0 ! Number of phase-specific empirical TT output files
   ttou = .false. ! Output of empirical TTs for specific phases
   ttou_phase = ' ' ! List of phases for empirical TT output files

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calibration

   ! Calibration type
   direct_cal = .false.
   icaltype = 2
   indirect_cal = .false.

   ! Direct calibration
   dcal_pr = 'off'

   ! Indirect calibration
   cal_level = 99
   icaltype_pr = 'systematic  '
   ncal = 0
   rdbt = 0. ! radius of doubt

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Output files

   ! Tomography output files
   nitomo = 0 ! Number of TOMO commands given
   tomo_sub = .false. ! if subdividing is requested

   ! Print variables
   data_weight_pr = 'on '
   dflag_pr = 'on '
   pttt_pr = 'off'

   ! Output file for large residuals
   lres = 3.0 ! Threshold for large cluster residuals
   lresout = .true. ! .lres file for large cluster residuals

   ! Output file to pass for input to BAYESLOC
   blocout = .false. ! data file for import into BAYESLOC
   
   ! .datf file
   datfout = .false. ! ~.datf file of flagged, re-identified arrival time data.

   ! GCCEL files
   commentary_buffer = ' '
   commentary_fname = ' '
   gcat_folder = ' '
   gccelout = .false. ! data file for GCCEL

   ! Output file for limited distance range (usually used for Pg/Pn cross-over range)
   dist1 = 0. ! Minimum distance for oldr
   dist2 = 180. ! Max distance for oldr
   oldr_out = .false. ! Output for limited distance range (oldr)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Plotting

   ! Plot stars
   star_auto = .false.
   star_plot = .false.

   ! Plot stations
   n_stat = 0
   stat_plot = .false.

   ! Plot ellipses
   ellipse_plot = .false.
   n_ellipse = 0

   ! Digital fault map
   fault_map = .false.
   fault_map_filename = ' '
   fault_plot_style = '-Wthick,red'
   n_fault_map = 0

   ! Plot topography in GMT
   dem2_filename = ' '
   plot_dem1 = .false.
   plot_dem2 = .false.

   ! Cross sections
   n_xsec = 0

   ! Single-station tt5 plots
   n_tt5s = 0
   tt5s = .false.
   tt5s_sta = ' '
   
   ! Single-event tt5 plots
   n_tt5e = 0
   tt5e = .false.
   tt5e_all = .false.

   ! Phase residual plots
   n_phrp = 0
   phrp = .false.
   phrp_delt = 0.
   phrp_phase = ' '

   ! Empirical path anomaly plots
   epa_plot = .false.
   n_epa_plot = 0

   ! Event residual plots
   n_evrp = 0
   evrp_delta_max_all = 180.
   evrp_delta_min_all = 0.
   evrp = .false. ! Residual plots for individual events
   evrp_all = .false. ! Residual plots for all events

   ! Plot types and options
   eplt = .false. ! confidence ellipse plot
   fdhp = .false. ! Focal depth histogram plot
   gmtp = .true. ! Enable GMT plotting
   reduced = .true. ! Use reduced velocities in type 6 and 7 TT plots
   rpth = .true. ! Plot direct calibration raypaths
   splt = .false. ! seismicity plot
   tt1 = .false.
   tt2 = .false.
   tt3 = .false.
   tt4 = .false.
   tt5 = .false.
   tt6 = .false.
   tt7 = .false.
   tt8 = .false.
   tt9 = .false.
   vectors = .true.

   ! Plot slip data (PSLP)
   n_pslp = 0
   pslp = .false.
   
   return
   
end subroutine initialize
   
   
end module mloc_init


