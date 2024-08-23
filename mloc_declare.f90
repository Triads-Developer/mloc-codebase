!> Declaration of most of the variables that used to be in the "mloc.inc" include file.
!> 'Simple' variables are initialized here. Allocatable arrays are initialized in module
!> 'mloc_init' after their dimensions have been calculated and they've been allocated.
!> Also carries the procedures used for messaging in mloc and the procedure for selecting
!> an unused logical unit number.

module mloc_declare

   implicit none
   save
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Program limits
   
   ! Parameters
   integer, parameter :: n_bdp_max = 60 ! Maximum number of stations in the 'bad depth phase' list
   integer, parameter :: n_bp_lat = 2161 ! Number of latitude values of bouncepoint topography data file
   integer, parameter :: n_bp_lon = 4321 ! Number of longitude values of bouncepoint topography data file
   integer, parameter :: n_cell_max = 100 ! Maximum number of cells in the grid for subdivided type 3 tomographic output
   integer, parameter :: n_diff_max = 500 ! Maximum number of differential times
   integer, parameter :: n_ellipse_nmax = 20 ! Maximum number of individual ellipses that can plotted
   integer, parameter :: n_epa_max = 6 ! Maximum number of empirical path anomaly plots
   integer, parameter :: n_fault_map_max = 4 ! Maximum number of digital fault map files that can be read
   integer, parameter :: n_hrdp_max = 760 ! Maximum depth for testing and plotting depth from relative depth phases
   integer, parameter :: n_ifsm_max = 1000 ! Maximum number of matched entries
   integer, parameter :: n_it_max = 4 ! Maximum number of iterations
   integer, parameter :: n_junk_max = 100 ! Maximum number of entries in the junk phase list
   integer, parameter :: n_miss_sta_max = 400 ! Maximum number of missing stations for any one event
    ! or total number of unique missing stations over all events
   integer, parameter :: n_no_phreid_max = 6 ! Maximum number of phases that can be prevented from phase re-identification
   integer, parameter :: n_nsmd_max = 2000 ! Maximum number of matched entries
   integer, parameter :: n_phrp_max = 10 ! Maximum number of phase residual plots
   integer, parameter :: n_psdre_max = 20 ! Maximum number of phase-specific default reading error values
   integer, parameter :: n_q_max = 7000 ! max number of distinct station-phases used
   integer, parameter :: n_qres_max = 400 ! Maximum number of samples allowed for estimating empirical reading error
   integer, parameter :: n_radf_max = 40 ! Maximum number of stations for which agency and deployment fields can be read
   integer, parameter :: n_selected_max = 9 ! Maximum number of selected-event plots
   integer, parameter :: n_skip_max = 10 ! Maximum number of stations or phases that can be skipped
   integer, parameter :: n_sprd_max = 60 ! Maximum number of entries in the ttsprd file
   integer, parameter :: n_sstn_max = 8 ! Maximum number of supplemental station files that can be read
   integer, parameter :: n_sta_max = 30000 ! Maximum number of stations in the station list (plus any auxilliary list)
   integer, parameter :: n_star_table = 5 ! Number of entries in the magnitude-size table for plotting stars
   integer, parameter :: n_stat_max = 30 ! Maximum number of station locations that can plotted
   integer, parameter :: n_terr_max = 4 ! Maximum number of stations for which a time correction can be made
   integer, parameter :: n_tomo_max = 10 ! Maximum number of tomography output files
   integer, parameter :: n_tt5s_max = 20 ! Maximum number of single-station tt5 plots
   integer, parameter :: n_ttou_max = 6 ! Maximum number of output files with TT data for specific phases
   integer, parameter :: n_xsec_max = 4 ! Maximum number of cross sections that can plotted
   
   ! Derived limits, based on n_it_max
   integer, parameter :: n_it_max1 = n_it_max + 1
   integer, parameter :: n_it_max2 = n_it_max + 2

   ! Limits updated during runtime
   integer :: n_arrtim_max = 4 ! maximum number of unfiltered arrival times in any event of the cluster (proxy for ntmax0)
   integer :: n_event = 1 ! number of events
   integer :: n_hdf_max = 1 ! Maximum number of events in an HDF file for input as starting locations (nominally, n_event).
   integer :: n_star_max = 1 ! Maximum number of individual stars that can plotted
   
   ! Derived limits, based on n_event
   integer :: n_fp_max = 4 ! max number of free parameters (n_event*4)
   integer :: n_qi_max = 4 ! max number of instances of a particular station-phase (n_event*4)
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Constants

   real, parameter :: pi = 3.1415926536
   real, parameter :: radius = 6371. ! global-average radius of the Earth
   real, parameter :: rpd = 1.7453293e-2 ! radians per degree
   real, parameter :: dpr = 57.29577951  ! degrees per radian
   real, parameter :: dgkmla = dpr/radius ! degrees/km of latitude

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Environment

   integer :: gmt_version
   integer :: ishell ! shell ( 1 for bash, 2 for csh or 3 for zsh) used in GMT scripts
   character(len=32) :: basename ! name for the current run
   character(len=80) :: datadir ! directory containing the data for this run
   character(len=1) :: dirsym ! symbol used to separate directories
   character(len=8) :: mloc_author ! Author
   character(len=40) :: mloc_version
   character(len=12) :: mloc_status
   
   ! File paths
   character(len=16) :: cpt_path
   character(len=16) :: dem_path
   character(len=20) :: ellip_path
   character(len=132) :: gmt_script_dir
   character(len=132) :: mloc_path
   character(len=12) :: star_auto_path
   character(len=12) :: station_path
   character(len=12) :: taup_path

   ! File names
   character(len=40) :: cpt_file
   character(len=80) :: outfile ! Base-name for all output files
   character(len=20) :: star_auto_file
   character(len=16) :: station_master

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Configuration

   ! Logging and verbosity
   logical :: debug
   logical :: verbose_log
   logical :: verbose_screen
   
   ! Longitude basis
   integer :: longitude_range ! Longitudes will be converted to the range -180 to 180°

   ! Windowing residuals
   real :: wind1
   real :: wind2
   real :: windloclim ! Windows are expanded at epicentral distances less than this

   ! Epicentral distance ranges for hypocentroid and cluster vectors
   integer :: nclim ! Number of distance windows in use for cluster vectors
   integer :: nhlim ! Number of distance windows in use for hypocentroid
   real, dimension(3,2) :: clim
   real, dimension(3,2) :: hlim
   
   ! Data weighting and bias correction (in the Jordan & Sverdrup sense)
   logical :: bias_corr ! Bias correction of hypocentroid from events with many more readings
   logical :: data_weight ! Weighting data by inverse of error
   logical :: pttt ! perfect theoretical travel times assumed for hypocentroid

    ! Test mode
    logical :: test_mode ! activation of alternative code segements
   
   ! Selection of a direct calibration subcluster
   integer :: subc_nmin ! minimum number of readings within default distance for selection of subclusters
   integer :: subc_nconnect ! minimum number of connected readings within default distance for selection of subclusters
   real :: subc_delt ! Default distance for selection of subclusters
   logical :: subc_set ! Criteria for selection of subclusters
   
   ! Convergence limits
   integer :: convergence_test_index
   real :: cl_epi_h ! Hypocentroid lat-lon, degrees
   real :: cl_epi_c ! Cluster vector lat-lon, km
   real :: cl_dep_h ! Hypocentroid depth, km
   real :: cl_dep_c ! Cluster vector depth, km
   real :: cl_ot_h ! Hypocentroid OT, sec
   real :: cl_ot_c ! Cluster vector OT, sec
   
   ! Kill events
   integer :: kill_count
   logical :: kill_all
   logical :: kill_one
   character(len=76) :: kill_reason
   
   ! Skip specified stations, phases and authors
   integer :: n_skip ! Number of SKIP commands
   character(len=8), dimension(n_skip_max,3) :: skip_params
   logical :: skip
      
   ! Take starting locations from HDF file of a previous run
   integer, allocatable, dimension(:) :: hdf_day ! read from an hdf file
   integer, allocatable, dimension(:) :: hdf_mon ! read from an hdf file
   integer, allocatable, dimension(:) :: hdf_yr4 ! read from an hdf file
   integer :: nrhdf
   real, allocatable, dimension(:) :: hdf_dep ! read from an hdf file
   real, allocatable, dimension(:) :: hdf_dep_minus ! read from an hdf file
   real, allocatable, dimension(:) :: hdf_dep_plus ! read from an hdf file
   real, allocatable, dimension(:) :: hdf_lat ! read from an hdf file
   real, allocatable, dimension(:) :: hdf_lon ! read from an hdf file
   real, allocatable, dimension(:) :: hdf_time ! read from an hdf file
   character(len=2), allocatable, dimension(:) :: hdf_dep_code ! read from an hdf file
   character(len=10), allocatable, dimension(:) :: hdf_evid ! read from an hdf file
   character(len=60) :: rhdf_filnam
   logical :: read_hdf

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
!  Inversion

   integer :: nstep = n_it_max ! Number of iterations to perform
   real :: tikhonov_factor ! Tikhonov regularization
   double precision, allocatable, dimension(:) :: e_dsvd ! used in dsvd
   logical :: damping = .false.
   logical :: no_xflags = .false. ! Ignore 'x' flags

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Command arguments

   integer :: nvar
   real, dimension(12) :: value
   character(len=4) :: comd
   character(len=76) :: params
   logical :: cmndfil ! reading from a command_file
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Hypocentroidal decomposition

   integer, allocatable, dimension(:) :: jb
   integer, allocatable, dimension(:,:) :: ntq
   integer, allocatable, dimension(:) :: ntqi
   real, allocatable, dimension(:,:,:) :: a
   real, dimension(0:n_it_max1) :: otsh
   real, allocatable, dimension(:,:) :: otsp
   real, allocatable, dimension(:) :: sighatj
   real, allocatable, dimension(:) :: vnhat
   real, allocatable, dimension(:) :: wq
   real, allocatable, dimension(:) :: wq2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Hypocentroid

   integer, dimension(0:n_it_max1) :: hourh
   integer, dimension(0:n_it_max1) :: minh
   integer, dimension(0:n_it_max1) :: nqmth
   real, dimension(4) :: bcorr
   real, dimension(4,0:n_it_max1) :: delx0
   real, dimension(0:n_it_max1) :: depthh
   real, dimension(12) :: dimph
   real, allocatable, dimension(:,:) :: dtmph
   real, dimension(0:n_it_max1) :: ehatsqh
   real, allocatable, dimension(:,:) :: ehiev
   real :: hlatshift
   real :: hlonshift
   real :: hdepthshift
   real :: htimeshift
   real, dimension(0:n_it_max1) :: lath
   real, dimension(0:n_it_max1) :: lathgc ! hypocentroid latitude in geocentric coordinates.
   real, dimension(0:n_it_max1) :: lonh
   real, dimension(4) :: sdxhath
   real, dimension(0:n_it_max1) :: sech
   real, dimension(0:n_it_max1) :: shath
   character(len=5), dimension(4) :: fixprh
   logical :: depthfh ! depth as a free parameter for the hypocentroid
   logical :: latfh ! latitude as a free parameter for the hypocentroid
   logical :: lonfh ! longitude as a free parameter for the hypocentroid
   logical :: ponly ! Use P only or all phase types for hypocentroid
   logical :: timefh ! origin time as a free parameter for the hypocentroid

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Cluster Vectors

   integer, allocatable, dimension(:,:) :: hourp
   integer, allocatable, dimension(:,:) :: minp
   integer :: ntc
   integer :: nqc
   real, allocatable, dimension(:,:) :: depthp
   real, dimension(12) :: dimpc
   real, allocatable, dimension(:,:) :: dimpciev
   real, allocatable, dimension(:,:) :: dtmpc
   real, allocatable, dimension(:,:,:) :: dxp
   real, allocatable, dimension(:,:) :: eci
   real, allocatable, dimension(:,:) :: eciev
   real, dimension(0:n_it_max1) :: ehatsqc
   real, allocatable, dimension(:) :: lat_cf ! latitude for an event declared in the command_file
   real, allocatable, dimension(:) :: lat_hdf ! latitude for an event read from an hdf file
   real, allocatable, dimension(:) :: lon_cf ! longitude for an event declared in the command_file
   real, allocatable, dimension(:) :: lon_hdf ! longitude for an event read from an hdf file
   real :: radius_cvff ! cluster vector fudge factor, km
   real, allocatable, dimension(:,:) :: sdxhatc
   real, dimension(0:n_it_max1) :: shatc
   real, allocatable, dimension(:) :: shatsqci
   real, allocatable, dimension(:) :: time_cf
   real, allocatable, dimension(:) :: time_hdf
   real, allocatable, dimension(:,:) :: latp
   real, allocatable, dimension(:,:) :: latpgc ! event latitude in geocentric coordinates.
   real, allocatable, dimension(:,:) :: lonp
   real, allocatable, dimension(:,:) :: secp
   double precision :: shatsqc
   double precision, allocatable, dimension(:,:) :: vhatc
   character(len=1), allocatable, dimension(:) :: depset_pr ! code for how depth was set for individual events
   character(len=5), allocatable, dimension(:,:) :: fixpr ! 'free' or 'fixed' label for parameters
   logical, allocatable, dimension(:) :: depthf ! depth as a free parameter for individual events
   logical, allocatable, dimension(:) :: latf ! latitude as a free parameter for individual events
   logical, allocatable, dimension(:) :: lonf ! longitude as a free parameter for individual events
   logical, allocatable, dimension(:) :: timef ! origin time as a free parameter for individual events
   
   integer, allocatable, dimension(:) :: mtiev ! number of free parameters for each event
   integer, allocatable, dimension(:,:) :: mindx
   character(len=76), allocatable, dimension(:) :: event_comment
   character(len=10), allocatable, dimension(:) :: evid ! event ID
   character(len=6), allocatable, dimension(:) :: evid_source ! source of an event ID
   character(len=30), allocatable, dimension(:) :: evtnam ! event name, based on date and origin time
   character(len=100), allocatable, dimension(:) :: infile ! event data file pathname
   character(len=20), allocatable, dimension(:) :: infile20 ! event data filename

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Events

   ! Depth-related
   integer :: n_cfd
   real, allocatable, dimension(:,:) :: depth_cf
   real :: depth_default
   real :: depth_default_minus
   real :: depth_default_plus
   real, allocatable, dimension(:,:) :: depth_hdf
   real, allocatable, dimension(:,:) :: depth_inp
   real, allocatable, dimension(:) :: depthp_minus ! negative (shallower) depth uncertainty for each event
   real, allocatable, dimension(:) :: depthp_plus ! positive (deeper) depth uncertainty for each event
   real :: median_constrained_depths
   real :: rescfd
   character(len=1), allocatable, dimension(:) :: depth_inp_c ! constrained depth from input data file
   character(len=1), allocatable, dimension(:) :: depth_hdf_c
   character(len=1), allocatable, dimension(:) :: depth_cf_c
   character(len=1) :: depth_default_c ! Depth code used for default cluster depth
   
   ! Ellipse-related
   real, allocatable, dimension(:) :: alic
   real, allocatable, dimension(:) :: alphac
   real, allocatable, dimension(:) :: alphacg
   real :: alphah
   real, allocatable, dimension(:) :: blic
   real, allocatable, dimension(:,:,:) :: ccv ! cluster vector covariance matrix
   real, allocatable, dimension(:) :: kcritc
   real, allocatable, dimension(:) :: xl1c
   real, allocatable, dimension(:) :: xl1cg
   real :: xl1h
   real, allocatable, dimension(:) :: xl2c
   real, allocatable, dimension(:) :: xl2cg
   real :: xl2h

   ! Hypocenters, magnitudes
   integer, allocatable, dimension(:) :: hour_inp
   integer, allocatable, dimension(:) :: idye
   integer, allocatable, dimension(:) :: iyre
   integer, allocatable, dimension(:) :: min_inp
   integer, allocatable, dimension(:) :: mone
   integer, allocatable, dimension(:) :: nevt
   real, allocatable, dimension(:) :: lat_inp
   real, allocatable, dimension(:) :: lon_inp
   real, allocatable, dimension(:) :: rmag
   real, allocatable, dimension(:) :: sec_inp
   real, allocatable, dimension(:) :: time_inp
   character(len=4), allocatable, dimension(:) :: cal_code_input
   character(len=9), allocatable, dimension(:) :: hypo_author ! author of a hypocenter record in an MNF file
   character(len=9), allocatable, dimension(:) :: magnitude_author ! author of a magnitude record in an MNF file
   character(len=2), allocatable, dimension(:) :: mmag


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Phase data

   integer, allocatable, dimension(:,:) :: ipad
   integer, allocatable, dimension(:,:) :: ipah
   integer, allocatable, dimension(:,:) :: ipam
   integer, allocatable, dimension(:,:) :: ipamo
   integer, allocatable, dimension(:,:) :: ipay
   integer, allocatable, dimension(:,:) :: iptim
   integer, allocatable, dimension(:,:) :: mnf_line
   integer, allocatable, dimension(:) :: nst
   real, allocatable, dimension(:,:) :: pas
   real, allocatable, dimension(:,:) :: resisc
   real, allocatable, dimension(:,:) :: sdread
   real, allocatable, dimension(:,:) :: ttoff
   real, allocatable, dimension(:,:) :: ttsprd
   real, allocatable, dimension(:,:) :: weight
   character(len=52), allocatable, dimension(:,:) :: adslcaa
   character(len=5), allocatable, dimension(:,:) :: agency
   character(len=3), allocatable, dimension(:,:) :: channel
   character(len=8), allocatable, dimension(:,:) :: deployment
   character(len=10), allocatable, dimension(:,:) :: dist_az ! distance-azimuth fields in a phase record
   character(len=1), allocatable, dimension(:,:) :: fcode
   character(len=2), allocatable, dimension(:,:) :: location
   character(len=8), allocatable, dimension(:,:) :: phase
   character(len=8), allocatable, dimension(:,:) :: phase0
   character(len=8), allocatable, dimension(:,:) :: readsrc
   character(len=20), allocatable, dimension(:,:) :: sad
   character(len=5), allocatable, dimension(:,:) :: stname
   logical, allocatable, dimension(:,:) :: connected
   logical :: dflag ! use data flags (x, p, d, k, l)
   logical, allocatable, dimension(:,:) :: fltrcdelt
   logical, allocatable, dimension(:,:) :: fltrcres
   logical, allocatable, dimension(:,:) :: fltrcflag
   logical, allocatable, dimension(:,:) :: fltrc
   logical, allocatable, dimension(:,:) :: fltrhdelt
   logical, allocatable, dimension(:,:) :: fltrhres
   logical, allocatable, dimension(:,:) :: fltrhflag
   logical, allocatable, dimension(:,:) :: fltrh
   logical, allocatable, dimension(:,:) :: rel_depth_phase ! Relative depth phase flag
   logical, allocatable, dimension(:,:) :: rel_phase ! Relative phase flag

   ! Timing error corrections
   integer :: n_terr ! Number of stations with timing corrections
   real, dimension(n_terr_max) :: terr_corr ! Timing corrections
   character(len=5), dimension(n_terr_max) :: terr_stn ! Station codes with timing corrections
   logical :: timing_error_correction ! Timing correction
   
   ! Phase-specific default reading error values
   integer :: n_psdre
   real, dimension(n_psdre_max) :: psdre
   character(len=8), dimension(n_psdre_max) :: psdre_phase
   
   ! Relative depth phases
   integer, allocatable, dimension(:,:) :: h_rdp
   integer :: h_rdp1
   integer :: h_rdp2
   integer :: n_rdpp
   integer, allocatable, dimension(:) :: n_samples 
   real, allocatable, dimension(:,:) :: hrmsi
   real, allocatable, dimension(:,:) :: z_hrdp
   real, allocatable, dimension(:) :: z_test
   character(len=30), allocatable, dimension(:) :: rdpp_evt
   logical :: rdpp ! Relative depth phase plot for a single event
   logical :: rdpp_all ! Relative depth phase plot for all events with data
   
   ! Differential time data
   integer, allocatable, dimension(:,:) :: diff_line ! line number of the data file on which the reading occurs
   integer, allocatable, dimension(:,:) :: idiff ! index of differential time readings
   integer, allocatable, dimension(:) :: idiff0
   integer :: ndiff
   character(len=80) :: diffdatfilnam
   logical :: diffdat
   
   ! Reading errors
   integer, allocatable, dimension(:) :: indexq
   integer, allocatable, dimension(:,:) :: iqiev
   integer :: nqname ! number of station-phase reading errors read from ~.rderr file
   real, allocatable, dimension(:) :: epa
   real, allocatable, dimension(:) :: ere
   real, allocatable, dimension(:) :: qelev
   real, allocatable, dimension(:) :: qlat
   real, allocatable, dimension(:) :: qlon
   real, allocatable, dimension(:,:) :: qres
   real :: rderr_loc_delt ! Distance limit for applying comand RELS
   real :: rderr_loc_p ! Reading error for P at local stations (RELS)
   real :: rderr_loc_s ! Reading error for S at local stations (RELS)
   real :: rderr_min ! Minimum value allowed for reading errors, for phases beyond local distance, except depth phases (MARE)
   real :: rderr_min_depth ! Minimum value allowed for reading errors of teleseismic depth phases (MARE)
   real :: rderr_min_loc ! Minimum value allowed for reading errors at local distance (MARE)
   real, allocatable, dimension(:) :: rderr0
   real, allocatable, dimension(:) :: rdsigma
   character(len=21), allocatable, dimension(:) :: qname
   character(len=21), allocatable, dimension(:) :: qname1
   character(len=100) :: rderrfname
   logical :: read_rderr ! Use a file of empirical reading errors for individual station-phases
   logical :: rels_set ! Status of setting of reading errors for local stations (RELS)

   ! Junk Phases
   integer :: n_junk
   character(len=8), dimension(n_junk_max) :: junk_phase ! phases to ignore
   
   ! Phase re-identification
   integer :: n_no_phreid ! Number of phases for which phase re-identification is prevented (PPRI)
   character(len=1), allocatable, dimension(:,:) :: no_phid ! Flag to prevent phase reidentification for particular readings
   character(len=8), dimension(n_no_phreid_max) :: no_phreid
   logical :: phid ! phase identification
   logical :: plogout ! ~.plog file output, phase identification debugging
   logical, allocatable, dimension(:,:) :: phidird ! Individual reading flag to allow change of phase ID if .true.
   
   ! Spread and baseline offset of specific phases
   integer :: nsprd ! number of phase spread entries read from a ~.ttsprd file
   real, dimension(n_sprd_max) :: sprd
   real, dimension(n_sprd_max) :: ttoffset
   logical :: read_ttsprd ! Use a file of empirical travel-time spreads for different phases
   character(len=8), dimension(n_sprd_max) :: sprdph
   character(len=100) :: ttsprdfname
   
   ! Crustal heterogeneity factor for differential-time data
   real :: cscale ! Spatial scale over which crustal velocities may vary
   real :: vmr ! Percentage variation in crustal velocities

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Stations

   integer, allocatable, dimension(:) :: jdate_off
   integer, allocatable, dimension(:) :: jdate_on
   integer, allocatable, dimension(:,:) :: kcode
   integer, allocatable, dimension(:) :: kstat
   integer :: n_failed_date_range
   integer :: n_radf
   integer :: n_supp_stn ! total number of stations read from supplementary station files
   integer :: n_supp_stn_file ! Number of supplementary station files
   integer :: nkstat
   integer :: nstat1 ! total number of stations read (master + supplementary)
   integer :: nstn_used ! number of stations used in direct calibration
   real, allocatable, dimension(:) :: ahgtr
   real, allocatable, dimension(:,:) :: ahgts
   real, allocatable, dimension(:,:) :: azes
   real, allocatable, dimension(:,:) :: azse
   real, allocatable, dimension(:,:) :: delt
   real, allocatable, dimension(:,:,:) :: dt
   real, allocatable, dimension(:,:) :: elcr
   real, allocatable, dimension(:,:) :: psd
   real, allocatable, dimension(:,:,:) :: s
   real, allocatable, dimension(:) :: sca0
   real, allocatable, dimension(:,:) :: sca0s
   real, allocatable, dimension(:,:) :: sca1s
   real, allocatable, dimension(:,:) :: sca2s
   real, allocatable, dimension(:,:) :: scb1s
   real, allocatable, dimension(:,:) :: scb2s
   real, allocatable, dimension(:) :: sd1
   real, allocatable, dimension(:,:) :: sd1s
   real, allocatable, dimension(:,:) :: sd2s
   real, allocatable, dimension(:,:) :: sdstcs
   real, allocatable, dimension(:,:) :: stacs
   real, allocatable, dimension(:) :: stalat
   real, allocatable, dimension(:) :: stalon
   real, allocatable, dimension(:,:) :: stladg
   real, allocatable, dimension(:,:) :: stlndg
   real, allocatable, dimension(:,:) :: ttcomp
   real, allocatable, dimension(:,:) :: tto
   logical :: read_ad ! Read agency and deployment fields to resolve station code conflicts
   logical :: suppstn ! Use upplementary station code files
   character(len=5), allocatable, dimension(:) :: nstr1
   character(len=5), dimension(n_radf_max) :: radf_stn
   character(len=5), allocatable, dimension(:) :: sta_agency
   character(len=8), allocatable, dimension(:) :: sta_author
   character(len=8), allocatable, dimension(:) :: sta_deployment
   character(len=3), allocatable, dimension(:) :: sta_cha ! station channel code from ADSLC
   character(len=2), allocatable, dimension(:) :: sta_loc ! station location code from ADSLC
   character(len=13), allocatable, dimension(:) :: stn_dcal_used 
   character(len=100), dimension(n_sstn_max) :: suppfilnam

   ! List of stations suspected of reporting bogus depth phases
   integer :: n_bdp ! Number of entries in the 'bad depth phases' list
   character(len=100) :: bdp_filnam
   character(len=108), dimension(n_bdp_max) :: bdp_station
   logical :: bdp_list ! Use list of stations known to report bad depth phases
   
   ! Duplicate entries in station lists
   integer :: n_dupe
   integer :: n_dupe_conflict
   integer :: n_dupe_minor
   integer :: n_dupe_pure
   integer :: n_dupe_significant
   character(len=40), dimension(n_sta_max) :: duplication

   ! Missing station codes and other station statistics
   integer, dimension(n_miss_sta_max) :: n_miss_sta
   integer :: n_miss_sta_total ! number of stations missing from the station files
   integer, allocatable, dimension(:,:) :: ndat
   integer, allocatable, dimension(:,:) :: ndatc
   integer, allocatable, dimension(:,:) :: ndatdl
   integer, allocatable, dimension(:) :: ndatfl
   integer, allocatable, dimension(:,:) :: ndatpr
   integer, allocatable, dimension(:) :: nmiss ! # missing station codes for each event
   character(len=20), allocatable, dimension(:,:) :: missta ! list of missing station codes for each event
   character(len=20), dimension(n_miss_sta_max) :: n_miss_sta_list ! List of unique missing station codes over all events

   ! NEIC station metadata
    logical :: nsmd ! Search of NEIC station metadata for missing station codes

   ! ISF-FDSN station metadata
   logical :: ifsm ! Search of ISC-FDSN station metadata for missing station codes

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Travel time

   ! Local velocity model
   integer, parameter :: indph = 30000 ! Index controlling what phases are returned from hyposat
   real :: dlimlocmod ! Epicentral distance limit for local model (LMOD)
   real :: zlimlocmod ! Focal depth limit for local model (LMOD)
   character(len=80) :: locmodfname ! Pathname to local crustal model file (LMOD)
   logical :: locmod ! Use a local velocity model (LMOD)
      
   ! Phases not in ak135
   real :: lg_a ! Intercept, in sec (LGTT)
   real :: lg_b ! Slope, in sec/degree (LGTT)
   real :: lg_min ! Minimum epicentral distance (degrees) at which Lg travel times will be calculated (LGTT)
   real :: p__a ! Intercept, in sec (P_TT)
   real :: p__b ! Slope, in sec/degree (P_TT)
   real :: p__min ! Minimum epicentral distance (degrees) at which P_ travel times will be calculated (P_TT)
   real :: rg_a ! Intercept, in sec (RGTT)
   real :: rg_b ! Slope, in sec/degree (RGTT)
   real :: rg_min ! Minimum epicentral distance (degrees) at which Rg travel times will be calculated (RGTT)
   real :: tphase_a ! Intercept, in sec (TPTT)
   real :: tphase_b ! Slope, in sec/degree (TPTT)
   
   ! Station elevation corrections
   integer :: tt_corr ! Station elevation correction (CORR)
   real :: pcrvel ! P velocity for station elevation correction
   real :: scrvel ! S velocity for station elevation correction
   
   ! tau-P
   character(len=5) :: taup_model ! (TAUP)
   integer, parameter :: jsrc = 150 ! Maximum number of discrete model slowness samples above the maximum source depth of interest
   integer, parameter :: jseg = 30 ! Maximum number of different types of travel-times considered
   integer, parameter :: jbrn = 100 ! Maximum number of different travel-time branches to be searched
   integer, parameter :: jout = 2500 ! Maximum length of all travel-time branches strung together
   integer, parameter :: jtsm = 350 ! Maximum length of the tau depth increments
   ! Derived parameters:
   integer, parameter :: jxsm = jbrn ! Maximum number of x-values needed for the depth increments
   integer, parameter :: jbrnu = jbrn ! Maximum length of the up-going branches
   integer, parameter :: jbrna = jbrn ! Maximum length of branches which may need re-interpolation
   integer, parameter :: jrec = jtsm + jxsm
   integer, parameter :: jtsm0 = jtsm + 1
   logical :: log_tt ! Log file for intermediate tau-p results
   logical, dimension(2) :: prnt ! common /prtflc/
   logical, dimension(jseg) :: segmsk ! (common /prtflc/)
   
   ! Main phase list (mpl), travel times and derivatives
   integer, parameter :: max_mpl = 60
!    integer :: nphase ! Number of entries in the main phase list
!    real, dimension(max_mpl) :: dddp
!    real, dimension(max_mpl) :: dtdd
!    real, dimension(max_mpl) :: dtdh
!    real, dimension(max_mpl) :: tt ! travel times
!    character(len=8), dimension(max_mpl) :: phcd ! phase codes
   
   ! Bounce point topography
   integer, dimension(n_bp_lon,n_bp_lat) :: bp_topo
   logical :: bptc ! Bounce point topography correction, for pP, sP and pwP
   
   ! Empirical TT output files
   integer :: n_ttou ! Number of phase-specific empirical TT output files
   character(len=8), dimension(n_ttou_max) :: ttou_phase ! List of phases for empirical TT output files
   logical :: ttou ! Output of empirical TTs for specific phases
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Calibration

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
   real, allocatable, dimension(:,:,:,:) :: rcv
   real :: rdbt ! radius of doubt
   real, allocatable, dimension(:,:) :: s_cal
   real :: s1gt
   real :: s2gt
   real, allocatable, dimension(:,:) :: sdstcs_cal
   real, allocatable, dimension(:,:) :: stacs_cal
   real, allocatable, dimension(:) :: w12
   real, allocatable, dimension(:) :: w3
   real, allocatable, dimension(:) :: w4
   logical, allocatable, dimension(:,:) :: cal_event
   logical :: ot_cal = .false.
   character(len=1), allocatable, dimension(:) :: cal_par
   character(len=12) :: icaltype_pr
   
   ! Direct and indirect calibration
   real, allocatable, dimension(:) :: eqr


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Output files

   ! Tomography output files
   integer, dimension(n_tomo_max) :: itomo ! index of type of output
   integer, dimension(n_tomo_max) :: n_cell_lat ! # of cells in latitude
   integer, dimension(n_tomo_max) :: n_cell_lon ! # of cells in longitude
   integer :: nitomo ! Number of TOMO commands given
   real, dimension(n_tomo_max) :: cell_lat_inc ! cell dimension for latitude in degrees
   real, dimension(n_tomo_max) :: cell_lon_inc ! cell dimension for longitude in degrees
   real, dimension(n_tomo_max) :: grid_origin_lat ! latitude of grid origin
   real, dimension(n_tomo_max) :: grid_origin_lon ! longitude of grid_origin
   character(len=8), dimension(n_tomo_max) :: tomo_phase ! seismic phase name
   logical, dimension(n_tomo_max) :: tomo_sub ! if subdividing is requested
   
   character(len=20), allocatable, dimension(:) :: annotation
   character(len=64), allocatable, dimension(:) :: hypocenter_list
   character(len=180), allocatable, dimension(:) :: hdfline

   ! Print variables
   integer :: hourpr
   integer :: minpr
   real :: depout
   real :: latout
   real :: lonout
   real :: secpr
   character(len=3), dimension(0:1) :: biaspr
   character(len=3) :: data_weight_pr
   character(len=5), allocatable, dimension(:) :: depthpr ! labels for free parameters for each event
   character(len=3) :: dflag_pr
   character(len=3), allocatable, dimension(:) :: latpr
   character(len=3), allocatable, dimension(:) :: lonpr
   character(len=3) :: pttt_pr
   character(len=80), dimension(2) :: stacpr
   character(len=20), dimension(2) :: tablpr
   character(len=11), allocatable, dimension(:) :: timepr
   character(len=3), dimension(0:1) :: weigpr

   ! Output file for large residuals
   real :: lres ! Threshold for large cluster residuals
   logical :: lresout ! .lres file for large cluster residuals

   ! Output file to pass for input to BAYESLOC
   logical :: blocout ! data file for import into BAYESLOC
   
   ! .datf file
   logical :: datfout ! ~.datf file of flagged, re-identified arrival time data.
   
   ! GCCEL files
   character(len=80), dimension(9) :: commentary_buffer
   character(len=80) :: commentary_fname
   character(len=132) :: gcat_folder
   logical :: gccelout ! data file for GCCEL
   
   ! Output file for limited distance range (usually used for Pg/Pn cross-over range)
   real :: dist1 ! Minimum distance for oldr
   real :: dist2 ! Max distance for oldr
   logical :: oldr_out ! Output for limited distance range (oldr)
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Plotting

   ! Plot stars
   integer :: n_star
   integer, allocatable, dimension(:) :: iev_star
   real, allocatable, dimension(:) :: star_size
   real, dimension(n_star_table,3) :: star_table
   logical :: star_auto
   logical :: star_plot
   
   ! Plot stations
   integer :: n_stat
   real, dimension(n_stat_max,3) :: stat_data
   logical :: stat_plot
   
   ! Plot ellipses
   integer :: n_ellipse
   real, dimension(n_ellipse_nmax,5) :: ellipse_data
   logical :: ellipse_plot
   
   ! Digital fault map
   integer :: n_fault_map
   character(len=100), dimension(n_fault_map_max) :: fault_map_filename
   character(len=30) :: fault_plot_style
   logical :: fault_map
   
   ! Plot topography in GMT
   character(len=132) :: dem2_filename = ' '
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
   logical :: phrp = .false.

   ! Empirical path anomaly plots
   integer :: n_epa_plot
   real, dimension(n_epa_max) :: epa_plot_distance
   character(len=8), dimension(n_epa_max) :: epa_plot_phase
   logical :: epa_plot
   logical, dimension(n_epa_max) :: epa_plot_raypath

   ! Event residual plots
   integer :: n_evrp
   real, allocatable, dimension(:) :: evrp_delta_min
   real, allocatable, dimension(:) :: evrp_delta_max
   real :: evrp_delta_min_all
   real :: evrp_delta_max_all
   character(len=30), allocatable, dimension(:) :: evrp_event
   logical :: evrp ! Residual plots for individual events
   logical :: evrp_all ! Residual plots for all events
    
   ! Plot types and options
   logical :: eplt ! confidence ellipse plot
   logical :: fdhp ! Focal depth histogram plot
   logical :: gmtp ! Enable GMT plotting
   logical, allocatable, dimension(:,:) :: plot
   logical :: reduced ! Use reduced velocities in type 6 and 7 TT plots
   logical :: rpth ! Plot direct calibration raypaths
   logical :: splt ! seismicity plot
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

   ! Plot slip data (PSLP)
   integer :: n_pslp
   character(len=8), allocatable, dimension(:) :: pslp_color
   character(len=80), allocatable, dimension(:) :: pslp_filnam
   logical :: pslp
   
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  Logical unit numbers

   integer, parameter :: io_alloc_log = 43 ! array allocation log
   integer, parameter :: io_bayes = 11 ! BAYESLOC output
   integer, parameter :: io_bdp = 12 ! Bad depth phase station list
   integer, parameter :: io_bptc_log = 13 ! Bounce point correction log file
   integer, parameter :: io_cal = 14 ! .cal file
   integer, parameter :: io_cfil = 15 ! command file
   integer, parameter :: io_dat0 = 16 ! .dat0 file
   integer, parameter :: io_datf = 17 ! .datf file
   integer, parameter :: io_depth_phase = 18 ! depth phases output
   integer, parameter :: io_diffdat = 19 ! Differential time data
   integer, parameter :: io_dpnc = 20 ! Output file for depth phase name changes
   integer, parameter :: io_gccel = 21 ! GCCEL output
   integer, parameter :: io_gmt = 22 ! GMT scripts
   integer, parameter :: io_ifsm = 23 ! ISC-FDSN station metadata file
   integer, parameter :: io_in = 24 ! Event data files
   integer, parameter :: io_junk = 25 ! Input file of junk phase names that will not be read
   integer, parameter :: io_locmod = 26 ! Custom crustal velocity model
   integer, parameter :: io_log = 27 ! log file
   integer, parameter :: io_lres = 28 ! .lres file (large cluster residuals)
   integer, parameter :: io_nsmd = 29 ! NEIC station metadata file
   integer, parameter :: io_oldr = 30 ! Output for limited distance range
   integer, parameter :: io_out = 31 ! Used for several output files
   integer, parameter :: io_pdf = 32 ! Probability density function file
   integer, parameter :: io_plog = 33 ! phase re-identification log
   integer, parameter :: io_rderr = 34 ! .rderr (reading errors file)
   integer, parameter :: io_rhdf = 35 ! HDF file to set starting locations
   integer, parameter :: io_shifted = 36 ! Shifted phase_data after indirect calibration
   integer, parameter :: io_stn_log = 37 ! station data log
   integer, parameter :: io_taup = 38 ! .hed and .tbl files (ak135)
   integer, parameter :: io_tt_log = 39 ! Log file for debugging travel-time calculations
   integer, parameter :: io_ttou = 40 ! Empirical TT data for specific phases
   integer, parameter :: io_ttsprd = 41 ! .ttsprd file
   integer, parameter :: io_xdat = 42 ! .xdat file

contains

!*****************************************************************************************
   subroutine oops (msg)
   
   ! Error report and stop
   
      character(len=*) :: msg
   
      write (*,'(/a)') 'Oops! Fatal error in '//trim(msg)
   
      stop
      
   end subroutine oops

      
!*****************************************************************************************
   subroutine warnings (msg)
   
   ! Warnings
   
      character(len=*) :: msg
         
      write (*,'(t3,a)') 'Warning from '//trim(msg)
         
      return
      
   end subroutine warnings

      
!*****************************************************************************************
   subroutine fyi (msg)
   
   ! Informational messages
   
      character(len=*) :: msg
         
      write (*,'(t3,a)') 'FYI from '//trim(msg)
         
      return
      
   end subroutine fyi
      
!*****************************************************************************************
   subroutine allocation_error (procedure_name, array_name, error)
   
      integer, intent(in) :: error
      character(len=*), intent(in) :: array_name
      character(len=*), intent(in) :: procedure_name
      character(len=132) :: msg
      
      write (msg,'(a,i8,3a)') trim(procedure_name)//': allocation error (', error,&
       ') for array "', trim(array_name), '"'
      call oops (trim(msg))
      
      return
   
   end subroutine allocation_error

   
!*****************************************************************************************
   subroutine deallocation_error (procedure_name, array_name, error)
   
      integer, intent(in) :: error
      character(len=*), intent(in) :: array_name
      character(len=*), intent(in) :: procedure_name
      character(len=132) :: msg
      
      write (msg,'(a,i8,3a)') trim(procedure_name)//': deallocation error (', error,&
       ') for array "', trim(array_name), '"'
      call oops (trim(msg))
      
      return
   
   end subroutine deallocation_error

   
!*****************************************************************************************
   integer function lunit ()
   
   ! find an unopened fortran unit number above the ones defined in this module.
         
      logical :: lopen
      
      do lunit = 43,99
         inquire (unit=lunit,opened=lopen)
         if (.not.lopen) return
      end do
      lunit = 0
      write (*,'(t3,a)') 'Oops! Fatal error in lunit: failed to find an unopened unit between 43 and 99'
      stop
      
   end function lunit

      
!*****************************************************************************************
end module mloc_declare

