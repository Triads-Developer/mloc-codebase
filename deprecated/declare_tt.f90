!> Declaration of variables related to travel time calculations

module declare_tt

   use declare_limits

   implicit none
   save
   
   ! Local Model
   integer :: indph
   real :: dlimlocmod = 0. ! (LMOD)
   real :: zlimlocmod = 0. ! (LMOD)
   character(len=80) :: locmodfname = ' ' ! (LMOD)
   logical :: locmod ! (LMOD)
      
   ! Phases not in ak135
   real :: lg_a ! (LGTT)
   real :: lg_b ! (LGTT)
   real :: lg_min ! (LGTT)
   real :: p__a ! (P_TT)
   real :: p__b ! (P_TT)
   real :: p__min ! (P_TT)
   real :: rg_a ! (RGTT)
   real :: rg_b  ! (RGTT)
   real :: rg_min ! (RGTT)
   real :: tphase_a ! (TPTT)
   real :: tphase_b ! (TPTT)
   
   ! Station elevation corrections
   integer :: tt_corr ! (CORR)
   real :: secv_p ! (SECV)
   real :: secv_s ! (SECV)
   
   ! Tau-P
   character(len=5) :: taup_model ! (TAUP)
   ! User-settable parameters:
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
   logical :: log_taup
   logical :: verbose_taup
   
   ! Main phase list (mpl), travel times and derivatives
   integer, parameter :: max_mpl = 60
   integer :: nphase ! Number of entries in the main phase list
   real, dimension(max_mpl) :: dddp
   real, dimension(max_mpl) :: dtdd
   real, dimension(max_mpl) :: dtdh
   real, dimension(max_mpl) :: tt ! travel times
   character(len=8), dimension(max_mpl) :: phcd ! phase codes
   
   ! Bounce point topography
   integer, dimension(bp_nlon,bp_nlat) :: bp_topo
   logical :: bptc
   
   ! Empirical TT output files
   integer :: n_ttou
   character(len=8), dimension(n_ttou_max) :: ttou_phase
   logical :: ttou
   
   ! Station elevation correction
   real :: pcrvel
   real :: scrvel
   
end module declare_tt
