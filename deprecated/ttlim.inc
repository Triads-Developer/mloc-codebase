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

   integer :: io_taup_log ! Log file unit number for tau-p debugging
   logical :: verbose_taup, log_taup
   common /taup_print/ verbose_taup, log_taup, io_taup_log
