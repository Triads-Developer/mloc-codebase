   ! nevmax and nmax are defined in mloc.inc
   integer, parameter :: nmax = nevmax
   integer, parameter :: nmax2 = nmax*2

   real :: rsvx, rsvy
   common /rsv/ rsvx(nevmax), rsvy(nevmax)

   real :: covs, cvs
   common /tst/ covs(nmax2,nmax2), cvs(nmax,2)
