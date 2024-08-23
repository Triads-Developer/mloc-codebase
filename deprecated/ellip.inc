   integer, parameter :: numph = 57
   integer, parameter :: mdel = 50

   character(len=8) :: phcod
   common /tau1/ phcod(numph)

   integer :: phnch, np
   common /tau2/ phnch(numph), np(numph)

   real :: dpth, delta, di1, di2, t0, t1, t2
   common /tau3/ dpth(6), delta(numph,mdel), di1(numph), di2(numph), t0(numph,mdel,6),&
    t1(numph,mdel,6), t2(numph,mdel,6)

