module declare_constants

   implicit none
   save
   
   real, parameter :: pi = 3.1415926536
   real, parameter :: radius = 6371. ! global-average radius of the Earth
   real, parameter :: rpd = 1.7453293e-2 ! radians per degree
   real, parameter :: dpr = 57.29577951  ! degrees per radian
   real, parameter :: dgkmla = dpr/radius ! degrees/km of latitude

end module declare_constants

