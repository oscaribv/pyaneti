module constants

implicit none

  integer, parameter :: mireal = selected_real_kind(8)
  integer, parameter :: mkernel = 4
  real(kind=mireal), parameter :: pi = 3.14159265358979311599796346854418516 ! acos(-1.d0)
  real(kind=mireal), parameter :: two_pi = 2.0*pi
  real(kind=mireal), parameter :: log_two_pi = log(2.0*pi)
  real(kind=mireal), parameter :: uno = 1.d0
  real(kind=mireal), parameter :: fmin=1.d-8
  real(kind=mireal), parameter :: small = 1.d-5
  real(kind=mireal), parameter :: a_factor = 2.d0
  integer, parameter :: imax = int(1e8)
  !Physical parameters
  real(kind=mireal), parameter :: S_radius_SI = 6.957d8 !R_sun
  real(kind=mireal), parameter :: S_GM_SI = 1.3271244d20 ! G M_sun
  real(kind=mireal), parameter :: G_cgs = 6.67508d-11*1.d3 !Gravitational constant cgs
  real(kind=mireal), parameter :: sind =  86400. !Secs in a day
  real(kind=mireal), parameter :: sind2 =  7464960000.d0 !(Secs in a day)**2

contains

end module constants