!This suborutine calulates the stellar flux based on Maxted & Gill (2018)
!qpower2 routine
!This works for cases where p (Rp/R*) < 0.2 R*
subroutine qpower2(z,p,c,alpha,flux,nz)
use constants
implicit none

!In/Out variables
integer, intent(in) :: nz
real(kind=mireal), intent(in) :: z(nz),p,c,alpha
real(kind=mireal),intent(out) :: flux(nz)
!
integer :: i

  flux = 1.d0
  if( all(z > 1.d0 - p) .and. all (z < 1.d0 + p) ) then

     call q2(z,p,c,alpha,flux,nz)

  else if( all ( z <= 1.d0 - p ) ) then

     call q1(z,p,c,alpha,flux,nz)

  else

    do i = 1, nz
      if( z(i) > 1.d0 - p .and. z(i) < 1.d0 + p ) then
        call q2(z(i),p,c,alpha,flux(i),1)
      else if( z(i) <= 1.d0 - p ) then
        call q1(z(i),p,c,alpha,flux(i),1)
      end if
    end do

  end if

end subroutine qpower2


subroutine q1(z,p,c,alpha,flux,nz)
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: nz
  real(kind=mireal), intent(in) :: z(nz),p,c,alpha
  real(kind=mireal), intent(out) :: flux(nz)
  !Local variables
  real(kind=mireal) :: I0, g, p2
  real(kind=mireal), dimension(nz) :: zt, s, c0, c2

  I0 = alpha + 2. / (pi * ( alpha - c * alpha + 2.) )
  g  = 0.5 * alpha

  p2 = p*p

  call clip(abs(z), 0.d0, 1.d0 - p, zt, nz)

  s  = 1. - zt * zt
  c0 = 1. - c + c*s**g
  c2 = 0.5*alpha*c*s**(g-2.)*( (alpha-1.)*zt*zt - 1. )

  flux = ( c0 + 0.25*p2*c2 - 0.125*alpha*c*p2*s**(g-1.) )
  flux = 1. - I0*pi*p2*flux

end subroutine


subroutine q2(z,p,c,alpha,flux,nz)
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: nz
  real(kind=mireal), intent(in) :: z(nz),p,c,alpha
  real(kind=mireal), intent(out) :: flux(nz)
  !Local variables
  real(kind=mireal) :: I0, g, p2, finfo
  real(kind=mireal), dimension(nz) :: zt, d, ra, rb, sa, sb, q
  real(kind=mireal), dimension(nz) :: w2, w, b0, b1, b2, a0, a1, aq
  real(kind=mireal), dimension(nz) :: J1, J2, K1, K2, d0, d1
  real(kind=mireal), dimension(nz) :: ztra, dzt

  finfo = tiny(0.d0)

  I0 = alpha + 2. / (pi * ( alpha - c * alpha + 2.) )
  g  = 0.5 * alpha

  p2 = p*p

  call clip(abs(z), 1.d0 - p , 1.d0 + p, zt, nz)
  call clip( (zt*zt - p2 + 1.d0)/(2*zt),0.d0,1.d0,d,nz)

  ra = 0.5 * ( zt - p + d )
  rb = 0.5 * ( 1. + d )

  call clip(1.d0-ra*ra,finfo,1.d0,sa,nz)
  call clip(1.d0-rb*rb,finfo,1.d0,sb,nz)
  call clip((zt - d)/p,-1.d0,1.d0,q,nz)

  dzt = d-zt
  w2 = p2 - dzt*dzt

  call clip(w2,finfo,1.d0,w,nz)
  w = sqrt(w)

  b0 = 1. - c + c*sa**g
  b1 = - alpha * c * ra * sa**(g-1.)
  b2 = 0.5*alpha*c*sa**(g-2.)*( (alpha-1.d0)*ra*ra - 1.d0 )

  ztra = zt-ra
  a0 = b0 + b1*ztra + b2*ztra*ztra
  a1 = b1+2*b2*ztra
  aq = acos(q)

  J1 = ( a0*dzt - (2.d0/3.)*a1*w2 + 0.25*b2*dzt*( 2*dzt*dzt-p2 ) ) * w &
       + ( a0 *p2 + 0.25*b2*p2*p2)*aq

  call clip(1.d0-q*q,0.d0,1.d0,J2,nz)
  J2 = alpha*c*sa**(g-1.)*p2*p2*(0.125*aq + (1.d0/12.)*q*(q*q-2.5d0)*sqrt(J2))

  d0 = 1.d0 - c + c * sb**g
  d1 = -alpha*c*rb*sb**(g-1)

  call clip(1.d0-d*d,0.d0,1.d0,K1,nz)
  K1 = (d0-rb*d1)*acos(d) + ( (rb*d + (2.d0/3.)*(1.-d*d))*d1 - d*d0 )*sqrt(K1)

  K2 = (1./3.)*c*alpha*sb**(g + 0.5)*(1.d0 - d )

  flux = 1.d0 - I0 * (J1 - J2 + K1 - K2)

end subroutine

subroutine clip(a,amin,amax,aout,na)
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: na
  real(kind=mireal), intent(in) :: a(na), amin, amax
  real(kind=mireal), intent(out) :: aout(na)

  aout = a

  where (a < amin)
    aout = amin
  else where (a > amax)
    aout = amax
  end where

 end subroutine