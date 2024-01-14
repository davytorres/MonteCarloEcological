subroutine mct(nnn,rho,sdnor,sdnor2,expectation,var)
 implicit none
! Input
  integer nnn
  double precision rho
  double precision sdnor,sdnor2
! Output
  double precision expectation,var
! Local
  integer ii,nnnm1
  double precision hyp,dhyp,expr2
  double precision fac,s,ssz,sstop,ssbot
      
  ii = 2
  fac = 1.d0
  s = dble(nnn)
  ssz = s/2.d0

  do while (ssz .gt. 3)
     fac = fac*(s-dble(ii))/(s-dble(ii+1))
     ssz = (s-dble(ii))/2.d0
     ii = ii + 2
  end do

  call GAMMA(ssz,sstop)
  call GAMMA(ssz-.5d0,ssbot)      
  fac = fac*sstop/ssbot
      
  expectation = 2.d0*(fac**2)/dble(nnn-1)      

  fac = fac*((1.d0-rho**2)**((s-1.d0)/2.d0))
  fac = fac/sqrt(acos(-1.d0))
  fac = fac*sdnor/sdnor2

  expectation = expectation*rho

  call hypser(0.5d0,0.5d0,dble(nnn+1)/2.d0,rho**2,hyp,dhyp)
  expectation = expectation*hyp

  nnnm1 = nnn - 1
  call hypser(1.d0,1.d0,dble(nnnm1)/2.d0+1.d0,rho**2,hyp,dhyp)
  expr2 = 1.d0 - hyp*(1.d0-rho**2)*dble(nnnm1-1)/dble(nnnm1)

  var = expr2 - expectation**2

  return
end
