      double precision function gcoeff7()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff7 =  + UBAR**(-2) * (  - 12.D0*b*t1*vlwm*xnsq**2 + 12.D0*b*
     &    t1*vlwm - 12.D0*b*t1**2*vlwm*xnsq + 12.D0*b*t1**2*vlwm*
     &    xnsq**2 + 12.D0*b*vlwm*xnsq - 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + UBAR**(-1) * (  - 36.D0*b*t1*vlwm*xnsq + 84.D0
     &    *b*t1*vlwm*xnsq**2 - 12.D0*b*t1*xnsq + 12.D0*b*t1*xnsq**2 - 
     &    84.D0*b*vlwm*xnsq + 36.D0*b*vlwm - 12.D0*b*xnsq + 12.D0*b )
      gcoeff7 = gcoeff7 + INVG * ( 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2*xnsq
     &     + 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2 + 3.D0*b*t1**(-1)*ro*vlsm*
     &    xnsq**2 + 6.D0*b*t1**(-1)*ro*vltm**2*xnsq + 6.D0*b*t1**(-1)*
     &    ro*vltm**2 - 18.D0*b*t1**(-1)*ro*vlwm*xnsq + 12.D0*b*t1**(-1)
     &    *ro*vlwm*xnsq**2 + 6.D0*b*t1**(-1)*ro*vlwm + 6.D0*b*t1**(-1)*
     &    ro*vdt*xnsq + 6.D0*b*t1**(-1)*ro*vdt - 1.D0/2.D0*b*t1**(-1)*
     &    ro*pi**2*xnsq - 1.D0/2.D0*b*t1**(-1)*ro*pi**2 - 15.D0/8.D0*b*
     &    t1**(-1)*ro**2*vlpm**2*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**2*
     &    vlpm**2 - 9.D0/2.D0*b*t1**(-1)*ro**2*vlsm*vltm*xnsq**2 - 9.D0/
     &    2.D0*b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 3.D0*b*t1**(-1)*
     &    ro**2*vlsm*xnsq**2 + 9.D0/4.D0*b*t1**(-1)*ro**2*vlsm**2*
     &    xnsq**2 + 9.D0/2.D0*b*t1**(-1)*ro**2*vltm*vlwm*xnsq + 3.D0/2.D
     &    0*b*t1**(-1)*ro**2*vltm*xnsq - 3.D0*b*t1**(-1)*ro**2*vltm*
     &    xnsq**2 - 21.D0/4.D0*b*t1**(-1)*ro**2*vltm**2*xnsq - 3.D0*b*
     &    t1**(-1)*ro**2*vltm**2 + 3.D0/2.D0*b*t1**(-1)*ro**2*vlwm*xnsq
     &     - 3.D0*b*t1**(-1)*ro**2*vlwm*xnsq**2 )
      gcoeff7 = gcoeff7 + INVG * ( 9.D0/4.D0*b*t1**(-1)*ro**2*vlwm**2*
     &    xnsq + 9.D0/2.D0*b*t1**(-1)*ro**2*vdw*xnsq - 9.D0/2.D0*b*
     &    t1**(-1)*ro**2*vdw*xnsq**2 - 3.D0*b*t1**(-1)*ro**2*vdt*xnsq
     &     - 9.D0/2.D0*b*t1**(-1)*ro**2*vdt*xnsq**2 - 3.D0*b*t1**(-1)*
     &    ro**2*vdt - 1.D0/8.D0*b*t1**(-1)*ro**2*pi**2*xnsq + 3.D0/4.D0
     &    *b*t1**(-1)*ro**2*pi**2*xnsq**2 + 1.D0/4.D0*b*t1**(-1)*ro**2*
     &    pi**2 - 3.D0*b*t1**(-1)*vlsm*omro**(-1)*xnsq**2 + 3.D0*b*
     &    t1**(-1)*vlsm*xnsq**2 + 3.D0*b*t2**(-1)*ro*vlsm*xnsq**2 + 3.D0
     &    /4.D0*b*t2**(-1)*ro**2*vlpm**2 + 3.D0*b*t2**(-1)*ro**2*
     &    vlwm**2 + 3.D0*b*t2**(-1)*ro**2*vdw - 1.D0/4.D0*b*t2**(-1)*
     &    ro**2*pi**2 - 3.D0*b*t2**(-1)*vlsm*omro**(-1)*xnsq**2 + 3.D0*
     &    b*t2**(-1)*vlsm*xnsq**2 + 9.D0/2.D0*b*ro*vlpm**2*xnsq + 9.D0/
     &    2.D0*b*ro*vlpm**2 + 12.D0*b*ro*vlsm*vltm*xnsq**2 + 6.D0*b*ro*
     &    vlsm*vlwm*xnsq**2 - 12.D0*b*ro*vlsm*xnsq**2 - 9.D0/2.D0*b*ro*
     &    vlsm**2*xnsq**2 - 12.D0*b*ro*vltm*vlwm*xnsq + 12.D0*b*ro*vltm
     &    *xnsq )
      gcoeff7 = gcoeff7 + INVG * (  - 6.D0*b*ro*vltm + 12.D0*b*ro*
     &    vltm**2*xnsq + 12.D0*b*ro*vltm**2 - 24.D0*b*ro*vlwm*xnsq + 24.
     &    D0*b*ro*vlwm*xnsq**2 + 6.D0*b*ro*vlwm - 6.D0*b*ro*vlwm**2*
     &    xnsq + 6.D0*b*ro*vlwm**2 - 12.D0*b*ro*vdw*xnsq + 6.D0*b*ro*
     &    vdw*xnsq**2 + 6.D0*b*ro*vdw + 6.D0*b*ro*vdt*xnsq + 12.D0*b*ro
     &    *vdt*xnsq**2 + 12.D0*b*ro*vdt + 1.D0/2.D0*b*ro*pi**2*xnsq - 3.
     &    D0/2.D0*b*ro*pi**2*xnsq**2 - 3.D0/2.D0*b*ro*pi**2 - 3.D0*b*
     &    ro**2*vltm**2*xnsq + 3.D0*b*ro**2*vlwm**2*xnsq + 3.D0*b*ro**2
     &    *vdw*xnsq - 3.D0*b*ro**2*vdt*xnsq + 12.D0*b*vlsm*omro**(-1)*
     &    xnsq**2 - 12.D0*b*vlsm*xnsq**2 + 72.D0*b*vlwm*xnsq - 48.D0*b*
     &    vlwm*xnsq**2 - 24.D0*b*vlwm + 3.D0*t1**(-1)*ro*vlpm*vlsm*xnsq
     &     + 3.D0*t1**(-1)*ro*vlpm*vlsm - 6.D0*t1**(-1)*ro*vlpm*vltm*
     &    xnsq - 6.D0*t1**(-1)*ro*vlpm*vltm - 3.D0/2.D0*t1**(-1)*ro*
     &    vlpm**2*xnsq - 3.D0/4.D0*t1**(-1)*ro*vlpm**2*xnsq**2 - 3.D0/2.
     &    D0*t1**(-1)*ro*vlpm**2 - 6.D0*t1**(-1)*ro*vdmp*xnsq - 3.D0*
     &    t1**(-1)*ro*vdmp*xnsq**2 )
      gcoeff7 = gcoeff7 + INVG * (  - 6.D0*t1**(-1)*ro*vdmp - 1.D0/2.D0
     &    *t1**(-1)*ro*pi**2*xnsq - 1.D0/4.D0*t1**(-1)*ro*pi**2*xnsq**2
     &     - 1.D0/2.D0*t1**(-1)*ro*pi**2 - 21.D0/4.D0*t1**(-1)*ro**2*
     &    vlpm*vlsm*xnsq - 9.D0/2.D0*t1**(-1)*ro**2*vlpm*vlsm + 21.D0/2.
     &    D0*t1**(-1)*ro**2*vlpm*vltm*xnsq + 9.D0*t1**(-1)*ro**2*vlpm*
     &    vltm - 3.D0/2.D0*t1**(-1)*ro**2*vlpm*xnsq - 3.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm + 21.D0/8.D0*t1**(-1)*ro**2*vlpm**2*xnsq
     &     + 9.D0/4.D0*t1**(-1)*ro**2*vlpm**2*xnsq**2 + 9.D0/4.D0*
     &    t1**(-1)*ro**2*vlpm**2 + 21.D0/2.D0*t1**(-1)*ro**2*vdmp*xnsq
     &     + 9.D0*t1**(-1)*ro**2*vdmp*xnsq**2 + 9.D0*t1**(-1)*ro**2*
     &    vdmp + 7.D0/8.D0*t1**(-1)*ro**2*pi**2*xnsq + 3.D0/4.D0*
     &    t1**(-1)*ro**2*pi**2*xnsq**2 + 3.D0/4.D0*t1**(-1)*ro**2*pi**2
     &     + 9.D0/4.D0*t1**(-1)*ro**3*vlpm*vlsm*xnsq - 9.D0/4.D0*
     &    t1**(-1)*ro**3*vlpm*vltm*xnsq - 9.D0/4.D0*t1**(-1)*ro**3*vlpm
     &    *vlwm*xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm*xnsq - 9.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm**2*xnsq )
      gcoeff7 = gcoeff7 + INVG * (  - 9.D0/2.D0*t1**(-1)*ro**3*vdmp*
     &    xnsq - 3.D0/8.D0*t1**(-1)*ro**3*pi**2*xnsq + 3.D0/4.D0*
     &    t1**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0/4.D0*t1**(-1)*
     &    vlpm**2*xnsq**2 + 3.D0*t1**(-1)*vdmp*omro**(-1)*xnsq**2 - 3.D0
     &    *t1**(-1)*vdmp*xnsq**2 + 1.D0/4.D0*t1**(-1)*omro**(-1)*pi**2*
     &    xnsq**2 - 1.D0/4.D0*t1**(-1)*pi**2*xnsq**2 - 3.D0/4.D0*
     &    t2**(-1)*ro*vlpm**2*xnsq**2 - 3.D0*t2**(-1)*ro*vdmp*xnsq**2
     &     - 1.D0/4.D0*t2**(-1)*ro*pi**2*xnsq**2 - 3.D0/2.D0*t2**(-1)*
     &    ro**2*vlpm*vlsm + 3.D0*t2**(-1)*ro**2*vlpm*vlwm - 3.D0/2.D0*
     &    t2**(-1)*ro**2*vlpm + 3.D0/2.D0*t2**(-1)*ro**2*vlpm**2*
     &    xnsq**2 + 3.D0/4.D0*t2**(-1)*ro**2*vlpm**2 + 6.D0*t2**(-1)*
     &    ro**2*vdmp*xnsq**2 + 3.D0*t2**(-1)*ro**2*vdmp + 1.D0/2.D0*
     &    t2**(-1)*ro**2*pi**2*xnsq**2 + 1.D0/4.D0*t2**(-1)*ro**2*pi**2
     &     + 3.D0/4.D0*t2**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0/4.D0*
     &    t2**(-1)*vlpm**2*xnsq**2 + 3.D0*t2**(-1)*vdmp*omro**(-1)*
     &    xnsq**2 )
      gcoeff7 = gcoeff7 + INVG * (  - 3.D0*t2**(-1)*vdmp*xnsq**2 + 1.D0/
     &    4.D0*t2**(-1)*omro**(-1)*pi**2*xnsq**2 - 1.D0/4.D0*t2**(-1)*
     &    pi**2*xnsq**2 + 9.D0*ro*vlpm*vlsm*xnsq + 9.D0*ro*vlpm*vlsm - 
     &    18.D0*ro*vlpm*vltm*xnsq - 12.D0*ro*vlpm*vltm - 6.D0*ro*vlpm*
     &    vlwm + 6.D0*ro*vlpm*xnsq + 6.D0*ro*vlpm - 9.D0/2.D0*ro*
     &    vlpm**2*xnsq - 3.D0/2.D0*ro*vlpm**2*xnsq**2 - 9.D0/2.D0*ro*
     &    vlpm**2 - 18.D0*ro*vdmp*xnsq - 6.D0*ro*vdmp*xnsq**2 - 18.D0*
     &    ro*vdmp - 3.D0/2.D0*ro*pi**2*xnsq - 1.D0/2.D0*ro*pi**2*
     &    xnsq**2 - 3.D0/2.D0*ro*pi**2 - 15.D0/2.D0*ro**2*vlpm*vlsm*
     &    xnsq + 12.D0*ro**2*vlpm*vltm*xnsq + 3.D0*ro**2*vlpm*vltm + 3.D
     &    0*ro**2*vlpm*vlwm*xnsq - 3.D0*ro**2*vlpm*vlwm - 6.D0*ro**2*
     &    vlpm*xnsq + 15.D0/4.D0*ro**2*vlpm**2*xnsq + 15.D0*ro**2*vdmp*
     &    xnsq + 5.D0/4.D0*ro**2*pi**2*xnsq - 3.D0*vlpm**2*omro**(-1)*
     &    xnsq**2 + 3.D0*vlpm**2*xnsq**2 - 12.D0*vdmp*omro**(-1)*
     &    xnsq**2 + 12.D0*vdmp*xnsq**2 - omro**(-1)*pi**2*xnsq**2 + 
     &    pi**2*xnsq**2 )
      gcoeff7 = gcoeff7 + INVG**2 * ( 3.D0/8.D0*b*t1**(-1)*ro**2*
     &    vlpm**2*xnsq + 3.D0/8.D0*b*t1**(-1)*ro**2*vlpm**2 + 3.D0/2.D0
     &    *b*t1**(-1)*ro**2*vltm**2*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**2*
     &    vltm**2 + 3.D0/2.D0*b*t1**(-1)*ro**2*vdt*xnsq + 3.D0/2.D0*b*
     &    t1**(-1)*ro**2*vdt - 1.D0/8.D0*b*t1**(-1)*ro**2*pi**2*xnsq - 
     &    1.D0/8.D0*b*t1**(-1)*ro**2*pi**2 - 3.D0/8.D0*b*t1**(-1)*ro**3
     &    *vlpm**2*xnsq - 3.D0/16.D0*b*t1**(-1)*ro**3*vlpm**2 - 3.D0/4.D
     &    0*b*t1**(-1)*ro**3*vlsm*vltm*xnsq**2 - 3.D0/4.D0*b*t1**(-1)*
     &    ro**3*vlsm*vlwm*xnsq**2 + 3.D0/8.D0*b*t1**(-1)*ro**3*vlsm**2*
     &    xnsq**2 + 9.D0/8.D0*b*t1**(-1)*ro**3*vltm*vlwm*xnsq - 9.D0/16.
     &    D0*b*t1**(-1)*ro**3*vltm**2*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**3
     &    *vltm**2 + 3.D0/16.D0*b*t1**(-1)*ro**3*vlwm**2*xnsq + 3.D0/4.D
     &    0*b*t1**(-1)*ro**3*vdw*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**3*vdw*
     &    xnsq**2 - 3.D0/4.D0*b*t1**(-1)*ro**3*vdt*xnsq**2 - 3.D0/4.D0*
     &    b*t1**(-1)*ro**3*vdt - 1.D0/16.D0*b*t1**(-1)*ro**3*pi**2*xnsq
     &     + 1.D0/8.D0*b*t1**(-1)*ro**3*pi**2*xnsq**2 )
      gcoeff7 = gcoeff7 + INVG**2 * ( 1.D0/16.D0*b*t1**(-1)*ro**3*pi**2
     &     + 3.D0/16.D0*b*t2**(-1)*ro**3*vlpm**2 + 3.D0/4.D0*b*t2**(-1)
     &    *ro**3*vlwm**2 + 3.D0/4.D0*b*t2**(-1)*ro**3*vdw - 1.D0/16.D0*
     &    b*t2**(-1)*ro**3*pi**2 - 3.D0/2.D0*b*ro*vlpm**2*xnsq - 3.D0/2.
     &    D0*b*ro*vlpm**2 - 6.D0*b*ro*vltm**2*xnsq - 6.D0*b*ro*vltm**2
     &     - 6.D0*b*ro*vdt*xnsq - 6.D0*b*ro*vdt + 1.D0/2.D0*b*ro*pi**2*
     &    xnsq + 1.D0/2.D0*b*ro*pi**2 + 3.D0/2.D0*b*ro**2*vlpm**2*xnsq
     &     + 3.D0/4.D0*b*ro**2*vlpm**2 + 3.D0*b*ro**2*vlsm*vltm*xnsq**2
     &     + 3.D0*b*ro**2*vlsm*vlwm*xnsq**2 - 3.D0/2.D0*b*ro**2*vlsm**2
     &    *xnsq**2 - 9.D0/2.D0*b*ro**2*vltm*vlwm*xnsq + 15.D0/4.D0*b*
     &    ro**2*vltm**2*xnsq + 9.D0/2.D0*b*ro**2*vltm**2 - 9.D0/4.D0*b*
     &    ro**2*vlwm**2*xnsq - 3.D0/2.D0*b*ro**2*vlwm**2 - 9.D0/2.D0*b*
     &    ro**2*vdw*xnsq + 3.D0*b*ro**2*vdw*xnsq**2 - 3.D0/2.D0*b*ro**2
     &    *vdw + 3.D0/2.D0*b*ro**2*vdt*xnsq + 3.D0*b*ro**2*vdt*xnsq**2
     &     + 9.D0/2.D0*b*ro**2*vdt + 1.D0/4.D0*b*ro**2*pi**2*xnsq - 1.D0
     &    /2.D0*b*ro**2*pi**2*xnsq**2 )
      gcoeff7 = gcoeff7 + INVG**2 * (  - 1.D0/4.D0*b*ro**2*pi**2 - 3.D0/
     &    4.D0*b*ro**3*vltm**2*xnsq + 3.D0/4.D0*b*ro**3*vlwm**2*xnsq + 
     &    3.D0/4.D0*b*ro**3*vdw*xnsq - 3.D0/4.D0*b*ro**3*vdt*xnsq + 3.D0
     &    /4.D0*t1**(-1)*ro**2*vlpm*vlsm*xnsq + 3.D0/4.D0*t1**(-1)*
     &    ro**2*vlpm*vlsm - 3.D0/2.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 3.
     &    D0/2.D0*t1**(-1)*ro**2*vlpm*vltm - 3.D0/8.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq - 3.D0/8.D0*t1**(-1)*ro**2*vlpm**2 - 3.D0/2.D0*
     &    t1**(-1)*ro**2*vdmp*xnsq - 3.D0/2.D0*t1**(-1)*ro**2*vdmp - 1.D
     &    0/8.D0*t1**(-1)*ro**2*pi**2*xnsq - 1.D0/8.D0*t1**(-1)*ro**2*
     &    pi**2 - 9.D0/8.D0*t1**(-1)*ro**3*vlpm*vlsm*xnsq - 3.D0/4.D0*
     &    t1**(-1)*ro**3*vlpm*vlsm + 15.D0/8.D0*t1**(-1)*ro**3*vlpm*
     &    vltm*xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm*vltm + 3.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm*vlwm*xnsq + 9.D0/16.D0*t1**(-1)*ro**3*
     &    vlpm**2*xnsq + 3.D0/8.D0*t1**(-1)*ro**3*vlpm**2*xnsq**2 + 3.D0
     &    /8.D0*t1**(-1)*ro**3*vlpm**2 + 9.D0/4.D0*t1**(-1)*ro**3*vdmp*
     &    xnsq )
      gcoeff7 = gcoeff7 + INVG**2 * ( 3.D0/2.D0*t1**(-1)*ro**3*vdmp*
     &    xnsq**2 + 3.D0/2.D0*t1**(-1)*ro**3*vdmp + 3.D0/16.D0*t1**(-1)
     &    *ro**3*pi**2*xnsq + 1.D0/8.D0*t1**(-1)*ro**3*pi**2*xnsq**2 + 
     &    1.D0/8.D0*t1**(-1)*ro**3*pi**2 + 3.D0/8.D0*t1**(-1)*ro**4*
     &    vlpm*vlsm*xnsq - 3.D0/8.D0*t1**(-1)*ro**4*vlpm*vltm*xnsq - 3.D
     &    0/8.D0*t1**(-1)*ro**4*vlpm*vlwm*xnsq - 3.D0/16.D0*t1**(-1)*
     &    ro**4*vlpm**2*xnsq - 3.D0/4.D0*t1**(-1)*ro**4*vdmp*xnsq - 1.D0
     &    /16.D0*t1**(-1)*ro**4*pi**2*xnsq + 3.D0/8.D0*t2**(-1)*ro**3*
     &    vlpm**2*xnsq**2 + 3.D0/2.D0*t2**(-1)*ro**3*vdmp*xnsq**2 + 1.D0
     &    /8.D0*t2**(-1)*ro**3*pi**2*xnsq**2 - 3.D0*ro*vlpm*vlsm*xnsq
     &     - 3.D0*ro*vlpm*vlsm + 6.D0*ro*vlpm*vltm*xnsq + 6.D0*ro*vlpm*
     &    vltm + 3.D0/2.D0*ro*vlpm**2*xnsq + 3.D0/2.D0*ro*vlpm**2 + 6.D0
     &    *ro*vdmp*xnsq + 6.D0*ro*vdmp + 1.D0/2.D0*ro*pi**2*xnsq + 1.D0/
     &    2.D0*ro*pi**2 + 9.D0/2.D0*ro**2*vlpm*vlsm*xnsq + 3.D0*ro**2*
     &    vlpm*vlsm - 9.D0*ro**2*vlpm*vltm*xnsq - 15.D0/2.D0*ro**2*vlpm
     &    *vltm )
      gcoeff7 = gcoeff7 + INVG**2 * ( 3.D0/2.D0*ro**2*vlpm*vlwm - 9.D0/
     &    4.D0*ro**2*vlpm**2*xnsq - 3.D0/2.D0*ro**2*vlpm**2*xnsq**2 - 3.
     &    D0/2.D0*ro**2*vlpm**2 - 9.D0*ro**2*vdmp*xnsq - 6.D0*ro**2*
     &    vdmp*xnsq**2 - 6.D0*ro**2*vdmp - 3.D0/4.D0*ro**2*pi**2*xnsq
     &     - 1.D0/2.D0*ro**2*pi**2*xnsq**2 - 1.D0/2.D0*ro**2*pi**2 - 3.D
     &    0/2.D0*ro**3*vlpm*vlsm*xnsq + 3.D0*ro**3*vlpm*vltm*xnsq + 3.D0
     &    /2.D0*ro**3*vlpm*vltm - 3.D0/2.D0*ro**3*vlpm*vlwm + 3.D0/4.D0
     &    *ro**3*vlpm**2*xnsq + 3.D0*ro**3*vdmp*xnsq + 1.D0/4.D0*ro**3*
     &    pi**2*xnsq )
      gcoeff7 = gcoeff7 - 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2*xnsq - 6.D0*b
     &    *t1**(-1)*ro*vlsm*vltm*xnsq**2 - 6.D0*b*t1**(-1)*ro*vlsm*vlwm
     &    *xnsq**2 + 12.D0*b*t1**(-1)*ro*vlsm*xnsq**2 + 3.D0*b*t1**(-1)
     &    *ro*vlsm**2*xnsq**2 + 6.D0*b*t1**(-1)*ro*vltm*xnsq - 12.D0*b*
     &    t1**(-1)*ro*vltm*xnsq**2 - 12.D0*b*t1**(-1)*ro*vltm**2*xnsq
     &     + 6.D0*b*t1**(-1)*ro*vlwm*xnsq - 12.D0*b*t1**(-1)*ro*vlwm*
     &    xnsq**2 + 6.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 6.D0*b*t1**(-1)*
     &    ro*vdw*xnsq - 6.D0*b*t1**(-1)*ro*vdw*xnsq**2 - 12.D0*b*
     &    t1**(-1)*ro*vdt*xnsq - 6.D0*b*t1**(-1)*ro*vdt*xnsq**2 + 1.D0/
     &    2.D0*b*t1**(-1)*ro*pi**2*xnsq + b*t1**(-1)*ro*pi**2*xnsq**2
     &     - 12.D0*b*t1**(-1)*vlsm*omro**(-1)*xnsq**2 + 12.D0*b*
     &    t1**(-1)*vlsm*xnsq**2 - 72.D0*b*t1**(-1)*vlwm*xnsq + 48.D0*b*
     &    t1**(-1)*vlwm*xnsq**2 + 24.D0*b*t1**(-1)*vlwm - 12.D0*b*
     &    t2**(-1)*vlsm*omro**(-1)*xnsq**2 + 12.D0*b*t2**(-1)*vlsm*
     &    xnsq**2 + 24.D0*b*t2**(-1)*vlwm*xnsq - 24.D0*b*t2**(-1)*vlwm*
     &    xnsq**2
      gcoeff7 = gcoeff7 + 24.D0*b*t2**(-1)*xnsq - 12.D0*b*t2**(-1)*
     &    xnsq**2 - 12.D0*b*t2**(-1) + 24.D0*b*vlsm*omro**(-1)*xnsq**2
     &     - 24.D0*b*vlsm*xnsq**2 - 24.D0*b*vlwm*xnsq + 72.D0*b*vlwm*
     &    xnsq**2 - 12.D0*b*xnsq + 12.D0*b*xnsq**2 - 3.D0*t1**(-1)*ro*
     &    vlpm*vlsm*xnsq - 6.D0*t1**(-1)*ro*vlpm*vlsm + 12.D0*t1**(-1)*
     &    ro*vlpm*vltm*xnsq + 12.D0*t1**(-1)*ro*vlpm*vltm - 6.D0*
     &    t1**(-1)*ro*vlpm*vlwm*xnsq - 6.D0*t1**(-1)*ro*vlpm*xnsq - 6.D0
     &    *t1**(-1)*ro*vlpm + 3.D0/2.D0*t1**(-1)*ro*vlpm**2*xnsq + 3.D0
     &    *t1**(-1)*ro*vlpm**2*xnsq**2 + 3.D0*t1**(-1)*ro*vlpm**2 + 6.D0
     &    *t1**(-1)*ro*vdmp*xnsq + 12.D0*t1**(-1)*ro*vdmp*xnsq**2 + 12.D
     &    0*t1**(-1)*ro*vdmp + 1.D0/2.D0*t1**(-1)*ro*pi**2*xnsq + 
     &    t1**(-1)*ro*pi**2*xnsq**2 + t1**(-1)*ro*pi**2 + 3.D0*t1**(-1)
     &    *ro**2*vlpm*vlsm*xnsq - 3.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 
     &    3.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 6.D0*t1**(-1)*ro**2*vlpm
     &    *xnsq - 3.D0/2.D0*t1**(-1)*ro**2*vlpm**2*xnsq - 6.D0*t1**(-1)
     &    *ro**2*vdmp*xnsq
      gcoeff7 = gcoeff7 - 1.D0/2.D0*t1**(-1)*ro**2*pi**2*xnsq + 3.D0*
     &    t1**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0*t1**(-1)*vlpm**2*
     &    xnsq**2 + 12.D0*t1**(-1)*vdmp*omro**(-1)*xnsq**2 - 12.D0*
     &    t1**(-1)*vdmp*xnsq**2 + t1**(-1)*omro**(-1)*pi**2*xnsq**2 - 
     &    t1**(-1)*pi**2*xnsq**2 - 6.D0*t2**(-1)*ro*vlpm*vlsm + 12.D0*
     &    t2**(-1)*ro*vlpm*vlwm - 6.D0*t2**(-1)*ro*vlpm + 3.D0*t2**(-1)
     &    *ro*vlpm**2 + 12.D0*t2**(-1)*ro*vdmp + t2**(-1)*ro*pi**2 + 3.D
     &    0*t2**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0*t2**(-1)*vlpm**2
     &    *xnsq**2 + 12.D0*t2**(-1)*vdmp*omro**(-1)*xnsq**2 - 12.D0*
     &    t2**(-1)*vdmp*xnsq**2 + t2**(-1)*omro**(-1)*pi**2*xnsq**2 - 
     &    t2**(-1)*pi**2*xnsq**2 - 12.D0*ro*vlpm*vlsm*xnsq + 12.D0*ro*
     &    vlpm*vltm*xnsq + 12.D0*ro*vlpm*vlwm*xnsq - 12.D0*ro*vlpm*xnsq
     &     + 6.D0*ro*vlpm**2*xnsq + 24.D0*ro*vdmp*xnsq + 2.D0*ro*pi**2*
     &    xnsq - 6.D0*vlpm**2*omro**(-1)*xnsq**2 + 6.D0*vlpm**2*xnsq**2
     &     - 24.D0*vdmp*omro**(-1)*xnsq**2 + 24.D0*vdmp*xnsq**2 - 2.D0*
     &    omro**(-1)*pi**2*xnsq**2
      gcoeff7 = gcoeff7 + 2.D0*pi**2*xnsq**2

      return
      end
