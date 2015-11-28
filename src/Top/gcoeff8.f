      double precision function gcoeff8()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff8 =  + TBAR**(-1) * (  - 24.D0*b*t1*vltm*xnsq + 24.D0*b*t1*
     &    vltm*xnsq**2 + 48.D0*b*vltm*xnsq - 24.D0*b*vltm*xnsq**2 - 24.D
     &    0*b*vltm )
      gcoeff8 = gcoeff8 + UBAR**(-1) * ( 24.D0*b*t1*vlwm*xnsq - 24.D0*b
     &    *t1*vlwm*xnsq**2 + 24.D0*b*vlwm*xnsq - 24.D0*b*vlwm )
      gcoeff8 = gcoeff8 + INVG * (  - 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2
     &     - 6.D0*b*t1**(-1)*ro*vlsm*vlwm*xnsq**2 + 15.D0/4.D0*b*
     &    t1**(-1)*ro*vlsm*xnsq**2 + 3.D0/2.D0*b*t1**(-1)*ro*vlsm**2*
     &    xnsq**2 - 12.D0*b*t1**(-1)*ro*vltm*xnsq + 6.D0*b*t1**(-1)*ro*
     &    vltm*xnsq**2 + 6.D0*b*t1**(-1)*ro*vltm - 6.D0*b*t1**(-1)*ro*
     &    vltm**2*xnsq - 6.D0*b*t1**(-1)*ro*vltm**2 + 12.D0*b*t1**(-1)*
     &    ro*vlwm*xnsq - 24.D0*b*t1**(-1)*ro*vlwm*xnsq**2 + 12.D0*b*
     &    t1**(-1)*ro*vlwm + 6.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 6.D0*b*
     &    t1**(-1)*ro*vdw*xnsq - 6.D0*b*t1**(-1)*ro*vdw*xnsq**2 - 6.D0*
     &    b*t1**(-1)*ro*vdt*xnsq - 6.D0*b*t1**(-1)*ro*vdt + 1.D0/2.D0*b
     &    *t1**(-1)*ro*pi**2*xnsq**2 + 1.D0/2.D0*b*t1**(-1)*ro*pi**2 - 
     &    3.D0/2.D0*b*t1**(-1)*ro*xnsq**2 + 3.D0*b*t1**(-1)*ro**2*
     &    vlpm**2 - 3.D0*b*t1**(-1)*ro**2*vlsm*vltm*xnsq**2 + 3.D0/2.D0
     &    *b*t1**(-1)*ro**2*vlsm*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*ro**2*
     &    vlsm**2*xnsq**2 + 12.D0*b*t1**(-1)*ro**2*vltm*xnsq - 9.D0*b*
     &    t1**(-1)*ro**2*vltm*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG * (  - 12.D0*b*t1**(-1)*ro**2*vltm + 15.D
     &    0*b*t1**(-1)*ro**2*vltm**2*xnsq + 12.D0*b*t1**(-1)*ro**2*
     &    vltm**2 - 12.D0*b*t1**(-1)*ro**2*vlwm*xnsq + 6.D0*b*t1**(-1)*
     &    ro**2*vlwm*xnsq**2 - 15.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq - 15.
     &    D0*b*t1**(-1)*ro**2*vdw*xnsq + 15.D0*b*t1**(-1)*ro**2*vdt*
     &    xnsq - 3.D0*b*t1**(-1)*ro**2*vdt*xnsq**2 + 12.D0*b*t1**(-1)*
     &    ro**2*vdt + 1.D0/4.D0*b*t1**(-1)*ro**2*pi**2*xnsq**2 - b*
     &    t1**(-1)*ro**2*pi**2 + 3.D0/2.D0*b*t1**(-1)*ro**2 - 21.D0/4.D0
     &    *b*t1**(-1)*ro**3*vltm**2*xnsq + 21.D0/4.D0*b*t1**(-1)*ro**3*
     &    vlwm**2*xnsq + 21.D0/4.D0*b*t1**(-1)*ro**3*vdw*xnsq - 21.D0/4.
     &    D0*b*t1**(-1)*ro**3*vdt*xnsq - 9.D0/4.D0*b*t1**(-1)*vlsm*
     &    omro**(-2)*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*vlsm*omro**(-1)*
     &    xnsq**2 + 3.D0/2.D0*b*t1**(-1)*vlsm*xnsq**2 + 3.D0/2.D0*b*
     &    t1**(-1)*omro**(-1)*xnsq**2 - 3.D0/2.D0*b*t1**(-1)*xnsq**2 - 
     &    3.D0/2.D0*b*t2**(-1)*ro*vlpm**2 - 6.D0*b*t2**(-1)*ro*vlsm*
     &    vlwm*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG * ( 15.D0/4.D0*b*t2**(-1)*ro*vlsm*
     &    xnsq**2 + 3.D0/2.D0*b*t2**(-1)*ro*vlsm**2*xnsq**2 - 18.D0*b*
     &    t2**(-1)*ro*vlwm*xnsq**2 + 18.D0*b*t2**(-1)*ro*vlwm - 6.D0*b*
     &    t2**(-1)*ro*vlwm**2 - 6.D0*b*t2**(-1)*ro*vdw*xnsq**2 - 6.D0*b
     &    *t2**(-1)*ro*vdw + 1.D0/2.D0*b*t2**(-1)*ro*pi**2*xnsq**2 + 1.D
     &    0/2.D0*b*t2**(-1)*ro*pi**2 - 3.D0/2.D0*b*t2**(-1)*ro*xnsq**2
     &     + 3.D0*b*t2**(-1)*ro**2*vlpm**2 - 3.D0*b*t2**(-1)*ro**2*vlsm
     &    *vlwm*xnsq**2 + 3.D0/2.D0*b*t2**(-1)*ro**2*vlsm*xnsq**2 + 3.D0
     &    /4.D0*b*t2**(-1)*ro**2*vlsm**2*xnsq**2 - 3.D0*b*t2**(-1)*
     &    ro**2*vlwm*xnsq**2 - 12.D0*b*t2**(-1)*ro**2*vlwm + 12.D0*b*
     &    t2**(-1)*ro**2*vlwm**2 - 3.D0*b*t2**(-1)*ro**2*vdw*xnsq**2 + 
     &    12.D0*b*t2**(-1)*ro**2*vdw + 1.D0/4.D0*b*t2**(-1)*ro**2*pi**2
     &    *xnsq**2 - b*t2**(-1)*ro**2*pi**2 + 3.D0/2.D0*b*t2**(-1)*
     &    ro**2 - 9.D0/4.D0*b*t2**(-1)*vlsm*omro**(-2)*xnsq**2 + 3.D0/4.
     &    D0*b*t2**(-1)*vlsm*omro**(-1)*xnsq**2 + 3.D0/2.D0*b*t2**(-1)*
     &    vlsm*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG * ( 3.D0/2.D0*b*t2**(-1)*omro**(-1)*
     &    xnsq**2 - 3.D0/2.D0*b*t2**(-1)*xnsq**2 - 3.D0/2.D0*b*ro*
     &    vlpm**2*xnsq + 3.D0*b*ro*vlpm**2 + 6.D0*b*ro*vlsm*vltm*
     &    xnsq**2 + 12.D0*b*ro*vlsm*vlwm*xnsq**2 - 9.D0/2.D0*b*ro*vlsm*
     &    xnsq**2 - 9.D0/2.D0*b*ro*vlsm**2*xnsq**2 - 24.D0*b*ro*vltm*
     &    vlwm*xnsq + 18.D0*b*ro*vltm*xnsq - 6.D0*b*ro*vltm*xnsq**2 - 
     &    12.D0*b*ro*vltm**2*xnsq + 6.D0*b*ro*vltm**2 + 54.D0*b*ro*vlwm
     &    *xnsq - 18.D0*b*ro*vlwm*xnsq**2 - 18.D0*b*ro*vlwm**2*xnsq + 6.
     &    D0*b*ro*vlwm**2 - 30.D0*b*ro*vdw*xnsq + 12.D0*b*ro*vdw*
     &    xnsq**2 + 6.D0*b*ro*vdw - 24.D0*b*ro*vdt*xnsq + 6.D0*b*ro*vdt
     &    *xnsq**2 + 6.D0*b*ro*vdt + 9.D0/2.D0*b*ro*pi**2*xnsq - 3.D0/2.
     &    D0*b*ro*pi**2*xnsq**2 - b*ro*pi**2 - 12.D0*b*ro*xnsq + 6.D0*b
     &    *ro*xnsq**2 + 6.D0*b*ro**2*vlpm**2*xnsq + 3.D0*b*ro**2*vlsm*
     &    vltm*xnsq**2 + 3.D0*b*ro**2*vlsm*vlwm*xnsq**2 - 3.D0*b*ro**2*
     &    vlsm*xnsq**2 - 3.D0/2.D0*b*ro**2*vlsm**2*xnsq**2 - 12.D0*b*
     &    ro**2*vltm*xnsq )
      gcoeff8 = gcoeff8 + INVG * ( 3.D0*b*ro**2*vltm*xnsq**2 + 12.D0*b*
     &    ro**2*vltm**2*xnsq - 12.D0*b*ro**2*vlwm*xnsq + 3.D0*b*ro**2*
     &    vlwm*xnsq**2 + 12.D0*b*ro**2*vlwm**2*xnsq + 12.D0*b*ro**2*vdw
     &    *xnsq + 3.D0*b*ro**2*vdw*xnsq**2 + 12.D0*b*ro**2*vdt*xnsq + 3.
     &    D0*b*ro**2*vdt*xnsq**2 - 2.D0*b*ro**2*pi**2*xnsq - 1.D0/2.D0*
     &    b*ro**2*pi**2*xnsq**2 + 3.D0*b*ro**2*xnsq + 9.D0*b*vlsm*
     &    omro**(-2)*xnsq**2 - 15.D0/2.D0*b*vlsm*omro**(-1)*xnsq**2 - 3.
     &    D0/2.D0*b*vlsm*xnsq**2 - 48.D0*b*vlwm*xnsq + 48.D0*b*vlwm*
     &    xnsq**2 - 6.D0*b*omro**(-1)*xnsq**2 + 6.D0*b*xnsq**2 - 3.D0*
     &    t1**(-1)*ro*vlpm*vlsm + 6.D0*t1**(-1)*ro*vlpm*vltm*xnsq + 6.D0
     &    *t1**(-1)*ro*vlpm*vltm - 6.D0*t1**(-1)*ro*vlpm*vlwm*xnsq + 3.D
     &    0/4.D0*t1**(-1)*ro*vlpm**2*xnsq**2 + 3.D0/2.D0*t1**(-1)*ro*
     &    vlpm**2 + 3.D0*t1**(-1)*ro*vdmp*xnsq**2 + 6.D0*t1**(-1)*ro*
     &    vdmp + 1.D0/4.D0*t1**(-1)*ro*pi**2*xnsq**2 + 1.D0/2.D0*
     &    t1**(-1)*ro*pi**2 + 9.D0/2.D0*t1**(-1)*ro**2*vlpm*vlsm - 18.D0
     &    *t1**(-1)*ro**2*vlpm*vltm*xnsq )
      gcoeff8 = gcoeff8 + INVG * (  - 9.D0*t1**(-1)*ro**2*vlpm*vltm + 
     &    18.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 9.D0/2.D0*t1**(-1)*
     &    ro**2*vlpm + 3.D0/16.D0*t1**(-1)*ro**2*vlpm**2*xnsq**2 - 9.D0/
     &    4.D0*t1**(-1)*ro**2*vlpm**2 + 3.D0/4.D0*t1**(-1)*ro**2*vdmp*
     &    xnsq**2 - 9.D0*t1**(-1)*ro**2*vdmp + 1.D0/16.D0*t1**(-1)*
     &    ro**2*pi**2*xnsq**2 - 3.D0/4.D0*t1**(-1)*ro**2*pi**2 - 9.D0/4.
     &    D0*t1**(-1)*ro**3*vlpm*vlsm + 27.D0/2.D0*t1**(-1)*ro**3*vlpm*
     &    vltm*xnsq + 9.D0/2.D0*t1**(-1)*ro**3*vlpm*vltm - 27.D0/2.D0*
     &    t1**(-1)*ro**3*vlpm*vlwm*xnsq - 21.D0/4.D0*t1**(-1)*ro**3*
     &    vlpm + 9.D0/8.D0*t1**(-1)*ro**3*vlpm**2 + 9.D0/2.D0*t1**(-1)*
     &    ro**3*vdmp + 3.D0/8.D0*t1**(-1)*ro**3*pi**2 + 9.D0/16.D0*
     &    t1**(-1)*vlpm**2*omro**(-2)*xnsq**2 - 3.D0/8.D0*t1**(-1)*
     &    vlpm**2*omro**(-1)*xnsq**2 - 3.D0/16.D0*t1**(-1)*vlpm**2*
     &    xnsq**2 + 9.D0/4.D0*t1**(-1)*vdmp*omro**(-2)*xnsq**2 - 3.D0/2.
     &    D0*t1**(-1)*vdmp*omro**(-1)*xnsq**2 - 3.D0/4.D0*t1**(-1)*vdmp
     &    *xnsq**2 )
      gcoeff8 = gcoeff8 + INVG * ( 3.D0/16.D0*t1**(-1)*omro**(-2)*pi**2
     &    *xnsq**2 - 1.D0/8.D0*t1**(-1)*omro**(-1)*pi**2*xnsq**2 - 1.D0/
     &    16.D0*t1**(-1)*pi**2*xnsq**2 - 3.D0*t2**(-1)*ro*vlpm*vlsm + 6.
     &    D0*t2**(-1)*ro*vlpm*vlwm + 3.D0/4.D0*t2**(-1)*ro*vlpm**2*
     &    xnsq**2 + 3.D0/2.D0*t2**(-1)*ro*vlpm**2 + 3.D0*t2**(-1)*ro*
     &    vdmp*xnsq**2 + 6.D0*t2**(-1)*ro*vdmp + 1.D0/4.D0*t2**(-1)*ro*
     &    pi**2*xnsq**2 + 1.D0/2.D0*t2**(-1)*ro*pi**2 + 9.D0/2.D0*
     &    t2**(-1)*ro**2*vlpm*vlsm - 9.D0*t2**(-1)*ro**2*vlpm*vlwm + 9.D
     &    0/2.D0*t2**(-1)*ro**2*vlpm + 3.D0/16.D0*t2**(-1)*ro**2*
     &    vlpm**2*xnsq**2 - 9.D0/4.D0*t2**(-1)*ro**2*vlpm**2 + 3.D0/4.D0
     &    *t2**(-1)*ro**2*vdmp*xnsq**2 - 9.D0*t2**(-1)*ro**2*vdmp + 1.D0
     &    /16.D0*t2**(-1)*ro**2*pi**2*xnsq**2 - 3.D0/4.D0*t2**(-1)*
     &    ro**2*pi**2 - 9.D0/4.D0*t2**(-1)*ro**3*vlpm*vlsm + 9.D0/2.D0*
     &    t2**(-1)*ro**3*vlpm*vlwm - 21.D0/4.D0*t2**(-1)*ro**3*vlpm + 9.
     &    D0/8.D0*t2**(-1)*ro**3*vlpm**2 + 9.D0/2.D0*t2**(-1)*ro**3*
     &    vdmp )
      gcoeff8 = gcoeff8 + INVG * ( 3.D0/8.D0*t2**(-1)*ro**3*pi**2 + 9.D0
     &    /16.D0*t2**(-1)*vlpm**2*omro**(-2)*xnsq**2 - 3.D0/8.D0*
     &    t2**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0/16.D0*t2**(-1)*
     &    vlpm**2*xnsq**2 + 9.D0/4.D0*t2**(-1)*vdmp*omro**(-2)*xnsq**2
     &     - 3.D0/2.D0*t2**(-1)*vdmp*omro**(-1)*xnsq**2 - 3.D0/4.D0*
     &    t2**(-1)*vdmp*xnsq**2 + 3.D0/16.D0*t2**(-1)*omro**(-2)*pi**2*
     &    xnsq**2 - 1.D0/8.D0*t2**(-1)*omro**(-1)*pi**2*xnsq**2 - 1.D0/
     &    16.D0*t2**(-1)*pi**2*xnsq**2 - 3.D0*ro*vlpm*vlsm*xnsq + 6.D0*
     &    ro*vlpm*vlsm - 6.D0*ro*vlpm*vltm + 6.D0*ro*vlpm*vlwm*xnsq - 6.
     &    D0*ro*vlpm*vlwm - 6.D0*ro*vlpm*xnsq - 12.D0*ro*vlpm + 3.D0/2.D
     &    0*ro*vlpm**2*xnsq - 21.D0/8.D0*ro*vlpm**2*xnsq**2 - 3.D0*ro*
     &    vlpm**2 + 6.D0*ro*vdmp*xnsq - 21.D0/2.D0*ro*vdmp*xnsq**2 - 12.
     &    D0*ro*vdmp + 1.D0/2.D0*ro*pi**2*xnsq - 7.D0/8.D0*ro*pi**2*
     &    xnsq**2 - ro*pi**2 + 21.D0/2.D0*ro**2*vlpm*vlsm*xnsq + 3.D0*
     &    ro**2*vlpm*vlsm - 12.D0*ro**2*vlpm*vltm*xnsq - 3.D0*ro**2*
     &    vlpm*vltm )
      gcoeff8 = gcoeff8 + INVG * (  - 9.D0*ro**2*vlpm*vlwm*xnsq - 3.D0*
     &    ro**2*vlpm*vlwm + 39.D0/2.D0*ro**2*vlpm*xnsq + 21.D0*ro**2*
     &    vlpm - 21.D0/4.D0*ro**2*vlpm**2*xnsq + 3.D0/8.D0*ro**2*
     &    vlpm**2*xnsq**2 - 3.D0/2.D0*ro**2*vlpm**2 - 21.D0*ro**2*vdmp*
     &    xnsq + 3.D0/2.D0*ro**2*vdmp*xnsq**2 - 6.D0*ro**2*vdmp - 7.D0/
     &    4.D0*ro**2*pi**2*xnsq + 1.D0/8.D0*ro**2*pi**2*xnsq**2 - 1.D0/
     &    2.D0*ro**2*pi**2 - 9.D0/2.D0*ro**3*vlpm*vlsm*xnsq + 9.D0/2.D0
     &    *ro**3*vlpm*vltm*xnsq + 9.D0/2.D0*ro**3*vlpm*vlwm*xnsq - 21.D0
     &    /2.D0*ro**3*vlpm*xnsq + 9.D0/4.D0*ro**3*vlpm**2*xnsq + 9.D0*
     &    ro**3*vdmp*xnsq + 3.D0/4.D0*ro**3*pi**2*xnsq - 9.D0/4.D0*
     &    vlpm**2*omro**(-2)*xnsq**2 + 21.D0/8.D0*vlpm**2*omro**(-1)*
     &    xnsq**2 - 3.D0/8.D0*vlpm**2*xnsq**2 - 9.D0*vdmp*omro**(-2)*
     &    xnsq**2 + 21.D0/2.D0*vdmp*omro**(-1)*xnsq**2 - 3.D0/2.D0*vdmp
     &    *xnsq**2 - 3.D0/4.D0*omro**(-2)*pi**2*xnsq**2 + 7.D0/8.D0*
     &    omro**(-1)*pi**2*xnsq**2 - 1.D0/8.D0*pi**2*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**2 * ( 3.D0/4.D0*b*t1**(-1)*ro*vlsm*
     &    xnsq**2 - 9.D0/8.D0*b*t1**(-1)*ro**2*vlpm**2 - 9.D0/2.D0*b*
     &    t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*ro**2
     &    *vlsm*xnsq**2 + 9.D0/8.D0*b*t1**(-1)*ro**2*vlsm**2*xnsq**2 - 
     &    9.D0/2.D0*b*t1**(-1)*ro**2*vltm**2*xnsq - 9.D0/2.D0*b*
     &    t1**(-1)*ro**2*vltm**2 - 3.D0*b*t1**(-1)*ro**2*vlwm*xnsq**2
     &     + 3.D0*b*t1**(-1)*ro**2*vlwm + 9.D0/2.D0*b*t1**(-1)*ro**2*
     &    vlwm**2*xnsq + 9.D0/2.D0*b*t1**(-1)*ro**2*vdw*xnsq - 9.D0/2.D0
     &    *b*t1**(-1)*ro**2*vdw*xnsq**2 - 9.D0/2.D0*b*t1**(-1)*ro**2*
     &    vdt*xnsq - 9.D0/2.D0*b*t1**(-1)*ro**2*vdt + 3.D0/8.D0*b*
     &    t1**(-1)*ro**2*pi**2*xnsq**2 + 3.D0/8.D0*b*t1**(-1)*ro**2*
     &    pi**2 + 3.D0/2.D0*b*t1**(-1)*ro**3*vlpm**2 - 15.D0/8.D0*b*
     &    t1**(-1)*ro**3*vlsm*vltm*xnsq**2 + 9.D0/8.D0*b*t1**(-1)*ro**3
     &    *vlsm*vlwm*xnsq**2 + 3.D0/16.D0*b*t1**(-1)*ro**3*vlsm**2*
     &    xnsq**2 + 3.D0/4.D0*b*t1**(-1)*ro**3*vltm*xnsq - 3.D0/4.D0*b*
     &    t1**(-1)*ro**3*vltm*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**2 * (  - 3.D0/2.D0*b*t1**(-1)*ro**3*
     &    vltm + 63.D0/8.D0*b*t1**(-1)*ro**3*vltm**2*xnsq + 6.D0*b*
     &    t1**(-1)*ro**3*vltm**2 - 3.D0/4.D0*b*t1**(-1)*ro**3*vlwm*xnsq
     &     + 3.D0/4.D0*b*t1**(-1)*ro**3*vlwm*xnsq**2 - 63.D0/8.D0*b*
     &    t1**(-1)*ro**3*vlwm**2*xnsq - 63.D0/8.D0*b*t1**(-1)*ro**3*vdw
     &    *xnsq + 9.D0/8.D0*b*t1**(-1)*ro**3*vdw*xnsq**2 + 63.D0/8.D0*b
     &    *t1**(-1)*ro**3*vdt*xnsq - 15.D0/8.D0*b*t1**(-1)*ro**3*vdt*
     &    xnsq**2 + 6.D0*b*t1**(-1)*ro**3*vdt + 1.D0/16.D0*b*t1**(-1)*
     &    ro**3*pi**2*xnsq**2 - 1.D0/2.D0*b*t1**(-1)*ro**3*pi**2 - 33.D0
     &    /16.D0*b*t1**(-1)*ro**4*vltm**2*xnsq + 33.D0/16.D0*b*t1**(-1)
     &    *ro**4*vlwm**2*xnsq + 33.D0/16.D0*b*t1**(-1)*ro**4*vdw*xnsq
     &     - 33.D0/16.D0*b*t1**(-1)*ro**4*vdt*xnsq - 3.D0/4.D0*b*
     &    t1**(-1)*vlsm*omro**(-1)*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*vlsm*
     &    xnsq**2 + 3.D0/4.D0*b*t2**(-1)*ro*vlsm*xnsq**2 - 9.D0/8.D0*b*
     &    t2**(-1)*ro**2*vlpm**2 - 9.D0/2.D0*b*t2**(-1)*ro**2*vlsm*vlwm
     &    *xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**2 * ( 3.D0/4.D0*b*t2**(-1)*ro**2*vlsm*
     &    xnsq**2 + 9.D0/8.D0*b*t2**(-1)*ro**2*vlsm**2*xnsq**2 - 3.D0*b
     &    *t2**(-1)*ro**2*vlwm*xnsq**2 + 3.D0*b*t2**(-1)*ro**2*vlwm - 9.
     &    D0/2.D0*b*t2**(-1)*ro**2*vlwm**2 - 9.D0/2.D0*b*t2**(-1)*ro**2
     &    *vdw*xnsq**2 - 9.D0/2.D0*b*t2**(-1)*ro**2*vdw + 3.D0/8.D0*b*
     &    t2**(-1)*ro**2*pi**2*xnsq**2 + 3.D0/8.D0*b*t2**(-1)*ro**2*
     &    pi**2 + 3.D0/2.D0*b*t2**(-1)*ro**3*vlpm**2 - 3.D0/4.D0*b*
     &    t2**(-1)*ro**3*vlsm*vlwm*xnsq**2 + 3.D0/16.D0*b*t2**(-1)*
     &    ro**3*vlsm**2*xnsq**2 - 3.D0/2.D0*b*t2**(-1)*ro**3*vlwm + 6.D0
     &    *b*t2**(-1)*ro**3*vlwm**2 - 3.D0/4.D0*b*t2**(-1)*ro**3*vdw*
     &    xnsq**2 + 6.D0*b*t2**(-1)*ro**3*vdw + 1.D0/16.D0*b*t2**(-1)*
     &    ro**3*pi**2*xnsq**2 - 1.D0/2.D0*b*t2**(-1)*ro**3*pi**2 - 3.D0/
     &    4.D0*b*t2**(-1)*vlsm*omro**(-1)*xnsq**2 + 3.D0/4.D0*b*
     &    t2**(-1)*vlsm*xnsq**2 + 3.D0/2.D0*b*ro*vlpm**2*xnsq + 3.D0*b*
     &    ro*vlpm**2 + 6.D0*b*ro*vlsm*vlwm*xnsq**2 - 3.D0*b*ro*vlsm*
     &    xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**2 * (  - 3.D0/2.D0*b*ro*vlsm**2*xnsq**2
     &     + 6.D0*b*ro*vltm**2*xnsq + 6.D0*b*ro*vltm**2 + 12.D0*b*ro*
     &    vlwm*xnsq**2 - 12.D0*b*ro*vlwm + 6.D0*b*ro*vlwm**2 + 6.D0*b*
     &    ro*vdw*xnsq**2 + 6.D0*b*ro*vdw + 6.D0*b*ro*vdt*xnsq + 6.D0*b*
     &    ro*vdt - 1.D0/2.D0*b*ro*pi**2*xnsq - 1.D0/2.D0*b*ro*pi**2*
     &    xnsq**2 - b*ro*pi**2 - 51.D0/8.D0*b*ro**2*vlpm**2*xnsq - 33.D0
     &    /4.D0*b*ro**2*vlpm**2 + 9.D0/2.D0*b*ro**2*vlsm*vltm*xnsq**2
     &     + 3.D0*b*ro**2*vlsm*vlwm*xnsq**2 - 15.D0/8.D0*b*ro**2*
     &    vlsm**2*xnsq**2 + 3.D0*b*ro**2*vltm*xnsq + 3.D0*b*ro**2*vltm
     &     - 21.D0*b*ro**2*vltm**2*xnsq - 33.D0/2.D0*b*ro**2*vltm**2 + 
     &    6.D0*b*ro**2*vlwm*xnsq - 3.D0*b*ro**2*vlwm*xnsq**2 + 3.D0*b*
     &    ro**2*vlwm - 9.D0/2.D0*b*ro**2*vlwm**2*xnsq - 33.D0/2.D0*b*
     &    ro**2*vlwm**2 - 9.D0/2.D0*b*ro**2*vdw*xnsq + 3.D0*b*ro**2*vdw
     &    *xnsq**2 - 33.D0/2.D0*b*ro**2*vdw - 21.D0*b*ro**2*vdt*xnsq + 
     &    9.D0/2.D0*b*ro**2*vdt*xnsq**2 - 33.D0/2.D0*b*ro**2*vdt + 17.D0
     &    /8.D0*b*ro**2*pi**2*xnsq )
      gcoeff8 = gcoeff8 + INVG**2 * (  - 5.D0/8.D0*b*ro**2*pi**2*
     &    xnsq**2 + 11.D0/4.D0*b*ro**2*pi**2 + 69.D0/16.D0*b*ro**3*
     &    vlpm**2*xnsq + 21.D0/8.D0*b*ro**3*vlpm**2 + 3.D0/4.D0*b*ro**3
     &    *vlsm*vltm*xnsq**2 + 3.D0/4.D0*b*ro**3*vlsm*vlwm*xnsq**2 - 3.D
     &    0/8.D0*b*ro**3*vlsm**2*xnsq**2 - 15.D0/4.D0*b*ro**3*vltm*vlwm
     &    *xnsq - 3.D0/2.D0*b*ro**3*vltm*xnsq + 75.D0/8.D0*b*ro**3*
     &    vltm**2*xnsq + 21.D0/4.D0*b*ro**3*vltm**2 - 3.D0/2.D0*b*ro**3
     &    *vlwm*xnsq + 33.D0/8.D0*b*ro**3*vlwm**2*xnsq + 21.D0/4.D0*b*
     &    ro**3*vlwm**2 + 9.D0/4.D0*b*ro**3*vdw*xnsq + 3.D0/4.D0*b*
     &    ro**3*vdw*xnsq**2 + 21.D0/4.D0*b*ro**3*vdw + 15.D0/2.D0*b*
     &    ro**3*vdt*xnsq + 3.D0/4.D0*b*ro**3*vdt*xnsq**2 + 21.D0/4.D0*b
     &    *ro**3*vdt - 13.D0/16.D0*b*ro**3*pi**2*xnsq - 1.D0/8.D0*b*
     &    ro**3*pi**2*xnsq**2 - 7.D0/8.D0*b*ro**3*pi**2 + 3.D0*b*vlsm*
     &    omro**(-1)*xnsq**2 - 3.D0*b*vlsm*xnsq**2 - 3.D0/16.D0*
     &    t1**(-1)*ro*vlpm**2*xnsq**2 - 3.D0/4.D0*t1**(-1)*ro*vdmp*
     &    xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**2 * (  - 1.D0/16.D0*t1**(-1)*ro*pi**2*
     &    xnsq**2 - 9.D0/4.D0*t1**(-1)*ro**2*vlpm*vlsm + 9.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm*vltm*xnsq + 9.D0/2.D0*t1**(-1)*ro**2*vlpm
     &    *vltm - 9.D0/2.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 15.D0/16.D0
     &    *t1**(-1)*ro**2*vlpm**2*xnsq**2 + 9.D0/8.D0*t1**(-1)*ro**2*
     &    vlpm**2 + 15.D0/4.D0*t1**(-1)*ro**2*vdmp*xnsq**2 + 9.D0/2.D0*
     &    t1**(-1)*ro**2*vdmp + 5.D0/16.D0*t1**(-1)*ro**2*pi**2*xnsq**2
     &     + 3.D0/8.D0*t1**(-1)*ro**2*pi**2 + 33.D0/8.D0*t1**(-1)*ro**3
     &    *vlpm*vlsm - 81.D0/8.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq - 33.D0/
     &    4.D0*t1**(-1)*ro**3*vlpm*vltm + 81.D0/8.D0*t1**(-1)*ro**3*
     &    vlpm*vlwm*xnsq + 3.D0/4.D0*t1**(-1)*ro**3*vlpm - 9.D0/16.D0*
     &    t1**(-1)*ro**3*vlpm**2*xnsq**2 - 33.D0/16.D0*t1**(-1)*ro**3*
     &    vlpm**2 - 9.D0/4.D0*t1**(-1)*ro**3*vdmp*xnsq**2 - 33.D0/4.D0*
     &    t1**(-1)*ro**3*vdmp - 3.D0/16.D0*t1**(-1)*ro**3*pi**2*xnsq**2
     &     - 11.D0/16.D0*t1**(-1)*ro**3*pi**2 - 27.D0/16.D0*t1**(-1)*
     &    ro**4*vlpm*vlsm )
      gcoeff8 = gcoeff8 + INVG**2 * ( 87.D0/16.D0*t1**(-1)*ro**4*vlpm*
     &    vltm*xnsq + 27.D0/8.D0*t1**(-1)*ro**4*vlpm*vltm - 87.D0/16.D0
     &    *t1**(-1)*ro**4*vlpm*vlwm*xnsq - 3.D0/4.D0*t1**(-1)*ro**4*
     &    vlpm + 27.D0/32.D0*t1**(-1)*ro**4*vlpm**2 + 27.D0/8.D0*
     &    t1**(-1)*ro**4*vdmp + 9.D0/32.D0*t1**(-1)*ro**4*pi**2 + 3.D0/
     &    16.D0*t1**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0/16.D0*
     &    t1**(-1)*vlpm**2*xnsq**2 + 3.D0/4.D0*t1**(-1)*vdmp*omro**(-1)
     &    *xnsq**2 - 3.D0/4.D0*t1**(-1)*vdmp*xnsq**2 + 1.D0/16.D0*
     &    t1**(-1)*omro**(-1)*pi**2*xnsq**2 - 1.D0/16.D0*t1**(-1)*pi**2
     &    *xnsq**2 - 3.D0/16.D0*t2**(-1)*ro*vlpm**2*xnsq**2 - 3.D0/4.D0
     &    *t2**(-1)*ro*vdmp*xnsq**2 - 1.D0/16.D0*t2**(-1)*ro*pi**2*
     &    xnsq**2 - 9.D0/4.D0*t2**(-1)*ro**2*vlpm*vlsm + 9.D0/2.D0*
     &    t2**(-1)*ro**2*vlpm*vlwm + 15.D0/16.D0*t2**(-1)*ro**2*vlpm**2
     &    *xnsq**2 + 9.D0/8.D0*t2**(-1)*ro**2*vlpm**2 + 15.D0/4.D0*
     &    t2**(-1)*ro**2*vdmp*xnsq**2 + 9.D0/2.D0*t2**(-1)*ro**2*vdmp
     &     + 5.D0/16.D0*t2**(-1)*ro**2*pi**2*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**2 * ( 3.D0/8.D0*t2**(-1)*ro**2*pi**2 + 
     &    33.D0/8.D0*t2**(-1)*ro**3*vlpm*vlsm - 33.D0/4.D0*t2**(-1)*
     &    ro**3*vlpm*vlwm + 3.D0/4.D0*t2**(-1)*ro**3*vlpm - 9.D0/16.D0*
     &    t2**(-1)*ro**3*vlpm**2*xnsq**2 - 33.D0/16.D0*t2**(-1)*ro**3*
     &    vlpm**2 - 9.D0/4.D0*t2**(-1)*ro**3*vdmp*xnsq**2 - 33.D0/4.D0*
     &    t2**(-1)*ro**3*vdmp - 3.D0/16.D0*t2**(-1)*ro**3*pi**2*xnsq**2
     &     - 11.D0/16.D0*t2**(-1)*ro**3*pi**2 - 27.D0/16.D0*t2**(-1)*
     &    ro**4*vlpm*vlsm + 27.D0/8.D0*t2**(-1)*ro**4*vlpm*vlwm - 3.D0/
     &    4.D0*t2**(-1)*ro**4*vlpm + 27.D0/32.D0*t2**(-1)*ro**4*vlpm**2
     &     + 27.D0/8.D0*t2**(-1)*ro**4*vdmp + 9.D0/32.D0*t2**(-1)*ro**4
     &    *pi**2 + 3.D0/16.D0*t2**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0
     &    /16.D0*t2**(-1)*vlpm**2*xnsq**2 + 3.D0/4.D0*t2**(-1)*vdmp*
     &    omro**(-1)*xnsq**2 - 3.D0/4.D0*t2**(-1)*vdmp*xnsq**2 + 1.D0/
     &    16.D0*t2**(-1)*omro**(-1)*pi**2*xnsq**2 - 1.D0/16.D0*t2**(-1)
     &    *pi**2*xnsq**2 + 3.D0*ro*vlpm*vlsm*xnsq + 6.D0*ro*vlpm*vlsm
     &     - 6.D0*ro*vlpm*vltm*xnsq )
      gcoeff8 = gcoeff8 + INVG**2 * (  - 6.D0*ro*vlpm*vltm - 6.D0*ro*
     &    vlpm*vlwm - 3.D0/2.D0*ro*vlpm**2*xnsq - 3.D0/4.D0*ro*vlpm**2*
     &    xnsq**2 - 3.D0*ro*vlpm**2 - 6.D0*ro*vdmp*xnsq - 3.D0*ro*vdmp*
     &    xnsq**2 - 12.D0*ro*vdmp - 1.D0/2.D0*ro*pi**2*xnsq - 1.D0/4.D0
     &    *ro*pi**2*xnsq**2 - ro*pi**2 - 57.D0/4.D0*ro**2*vlpm*vlsm*
     &    xnsq - 39.D0/2.D0*ro**2*vlpm*vlsm + 24.D0*ro**2*vlpm*vltm*
     &    xnsq + 39.D0/2.D0*ro**2*vlpm*vltm + 9.D0/2.D0*ro**2*vlpm*vlwm
     &    *xnsq + 39.D0/2.D0*ro**2*vlpm*vlwm - 3.D0*ro**2*vlpm*xnsq - 6.
     &    D0*ro**2*vlpm + 57.D0/8.D0*ro**2*vlpm**2*xnsq - 3.D0/8.D0*
     &    ro**2*vlpm**2*xnsq**2 + 39.D0/4.D0*ro**2*vlpm**2 + 57.D0/2.D0
     &    *ro**2*vdmp*xnsq - 3.D0/2.D0*ro**2*vdmp*xnsq**2 + 39.D0*ro**2
     &    *vdmp + 19.D0/8.D0*ro**2*pi**2*xnsq - 1.D0/8.D0*ro**2*pi**2*
     &    xnsq**2 + 13.D0/4.D0*ro**2*pi**2 + 117.D0/8.D0*ro**3*vlpm*
     &    vlsm*xnsq + 51.D0/4.D0*ro**3*vlpm*vlsm - 21.D0*ro**3*vlpm*
     &    vltm*xnsq - 51.D0/4.D0*ro**3*vlpm*vltm - 33.D0/4.D0*ro**3*
     &    vlpm*vlwm*xnsq )
      gcoeff8 = gcoeff8 + INVG**2 * (  - 51.D0/4.D0*ro**3*vlpm*vlwm + 9.
     &    D0/2.D0*ro**3*vlpm*xnsq + 6.D0*ro**3*vlpm - 117.D0/16.D0*
     &    ro**3*vlpm**2*xnsq + 3.D0/2.D0*ro**3*vlpm**2*xnsq**2 - 51.D0/
     &    8.D0*ro**3*vlpm**2 - 117.D0/4.D0*ro**3*vdmp*xnsq + 6.D0*ro**3
     &    *vdmp*xnsq**2 - 51.D0/2.D0*ro**3*vdmp - 39.D0/16.D0*ro**3*
     &    pi**2*xnsq + 1.D0/2.D0*ro**3*pi**2*xnsq**2 - 17.D0/8.D0*ro**3
     &    *pi**2 - 27.D0/8.D0*ro**4*vlpm*vlsm*xnsq + 27.D0/8.D0*ro**4*
     &    vlpm*vltm*xnsq + 27.D0/8.D0*ro**4*vlpm*vlwm*xnsq - 3.D0/2.D0*
     &    ro**4*vlpm*xnsq + 27.D0/16.D0*ro**4*vlpm**2*xnsq + 27.D0/4.D0
     &    *ro**4*vdmp*xnsq + 9.D0/16.D0*ro**4*pi**2*xnsq - 3.D0/4.D0*
     &    vlpm**2*omro**(-1)*xnsq**2 + 3.D0/4.D0*vlpm**2*xnsq**2 - 3.D0
     &    *vdmp*omro**(-1)*xnsq**2 + 3.D0*vdmp*xnsq**2 - 1.D0/4.D0*
     &    omro**(-1)*pi**2*xnsq**2 + 1.D0/4.D0*pi**2*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**3 * (  - 3.D0/16.D0*b*t1**(-1)*ro**3*
     &    vlpm**2 - 3.D0/4.D0*b*t1**(-1)*ro**3*vlsm*vlwm*xnsq**2 + 3.D0/
     &    16.D0*b*t1**(-1)*ro**3*vlsm**2*xnsq**2 - 3.D0/4.D0*b*t1**(-1)
     &    *ro**3*vltm**2*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**3*vltm**2 + 3.D
     &    0/4.D0*b*t1**(-1)*ro**3*vlwm**2*xnsq + 3.D0/4.D0*b*t1**(-1)*
     &    ro**3*vdw*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**3*vdw*xnsq**2 - 3.D0
     &    /4.D0*b*t1**(-1)*ro**3*vdt*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**3*
     &    vdt + 1.D0/16.D0*b*t1**(-1)*ro**3*pi**2*xnsq**2 + 1.D0/16.D0*
     &    b*t1**(-1)*ro**3*pi**2 + 3.D0/16.D0*b*t1**(-1)*ro**4*vlpm**2
     &     - 3.D0/16.D0*b*t1**(-1)*ro**4*vlsm*vltm*xnsq**2 + 3.D0/16.D0
     &    *b*t1**(-1)*ro**4*vlsm*vlwm*xnsq**2 + 15.D0/16.D0*b*t1**(-1)*
     &    ro**4*vltm**2*xnsq + 3.D0/4.D0*b*t1**(-1)*ro**4*vltm**2 - 15.D
     &    0/16.D0*b*t1**(-1)*ro**4*vlwm**2*xnsq - 15.D0/16.D0*b*
     &    t1**(-1)*ro**4*vdw*xnsq + 3.D0/16.D0*b*t1**(-1)*ro**4*vdw*
     &    xnsq**2 + 15.D0/16.D0*b*t1**(-1)*ro**4*vdt*xnsq - 3.D0/16.D0*
     &    b*t1**(-1)*ro**4*vdt*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**3 * ( 3.D0/4.D0*b*t1**(-1)*ro**4*vdt - 
     &    1.D0/16.D0*b*t1**(-1)*ro**4*pi**2 - 3.D0/16.D0*b*t1**(-1)*
     &    ro**5*vltm**2*xnsq + 3.D0/16.D0*b*t1**(-1)*ro**5*vlwm**2*xnsq
     &     + 3.D0/16.D0*b*t1**(-1)*ro**5*vdw*xnsq - 3.D0/16.D0*b*
     &    t1**(-1)*ro**5*vdt*xnsq - 3.D0/16.D0*b*t2**(-1)*ro**3*vlpm**2
     &     - 3.D0/4.D0*b*t2**(-1)*ro**3*vlsm*vlwm*xnsq**2 + 3.D0/16.D0*
     &    b*t2**(-1)*ro**3*vlsm**2*xnsq**2 - 3.D0/4.D0*b*t2**(-1)*ro**3
     &    *vlwm**2 - 3.D0/4.D0*b*t2**(-1)*ro**3*vdw*xnsq**2 - 3.D0/4.D0
     &    *b*t2**(-1)*ro**3*vdw + 1.D0/16.D0*b*t2**(-1)*ro**3*pi**2*
     &    xnsq**2 + 1.D0/16.D0*b*t2**(-1)*ro**3*pi**2 + 3.D0/16.D0*b*
     &    t2**(-1)*ro**4*vlpm**2 + 3.D0/4.D0*b*t2**(-1)*ro**4*vlwm**2
     &     + 3.D0/4.D0*b*t2**(-1)*ro**4*vdw - 1.D0/16.D0*b*t2**(-1)*
     &    ro**4*pi**2 + 3.D0/4.D0*b*ro**2*vlpm**2*xnsq + 3.D0/2.D0*b*
     &    ro**2*vlpm**2 + 3.D0*b*ro**2*vlsm*vlwm*xnsq**2 - 3.D0/4.D0*b*
     &    ro**2*vlsm**2*xnsq**2 + 3.D0*b*ro**2*vltm**2*xnsq + 3.D0*b*
     &    ro**2*vltm**2 )
      gcoeff8 = gcoeff8 + INVG**3 * ( 3.D0*b*ro**2*vlwm**2 + 3.D0*b*
     &    ro**2*vdw*xnsq**2 + 3.D0*b*ro**2*vdw + 3.D0*b*ro**2*vdt*xnsq
     &     + 3.D0*b*ro**2*vdt - 1.D0/4.D0*b*ro**2*pi**2*xnsq - 1.D0/4.D0
     &    *b*ro**2*pi**2*xnsq**2 - 1.D0/2.D0*b*ro**2*pi**2 - 21.D0/16.D0
     &    *b*ro**3*vlpm**2*xnsq - 15.D0/8.D0*b*ro**3*vlpm**2 + 3.D0/4.D0
     &    *b*ro**3*vlsm*vltm*xnsq**2 - 3.D0/16.D0*b*ro**3*vlsm**2*
     &    xnsq**2 - 9.D0/2.D0*b*ro**3*vltm**2*xnsq - 15.D0/4.D0*b*ro**3
     &    *vltm**2 - 3.D0/4.D0*b*ro**3*vlwm**2*xnsq - 15.D0/4.D0*b*
     &    ro**3*vlwm**2 - 3.D0/4.D0*b*ro**3*vdw*xnsq - 15.D0/4.D0*b*
     &    ro**3*vdw - 9.D0/2.D0*b*ro**3*vdt*xnsq + 3.D0/4.D0*b*ro**3*
     &    vdt*xnsq**2 - 15.D0/4.D0*b*ro**3*vdt + 7.D0/16.D0*b*ro**3*
     &    pi**2*xnsq - 1.D0/16.D0*b*ro**3*pi**2*xnsq**2 + 5.D0/8.D0*b*
     &    ro**3*pi**2 + 9.D0/16.D0*b*ro**4*vlpm**2*xnsq + 3.D0/8.D0*b*
     &    ro**4*vlpm**2 - 3.D0/8.D0*b*ro**4*vltm*vlwm*xnsq + 21.D0/16.D0
     &    *b*ro**4*vltm**2*xnsq + 3.D0/4.D0*b*ro**4*vltm**2 + 9.D0/16.D0
     &    *b*ro**4*vlwm**2*xnsq )
      gcoeff8 = gcoeff8 + INVG**3 * ( 3.D0/4.D0*b*ro**4*vlwm**2 + 3.D0/
     &    8.D0*b*ro**4*vdw*xnsq + 3.D0/4.D0*b*ro**4*vdw + 9.D0/8.D0*b*
     &    ro**4*vdt*xnsq + 3.D0/4.D0*b*ro**4*vdt - 1.D0/8.D0*b*ro**4*
     &    pi**2*xnsq - 1.D0/8.D0*b*ro**4*pi**2 - 3.D0/8.D0*t1**(-1)*
     &    ro**3*vlpm*vlsm + 3.D0/4.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq + 3.
     &    D0/4.D0*t1**(-1)*ro**3*vlpm*vltm - 3.D0/4.D0*t1**(-1)*ro**3*
     &    vlpm*vlwm*xnsq + 3.D0/16.D0*t1**(-1)*ro**3*vlpm**2*xnsq**2 + 
     &    3.D0/16.D0*t1**(-1)*ro**3*vlpm**2 + 3.D0/4.D0*t1**(-1)*ro**3*
     &    vdmp*xnsq**2 + 3.D0/4.D0*t1**(-1)*ro**3*vdmp + 1.D0/16.D0*
     &    t1**(-1)*ro**3*pi**2*xnsq**2 + 1.D0/16.D0*t1**(-1)*ro**3*
     &    pi**2 + 9.D0/16.D0*t1**(-1)*ro**4*vlpm*vlsm - 21.D0/16.D0*
     &    t1**(-1)*ro**4*vlpm*vltm*xnsq - 9.D0/8.D0*t1**(-1)*ro**4*vlpm
     &    *vltm + 21.D0/16.D0*t1**(-1)*ro**4*vlpm*vlwm*xnsq - 3.D0/32.D0
     &    *t1**(-1)*ro**4*vlpm**2*xnsq**2 - 9.D0/32.D0*t1**(-1)*ro**4*
     &    vlpm**2 - 3.D0/8.D0*t1**(-1)*ro**4*vdmp*xnsq**2 - 9.D0/8.D0*
     &    t1**(-1)*ro**4*vdmp )
      gcoeff8 = gcoeff8 + INVG**3 * (  - 1.D0/32.D0*t1**(-1)*ro**4*
     &    pi**2*xnsq**2 - 3.D0/32.D0*t1**(-1)*ro**4*pi**2 - 3.D0/16.D0*
     &    t1**(-1)*ro**5*vlpm*vlsm + 9.D0/16.D0*t1**(-1)*ro**5*vlpm*
     &    vltm*xnsq + 3.D0/8.D0*t1**(-1)*ro**5*vlpm*vltm - 9.D0/16.D0*
     &    t1**(-1)*ro**5*vlpm*vlwm*xnsq + 3.D0/32.D0*t1**(-1)*ro**5*
     &    vlpm**2 + 3.D0/8.D0*t1**(-1)*ro**5*vdmp + 1.D0/32.D0*t1**(-1)
     &    *ro**5*pi**2 - 3.D0/8.D0*t2**(-1)*ro**3*vlpm*vlsm + 3.D0/4.D0
     &    *t2**(-1)*ro**3*vlpm*vlwm + 3.D0/16.D0*t2**(-1)*ro**3*vlpm**2
     &    *xnsq**2 + 3.D0/16.D0*t2**(-1)*ro**3*vlpm**2 + 3.D0/4.D0*
     &    t2**(-1)*ro**3*vdmp*xnsq**2 + 3.D0/4.D0*t2**(-1)*ro**3*vdmp
     &     + 1.D0/16.D0*t2**(-1)*ro**3*pi**2*xnsq**2 + 1.D0/16.D0*
     &    t2**(-1)*ro**3*pi**2 + 9.D0/16.D0*t2**(-1)*ro**4*vlpm*vlsm - 
     &    9.D0/8.D0*t2**(-1)*ro**4*vlpm*vlwm - 3.D0/32.D0*t2**(-1)*
     &    ro**4*vlpm**2*xnsq**2 - 9.D0/32.D0*t2**(-1)*ro**4*vlpm**2 - 3.
     &    D0/8.D0*t2**(-1)*ro**4*vdmp*xnsq**2 - 9.D0/8.D0*t2**(-1)*
     &    ro**4*vdmp )
      gcoeff8 = gcoeff8 + INVG**3 * (  - 1.D0/32.D0*t2**(-1)*ro**4*
     &    pi**2*xnsq**2 - 3.D0/32.D0*t2**(-1)*ro**4*pi**2 - 3.D0/16.D0*
     &    t2**(-1)*ro**5*vlpm*vlsm + 3.D0/8.D0*t2**(-1)*ro**5*vlpm*vlwm
     &     + 3.D0/32.D0*t2**(-1)*ro**5*vlpm**2 + 3.D0/8.D0*t2**(-1)*
     &    ro**5*vdmp + 1.D0/32.D0*t2**(-1)*ro**5*pi**2 + 3.D0/2.D0*
     &    ro**2*vlpm*vlsm*xnsq + 3.D0*ro**2*vlpm*vlsm - 3.D0*ro**2*vlpm
     &    *vltm*xnsq - 3.D0*ro**2*vlpm*vltm - 3.D0*ro**2*vlpm*vlwm - 3.D
     &    0/4.D0*ro**2*vlpm**2*xnsq - 3.D0/4.D0*ro**2*vlpm**2*xnsq**2
     &     - 3.D0/2.D0*ro**2*vlpm**2 - 3.D0*ro**2*vdmp*xnsq - 3.D0*
     &    ro**2*vdmp*xnsq**2 - 6.D0*ro**2*vdmp - 1.D0/4.D0*ro**2*pi**2*
     &    xnsq - 1.D0/4.D0*ro**2*pi**2*xnsq**2 - 1.D0/2.D0*ro**2*pi**2
     &     - 27.D0/8.D0*ro**3*vlpm*vlsm*xnsq - 21.D0/4.D0*ro**3*vlpm*
     &    vlsm + 6.D0*ro**3*vlpm*vltm*xnsq + 21.D0/4.D0*ro**3*vlpm*vltm
     &     + 3.D0/4.D0*ro**3*vlpm*vlwm*xnsq + 21.D0/4.D0*ro**3*vlpm*
     &    vlwm + 27.D0/16.D0*ro**3*vlpm**2*xnsq + 3.D0/16.D0*ro**3*
     &    vlpm**2*xnsq**2 )
      gcoeff8 = gcoeff8 + INVG**3 * ( 21.D0/8.D0*ro**3*vlpm**2 + 27.D0/
     &    4.D0*ro**3*vdmp*xnsq + 3.D0/4.D0*ro**3*vdmp*xnsq**2 + 21.D0/2.
     &    D0*ro**3*vdmp + 9.D0/16.D0*ro**3*pi**2*xnsq + 1.D0/16.D0*
     &    ro**3*pi**2*xnsq**2 + 7.D0/8.D0*ro**3*pi**2 + 9.D0/4.D0*ro**4
     &    *vlpm*vlsm*xnsq + 9.D0/4.D0*ro**4*vlpm*vlsm - 27.D0/8.D0*
     &    ro**4*vlpm*vltm*xnsq - 9.D0/4.D0*ro**4*vlpm*vltm - 9.D0/8.D0*
     &    ro**4*vlpm*vlwm*xnsq - 9.D0/4.D0*ro**4*vlpm*vlwm - 9.D0/8.D0*
     &    ro**4*vlpm**2*xnsq + 3.D0/16.D0*ro**4*vlpm**2*xnsq**2 - 9.D0/
     &    8.D0*ro**4*vlpm**2 - 9.D0/2.D0*ro**4*vdmp*xnsq + 3.D0/4.D0*
     &    ro**4*vdmp*xnsq**2 - 9.D0/2.D0*ro**4*vdmp - 3.D0/8.D0*ro**4*
     &    pi**2*xnsq + 1.D0/16.D0*ro**4*pi**2*xnsq**2 - 3.D0/8.D0*ro**4
     &    *pi**2 - 3.D0/8.D0*ro**5*vlpm*vlsm*xnsq + 3.D0/8.D0*ro**5*
     &    vlpm*vltm*xnsq + 3.D0/8.D0*ro**5*vlpm*vlwm*xnsq + 3.D0/16.D0*
     &    ro**5*vlpm**2*xnsq + 3.D0/4.D0*ro**5*vdmp*xnsq + 1.D0/16.D0*
     &    ro**5*pi**2*xnsq )
      gcoeff8 = gcoeff8 + 6.D0*b*t1**(-1)*ro*vlsm*vltm*xnsq**2 - 6.D0*b
     &    *t1**(-1)*ro*vlsm*vlwm*xnsq**2 + 6.D0*b*t1**(-1)*ro*vlsm*
     &    xnsq**2 + 36.D0*b*t1**(-1)*ro*vltm*xnsq - 24.D0*b*t1**(-1)*ro
     &    *vltm*xnsq**2 - 24.D0*b*t1**(-1)*ro*vltm - 6.D0*b*t1**(-1)*ro
     &    *vltm**2*xnsq - 36.D0*b*t1**(-1)*ro*vlwm*xnsq + 12.D0*b*
     &    t1**(-1)*ro*vlwm*xnsq**2 + 6.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 
     &    6.D0*b*t1**(-1)*ro*vdw*xnsq - 6.D0*b*t1**(-1)*ro*vdw*xnsq**2
     &     - 6.D0*b*t1**(-1)*ro*vdt*xnsq + 6.D0*b*t1**(-1)*ro*vdt*
     &    xnsq**2 + 6.D0*b*t1**(-1)*ro - 9.D0*b*t1**(-1)*vlsm*
     &    omro**(-2)*xnsq**2 + 6.D0*b*t1**(-1)*vlsm*omro**(-1)*xnsq**2
     &     + 3.D0*b*t1**(-1)*vlsm*xnsq**2 - 48.D0*b*t1**(-1)*vltm*xnsq
     &     + 24.D0*b*t1**(-1)*vltm*xnsq**2 + 24.D0*b*t1**(-1)*vltm + 48.
     &    D0*b*t1**(-1)*vlwm*xnsq - 48.D0*b*t1**(-1)*vlwm*xnsq**2 + 6.D0
     &    *b*t1**(-1)*omro**(-1)*xnsq**2 - 6.D0*b*t1**(-1)*xnsq**2 + 6.D
     &    0*b*t2**(-1)*ro*vlsm*xnsq**2 - 12.D0*b*t2**(-1)*ro*vlwm*
     &    xnsq**2
      gcoeff8 = gcoeff8 - 24.D0*b*t2**(-1)*ro*vlwm + 6.D0*b*t2**(-1)*ro
     &     - 9.D0*b*t2**(-1)*vlsm*omro**(-2)*xnsq**2 + 6.D0*b*t2**(-1)*
     &    vlsm*omro**(-1)*xnsq**2 + 3.D0*b*t2**(-1)*vlsm*xnsq**2 - 24.D0
     &    *b*t2**(-1)*vlwm*xnsq**2 + 24.D0*b*t2**(-1)*vlwm + 6.D0*b*
     &    t2**(-1)*omro**(-1)*xnsq**2 - 6.D0*b*t2**(-1)*xnsq**2 - 12.D0
     &    *b*ro*vlsm*xnsq**2 + 36.D0*b*ro*vltm*vlwm*xnsq - 24.D0*b*ro*
     &    vltm*xnsq + 12.D0*b*ro*vltm*xnsq**2 - 18.D0*b*ro*vltm**2*xnsq
     &     - 24.D0*b*ro*vlwm*xnsq + 12.D0*b*ro*vlwm*xnsq**2 - 18.D0*b*
     &    ro*vlwm**2*xnsq - 18.D0*b*ro*pi**2*xnsq + 12.D0*b*ro*xnsq + 
     &    18.D0*b*vlsm*omro**(-2)*xnsq**2 - 12.D0*b*vlsm*omro**(-1)*
     &    xnsq**2 - 6.D0*b*vlsm*xnsq**2 + 24.D0*b*vltm*xnsq - 24.D0*b*
     &    vltm*xnsq**2 + 24.D0*b*vlwm*xnsq - 24.D0*b*vlwm*xnsq**2 - 12.D
     &    0*b*omro**(-1)*xnsq**2 + 12.D0*b*xnsq**2 - 12.D0*t1**(-1)*ro*
     &    vlpm*vlsm + 6.D0*t1**(-1)*ro*vlpm*vltm*xnsq + 24.D0*t1**(-1)*
     &    ro*vlpm*vltm - 6.D0*t1**(-1)*ro*vlpm*vlwm*xnsq + 6.D0*
     &    t1**(-1)*ro*vlpm
      gcoeff8 = gcoeff8 + 15.D0/4.D0*t1**(-1)*ro*vlpm**2*xnsq**2 + 6.D0
     &    *t1**(-1)*ro*vlpm**2 + 15.D0*t1**(-1)*ro*vdmp*xnsq**2 + 24.D0
     &    *t1**(-1)*ro*vdmp + 5.D0/4.D0*t1**(-1)*ro*pi**2*xnsq**2 + 2.D0
     &    *t1**(-1)*ro*pi**2 + 6.D0*t1**(-1)*ro**2*vlpm*vlsm + 3.D0*
     &    t1**(-1)*ro**2*vlpm*vltm*xnsq - 12.D0*t1**(-1)*ro**2*vlpm*
     &    vltm - 3.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq - 9.D0*t1**(-1)*
     &    ro**2*vlpm - 3.D0*t1**(-1)*ro**2*vlpm**2 - 12.D0*t1**(-1)*
     &    ro**2*vdmp - t1**(-1)*ro**2*pi**2 + 9.D0/4.D0*t1**(-1)*
     &    vlpm**2*omro**(-2)*xnsq**2 - 9.D0/4.D0*t1**(-1)*vlpm**2*
     &    omro**(-1)*xnsq**2 + 9.D0*t1**(-1)*vdmp*omro**(-2)*xnsq**2 - 
     &    9.D0*t1**(-1)*vdmp*omro**(-1)*xnsq**2 + 3.D0/4.D0*t1**(-1)*
     &    omro**(-2)*pi**2*xnsq**2 - 3.D0/4.D0*t1**(-1)*omro**(-1)*
     &    pi**2*xnsq**2 - 12.D0*t2**(-1)*ro*vlpm*vlsm + 24.D0*t2**(-1)*
     &    ro*vlpm*vlwm + 6.D0*t2**(-1)*ro*vlpm + 15.D0/4.D0*t2**(-1)*ro
     &    *vlpm**2*xnsq**2 + 6.D0*t2**(-1)*ro*vlpm**2 + 15.D0*t2**(-1)*
     &    ro*vdmp*xnsq**2
      gcoeff8 = gcoeff8 + 24.D0*t2**(-1)*ro*vdmp + 5.D0/4.D0*t2**(-1)*
     &    ro*pi**2*xnsq**2 + 2.D0*t2**(-1)*ro*pi**2 + 6.D0*t2**(-1)*
     &    ro**2*vlpm*vlsm - 12.D0*t2**(-1)*ro**2*vlpm*vlwm - 9.D0*
     &    t2**(-1)*ro**2*vlpm - 3.D0*t2**(-1)*ro**2*vlpm**2 - 12.D0*
     &    t2**(-1)*ro**2*vdmp - t2**(-1)*ro**2*pi**2 + 9.D0/4.D0*
     &    t2**(-1)*vlpm**2*omro**(-2)*xnsq**2 - 9.D0/4.D0*t2**(-1)*
     &    vlpm**2*omro**(-1)*xnsq**2 + 9.D0*t2**(-1)*vdmp*omro**(-2)*
     &    xnsq**2 - 9.D0*t2**(-1)*vdmp*omro**(-1)*xnsq**2 + 3.D0/4.D0*
     &    t2**(-1)*omro**(-2)*pi**2*xnsq**2 - 3.D0/4.D0*t2**(-1)*
     &    omro**(-1)*pi**2*xnsq**2 - 24.D0*ro*vlpm*vlsm*xnsq + 24.D0*ro
     &    *vlpm*vltm*xnsq + 24.D0*ro*vlpm*vlwm*xnsq + 12.D0*ro*vlpm*
     &    xnsq + 12.D0*ro*vlpm**2*xnsq - 15.D0/2.D0*ro*vlpm**2*xnsq**2
     &     + 48.D0*ro*vdmp*xnsq - 30.D0*ro*vdmp*xnsq**2 + 4.D0*ro*pi**2
     &    *xnsq - 5.D0/2.D0*ro*pi**2*xnsq**2 + 12.D0*ro**2*vlpm*vlsm*
     &    xnsq - 12.D0*ro**2*vlpm*vltm*xnsq - 12.D0*ro**2*vlpm*vlwm*
     &    xnsq
      gcoeff8 = gcoeff8 - 18.D0*ro**2*vlpm*xnsq - 6.D0*ro**2*vlpm**2*
     &    xnsq - 24.D0*ro**2*vdmp*xnsq - 2.D0*ro**2*pi**2*xnsq - 9.D0/2.
     &    D0*vlpm**2*omro**(-2)*xnsq**2 + 9.D0/2.D0*vlpm**2*omro**(-1)*
     &    xnsq**2 - 18.D0*vdmp*omro**(-2)*xnsq**2 + 18.D0*vdmp*
     &    omro**(-1)*xnsq**2 - 3.D0/2.D0*omro**(-2)*pi**2*xnsq**2 + 3.D0
     &    /2.D0*omro**(-1)*pi**2*xnsq**2

      return
      end
