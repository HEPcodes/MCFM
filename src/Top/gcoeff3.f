      double precision function gcoeff3()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff3 =  + UBAR**(-2) * (  - 12.D0*b*t1*vlwm*xnsq + 12.D0*b*t1*
     &    vlwm*xnsq**2 - 12.D0*b*vlwm*xnsq + 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + UBAR**(-1)*INVG * ( 12.D0*b*t1**(-1)*vlwm*
     &    xnsq - 12.D0*b*t1**(-1)*vlwm - 24.D0*b*t1**(-1)*xnsq + 12.D0*
     &    b*t1**(-1)*xnsq**2 + 12.D0*b*t1**(-1) - 36.D0*b*vlwm*xnsq + 
     &    12.D0*b*vlwm*xnsq**2 + 24.D0*b*vlwm + 24.D0*b*xnsq - 12.D0*b*
     &    xnsq**2 - 12.D0*b )
      gcoeff3 = gcoeff3 + UBAR**(-1) * ( 12.D0*b*t1**(-1)*vlwm*xnsq - 
     &    12.D0*b*t1**(-1)*vlwm - 24.D0*b*t1**(-1)*xnsq + 12.D0*b*
     &    t1**(-1)*xnsq**2 + 12.D0*b*t1**(-1) - 12.D0*b*t2**(-1)*vlwm*
     &    xnsq**2 + 12.D0*b*t2**(-1)*vlwm + 24.D0*b*vlwm*xnsq + 24.D0*b
     &    *vlwm*xnsq**2 - 12.D0*b*xnsq + 12.D0*b*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG * (  - 9.D0*b*t1**(-1)*ro*vlpm**2*xnsq
     &     - 9.D0/2.D0*b*t1**(-1)*ro*vlpm**2 + 48.D0*b*t1**(-1)*ro*vlsm
     &    *vltm*xnsq**2 - 6.D0*b*t1**(-1)*ro*vlsm*vlwm*xnsq**2 - 27.D0/
     &    2.D0*b*t1**(-1)*ro*vlsm*xnsq**2 - 21.D0/2.D0*b*t1**(-1)*ro*
     &    vlsm**2*xnsq**2 - 12.D0*b*t1**(-1)*ro*vltm*vlwm*xnsq + 36.D0*
     &    b*t1**(-1)*ro*vltm*xnsq + 12.D0*b*t1**(-1)*ro*vltm*xnsq**2 + 
     &    18.D0*b*t1**(-1)*ro*vltm - 60.D0*b*t1**(-1)*ro*vltm**2*xnsq
     &     - 18.D0*b*t1**(-1)*ro*vltm**2 + 42.D0*b*t1**(-1)*ro*vlwm*
     &    xnsq - 6.D0*b*t1**(-1)*ro*vlwm*xnsq**2 - 36.D0*b*t1**(-1)*ro*
     &    vlwm + 12.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 6.D0*b*t1**(-1)*ro*
     &    vdw*xnsq - 6.D0*b*t1**(-1)*ro*vdw*xnsq**2 - 66.D0*b*t1**(-1)*
     &    ro*vdt*xnsq + 48.D0*b*t1**(-1)*ro*vdt*xnsq**2 - 18.D0*b*
     &    t1**(-1)*ro*vdt + 5.D0*b*t1**(-1)*ro*pi**2*xnsq - 7.D0/2.D0*b
     &    *t1**(-1)*ro*pi**2*xnsq**2 + 3.D0/2.D0*b*t1**(-1)*ro*pi**2 - 
     &    18.D0*b*t1**(-1)*ro*xnsq + 6.D0*b*t1**(-1)*ro*xnsq**2 + 9.D0/
     &    2.D0*b*t1**(-1)*ro**2*vlpm**2*xnsq )
      gcoeff3 = gcoeff3 + INVG * (  - 9.D0*b*t1**(-1)*ro**2*vltm*xnsq
     &     + 39.D0/2.D0*b*t1**(-1)*ro**2*vltm**2*xnsq - 9.D0*b*t1**(-1)
     &    *ro**2*vlwm*xnsq + 9.D0/2.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq + 
     &    9.D0/2.D0*b*t1**(-1)*ro**2*vdw*xnsq + 39.D0/2.D0*b*t1**(-1)*
     &    ro**2*vdt*xnsq - 1.D0/2.D0*b*t1**(-1)*ro**2*pi**2*xnsq + 3.D0
     &    *b*t1**(-1)*ro**2*xnsq + 6.D0*b*t1**(-1)*vlpm**2*xnsq + 3.D0*
     &    b*t1**(-1)*vlpm**2 - 24.D0*b*t1**(-1)*vlsm*vltm*xnsq**2 + 12.D
     &    0*b*t1**(-1)*vlsm*vlwm*xnsq**2 + 9.D0/2.D0*b*t1**(-1)*vlsm*
     &    omro**(-1)*xnsq**2 - 9.D0/2.D0*b*t1**(-1)*vlsm*xnsq**2 + 3.D0
     &    *b*t1**(-1)*vlsm**2*xnsq**2 - 24.D0*b*t1**(-1)*vltm*xnsq - 24.
     &    D0*b*t1**(-1)*vltm + 36.D0*b*t1**(-1)*vltm**2*xnsq + 12.D0*b*
     &    t1**(-1)*vltm**2 - 60.D0*b*t1**(-1)*vlwm*xnsq + 24.D0*b*
     &    t1**(-1)*vlwm*xnsq**2 + 36.D0*b*t1**(-1)*vlwm - 12.D0*b*
     &    t1**(-1)*vlwm**2*xnsq - 12.D0*b*t1**(-1)*vdw*xnsq + 12.D0*b*
     &    t1**(-1)*vdw*xnsq**2 + 36.D0*b*t1**(-1)*vdt*xnsq - 24.D0*b*
     &    t1**(-1)*vdt*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG * ( 12.D0*b*t1**(-1)*vdt - 2.D0*b*
     &    t1**(-1)*pi**2*xnsq + b*t1**(-1)*pi**2*xnsq**2 - b*t1**(-1)*
     &    pi**2 + 24.D0*b*t1**(-1)*xnsq - 12.D0*b*t1**(-1)*xnsq**2 - 12.
     &    D0*b*t1**(-1) + 6.D0*b*t2**(-2)*ro**2*vlwm**2 + 6.D0*b*
     &    t2**(-2)*ro**2*vdw + b*t2**(-2)*ro**2*pi**2 - 9.D0/2.D0*b*
     &    t2**(-1)*ro*vlpm**2 - 6.D0*b*t2**(-1)*ro*vlsm*vlwm*xnsq**2 + 
     &    3.D0/2.D0*b*t2**(-1)*ro*vlsm**2*xnsq**2 - 6.D0*b*t2**(-1)*ro*
     &    vlwm - 30.D0*b*t2**(-1)*ro*vlwm**2 - 6.D0*b*t2**(-1)*ro*vdw*
     &    xnsq**2 - 30.D0*b*t2**(-1)*ro*vdw + 1.D0/2.D0*b*t2**(-1)*ro*
     &    pi**2*xnsq**2 - 1.D0/2.D0*b*t2**(-1)*ro*pi**2 - 6.D0*b*
     &    t2**(-1)*ro + 6.D0*b*t2**(-1)*ro**2*vlwm**2*xnsq + 6.D0*b*
     &    t2**(-1)*ro**2*vdw*xnsq + b*t2**(-1)*ro**2*pi**2*xnsq + 3.D0*
     &    b*t2**(-1)*vlpm**2 + 12.D0*b*t2**(-1)*vlsm*vlwm*xnsq**2 + 6.D0
     &    *b*t2**(-1)*vlsm*omro**(-1)*xnsq**2 - 6.D0*b*t2**(-1)*vlsm*
     &    xnsq**2 - 3.D0*b*t2**(-1)*vlsm**2*xnsq**2 + 12.D0*b*t2**(-1)*
     &    vlwm*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG * (  - 12.D0*b*t2**(-1)*vlwm + 12.D0*b*
     &    t2**(-1)*vlwm**2 + 12.D0*b*t2**(-1)*vdw*xnsq**2 + 12.D0*b*
     &    t2**(-1)*vdw - b*t2**(-1)*pi**2*xnsq**2 - b*t2**(-1)*pi**2 - 
     &    9.D0*b*ro*vlpm**2*xnsq - 3.D0*b*ro*vlpm**2 + 18.D0*b*ro*vlsm*
     &    vltm*xnsq**2 + 6.D0*b*ro*vlsm*vlwm*xnsq**2 - 12.D0*b*ro*vlsm*
     &    xnsq**2 - 6.D0*b*ro*vlsm**2*xnsq**2 + 18.D0*b*ro*vltm*xnsq + 
     &    24.D0*b*ro*vltm*xnsq**2 + 12.D0*b*ro*vltm - 18.D0*b*ro*
     &    vltm**2*xnsq + 42.D0*b*ro*vlwm*xnsq - 12.D0*b*ro*vlwm - 30.D0
     &    *b*ro*vlwm**2*xnsq - 30.D0*b*ro*vdw*xnsq + 6.D0*b*ro*vdw*
     &    xnsq**2 - 18.D0*b*ro*vdt*xnsq + 18.D0*b*ro*vdt*xnsq**2 + b*ro
     &    *pi**2*xnsq - 2.D0*b*ro*pi**2*xnsq**2 + 3.D0*b*ro*pi**2 - 18.D
     &    0*b*ro*xnsq + 3.D0*b*vlpm**2*xnsq - 84.D0*b*vlsm*vltm*xnsq**2
     &     - 21.D0*b*vlsm*omro**(-1)*xnsq**2 + 33.D0*b*vlsm*xnsq**2 + 
     &    21.D0*b*vlsm**2*xnsq**2 - 96.D0*b*vltm*xnsq - 24.D0*b*vltm*
     &    xnsq**2 - 12.D0*b*vltm**2 - 48.D0*b*vlwm*xnsq + 72.D0*b*vlwm
     &     + 12.D0*b*vlwm**2*xnsq )
      gcoeff3 = gcoeff3 + INVG * ( 12.D0*b*vlwm**2 + 12.D0*b*vdw*xnsq
     &     + 12.D0*b*vdw - 84.D0*b*vdt*xnsq**2 - 12.D0*b*vdt - b*pi**2*
     &    xnsq + 7.D0*b*pi**2*xnsq**2 + 96.D0*b*xnsq - 36.D0*b*xnsq**2
     &     - 12.D0*b - 24.D0*t1**(-1)*ro*vlpm*vlsm*xnsq - 12.D0*
     &    t1**(-1)*ro*vlpm*vlsm + 72.D0*t1**(-1)*ro*vlpm*vltm*xnsq + 24.
     &    D0*t1**(-1)*ro*vlpm*vltm - 24.D0*t1**(-1)*ro*vlpm*vlwm*xnsq
     &     - 45.D0*t1**(-1)*ro*vlpm*xnsq - 36.D0*t1**(-1)*ro*vlpm + 12.D
     &    0*t1**(-1)*ro*vlpm**2*xnsq - 87.D0/8.D0*t1**(-1)*ro*vlpm**2*
     &    xnsq**2 + 6.D0*t1**(-1)*ro*vlpm**2 + 48.D0*t1**(-1)*ro*vdmp*
     &    xnsq - 87.D0/2.D0*t1**(-1)*ro*vdmp*xnsq**2 + 24.D0*t1**(-1)*
     &    ro*vdmp + 4.D0*t1**(-1)*ro*pi**2*xnsq - 29.D0/8.D0*t1**(-1)*
     &    ro*pi**2*xnsq**2 + 2.D0*t1**(-1)*ro*pi**2 + 15.D0*t1**(-1)*
     &    ro**2*vlpm*vlsm*xnsq + 15.D0/2.D0*t1**(-1)*ro**2*vlpm*vlsm - 
     &    42.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 15.D0*t1**(-1)*ro**2*
     &    vlpm*vltm + 12.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 57.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm*xnsq )
      gcoeff3 = gcoeff3 + INVG * ( 15.D0*t1**(-1)*ro**2*vlpm - 15.D0/2.D
     &    0*t1**(-1)*ro**2*vlpm**2*xnsq + 27.D0/8.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq**2 - 15.D0/4.D0*t1**(-1)*ro**2*vlpm**2 - 30.D0*
     &    t1**(-1)*ro**2*vdmp*xnsq + 27.D0/2.D0*t1**(-1)*ro**2*vdmp*
     &    xnsq**2 - 15.D0*t1**(-1)*ro**2*vdmp - 5.D0/2.D0*t1**(-1)*
     &    ro**2*pi**2*xnsq + 9.D0/8.D0*t1**(-1)*ro**2*pi**2*xnsq**2 - 5.
     &    D0/4.D0*t1**(-1)*ro**2*pi**2 - 3.D0*t1**(-1)*ro**3*vlpm*vlsm*
     &    xnsq + 3.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq + 3.D0*t1**(-1)*
     &    ro**3*vlpm*vlwm*xnsq - 15.D0/2.D0*t1**(-1)*ro**3*vlpm*xnsq + 
     &    3.D0/2.D0*t1**(-1)*ro**3*vlpm**2*xnsq + 6.D0*t1**(-1)*ro**3*
     &    vdmp*xnsq + 1.D0/2.D0*t1**(-1)*ro**3*pi**2*xnsq + 12.D0*
     &    t1**(-1)*vlpm*vlsm*xnsq + 6.D0*t1**(-1)*vlpm*vlsm - 36.D0*
     &    t1**(-1)*vlpm*vltm*xnsq - 12.D0*t1**(-1)*vlpm*vltm + 12.D0*
     &    t1**(-1)*vlpm*vlwm*xnsq + 24.D0*t1**(-1)*vlpm*xnsq + 24.D0*
     &    t1**(-1)*vlpm - 9.D0/8.D0*t1**(-1)*vlpm**2*omro**(-1)*xnsq**2
     &     - 6.D0*t1**(-1)*vlpm**2*xnsq )
      gcoeff3 = gcoeff3 + INVG * ( 33.D0/8.D0*t1**(-1)*vlpm**2*xnsq**2
     &     - 3.D0*t1**(-1)*vlpm**2 - 9.D0/2.D0*t1**(-1)*vdmp*omro**(-1)
     &    *xnsq**2 - 24.D0*t1**(-1)*vdmp*xnsq + 33.D0/2.D0*t1**(-1)*
     &    vdmp*xnsq**2 - 12.D0*t1**(-1)*vdmp - 3.D0/8.D0*t1**(-1)*
     &    omro**(-1)*pi**2*xnsq**2 - 2.D0*t1**(-1)*pi**2*xnsq + 11.D0/8.
     &    D0*t1**(-1)*pi**2*xnsq**2 - t1**(-1)*pi**2 - 6.D0*t2**(-1)*ro
     &    *vlpm*vlsm + 12.D0*t2**(-1)*ro*vlpm*vlwm + 3.D0/2.D0*t2**(-1)
     &    *ro*vlpm**2*xnsq**2 + 3.D0*t2**(-1)*ro*vlpm**2 + 6.D0*
     &    t2**(-1)*ro*vdmp*xnsq**2 + 12.D0*t2**(-1)*ro*vdmp + 1.D0/2.D0
     &    *t2**(-1)*ro*pi**2*xnsq**2 + t2**(-1)*ro*pi**2 + 3.D0/2.D0*
     &    t2**(-1)*ro**2*vlpm*vlsm - 3.D0*t2**(-1)*ro**2*vlpm*vlwm + 3.D
     &    0*t2**(-1)*ro**2*vlpm - 3.D0/4.D0*t2**(-1)*ro**2*vlpm**2 - 3.D
     &    0*t2**(-1)*ro**2*vdmp - 1.D0/4.D0*t2**(-1)*ro**2*pi**2 + 6.D0
     &    *t2**(-1)*vlpm*vlsm - 12.D0*t2**(-1)*vlpm*vlwm - 3.D0/2.D0*
     &    t2**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.D0/2.D0*t2**(-1)*
     &    vlpm**2*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG * (  - 3.D0*t2**(-1)*vlpm**2 - 6.D0*
     &    t2**(-1)*vdmp*omro**(-1)*xnsq**2 - 6.D0*t2**(-1)*vdmp*xnsq**2
     &     - 12.D0*t2**(-1)*vdmp - 1.D0/2.D0*t2**(-1)*omro**(-1)*pi**2*
     &    xnsq**2 - 1.D0/2.D0*t2**(-1)*pi**2*xnsq**2 - t2**(-1)*pi**2
     &     - 15.D0*ro*vlpm*vlsm*xnsq - 12.D0*ro*vlpm*vlsm + 18.D0*ro*
     &    vlpm*vltm*xnsq - 6.D0*ro*vlpm*vltm + 12.D0*ro*vlpm*vlwm*xnsq
     &     + 30.D0*ro*vlpm*vlwm - 45.D0*ro*vlpm*xnsq - 24.D0*ro*vlpm + 
     &    15.D0/2.D0*ro*vlpm**2*xnsq - 63.D0/4.D0*ro*vlpm**2*xnsq**2 + 
     &    6.D0*ro*vlpm**2 + 30.D0*ro*vdmp*xnsq - 63.D0*ro*vdmp*xnsq**2
     &     + 24.D0*ro*vdmp + 5.D0/2.D0*ro*pi**2*xnsq - 21.D0/4.D0*ro*
     &    pi**2*xnsq**2 + 2.D0*ro*pi**2 + 6.D0*ro**2*vlpm*vlsm*xnsq - 9.
     &    D0*ro**2*vlpm*vltm*xnsq + 6.D0*ro**2*vlpm*vltm - 3.D0*ro**2*
     &    vlpm*vlwm*xnsq - 6.D0*ro**2*vlpm*vlwm + 21.D0*ro**2*vlpm*xnsq
     &     - 3.D0*ro**2*vlpm**2*xnsq - 12.D0*ro**2*vdmp*xnsq - ro**2*
     &    pi**2*xnsq + 6.D0*vlpm*vlsm*xnsq + 12.D0*vlpm*vltm - 12.D0*
     &    vlpm*vlwm*xnsq )
      gcoeff3 = gcoeff3 + INVG * (  - 12.D0*vlpm*vlwm + 24.D0*vlpm*xnsq
     &     + 12.D0*vlpm + 21.D0/4.D0*vlpm**2*omro**(-1)*xnsq**2 - 3.D0*
     &    vlpm**2*xnsq + 63.D0/4.D0*vlpm**2*xnsq**2 + 21.D0*vdmp*
     &    omro**(-1)*xnsq**2 - 12.D0*vdmp*xnsq + 63.D0*vdmp*xnsq**2 + 7.
     &    D0/4.D0*omro**(-1)*pi**2*xnsq**2 - pi**2*xnsq + 21.D0/4.D0*
     &    pi**2*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**2 * ( 15.D0*b*t1**(-1)*ro*vlpm**2*xnsq
     &     + 51.D0/4.D0*b*t1**(-1)*ro*vlpm**2 - 12.D0*b*t1**(-1)*ro*
     &    vlsm*vltm*xnsq**2 + 3.D0*b*t1**(-1)*ro*vlsm*vlwm*xnsq**2 + 9.D
     &    0/4.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2 - 6.D0*b*t1**(-1)*ro*
     &    vltm*xnsq - 6.D0*b*t1**(-1)*ro*vltm + 63.D0*b*t1**(-1)*ro*
     &    vltm**2*xnsq + 51.D0*b*t1**(-1)*ro*vltm**2 - 12.D0*b*t1**(-1)
     &    *ro*vlwm*xnsq + 6.D0*b*t1**(-1)*ro*vlwm*xnsq**2 + 6.D0*b*
     &    t1**(-1)*ro*vlwm - 3.D0*b*t1**(-1)*ro*vlwm**2*xnsq - 3.D0*b*
     &    t1**(-1)*ro*vdw*xnsq + 3.D0*b*t1**(-1)*ro*vdw*xnsq**2 + 63.D0
     &    *b*t1**(-1)*ro*vdt*xnsq - 12.D0*b*t1**(-1)*ro*vdt*xnsq**2 + 
     &    51.D0*b*t1**(-1)*ro*vdt - 5.D0*b*t1**(-1)*ro*pi**2*xnsq + 3.D0
     &    /4.D0*b*t1**(-1)*ro*pi**2*xnsq**2 - 17.D0/4.D0*b*t1**(-1)*ro*
     &    pi**2 - 99.D0/8.D0*b*t1**(-1)*ro**2*vlpm**2*xnsq - 63.D0/8.D0
     &    *b*t1**(-1)*ro**2*vlpm**2 + 51.D0/4.D0*b*t1**(-1)*ro**2*vlsm*
     &    vltm*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2
     &     - 3.D0/2.D0*b*t1**(-1)*ro**2*vlsm*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**2 * (  - 27.D0/8.D0*b*t1**(-1)*ro**2*
     &    vlsm**2*xnsq**2 + 9.D0/2.D0*b*t1**(-1)*ro**2*vltm*xnsq + 3.D0/
     &    2.D0*b*t1**(-1)*ro**2*vltm*xnsq**2 + 6.D0*b*t1**(-1)*ro**2*
     &    vltm - 207.D0/4.D0*b*t1**(-1)*ro**2*vltm**2*xnsq - 63.D0/2.D0
     &    *b*t1**(-1)*ro**2*vltm**2 + 21.D0/2.D0*b*t1**(-1)*ro**2*vlwm*
     &    xnsq - 3.D0/2.D0*b*t1**(-1)*ro**2*vlwm*xnsq**2 - 6.D0*b*
     &    t1**(-1)*ro**2*vlwm + 9.D0/4.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq
     &     + 9.D0/4.D0*b*t1**(-1)*ro**2*vdw*xnsq + 3.D0/4.D0*b*t1**(-1)
     &    *ro**2*vdw*xnsq**2 - 207.D0/4.D0*b*t1**(-1)*ro**2*vdt*xnsq + 
     &    51.D0/4.D0*b*t1**(-1)*ro**2*vdt*xnsq**2 - 63.D0/2.D0*b*
     &    t1**(-1)*ro**2*vdt + 33.D0/8.D0*b*t1**(-1)*ro**2*pi**2*xnsq
     &     - 9.D0/8.D0*b*t1**(-1)*ro**2*pi**2*xnsq**2 + 21.D0/8.D0*b*
     &    t1**(-1)*ro**2*pi**2 + 57.D0/16.D0*b*t1**(-1)*ro**3*vlpm**2*
     &    xnsq + 15.D0/16.D0*b*t1**(-1)*ro**3*vlpm**2 - 9.D0/4.D0*b*
     &    t1**(-1)*ro**3*vltm*vlwm*xnsq - 3.D0/2.D0*b*t1**(-1)*ro**3*
     &    vltm*xnsq )
      gcoeff3 = gcoeff3 + INVG**2 * ( 45.D0/4.D0*b*t1**(-1)*ro**3*
     &    vltm**2*xnsq + 15.D0/4.D0*b*t1**(-1)*ro**3*vltm**2 - 3.D0/2.D0
     &    *b*t1**(-1)*ro**3*vlwm*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**3*
     &    vlwm**2*xnsq + 3.D0/8.D0*b*t1**(-1)*ro**3*vdw*xnsq + 81.D0/8.D
     &    0*b*t1**(-1)*ro**3*vdt*xnsq + 15.D0/4.D0*b*t1**(-1)*ro**3*vdt
     &     - 11.D0/16.D0*b*t1**(-1)*ro**3*pi**2*xnsq - 5.D0/16.D0*b*
     &    t1**(-1)*ro**3*pi**2 - 6.D0*b*t1**(-1)*vlpm**2*xnsq - 6.D0*b*
     &    t1**(-1)*vlpm**2 - 24.D0*b*t1**(-1)*vltm**2*xnsq - 24.D0*b*
     &    t1**(-1)*vltm**2 - 24.D0*b*t1**(-1)*vdt*xnsq - 24.D0*b*
     &    t1**(-1)*vdt + 2.D0*b*t1**(-1)*pi**2*xnsq + 2.D0*b*t1**(-1)*
     &    pi**2 + 3.D0/4.D0*b*t2**(-2)*ro**3*vlwm**2 + 3.D0/4.D0*b*
     &    t2**(-2)*ro**3*vdw + 1.D0/8.D0*b*t2**(-2)*ro**3*pi**2 + 3.D0/
     &    4.D0*b*t2**(-1)*ro*vlpm**2 + 3.D0*b*t2**(-1)*ro*vlsm*vlwm*
     &    xnsq**2 - 3.D0/4.D0*b*t2**(-1)*ro*vlsm**2*xnsq**2 + 3.D0*b*
     &    t2**(-1)*ro*vlwm**2 + 3.D0*b*t2**(-1)*ro*vdw*xnsq**2 + 3.D0*b
     &    *t2**(-1)*ro*vdw )
      gcoeff3 = gcoeff3 + INVG**2 * (  - 1.D0/4.D0*b*t2**(-1)*ro*pi**2*
     &    xnsq**2 - 1.D0/4.D0*b*t2**(-1)*ro*pi**2 - 3.D0/8.D0*b*
     &    t2**(-1)*ro**2*vlpm**2 + 3.D0/2.D0*b*t2**(-1)*ro**2*vlsm*vlwm
     &    *xnsq**2 - 3.D0/2.D0*b*t2**(-1)*ro**2*vlsm*xnsq**2 - 3.D0/8.D0
     &    *b*t2**(-1)*ro**2*vlsm**2*xnsq**2 + 3.D0*b*t2**(-1)*ro**2*
     &    vlwm*xnsq**2 - 9.D0/2.D0*b*t2**(-1)*ro**2*vlwm**2 + 3.D0/2.D0
     &    *b*t2**(-1)*ro**2*vdw*xnsq**2 - 9.D0/2.D0*b*t2**(-1)*ro**2*
     &    vdw - 1.D0/8.D0*b*t2**(-1)*ro**2*pi**2*xnsq**2 - 3.D0/8.D0*b*
     &    t2**(-1)*ro**2*pi**2 - 15.D0/16.D0*b*t2**(-1)*ro**3*vlpm**2
     &     + 3.D0/4.D0*b*t2**(-1)*ro**3*vlwm**2*xnsq - 15.D0/4.D0*b*
     &    t2**(-1)*ro**3*vlwm**2 + 3.D0/4.D0*b*t2**(-1)*ro**3*vdw*xnsq
     &     - 15.D0/4.D0*b*t2**(-1)*ro**3*vdw + 1.D0/8.D0*b*t2**(-1)*
     &    ro**3*pi**2*xnsq + 5.D0/16.D0*b*t2**(-1)*ro**3*pi**2 + 111.D0/
     &    4.D0*b*ro*vlpm**2*xnsq + 18.D0*b*ro*vlpm**2 - 45.D0*b*ro*vlsm
     &    *vltm*xnsq**2 + 6.D0*b*ro*vlsm*xnsq**2 + 45.D0/4.D0*b*ro*
     &    vlsm**2*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**2 * ( 12.D0*b*ro*vltm*vlwm*xnsq - 24.D0
     &    *b*ro*vltm*xnsq - 6.D0*b*ro*vltm*xnsq**2 - 30.D0*b*ro*vltm + 
     &    114.D0*b*ro*vltm**2*xnsq + 57.D0*b*ro*vltm**2 - 54.D0*b*ro*
     &    vlwm*xnsq + 12.D0*b*ro*vlwm*xnsq**2 + 30.D0*b*ro*vlwm + 9.D0*
     &    b*ro*vlwm**2*xnsq + 15.D0*b*ro*vlwm**2 + 15.D0*b*ro*vdw*xnsq
     &     + 15.D0*b*ro*vdw + 120.D0*b*ro*vdt*xnsq - 45.D0*b*ro*vdt*
     &    xnsq**2 + 57.D0*b*ro*vdt - 45.D0/4.D0*b*ro*pi**2*xnsq + 15.D0/
     &    4.D0*b*ro*pi**2*xnsq**2 - 6.D0*b*ro*pi**2 - 105.D0/8.D0*b*
     &    ro**2*vlpm**2*xnsq - 9.D0/2.D0*b*ro**2*vlpm**2 + 15.D0/2.D0*b
     &    *ro**2*vlsm*vltm*xnsq**2 - 3.D0/2.D0*b*ro**2*vlsm*vlwm*
     &    xnsq**2 - 3.D0/2.D0*b*ro**2*vlsm**2*xnsq**2 + 9.D0/2.D0*b*
     &    ro**2*vltm*vlwm*xnsq + 6.D0*b*ro**2*vltm*xnsq + 3.D0*b*ro**2*
     &    vltm*xnsq**2 + 6.D0*b*ro**2*vltm - 195.D0/4.D0*b*ro**2*
     &    vltm**2*xnsq - 39.D0/2.D0*b*ro**2*vltm**2 + 12.D0*b*ro**2*
     &    vlwm*xnsq - 3.D0*b*ro**2*vlwm*xnsq**2 - 6.D0*b*ro**2*vlwm - 9.
     &    D0/4.D0*b*ro**2*vlwm**2*xnsq )
      gcoeff3 = gcoeff3 + INVG**2 * ( 9.D0/2.D0*b*ro**2*vlwm**2 - 3.D0/
     &    2.D0*b*ro**2*vdw*xnsq**2 + 9.D0/2.D0*b*ro**2*vdw - 93.D0/2.D0
     &    *b*ro**2*vdt*xnsq + 15.D0/2.D0*b*ro**2*vdt*xnsq**2 - 39.D0/2.D
     &    0*b*ro**2*vdt + 25.D0/8.D0*b*ro**2*pi**2*xnsq - 1.D0/2.D0*b*
     &    ro**2*pi**2*xnsq**2 + 2.D0*b*ro**2*pi**2 + 15.D0/4.D0*b*ro**3
     &    *vltm**2*xnsq - 15.D0/4.D0*b*ro**3*vlwm**2*xnsq - 15.D0/4.D0*
     &    b*ro**3*vdw*xnsq + 15.D0/4.D0*b*ro**3*vdt*xnsq - 15.D0*b*
     &    vlpm**2*xnsq - 12.D0*b*vlpm**2 + 24.D0*b*vlsm*vltm*xnsq**2 - 
     &    12.D0*b*vlsm*vlwm*xnsq**2 - 3.D0*b*vlsm**2*xnsq**2 + 24.D0*b*
     &    vltm*xnsq + 24.D0*b*vltm - 60.D0*b*vltm**2*xnsq - 36.D0*b*
     &    vltm**2 + 48.D0*b*vlwm*xnsq - 24.D0*b*vlwm*xnsq**2 - 24.D0*b*
     &    vlwm - 12.D0*b*vlwm**2 - 12.D0*b*vdw*xnsq**2 - 12.D0*b*vdw - 
     &    60.D0*b*vdt*xnsq + 24.D0*b*vdt*xnsq**2 - 36.D0*b*vdt + 5.D0*b
     &    *pi**2*xnsq - b*pi**2*xnsq**2 + 4.D0*b*pi**2 + 36.D0*t1**(-1)
     &    *ro*vlpm*vlsm*xnsq + 63.D0/2.D0*t1**(-1)*ro*vlpm*vlsm - 75.D0
     &    *t1**(-1)*ro*vlpm*vltm*xnsq )
      gcoeff3 = gcoeff3 + INVG**2 * (  - 63.D0*t1**(-1)*ro*vlpm*vltm + 
     &    3.D0*t1**(-1)*ro*vlpm*vlwm*xnsq + 6.D0*t1**(-1)*ro*vlpm*xnsq
     &     + 6.D0*t1**(-1)*ro*vlpm - 18.D0*t1**(-1)*ro*vlpm**2*xnsq + 9.
     &    D0/4.D0*t1**(-1)*ro*vlpm**2*xnsq**2 - 63.D0/4.D0*t1**(-1)*ro*
     &    vlpm**2 - 72.D0*t1**(-1)*ro*vdmp*xnsq + 9.D0*t1**(-1)*ro*vdmp
     &    *xnsq**2 - 63.D0*t1**(-1)*ro*vdmp - 6.D0*t1**(-1)*ro*pi**2*
     &    xnsq + 3.D0/4.D0*t1**(-1)*ro*pi**2*xnsq**2 - 21.D0/4.D0*
     &    t1**(-1)*ro*pi**2 - 153.D0/4.D0*t1**(-1)*ro**2*vlpm*vlsm*xnsq
     &     - 27.D0*t1**(-1)*ro**2*vlpm*vlsm + 321.D0/4.D0*t1**(-1)*
     &    ro**2*vlpm*vltm*xnsq + 54.D0*t1**(-1)*ro**2*vlpm*vltm - 15.D0/
     &    4.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq - 12.D0*t1**(-1)*ro**2*
     &    vlpm*xnsq - 21.D0/2.D0*t1**(-1)*ro**2*vlpm + 153.D0/8.D0*
     &    t1**(-1)*ro**2*vlpm**2*xnsq - 9.D0/2.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq**2 + 27.D0/2.D0*t1**(-1)*ro**2*vlpm**2 + 153.D0/
     &    2.D0*t1**(-1)*ro**2*vdmp*xnsq - 18.D0*t1**(-1)*ro**2*vdmp*
     &    xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**2 * ( 54.D0*t1**(-1)*ro**2*vdmp + 51.D0/
     &    8.D0*t1**(-1)*ro**2*pi**2*xnsq - 3.D0/2.D0*t1**(-1)*ro**2*
     &    pi**2*xnsq**2 + 9.D0/2.D0*t1**(-1)*ro**2*pi**2 + 135.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm*vlsm*xnsq + 15.D0/2.D0*t1**(-1)*ro**3*
     &    vlpm*vlsm - 255.D0/8.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq - 15.D0
     &    *t1**(-1)*ro**3*vlpm*vltm - 15.D0/8.D0*t1**(-1)*ro**3*vlpm*
     &    vlwm*xnsq + 15.D0/2.D0*t1**(-1)*ro**3*vlpm*xnsq + 9.D0/2.D0*
     &    t1**(-1)*ro**3*vlpm - 135.D0/16.D0*t1**(-1)*ro**3*vlpm**2*
     &    xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm**2*xnsq**2 - 15.D0/4.D0*
     &    t1**(-1)*ro**3*vlpm**2 - 135.D0/4.D0*t1**(-1)*ro**3*vdmp*xnsq
     &     + 6.D0*t1**(-1)*ro**3*vdmp*xnsq**2 - 15.D0*t1**(-1)*ro**3*
     &    vdmp - 45.D0/16.D0*t1**(-1)*ro**3*pi**2*xnsq + 1.D0/2.D0*
     &    t1**(-1)*ro**3*pi**2*xnsq**2 - 5.D0/4.D0*t1**(-1)*ro**3*pi**2
     &     - 21.D0/8.D0*t1**(-1)*ro**4*vlpm*vlsm*xnsq + 21.D0/8.D0*
     &    t1**(-1)*ro**4*vlpm*vltm*xnsq + 21.D0/8.D0*t1**(-1)*ro**4*
     &    vlpm*vlwm*xnsq )
      gcoeff3 = gcoeff3 + INVG**2 * (  - 3.D0/2.D0*t1**(-1)*ro**4*vlpm*
     &    xnsq + 21.D0/16.D0*t1**(-1)*ro**4*vlpm**2*xnsq + 21.D0/4.D0*
     &    t1**(-1)*ro**4*vdmp*xnsq + 7.D0/16.D0*t1**(-1)*ro**4*pi**2*
     &    xnsq - 12.D0*t1**(-1)*vlpm*vlsm*xnsq - 12.D0*t1**(-1)*vlpm*
     &    vlsm + 24.D0*t1**(-1)*vlpm*vltm*xnsq + 24.D0*t1**(-1)*vlpm*
     &    vltm + 6.D0*t1**(-1)*vlpm**2*xnsq + 6.D0*t1**(-1)*vlpm**2 + 
     &    24.D0*t1**(-1)*vdmp*xnsq + 24.D0*t1**(-1)*vdmp + 2.D0*
     &    t1**(-1)*pi**2*xnsq + 2.D0*t1**(-1)*pi**2 + 3.D0/2.D0*
     &    t2**(-1)*ro*vlpm*vlsm - 3.D0*t2**(-1)*ro*vlpm*vlwm - 3.D0/4.D0
     &    *t2**(-1)*ro*vlpm**2*xnsq**2 - 3.D0/4.D0*t2**(-1)*ro*vlpm**2
     &     - 3.D0*t2**(-1)*ro*vdmp*xnsq**2 - 3.D0*t2**(-1)*ro*vdmp - 1.D
     &    0/4.D0*t2**(-1)*ro*pi**2*xnsq**2 - 1.D0/4.D0*t2**(-1)*ro*
     &    pi**2 - 3.D0/2.D0*t2**(-1)*ro**2*vlpm*vlsm + 3.D0*t2**(-1)*
     &    ro**2*vlpm*vlwm - 3.D0/2.D0*t2**(-1)*ro**2*vlpm + 3.D0/4.D0*
     &    t2**(-1)*ro**2*vlpm**2 + 3.D0*t2**(-1)*ro**2*vdmp + 1.D0/4.D0
     &    *t2**(-1)*ro**2*pi**2 )
      gcoeff3 = gcoeff3 + INVG**2 * ( 3.D0/2.D0*t2**(-1)*ro**3*vlpm + 
     &    141.D0/2.D0*ro*vlpm*vlsm*xnsq + 48.D0*ro*vlpm*vlsm - 138.D0*
     &    ro*vlpm*vltm*xnsq - 75.D0*ro*vlpm*vltm - 3.D0*ro*vlpm*vlwm*
     &    xnsq - 21.D0*ro*vlpm*vlwm + 54.D0*ro*vlpm*xnsq + 48.D0*ro*
     &    vlpm - 141.D0/4.D0*ro*vlpm**2*xnsq + 51.D0/4.D0*ro*vlpm**2*
     &    xnsq**2 - 24.D0*ro*vlpm**2 - 141.D0*ro*vdmp*xnsq + 51.D0*ro*
     &    vdmp*xnsq**2 - 96.D0*ro*vdmp - 47.D0/4.D0*ro*pi**2*xnsq + 17.D
     &    0/4.D0*ro*pi**2*xnsq**2 - 8.D0*ro*pi**2 - 201.D0/4.D0*ro**2*
     &    vlpm*vlsm*xnsq - 24.D0*ro**2*vlpm*vlsm + 195.D0/2.D0*ro**2*
     &    vlpm*vltm*xnsq + 87.D0/2.D0*ro**2*vlpm*vltm + 3.D0*ro**2*vlpm
     &    *vlwm*xnsq + 9.D0/2.D0*ro**2*vlpm*vlwm - 39.D0*ro**2*vlpm*
     &    xnsq - 24.D0*ro**2*vlpm + 201.D0/8.D0*ro**2*vlpm**2*xnsq - 27.
     &    D0/4.D0*ro**2*vlpm**2*xnsq**2 + 12.D0*ro**2*vlpm**2 + 201.D0/
     &    2.D0*ro**2*vdmp*xnsq - 27.D0*ro**2*vdmp*xnsq**2 + 48.D0*ro**2
     &    *vdmp + 67.D0/8.D0*ro**2*pi**2*xnsq - 9.D0/4.D0*ro**2*pi**2*
     &    xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**2 * ( 4.D0*ro**2*pi**2 + 39.D0/4.D0*
     &    ro**3*vlpm*vlsm*xnsq - 39.D0/2.D0*ro**3*vlpm*vltm*xnsq - 9.D0/
     &    2.D0*ro**3*vlpm*vltm + 9.D0/2.D0*ro**3*vlpm*vlwm + 9.D0*ro**3
     &    *vlpm*xnsq - 39.D0/8.D0*ro**3*vlpm**2*xnsq - 39.D0/2.D0*ro**3
     &    *vdmp*xnsq - 13.D0/8.D0*ro**3*pi**2*xnsq - 30.D0*vlpm*vlsm*
     &    xnsq - 24.D0*vlpm*vlsm + 60.D0*vlpm*vltm*xnsq + 36.D0*vlpm*
     &    vltm + 12.D0*vlpm*vlwm - 24.D0*vlpm*xnsq - 24.D0*vlpm + 15.D0
     &    *vlpm**2*xnsq - 3.D0*vlpm**2*xnsq**2 + 12.D0*vlpm**2 + 60.D0*
     &    vdmp*xnsq - 12.D0*vdmp*xnsq**2 + 48.D0*vdmp + 5.D0*pi**2*xnsq
     &     - pi**2*xnsq**2 + 4.D0*pi**2 )
      gcoeff3 = gcoeff3 + INVG**3 * (  - 3.D0/2.D0*b*t1**(-1)*ro*
     &    vlpm**2*xnsq - 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2 - 6.D0*b*
     &    t1**(-1)*ro*vltm**2*xnsq - 6.D0*b*t1**(-1)*ro*vltm**2 - 6.D0*
     &    b*t1**(-1)*ro*vdt*xnsq - 6.D0*b*t1**(-1)*ro*vdt + 1.D0/2.D0*b
     &    *t1**(-1)*ro*pi**2*xnsq + 1.D0/2.D0*b*t1**(-1)*ro*pi**2 + 27.D
     &    0/8.D0*b*t1**(-1)*ro**2*vlpm**2*xnsq + 3.D0*b*t1**(-1)*ro**2*
     &    vlpm**2 - 3.D0/2.D0*b*t1**(-1)*ro**2*vlsm*vltm*xnsq**2 + 3.D0/
     &    8.D0*b*t1**(-1)*ro**2*vlsm**2*xnsq**2 + 27.D0/2.D0*b*t1**(-1)
     &    *ro**2*vltm**2*xnsq + 12.D0*b*t1**(-1)*ro**2*vltm**2 + 27.D0/
     &    2.D0*b*t1**(-1)*ro**2*vdt*xnsq - 3.D0/2.D0*b*t1**(-1)*ro**2*
     &    vdt*xnsq**2 + 12.D0*b*t1**(-1)*ro**2*vdt - 9.D0/8.D0*b*
     &    t1**(-1)*ro**2*pi**2*xnsq + 1.D0/8.D0*b*t1**(-1)*ro**2*pi**2*
     &    xnsq**2 - b*t1**(-1)*ro**2*pi**2 - 39.D0/16.D0*b*t1**(-1)*
     &    ro**3*vlpm**2*xnsq - 27.D0/16.D0*b*t1**(-1)*ro**3*vlpm**2 + 9.
     &    D0/8.D0*b*t1**(-1)*ro**3*vlsm*vltm*xnsq**2 + 3.D0/8.D0*b*
     &    t1**(-1)*ro**3*vlsm*vlwm*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**3 * (  - 3.D0/8.D0*b*t1**(-1)*ro**3*
     &    vlsm**2*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*ro**3*vltm*vlwm*xnsq
     &     - 9.D0*b*t1**(-1)*ro**3*vltm**2*xnsq - 27.D0/4.D0*b*t1**(-1)
     &    *ro**3*vltm**2 + 3.D0/8.D0*b*t1**(-1)*ro**3*vdw*xnsq + 3.D0/8.
     &    D0*b*t1**(-1)*ro**3*vdw*xnsq**2 - 69.D0/8.D0*b*t1**(-1)*ro**3
     &    *vdt*xnsq + 9.D0/8.D0*b*t1**(-1)*ro**3*vdt*xnsq**2 - 27.D0/4.D
     &    0*b*t1**(-1)*ro**3*vdt + 11.D0/16.D0*b*t1**(-1)*ro**3*pi**2*
     &    xnsq - 1.D0/8.D0*b*t1**(-1)*ro**3*pi**2*xnsq**2 + 9.D0/16.D0*
     &    b*t1**(-1)*ro**3*pi**2 + 9.D0/16.D0*b*t1**(-1)*ro**4*vlpm**2*
     &    xnsq + 3.D0/16.D0*b*t1**(-1)*ro**4*vlpm**2 - 3.D0/8.D0*b*
     &    t1**(-1)*ro**4*vltm*vlwm*xnsq + 27.D0/16.D0*b*t1**(-1)*ro**4*
     &    vltm**2*xnsq + 3.D0/4.D0*b*t1**(-1)*ro**4*vltm**2 + 3.D0/16.D0
     &    *b*t1**(-1)*ro**4*vlwm**2*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**4*
     &    vdt*xnsq + 3.D0/4.D0*b*t1**(-1)*ro**4*vdt - 1.D0/8.D0*b*
     &    t1**(-1)*ro**4*pi**2*xnsq - 1.D0/16.D0*b*t1**(-1)*ro**4*pi**2
     &     + 3.D0/16.D0*b*t2**(-1)*ro**3*vlpm**2 )
      gcoeff3 = gcoeff3 + INVG**3 * ( 3.D0/4.D0*b*t2**(-1)*ro**3*vlsm*
     &    vlwm*xnsq**2 - 3.D0/16.D0*b*t2**(-1)*ro**3*vlsm**2*xnsq**2 + 
     &    3.D0/4.D0*b*t2**(-1)*ro**3*vlwm**2 + 3.D0/4.D0*b*t2**(-1)*
     &    ro**3*vdw*xnsq**2 + 3.D0/4.D0*b*t2**(-1)*ro**3*vdw - 1.D0/16.D
     &    0*b*t2**(-1)*ro**3*pi**2*xnsq**2 - 1.D0/16.D0*b*t2**(-1)*
     &    ro**3*pi**2 - 3.D0/16.D0*b*t2**(-1)*ro**4*vlpm**2 - 3.D0/4.D0
     &    *b*t2**(-1)*ro**4*vlwm**2 - 3.D0/4.D0*b*t2**(-1)*ro**4*vdw + 
     &    1.D0/16.D0*b*t2**(-1)*ro**4*pi**2 - 15.D0*b*ro*vlpm**2*xnsq
     &     - 27.D0/2.D0*b*ro*vlpm**2 + 6.D0*b*ro*vlsm*vltm*xnsq**2 - 3.D
     &    0/2.D0*b*ro*vlsm**2*xnsq**2 - 60.D0*b*ro*vltm**2*xnsq - 54.D0
     &    *b*ro*vltm**2 - 60.D0*b*ro*vdt*xnsq + 6.D0*b*ro*vdt*xnsq**2
     &     - 54.D0*b*ro*vdt + 5.D0*b*ro*pi**2*xnsq - 1.D0/2.D0*b*ro*
     &    pi**2*xnsq**2 + 9.D0/2.D0*b*ro*pi**2 + 99.D0/8.D0*b*ro**2*
     &    vlpm**2*xnsq + 9.D0*b*ro**2*vlpm**2 - 6.D0*b*ro**2*vlsm*vltm*
     &    xnsq**2 - 3.D0/2.D0*b*ro**2*vlsm*vlwm*xnsq**2 + 15.D0/8.D0*b*
     &    ro**2*vlsm**2*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**3 * (  - 3.D0*b*ro**2*vltm*vlwm*xnsq + 
     &    48.D0*b*ro**2*vltm**2*xnsq + 75.D0/2.D0*b*ro**2*vltm**2 - 3.D0
     &    /2.D0*b*ro**2*vlwm**2*xnsq - 3.D0/2.D0*b*ro**2*vlwm**2 - 3.D0
     &    *b*ro**2*vdw*xnsq - 3.D0/2.D0*b*ro**2*vdw*xnsq**2 - 3.D0/2.D0
     &    *b*ro**2*vdw + 93.D0/2.D0*b*ro**2*vdt*xnsq - 6.D0*b*ro**2*vdt
     &    *xnsq**2 + 75.D0/2.D0*b*ro**2*vdt - 29.D0/8.D0*b*ro**2*pi**2*
     &    xnsq + 5.D0/8.D0*b*ro**2*pi**2*xnsq**2 - 3.D0*b*ro**2*pi**2
     &     - 27.D0/8.D0*b*ro**3*vlpm**2*xnsq - 3.D0/2.D0*b*ro**3*
     &    vlpm**2 + 3.D0/4.D0*b*ro**3*vlsm*vltm*xnsq**2 - 3.D0/4.D0*b*
     &    ro**3*vlsm*vlwm*xnsq**2 + 9.D0/4.D0*b*ro**3*vltm*vlwm*xnsq - 
     &    105.D0/8.D0*b*ro**3*vltm**2*xnsq - 15.D0/2.D0*b*ro**3*vltm**2
     &     + 15.D0/8.D0*b*ro**3*vlwm**2*xnsq + 3.D0/2.D0*b*ro**3*
     &    vlwm**2 + 3.D0*b*ro**3*vdw*xnsq - 3.D0/4.D0*b*ro**3*vdw*
     &    xnsq**2 + 3.D0/2.D0*b*ro**3*vdw - 12.D0*b*ro**3*vdt*xnsq + 3.D
     &    0/4.D0*b*ro**3*vdt*xnsq**2 - 15.D0/2.D0*b*ro**3*vdt + 3.D0/4.D
     &    0*b*ro**3*pi**2*xnsq )
      gcoeff3 = gcoeff3 + INVG**3 * ( 1.D0/2.D0*b*ro**3*pi**2 + 3.D0/4.D
     &    0*b*ro**4*vltm**2*xnsq - 3.D0/4.D0*b*ro**4*vlwm**2*xnsq - 3.D0
     &    /4.D0*b*ro**4*vdw*xnsq + 3.D0/4.D0*b*ro**4*vdt*xnsq + 6.D0*b*
     &    vlpm**2*xnsq + 6.D0*b*vlpm**2 + 24.D0*b*vltm**2*xnsq + 24.D0*
     &    b*vltm**2 + 24.D0*b*vdt*xnsq + 24.D0*b*vdt - 2.D0*b*pi**2*
     &    xnsq - 2.D0*b*pi**2 - 3.D0*t1**(-1)*ro*vlpm*vlsm*xnsq - 3.D0*
     &    t1**(-1)*ro*vlpm*vlsm + 6.D0*t1**(-1)*ro*vlpm*vltm*xnsq + 6.D0
     &    *t1**(-1)*ro*vlpm*vltm + 3.D0/2.D0*t1**(-1)*ro*vlpm**2*xnsq
     &     + 3.D0/2.D0*t1**(-1)*ro*vlpm**2 + 6.D0*t1**(-1)*ro*vdmp*xnsq
     &     + 6.D0*t1**(-1)*ro*vdmp + 1.D0/2.D0*t1**(-1)*ro*pi**2*xnsq
     &     + 1.D0/2.D0*t1**(-1)*ro*pi**2 + 33.D0/4.D0*t1**(-1)*ro**2*
     &    vlpm*vlsm*xnsq + 15.D0/2.D0*t1**(-1)*ro**2*vlpm*vlsm - 33.D0/
     &    2.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 15.D0*t1**(-1)*ro**2*
     &    vlpm*vltm - 33.D0/8.D0*t1**(-1)*ro**2*vlpm**2*xnsq + 3.D0/8.D0
     &    *t1**(-1)*ro**2*vlpm**2*xnsq**2 - 15.D0/4.D0*t1**(-1)*ro**2*
     &    vlpm**2 )
      gcoeff3 = gcoeff3 + INVG**3 * (  - 33.D0/2.D0*t1**(-1)*ro**2*vdmp
     &    *xnsq + 3.D0/2.D0*t1**(-1)*ro**2*vdmp*xnsq**2 - 15.D0*
     &    t1**(-1)*ro**2*vdmp - 11.D0/8.D0*t1**(-1)*ro**2*pi**2*xnsq + 
     &    1.D0/8.D0*t1**(-1)*ro**2*pi**2*xnsq**2 - 5.D0/4.D0*t1**(-1)*
     &    ro**2*pi**2 - 63.D0/8.D0*t1**(-1)*ro**3*vlpm*vlsm*xnsq - 6.D0
     &    *t1**(-1)*ro**3*vlpm*vlsm + 123.D0/8.D0*t1**(-1)*ro**3*vlpm*
     &    vltm*xnsq + 12.D0*t1**(-1)*ro**3*vlpm*vltm + 3.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm*vlwm*xnsq + 63.D0/16.D0*t1**(-1)*ro**3*
     &    vlpm**2*xnsq - 9.D0/16.D0*t1**(-1)*ro**3*vlpm**2*xnsq**2 + 3.D
     &    0*t1**(-1)*ro**3*vlpm**2 + 63.D0/4.D0*t1**(-1)*ro**3*vdmp*
     &    xnsq - 9.D0/4.D0*t1**(-1)*ro**3*vdmp*xnsq**2 + 12.D0*t1**(-1)
     &    *ro**3*vdmp + 21.D0/16.D0*t1**(-1)*ro**3*pi**2*xnsq - 3.D0/16.
     &    D0*t1**(-1)*ro**3*pi**2*xnsq**2 + t1**(-1)*ro**3*pi**2 + 3.D0
     &    *t1**(-1)*ro**4*vlpm*vlsm*xnsq + 3.D0/2.D0*t1**(-1)*ro**4*
     &    vlpm*vlsm - 21.D0/4.D0*t1**(-1)*ro**4*vlpm*vltm*xnsq - 3.D0*
     &    t1**(-1)*ro**4*vlpm*vltm )
      gcoeff3 = gcoeff3 + INVG**3 * (  - 3.D0/4.D0*t1**(-1)*ro**4*vlpm*
     &    vlwm*xnsq - 3.D0/2.D0*t1**(-1)*ro**4*vlpm**2*xnsq + 3.D0/16.D0
     &    *t1**(-1)*ro**4*vlpm**2*xnsq**2 - 3.D0/4.D0*t1**(-1)*ro**4*
     &    vlpm**2 - 6.D0*t1**(-1)*ro**4*vdmp*xnsq + 3.D0/4.D0*t1**(-1)*
     &    ro**4*vdmp*xnsq**2 - 3.D0*t1**(-1)*ro**4*vdmp - 1.D0/2.D0*
     &    t1**(-1)*ro**4*pi**2*xnsq + 1.D0/16.D0*t1**(-1)*ro**4*pi**2*
     &    xnsq**2 - 1.D0/4.D0*t1**(-1)*ro**4*pi**2 - 3.D0/8.D0*t1**(-1)
     &    *ro**5*vlpm*vlsm*xnsq + 3.D0/8.D0*t1**(-1)*ro**5*vlpm*vltm*
     &    xnsq + 3.D0/8.D0*t1**(-1)*ro**5*vlpm*vlwm*xnsq + 3.D0/16.D0*
     &    t1**(-1)*ro**5*vlpm**2*xnsq + 3.D0/4.D0*t1**(-1)*ro**5*vdmp*
     &    xnsq + 1.D0/16.D0*t1**(-1)*ro**5*pi**2*xnsq - 36.D0*ro*vlpm*
     &    vlsm*xnsq - 33.D0*ro*vlpm*vlsm + 72.D0*ro*vlpm*vltm*xnsq + 66.
     &    D0*ro*vlpm*vltm + 18.D0*ro*vlpm**2*xnsq - 3.D0/2.D0*ro*
     &    vlpm**2*xnsq**2 + 33.D0/2.D0*ro*vlpm**2 + 72.D0*ro*vdmp*xnsq
     &     - 6.D0*ro*vdmp*xnsq**2 + 66.D0*ro*vdmp + 6.D0*ro*pi**2*xnsq
     &     - 1.D0/2.D0*ro*pi**2*xnsq**2 )
      gcoeff3 = gcoeff3 + INVG**3 * ( 11.D0/2.D0*ro*pi**2 + 153.D0/4.D0
     &    *ro**2*vlpm*vlsm*xnsq + 30.D0*ro**2*vlpm*vlsm - 153.D0/2.D0*
     &    ro**2*vlpm*vltm*xnsq - 123.D0/2.D0*ro**2*vlpm*vltm + 3.D0/2.D0
     &    *ro**2*vlpm*vlwm - 153.D0/8.D0*ro**2*vlpm**2*xnsq + 21.D0/8.D0
     &    *ro**2*vlpm**2*xnsq**2 - 15.D0*ro**2*vlpm**2 - 153.D0/2.D0*
     &    ro**2*vdmp*xnsq + 21.D0/2.D0*ro**2*vdmp*xnsq**2 - 60.D0*ro**2
     &    *vdmp - 51.D0/8.D0*ro**2*pi**2*xnsq + 7.D0/8.D0*ro**2*pi**2*
     &    xnsq**2 - 5.D0*ro**2*pi**2 - 33.D0/2.D0*ro**3*vlpm*vlsm*xnsq
     &     - 9.D0*ro**3*vlpm*vlsm + 33.D0*ro**3*vlpm*vltm*xnsq + 21.D0*
     &    ro**3*vlpm*vltm - 3.D0*ro**3*vlpm*vlwm + 33.D0/4.D0*ro**3*
     &    vlpm**2*xnsq - 9.D0/8.D0*ro**3*vlpm**2*xnsq**2 + 9.D0/2.D0*
     &    ro**3*vlpm**2 + 33.D0*ro**3*vdmp*xnsq - 9.D0/2.D0*ro**3*vdmp*
     &    xnsq**2 + 18.D0*ro**3*vdmp + 11.D0/4.D0*ro**3*pi**2*xnsq - 3.D
     &    0/8.D0*ro**3*pi**2*xnsq**2 + 3.D0/2.D0*ro**3*pi**2 + 9.D0/4.D0
     &    *ro**4*vlpm*vlsm*xnsq - 9.D0/2.D0*ro**4*vlpm*vltm*xnsq - 3.D0/
     &    2.D0*ro**4*vlpm*vltm )
      gcoeff3 = gcoeff3 + INVG**3 * ( 3.D0/2.D0*ro**4*vlpm*vlwm - 9.D0/
     &    8.D0*ro**4*vlpm**2*xnsq - 9.D0/2.D0*ro**4*vdmp*xnsq - 3.D0/8.D
     &    0*ro**4*pi**2*xnsq + 12.D0*vlpm*vlsm*xnsq + 12.D0*vlpm*vlsm
     &     - 24.D0*vlpm*vltm*xnsq - 24.D0*vlpm*vltm - 6.D0*vlpm**2*xnsq
     &     - 6.D0*vlpm**2 - 24.D0*vdmp*xnsq - 24.D0*vdmp - 2.D0*pi**2*
     &    xnsq - 2.D0*pi**2 )
      gcoeff3 = gcoeff3 - 32.D0*XLF*b*t1*TR*xn*xnsq + 16.D0*XLF*b*TR*xn
     &    *xnsq - 3.D0*b*t1**(-1)*ro*vlpm**2*xnsq - 3.D0*b*t1**(-1)*ro*
     &    vlpm**2 + 12.D0*b*t1**(-1)*ro*vltm*vlwm*xnsq - 12.D0*b*
     &    t1**(-1)*ro*vltm*xnsq - 6.D0*b*t1**(-1)*ro*vltm**2*xnsq - 12.D
     &    0*b*t1**(-1)*ro*vlwm*xnsq + 6.D0*b*t1**(-1)*ro*vlwm**2*xnsq
     &     + 12.D0*b*t1**(-1)*ro*vdw*xnsq - b*t1**(-1)*ro*pi**2*xnsq + 
     &    3.D0*b*t1**(-1)*ro*pi**2 + 12.D0*b*t1**(-1)*ro*xnsq + 6.D0*b*
     &    t1**(-1)*vlpm**2*xnsq + 60.D0*b*t1**(-1)*vlsm*vltm*xnsq**2 - 
     &    12.D0*b*t1**(-1)*vlsm*vlwm*xnsq**2 + 18.D0*b*t1**(-1)*vlsm*
     &    omro**(-1)*xnsq**2 - 30.D0*b*t1**(-1)*vlsm*xnsq**2 - 12.D0*b*
     &    t1**(-1)*vlsm**2*xnsq**2 + 72.D0*b*t1**(-1)*vltm*xnsq + 24.D0
     &    *b*t1**(-1)*vltm*xnsq**2 - 24.D0*b*t1**(-1)*vltm + 12.D0*b*
     &    t1**(-1)*vltm**2*xnsq - 48.D0*b*t1**(-1)*vlwm + 12.D0*b*
     &    t1**(-1)*vlwm**2*xnsq + 12.D0*b*t1**(-1)*vdw*xnsq - 12.D0*b*
     &    t1**(-1)*vdw*xnsq**2 + 12.D0*b*t1**(-1)*vdt*xnsq + 60.D0*b*
     &    t1**(-1)*vdt*xnsq**2
      gcoeff3 = gcoeff3 - 2.D0*b*t1**(-1)*pi**2*xnsq - 4.D0*b*t1**(-1)*
     &    pi**2*xnsq**2 - 72.D0*b*t1**(-1)*xnsq + 24.D0*b*t1**(-1)*
     &    xnsq**2 - 24.D0*b*t1*ro*vlpm**2*TR*xn*xnsq + 24.D0*b*t1*ro*TR
     &    *pi**2*xn*xnsq - 192.D0*b*t1*ro*TR*xn*xnsq - 32.D0*b*t1*TR*xn
     &    *xnsq + 16.D0*b*t1*xnsq**2 + 12.D0*b*t2**(-2)*ro*vlwm**2*xnsq
     &     + 12.D0*b*t2**(-2)*ro*vdw*xnsq + 2.D0*b*t2**(-2)*ro*pi**2*
     &    xnsq + 3.D0*b*t2**(-1)*ro*vlpm**2 + 12.D0*b*t2**(-1)*ro*
     &    vlwm**2*xnsq + 12.D0*b*t2**(-1)*ro*vdw*xnsq + 2.D0*b*t2**(-1)
     &    *ro*pi**2*xnsq - 3.D0*b*t2**(-1)*ro*pi**2 + 24.D0*b*t2**(-1)*
     &    vlsm*omro**(-1)*xnsq**2 + 24.D0*b*t2**(-1)*vlsm*xnsq**2 + 48.D
     &    0*b*t2**(-1)*vlwm*xnsq - 48.D0*b*t2**(-1)*vlwm*xnsq**2 - 72.D0
     &    *b*t2**(-1)*vlwm - 48.D0*b*t2**(-1)*xnsq + 24.D0*b*t2**(-1)
     &     + 12.D0*b*ro*vlpm**2*TR*xn*xnsq - 12.D0*b*ro*TR*pi**2*xn*
     &    xnsq + 96.D0*b*ro*TR*xn*xnsq - 48.D0*b*vlsm*omro**(-1)*
     &    xnsq**2 - 48.D0*b*vlsm*xnsq**2 + 72.D0*b*vltm*vlwm*xnsq - 24.D
     &    0*b*vltm*xnsq
      gcoeff3 = gcoeff3 + 48.D0*b*vltm*xnsq**2 - 36.D0*b*vltm**2*xnsq
     &     - 24.D0*b*vlwm*xnsq + 48.D0*b*vlwm*xnsq**2 - 36.D0*b*vlwm**2
     &    *xnsq + 16.D0*b*TR*xn*xnsq - 36.D0*b*pi**2*xnsq - 8.D0*b*
     &    xnsq**2 - 18.D0*t1**(-1)*ro*vlpm*vlsm*xnsq + 6.D0*t1**(-1)*ro
     &    *vlpm*vlsm + 6.D0*t1**(-1)*ro*vlpm*vltm*xnsq - 12.D0*t1**(-1)
     &    *ro*vlpm*vltm + 30.D0*t1**(-1)*ro*vlpm*vlwm*xnsq - 6.D0*
     &    t1**(-1)*ro*vlpm*xnsq - 12.D0*t1**(-1)*ro*vlpm + 9.D0*
     &    t1**(-1)*ro*vlpm**2*xnsq + 3.D0/2.D0*t1**(-1)*ro*vlpm**2*
     &    xnsq**2 - 3.D0*t1**(-1)*ro*vlpm**2 + 36.D0*t1**(-1)*ro*vdmp*
     &    xnsq + 6.D0*t1**(-1)*ro*vdmp*xnsq**2 - 12.D0*t1**(-1)*ro*vdmp
     &     + 3.D0*t1**(-1)*ro*pi**2*xnsq + 1.D0/2.D0*t1**(-1)*ro*pi**2*
     &    xnsq**2 - t1**(-1)*ro*pi**2 + 6.D0*t1**(-1)*ro**2*vlpm*vlsm*
     &    xnsq - 6.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 6.D0*t1**(-1)*
     &    ro**2*vlpm*vlwm*xnsq - 6.D0*t1**(-1)*ro**2*vlpm*xnsq - 3.D0*
     &    t1**(-1)*ro**2*vlpm**2*xnsq - 12.D0*t1**(-1)*ro**2*vdmp*xnsq
     &     - t1**(-1)*ro**2*pi**2*xnsq
      gcoeff3 = gcoeff3 + 12.D0*t1**(-1)*vlpm*vlsm*xnsq - 12.D0*
     &    t1**(-1)*vlpm*vltm*xnsq - 12.D0*t1**(-1)*vlpm*vlwm*xnsq + 12.D
     &    0*t1**(-1)*vlpm*xnsq + 24.D0*t1**(-1)*vlpm - 9.D0/2.D0*
     &    t1**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 6.D0*t1**(-1)*vlpm**2*
     &    xnsq - 15.D0/2.D0*t1**(-1)*vlpm**2*xnsq**2 - 18.D0*t1**(-1)*
     &    vdmp*omro**(-1)*xnsq**2 - 24.D0*t1**(-1)*vdmp*xnsq - 30.D0*
     &    t1**(-1)*vdmp*xnsq**2 - 3.D0/2.D0*t1**(-1)*omro**(-1)*pi**2*
     &    xnsq**2 - 2.D0*t1**(-1)*pi**2*xnsq - 5.D0/2.D0*t1**(-1)*pi**2
     &    *xnsq**2 + 96.D0*t1*ro*vlpm*TR*xn*xnsq - 96.D0*t1*ro**2*vlpm*
     &    TR*xn*xnsq + 6.D0*t2**(-1)*ro*vlpm*vlsm - 12.D0*t2**(-1)*ro*
     &    vlpm*vlwm - 12.D0*t2**(-1)*ro*vlpm - 3.D0*t2**(-1)*ro*vlpm**2
     &     - 12.D0*t2**(-1)*ro*vdmp - t2**(-1)*ro*pi**2 + 24.D0*
     &    t2**(-1)*vlpm - 6.D0*t2**(-1)*vlpm**2*omro**(-1)*xnsq**2 + 6.D
     &    0*t2**(-1)*vlpm**2*xnsq**2 - 24.D0*t2**(-1)*vdmp*omro**(-1)*
     &    xnsq**2 + 24.D0*t2**(-1)*vdmp*xnsq**2 - 2.D0*t2**(-1)*
     &    omro**(-1)*pi**2*xnsq**2
      gcoeff3 = gcoeff3 + 2.D0*t2**(-1)*pi**2*xnsq**2 + 12.D0*ro*vlpm*
     &    vlsm*xnsq - 12.D0*ro*vlpm*vltm*xnsq - 12.D0*ro*vlpm*vlwm*xnsq
     &     - 48.D0*ro*vlpm*TR*xn*xnsq - 24.D0*ro*vlpm*xnsq - 6.D0*ro*
     &    vlpm**2*xnsq - 24.D0*ro*vdmp*xnsq - 2.D0*ro*pi**2*xnsq + 48.D0
     &    *ro**2*vlpm*TR*xn*xnsq + 48.D0*vlpm*xnsq + 12.D0*vlpm**2*
     &    omro**(-1)*xnsq**2 - 12.D0*vlpm**2*xnsq**2 + 48.D0*vdmp*
     &    omro**(-1)*xnsq**2 - 48.D0*vdmp*xnsq**2 + 4.D0*omro**(-1)*
     &    pi**2*xnsq**2 - 4.D0*pi**2*xnsq**2

      return
      end
