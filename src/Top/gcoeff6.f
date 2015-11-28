      double precision function gcoeff6()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff6 =  + TBAR**(-2) * ( 24.D0*b*t1*vltm*xnsq - 12.D0*b*t1*
     &    vltm*xnsq**2 - 12.D0*b*t1*vltm - 12.D0*b*t1**2*vltm*xnsq + 12.
     &    D0*b*t1**2*vltm*xnsq**2 )
      gcoeff6 = gcoeff6 + TBAR**(-1) * ( 72.D0*b*t1*vltm*xnsq - 120.D0*
     &    b*t1*vltm*xnsq**2 + 12.D0*b*t1*xnsq - 12.D0*b*t1*xnsq**2 - 
     &    192.D0*b*vltm*xnsq + 120.D0*b*vltm*xnsq**2 + 72.D0*b*vltm - 
     &    24.D0*b*xnsq + 12.D0*b*xnsq**2 + 12.D0*b )
      gcoeff6 = gcoeff6 + UBAR**(-2) * (  - 12.D0*b*t1*vlwm*xnsq**2 + 
     &    12.D0*b*t1*vlwm - 12.D0*b*t1**2*vlwm*xnsq + 12.D0*b*t1**2*
     &    vlwm*xnsq**2 + 12.D0*b*vlwm*xnsq - 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + UBAR**(-1)*INVG * (  - 30.D0*b*t1**(-1)*vlwm*
     &    xnsq**2 - 36.D0*b*t1**(-1)*vlwm + 30.D0*b*vlwm*xnsq**2 + 36.D0
     &    *b*vlwm )
      gcoeff6 = gcoeff6 + UBAR**(-1) * (  - 30.D0*b*t1**(-1)*vlwm*
     &    xnsq**2 - 36.D0*b*t1**(-1)*vlwm - 72.D0*b*t1*vlwm*xnsq + 120.D
     &    0*b*t1*vlwm*xnsq**2 - 12.D0*b*t1*xnsq + 12.D0*b*t1*xnsq**2 - 
     &    30.D0*b*t2**(-1)*vlwm*xnsq**2 - 36.D0*b*t2**(-1)*vlwm - 120.D0
     &    *b*vlwm*xnsq + 72.D0*b*vlwm - 12.D0*b*xnsq + 12.D0*b )
      gcoeff6 = gcoeff6 + INVG * (  - 15.D0/2.D0*b*t1**(-1)*ro*vlpm**2*
     &    xnsq + 15.D0/2.D0*b*t1**(-1)*ro*vlpm**2 - 15.D0*b*t1**(-1)*ro
     &    *vlsm*vlwm*xnsq**2 + 9.D0*b*t1**(-1)*ro*vlsm*xnsq**2 + 15.D0/
     &    4.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2 - 12.D0*b*t1**(-1)*ro*vltm
     &    *vlwm*xnsq + 24.D0*b*t1**(-1)*ro*vltm*xnsq - 15.D0*b*t1**(-1)
     &    *ro*vltm*xnsq**2 - 9.D0*b*t1**(-1)*ro*vltm - 12.D0*b*t1**(-1)
     &    *ro*vltm**2*xnsq - 6.D0*b*t1**(-1)*ro*vltm**2 + 3.D0*b*
     &    t1**(-1)*ro*vlwm*xnsq - 117.D0/4.D0*b*t1**(-1)*ro*vlwm*
     &    xnsq**2 - 30.D0*b*t1**(-1)*ro*vlwm - 30.D0*b*t1**(-1)*ro*
     &    vlwm**2*xnsq + 36.D0*b*t1**(-1)*ro*vlwm**2 - 36.D0*b*t1**(-1)
     &    *ro*vdw*xnsq - 15.D0*b*t1**(-1)*ro*vdw*xnsq**2 + 36.D0*b*
     &    t1**(-1)*ro*vdw - 18.D0*b*t1**(-1)*ro*vdt*xnsq - 6.D0*b*
     &    t1**(-1)*ro*vdt + 9.D0/2.D0*b*t1**(-1)*ro*pi**2*xnsq + 5.D0/4.
     &    D0*b*t1**(-1)*ro*pi**2*xnsq**2 - 5.D0/2.D0*b*t1**(-1)*ro*
     &    pi**2 + 3.D0*b*t1**(-1)*ro*xnsq**2 - 3.D0*b*t1**(-1)*ro + 51.D
     &    0/8.D0*b*t1**(-1)*ro**2*vlpm**2*xnsq )
      gcoeff6 = gcoeff6 + INVG * ( 9.D0/8.D0*b*t1**(-1)*ro**2*vlpm**2
     &     + 6.D0*b*t1**(-1)*ro**2*vlsm*vltm*xnsq**2 + 6.D0*b*t1**(-1)*
     &    ro**2*vlsm*vlwm*xnsq**2 - 3.D0*b*t1**(-1)*ro**2*vlsm*xnsq**2
     &     - 3.D0*b*t1**(-1)*ro**2*vlsm**2*xnsq**2 - 15.D0/2.D0*b*
     &    t1**(-1)*ro**2*vltm*vlwm*xnsq - 3.D0*b*t1**(-1)*ro**2*vltm*
     &    xnsq + 3.D0*b*t1**(-1)*ro**2*vltm*xnsq**2 + 9.D0/4.D0*b*
     &    t1**(-1)*ro**2*vltm**2*xnsq + 15.D0/2.D0*b*t1**(-1)*ro**2*
     &    vltm**2 - 3.D0*b*t1**(-1)*ro**2*vlwm*xnsq + 3.D0*b*t1**(-1)*
     &    ro**2*vlwm*xnsq**2 + 75.D0/4.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq
     &     - 3.D0*b*t1**(-1)*ro**2*vlwm**2 + 15.D0*b*t1**(-1)*ro**2*vdw
     &    *xnsq + 6.D0*b*t1**(-1)*ro**2*vdw*xnsq**2 - 3.D0*b*t1**(-1)*
     &    ro**2*vdw - 3.D0/2.D0*b*t1**(-1)*ro**2*vdt*xnsq + 6.D0*b*
     &    t1**(-1)*ro**2*vdt*xnsq**2 + 15.D0/2.D0*b*t1**(-1)*ro**2*vdt
     &     - 3.D0/8.D0*b*t1**(-1)*ro**2*pi**2*xnsq - b*t1**(-1)*ro**2*
     &    pi**2*xnsq**2 - 3.D0/8.D0*b*t1**(-1)*ro**2*pi**2 + 3.D0*b*
     &    t1**(-1)*vlpm**2*xnsq )
      gcoeff6 = gcoeff6 + INVG * (  - 6.D0*b*t1**(-1)*vlpm**2 + 30.D0*b
     &    *t1**(-1)*vlsm*vlwm*xnsq**2 - 3.D0/2.D0*b*t1**(-1)*vlsm*
     &    omro**(-1)*xnsq**2 + 3.D0/2.D0*b*t1**(-1)*vlsm*xnsq**2 - 15.D0
     &    /2.D0*b*t1**(-1)*vlsm**2*xnsq**2 + 30.D0*b*t1**(-1)*vlwm*
     &    xnsq**2 + 36.D0*b*t1**(-1)*vlwm + 12.D0*b*t1**(-1)*vlwm**2*
     &    xnsq - 24.D0*b*t1**(-1)*vlwm**2 + 12.D0*b*t1**(-1)*vdw*xnsq
     &     + 30.D0*b*t1**(-1)*vdw*xnsq**2 - 24.D0*b*t1**(-1)*vdw - b*
     &    t1**(-1)*pi**2*xnsq - 5.D0/2.D0*b*t1**(-1)*pi**2*xnsq**2 + 2.D
     &    0*b*t1**(-1)*pi**2 + 3.D0*b*t2**(-2)*ro**2*vlwm**2 + 3.D0*b*
     &    t2**(-2)*ro**2*vdw + 1.D0/2.D0*b*t2**(-2)*ro**2*pi**2 + 12.D0
     &    *b*t2**(-1)*ro*vlpm**2 - 27.D0*b*t2**(-1)*ro*vlsm*vlwm*
     &    xnsq**2 + 9.D0*b*t2**(-1)*ro*vlsm*xnsq**2 + 27.D0/4.D0*b*
     &    t2**(-1)*ro*vlsm**2*xnsq**2 - 117.D0/4.D0*b*t2**(-1)*ro*vlwm*
     &    xnsq**2 - 33.D0*b*t2**(-1)*ro*vlwm + 48.D0*b*t2**(-1)*ro*
     &    vlwm**2 - 27.D0*b*t2**(-1)*ro*vdw*xnsq**2 + 48.D0*b*t2**(-1)*
     &    ro*vdw )
      gcoeff6 = gcoeff6 + INVG * ( 9.D0/4.D0*b*t2**(-1)*ro*pi**2*
     &    xnsq**2 - 4.D0*b*t2**(-1)*ro*pi**2 + 3.D0*b*t2**(-1)*ro*
     &    xnsq**2 - 3.D0*b*t2**(-1)*ro - 21.D0/8.D0*b*t2**(-1)*ro**2*
     &    vlpm**2 + 3.D0*b*t2**(-1)*ro**2*vlwm**2*xnsq - 21.D0/2.D0*b*
     &    t2**(-1)*ro**2*vlwm**2 + 3.D0*b*t2**(-1)*ro**2*vdw*xnsq - 21.D
     &    0/2.D0*b*t2**(-1)*ro**2*vdw + 1.D0/2.D0*b*t2**(-1)*ro**2*
     &    pi**2*xnsq + 7.D0/8.D0*b*t2**(-1)*ro**2*pi**2 - 9.D0*b*
     &    t2**(-1)*vlpm**2 + 30.D0*b*t2**(-1)*vlsm*vlwm*xnsq**2 - 3.D0/
     &    2.D0*b*t2**(-1)*vlsm*omro**(-1)*xnsq**2 + 3.D0/2.D0*b*
     &    t2**(-1)*vlsm*xnsq**2 - 15.D0/2.D0*b*t2**(-1)*vlsm**2*xnsq**2
     &     + 30.D0*b*t2**(-1)*vlwm*xnsq**2 + 36.D0*b*t2**(-1)*vlwm - 36.
     &    D0*b*t2**(-1)*vlwm**2 + 30.D0*b*t2**(-1)*vdw*xnsq**2 - 36.D0*
     &    b*t2**(-1)*vdw - 5.D0/2.D0*b*t2**(-1)*pi**2*xnsq**2 + 3.D0*b*
     &    t2**(-1)*pi**2 + 15.D0/2.D0*b*ro*vlpm**2 - 12.D0*b*ro*vlsm*
     &    vltm*xnsq**2 - 12.D0*b*ro*vlsm*vlwm*xnsq**2 + 6.D0*b*ro*
     &    vlsm**2*xnsq**2 )
      gcoeff6 = gcoeff6 + INVG * (  - 12.D0*b*ro*vltm*vlwm*xnsq - 21.D0
     &    *b*ro*vltm*xnsq + 15.D0*b*ro*vltm*xnsq**2 + 12.D0*b*ro*vltm
     &     - 18.D0*b*ro*vltm**2*xnsq - 6.D0*b*ro*vltm**2 + 21.D0*b*ro*
     &    vlwm*xnsq - 15.D0*b*ro*vlwm*xnsq**2 - 12.D0*b*ro*vlwm + 6.D0*
     &    b*ro*vlwm**2*xnsq + 36.D0*b*ro*vlwm**2 - 12.D0*b*ro*vdw*
     &    xnsq**2 + 36.D0*b*ro*vdw - 24.D0*b*ro*vdt*xnsq - 12.D0*b*ro*
     &    vdt*xnsq**2 - 6.D0*b*ro*vdt + 2.D0*b*ro*pi**2*xnsq + 2.D0*b*
     &    ro*pi**2*xnsq**2 - 5.D0/2.D0*b*ro*pi**2 + 15.D0/2.D0*b*ro**2*
     &    vltm**2*xnsq - 15.D0/2.D0*b*ro**2*vlwm**2*xnsq - 15.D0/2.D0*b
     &    *ro**2*vdw*xnsq + 15.D0/2.D0*b*ro**2*vdt*xnsq - 9.D0*b*
     &    vlpm**2 - 30.D0*b*vlsm*vlwm*xnsq**2 + 6.D0*b*vlsm*omro**(-1)*
     &    xnsq**2 - 36.D0*b*vlsm*xnsq**2 + 15.D0/2.D0*b*vlsm**2*xnsq**2
     &     + 48.D0*b*vltm*vlwm*xnsq + 36.D0*b*vltm**2*xnsq + 12.D0*b*
     &    vltm**2 + 117.D0*b*vlwm*xnsq**2 + 132.D0*b*vlwm + 12.D0*b*
     &    vlwm**2*xnsq - 48.D0*b*vlwm**2 + 36.D0*b*vdw*xnsq - 30.D0*b*
     &    vdw*xnsq**2 )
      gcoeff6 = gcoeff6 + INVG * (  - 48.D0*b*vdw + 60.D0*b*vdt*xnsq + 
     &    12.D0*b*vdt - 8.D0*b*pi**2*xnsq + 5.D0/2.D0*b*pi**2*xnsq**2
     &     + 3.D0*b*pi**2 - 12.D0*b*xnsq**2 + 12.D0*b - 18.D0*t1**(-1)*
     &    ro*vlpm*vlsm*xnsq + 21.D0*t1**(-1)*ro*vlpm*vlsm + 6.D0*
     &    t1**(-1)*ro*vlpm*vltm*xnsq + 6.D0*t1**(-1)*ro*vlpm*vltm + 30.D
     &    0*t1**(-1)*ro*vlpm*vlwm*xnsq - 48.D0*t1**(-1)*ro*vlpm*vlwm - 
     &    3.D0*t1**(-1)*ro*vlpm*xnsq + 6.D0*t1**(-1)*ro*vlpm + 9.D0*
     &    t1**(-1)*ro*vlpm**2*xnsq + 57.D0/8.D0*t1**(-1)*ro*vlpm**2*
     &    xnsq**2 - 21.D0/2.D0*t1**(-1)*ro*vlpm**2 + 36.D0*t1**(-1)*ro*
     &    vdmp*xnsq + 57.D0/2.D0*t1**(-1)*ro*vdmp*xnsq**2 - 42.D0*
     &    t1**(-1)*ro*vdmp + 3.D0*t1**(-1)*ro*pi**2*xnsq + 19.D0/8.D0*
     &    t1**(-1)*ro*pi**2*xnsq**2 - 7.D0/2.D0*t1**(-1)*ro*pi**2 + 18.D
     &    0*t1**(-1)*ro**2*vlpm*vlsm*xnsq - 6.D0*t1**(-1)*ro**2*vlpm*
     &    vlsm - 9.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 3.D0*t1**(-1)*
     &    ro**2*vlpm*vltm - 27.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 15.D0
     &    *t1**(-1)*ro**2*vlpm*vlwm )
      gcoeff6 = gcoeff6 + INVG * ( 6.D0*t1**(-1)*ro**2*vlpm*xnsq - 6.D0
     &    *t1**(-1)*ro**2*vlpm - 9.D0*t1**(-1)*ro**2*vlpm**2*xnsq - 57.D
     &    0/16.D0*t1**(-1)*ro**2*vlpm**2*xnsq**2 + 3.D0*t1**(-1)*ro**2*
     &    vlpm**2 - 36.D0*t1**(-1)*ro**2*vdmp*xnsq - 57.D0/4.D0*
     &    t1**(-1)*ro**2*vdmp*xnsq**2 + 12.D0*t1**(-1)*ro**2*vdmp - 3.D0
     &    *t1**(-1)*ro**2*pi**2*xnsq - 19.D0/16.D0*t1**(-1)*ro**2*pi**2
     &    *xnsq**2 + t1**(-1)*ro**2*pi**2 - 21.D0/4.D0*t1**(-1)*ro**3*
     &    vlpm*vlsm*xnsq + 21.D0/4.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq + 
     &    21.D0/4.D0*t1**(-1)*ro**3*vlpm*vlwm*xnsq - 3.D0*t1**(-1)*
     &    ro**3*vlpm*xnsq + 21.D0/8.D0*t1**(-1)*ro**3*vlpm**2*xnsq + 21.
     &    D0/2.D0*t1**(-1)*ro**3*vdmp*xnsq + 7.D0/8.D0*t1**(-1)*ro**3*
     &    pi**2*xnsq + 6.D0*t1**(-1)*vlpm*vlsm*xnsq - 12.D0*t1**(-1)*
     &    vlpm*vlsm - 12.D0*t1**(-1)*vlpm*vlwm*xnsq + 24.D0*t1**(-1)*
     &    vlpm*vlwm + 3.D0/8.D0*t1**(-1)*vlpm**2*omro**(-1)*xnsq**2 - 3.
     &    D0*t1**(-1)*vlpm**2*xnsq - 63.D0/8.D0*t1**(-1)*vlpm**2*
     &    xnsq**2 )
      gcoeff6 = gcoeff6 + INVG * ( 6.D0*t1**(-1)*vlpm**2 + 3.D0/2.D0*
     &    t1**(-1)*vdmp*omro**(-1)*xnsq**2 - 12.D0*t1**(-1)*vdmp*xnsq
     &     - 63.D0/2.D0*t1**(-1)*vdmp*xnsq**2 + 24.D0*t1**(-1)*vdmp + 1.
     &    D0/8.D0*t1**(-1)*omro**(-1)*pi**2*xnsq**2 - t1**(-1)*pi**2*
     &    xnsq - 21.D0/8.D0*t1**(-1)*pi**2*xnsq**2 + 2.D0*t1**(-1)*
     &    pi**2 + 39.D0*t2**(-1)*ro*vlpm*vlsm - 78.D0*t2**(-1)*ro*vlpm*
     &    vlwm + 9.D0*t2**(-1)*ro*vlpm + 57.D0/8.D0*t2**(-1)*ro*vlpm**2
     &    *xnsq**2 - 39.D0/2.D0*t2**(-1)*ro*vlpm**2 + 57.D0/2.D0*
     &    t2**(-1)*ro*vdmp*xnsq**2 - 78.D0*t2**(-1)*ro*vdmp + 19.D0/8.D0
     &    *t2**(-1)*ro*pi**2*xnsq**2 - 13.D0/2.D0*t2**(-1)*ro*pi**2 - 
     &    39.D0/2.D0*t2**(-1)*ro**2*vlpm*vlsm + 39.D0*t2**(-1)*ro**2*
     &    vlpm*vlwm - 9.D0*t2**(-1)*ro**2*vlpm - 63.D0/16.D0*t2**(-1)*
     &    ro**2*vlpm**2*xnsq**2 + 39.D0/4.D0*t2**(-1)*ro**2*vlpm**2 - 
     &    63.D0/4.D0*t2**(-1)*ro**2*vdmp*xnsq**2 + 39.D0*t2**(-1)*ro**2
     &    *vdmp - 21.D0/16.D0*t2**(-1)*ro**2*pi**2*xnsq**2 + 13.D0/4.D0
     &    *t2**(-1)*ro**2*pi**2 )
      gcoeff6 = gcoeff6 + INVG * (  - 18.D0*t2**(-1)*vlpm*vlsm + 36.D0*
     &    t2**(-1)*vlpm*vlwm + 3.D0/8.D0*t2**(-1)*vlpm**2*omro**(-1)*
     &    xnsq**2 - 63.D0/8.D0*t2**(-1)*vlpm**2*xnsq**2 + 9.D0*t2**(-1)
     &    *vlpm**2 + 3.D0/2.D0*t2**(-1)*vdmp*omro**(-1)*xnsq**2 - 63.D0/
     &    2.D0*t2**(-1)*vdmp*xnsq**2 + 36.D0*t2**(-1)*vdmp + 1.D0/8.D0*
     &    t2**(-1)*omro**(-1)*pi**2*xnsq**2 - 21.D0/8.D0*t2**(-1)*pi**2
     &    *xnsq**2 + 3.D0*t2**(-1)*pi**2 + 6.D0*ro*vlpm*vlsm*xnsq + 12.D
     &    0*ro*vlpm*vlsm + 18.D0*ro*vlpm*vltm*xnsq + 12.D0*ro*vlpm*vltm
     &     - 30.D0*ro*vlpm*vlwm*xnsq - 36.D0*ro*vlpm*vlwm + 36.D0*ro*
     &    vlpm - 3.D0*ro*vlpm**2*xnsq + 45.D0/4.D0*ro*vlpm**2*xnsq**2
     &     - 6.D0*ro*vlpm**2 - 12.D0*ro*vdmp*xnsq + 45.D0*ro*vdmp*
     &    xnsq**2 - 24.D0*ro*vdmp - ro*pi**2*xnsq + 15.D0/4.D0*ro*pi**2
     &    *xnsq**2 - 2.D0*ro*pi**2 - 6.D0*ro**2*vlpm*vlsm*xnsq - 12.D0*
     &    ro**2*vlpm*vltm*xnsq - 9.D0*ro**2*vlpm*vltm + 24.D0*ro**2*
     &    vlpm*vlwm*xnsq + 9.D0*ro**2*vlpm*vlwm + 3.D0*ro**2*vlpm**2*
     &    xnsq )
      gcoeff6 = gcoeff6 + INVG * ( 12.D0*ro**2*vdmp*xnsq + ro**2*pi**2*
     &    xnsq - 18.D0*vlpm*vlsm - 12.D0*vlpm*vltm*xnsq - 12.D0*vlpm*
     &    vltm + 12.D0*vlpm*vlwm*xnsq + 48.D0*vlpm*vlwm - 36.D0*vlpm - 
     &    3.D0/2.D0*vlpm**2*omro**(-1)*xnsq**2 + 9.D0*vlpm**2*xnsq**2
     &     + 9.D0*vlpm**2 - 6.D0*vdmp*omro**(-1)*xnsq**2 + 36.D0*vdmp*
     &    xnsq**2 + 36.D0*vdmp - 1.D0/2.D0*omro**(-1)*pi**2*xnsq**2 + 3.
     &    D0*pi**2*xnsq**2 + 3.D0*pi**2 )
      gcoeff6 = gcoeff6 + INVG**2 * ( 3.D0/4.D0*b*t1**(-1)*ro*vlpm**2*
     &    xnsq - 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2 + 15.D0/2.D0*b*
     &    t1**(-1)*ro*vlsm*vlwm*xnsq**2 - 15.D0/8.D0*b*t1**(-1)*ro*
     &    vlsm**2*xnsq**2 + 3.D0*b*t1**(-1)*ro*vlwm**2*xnsq - 6.D0*b*
     &    t1**(-1)*ro*vlwm**2 + 3.D0*b*t1**(-1)*ro*vdw*xnsq + 15.D0/2.D0
     &    *b*t1**(-1)*ro*vdw*xnsq**2 - 6.D0*b*t1**(-1)*ro*vdw - 1.D0/4.D
     &    0*b*t1**(-1)*ro*pi**2*xnsq - 5.D0/8.D0*b*t1**(-1)*ro*pi**2*
     &    xnsq**2 + 1.D0/2.D0*b*t1**(-1)*ro*pi**2 - 27.D0/16.D0*b*
     &    t1**(-1)*ro**2*vlpm**2*xnsq + 21.D0/16.D0*b*t1**(-1)*ro**2*
     &    vlpm**2 - 51.D0/8.D0*b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 51.D
     &    0/32.D0*b*t1**(-1)*ro**2*vlsm**2*xnsq**2 - 3.D0/4.D0*b*
     &    t1**(-1)*ro**2*vltm**2*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**2*
     &    vltm**2 - 6.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq + 6.D0*b*
     &    t1**(-1)*ro**2*vlwm**2 - 6.D0*b*t1**(-1)*ro**2*vdw*xnsq - 51.D
     &    0/8.D0*b*t1**(-1)*ro**2*vdw*xnsq**2 + 6.D0*b*t1**(-1)*ro**2*
     &    vdw )
      gcoeff6 = gcoeff6 + INVG**2 * (  - 3.D0/4.D0*b*t1**(-1)*ro**2*vdt
     &    *xnsq - 3.D0/4.D0*b*t1**(-1)*ro**2*vdt + 9.D0/16.D0*b*
     &    t1**(-1)*ro**2*pi**2*xnsq + 17.D0/32.D0*b*t1**(-1)*ro**2*
     &    pi**2*xnsq**2 - 7.D0/16.D0*b*t1**(-1)*ro**2*pi**2 + 15.D0/16.D
     &    0*b*t1**(-1)*ro**3*vlpm**2*xnsq + 3.D0/16.D0*b*t1**(-1)*ro**3
     &    *vlpm**2 + 3.D0/4.D0*b*t1**(-1)*ro**3*vlsm*vltm*xnsq**2 + 3.D0
     &    /4.D0*b*t1**(-1)*ro**3*vlsm*vlwm*xnsq**2 - 3.D0/8.D0*b*
     &    t1**(-1)*ro**3*vlsm**2*xnsq**2 - 3.D0/2.D0*b*t1**(-1)*ro**3*
     &    vltm*vlwm*xnsq + 3.D0/8.D0*b*t1**(-1)*ro**3*vltm**2*xnsq + 3.D
     &    0/2.D0*b*t1**(-1)*ro**3*vltm**2 + 21.D0/8.D0*b*t1**(-1)*ro**3
     &    *vlwm**2*xnsq - 3.D0/4.D0*b*t1**(-1)*ro**3*vlwm**2 + 15.D0/8.D
     &    0*b*t1**(-1)*ro**3*vdw*xnsq + 3.D0/4.D0*b*t1**(-1)*ro**3*vdw*
     &    xnsq**2 - 3.D0/4.D0*b*t1**(-1)*ro**3*vdw - 3.D0/8.D0*b*
     &    t1**(-1)*ro**3*vdt*xnsq + 3.D0/4.D0*b*t1**(-1)*ro**3*vdt*
     &    xnsq**2 + 3.D0/2.D0*b*t1**(-1)*ro**3*vdt + 1.D0/16.D0*b*
     &    t1**(-1)*ro**3*pi**2*xnsq )
      gcoeff6 = gcoeff6 + INVG**2 * (  - 1.D0/8.D0*b*t1**(-1)*ro**3*
     &    pi**2*xnsq**2 - 1.D0/16.D0*b*t1**(-1)*ro**3*pi**2 + 3.D0/4.D0
     &    *b*t2**(-2)*ro**3*vlwm**2 + 3.D0/4.D0*b*t2**(-2)*ro**3*vdw + 
     &    1.D0/8.D0*b*t2**(-2)*ro**3*pi**2 - 9.D0/4.D0*b*t2**(-1)*ro*
     &    vlpm**2 + 15.D0/2.D0*b*t2**(-1)*ro*vlsm*vlwm*xnsq**2 - 15.D0/
     &    8.D0*b*t2**(-1)*ro*vlsm**2*xnsq**2 - 9.D0*b*t2**(-1)*ro*
     &    vlwm**2 + 15.D0/2.D0*b*t2**(-1)*ro*vdw*xnsq**2 - 9.D0*b*
     &    t2**(-1)*ro*vdw - 5.D0/8.D0*b*t2**(-1)*ro*pi**2*xnsq**2 + 3.D0
     &    /4.D0*b*t2**(-1)*ro*pi**2 + 45.D0/16.D0*b*t2**(-1)*ro**2*
     &    vlpm**2 - 57.D0/8.D0*b*t2**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 57.D
     &    0/32.D0*b*t2**(-1)*ro**2*vlsm**2*xnsq**2 + 33.D0/4.D0*b*
     &    t2**(-1)*ro**2*vlwm**2 - 57.D0/8.D0*b*t2**(-1)*ro**2*vdw*
     &    xnsq**2 + 33.D0/4.D0*b*t2**(-1)*ro**2*vdw + 19.D0/32.D0*b*
     &    t2**(-1)*ro**2*pi**2*xnsq**2 - 23.D0/16.D0*b*t2**(-1)*ro**2*
     &    pi**2 - 9.D0/16.D0*b*t2**(-1)*ro**3*vlpm**2 + 3.D0/4.D0*b*
     &    t2**(-1)*ro**3*vlwm**2*xnsq )
      gcoeff6 = gcoeff6 + INVG**2 * (  - 9.D0/4.D0*b*t2**(-1)*ro**3*
     &    vlwm**2 + 3.D0/4.D0*b*t2**(-1)*ro**3*vdw*xnsq - 9.D0/4.D0*b*
     &    t2**(-1)*ro**3*vdw + 1.D0/8.D0*b*t2**(-1)*ro**3*pi**2*xnsq + 
     &    3.D0/16.D0*b*t2**(-1)*ro**3*pi**2 - 45.D0/4.D0*b*ro*vlpm**2
     &     + 57.D0/2.D0*b*ro*vlsm*vlwm*xnsq**2 - 57.D0/8.D0*b*ro*
     &    vlsm**2*xnsq**2 + 3.D0*b*ro*vltm**2*xnsq + 3.D0*b*ro*vltm**2
     &     - 3.D0*b*ro*vlwm**2*xnsq - 48.D0*b*ro*vlwm**2 - 3.D0*b*ro*
     &    vdw*xnsq + 57.D0/2.D0*b*ro*vdw*xnsq**2 - 48.D0*b*ro*vdw + 3.D0
     &    *b*ro*vdt*xnsq + 3.D0*b*ro*vdt - 19.D0/8.D0*b*ro*pi**2*
     &    xnsq**2 + 15.D0/4.D0*b*ro*pi**2 + 9.D0/4.D0*b*ro**2*vlpm**2
     &     + 3.D0/4.D0*b*ro**2*vlsm*vltm*xnsq**2 - 3.D0/4.D0*b*ro**2*
     &    vlsm*vlwm*xnsq**2 - 21.D0/4.D0*b*ro**2*vltm**2*xnsq - 9.D0/2.D
     &    0*b*ro**2*vltm**2 + 9.D0/4.D0*b*ro**2*vlwm**2*xnsq + 33.D0/2.D
     &    0*b*ro**2*vlwm**2 + 9.D0/4.D0*b*ro**2*vdw*xnsq - 3.D0/4.D0*b*
     &    ro**2*vdw*xnsq**2 + 33.D0/2.D0*b*ro**2*vdw - 21.D0/4.D0*b*
     &    ro**2*vdt*xnsq )
      gcoeff6 = gcoeff6 + INVG**2 * ( 3.D0/4.D0*b*ro**2*vdt*xnsq**2 - 9.
     &    D0/2.D0*b*ro**2*vdt - 1.D0/2.D0*b*ro**2*pi**2*xnsq - 1.D0/4.D0
     &    *b*ro**2*pi**2 + 3.D0/2.D0*b*ro**3*vltm**2*xnsq - 3.D0/2.D0*b
     &    *ro**3*vlwm**2*xnsq - 3.D0/2.D0*b*ro**3*vdw*xnsq + 3.D0/2.D0*
     &    b*ro**3*vdt*xnsq + 9.D0*b*vlpm**2 - 30.D0*b*vlsm*vlwm*xnsq**2
     &     + 15.D0/2.D0*b*vlsm**2*xnsq**2 + 36.D0*b*vlwm**2 - 30.D0*b*
     &    vdw*xnsq**2 + 36.D0*b*vdw + 5.D0/2.D0*b*pi**2*xnsq**2 - 3.D0*
     &    b*pi**2 + 3.D0/2.D0*t1**(-1)*ro*vlpm*vlsm*xnsq - 3.D0*
     &    t1**(-1)*ro*vlpm*vlsm - 3.D0*t1**(-1)*ro*vlpm*vlwm*xnsq + 6.D0
     &    *t1**(-1)*ro*vlpm*vlwm - 3.D0/4.D0*t1**(-1)*ro*vlpm**2*xnsq
     &     - 15.D0/8.D0*t1**(-1)*ro*vlpm**2*xnsq**2 + 3.D0/2.D0*
     &    t1**(-1)*ro*vlpm**2 - 3.D0*t1**(-1)*ro*vdmp*xnsq - 15.D0/2.D0
     &    *t1**(-1)*ro*vdmp*xnsq**2 + 6.D0*t1**(-1)*ro*vdmp - 1.D0/4.D0
     &    *t1**(-1)*ro*pi**2*xnsq - 5.D0/8.D0*t1**(-1)*ro*pi**2*xnsq**2
     &     + 1.D0/2.D0*t1**(-1)*ro*pi**2 - 33.D0/8.D0*t1**(-1)*ro**2*
     &    vlpm*vlsm*xnsq )
      gcoeff6 = gcoeff6 + INVG**2 * ( 33.D0/8.D0*t1**(-1)*ro**2*vlpm*
     &    vlsm + 3.D0/4.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq + 3.D0/4.D0*
     &    t1**(-1)*ro**2*vlpm*vltm + 15.D0/2.D0*t1**(-1)*ro**2*vlpm*
     &    vlwm*xnsq - 9.D0*t1**(-1)*ro**2*vlpm*vlwm + 33.D0/16.D0*
     &    t1**(-1)*ro**2*vlpm**2*xnsq + 81.D0/32.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq**2 - 33.D0/16.D0*t1**(-1)*ro**2*vlpm**2 + 33.D0/
     &    4.D0*t1**(-1)*ro**2*vdmp*xnsq + 81.D0/8.D0*t1**(-1)*ro**2*
     &    vdmp*xnsq**2 - 33.D0/4.D0*t1**(-1)*ro**2*vdmp + 11.D0/16.D0*
     &    t1**(-1)*ro**2*pi**2*xnsq + 27.D0/32.D0*t1**(-1)*ro**2*pi**2*
     &    xnsq**2 - 11.D0/16.D0*t1**(-1)*ro**2*pi**2 + 27.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm*vlsm*xnsq - 9.D0/8.D0*t1**(-1)*ro**3*vlpm
     &    *vlsm - 3.D0/2.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq - 3.D0/4.D0*
     &    t1**(-1)*ro**3*vlpm*vltm - 21.D0/4.D0*t1**(-1)*ro**3*vlpm*
     &    vlwm*xnsq + 3.D0*t1**(-1)*ro**3*vlpm*vlwm - 27.D0/16.D0*
     &    t1**(-1)*ro**3*vlpm**2*xnsq - 15.D0/32.D0*t1**(-1)*ro**3*
     &    vlpm**2*xnsq**2 )
      gcoeff6 = gcoeff6 + INVG**2 * ( 9.D0/16.D0*t1**(-1)*ro**3*vlpm**2
     &     - 27.D0/4.D0*t1**(-1)*ro**3*vdmp*xnsq - 15.D0/8.D0*t1**(-1)*
     &    ro**3*vdmp*xnsq**2 + 9.D0/4.D0*t1**(-1)*ro**3*vdmp - 9.D0/16.D
     &    0*t1**(-1)*ro**3*pi**2*xnsq - 5.D0/32.D0*t1**(-1)*ro**3*pi**2
     &    *xnsq**2 + 3.D0/16.D0*t1**(-1)*ro**3*pi**2 - 3.D0/4.D0*
     &    t1**(-1)*ro**4*vlpm*vlsm*xnsq + 3.D0/4.D0*t1**(-1)*ro**4*vlpm
     &    *vltm*xnsq + 3.D0/4.D0*t1**(-1)*ro**4*vlpm*vlwm*xnsq + 3.D0/8.
     &    D0*t1**(-1)*ro**4*vlpm**2*xnsq + 3.D0/2.D0*t1**(-1)*ro**4*
     &    vdmp*xnsq + 1.D0/8.D0*t1**(-1)*ro**4*pi**2*xnsq - 9.D0/2.D0*
     &    t2**(-1)*ro*vlpm*vlsm + 9.D0*t2**(-1)*ro*vlpm*vlwm - 15.D0/8.D
     &    0*t2**(-1)*ro*vlpm**2*xnsq**2 + 9.D0/4.D0*t2**(-1)*ro*vlpm**2
     &     - 15.D0/2.D0*t2**(-1)*ro*vdmp*xnsq**2 + 9.D0*t2**(-1)*ro*
     &    vdmp - 5.D0/8.D0*t2**(-1)*ro*pi**2*xnsq**2 + 3.D0/4.D0*
     &    t2**(-1)*ro*pi**2 + 63.D0/8.D0*t2**(-1)*ro**2*vlpm*vlsm - 63.D
     &    0/4.D0*t2**(-1)*ro**2*vlpm*vlwm + 87.D0/32.D0*t2**(-1)*ro**2*
     &    vlpm**2*xnsq**2 )
      gcoeff6 = gcoeff6 + INVG**2 * (  - 63.D0/16.D0*t2**(-1)*ro**2*
     &    vlpm**2 + 87.D0/8.D0*t2**(-1)*ro**2*vdmp*xnsq**2 - 63.D0/4.D0
     &    *t2**(-1)*ro**2*vdmp + 29.D0/32.D0*t2**(-1)*ro**2*pi**2*
     &    xnsq**2 - 21.D0/16.D0*t2**(-1)*ro**2*pi**2 - 27.D0/8.D0*
     &    t2**(-1)*ro**3*vlpm*vlsm + 27.D0/4.D0*t2**(-1)*ro**3*vlpm*
     &    vlwm - 21.D0/32.D0*t2**(-1)*ro**3*vlpm**2*xnsq**2 + 27.D0/16.D
     &    0*t2**(-1)*ro**3*vlpm**2 - 21.D0/8.D0*t2**(-1)*ro**3*vdmp*
     &    xnsq**2 + 27.D0/4.D0*t2**(-1)*ro**3*vdmp - 7.D0/32.D0*
     &    t2**(-1)*ro**3*pi**2*xnsq**2 + 9.D0/16.D0*t2**(-1)*ro**3*
     &    pi**2 - 63.D0/2.D0*ro*vlpm*vlsm - 3.D0*ro*vlpm*vltm*xnsq - 3.D
     &    0*ro*vlpm*vltm + 3.D0*ro*vlpm*vlwm*xnsq + 66.D0*ro*vlpm*vlwm
     &     - 87.D0/8.D0*ro*vlpm**2*xnsq**2 + 63.D0/4.D0*ro*vlpm**2 - 87.
     &    D0/2.D0*ro*vdmp*xnsq**2 + 63.D0*ro*vdmp - 29.D0/8.D0*ro*pi**2
     &    *xnsq**2 + 21.D0/4.D0*ro*pi**2 + 27.D0/2.D0*ro**2*vlpm*vlsm
     &     + 27.D0/4.D0*ro**2*vlpm*vltm*xnsq + 6.D0*ro**2*vlpm*vltm - 
     &    27.D0/4.D0*ro**2*vlpm*vlwm*xnsq )
      gcoeff6 = gcoeff6 + INVG**2 * (  - 33.D0*ro**2*vlpm*vlwm + 21.D0/
     &    8.D0*ro**2*vlpm**2*xnsq**2 - 27.D0/4.D0*ro**2*vlpm**2 + 21.D0/
     &    2.D0*ro**2*vdmp*xnsq**2 - 27.D0*ro**2*vdmp + 7.D0/8.D0*ro**2*
     &    pi**2*xnsq**2 - 9.D0/4.D0*ro**2*pi**2 - 15.D0/4.D0*ro**3*vlpm
     &    *vltm*xnsq - 3.D0*ro**3*vlpm*vltm + 15.D0/4.D0*ro**3*vlpm*
     &    vlwm*xnsq + 3.D0*ro**3*vlpm*vlwm + 18.D0*vlpm*vlsm - 36.D0*
     &    vlpm*vlwm + 15.D0/2.D0*vlpm**2*xnsq**2 - 9.D0*vlpm**2 + 30.D0
     &    *vdmp*xnsq**2 - 36.D0*vdmp + 5.D0/2.D0*pi**2*xnsq**2 - 3.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 - 64.D0*XLF*b*t1*TR*rmuom2*xn*xnsq + 32.D0*XLF*
     &    b*t2**(-1)*TR*rmuom2*xn*xnsq - 32.D0*XLF*b*t2**(-1)*TR*rmuom2
     &    *xn + 6.D0*b*t1**(-2)*ro*vltm**2*xnsq - 6.D0*b*t1**(-2)*ro*
     &    vltm**2 + 6.D0*b*t1**(-2)*ro*vdt*xnsq - 6.D0*b*t1**(-2)*ro*
     &    vdt + b*t1**(-2)*ro*pi**2*xnsq - b*t1**(-2)*ro*pi**2 + 21.D0/
     &    2.D0*b*t1**(-1)*ro*vlpm**2*xnsq + 3.D0/2.D0*b*t1**(-1)*ro*
     &    vlpm**2 + 12.D0*b*t1**(-1)*ro*vlsm*vltm*xnsq**2 + 12.D0*b*
     &    t1**(-1)*ro*vlsm*vlwm*xnsq**2 - 12.D0*b*t1**(-1)*ro*vlsm*
     &    xnsq**2 - 6.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2 - 6.D0*b*
     &    t1**(-1)*ro*vltm*vlwm*xnsq - 12.D0*b*t1**(-1)*ro*vltm*xnsq + 
     &    12.D0*b*t1**(-1)*ro*vltm*xnsq**2 + 3.D0*b*t1**(-1)*ro*vltm**2
     &    *xnsq - 12.D0*b*t1**(-1)*ro*vlwm*xnsq + 12.D0*b*t1**(-1)*ro*
     &    vlwm*xnsq**2 + 33.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 30.D0*b*
     &    t1**(-1)*ro*vdw*xnsq + 12.D0*b*t1**(-1)*ro*vdw*xnsq**2 + 12.D0
     &    *b*t1**(-1)*ro*vdt*xnsq**2 - 5.D0/2.D0*b*t1**(-1)*ro*pi**2*
     &    xnsq
      gcoeff6 = gcoeff6 - 2.D0*b*t1**(-1)*ro*pi**2*xnsq**2 - 3.D0/2.D0*
     &    b*t1**(-1)*ro*pi**2 - 3.D0*b*t1**(-1)*vlpm**2*xnsq + 9.D0*b*
     &    t1**(-1)*vlpm**2 + 42.D0*b*t1**(-1)*vlsm*vlwm*xnsq**2 - 6.D0*
     &    b*t1**(-1)*vlsm*omro**(-1)*xnsq**2 + 36.D0*b*t1**(-1)*vlsm*
     &    xnsq**2 - 21.D0/2.D0*b*t1**(-1)*vlsm**2*xnsq**2 - 48.D0*b*
     &    t1**(-1)*vltm*vlwm*xnsq + 168.D0*b*t1**(-1)*vltm*xnsq - 108.D0
     &    *b*t1**(-1)*vltm*xnsq**2 - 60.D0*b*t1**(-1)*vltm - 36.D0*b*
     &    t1**(-1)*vltm**2*xnsq - 12.D0*b*t1**(-1)*vltm**2 + 12.D0*b*
     &    t1**(-1)*vlwm*xnsq - 117.D0*b*t1**(-1)*vlwm*xnsq**2 - 120.D0*
     &    b*t1**(-1)*vlwm - 24.D0*b*t1**(-1)*vlwm**2*xnsq + 48.D0*b*
     &    t1**(-1)*vlwm**2 - 48.D0*b*t1**(-1)*vdw*xnsq + 42.D0*b*
     &    t1**(-1)*vdw*xnsq**2 + 48.D0*b*t1**(-1)*vdw - 60.D0*b*
     &    t1**(-1)*vdt*xnsq - 12.D0*b*t1**(-1)*vdt + 9.D0*b*t1**(-1)*
     &    pi**2*xnsq - 7.D0/2.D0*b*t1**(-1)*pi**2*xnsq**2 - 3.D0*b*
     &    t1**(-1)*pi**2 + 24.D0*b*t1**(-1)*xnsq - 12.D0*b*t1**(-1)*
     &    xnsq**2
      gcoeff6 = gcoeff6 - 12.D0*b*t1**(-1) - 48.D0*b*t1*vlsm*vltm*
     &    xnsq**2 - 48.D0*b*t1*vlsm*vlwm*xnsq**2 - 72.D0*b*t1*vlsm*
     &    omro**(-1)*xnsq**2 + 48.D0*b*t1*vlsm*xnsq**2 + 72.D0*b*t1*
     &    vltm*vlwm*xnsq - 24.D0*b*t1*vltm*xnsq + 48.D0*b*t1*vltm*
     &    xnsq**2 - 36.D0*b*t1*vltm**2*xnsq - 24.D0*b*t1*vlwm*xnsq + 48.
     &    D0*b*t1*vlwm*xnsq**2 - 36.D0*b*t1*vlwm**2*xnsq - 36.D0*b*t1*
     &    pi**2*xnsq + 24.D0*b*t1*pi**2*xnsq**2 + 176.D0*b*t1*rmuom2*
     &    xnsq**2 + 96.D0*b*t1*xnsq - 48.D0*b*t1*xnsq**2 + 12.D0*b*
     &    t2**(-2)*ro*vlwm**2*xnsq - 12.D0*b*t2**(-2)*ro*vlwm**2 + 12.D0
     &    *b*t2**(-2)*ro*vdw*xnsq - 12.D0*b*t2**(-2)*ro*vdw + 2.D0*b*
     &    t2**(-2)*ro*pi**2*xnsq - 2.D0*b*t2**(-2)*ro*pi**2 - 3.D0/2.D0
     &    *b*t2**(-1)*ro*vlpm**2 - 18.D0*b*t2**(-1)*ro*vlwm**2*xnsq - 
     &    18.D0*b*t2**(-1)*ro*vdw*xnsq - 3.D0*b*t2**(-1)*ro*pi**2*xnsq
     &     + 3.D0/2.D0*b*t2**(-1)*ro*pi**2 + 3.D0*b*t2**(-1)*vlpm**2 + 
     &    54.D0*b*t2**(-1)*vlsm*vlwm*xnsq**2 - 6.D0*b*t2**(-1)*vlsm*
     &    omro**(-1)*xnsq**2
      gcoeff6 = gcoeff6 + 36.D0*b*t2**(-1)*vlsm*xnsq**2 - 3.D0/2.D0*b*
     &    t2**(-1)*vlsm**2*xnsq**2 - 96.D0*b*t2**(-1)*vltm*vlwm*xnsq + 
     &    240.D0*b*t2**(-1)*vlwm*xnsq - 261.D0*b*t2**(-1)*vlwm*xnsq**2
     &     - 228.D0*b*t2**(-1)*vlwm + 12.D0*b*t2**(-1)*vlwm**2 + 6.D0*b
     &    *t2**(-1)*vdw*xnsq**2 + 12.D0*b*t2**(-1)*vdw + 24.D0*b*
     &    t2**(-1)*pi**2*xnsq - 25.D0/2.D0*b*t2**(-1)*pi**2*xnsq**2 - b
     &    *t2**(-1)*pi**2 + 88.D0*b*t2**(-1)*rmuom2*xnsq - 88.D0*b*
     &    t2**(-1)*rmuom2*xnsq**2 - 48.D0*b*t2**(-1)*xnsq + 12.D0*b*
     &    t2**(-1)*xnsq**2 + 36.D0*b*t2**(-1) - 12.D0*b*vlpm**2*xnsq - 
     &    12.D0*b*vlsm*vltm*xnsq**2 - 84.D0*b*vlsm*vlwm*xnsq**2 + 72.D0
     &    *b*vlsm*omro**(-1)*xnsq**2 - 48.D0*b*vlsm*xnsq**2 + 24.D0*b*
     &    vlsm**2*xnsq**2 + 72.D0*b*vltm*vlwm*xnsq + 12.D0*b*vltm*xnsq
     &     + 60.D0*b*vltm*xnsq**2 - 24.D0*b*vltm - 48.D0*b*vltm**2*xnsq
     &     - 156.D0*b*vlwm*xnsq + 132.D0*b*vlwm*xnsq**2 + 24.D0*b*vlwm
     &     - 72.D0*b*vlwm**2*xnsq - 36.D0*b*vdw*xnsq - 36.D0*b*vdw*
     &    xnsq**2
      gcoeff6 = gcoeff6 - 12.D0*b*vdt*xnsq - 60.D0*b*vdt*xnsq**2 - 32.D0
     &    *b*pi**2*xnsq + 8.D0*b*pi**2*xnsq**2 + 48.D0*b*xnsq**2 + 18.D0
     &    *t1**(-1)*ro*vlpm*vlsm*xnsq - 6.D0*t1**(-1)*ro*vlpm*vlsm - 12.
     &    D0*t1**(-1)*ro*vlpm*vltm*xnsq - 24.D0*t1**(-1)*ro*vlpm*vlwm*
     &    xnsq + 12.D0*t1**(-1)*ro*vlpm*vlwm + 24.D0*t1**(-1)*ro*vlpm*
     &    xnsq - 24.D0*t1**(-1)*ro*vlpm - 9.D0*t1**(-1)*ro*vlpm**2*xnsq
     &     - 27.D0/4.D0*t1**(-1)*ro*vlpm**2*xnsq**2 + 3.D0*t1**(-1)*ro*
     &    vlpm**2 - 36.D0*t1**(-1)*ro*vdmp*xnsq - 27.D0*t1**(-1)*ro*
     &    vdmp*xnsq**2 + 12.D0*t1**(-1)*ro*vdmp - 3.D0*t1**(-1)*ro*
     &    pi**2*xnsq - 9.D0/4.D0*t1**(-1)*ro*pi**2*xnsq**2 + t1**(-1)*
     &    ro*pi**2 - 9.D0*t1**(-1)*ro**2*vlpm*vlsm*xnsq + 9.D0*t1**(-1)
     &    *ro**2*vlpm*vltm*xnsq + 9.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq - 
     &    12.D0*t1**(-1)*ro**2*vlpm*xnsq + 9.D0/2.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq + 18.D0*t1**(-1)*ro**2*vdmp*xnsq + 3.D0/2.D0*
     &    t1**(-1)*ro**2*pi**2*xnsq - 6.D0*t1**(-1)*vlpm*vlsm*xnsq + 18.
     &    D0*t1**(-1)*vlpm*vlsm
      gcoeff6 = gcoeff6 + 12.D0*t1**(-1)*vlpm*vltm*xnsq + 12.D0*
     &    t1**(-1)*vlpm*vltm - 48.D0*t1**(-1)*vlpm*vlwm - 12.D0*
     &    t1**(-1)*vlpm*xnsq + 24.D0*t1**(-1)*vlpm + 3.D0/2.D0*t1**(-1)
     &    *vlpm**2*omro**(-1)*xnsq**2 + 3.D0*t1**(-1)*vlpm**2*xnsq - 12.
     &    D0*t1**(-1)*vlpm**2*xnsq**2 - 9.D0*t1**(-1)*vlpm**2 + 6.D0*
     &    t1**(-1)*vdmp*omro**(-1)*xnsq**2 + 12.D0*t1**(-1)*vdmp*xnsq
     &     - 48.D0*t1**(-1)*vdmp*xnsq**2 - 36.D0*t1**(-1)*vdmp + 1.D0/2.
     &    D0*t1**(-1)*omro**(-1)*pi**2*xnsq**2 + t1**(-1)*pi**2*xnsq - 
     &    4.D0*t1**(-1)*pi**2*xnsq**2 - 3.D0*t1**(-1)*pi**2 + 12.D0*t1*
     &    ro*vlpm*vlsm*xnsq - 12.D0*t1*ro*vlpm*vltm*xnsq - 12.D0*t1*ro*
     &    vlpm*vlwm*xnsq + 48.D0*t1*ro*vlpm*xnsq + 6.D0*t1*ro*vlpm**2*
     &    xnsq - 24.D0*t1*ro*vdmp*xnsq - 48.D0*t1*ro*vdmb*xnsq + 14.D0*
     &    t1*ro*pi**2*xnsq - 24.D0*t1*vlpm*xnsq + 18.D0*t1*vlpm**2*
     &    omro**(-1)*xnsq**2 - 24.D0*t1*vlpm**2*xnsq + 6.D0*t1*vlpm**2*
     &    xnsq**2 + 72.D0*t1*vdmp*omro**(-1)*xnsq**2 + 24.D0*t1*vdmp*
     &    xnsq**2
      gcoeff6 = gcoeff6 + 96.D0*t1*vdmb*xnsq + 6.D0*t1*omro**(-1)*pi**2
     &    *xnsq**2 - 32.D0*t1*pi**2*xnsq + 2.D0*t1*pi**2*xnsq**2 - 12.D0
     &    *t2**(-1)*ro*vlpm*vlsm + 24.D0*t2**(-1)*ro*vlpm*vlwm - 36.D0*
     &    t2**(-1)*ro*vlpm - 21.D0/4.D0*t2**(-1)*ro*vlpm**2*xnsq**2 + 
     &    12.D0*t2**(-1)*ro*vlpm**2 - 21.D0*t2**(-1)*ro*vdmp*xnsq**2 + 
     &    24.D0*t2**(-1)*ro*vdmp - 24.D0*t2**(-1)*ro*vdmb - 7.D0/4.D0*
     &    t2**(-1)*ro*pi**2*xnsq**2 + 10.D0*t2**(-1)*ro*pi**2 + 6.D0*
     &    t2**(-1)*vlpm*vlsm - 12.D0*t2**(-1)*vlpm*vlwm + 36.D0*
     &    t2**(-1)*vlpm + 3.D0/2.D0*t2**(-1)*vlpm**2*omro**(-1)*xnsq**2
     &     - 15.D0*t2**(-1)*vlpm**2*xnsq**2 - 15.D0*t2**(-1)*vlpm**2 + 
     &    6.D0*t2**(-1)*vdmp*omro**(-1)*xnsq**2 - 60.D0*t2**(-1)*vdmp*
     &    xnsq**2 - 12.D0*t2**(-1)*vdmp + 48.D0*t2**(-1)*vdmb + 1.D0/2.D
     &    0*t2**(-1)*omro**(-1)*pi**2*xnsq**2 - 5.D0*t2**(-1)*pi**2*
     &    xnsq**2 - 17.D0*t2**(-1)*pi**2 - 12.D0*ro*vlpm*vlsm*xnsq - 12.
     &    D0*ro*vlpm*vltm + 24.D0*ro*vlpm*vlwm*xnsq + 12.D0*ro*vlpm*
     &    vlwm
      gcoeff6 = gcoeff6 - 48.D0*ro*vlpm*xnsq + 6.D0*ro*vlpm**2*xnsq + 
     &    24.D0*ro*vdmp*xnsq + 2.D0*ro*pi**2*xnsq - 24.D0*vlpm*vlsm*
     &    xnsq + 12.D0*vlpm*vltm*xnsq + 36.D0*vlpm*vlwm*xnsq + 24.D0*
     &    vlpm*xnsq - 18.D0*vlpm**2*omro**(-1)*xnsq**2 + 12.D0*vlpm**2*
     &    xnsq + 42.D0*vlpm**2*xnsq**2 - 72.D0*vdmp*omro**(-1)*xnsq**2
     &     + 48.D0*vdmp*xnsq + 168.D0*vdmp*xnsq**2 - 6.D0*omro**(-1)*
     &    pi**2*xnsq**2 + 4.D0*pi**2*xnsq + 14.D0*pi**2*xnsq**2

      return
      end
