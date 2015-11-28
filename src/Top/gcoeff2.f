      double precision function gcoeff2()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff2 =  + TBAR**(-1) * ( 24.D0*b*t1**(-1)*vltm*xnsq - 12.D0*b*
     &    t1**(-1)*vltm*xnsq**2 - 12.D0*b*t1**(-1)*vltm - 12.D0*b*vltm*
     &    xnsq + 12.D0*b*vltm*xnsq**2 )
      gcoeff2 = gcoeff2 + UBAR**(-1)*INVG * ( 24.D0*b*t1**(-1)*vlwm*
     &    xnsq - 12.D0*b*t1**(-1)*vlwm*xnsq**2 - 12.D0*b*t1**(-1)*vlwm
     &     - 24.D0*b*vlwm*xnsq + 12.D0*b*vlwm*xnsq**2 + 12.D0*b*vlwm )
      gcoeff2 = gcoeff2 + UBAR**(-1) * ( 24.D0*b*t1**(-1)*vlwm*xnsq - 
     &    12.D0*b*t1**(-1)*vlwm*xnsq**2 - 12.D0*b*t1**(-1)*vlwm + 12.D0
     &    *b*vlwm*xnsq - 12.D0*b*vlwm*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG * (  - 6.D0*b*t1**(-1)*ro*vlpm**2*xnsq
     &     - 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2 + 6.D0*b*t1**(-1)*ro*vlsm*
     &    vltm*xnsq**2 + 12.D0*b*t1**(-1)*ro*vlsm*vlwm*xnsq**2 - 15.D0/
     &    2.D0*b*t1**(-1)*ro*vlsm*xnsq**2 - 9.D0/2.D0*b*t1**(-1)*ro*
     &    vlsm**2*xnsq**2 + 27.D0*b*t1**(-1)*ro*vltm*xnsq - 9.D0*b*
     &    t1**(-1)*ro*vltm*xnsq**2 - 6.D0*b*t1**(-1)*ro*vltm - 12.D0*b*
     &    t1**(-1)*ro*vltm**2*xnsq - 6.D0*b*t1**(-1)*ro*vltm**2 + 27.D0
     &    *b*t1**(-1)*ro*vlwm*xnsq + 3.D0*b*t1**(-1)*ro*vlwm*xnsq**2 - 
     &    36.D0*b*t1**(-1)*ro*vlwm - 12.D0*b*t1**(-1)*ro*vlwm**2*xnsq
     &     - 12.D0*b*t1**(-1)*ro*vdw*xnsq + 12.D0*b*t1**(-1)*ro*vdw*
     &    xnsq**2 - 12.D0*b*t1**(-1)*ro*vdt*xnsq + 6.D0*b*t1**(-1)*ro*
     &    vdt*xnsq**2 - 6.D0*b*t1**(-1)*ro*vdt + 2.D0*b*t1**(-1)*ro*
     &    pi**2*xnsq - 3.D0/2.D0*b*t1**(-1)*ro*pi**2*xnsq**2 + 1.D0/2.D0
     &    *b*t1**(-1)*ro*pi**2 - 12.D0*b*t1**(-1)*ro*xnsq + 6.D0*b*
     &    t1**(-1)*ro*xnsq**2 + 6.D0*b*t1**(-1)*ro + 9.D0/2.D0*b*
     &    t1**(-1)*ro**2*vlpm**2*xnsq )
      gcoeff2 = gcoeff2 + INVG * (  - 9.D0*b*t1**(-1)*ro**2*vltm*xnsq
     &     + 9.D0*b*t1**(-1)*ro**2*vltm**2*xnsq - 9.D0*b*t1**(-1)*ro**2
     &    *vlwm*xnsq + 21.D0/2.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq + 21.D0/
     &    2.D0*b*t1**(-1)*ro**2*vdw*xnsq + 9.D0*b*t1**(-1)*ro**2*vdt*
     &    xnsq - 5.D0/4.D0*b*t1**(-1)*ro**2*pi**2*xnsq + 3.D0*b*
     &    t1**(-1)*ro**2*xnsq + 3.D0/2.D0*b*t1**(-1)*vlsm*omro**(-1)*
     &    xnsq**2 - 3.D0/2.D0*b*t1**(-1)*vlsm*xnsq**2 - 24.D0*b*
     &    t1**(-1)*vltm*xnsq + 12.D0*b*t1**(-1)*vltm*xnsq**2 + 12.D0*b*
     &    t1**(-1)*vltm - 24.D0*b*t1**(-1)*vlwm*xnsq - 12.D0*b*t1**(-1)
     &    *vlwm*xnsq**2 + 36.D0*b*t1**(-1)*vlwm + 3.D0/2.D0*b*t2**(-2)*
     &    ro**2*vlwm**2 + 3.D0/2.D0*b*t2**(-2)*ro**2*vdw + 1.D0/4.D0*b*
     &    t2**(-2)*ro**2*pi**2 + 3.D0/2.D0*b*t2**(-1)*ro*vlpm**2 + 6.D0
     &    *b*t2**(-1)*ro*vlsm*vlwm*xnsq**2 - 6.D0*b*t2**(-1)*ro*vlsm*
     &    xnsq**2 - 3.D0/2.D0*b*t2**(-1)*ro*vlsm**2*xnsq**2 + 12.D0*b*
     &    t2**(-1)*ro*vlwm*xnsq**2 - 30.D0*b*t2**(-1)*ro*vlwm + 6.D0*b*
     &    t2**(-1)*ro*vdw*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG * (  - 1.D0/2.D0*b*t2**(-1)*ro*pi**2*
     &    xnsq**2 - 3.D0/2.D0*b*t2**(-1)*ro*pi**2 + 3.D0/2.D0*b*
     &    t2**(-1)*ro**2*vlwm**2*xnsq + 3.D0/2.D0*b*t2**(-1)*ro**2*vdw*
     &    xnsq + 1.D0/4.D0*b*t2**(-1)*ro**2*pi**2*xnsq - 24.D0*b*
     &    t2**(-1)*vlwm*xnsq**2 + 24.D0*b*t2**(-1)*vlwm + 3.D0/2.D0*b*
     &    ro*vlpm**2*xnsq + 6.D0*b*ro*vlsm*vltm*xnsq**2 - 6.D0*b*ro*
     &    vlsm*vlwm*xnsq**2 - 6.D0*b*ro*vltm*vlwm*xnsq - 6.D0*b*ro*vltm
     &    *xnsq + 12.D0*b*ro*vltm*xnsq**2 + 12.D0*b*ro*vltm - 3.D0*b*ro
     &    *vltm**2*xnsq + 6.D0*b*ro*vltm**2 + 18.D0*b*ro*vlwm*xnsq - 12.
     &    D0*b*ro*vlwm*xnsq**2 - 12.D0*b*ro*vlwm - 3.D0*b*ro*vlwm**2*
     &    xnsq - 6.D0*b*ro*vdw*xnsq - 6.D0*b*ro*vdw*xnsq**2 - 6.D0*b*ro
     &    *vdt*xnsq + 6.D0*b*ro*vdt*xnsq**2 + 6.D0*b*ro*vdt - 1.D0/2.D0
     &    *b*ro*pi**2*xnsq + b*ro*pi**2 - 6.D0*b*ro*xnsq - 3.D0*b*vlsm*
     &    omro**(-1)*xnsq**2 + 3.D0*b*vlsm*xnsq**2 + 24.D0*b*vltm*xnsq
     &     - 12.D0*b*vltm*xnsq**2 - 12.D0*b*vltm - 60.D0*b*vlwm*xnsq + 
     &    24.D0*b*vlwm*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG * ( 60.D0*b*vlwm + 24.D0*b*xnsq - 12.D0*
     &    b*xnsq**2 - 12.D0*b - 12.D0*t1**(-1)*ro*vlpm*vlsm*xnsq - 3.D0
     &    *t1**(-1)*ro*vlpm*vlsm + 12.D0*t1**(-1)*ro*vlpm*vltm*xnsq + 6.
     &    D0*t1**(-1)*ro*vlpm*vltm + 12.D0*t1**(-1)*ro*vlpm*vlwm*xnsq
     &     - 12.D0*t1**(-1)*ro*vlpm*xnsq - 6.D0*t1**(-1)*ro*vlpm + 6.D0
     &    *t1**(-1)*ro*vlpm**2*xnsq - 33.D0/8.D0*t1**(-1)*ro*vlpm**2*
     &    xnsq**2 + 3.D0/2.D0*t1**(-1)*ro*vlpm**2 + 24.D0*t1**(-1)*ro*
     &    vdmp*xnsq - 33.D0/2.D0*t1**(-1)*ro*vdmp*xnsq**2 + 6.D0*
     &    t1**(-1)*ro*vdmp + 2.D0*t1**(-1)*ro*pi**2*xnsq - 11.D0/8.D0*
     &    t1**(-1)*ro*pi**2*xnsq**2 + 1.D0/2.D0*t1**(-1)*ro*pi**2 + 15.D
     &    0*t1**(-1)*ro**2*vlpm*vlsm*xnsq + 3.D0*t1**(-1)*ro**2*vlpm*
     &    vlsm - 15.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 6.D0*t1**(-1)*
     &    ro**2*vlpm*vltm - 15.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 39.D0/
     &    2.D0*t1**(-1)*ro**2*vlpm*xnsq + 6.D0*t1**(-1)*ro**2*vlpm - 15.
     &    D0/2.D0*t1**(-1)*ro**2*vlpm**2*xnsq + 27.D0/8.D0*t1**(-1)*
     &    ro**2*vlpm**2*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG * (  - 3.D0/2.D0*t1**(-1)*ro**2*vlpm**2
     &     - 30.D0*t1**(-1)*ro**2*vdmp*xnsq + 27.D0/2.D0*t1**(-1)*ro**2
     &    *vdmp*xnsq**2 - 6.D0*t1**(-1)*ro**2*vdmp - 5.D0/2.D0*t1**(-1)
     &    *ro**2*pi**2*xnsq + 9.D0/8.D0*t1**(-1)*ro**2*pi**2*xnsq**2 - 
     &    1.D0/2.D0*t1**(-1)*ro**2*pi**2 - 3.D0*t1**(-1)*ro**3*vlpm*
     &    vlsm*xnsq + 3.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq + 3.D0*
     &    t1**(-1)*ro**3*vlpm*vlwm*xnsq - 15.D0/2.D0*t1**(-1)*ro**3*
     &    vlpm*xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm**2*xnsq + 6.D0*
     &    t1**(-1)*ro**3*vdmp*xnsq + 1.D0/2.D0*t1**(-1)*ro**3*pi**2*
     &    xnsq - 3.D0/8.D0*t1**(-1)*vlpm**2*omro**(-1)*xnsq**2 + 3.D0/8.
     &    D0*t1**(-1)*vlpm**2*xnsq**2 - 3.D0/2.D0*t1**(-1)*vdmp*
     &    omro**(-1)*xnsq**2 + 3.D0/2.D0*t1**(-1)*vdmp*xnsq**2 - 1.D0/8.
     &    D0*t1**(-1)*omro**(-1)*pi**2*xnsq**2 + 1.D0/8.D0*t1**(-1)*
     &    pi**2*xnsq**2 + 3.D0*t2**(-1)*ro*vlpm*vlsm - 6.D0*t2**(-1)*ro
     &    *vlpm*vlwm + 6.D0*t2**(-1)*ro*vlpm - 3.D0/2.D0*t2**(-1)*ro*
     &    vlpm**2*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG * (  - 3.D0/2.D0*t2**(-1)*ro*vlpm**2 - 6.
     &    D0*t2**(-1)*ro*vdmp*xnsq**2 - 6.D0*t2**(-1)*ro*vdmp - 1.D0/2.D
     &    0*t2**(-1)*ro*pi**2*xnsq**2 - 1.D0/2.D0*t2**(-1)*ro*pi**2 - 3.
     &    D0*t2**(-1)*ro**2*vlpm*vlsm + 6.D0*t2**(-1)*ro**2*vlpm*vlwm
     &     - 6.D0*t2**(-1)*ro**2*vlpm + 3.D0/2.D0*t2**(-1)*ro**2*
     &    vlpm**2 + 6.D0*t2**(-1)*ro**2*vdmp + 1.D0/2.D0*t2**(-1)*ro**2
     &    *pi**2 + 3.D0*ro*vlpm*vlsm*xnsq - 6.D0*ro*vlpm*vltm - 6.D0*ro
     &    *vlpm*vlwm*xnsq + 6.D0*ro*vlpm*vlwm - 3.D0*ro*vlpm*xnsq - 3.D0
     &    /2.D0*ro*vlpm**2*xnsq - 3.D0/4.D0*ro*vlpm**2*xnsq**2 - 6.D0*
     &    ro*vdmp*xnsq - 3.D0*ro*vdmp*xnsq**2 - 1.D0/2.D0*ro*pi**2*xnsq
     &     - 1.D0/4.D0*ro*pi**2*xnsq**2 - 3.D0*ro**2*vlpm*vlsm*xnsq + 6.
     &    D0*ro**2*vlpm*vltm + 6.D0*ro**2*vlpm*vlwm*xnsq - 6.D0*ro**2*
     &    vlpm*vlwm + 3.D0*ro**2*vlpm*xnsq + 3.D0/2.D0*ro**2*vlpm**2*
     &    xnsq + 6.D0*ro**2*vdmp*xnsq + 1.D0/2.D0*ro**2*pi**2*xnsq + 3.D
     &    0/4.D0*vlpm**2*omro**(-1)*xnsq**2 - 3.D0/4.D0*vlpm**2*xnsq**2
     &     + 3.D0*vdmp*omro**(-1)*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG * (  - 3.D0*vdmp*xnsq**2 + 1.D0/4.D0*
     &    omro**(-1)*pi**2*xnsq**2 - 1.D0/4.D0*pi**2*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG**2 * ( 3.D0*b*t1**(-1)*ro*vlpm**2*xnsq
     &     + 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2 - 6.D0*b*t1**(-1)*ro*vlsm*
     &    vlwm*xnsq**2 + 3.D0/2.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2 + 6.D0
     &    *b*t1**(-1)*ro*vltm**2*xnsq + 6.D0*b*t1**(-1)*ro*vltm**2 - 6.D
     &    0*b*t1**(-1)*ro*vlwm*xnsq**2 + 6.D0*b*t1**(-1)*ro*vlwm + 6.D0
     &    *b*t1**(-1)*ro*vlwm**2*xnsq + 6.D0*b*t1**(-1)*ro*vdw*xnsq - 6.
     &    D0*b*t1**(-1)*ro*vdw*xnsq**2 + 6.D0*b*t1**(-1)*ro*vdt*xnsq + 
     &    6.D0*b*t1**(-1)*ro*vdt - b*t1**(-1)*ro*pi**2*xnsq + 1.D0/2.D0
     &    *b*t1**(-1)*ro*pi**2*xnsq**2 - 1.D0/2.D0*b*t1**(-1)*ro*pi**2
     &     - 27.D0/4.D0*b*t1**(-1)*ro**2*vlpm**2*xnsq - 21.D0/8.D0*b*
     &    t1**(-1)*ro**2*vlpm**2 + 3.D0*b*t1**(-1)*ro**2*vlsm*vltm*
     &    xnsq**2 + 15.D0/2.D0*b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 - 3.D0
     &    /2.D0*b*t1**(-1)*ro**2*vlsm*xnsq**2 - 21.D0/8.D0*b*t1**(-1)*
     &    ro**2*vlsm**2*xnsq**2 + 9.D0/2.D0*b*t1**(-1)*ro**2*vltm*xnsq
     &     - 3.D0/2.D0*b*t1**(-1)*ro**2*vltm*xnsq**2 - 27.D0/2.D0*b*
     &    t1**(-1)*ro**2*vltm**2*xnsq )
      gcoeff2 = gcoeff2 + INVG**2 * (  - 21.D0/2.D0*b*t1**(-1)*ro**2*
     &    vltm**2 + 9.D0/2.D0*b*t1**(-1)*ro**2*vlwm*xnsq + 3.D0/2.D0*b*
     &    t1**(-1)*ro**2*vlwm*xnsq**2 - 6.D0*b*t1**(-1)*ro**2*vlwm - 27.
     &    D0/2.D0*b*t1**(-1)*ro**2*vlwm**2*xnsq - 27.D0/2.D0*b*t1**(-1)
     &    *ro**2*vdw*xnsq + 15.D0/2.D0*b*t1**(-1)*ro**2*vdw*xnsq**2 - 
     &    27.D0/2.D0*b*t1**(-1)*ro**2*vdt*xnsq + 3.D0*b*t1**(-1)*ro**2*
     &    vdt*xnsq**2 - 21.D0/2.D0*b*t1**(-1)*ro**2*vdt + 9.D0/4.D0*b*
     &    t1**(-1)*ro**2*pi**2*xnsq - 7.D0/8.D0*b*t1**(-1)*ro**2*pi**2*
     &    xnsq**2 + 7.D0/8.D0*b*t1**(-1)*ro**2*pi**2 + 57.D0/16.D0*b*
     &    t1**(-1)*ro**3*vlpm**2*xnsq + 15.D0/16.D0*b*t1**(-1)*ro**3*
     &    vlpm**2 - 9.D0/4.D0*b*t1**(-1)*ro**3*vltm*vlwm*xnsq - 3.D0/2.D
     &    0*b*t1**(-1)*ro**3*vltm*xnsq + 6.D0*b*t1**(-1)*ro**3*vltm**2*
     &    xnsq + 15.D0/4.D0*b*t1**(-1)*ro**3*vltm**2 - 3.D0/2.D0*b*
     &    t1**(-1)*ro**3*vlwm*xnsq + 6.D0*b*t1**(-1)*ro**3*vlwm**2*xnsq
     &     + 39.D0/8.D0*b*t1**(-1)*ro**3*vdw*xnsq + 39.D0/8.D0*b*
     &    t1**(-1)*ro**3*vdt*xnsq )
      gcoeff2 = gcoeff2 + INVG**2 * ( 15.D0/4.D0*b*t1**(-1)*ro**3*vdt
     &     - 13.D0/16.D0*b*t1**(-1)*ro**3*pi**2*xnsq - 5.D0/16.D0*b*
     &    t1**(-1)*ro**3*pi**2 - 3.D0/2.D0*b*t2**(-1)*ro*vlpm**2 - 6.D0
     &    *b*t2**(-1)*ro*vlsm*vlwm*xnsq**2 + 3.D0/2.D0*b*t2**(-1)*ro*
     &    vlsm**2*xnsq**2 - 6.D0*b*t2**(-1)*ro*vlwm*xnsq**2 + 6.D0*b*
     &    t2**(-1)*ro*vlwm - 6.D0*b*t2**(-1)*ro*vlwm**2 - 6.D0*b*
     &    t2**(-1)*ro*vdw*xnsq**2 - 6.D0*b*t2**(-1)*ro*vdw + 1.D0/2.D0*
     &    b*t2**(-1)*ro*pi**2*xnsq**2 + 1.D0/2.D0*b*t2**(-1)*ro*pi**2
     &     + 21.D0/8.D0*b*t2**(-1)*ro**2*vlpm**2 + 9.D0/2.D0*b*t2**(-1)
     &    *ro**2*vlsm*vlwm*xnsq**2 - 3.D0/2.D0*b*t2**(-1)*ro**2*vlsm*
     &    xnsq**2 - 9.D0/8.D0*b*t2**(-1)*ro**2*vlsm**2*xnsq**2 + 3.D0*b
     &    *t2**(-1)*ro**2*vlwm*xnsq**2 - 6.D0*b*t2**(-1)*ro**2*vlwm + 
     &    21.D0/2.D0*b*t2**(-1)*ro**2*vlwm**2 + 9.D0/2.D0*b*t2**(-1)*
     &    ro**2*vdw*xnsq**2 + 21.D0/2.D0*b*t2**(-1)*ro**2*vdw - 3.D0/8.D
     &    0*b*t2**(-1)*ro**2*pi**2*xnsq**2 - 7.D0/8.D0*b*t2**(-1)*ro**2
     &    *pi**2 )
      gcoeff2 = gcoeff2 + INVG**2 * (  - 15.D0/16.D0*b*t2**(-1)*ro**3*
     &    vlpm**2 - 15.D0/4.D0*b*t2**(-1)*ro**3*vlwm**2 - 15.D0/4.D0*b*
     &    t2**(-1)*ro**3*vdw + 5.D0/16.D0*b*t2**(-1)*ro**3*pi**2 + 3.D0
     &    *b*ro*vlpm**2*xnsq - 6.D0*b*ro*vlsm*vltm*xnsq**2 - 6.D0*b*ro*
     &    vlsm*vlwm*xnsq**2 + 6.D0*b*ro*vlsm*xnsq**2 + 3.D0*b*ro*
     &    vlsm**2*xnsq**2 - 6.D0*b*ro*vltm*xnsq - 6.D0*b*ro*vltm + 18.D0
     &    *b*ro*vltm**2*xnsq + 12.D0*b*ro*vltm**2 - 12.D0*b*ro*vlwm*
     &    xnsq - 6.D0*b*ro*vlwm*xnsq**2 + 30.D0*b*ro*vlwm - 6.D0*b*ro*
     &    vlwm**2*xnsq - 12.D0*b*ro*vlwm**2 - 6.D0*b*ro*vdw*xnsq - 6.D0
     &    *b*ro*vdw*xnsq**2 - 12.D0*b*ro*vdw + 18.D0*b*ro*vdt*xnsq - 6.D
     &    0*b*ro*vdt*xnsq**2 + 12.D0*b*ro*vdt - b*ro*pi**2*xnsq + b*ro*
     &    pi**2*xnsq**2 - 21.D0/8.D0*b*ro**2*vlpm**2*xnsq + 9.D0/2.D0*b
     &    *ro**2*vlsm*vltm*xnsq**2 - 9.D0/2.D0*b*ro**2*vlsm*vlwm*
     &    xnsq**2 + 3.D0/2.D0*b*ro**2*vltm*vlwm*xnsq + 3.D0*b*ro**2*
     &    vltm*xnsq**2 + 6.D0*b*ro**2*vltm - 81.D0/4.D0*b*ro**2*vltm**2
     &    *xnsq )
      gcoeff2 = gcoeff2 + INVG**2 * (  - 21.D0/2.D0*b*ro**2*vltm**2 + 6.
     &    D0*b*ro**2*vlwm*xnsq - 3.D0*b*ro**2*vlwm*xnsq**2 - 6.D0*b*
     &    ro**2*vlwm + 45.D0/4.D0*b*ro**2*vlwm**2*xnsq + 21.D0/2.D0*b*
     &    ro**2*vlwm**2 + 12.D0*b*ro**2*vdw*xnsq - 9.D0/2.D0*b*ro**2*
     &    vdw*xnsq**2 + 21.D0/2.D0*b*ro**2*vdw - 39.D0/2.D0*b*ro**2*vdt
     &    *xnsq + 9.D0/2.D0*b*ro**2*vdt*xnsq**2 - 21.D0/2.D0*b*ro**2*
     &    vdt + 5.D0/8.D0*b*ro**2*pi**2*xnsq + 15.D0/4.D0*b*ro**3*
     &    vltm**2*xnsq - 15.D0/4.D0*b*ro**3*vlwm**2*xnsq - 15.D0/4.D0*b
     &    *ro**3*vdw*xnsq + 15.D0/4.D0*b*ro**3*vdt*xnsq + 24.D0*b*vlwm*
     &    xnsq**2 - 24.D0*b*vlwm + 6.D0*t1**(-1)*ro*vlpm*vlsm*xnsq + 3.D
     &    0*t1**(-1)*ro*vlpm*vlsm - 6.D0*t1**(-1)*ro*vlpm*vltm*xnsq - 6.
     &    D0*t1**(-1)*ro*vlpm*vltm - 6.D0*t1**(-1)*ro*vlpm*vlwm*xnsq - 
     &    3.D0*t1**(-1)*ro*vlpm**2*xnsq + 3.D0/2.D0*t1**(-1)*ro*vlpm**2
     &    *xnsq**2 - 3.D0/2.D0*t1**(-1)*ro*vlpm**2 - 12.D0*t1**(-1)*ro*
     &    vdmp*xnsq + 6.D0*t1**(-1)*ro*vdmp*xnsq**2 - 6.D0*t1**(-1)*ro*
     &    vdmp )
      gcoeff2 = gcoeff2 + INVG**2 * (  - t1**(-1)*ro*pi**2*xnsq + 1.D0/
     &    2.D0*t1**(-1)*ro*pi**2*xnsq**2 - 1.D0/2.D0*t1**(-1)*ro*pi**2
     &     - 33.D0/2.D0*t1**(-1)*ro**2*vlpm*vlsm*xnsq - 27.D0/4.D0*
     &    t1**(-1)*ro**2*vlpm*vlsm + 33.D0/2.D0*t1**(-1)*ro**2*vlpm*
     &    vltm*xnsq + 27.D0/2.D0*t1**(-1)*ro**2*vlpm*vltm + 33.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm*vlwm*xnsq - 3.D0*t1**(-1)*ro**2*vlpm*xnsq
     &     - 3.D0/2.D0*t1**(-1)*ro**2*vlpm + 33.D0/4.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq - 27.D0/8.D0*t1**(-1)*ro**2*vlpm**2*xnsq**2 + 27.
     &    D0/8.D0*t1**(-1)*ro**2*vlpm**2 + 33.D0*t1**(-1)*ro**2*vdmp*
     &    xnsq - 27.D0/2.D0*t1**(-1)*ro**2*vdmp*xnsq**2 + 27.D0/2.D0*
     &    t1**(-1)*ro**2*vdmp + 11.D0/4.D0*t1**(-1)*ro**2*pi**2*xnsq - 
     &    9.D0/8.D0*t1**(-1)*ro**2*pi**2*xnsq**2 + 9.D0/8.D0*t1**(-1)*
     &    ro**2*pi**2 + 105.D0/8.D0*t1**(-1)*ro**3*vlpm*vlsm*xnsq + 15.D
     &    0/4.D0*t1**(-1)*ro**3*vlpm*vlsm - 105.D0/8.D0*t1**(-1)*ro**3*
     &    vlpm*vltm*xnsq - 15.D0/2.D0*t1**(-1)*ro**3*vlpm*vltm - 105.D0/
     &    8.D0*t1**(-1)*ro**3*vlpm*vlwm*xnsq )
      gcoeff2 = gcoeff2 + INVG**2 * ( 9.D0/2.D0*t1**(-1)*ro**3*vlpm*
     &    xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm - 105.D0/16.D0*t1**(-1)*
     &    ro**3*vlpm**2*xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm**2*xnsq**2
     &     - 15.D0/8.D0*t1**(-1)*ro**3*vlpm**2 - 105.D0/4.D0*t1**(-1)*
     &    ro**3*vdmp*xnsq + 6.D0*t1**(-1)*ro**3*vdmp*xnsq**2 - 15.D0/2.D
     &    0*t1**(-1)*ro**3*vdmp - 35.D0/16.D0*t1**(-1)*ro**3*pi**2*xnsq
     &     + 1.D0/2.D0*t1**(-1)*ro**3*pi**2*xnsq**2 - 5.D0/8.D0*
     &    t1**(-1)*ro**3*pi**2 - 21.D0/8.D0*t1**(-1)*ro**4*vlpm*vlsm*
     &    xnsq + 21.D0/8.D0*t1**(-1)*ro**4*vlpm*vltm*xnsq + 21.D0/8.D0*
     &    t1**(-1)*ro**4*vlpm*vlwm*xnsq - 3.D0/2.D0*t1**(-1)*ro**4*vlpm
     &    *xnsq + 21.D0/16.D0*t1**(-1)*ro**4*vlpm**2*xnsq + 21.D0/4.D0*
     &    t1**(-1)*ro**4*vdmp*xnsq + 7.D0/16.D0*t1**(-1)*ro**4*pi**2*
     &    xnsq - 3.D0*t2**(-1)*ro*vlpm*vlsm + 6.D0*t2**(-1)*ro*vlpm*
     &    vlwm + 3.D0/2.D0*t2**(-1)*ro*vlpm**2*xnsq**2 + 3.D0/2.D0*
     &    t2**(-1)*ro*vlpm**2 + 6.D0*t2**(-1)*ro*vdmp*xnsq**2 + 6.D0*
     &    t2**(-1)*ro*vdmp )
      gcoeff2 = gcoeff2 + INVG**2 * ( 1.D0/2.D0*t2**(-1)*ro*pi**2*
     &    xnsq**2 + 1.D0/2.D0*t2**(-1)*ro*pi**2 + 27.D0/4.D0*t2**(-1)*
     &    ro**2*vlpm*vlsm - 27.D0/2.D0*t2**(-1)*ro**2*vlpm*vlwm + 3.D0/
     &    2.D0*t2**(-1)*ro**2*vlpm - 15.D0/8.D0*t2**(-1)*ro**2*vlpm**2*
     &    xnsq**2 - 27.D0/8.D0*t2**(-1)*ro**2*vlpm**2 - 15.D0/2.D0*
     &    t2**(-1)*ro**2*vdmp*xnsq**2 - 27.D0/2.D0*t2**(-1)*ro**2*vdmp
     &     - 5.D0/8.D0*t2**(-1)*ro**2*pi**2*xnsq**2 - 9.D0/8.D0*
     &    t2**(-1)*ro**2*pi**2 - 15.D0/4.D0*t2**(-1)*ro**3*vlpm*vlsm + 
     &    15.D0/2.D0*t2**(-1)*ro**3*vlpm*vlwm - 3.D0/2.D0*t2**(-1)*
     &    ro**3*vlpm + 15.D0/8.D0*t2**(-1)*ro**3*vlpm**2 + 15.D0/2.D0*
     &    t2**(-1)*ro**3*vdmp + 5.D0/8.D0*t2**(-1)*ro**3*pi**2 + 6.D0*
     &    ro*vlpm*vlsm*xnsq - 18.D0*ro*vlpm*vltm*xnsq - 12.D0*ro*vlpm*
     &    vltm + 6.D0*ro*vlpm*vlwm*xnsq + 12.D0*ro*vlpm*vlwm + 6.D0*ro*
     &    vlpm*xnsq - 3.D0*ro*vlpm**2*xnsq + 3.D0*ro*vlpm**2*xnsq**2 - 
     &    12.D0*ro*vdmp*xnsq + 12.D0*ro*vdmp*xnsq**2 - ro*pi**2*xnsq + 
     &    ro*pi**2*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG**2 * (  - 33.D0/4.D0*ro**2*vlpm*vlsm*
     &    xnsq + 30.D0*ro**2*vlpm*vltm*xnsq + 33.D0/2.D0*ro**2*vlpm*
     &    vltm - 27.D0/2.D0*ro**2*vlpm*vlwm*xnsq - 33.D0/2.D0*ro**2*
     &    vlpm*vlwm - 9.D0*ro**2*vlpm*xnsq + 33.D0/8.D0*ro**2*vlpm**2*
     &    xnsq - 3.D0/2.D0*ro**2*vlpm**2*xnsq**2 + 33.D0/2.D0*ro**2*
     &    vdmp*xnsq - 6.D0*ro**2*vdmp*xnsq**2 + 11.D0/8.D0*ro**2*pi**2*
     &    xnsq - 1.D0/2.D0*ro**2*pi**2*xnsq**2 + 9.D0/4.D0*ro**3*vlpm*
     &    vlsm*xnsq - 12.D0*ro**3*vlpm*vltm*xnsq - 9.D0/2.D0*ro**3*vlpm
     &    *vltm + 15.D0/2.D0*ro**3*vlpm*vlwm*xnsq + 9.D0/2.D0*ro**3*
     &    vlpm*vlwm + 3.D0*ro**3*vlpm*xnsq - 9.D0/8.D0*ro**3*vlpm**2*
     &    xnsq - 9.D0/2.D0*ro**3*vdmp*xnsq - 3.D0/8.D0*ro**3*pi**2*xnsq
     &     )
      gcoeff2 = gcoeff2 + INVG**3 * ( 3.D0/4.D0*b*t1**(-1)*ro**2*
     &    vlpm**2*xnsq + 3.D0/8.D0*b*t1**(-1)*ro**2*vlpm**2 - 3.D0/2.D0
     &    *b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 3.D0/8.D0*b*t1**(-1)*
     &    ro**2*vlsm**2*xnsq**2 + 3.D0/2.D0*b*t1**(-1)*ro**2*vltm**2*
     &    xnsq + 3.D0/2.D0*b*t1**(-1)*ro**2*vltm**2 + 3.D0/2.D0*b*
     &    t1**(-1)*ro**2*vlwm**2*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**2*vdw*
     &    xnsq - 3.D0/2.D0*b*t1**(-1)*ro**2*vdw*xnsq**2 + 3.D0/2.D0*b*
     &    t1**(-1)*ro**2*vdt*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**2*vdt - 1.D
     &    0/4.D0*b*t1**(-1)*ro**2*pi**2*xnsq + 1.D0/8.D0*b*t1**(-1)*
     &    ro**2*pi**2*xnsq**2 - 1.D0/8.D0*b*t1**(-1)*ro**2*pi**2 - 21.D0
     &    /16.D0*b*t1**(-1)*ro**3*vlpm**2*xnsq - 9.D0/16.D0*b*t1**(-1)*
     &    ro**3*vlpm**2 + 3.D0/8.D0*b*t1**(-1)*ro**3*vlsm*vltm*xnsq**2
     &     + 9.D0/8.D0*b*t1**(-1)*ro**3*vlsm*vlwm*xnsq**2 - 3.D0/8.D0*b
     &    *t1**(-1)*ro**3*vlsm**2*xnsq**2 - 21.D0/8.D0*b*t1**(-1)*ro**3
     &    *vltm**2*xnsq - 9.D0/4.D0*b*t1**(-1)*ro**3*vltm**2 - 21.D0/8.D
     &    0*b*t1**(-1)*ro**3*vlwm**2*xnsq )
      gcoeff2 = gcoeff2 + INVG**3 * (  - 21.D0/8.D0*b*t1**(-1)*ro**3*
     &    vdw*xnsq + 9.D0/8.D0*b*t1**(-1)*ro**3*vdw*xnsq**2 - 21.D0/8.D0
     &    *b*t1**(-1)*ro**3*vdt*xnsq + 3.D0/8.D0*b*t1**(-1)*ro**3*vdt*
     &    xnsq**2 - 9.D0/4.D0*b*t1**(-1)*ro**3*vdt + 7.D0/16.D0*b*
     &    t1**(-1)*ro**3*pi**2*xnsq - 1.D0/8.D0*b*t1**(-1)*ro**3*pi**2*
     &    xnsq**2 + 3.D0/16.D0*b*t1**(-1)*ro**3*pi**2 + 9.D0/16.D0*b*
     &    t1**(-1)*ro**4*vlpm**2*xnsq + 3.D0/16.D0*b*t1**(-1)*ro**4*
     &    vlpm**2 - 3.D0/8.D0*b*t1**(-1)*ro**4*vltm*vlwm*xnsq + 15.D0/
     &    16.D0*b*t1**(-1)*ro**4*vltm**2*xnsq + 3.D0/4.D0*b*t1**(-1)*
     &    ro**4*vltm**2 + 15.D0/16.D0*b*t1**(-1)*ro**4*vlwm**2*xnsq + 3.
     &    D0/4.D0*b*t1**(-1)*ro**4*vdw*xnsq + 3.D0/4.D0*b*t1**(-1)*
     &    ro**4*vdt*xnsq + 3.D0/4.D0*b*t1**(-1)*ro**4*vdt - 1.D0/8.D0*b
     &    *t1**(-1)*ro**4*pi**2*xnsq - 1.D0/16.D0*b*t1**(-1)*ro**4*
     &    pi**2 - 3.D0/8.D0*b*t2**(-1)*ro**2*vlpm**2 - 3.D0/2.D0*b*
     &    t2**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 3.D0/8.D0*b*t2**(-1)*ro**2
     &    *vlsm**2*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG**3 * (  - 3.D0/2.D0*b*t2**(-1)*ro**2*
     &    vlwm**2 - 3.D0/2.D0*b*t2**(-1)*ro**2*vdw*xnsq**2 - 3.D0/2.D0*
     &    b*t2**(-1)*ro**2*vdw + 1.D0/8.D0*b*t2**(-1)*ro**2*pi**2*
     &    xnsq**2 + 1.D0/8.D0*b*t2**(-1)*ro**2*pi**2 + 9.D0/16.D0*b*
     &    t2**(-1)*ro**3*vlpm**2 + 3.D0/4.D0*b*t2**(-1)*ro**3*vlsm*vlwm
     &    *xnsq**2 - 3.D0/16.D0*b*t2**(-1)*ro**3*vlsm**2*xnsq**2 + 9.D0/
     &    4.D0*b*t2**(-1)*ro**3*vlwm**2 + 3.D0/4.D0*b*t2**(-1)*ro**3*
     &    vdw*xnsq**2 + 9.D0/4.D0*b*t2**(-1)*ro**3*vdw - 1.D0/16.D0*b*
     &    t2**(-1)*ro**3*pi**2*xnsq**2 - 3.D0/16.D0*b*t2**(-1)*ro**3*
     &    pi**2 - 3.D0/16.D0*b*t2**(-1)*ro**4*vlpm**2 - 3.D0/4.D0*b*
     &    t2**(-1)*ro**4*vlwm**2 - 3.D0/4.D0*b*t2**(-1)*ro**4*vdw + 1.D0
     &    /16.D0*b*t2**(-1)*ro**4*pi**2 - 3.D0/2.D0*b*ro*vlpm**2*xnsq
     &     + 6.D0*b*ro*vlsm*vlwm*xnsq**2 - 3.D0/2.D0*b*ro*vlsm**2*
     &    xnsq**2 - 6.D0*b*ro*vltm**2*xnsq - 6.D0*b*ro*vltm**2 + 6.D0*b
     &    *ro*vlwm**2 + 6.D0*b*ro*vdw*xnsq**2 + 6.D0*b*ro*vdw - 6.D0*b*
     &    ro*vdt*xnsq )
      gcoeff2 = gcoeff2 + INVG**3 * (  - 6.D0*b*ro*vdt + 1.D0/2.D0*b*ro
     &    *pi**2*xnsq - 1.D0/2.D0*b*ro*pi**2*xnsq**2 + 21.D0/8.D0*b*
     &    ro**2*vlpm**2*xnsq - 3.D0/2.D0*b*ro**2*vlsm*vltm*xnsq**2 - 3.D
     &    0*b*ro**2*vlsm*vlwm*xnsq**2 + 9.D0/8.D0*b*ro**2*vlsm**2*
     &    xnsq**2 + 12.D0*b*ro**2*vltm**2*xnsq + 21.D0/2.D0*b*ro**2*
     &    vltm**2 - 3.D0/2.D0*b*ro**2*vlwm**2*xnsq - 21.D0/2.D0*b*ro**2
     &    *vlwm**2 - 3.D0/2.D0*b*ro**2*vdw*xnsq - 3.D0*b*ro**2*vdw*
     &    xnsq**2 - 21.D0/2.D0*b*ro**2*vdw + 12.D0*b*ro**2*vdt*xnsq - 3.
     &    D0/2.D0*b*ro**2*vdt*xnsq**2 + 21.D0/2.D0*b*ro**2*vdt - 7.D0/8.
     &    D0*b*ro**2*pi**2*xnsq + 3.D0/8.D0*b*ro**2*pi**2*xnsq**2 - 9.D0
     &    /8.D0*b*ro**3*vlpm**2*xnsq + 3.D0/4.D0*b*ro**3*vlsm*vltm*
     &    xnsq**2 - 3.D0/4.D0*b*ro**3*vlsm*vlwm*xnsq**2 + 3.D0/4.D0*b*
     &    ro**3*vltm*vlwm*xnsq - 51.D0/8.D0*b*ro**3*vltm**2*xnsq - 9.D0/
     &    2.D0*b*ro**3*vltm**2 + 21.D0/8.D0*b*ro**3*vlwm**2*xnsq + 9.D0/
     &    2.D0*b*ro**3*vlwm**2 + 3.D0*b*ro**3*vdw*xnsq - 3.D0/4.D0*b*
     &    ro**3*vdw*xnsq**2 )
      gcoeff2 = gcoeff2 + INVG**3 * ( 9.D0/2.D0*b*ro**3*vdw - 6.D0*b*
     &    ro**3*vdt*xnsq + 3.D0/4.D0*b*ro**3*vdt*xnsq**2 - 9.D0/2.D0*b*
     &    ro**3*vdt + 1.D0/4.D0*b*ro**3*pi**2*xnsq + 3.D0/4.D0*b*ro**4*
     &    vltm**2*xnsq - 3.D0/4.D0*b*ro**4*vlwm**2*xnsq - 3.D0/4.D0*b*
     &    ro**4*vdw*xnsq + 3.D0/4.D0*b*ro**4*vdt*xnsq + 3.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm*vlsm*xnsq + 3.D0/4.D0*t1**(-1)*ro**2*vlpm
     &    *vlsm - 3.D0/2.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 3.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm*vltm - 3.D0/2.D0*t1**(-1)*ro**2*vlpm*vlwm
     &    *xnsq - 3.D0/4.D0*t1**(-1)*ro**2*vlpm**2*xnsq + 3.D0/8.D0*
     &    t1**(-1)*ro**2*vlpm**2*xnsq**2 - 3.D0/8.D0*t1**(-1)*ro**2*
     &    vlpm**2 - 3.D0*t1**(-1)*ro**2*vdmp*xnsq + 3.D0/2.D0*t1**(-1)*
     &    ro**2*vdmp*xnsq**2 - 3.D0/2.D0*t1**(-1)*ro**2*vdmp - 1.D0/4.D0
     &    *t1**(-1)*ro**2*pi**2*xnsq + 1.D0/8.D0*t1**(-1)*ro**2*pi**2*
     &    xnsq**2 - 1.D0/8.D0*t1**(-1)*ro**2*pi**2 - 27.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm*vlsm*xnsq - 3.D0/2.D0*t1**(-1)*ro**3*vlpm
     &    *vlsm )
      gcoeff2 = gcoeff2 + INVG**3 * ( 27.D0/8.D0*t1**(-1)*ro**3*vlpm*
     &    vltm*xnsq + 3.D0*t1**(-1)*ro**3*vlpm*vltm + 27.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm*vlwm*xnsq + 27.D0/16.D0*t1**(-1)*ro**3*
     &    vlpm**2*xnsq - 9.D0/16.D0*t1**(-1)*ro**3*vlpm**2*xnsq**2 + 3.D
     &    0/4.D0*t1**(-1)*ro**3*vlpm**2 + 27.D0/4.D0*t1**(-1)*ro**3*
     &    vdmp*xnsq - 9.D0/4.D0*t1**(-1)*ro**3*vdmp*xnsq**2 + 3.D0*
     &    t1**(-1)*ro**3*vdmp + 9.D0/16.D0*t1**(-1)*ro**3*pi**2*xnsq - 
     &    3.D0/16.D0*t1**(-1)*ro**3*pi**2*xnsq**2 + 1.D0/4.D0*t1**(-1)*
     &    ro**3*pi**2 + 9.D0/4.D0*t1**(-1)*ro**4*vlpm*vlsm*xnsq + 3.D0/
     &    4.D0*t1**(-1)*ro**4*vlpm*vlsm - 9.D0/4.D0*t1**(-1)*ro**4*vlpm
     &    *vltm*xnsq - 3.D0/2.D0*t1**(-1)*ro**4*vlpm*vltm - 9.D0/4.D0*
     &    t1**(-1)*ro**4*vlpm*vlwm*xnsq - 9.D0/8.D0*t1**(-1)*ro**4*
     &    vlpm**2*xnsq + 3.D0/16.D0*t1**(-1)*ro**4*vlpm**2*xnsq**2 - 3.D
     &    0/8.D0*t1**(-1)*ro**4*vlpm**2 - 9.D0/2.D0*t1**(-1)*ro**4*vdmp
     &    *xnsq + 3.D0/4.D0*t1**(-1)*ro**4*vdmp*xnsq**2 - 3.D0/2.D0*
     &    t1**(-1)*ro**4*vdmp )
      gcoeff2 = gcoeff2 + INVG**3 * (  - 3.D0/8.D0*t1**(-1)*ro**4*pi**2
     &    *xnsq + 1.D0/16.D0*t1**(-1)*ro**4*pi**2*xnsq**2 - 1.D0/8.D0*
     &    t1**(-1)*ro**4*pi**2 - 3.D0/8.D0*t1**(-1)*ro**5*vlpm*vlsm*
     &    xnsq + 3.D0/8.D0*t1**(-1)*ro**5*vlpm*vltm*xnsq + 3.D0/8.D0*
     &    t1**(-1)*ro**5*vlpm*vlwm*xnsq + 3.D0/16.D0*t1**(-1)*ro**5*
     &    vlpm**2*xnsq + 3.D0/4.D0*t1**(-1)*ro**5*vdmp*xnsq + 1.D0/16.D0
     &    *t1**(-1)*ro**5*pi**2*xnsq - 3.D0/4.D0*t2**(-1)*ro**2*vlpm*
     &    vlsm + 3.D0/2.D0*t2**(-1)*ro**2*vlpm*vlwm + 3.D0/8.D0*
     &    t2**(-1)*ro**2*vlpm**2*xnsq**2 + 3.D0/8.D0*t2**(-1)*ro**2*
     &    vlpm**2 + 3.D0/2.D0*t2**(-1)*ro**2*vdmp*xnsq**2 + 3.D0/2.D0*
     &    t2**(-1)*ro**2*vdmp + 1.D0/8.D0*t2**(-1)*ro**2*pi**2*xnsq**2
     &     + 1.D0/8.D0*t2**(-1)*ro**2*pi**2 + 3.D0/2.D0*t2**(-1)*ro**3*
     &    vlpm*vlsm - 3.D0*t2**(-1)*ro**3*vlpm*vlwm - 3.D0/8.D0*
     &    t2**(-1)*ro**3*vlpm**2*xnsq**2 - 3.D0/4.D0*t2**(-1)*ro**3*
     &    vlpm**2 - 3.D0/2.D0*t2**(-1)*ro**3*vdmp*xnsq**2 - 3.D0*
     &    t2**(-1)*ro**3*vdmp )
      gcoeff2 = gcoeff2 + INVG**3 * (  - 1.D0/8.D0*t2**(-1)*ro**3*pi**2
     &    *xnsq**2 - 1.D0/4.D0*t2**(-1)*ro**3*pi**2 - 3.D0/4.D0*
     &    t2**(-1)*ro**4*vlpm*vlsm + 3.D0/2.D0*t2**(-1)*ro**4*vlpm*vlwm
     &     + 3.D0/8.D0*t2**(-1)*ro**4*vlpm**2 + 3.D0/2.D0*t2**(-1)*
     &    ro**4*vdmp + 1.D0/8.D0*t2**(-1)*ro**4*pi**2 - 3.D0*ro*vlpm*
     &    vlsm*xnsq + 6.D0*ro*vlpm*vltm*xnsq + 6.D0*ro*vlpm*vltm - 6.D0
     &    *ro*vlpm*vlwm + 3.D0/2.D0*ro*vlpm**2*xnsq - 3.D0/2.D0*ro*
     &    vlpm**2*xnsq**2 + 6.D0*ro*vdmp*xnsq - 6.D0*ro*vdmp*xnsq**2 + 
     &    1.D0/2.D0*ro*pi**2*xnsq - 1.D0/2.D0*ro*pi**2*xnsq**2 + 27.D0/
     &    4.D0*ro**2*vlpm*vlsm*xnsq - 15.D0*ro**2*vlpm*vltm*xnsq - 27.D0
     &    /2.D0*ro**2*vlpm*vltm + 3.D0/2.D0*ro**2*vlpm*vlwm*xnsq + 27.D0
     &    /2.D0*ro**2*vlpm*vlwm - 27.D0/8.D0*ro**2*vlpm**2*xnsq + 15.D0/
     &    8.D0*ro**2*vlpm**2*xnsq**2 - 27.D0/2.D0*ro**2*vdmp*xnsq + 15.D
     &    0/2.D0*ro**2*vdmp*xnsq**2 - 9.D0/8.D0*ro**2*pi**2*xnsq + 5.D0/
     &    8.D0*ro**2*pi**2*xnsq**2 - 9.D0/2.D0*ro**3*vlpm*vlsm*xnsq + 
     &    12.D0*ro**3*vlpm*vltm*xnsq )
      gcoeff2 = gcoeff2 + INVG**3 * ( 9.D0*ro**3*vlpm*vltm - 3.D0*ro**3
     &    *vlpm*vlwm*xnsq - 9.D0*ro**3*vlpm*vlwm + 9.D0/4.D0*ro**3*
     &    vlpm**2*xnsq - 3.D0/8.D0*ro**3*vlpm**2*xnsq**2 + 9.D0*ro**3*
     &    vdmp*xnsq - 3.D0/2.D0*ro**3*vdmp*xnsq**2 + 3.D0/4.D0*ro**3*
     &    pi**2*xnsq - 1.D0/8.D0*ro**3*pi**2*xnsq**2 + 3.D0/4.D0*ro**4*
     &    vlpm*vlsm*xnsq - 3.D0*ro**4*vlpm*vltm*xnsq - 3.D0/2.D0*ro**4*
     &    vlpm*vltm + 3.D0/2.D0*ro**4*vlpm*vlwm*xnsq + 3.D0/2.D0*ro**4*
     &    vlpm*vlwm - 3.D0/8.D0*ro**4*vlpm**2*xnsq - 3.D0/2.D0*ro**4*
     &    vdmp*xnsq - 1.D0/8.D0*ro**4*pi**2*xnsq )
      gcoeff2 = gcoeff2 - 32.D0*XLF*b*t1*TR*xn*xnsq + 16.D0*XLF*b*TR*xn
     &    *xnsq - 6.D0*b*t1**(-2)*ro*vltm**2*xnsq + 6.D0*b*t1**(-2)*ro*
     &    vltm**2 - 6.D0*b*t1**(-2)*ro*vdt*xnsq + 6.D0*b*t1**(-2)*ro*
     &    vdt - b*t1**(-2)*ro*pi**2*xnsq + b*t1**(-2)*ro*pi**2 - 3.D0*b
     &    *t1**(-1)*ro*vlpm**2*xnsq - 3.D0*b*t1**(-1)*ro*vlpm**2 + 12.D0
     &    *b*t1**(-1)*ro*vltm*vlwm*xnsq - 12.D0*b*t1**(-1)*ro*vltm*xnsq
     &     - 6.D0*b*t1**(-1)*ro*vltm**2*xnsq - 12.D0*b*t1**(-1)*ro*vlwm
     &    *xnsq + 6.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 12.D0*b*t1**(-1)*ro
     &    *vdw*xnsq - b*t1**(-1)*ro*pi**2*xnsq + 3.D0*b*t1**(-1)*ro*
     &    pi**2 + 12.D0*b*t1**(-1)*ro*xnsq + 6.D0*b*t1**(-1)*vlsm*
     &    omro**(-1)*xnsq**2 - 6.D0*b*t1**(-1)*vlsm*xnsq**2 + 12.D0*b*
     &    t1**(-1)*vltm*xnsq - 12.D0*b*t1**(-1)*vltm*xnsq**2 + 36.D0*b*
     &    t1**(-1)*vlwm*xnsq - 12.D0*b*t1**(-1)*vlwm*xnsq**2 - 48.D0*b*
     &    t1**(-1)*vlwm - 24.D0*b*t1**(-1)*xnsq + 24.D0*b*t1**(-1)*
     &    xnsq**2 - 24.D0*b*t1*ro*vlpm**2*TR*xn*xnsq + 24.D0*b*t1*ro*TR
     &    *pi**2*xn*xnsq
      gcoeff2 = gcoeff2 - 192.D0*b*t1*ro*TR*xn*xnsq - 32.D0*b*t1*TR*xn*
     &    xnsq + 16.D0*b*t1*xnsq**2 + 6.D0*b*t2**(-2)*ro*vlwm**2*xnsq
     &     + 6.D0*b*t2**(-2)*ro*vdw*xnsq + b*t2**(-2)*ro*pi**2*xnsq + 3.
     &    D0*b*t2**(-1)*ro*vlpm**2 + 12.D0*b*t2**(-1)*ro*vlwm**2*xnsq
     &     + 12.D0*b*t2**(-1)*ro*vdw*xnsq + 2.D0*b*t2**(-1)*ro*pi**2*
     &    xnsq - 3.D0*b*t2**(-1)*ro*pi**2 + 24.D0*b*t2**(-1)*vlwm*xnsq
     &     - 48.D0*b*t2**(-1)*vlwm - 24.D0*b*t2**(-1)*xnsq + 24.D0*b*
     &    t2**(-1) + 12.D0*b*ro*vlpm**2*TR*xn*xnsq - 12.D0*b*ro*TR*
     &    pi**2*xn*xnsq + 96.D0*b*ro*TR*xn*xnsq + 16.D0*b*TR*xn*xnsq - 
     &    8.D0*b*xnsq**2 - 6.D0*t1**(-1)*ro*vlpm*vlsm*xnsq + 6.D0*
     &    t1**(-1)*ro*vlpm*vltm*xnsq + 6.D0*t1**(-1)*ro*vlpm*vlwm*xnsq
     &     + 6.D0*t1**(-1)*ro*vlpm*xnsq + 3.D0*t1**(-1)*ro*vlpm**2*xnsq
     &     + 3.D0/2.D0*t1**(-1)*ro*vlpm**2*xnsq**2 + 12.D0*t1**(-1)*ro*
     &    vdmp*xnsq + 6.D0*t1**(-1)*ro*vdmp*xnsq**2 + t1**(-1)*ro*pi**2
     &    *xnsq + 1.D0/2.D0*t1**(-1)*ro*pi**2*xnsq**2 + 6.D0*t1**(-1)*
     &    ro**2*vlpm*vlsm*xnsq
      gcoeff2 = gcoeff2 - 6.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 6.D0*
     &    t1**(-1)*ro**2*vlpm*vlwm*xnsq - 6.D0*t1**(-1)*ro**2*vlpm*xnsq
     &     - 3.D0*t1**(-1)*ro**2*vlpm**2*xnsq - 12.D0*t1**(-1)*ro**2*
     &    vdmp*xnsq - t1**(-1)*ro**2*pi**2*xnsq - 3.D0/2.D0*t1**(-1)*
     &    vlpm**2*omro**(-1)*xnsq**2 + 3.D0/2.D0*t1**(-1)*vlpm**2*
     &    xnsq**2 - 6.D0*t1**(-1)*vdmp*omro**(-1)*xnsq**2 + 6.D0*
     &    t1**(-1)*vdmp*xnsq**2 - 1.D0/2.D0*t1**(-1)*omro**(-1)*pi**2*
     &    xnsq**2 + 1.D0/2.D0*t1**(-1)*pi**2*xnsq**2 + 96.D0*t1*ro*vlpm
     &    *TR*xn*xnsq - 96.D0*t1*ro**2*vlpm*TR*xn*xnsq - 48.D0*ro*vlpm*
     &    TR*xn*xnsq + 48.D0*ro**2*vlpm*TR*xn*xnsq

      return
      end
