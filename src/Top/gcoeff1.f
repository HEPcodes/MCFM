      double precision function gcoeff1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff1 =  + TBAR**(-1) * ( 12.D0*b*t1*vltm*xnsq - 12.D0*b*t1*
     &    vltm*xnsq**2 - 24.D0*b*vltm*xnsq + 12.D0*b*vltm*xnsq**2 + 12.D
     &    0*b*vltm )
      gcoeff1 = gcoeff1 + UBAR**(-2) * ( 12.D0*b*t1*vlwm*xnsq**2 - 12.D0
     &    *b*t1*vlwm + 12.D0*b*t1**2*vlwm*xnsq - 12.D0*b*t1**2*vlwm*
     &    xnsq**2 - 12.D0*b*vlwm*xnsq + 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + UBAR**(-1)*INVG * ( 12.D0*b*t1**(-1)*vlwm*
     &    xnsq**2 - 12.D0*b*t1**(-1)*vlwm - 12.D0*b*vlwm*xnsq**2 + 12.D0
     &    *b*vlwm )
      gcoeff1 = gcoeff1 + UBAR**(-1) * ( 12.D0*b*t1**(-1)*vlwm*xnsq**2
     &     - 12.D0*b*t1**(-1)*vlwm + 84.D0*b*t1*vlwm*xnsq - 132.D0*b*t1
     &    *vlwm*xnsq**2 + 12.D0*b*t1*xnsq - 12.D0*b*t1*xnsq**2 + 12.D0*
     &    b*t2**(-1)*vlwm*xnsq**2 - 12.D0*b*t2**(-1)*vlwm + 132.D0*b*
     &    vlwm*xnsq - 84.D0*b*vlwm + 12.D0*b*xnsq - 12.D0*b )
      gcoeff1 = gcoeff1 + INVG * (  - 12.D0*b*t1**(-1)*ro*vlpm**2*xnsq
     &     - 9.D0/2.D0*b*t1**(-1)*ro*vlpm**2 + 6.D0*b*t1**(-1)*ro*vlsm*
     &    vltm*xnsq**2 + 12.D0*b*t1**(-1)*ro*vlsm*vlwm*xnsq**2 - 3.D0*b
     &    *t1**(-1)*ro*vlsm*xnsq**2 - 9.D0/2.D0*b*t1**(-1)*ro*vlsm**2*
     &    xnsq**2 + 12.D0*b*t1**(-1)*ro*vltm*vlwm*xnsq + 9.D0*b*
     &    t1**(-1)*ro*vltm*xnsq - 3.D0*b*t1**(-1)*ro*vltm*xnsq**2 - 18.D
     &    0*b*t1**(-1)*ro*vltm**2*xnsq - 18.D0*b*t1**(-1)*ro*vltm**2 + 
     &    9.D0*b*t1**(-1)*ro*vlwm*xnsq + 9.D0/2.D0*b*t1**(-1)*ro*vlwm*
     &    xnsq**2 - 12.D0*b*t1**(-1)*ro*vlwm - 18.D0*b*t1**(-1)*ro*
     &    vlwm**2*xnsq - 12.D0*b*t1**(-1)*ro*vdw*xnsq + 12.D0*b*
     &    t1**(-1)*ro*vdw*xnsq**2 - 12.D0*b*t1**(-1)*ro*vdt*xnsq + 6.D0
     &    *b*t1**(-1)*ro*vdt*xnsq**2 - 18.D0*b*t1**(-1)*ro*vdt + 2.D0*b
     &    *t1**(-1)*ro*pi**2*xnsq - 3.D0/2.D0*b*t1**(-1)*ro*pi**2*
     &    xnsq**2 + 3.D0/2.D0*b*t1**(-1)*ro*pi**2 + 3.D0/4.D0*b*
     &    t1**(-1)*ro*xnsq**2 + 6.D0*b*t1**(-1)*ro**2*vlpm**2*xnsq + 9.D
     &    0/4.D0*b*t1**(-1)*ro**2*vlpm**2 )
      gcoeff1 = gcoeff1 + INVG * (  - 3.D0*b*t1**(-1)*ro**2*vlsm*vltm*
     &    xnsq**2 - 3.D0*b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 + 3.D0/2.D0
     &    *b*t1**(-1)*ro**2*vlsm**2*xnsq**2 - 6.D0*b*t1**(-1)*ro**2*
     &    vltm*vlwm*xnsq - 3.D0*b*t1**(-1)*ro**2*vltm*xnsq + 21.D0/2.D0
     &    *b*t1**(-1)*ro**2*vltm**2*xnsq + 9.D0*b*t1**(-1)*ro**2*
     &    vltm**2 - 3.D0*b*t1**(-1)*ro**2*vlwm*xnsq + 15.D0/2.D0*b*
     &    t1**(-1)*ro**2*vlwm**2*xnsq + 9.D0/2.D0*b*t1**(-1)*ro**2*vdw*
     &    xnsq - 3.D0*b*t1**(-1)*ro**2*vdw*xnsq**2 + 15.D0/2.D0*b*
     &    t1**(-1)*ro**2*vdt*xnsq - 3.D0*b*t1**(-1)*ro**2*vdt*xnsq**2
     &     + 9.D0*b*t1**(-1)*ro**2*vdt - b*t1**(-1)*ro**2*pi**2*xnsq + 
     &    1.D0/2.D0*b*t1**(-1)*ro**2*pi**2*xnsq**2 - 3.D0/4.D0*b*
     &    t1**(-1)*ro**2*pi**2 + 6.D0*b*t1**(-1)*vlpm**2*xnsq + 3.D0*b*
     &    t1**(-1)*vlpm**2 - 12.D0*b*t1**(-1)*vlsm*vlwm*xnsq**2 + 3.D0*
     &    b*t1**(-1)*vlsm**2*xnsq**2 + 12.D0*b*t1**(-1)*vltm**2*xnsq + 
     &    12.D0*b*t1**(-1)*vltm**2 - 12.D0*b*t1**(-1)*vlwm*xnsq**2 + 12.
     &    D0*b*t1**(-1)*vlwm )
      gcoeff1 = gcoeff1 + INVG * ( 12.D0*b*t1**(-1)*vlwm**2*xnsq + 12.D0
     &    *b*t1**(-1)*vdw*xnsq - 12.D0*b*t1**(-1)*vdw*xnsq**2 + 12.D0*b
     &    *t1**(-1)*vdt*xnsq + 12.D0*b*t1**(-1)*vdt - 2.D0*b*t1**(-1)*
     &    pi**2*xnsq + b*t1**(-1)*pi**2*xnsq**2 - b*t1**(-1)*pi**2 + 6.D
     &    0*b*t2**(-1)*ro*vlpm**2 + 12.D0*b*t2**(-1)*ro*vlsm*vlwm*
     &    xnsq**2 - 3.D0*b*t2**(-1)*ro*vlsm*xnsq**2 - 3.D0*b*t2**(-1)*
     &    ro*vlsm**2*xnsq**2 + 15.D0/2.D0*b*t2**(-1)*ro*vlwm*xnsq**2 - 
     &    12.D0*b*t2**(-1)*ro*vlwm + 24.D0*b*t2**(-1)*ro*vlwm**2 + 12.D0
     &    *b*t2**(-1)*ro*vdw*xnsq**2 + 24.D0*b*t2**(-1)*ro*vdw - b*
     &    t2**(-1)*ro*pi**2*xnsq**2 - 2.D0*b*t2**(-1)*ro*pi**2 + 3.D0/4.
     &    D0*b*t2**(-1)*ro*xnsq**2 - 9.D0/4.D0*b*t2**(-1)*ro**2*vlpm**2
     &     - 9.D0*b*t2**(-1)*ro**2*vlwm**2 - 9.D0*b*t2**(-1)*ro**2*vdw
     &     + 3.D0/4.D0*b*t2**(-1)*ro**2*pi**2 - 3.D0*b*t2**(-1)*vlpm**2
     &     - 12.D0*b*t2**(-1)*vlsm*vlwm*xnsq**2 + 3.D0*b*t2**(-1)*
     &    vlsm**2*xnsq**2 - 12.D0*b*t2**(-1)*vlwm*xnsq**2 + 12.D0*b*
     &    t2**(-1)*vlwm )
      gcoeff1 = gcoeff1 + INVG * (  - 12.D0*b*t2**(-1)*vlwm**2 - 12.D0*
     &    b*t2**(-1)*vdw*xnsq**2 - 12.D0*b*t2**(-1)*vdw + b*t2**(-1)*
     &    pi**2*xnsq**2 + b*t2**(-1)*pi**2 - 3.D0*b*ro*vlpm**2*xnsq - 3.
     &    D0*b*ro*vlpm**2 + 24.D0*b*ro*vlsm*vltm*xnsq**2 - 6.D0*b*ro*
     &    vlsm**2*xnsq**2 + 12.D0*b*ro*vltm*vlwm*xnsq + 6.D0*b*ro*vltm*
     &    xnsq**2 + 12.D0*b*ro*vltm - 30.D0*b*ro*vltm**2*xnsq - 18.D0*b
     &    *ro*vltm**2 + 12.D0*b*ro*vlwm*xnsq - 6.D0*b*ro*vlwm*xnsq**2
     &     - 12.D0*b*ro*vlwm + 30.D0*b*ro*vlwm**2*xnsq + 6.D0*b*ro*
     &    vlwm**2 + 36.D0*b*ro*vdw*xnsq + 6.D0*b*ro*vdw - 24.D0*b*ro*
     &    vdt*xnsq + 24.D0*b*ro*vdt*xnsq**2 - 18.D0*b*ro*vdt - b*ro*
     &    pi**2*xnsq - 2.D0*b*ro*pi**2*xnsq**2 + b*ro*pi**2 + 9.D0*b*
     &    ro**2*vltm**2*xnsq - 9.D0*b*ro**2*vlwm**2*xnsq - 9.D0*b*ro**2
     &    *vdw*xnsq + 9.D0*b*ro**2*vdt*xnsq + 3.D0*b*vlpm**2*xnsq - 12.D
     &    0*b*vlsm*vltm*xnsq**2 + 12.D0*b*vlsm*xnsq**2 + 3.D0*b*vlsm**2
     &    *xnsq**2 - 48.D0*b*vltm*vlwm*xnsq - 12.D0*b*vltm*xnsq - 12.D0
     &    *b*vltm )
      gcoeff1 = gcoeff1 + INVG * ( 12.D0*b*vltm**2 - 24.D0*b*vlwm*xnsq
     &     - 18.D0*b*vlwm*xnsq**2 + 60.D0*b*vlwm - 36.D0*b*vlwm**2*xnsq
     &     - 12.D0*b*vlwm**2 - 60.D0*b*vdw*xnsq - 12.D0*b*vdw - 24.D0*b
     &    *vdt*xnsq - 12.D0*b*vdt*xnsq**2 + 12.D0*b*vdt + 7.D0*b*pi**2*
     &    xnsq + b*pi**2*xnsq**2 - 3.D0*b*xnsq**2 - 30.D0*t1**(-1)*ro*
     &    vlpm*vlsm*xnsq - 12.D0*t1**(-1)*ro*vlpm*vlsm + 30.D0*t1**(-1)
     &    *ro*vlpm*vltm*xnsq + 24.D0*t1**(-1)*ro*vlpm*vltm + 30.D0*
     &    t1**(-1)*ro*vlpm*vlwm*xnsq - 6.D0*t1**(-1)*ro*vlpm*xnsq - 3.D0
     &    *t1**(-1)*ro*vlpm + 15.D0*t1**(-1)*ro*vlpm**2*xnsq - 6.D0*
     &    t1**(-1)*ro*vlpm**2*xnsq**2 + 6.D0*t1**(-1)*ro*vlpm**2 + 60.D0
     &    *t1**(-1)*ro*vdmp*xnsq - 24.D0*t1**(-1)*ro*vdmp*xnsq**2 + 24.D
     &    0*t1**(-1)*ro*vdmp + 5.D0*t1**(-1)*ro*pi**2*xnsq - 2.D0*
     &    t1**(-1)*ro*pi**2*xnsq**2 + 2.D0*t1**(-1)*ro*pi**2 + 24.D0*
     &    t1**(-1)*ro**2*vlpm*vlsm*xnsq + 6.D0*t1**(-1)*ro**2*vlpm*vlsm
     &     - 24.D0*t1**(-1)*ro**2*vlpm*vltm*xnsq - 12.D0*t1**(-1)*ro**2
     &    *vlpm*vltm )
      gcoeff1 = gcoeff1 + INVG * (  - 24.D0*t1**(-1)*ro**2*vlpm*vlwm*
     &    xnsq + 9.D0*t1**(-1)*ro**2*vlpm*xnsq + 3.D0*t1**(-1)*ro**2*
     &    vlpm - 12.D0*t1**(-1)*ro**2*vlpm**2*xnsq + 21.D0/4.D0*
     &    t1**(-1)*ro**2*vlpm**2*xnsq**2 - 3.D0*t1**(-1)*ro**2*vlpm**2
     &     - 48.D0*t1**(-1)*ro**2*vdmp*xnsq + 21.D0*t1**(-1)*ro**2*vdmp
     &    *xnsq**2 - 12.D0*t1**(-1)*ro**2*vdmp - 4.D0*t1**(-1)*ro**2*
     &    pi**2*xnsq + 7.D0/4.D0*t1**(-1)*ro**2*pi**2*xnsq**2 - 
     &    t1**(-1)*ro**2*pi**2 - 6.D0*t1**(-1)*ro**3*vlpm*vlsm*xnsq + 6.
     &    D0*t1**(-1)*ro**3*vlpm*vltm*xnsq + 6.D0*t1**(-1)*ro**3*vlpm*
     &    vlwm*xnsq - 3.D0*t1**(-1)*ro**3*vlpm*xnsq + 3.D0*t1**(-1)*
     &    ro**3*vlpm**2*xnsq + 12.D0*t1**(-1)*ro**3*vdmp*xnsq + 
     &    t1**(-1)*ro**3*pi**2*xnsq + 12.D0*t1**(-1)*vlpm*vlsm*xnsq + 6.
     &    D0*t1**(-1)*vlpm*vlsm - 12.D0*t1**(-1)*vlpm*vltm*xnsq - 12.D0
     &    *t1**(-1)*vlpm*vltm - 12.D0*t1**(-1)*vlpm*vlwm*xnsq - 6.D0*
     &    t1**(-1)*vlpm**2*xnsq + 3.D0*t1**(-1)*vlpm**2*xnsq**2 - 3.D0*
     &    t1**(-1)*vlpm**2 )
      gcoeff1 = gcoeff1 + INVG * (  - 24.D0*t1**(-1)*vdmp*xnsq + 12.D0*
     &    t1**(-1)*vdmp*xnsq**2 - 12.D0*t1**(-1)*vdmp - 2.D0*t1**(-1)*
     &    pi**2*xnsq + t1**(-1)*pi**2*xnsq**2 - t1**(-1)*pi**2 + 9.D0*
     &    t2**(-1)*ro*vlpm*vlsm - 18.D0*t2**(-1)*ro*vlpm*vlwm + 3.D0*
     &    t2**(-1)*ro*vlpm - 3.D0/2.D0*t2**(-1)*ro*vlpm**2*xnsq**2 - 9.D
     &    0/2.D0*t2**(-1)*ro*vlpm**2 - 6.D0*t2**(-1)*ro*vdmp*xnsq**2 - 
     &    18.D0*t2**(-1)*ro*vdmp - 1.D0/2.D0*t2**(-1)*ro*pi**2*xnsq**2
     &     - 3.D0/2.D0*t2**(-1)*ro*pi**2 - 3.D0*t2**(-1)*ro**2*vlpm*
     &    vlsm + 6.D0*t2**(-1)*ro**2*vlpm*vlwm - 3.D0*t2**(-1)*ro**2*
     &    vlpm + 3.D0/2.D0*t2**(-1)*ro**2*vlpm**2*xnsq**2 + 3.D0/2.D0*
     &    t2**(-1)*ro**2*vlpm**2 + 6.D0*t2**(-1)*ro**2*vdmp*xnsq**2 + 6.
     &    D0*t2**(-1)*ro**2*vdmp + 1.D0/2.D0*t2**(-1)*ro**2*pi**2*
     &    xnsq**2 + 1.D0/2.D0*t2**(-1)*ro**2*pi**2 - 6.D0*t2**(-1)*vlpm
     &    *vlsm + 12.D0*t2**(-1)*vlpm*vlwm + 3.D0*t2**(-1)*vlpm**2*
     &    xnsq**2 + 3.D0*t2**(-1)*vlpm**2 + 12.D0*t2**(-1)*vdmp*xnsq**2
     &     + 12.D0*t2**(-1)*vdmp )
      gcoeff1 = gcoeff1 + INVG * ( t2**(-1)*pi**2*xnsq**2 + t2**(-1)*
     &    pi**2 - 15.D0*ro*vlpm*vlsm*xnsq + 48.D0*ro*vlpm*vltm*xnsq + 
     &    24.D0*ro*vlpm*vltm - 18.D0*ro*vlpm*vlwm*xnsq - 24.D0*ro*vlpm*
     &    vlwm - 18.D0*ro*vlpm*xnsq + 15.D0/2.D0*ro*vlpm**2*xnsq - 27.D0
     &    /2.D0*ro*vlpm**2*xnsq**2 + 30.D0*ro*vdmp*xnsq - 54.D0*ro*vdmp
     &    *xnsq**2 + 5.D0/2.D0*ro*pi**2*xnsq - 9.D0/2.D0*ro*pi**2*
     &    xnsq**2 + 9.D0*ro**2*vlpm*vlsm*xnsq - 24.D0*ro**2*vlpm*vltm*
     &    xnsq - 12.D0*ro**2*vlpm*vltm + 6.D0*ro**2*vlpm*vlwm*xnsq + 12.
     &    D0*ro**2*vlpm*vlwm + 6.D0*ro**2*vlpm*xnsq - 9.D0/2.D0*ro**2*
     &    vlpm**2*xnsq - 18.D0*ro**2*vdmp*xnsq - 3.D0/2.D0*ro**2*pi**2*
     &    xnsq + 6.D0*vlpm*vlsm*xnsq - 24.D0*vlpm*vltm*xnsq - 12.D0*
     &    vlpm*vltm + 12.D0*vlpm*vlwm*xnsq + 12.D0*vlpm*vlwm + 12.D0*
     &    vlpm*xnsq - 3.D0*vlpm**2*xnsq + 3.D0*vlpm**2*xnsq**2 - 12.D0*
     &    vdmp*xnsq + 12.D0*vdmp*xnsq**2 - pi**2*xnsq + pi**2*xnsq**2 )
      gcoeff1 = gcoeff1 + INVG**2 * ( 3.D0/2.D0*b*t1**(-1)*ro*vlpm**2*
     &    xnsq + 3.D0/4.D0*b*t1**(-1)*ro*vlpm**2 - 3.D0*b*t1**(-1)*ro*
     &    vlsm*vlwm*xnsq**2 + 3.D0/4.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2
     &     + 3.D0*b*t1**(-1)*ro*vltm**2*xnsq + 3.D0*b*t1**(-1)*ro*
     &    vltm**2 + 3.D0*b*t1**(-1)*ro*vlwm**2*xnsq + 3.D0*b*t1**(-1)*
     &    ro*vdw*xnsq - 3.D0*b*t1**(-1)*ro*vdw*xnsq**2 + 3.D0*b*
     &    t1**(-1)*ro*vdt*xnsq + 3.D0*b*t1**(-1)*ro*vdt - 1.D0/2.D0*b*
     &    t1**(-1)*ro*pi**2*xnsq + 1.D0/4.D0*b*t1**(-1)*ro*pi**2*
     &    xnsq**2 - 1.D0/4.D0*b*t1**(-1)*ro*pi**2 - 21.D0/8.D0*b*
     &    t1**(-1)*ro**2*vlpm**2*xnsq - 9.D0/8.D0*b*t1**(-1)*ro**2*
     &    vlpm**2 + 3.D0/4.D0*b*t1**(-1)*ro**2*vlsm*vltm*xnsq**2 + 9.D0/
     &    4.D0*b*t1**(-1)*ro**2*vlsm*vlwm*xnsq**2 - 3.D0/4.D0*b*
     &    t1**(-1)*ro**2*vlsm**2*xnsq**2 - 21.D0/4.D0*b*t1**(-1)*ro**2*
     &    vltm**2*xnsq - 9.D0/2.D0*b*t1**(-1)*ro**2*vltm**2 - 21.D0/4.D0
     &    *b*t1**(-1)*ro**2*vlwm**2*xnsq - 21.D0/4.D0*b*t1**(-1)*ro**2*
     &    vdw*xnsq )
      gcoeff1 = gcoeff1 + INVG**2 * ( 9.D0/4.D0*b*t1**(-1)*ro**2*vdw*
     &    xnsq**2 - 21.D0/4.D0*b*t1**(-1)*ro**2*vdt*xnsq + 3.D0/4.D0*b*
     &    t1**(-1)*ro**2*vdt*xnsq**2 - 9.D0/2.D0*b*t1**(-1)*ro**2*vdt
     &     + 7.D0/8.D0*b*t1**(-1)*ro**2*pi**2*xnsq - 1.D0/4.D0*b*
     &    t1**(-1)*ro**2*pi**2*xnsq**2 + 3.D0/8.D0*b*t1**(-1)*ro**2*
     &    pi**2 + 9.D0/8.D0*b*t1**(-1)*ro**3*vlpm**2*xnsq + 3.D0/8.D0*b
     &    *t1**(-1)*ro**3*vlpm**2 - 3.D0/4.D0*b*t1**(-1)*ro**3*vltm*
     &    vlwm*xnsq + 15.D0/8.D0*b*t1**(-1)*ro**3*vltm**2*xnsq + 3.D0/2.
     &    D0*b*t1**(-1)*ro**3*vltm**2 + 15.D0/8.D0*b*t1**(-1)*ro**3*
     &    vlwm**2*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**3*vdw*xnsq + 3.D0/2.D0
     &    *b*t1**(-1)*ro**3*vdt*xnsq + 3.D0/2.D0*b*t1**(-1)*ro**3*vdt
     &     - 1.D0/4.D0*b*t1**(-1)*ro**3*pi**2*xnsq - 1.D0/8.D0*b*
     &    t1**(-1)*ro**3*pi**2 - 3.D0/4.D0*b*t2**(-1)*ro*vlpm**2 - 3.D0
     &    *b*t2**(-1)*ro*vlsm*vlwm*xnsq**2 + 3.D0/4.D0*b*t2**(-1)*ro*
     &    vlsm**2*xnsq**2 - 3.D0*b*t2**(-1)*ro*vlwm**2 - 3.D0*b*
     &    t2**(-1)*ro*vdw*xnsq**2 )
      gcoeff1 = gcoeff1 + INVG**2 * (  - 3.D0*b*t2**(-1)*ro*vdw + 1.D0/
     &    4.D0*b*t2**(-1)*ro*pi**2*xnsq**2 + 1.D0/4.D0*b*t2**(-1)*ro*
     &    pi**2 + 9.D0/8.D0*b*t2**(-1)*ro**2*vlpm**2 + 3.D0/2.D0*b*
     &    t2**(-1)*ro**2*vlsm*vlwm*xnsq**2 - 3.D0/8.D0*b*t2**(-1)*ro**2
     &    *vlsm**2*xnsq**2 + 9.D0/2.D0*b*t2**(-1)*ro**2*vlwm**2 + 3.D0/
     &    2.D0*b*t2**(-1)*ro**2*vdw*xnsq**2 + 9.D0/2.D0*b*t2**(-1)*
     &    ro**2*vdw - 1.D0/8.D0*b*t2**(-1)*ro**2*pi**2*xnsq**2 - 3.D0/8.
     &    D0*b*t2**(-1)*ro**2*pi**2 - 3.D0/8.D0*b*t2**(-1)*ro**3*
     &    vlpm**2 - 3.D0/2.D0*b*t2**(-1)*ro**3*vlwm**2 - 3.D0/2.D0*b*
     &    t2**(-1)*ro**3*vdw + 1.D0/8.D0*b*t2**(-1)*ro**3*pi**2 + 21.D0/
     &    4.D0*b*ro*vlpm**2*xnsq - 3.D0*b*ro*vlsm*vltm*xnsq**2 - 6.D0*b
     &    *ro*vlsm*vlwm*xnsq**2 + 9.D0/4.D0*b*ro*vlsm**2*xnsq**2 + 24.D0
     &    *b*ro*vltm**2*xnsq + 21.D0*b*ro*vltm**2 - 3.D0*b*ro*vlwm**2*
     &    xnsq - 21.D0*b*ro*vlwm**2 - 3.D0*b*ro*vdw*xnsq - 6.D0*b*ro*
     &    vdw*xnsq**2 - 21.D0*b*ro*vdw + 24.D0*b*ro*vdt*xnsq - 3.D0*b*
     &    ro*vdt*xnsq**2 )
      gcoeff1 = gcoeff1 + INVG**2 * ( 21.D0*b*ro*vdt - 7.D0/4.D0*b*ro*
     &    pi**2*xnsq + 3.D0/4.D0*b*ro*pi**2*xnsq**2 - 9.D0/4.D0*b*ro**2
     &    *vlpm**2*xnsq + 3.D0/2.D0*b*ro**2*vlsm*vltm*xnsq**2 - 3.D0/2.D
     &    0*b*ro**2*vlsm*vlwm*xnsq**2 + 3.D0/2.D0*b*ro**2*vltm*vlwm*
     &    xnsq - 51.D0/4.D0*b*ro**2*vltm**2*xnsq - 9.D0*b*ro**2*vltm**2
     &     + 21.D0/4.D0*b*ro**2*vlwm**2*xnsq + 9.D0*b*ro**2*vlwm**2 + 6.
     &    D0*b*ro**2*vdw*xnsq - 3.D0/2.D0*b*ro**2*vdw*xnsq**2 + 9.D0*b*
     &    ro**2*vdw - 12.D0*b*ro**2*vdt*xnsq + 3.D0/2.D0*b*ro**2*vdt*
     &    xnsq**2 - 9.D0*b*ro**2*vdt + 1.D0/2.D0*b*ro**2*pi**2*xnsq + 3.
     &    D0/2.D0*b*ro**3*vltm**2*xnsq - 3.D0/2.D0*b*ro**3*vlwm**2*xnsq
     &     - 3.D0/2.D0*b*ro**3*vdw*xnsq + 3.D0/2.D0*b*ro**3*vdt*xnsq - 
     &    3.D0*b*vlpm**2*xnsq + 12.D0*b*vlsm*vlwm*xnsq**2 - 3.D0*b*
     &    vlsm**2*xnsq**2 - 12.D0*b*vltm**2*xnsq - 12.D0*b*vltm**2 + 12.
     &    D0*b*vlwm**2 + 12.D0*b*vdw*xnsq**2 + 12.D0*b*vdw - 12.D0*b*
     &    vdt*xnsq - 12.D0*b*vdt + b*pi**2*xnsq - b*pi**2*xnsq**2 + 3.D0
     &    *t1**(-1)*ro*vlpm*vlsm*xnsq )
      gcoeff1 = gcoeff1 + INVG**2 * ( 3.D0/2.D0*t1**(-1)*ro*vlpm*vlsm
     &     - 3.D0*t1**(-1)*ro*vlpm*vltm*xnsq - 3.D0*t1**(-1)*ro*vlpm*
     &    vltm - 3.D0*t1**(-1)*ro*vlpm*vlwm*xnsq - 3.D0/2.D0*t1**(-1)*
     &    ro*vlpm**2*xnsq + 3.D0/4.D0*t1**(-1)*ro*vlpm**2*xnsq**2 - 3.D0
     &    /4.D0*t1**(-1)*ro*vlpm**2 - 6.D0*t1**(-1)*ro*vdmp*xnsq + 3.D0
     &    *t1**(-1)*ro*vdmp*xnsq**2 - 3.D0*t1**(-1)*ro*vdmp - 1.D0/2.D0
     &    *t1**(-1)*ro*pi**2*xnsq + 1.D0/4.D0*t1**(-1)*ro*pi**2*xnsq**2
     &     - 1.D0/4.D0*t1**(-1)*ro*pi**2 - 27.D0/4.D0*t1**(-1)*ro**2*
     &    vlpm*vlsm*xnsq - 3.D0*t1**(-1)*ro**2*vlpm*vlsm + 27.D0/4.D0*
     &    t1**(-1)*ro**2*vlpm*vltm*xnsq + 6.D0*t1**(-1)*ro**2*vlpm*vltm
     &     + 27.D0/4.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq + 27.D0/8.D0*
     &    t1**(-1)*ro**2*vlpm**2*xnsq - 9.D0/8.D0*t1**(-1)*ro**2*
     &    vlpm**2*xnsq**2 + 3.D0/2.D0*t1**(-1)*ro**2*vlpm**2 + 27.D0/2.D
     &    0*t1**(-1)*ro**2*vdmp*xnsq - 9.D0/2.D0*t1**(-1)*ro**2*vdmp*
     &    xnsq**2 + 6.D0*t1**(-1)*ro**2*vdmp + 9.D0/8.D0*t1**(-1)*ro**2
     &    *pi**2*xnsq )
      gcoeff1 = gcoeff1 + INVG**2 * (  - 3.D0/8.D0*t1**(-1)*ro**2*pi**2
     &    *xnsq**2 + 1.D0/2.D0*t1**(-1)*ro**2*pi**2 + 9.D0/2.D0*
     &    t1**(-1)*ro**3*vlpm*vlsm*xnsq + 3.D0/2.D0*t1**(-1)*ro**3*vlpm
     &    *vlsm - 9.D0/2.D0*t1**(-1)*ro**3*vlpm*vltm*xnsq - 3.D0*
     &    t1**(-1)*ro**3*vlpm*vltm - 9.D0/2.D0*t1**(-1)*ro**3*vlpm*vlwm
     &    *xnsq - 9.D0/4.D0*t1**(-1)*ro**3*vlpm**2*xnsq + 3.D0/8.D0*
     &    t1**(-1)*ro**3*vlpm**2*xnsq**2 - 3.D0/4.D0*t1**(-1)*ro**3*
     &    vlpm**2 - 9.D0*t1**(-1)*ro**3*vdmp*xnsq + 3.D0/2.D0*t1**(-1)*
     &    ro**3*vdmp*xnsq**2 - 3.D0*t1**(-1)*ro**3*vdmp - 3.D0/4.D0*
     &    t1**(-1)*ro**3*pi**2*xnsq + 1.D0/8.D0*t1**(-1)*ro**3*pi**2*
     &    xnsq**2 - 1.D0/4.D0*t1**(-1)*ro**3*pi**2 - 3.D0/4.D0*t1**(-1)
     &    *ro**4*vlpm*vlsm*xnsq + 3.D0/4.D0*t1**(-1)*ro**4*vlpm*vltm*
     &    xnsq + 3.D0/4.D0*t1**(-1)*ro**4*vlpm*vlwm*xnsq + 3.D0/8.D0*
     &    t1**(-1)*ro**4*vlpm**2*xnsq + 3.D0/2.D0*t1**(-1)*ro**4*vdmp*
     &    xnsq + 1.D0/8.D0*t1**(-1)*ro**4*pi**2*xnsq - 3.D0/2.D0*
     &    t2**(-1)*ro*vlpm*vlsm )
      gcoeff1 = gcoeff1 + INVG**2 * ( 3.D0*t2**(-1)*ro*vlpm*vlwm + 3.D0/
     &    4.D0*t2**(-1)*ro*vlpm**2*xnsq**2 + 3.D0/4.D0*t2**(-1)*ro*
     &    vlpm**2 + 3.D0*t2**(-1)*ro*vdmp*xnsq**2 + 3.D0*t2**(-1)*ro*
     &    vdmp + 1.D0/4.D0*t2**(-1)*ro*pi**2*xnsq**2 + 1.D0/4.D0*
     &    t2**(-1)*ro*pi**2 + 3.D0*t2**(-1)*ro**2*vlpm*vlsm - 6.D0*
     &    t2**(-1)*ro**2*vlpm*vlwm - 3.D0/4.D0*t2**(-1)*ro**2*vlpm**2*
     &    xnsq**2 - 3.D0/2.D0*t2**(-1)*ro**2*vlpm**2 - 3.D0*t2**(-1)*
     &    ro**2*vdmp*xnsq**2 - 6.D0*t2**(-1)*ro**2*vdmp - 1.D0/4.D0*
     &    t2**(-1)*ro**2*pi**2*xnsq**2 - 1.D0/2.D0*t2**(-1)*ro**2*pi**2
     &     - 3.D0/2.D0*t2**(-1)*ro**3*vlpm*vlsm + 3.D0*t2**(-1)*ro**3*
     &    vlpm*vlwm + 3.D0/4.D0*t2**(-1)*ro**3*vlpm**2 + 3.D0*t2**(-1)*
     &    ro**3*vdmp + 1.D0/4.D0*t2**(-1)*ro**3*pi**2 + 27.D0/2.D0*ro*
     &    vlpm*vlsm*xnsq - 30.D0*ro*vlpm*vltm*xnsq - 27.D0*ro*vlpm*vltm
     &     + 3.D0*ro*vlpm*vlwm*xnsq + 27.D0*ro*vlpm*vlwm - 27.D0/4.D0*
     &    ro*vlpm**2*xnsq + 15.D0/4.D0*ro*vlpm**2*xnsq**2 - 27.D0*ro*
     &    vdmp*xnsq )
      gcoeff1 = gcoeff1 + INVG**2 * ( 15.D0*ro*vdmp*xnsq**2 - 9.D0/4.D0
     &    *ro*pi**2*xnsq + 5.D0/4.D0*ro*pi**2*xnsq**2 - 9.D0*ro**2*vlpm
     &    *vlsm*xnsq + 24.D0*ro**2*vlpm*vltm*xnsq + 18.D0*ro**2*vlpm*
     &    vltm - 6.D0*ro**2*vlpm*vlwm*xnsq - 18.D0*ro**2*vlpm*vlwm + 9.D
     &    0/2.D0*ro**2*vlpm**2*xnsq - 3.D0/4.D0*ro**2*vlpm**2*xnsq**2
     &     + 18.D0*ro**2*vdmp*xnsq - 3.D0*ro**2*vdmp*xnsq**2 + 3.D0/2.D0
     &    *ro**2*pi**2*xnsq - 1.D0/4.D0*ro**2*pi**2*xnsq**2 + 3.D0/2.D0
     &    *ro**3*vlpm*vlsm*xnsq - 6.D0*ro**3*vlpm*vltm*xnsq - 3.D0*
     &    ro**3*vlpm*vltm + 3.D0*ro**3*vlpm*vlwm*xnsq + 3.D0*ro**3*vlpm
     &    *vlwm - 3.D0/4.D0*ro**3*vlpm**2*xnsq - 3.D0*ro**3*vdmp*xnsq
     &     - 1.D0/4.D0*ro**3*pi**2*xnsq - 6.D0*vlpm*vlsm*xnsq + 12.D0*
     &    vlpm*vltm*xnsq + 12.D0*vlpm*vltm - 12.D0*vlpm*vlwm + 3.D0*
     &    vlpm**2*xnsq - 3.D0*vlpm**2*xnsq**2 + 12.D0*vdmp*xnsq - 12.D0
     &    *vdmp*xnsq**2 + pi**2*xnsq - pi**2*xnsq**2 )
      gcoeff1 = gcoeff1 + 64.D0*XLF*b*t1*TR*rmuom2*xn*xnsq + 32.D0*XLF*
     &    b*t1*TR*xn*xnsq - 32.D0*XLF*b*t2**(-1)*TR*rmuom2*xn*xnsq + 32.
     &    D0*XLF*b*t2**(-1)*TR*rmuom2*xn - 16.D0*XLF*b*TR*xn*xnsq + 6.D0
     &    *b*t1**(-1)*ro*vlpm**2*xnsq + 3.D0*b*t1**(-1)*ro*vlpm**2 - 12.
     &    D0*b*t1**(-1)*ro*vlsm*vltm*xnsq**2 - 12.D0*b*t1**(-1)*ro*vlsm
     &    *vlwm*xnsq**2 + 6.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2 - 12.D0*b*
     &    t1**(-1)*ro*vltm*vlwm*xnsq - 12.D0*b*t1**(-1)*ro*vltm*xnsq + 
     &    24.D0*b*t1**(-1)*ro*vltm**2*xnsq - 12.D0*b*t1**(-1)*ro*vlwm*
     &    xnsq - 6.D0*b*t1**(-1)*ro*vdw*xnsq - 12.D0*b*t1**(-1)*ro*vdw*
     &    xnsq**2 + 18.D0*b*t1**(-1)*ro*vdt*xnsq - 12.D0*b*t1**(-1)*ro*
     &    vdt*xnsq**2 + 2.D0*b*t1**(-1)*ro*pi**2*xnsq + 2.D0*b*t1**(-1)
     &    *ro*pi**2*xnsq**2 - 3.D0*b*t1**(-1)*ro*pi**2 - 6.D0*b*
     &    t1**(-1)*vlpm**2*xnsq + 12.D0*b*t1**(-1)*vlsm*vltm*xnsq**2 + 
     &    12.D0*b*t1**(-1)*vlsm*vlwm*xnsq**2 - 12.D0*b*t1**(-1)*vlsm*
     &    xnsq**2 - 6.D0*b*t1**(-1)*vlsm**2*xnsq**2 + 48.D0*b*t1**(-1)*
     &    vltm*vlwm*xnsq
      gcoeff1 = gcoeff1 + 36.D0*b*t1**(-1)*vltm*xnsq - 12.D0*b*t1**(-1)
     &    *vltm*xnsq**2 + 12.D0*b*t1**(-1)*vltm**2*xnsq + 36.D0*b*
     &    t1**(-1)*vlwm*xnsq + 18.D0*b*t1**(-1)*vlwm*xnsq**2 - 48.D0*b*
     &    t1**(-1)*vlwm + 12.D0*b*t1**(-1)*vlwm**2*xnsq + 36.D0*b*
     &    t1**(-1)*vdw*xnsq + 12.D0*b*t1**(-1)*vdw*xnsq**2 + 36.D0*b*
     &    t1**(-1)*vdt*xnsq + 12.D0*b*t1**(-1)*vdt*xnsq**2 - 6.D0*b*
     &    t1**(-1)*pi**2*xnsq - 2.D0*b*t1**(-1)*pi**2*xnsq**2 + 3.D0*b*
     &    t1**(-1)*xnsq**2 + 24.D0*b*t1*ro*vlpm**2*TR*xn*xnsq - 24.D0*b
     &    *t1*ro*TR*pi**2*xn*xnsq + 192.D0*b*t1*ro*TR*xn*xnsq + 48.D0*b
     &    *t1*vlsm*vltm*xnsq**2 + 48.D0*b*t1*vlsm*vlwm*xnsq**2 + 24.D0*
     &    b*t1*vlsm*omro**(-1)*xnsq**2 - 96.D0*b*t1*vlsm*xnsq**2 + 32.D0
     &    *b*t1*TR*xn*xnsq - 24.D0*b*t1*pi**2*xnsq**2 - 176.D0*b*t1*
     &    rmuom2*xnsq**2 - 96.D0*b*t1*xnsq + 32.D0*b*t1*xnsq**2 - 12.D0
     &    *b*t2**(-2)*ro*vlwm**2*xnsq + 12.D0*b*t2**(-2)*ro*vlwm**2 - 
     &    12.D0*b*t2**(-2)*ro*vdw*xnsq + 12.D0*b*t2**(-2)*ro*vdw - 2.D0
     &    *b*t2**(-2)*ro*pi**2*xnsq
      gcoeff1 = gcoeff1 + 2.D0*b*t2**(-2)*ro*pi**2 - 3.D0*b*t2**(-1)*ro
     &    *vlpm**2 + 3.D0*b*t2**(-1)*ro*pi**2 + 6.D0*b*t2**(-1)*vlpm**2
     &     - 24.D0*b*t2**(-1)*vlsm*vlwm*xnsq**2 - 12.D0*b*t2**(-1)*vlsm
     &    *xnsq**2 - 6.D0*b*t2**(-1)*vlsm**2*xnsq**2 + 96.D0*b*t2**(-1)
     &    *vltm*vlwm*xnsq - 240.D0*b*t2**(-1)*vlwm*xnsq + 174.D0*b*
     &    t2**(-1)*vlwm*xnsq**2 + 48.D0*b*t2**(-1)*vlwm + 24.D0*b*
     &    t2**(-1)*vlwm**2 + 24.D0*b*t2**(-1)*vdw*xnsq**2 + 24.D0*b*
     &    t2**(-1)*vdw - 24.D0*b*t2**(-1)*pi**2*xnsq + 10.D0*b*t2**(-1)
     &    *pi**2*xnsq**2 - 2.D0*b*t2**(-1)*pi**2 - 88.D0*b*t2**(-1)*
     &    rmuom2*xnsq + 88.D0*b*t2**(-1)*rmuom2*xnsq**2 + 48.D0*b*
     &    t2**(-1)*xnsq + 3.D0*b*t2**(-1)*xnsq**2 - 48.D0*b*t2**(-1) - 
     &    12.D0*b*ro*vlpm**2*TR*xn*xnsq + 12.D0*b*ro*TR*pi**2*xn*xnsq
     &     - 96.D0*b*ro*TR*xn*xnsq + 6.D0*b*vlpm**2*xnsq + 24.D0*b*vlsm
     &    *vltm*xnsq**2 + 24.D0*b*vlsm*vlwm*xnsq**2 - 12.D0*b*vlsm*
     &    omro**(-1)*xnsq**2 + 48.D0*b*vlsm*xnsq**2 - 12.D0*b*vlsm**2*
     &    xnsq**2
      gcoeff1 = gcoeff1 - 72.D0*b*vltm*vlwm*xnsq + 24.D0*b*vltm*xnsq**2
     &     + 36.D0*b*vltm**2*xnsq + 96.D0*b*vlwm*xnsq - 168.D0*b*vlwm*
     &    xnsq**2 + 60.D0*b*vlwm**2*xnsq + 24.D0*b*vdw*xnsq - 24.D0*b*
     &    vdw*xnsq**2 + 72.D0*b*vdt*xnsq**2 - 16.D0*b*TR*xn*xnsq + 34.D0
     &    *b*pi**2*xnsq - 4.D0*b*pi**2*xnsq**2 - 16.D0*b*xnsq**2 + 24.D0
     &    *t1**(-1)*ro*vlpm*vlsm*xnsq - 24.D0*t1**(-1)*ro*vlpm*vltm*
     &    xnsq - 24.D0*t1**(-1)*ro*vlpm*vlwm*xnsq + 36.D0*t1**(-1)*ro*
     &    vlpm*xnsq + 12.D0*t1**(-1)*ro*vlpm - 12.D0*t1**(-1)*ro*
     &    vlpm**2*xnsq + 15.D0*t1**(-1)*ro*vlpm**2*xnsq**2 - 48.D0*
     &    t1**(-1)*ro*vdmp*xnsq + 60.D0*t1**(-1)*ro*vdmp*xnsq**2 - 4.D0
     &    *t1**(-1)*ro*pi**2*xnsq + 5.D0*t1**(-1)*ro*pi**2*xnsq**2 - 12.
     &    D0*t1**(-1)*ro**2*vlpm*vlsm*xnsq + 12.D0*t1**(-1)*ro**2*vlpm*
     &    vltm*xnsq + 12.D0*t1**(-1)*ro**2*vlpm*vlwm*xnsq - 12.D0*
     &    t1**(-1)*ro**2*vlpm*xnsq + 6.D0*t1**(-1)*ro**2*vlpm**2*xnsq
     &     + 24.D0*t1**(-1)*ro**2*vdmp*xnsq + 2.D0*t1**(-1)*ro**2*pi**2
     &    *xnsq
      gcoeff1 = gcoeff1 - 12.D0*t1**(-1)*vlpm*vlsm*xnsq + 12.D0*
     &    t1**(-1)*vlpm*vltm*xnsq + 12.D0*t1**(-1)*vlpm*vlwm*xnsq - 24.D
     &    0*t1**(-1)*vlpm*xnsq - 12.D0*t1**(-1)*vlpm + 6.D0*t1**(-1)*
     &    vlpm**2*xnsq - 6.D0*t1**(-1)*vlpm**2*xnsq**2 + 24.D0*t1**(-1)
     &    *vdmp*xnsq - 24.D0*t1**(-1)*vdmp*xnsq**2 + 2.D0*t1**(-1)*
     &    pi**2*xnsq - 2.D0*t1**(-1)*pi**2*xnsq**2 - 96.D0*t1*ro*vlpm*
     &    TR*xn*xnsq - 72.D0*t1*ro*vlpm*xnsq - 12.D0*t1*ro*vlpm**2*xnsq
     &     + 48.D0*t1*ro*vdmb*xnsq - 16.D0*t1*ro*pi**2*xnsq + 96.D0*t1*
     &    ro**2*vlpm*TR*xn*xnsq + 72.D0*t1*vlpm*xnsq - 6.D0*t1*vlpm**2*
     &    omro**(-1)*xnsq**2 + 24.D0*t1*vlpm**2*xnsq - 18.D0*t1*vlpm**2
     &    *xnsq**2 - 24.D0*t1*vdmp*omro**(-1)*xnsq**2 - 72.D0*t1*vdmp*
     &    xnsq**2 - 96.D0*t1*vdmb*xnsq - 2.D0*t1*omro**(-1)*pi**2*
     &    xnsq**2 + 32.D0*t1*pi**2*xnsq - 6.D0*t1*pi**2*xnsq**2 - 12.D0
     &    *t2**(-1)*ro*vlpm + 6.D0*t2**(-1)*ro*vlpm**2*xnsq**2 - 6.D0*
     &    t2**(-1)*ro*vlpm**2 + 24.D0*t2**(-1)*ro*vdmp*xnsq**2 + 24.D0*
     &    t2**(-1)*ro*vdmb
      gcoeff1 = gcoeff1 + 2.D0*t2**(-1)*ro*pi**2*xnsq**2 - 8.D0*
     &    t2**(-1)*ro*pi**2 + 12.D0*t2**(-1)*vlpm*vlsm - 24.D0*t2**(-1)
     &    *vlpm*vlwm + 12.D0*t2**(-1)*vlpm + 6.D0*t2**(-1)*vlpm**2*
     &    xnsq**2 + 6.D0*t2**(-1)*vlpm**2 + 24.D0*t2**(-1)*vdmp*xnsq**2
     &     - 24.D0*t2**(-1)*vdmp - 48.D0*t2**(-1)*vdmb + 2.D0*t2**(-1)*
     &    pi**2*xnsq**2 + 14.D0*t2**(-1)*pi**2 + 48.D0*ro*vlpm*TR*xn*
     &    xnsq + 36.D0*ro*vlpm*xnsq - 48.D0*ro**2*vlpm*TR*xn*xnsq + 12.D
     &    0*vlpm*vlsm*xnsq - 24.D0*vlpm*vlwm*xnsq - 36.D0*vlpm*xnsq + 3.
     &    D0*vlpm**2*omro**(-1)*xnsq**2 - 6.D0*vlpm**2*xnsq - 15.D0*
     &    vlpm**2*xnsq**2 + 12.D0*vdmp*omro**(-1)*xnsq**2 - 24.D0*vdmp*
     &    xnsq - 60.D0*vdmp*xnsq**2 + omro**(-1)*pi**2*xnsq**2 - 2.D0*
     &    pi**2*xnsq - 5.D0*pi**2*xnsq**2

      return
      end
