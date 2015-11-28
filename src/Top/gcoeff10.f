      double precision function gcoeff10()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeff10 =  + TBAR**(-2) * ( 12.D0*b*t1*vltm*xnsq - 6.D0*b*t1*
     &    vltm*xnsq**2 - 6.D0*b*t1*vltm - 6.D0*b*t1**2*vltm*xnsq + 6.D0
     &    *b*t1**2*vltm*xnsq**2 )
      gcoeff10 = gcoeff10 + TBAR**(-1) * ( 36.D0*b*t1*vltm*xnsq - 60.D0
     &    *b*t1*vltm*xnsq**2 + 6.D0*b*t1*xnsq - 6.D0*b*t1*xnsq**2 - 96.D
     &    0*b*vltm*xnsq + 60.D0*b*vltm*xnsq**2 + 36.D0*b*vltm - 12.D0*b
     &    *xnsq + 6.D0*b*xnsq**2 + 6.D0*b )
      gcoeff10 = gcoeff10 + UBAR**(-2) * (  - 6.D0*b*t1*vlwm*xnsq**2 + 
     &    6.D0*b*t1*vlwm - 6.D0*b*t1**2*vlwm*xnsq + 6.D0*b*t1**2*vlwm*
     &    xnsq**2 + 6.D0*b*vlwm*xnsq - 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + UBAR**(-1) * (  - 36.D0*b*t1*vlwm*xnsq + 60.
     &    D0*b*t1*vlwm*xnsq**2 - 6.D0*b*t1*xnsq + 6.D0*b*t1*xnsq**2 - 
     &    60.D0*b*vlwm*xnsq + 36.D0*b*vlwm - 6.D0*b*xnsq + 6.D0*b )
      gcoeff10 = gcoeff10 + INVG * (  - 3.D0/4.D0*b*t1**(-1)*ro*vlpm**2
     &     - 3.D0*b*t1**(-1)*ro*vlsm*vltm*xnsq**2 + 9.D0*b*t1**(-1)*ro*
     &    vlsm*vlwm*xnsq**2 - 3.D0/2.D0*b*t1**(-1)*ro*vlsm**2*xnsq**2
     &     - 3.D0*b*t1**(-1)*ro*vltm**2 - 3.D0/2.D0*b*t1**(-1)*ro*vlwm*
     &    xnsq**2 - 3.D0*b*t1**(-1)*ro*vlwm + 9.D0*b*t1**(-1)*ro*vdw*
     &    xnsq**2 - 3.D0*b*t1**(-1)*ro*vdt*xnsq**2 - 3.D0*b*t1**(-1)*ro
     &    *vdt - 1.D0/2.D0*b*t1**(-1)*ro*pi**2*xnsq**2 + 1.D0/4.D0*b*
     &    t1**(-1)*ro*pi**2 - 3.D0/4.D0*b*t1**(-1)*ro*xnsq**2 - 3.D0/2.D
     &    0*b*t1**(-1)*ro - 3.D0/2.D0*b*t1**(-1)*ro**2*vltm**2*xnsq + 3.
     &    D0*b*t1**(-1)*ro**2*vlwm**2*xnsq + 3.D0*b*t1**(-1)*ro**2*vdw*
     &    xnsq - 3.D0/2.D0*b*t1**(-1)*ro**2*vdt*xnsq + 1.D0/4.D0*b*
     &    t1**(-1)*ro**2*pi**2*xnsq + 3.D0/2.D0*b*t2**(-2)*ro**2*
     &    vlwm**2 + 3.D0/2.D0*b*t2**(-2)*ro**2*vdw + 1.D0/4.D0*b*
     &    t2**(-2)*ro**2*pi**2 - 3.D0/4.D0*b*t2**(-1)*ro*vlpm**2 + 6.D0
     &    *b*t2**(-1)*ro*vlsm*vlwm*xnsq**2 - 3.D0/2.D0*b*t2**(-1)*ro*
     &    vlsm**2*xnsq**2 )
      gcoeff10 = gcoeff10 + INVG * (  - 3.D0/2.D0*b*t2**(-1)*ro*vlwm*
     &    xnsq**2 - 3.D0*b*t2**(-1)*ro*vlwm - 3.D0*b*t2**(-1)*ro*
     &    vlwm**2 + 6.D0*b*t2**(-1)*ro*vdw*xnsq**2 - 3.D0*b*t2**(-1)*ro
     &    *vdw - 1.D0/2.D0*b*t2**(-1)*ro*pi**2*xnsq**2 + 1.D0/4.D0*b*
     &    t2**(-1)*ro*pi**2 - 3.D0/4.D0*b*t2**(-1)*ro*xnsq**2 - 3.D0/2.D
     &    0*b*t2**(-1)*ro + 3.D0/2.D0*b*t2**(-1)*ro**2*vlwm**2*xnsq + 3.
     &    D0/2.D0*b*t2**(-1)*ro**2*vdw*xnsq + 1.D0/4.D0*b*t2**(-1)*
     &    ro**2*pi**2*xnsq + 3.D0*b*ro*vlpm**2 - 6.D0*b*ro*vlsm*vltm*
     &    xnsq**2 - 6.D0*b*ro*vlsm*vlwm*xnsq**2 + 3.D0*b*ro*vlsm**2*
     &    xnsq**2 - 6.D0*b*ro*vltm*vlwm*xnsq + 6.D0*b*ro*vltm**2 - 6.D0
     &    *b*ro*vlwm**2*xnsq + 6.D0*b*ro*vlwm**2 - 9.D0*b*ro*vdw*xnsq
     &     - 6.D0*b*ro*vdw*xnsq**2 + 6.D0*b*ro*vdw - 3.D0*b*ro*vdt*xnsq
     &     - 6.D0*b*ro*vdt*xnsq**2 + 6.D0*b*ro*vdt + b*ro*pi**2*xnsq + 
     &    b*ro*pi**2*xnsq**2 - b*ro*pi**2 - 36.D0*b*vlsm*vlwm*xnsq**2
     &     + 9.D0*b*vlsm**2*xnsq**2 + 24.D0*b*vltm*vlwm*xnsq + 12.D0*b*
     &    vltm**2*xnsq )
      gcoeff10 = gcoeff10 + INVG * ( 6.D0*b*vlwm*xnsq**2 + 12.D0*b*vlwm
     &     + 12.D0*b*vlwm**2*xnsq + 24.D0*b*vdw*xnsq - 36.D0*b*vdw*
     &    xnsq**2 + 24.D0*b*vdt*xnsq - 4.D0*b*pi**2*xnsq + 3.D0*b*pi**2
     &    *xnsq**2 + 3.D0*b*xnsq**2 + 6.D0*b + 3.D0/2.D0*t1**(-1)*ro*
     &    vlpm*vlsm - 3.D0*t1**(-1)*ro*vlpm*vltm - 3.D0*t1**(-1)*ro*
     &    vlpm**2*xnsq**2 - 3.D0/4.D0*t1**(-1)*ro*vlpm**2 - 12.D0*
     &    t1**(-1)*ro*vdmp*xnsq**2 - 3.D0*t1**(-1)*ro*vdmp - t1**(-1)*
     &    ro*pi**2*xnsq**2 - 1.D0/4.D0*t1**(-1)*ro*pi**2 - 3.D0/2.D0*
     &    t1**(-1)*ro**2*vlpm*vlsm + 3.D0*t1**(-1)*ro**2*vlpm*vltm + 3.D
     &    0/4.D0*t1**(-1)*ro**2*vlpm**2 + 3.D0*t1**(-1)*ro**2*vdmp + 1.D
     &    0/4.D0*t1**(-1)*ro**2*pi**2 + 3.D0/2.D0*t2**(-1)*ro*vlpm*vlsm
     &     - 3.D0*t2**(-1)*ro*vlpm*vlwm - 3.D0*t2**(-1)*ro*vlpm**2*
     &    xnsq**2 - 3.D0/4.D0*t2**(-1)*ro*vlpm**2 - 12.D0*t2**(-1)*ro*
     &    vdmp*xnsq**2 - 3.D0*t2**(-1)*ro*vdmp - t2**(-1)*ro*pi**2*
     &    xnsq**2 - 1.D0/4.D0*t2**(-1)*ro*pi**2 - 3.D0/2.D0*t2**(-1)*
     &    ro**2*vlpm*vlsm )
      gcoeff10 = gcoeff10 + INVG * ( 3.D0*t2**(-1)*ro**2*vlpm*vlwm + 3.D
     &    0/4.D0*t2**(-1)*ro**2*vlpm**2 + 3.D0*t2**(-1)*ro**2*vdmp + 1.D
     &    0/4.D0*t2**(-1)*ro**2*pi**2 + 3.D0*ro*vlpm*vlsm*xnsq - 3.D0*
     &    ro*vlpm*vltm*xnsq - 3.D0*ro*vlpm*vlwm*xnsq - 3.D0/2.D0*ro*
     &    vlpm**2*xnsq + 3.D0*ro*vlpm**2*xnsq**2 - 6.D0*ro*vdmp*xnsq + 
     &    12.D0*ro*vdmp*xnsq**2 - 1.D0/2.D0*ro*pi**2*xnsq + ro*pi**2*
     &    xnsq**2 - 3.D0*ro**2*vlpm*vlsm*xnsq + 3.D0*ro**2*vlpm*vltm*
     &    xnsq + 3.D0*ro**2*vlpm*vlwm*xnsq + 3.D0/2.D0*ro**2*vlpm**2*
     &    xnsq + 6.D0*ro**2*vdmp*xnsq + 1.D0/2.D0*ro**2*pi**2*xnsq + 9.D
     &    0*vlpm**2*xnsq**2 + 36.D0*vdmp*xnsq**2 + 3.D0*pi**2*xnsq**2 )
      gcoeff10 = gcoeff10 + INVG**2 * ( 3.D0/8.D0*b*t1**(-1)*ro**3*
     &    vlwm**2*xnsq + 3.D0/8.D0*b*t1**(-1)*ro**3*vdw*xnsq + 1.D0/16.D
     &    0*b*t1**(-1)*ro**3*pi**2*xnsq + 3.D0/8.D0*b*t2**(-2)*ro**3*
     &    vlwm**2 + 3.D0/8.D0*b*t2**(-2)*ro**3*vdw + 1.D0/16.D0*b*
     &    t2**(-2)*ro**3*pi**2 - 3.D0/2.D0*b*t2**(-1)*ro**2*vlwm**2 - 3.
     &    D0/2.D0*b*t2**(-1)*ro**2*vdw - 1.D0/4.D0*b*t2**(-1)*ro**2*
     &    pi**2 + 3.D0/8.D0*b*t2**(-1)*ro**3*vlwm**2*xnsq + 3.D0/8.D0*b
     &    *t2**(-1)*ro**3*vdw*xnsq + 1.D0/16.D0*b*t2**(-1)*ro**3*pi**2*
     &    xnsq - 3.D0/2.D0*b*ro**2*vlwm**2*xnsq + 3.D0/2.D0*b*ro**2*
     &    vlwm**2 - 3.D0/2.D0*b*ro**2*vdw*xnsq + 3.D0/2.D0*b*ro**2*vdw
     &     - 1.D0/4.D0*b*ro**2*pi**2*xnsq + 1.D0/4.D0*b*ro**2*pi**2 )
      gcoeff10 = gcoeff10 + 16.D0*XLF*b*t1**(-1)*TR*rmuom2*xn*xnsq - 16.
     &    D0*XLF*b*t1**(-1)*TR*rmuom2*xn + 16.D0*XLF*b*t2**(-1)*TR*
     &    rmuom2*xn*xnsq - 16.D0*XLF*b*t2**(-1)*TR*rmuom2*xn - 32.D0*
     &    XLF*b*TR*rmuom2*xn*xnsq + 6.D0*b*t1**(-2)*ro*vltm**2*xnsq - 6.
     &    D0*b*t1**(-2)*ro*vltm**2 + 6.D0*b*t1**(-2)*ro*vdt*xnsq - 6.D0
     &    *b*t1**(-2)*ro*vdt + b*t1**(-2)*ro*pi**2*xnsq - b*t1**(-2)*ro
     &    *pi**2 - 12.D0*b*t1**(-1)*ro*vltm**2*xnsq + 6.D0*b*t1**(-1)*
     &    ro*vlwm**2*xnsq + 6.D0*b*t1**(-1)*ro*vdw*xnsq - 12.D0*b*
     &    t1**(-1)*ro*vdt*xnsq - b*t1**(-1)*ro*pi**2*xnsq - 3.D0*b*
     &    t1**(-1)*vlpm**2 + 12.D0*b*t1**(-1)*vlsm*vltm*xnsq**2 + 36.D0
     &    *b*t1**(-1)*vlsm*vlwm*xnsq**2 - 6.D0*b*t1**(-1)*vlsm**2*
     &    xnsq**2 - 48.D0*b*t1**(-1)*vltm*vlwm*xnsq + 120.D0*b*t1**(-1)
     &    *vltm*xnsq - 72.D0*b*t1**(-1)*vltm*xnsq**2 - 48.D0*b*t1**(-1)
     &    *vltm - 12.D0*b*t1**(-1)*vltm**2 - 6.D0*b*t1**(-1)*vlwm*
     &    xnsq**2 - 12.D0*b*t1**(-1)*vlwm + 36.D0*b*t1**(-1)*vdw*
     &    xnsq**2
      gcoeff10 = gcoeff10 - 12.D0*b*t1**(-1)*vdt*xnsq**2 - 12.D0*b*
     &    t1**(-1)*vdt + 12.D0*b*t1**(-1)*pi**2*xnsq - 8.D0*b*t1**(-1)*
     &    pi**2*xnsq**2 + b*t1**(-1)*pi**2 + 44.D0*b*t1**(-1)*rmuom2*
     &    xnsq - 44.D0*b*t1**(-1)*rmuom2*xnsq**2 - 24.D0*b*t1**(-1)*
     &    xnsq - 3.D0*b*t1**(-1)*xnsq**2 + 18.D0*b*t1**(-1) + 6.D0*b*
     &    t2**(-2)*ro*vlwm**2*xnsq - 6.D0*b*t2**(-2)*ro*vlwm**2 + 6.D0*
     &    b*t2**(-2)*ro*vdw*xnsq - 6.D0*b*t2**(-2)*ro*vdw + b*t2**(-2)*
     &    ro*pi**2*xnsq - b*t2**(-2)*ro*pi**2 - 6.D0*b*t2**(-1)*ro*
     &    vlwm**2*xnsq - 6.D0*b*t2**(-1)*ro*vdw*xnsq - b*t2**(-1)*ro*
     &    pi**2*xnsq - 3.D0*b*t2**(-1)*vlpm**2 + 48.D0*b*t2**(-1)*vlsm*
     &    vlwm*xnsq**2 - 6.D0*b*t2**(-1)*vlsm**2*xnsq**2 - 48.D0*b*
     &    t2**(-1)*vltm*vlwm*xnsq + 120.D0*b*t2**(-1)*vlwm*xnsq - 78.D0
     &    *b*t2**(-1)*vlwm*xnsq**2 - 60.D0*b*t2**(-1)*vlwm - 12.D0*b*
     &    t2**(-1)*vlwm**2 + 24.D0*b*t2**(-1)*vdw*xnsq**2 - 12.D0*b*
     &    t2**(-1)*vdw + 12.D0*b*t2**(-1)*pi**2*xnsq - 8.D0*b*t2**(-1)*
     &    pi**2*xnsq**2
      gcoeff10 = gcoeff10 + b*t2**(-1)*pi**2 + 44.D0*b*t2**(-1)*rmuom2*
     &    xnsq - 44.D0*b*t2**(-1)*rmuom2*xnsq**2 - 24.D0*b*t2**(-1)*
     &    xnsq - 3.D0*b*t2**(-1)*xnsq**2 + 18.D0*b*t2**(-1) - 6.D0*b*
     &    vlpm**2*xnsq - 48.D0*b*vlsm*vltm*xnsq**2 - 48.D0*b*vlsm*vlwm*
     &    xnsq**2 + 12.D0*b*vlsm**2*xnsq**2 + 72.D0*b*vltm*vlwm*xnsq - 
     &    48.D0*b*vltm*xnsq + 72.D0*b*vltm*xnsq**2 - 48.D0*b*vltm**2*
     &    xnsq - 48.D0*b*vlwm*xnsq + 72.D0*b*vlwm*xnsq**2 - 48.D0*b*
     &    vlwm**2*xnsq - 12.D0*b*vdw*xnsq - 24.D0*b*vdw*xnsq**2 - 12.D0
     &    *b*vdt*xnsq - 24.D0*b*vdt*xnsq**2 - 34.D0*b*pi**2*xnsq + 16.D0
     &    *b*pi**2*xnsq**2 + 88.D0*b*rmuom2*xnsq**2 + 48.D0*b*xnsq + 3.D
     &    0*t1**(-1)*ro*vlpm**2 - 12.D0*t1**(-1)*ro*vdmb + 4.D0*
     &    t1**(-1)*ro*pi**2 - 6.D0*t1**(-1)*vlpm*vlsm + 12.D0*t1**(-1)*
     &    vlpm*vltm - 12.D0*t1**(-1)*vlpm**2*xnsq**2 - 3.D0*t1**(-1)*
     &    vlpm**2 - 48.D0*t1**(-1)*vdmp*xnsq**2 + 12.D0*t1**(-1)*vdmp
     &     + 24.D0*t1**(-1)*vdmb - 4.D0*t1**(-1)*pi**2*xnsq**2 - 7.D0*
     &    t1**(-1)*pi**2
      gcoeff10 = gcoeff10 + 3.D0*t2**(-1)*ro*vlpm**2 - 12.D0*t2**(-1)*
     &    ro*vdmb + 4.D0*t2**(-1)*ro*pi**2 - 6.D0*t2**(-1)*vlpm*vlsm + 
     &    12.D0*t2**(-1)*vlpm*vlwm - 12.D0*t2**(-1)*vlpm**2*xnsq**2 - 3.
     &    D0*t2**(-1)*vlpm**2 - 48.D0*t2**(-1)*vdmp*xnsq**2 + 12.D0*
     &    t2**(-1)*vdmp + 24.D0*t2**(-1)*vdmb - 4.D0*t2**(-1)*pi**2*
     &    xnsq**2 - 7.D0*t2**(-1)*pi**2 + 6.D0*ro*vlpm**2*xnsq - 24.D0*
     &    ro*vdmb*xnsq + 8.D0*ro*pi**2*xnsq - 12.D0*vlpm*vlsm*xnsq + 12.
     &    D0*vlpm*vltm*xnsq + 12.D0*vlpm*vlwm*xnsq - 6.D0*vlpm**2*xnsq
     &     + 24.D0*vlpm**2*xnsq**2 + 24.D0*vdmp*xnsq + 96.D0*vdmp*
     &    xnsq**2 + 48.D0*vdmb*xnsq - 14.D0*pi**2*xnsq + 8.D0*pi**2*
     &    xnsq**2

      return
      end
