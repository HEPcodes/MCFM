      subroutine fill_gcoeff(
     . gcoeff1,gcoeff2,gcoeff3,gcoeff4,gcoeff5,
     . gcoeff6,gcoeff7,gcoeff8,gcoeff9,gcoeff10)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      double precision
     . gcoeff1,gcoeff2,gcoeff3,gcoeff4,gcoeff5,
     . gcoeff6,gcoeff7,gcoeff8,gcoeff9,gcoeff10

      gcoeff1 =  + t1**(-1)*ro*xnsq * ( 12.D0*b*vltm**2 + 12.D0*b*vdt
     &     + 2.D0*b*pi**2 )
      gcoeff1 = gcoeff1 + t1**(-1)*ro * (  - 12.D0*b*vltm**2 - 12.D0*b*
     &    vdt - 2.D0*b*pi**2 )
      gcoeff1 = gcoeff1 + t1**(-1)*UBAR**(-1)*INVG*xnsq**2 * ( 12.D0*b*
     &    vlwm )
      gcoeff1 = gcoeff1 + t1**(-1)*UBAR**(-1)*INVG * (  - 12.D0*b*vlwm
     &     )
      gcoeff1 = gcoeff1 + t1**(-1)*UBAR**(-1)*xnsq**2 * ( 12.D0*b*vlwm
     &     )
      gcoeff1 = gcoeff1 + t1**(-1)*UBAR**(-1) * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t1**(-1)*INVG*xnsq**2 * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t1**(-1)*INVG * ( 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t1*ro*xnsq * ( 24.D0*b*vlpm**2*TR*xn - 24.D0*
     &    b*TR*pi**2*xn + 192.D0*b*TR*xn - 96.D0*vlpm*TR*xn - 120.D0*
     &    vlpm - 12.D0*vlpm**2 + 48.D0*vdmb - 16.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1*ro**2*xnsq * ( 96.D0*vlpm*TR*xn )
      gcoeff1 = gcoeff1 + t1*UBAR**(-2)*xnsq**2 * ( 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t1*UBAR**(-2) * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t1*UBAR**(-1)*xnsq * ( 84.D0*b*vlwm + 12.D0*b
     &     )
      gcoeff1 = gcoeff1 + t1*UBAR**(-1)*xnsq**2 * (  - 132.D0*b*vlwm - 
     &    12.D0*b )
      gcoeff1 = gcoeff1 + t1*INVG*xnsq * (  - 96.D0*b*vlpm**2 - 84.D0*b
     &    *vltm - 36.D0*b*vlwm - 384.D0*b*vlwm**2 - 384.D0*b*vdw + 32.D0
     &    *b*pi**2 - 348.D0*vlpm*vlsm + 696.D0*vlpm*vlwm + 96.D0*vlpm
     &     + 174.D0*vlpm**2 + 696.D0*vdmp + 58.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1*INVG*xnsq**2 * ( 24.D0*b*vlsm*vlwm - 6.D0*
     &    b*vlsm**2 + 36.D0*b*vltm - 12.D0*b*vlwm + 24.D0*b*vdw - 2.D0*
     &    b*pi**2 - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1*INVG * (  - 78.D0*b*vlpm**2 + 48.D0*b*vltm
     &     - 48.D0*b*vlwm - 312.D0*b*vlwm**2 - 312.D0*b*vdw + 26.D0*b*
     &    pi**2 - 300.D0*vlpm*vlsm + 600.D0*vlpm*vlwm + 72.D0*vlpm + 
     &    150.D0*vlpm**2 + 600.D0*vdmp + 50.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1*INVG**2*xnsq * (  - 27.D0*b*vlpm**2 - 108.D
     &    0*b*vlwm**2 - 108.D0*b*vdw + 9.D0*b*pi**2 - 66.D0*vlpm*vlsm
     &     + 132.D0*vlpm*vlwm + 33.D0*vlpm**2 + 132.D0*vdmp + 11.D0*
     &    pi**2 )
      gcoeff1 = gcoeff1 + t1*INVG**2*xnsq**2 * ( 12.D0*b*vlsm*vlwm - 3.D
     &    0*b*vlsm**2 + 12.D0*b*vdw - b*pi**2 - 3.D0*vlpm**2 - 12.D0*
     &    vdmp - pi**2 )
      gcoeff1 = gcoeff1 + t1*INVG**2 * (  - 24.D0*b*vlpm**2 - 96.D0*b*
     &    vlwm**2 - 96.D0*b*vdw + 8.D0*b*pi**2 - 60.D0*vlpm*vlsm + 120.D
     &    0*vlpm*vlwm + 30.D0*vlpm**2 + 120.D0*vdmp + 10.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1*xnsq * ( 64.D0*XLF*b*TR*rmuom2*xn + 32.D0*
     &    XLF*b*TR*xn - 48.D0*b*vlpm**2 - 48.D0*b*vltm - 24.D0*b*
     &    vltm**2 - 48.D0*b*vlwm - 168.D0*b*vlwm**2 - 168.D0*b*vdw - 24.
     &    D0*b*vdt + 32.D0*b*TR*xn + 16.D0*b*pi**2 - 96.D0*b - 336.D0*
     &    vlpm*vlsm + 672.D0*vlpm*vlwm + 312.D0*vlpm + 192.D0*vlpm**2
     &     + 672.D0*vdmp - 96.D0*vdmb + 88.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1*xnsq**2 * ( 24.D0*b*vlsm*omro**(-1) - 96.D0
     &    *b*vlsm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.D0*b*vdt - 16.D0
     &    *b*pi**2 - 176.D0*b*rmuom2 + 32.D0*b - 6.D0*vlpm**2*
     &    omro**(-1) - 6.D0*vlpm**2 - 24.D0*vdmp*omro**(-1) - 24.D0*
     &    vdmp - 2.D0*omro**(-1)*pi**2 - 2.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1 * (  - 24.D0*b*vlpm**2 - 48.D0*b*vltm**2
     &     - 48.D0*b*vlwm**2 - 48.D0*b*vdw - 48.D0*b*vdt + 8.D0*b*pi**2
     &     - 240.D0*vlpm*vlsm + 480.D0*vlpm*vlwm + 96.D0*vlpm + 120.D0*
     &    vlpm**2 + 480.D0*vdmp + 40.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2*UBAR**(-2)*xnsq * ( 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t1**2*UBAR**(-2)*xnsq**2 * (  - 12.D0*b*vlwm
     &     )
      gcoeff1 = gcoeff1 + t1**2*INVG*xnsq * ( 180.D0*b*vlpm**2 - 96.D0*
     &    b*vltm*vlwm + 96.D0*b*vltm - 48.D0*b*vltm**2 + 48.D0*b*vlwm
     &     + 672.D0*b*vlwm**2 + 624.D0*b*vdw - 96.D0*b*vdt - 44.D0*b*
     &    pi**2 + 972.D0*vlpm*vlsm - 1944.D0*vlpm*vlwm - 312.D0*vlpm - 
     &    486.D0*vlpm**2 - 1944.D0*vdmp - 162.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2*INVG*xnsq**2 * ( 48.D0*b*vlsm*vltm + 48.
     &    D0*b*vlsm*vlwm - 24.D0*b*vlsm**2 - 24.D0*b*vltm + 24.D0*b*
     &    vlwm + 48.D0*b*vdw + 48.D0*b*vdt - 8.D0*b*pi**2 + 30.D0*
     &    vlpm**2 + 120.D0*vdmp + 10.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2*INVG * ( 120.D0*b*vlpm**2 - 48.D0*b*
     &    vltm + 72.D0*b*vltm**2 + 48.D0*b*vlwm + 408.D0*b*vlwm**2 + 
     &    408.D0*b*vdw + 72.D0*b*vdt - 40.D0*b*pi**2 + 624.D0*vlpm*vlsm
     &     + 96.D0*vlpm*vltm - 1344.D0*vlpm*vlwm - 144.D0*vlpm - 312.D0
     &    *vlpm**2 - 1248.D0*vdmp - 104.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2*INVG**2*xnsq * ( 99.D0*b*vlpm**2 - 24.D0
     &    *b*vltm*vlwm - 12.D0*b*vltm**2 + 384.D0*b*vlwm**2 + 372.D0*b*
     &    vdw - 24.D0*b*vdt - 29.D0*b*pi**2 + 306.D0*vlpm*vlsm - 612.D0
     &    *vlpm*vlwm - 153.D0*vlpm**2 - 612.D0*vdmp - 51.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2*INVG**2*xnsq**2 * ( 12.D0*b*vlsm*vltm
     &     - 48.D0*b*vlsm*vlwm + 9.D0*b*vlsm**2 - 48.D0*b*vdw + 12.D0*b
     &    *vdt + 3.D0*b*pi**2 + 21.D0*vlpm**2 + 84.D0*vdmp + 7.D0*pi**2
     &     )
      gcoeff1 = gcoeff1 + t1**2*INVG**2 * ( 78.D0*b*vlpm**2 + 12.D0*b*
     &    vltm**2 + 300.D0*b*vlwm**2 + 300.D0*b*vdw + 12.D0*b*vdt - 26.D
     &    0*b*pi**2 + 240.D0*vlpm*vlsm + 12.D0*vlpm*vltm - 492.D0*vlpm*
     &    vlwm - 120.D0*vlpm**2 - 480.D0*vdmp - 40.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2*xnsq * (  - 48.D0*b*vltm**2 + 48.D0*b*
     &    vlwm**2 + 48.D0*b*vdw - 48.D0*b*vdt + 336.D0*vlpm*vlsm - 672.D
     &    0*vlpm*vlwm - 288.D0*vlpm - 168.D0*vlpm**2 - 672.D0*vdmp - 56.
     &    D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**2 * ( 192.D0*vlpm*vltm - 192.D0*vlpm*vlwm
     &     )
      gcoeff1 = gcoeff1 + t1**3*INVG*xnsq * (  - 120.D0*b*vlpm**2 + 48.D
     &    0*b*vltm*vlwm - 48.D0*b*vltm + 96.D0*b*vltm**2 - 48.D0*b*vlwm
     &     - 528.D0*b*vlwm**2 - 504.D0*b*vdw + 120.D0*b*vdt + 32.D0*b*
     &    pi**2 - 1344.D0*vlpm*vlsm + 96.D0*vlpm*vltm + 2592.D0*vlpm*
     &    vlwm + 528.D0*vlpm + 672.D0*vlpm**2 + 2688.D0*vdmp + 224.D0*
     &    pi**2 )
      gcoeff1 = gcoeff1 + t1**3*INVG*xnsq**2 * (  - 48.D0*b*vlsm*vltm
     &     - 48.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.D0*
     &    b*vdt + 8.D0*b*pi**2 - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*
     &    pi**2 )
      gcoeff1 = gcoeff1 + t1**3*INVG * (  - 72.D0*b*vlpm**2 - 144.D0*b*
     &    vltm**2 - 144.D0*b*vlwm**2 - 144.D0*b*vdw - 144.D0*b*vdt + 24.
     &    D0*b*pi**2 - 432.D0*vlpm*vlsm - 384.D0*vlpm*vltm + 1248.D0*
     &    vlpm*vlwm + 96.D0*vlpm + 216.D0*vlpm**2 + 864.D0*vdmp + 72.D0
     &    *pi**2 )
      gcoeff1 = gcoeff1 + t1**3*INVG**2*xnsq * (  - 186.D0*b*vlpm**2 + 
     &    96.D0*b*vltm*vlwm + 60.D0*b*vltm**2 - 708.D0*b*vlwm**2 - 660.D
     &    0*b*vdw + 108.D0*b*vdt + 46.D0*b*pi**2 - 780.D0*vlpm*vlsm + 
     &    12.D0*vlpm*vltm + 1548.D0*vlpm*vlwm + 390.D0*vlpm**2 + 1560.D0
     &    *vdmp + 130.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**3*INVG**2*xnsq**2 * (  - 36.D0*b*vlsm*
     &    vltm + 60.D0*b*vlsm*vlwm - 6.D0*b*vlsm**2 + 60.D0*b*vdw - 36.D
     &    0*b*vdt - 2.D0*b*pi**2 - 54.D0*vlpm**2 - 216.D0*vdmp - 18.D0*
     &    pi**2 )
      gcoeff1 = gcoeff1 + t1**3*INVG**2 * (  - 132.D0*b*vlpm**2 - 72.D0
     &    *b*vltm**2 - 456.D0*b*vlwm**2 - 456.D0*b*vdw - 72.D0*b*vdt + 
     &    44.D0*b*pi**2 - 480.D0*vlpm*vlsm - 96.D0*vlpm*vltm + 1056.D0*
     &    vlpm*vlwm + 240.D0*vlpm**2 + 960.D0*vdmp + 80.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**3*xnsq * (  - 192.D0*vlpm*vlsm + 192.D0*
     &    vlpm*vltm + 192.D0*vlpm*vlwm + 192.D0*vlpm + 96.D0*vlpm**2 + 
     &    384.D0*vdmp + 32.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**4*INVG*xnsq * (  - 144.D0*b*vltm**2 + 144.
     &    D0*b*vlwm**2 + 144.D0*b*vdw - 144.D0*b*vdt + 1008.D0*vlpm*
     &    vlsm - 384.D0*vlpm*vltm - 1632.D0*vlpm*vlwm - 480.D0*vlpm - 
     &    504.D0*vlpm**2 - 2016.D0*vdmp - 168.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**4*INVG * ( 384.D0*vlpm*vltm - 384.D0*vlpm
     &    *vlwm )
      gcoeff1 = gcoeff1 + t1**4*INVG**2*xnsq * ( 180.D0*b*vlpm**2 - 120.
     &    D0*b*vltm*vlwm - 132.D0*b*vltm**2 + 732.D0*b*vlwm**2 + 672.D0
     &    *b*vdw - 192.D0*b*vdt - 40.D0*b*pi**2 + 1200.D0*vlpm*vlsm - 
     &    96.D0*vlpm*vltm - 2304.D0*vlpm*vlwm - 600.D0*vlpm**2 - 2400.D0
     &    *vdmp - 200.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**4*INVG**2*xnsq**2 * ( 24.D0*b*vlsm*vltm
     &     - 24.D0*b*vlsm*vlwm - 24.D0*b*vdw + 24.D0*b*vdt + 60.D0*
     &    vlpm**2 + 240.D0*vdmp + 20.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**4*INVG**2 * ( 120.D0*b*vlpm**2 + 144.D0*b
     &    *vltm**2 + 336.D0*b*vlwm**2 + 336.D0*b*vdw + 144.D0*b*vdt - 
     &    40.D0*b*pi**2 + 480.D0*vlpm*vlsm + 288.D0*vlpm*vltm - 1248.D0
     &    *vlpm*vlwm - 240.D0*vlpm**2 - 960.D0*vdmp - 80.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**5*INVG*xnsq * (  - 384.D0*vlpm*vlsm + 384.
     &    D0*vlpm*vltm + 384.D0*vlpm*vlwm + 192.D0*vlpm + 192.D0*
     &    vlpm**2 + 768.D0*vdmp + 64.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**5*INVG**2*xnsq * (  - 72.D0*b*vlpm**2 + 
     &    48.D0*b*vltm*vlwm + 168.D0*b*vltm**2 - 408.D0*b*vlwm**2 - 384.
     &    D0*b*vdw + 192.D0*b*vdt + 16.D0*b*pi**2 - 1152.D0*vlpm*vlsm
     &     + 288.D0*vlpm*vltm + 2016.D0*vlpm*vlwm + 576.D0*vlpm**2 + 
     &    2304.D0*vdmp + 192.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**5*INVG**2*xnsq**2 * (  - 24.D0*vlpm**2 - 
     &    96.D0*vdmp - 8.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**5*INVG**2 * (  - 48.D0*b*vlpm**2 - 96.D0*
     &    b*vltm**2 - 96.D0*b*vlwm**2 - 96.D0*b*vdw - 96.D0*b*vdt + 16.D
     &    0*b*pi**2 - 192.D0*vlpm*vlsm - 384.D0*vlpm*vltm + 768.D0*vlpm
     &    *vlwm + 96.D0*vlpm**2 + 384.D0*vdmp + 32.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**6*INVG**2*xnsq * (  - 96.D0*b*vltm**2 + 
     &    96.D0*b*vlwm**2 + 96.D0*b*vdw - 96.D0*b*vdt + 672.D0*vlpm*
     &    vlsm - 384.D0*vlpm*vltm - 960.D0*vlpm*vlwm - 336.D0*vlpm**2
     &     - 1344.D0*vdmp - 112.D0*pi**2 )
      gcoeff1 = gcoeff1 + t1**6*INVG**2 * ( 192.D0*vlpm*vltm - 192.D0*
     &    vlpm*vlwm )
      gcoeff1 = gcoeff1 + t1**7*INVG**2*xnsq * (  - 192.D0*vlpm*vlsm + 
     &    192.D0*vlpm*vltm + 192.D0*vlpm*vlwm + 96.D0*vlpm**2 + 384.D0*
     &    vdmp + 32.D0*pi**2 )
      gcoeff1 = gcoeff1 + t2**(-2)*ro*xnsq * (  - 12.D0*b*vlwm**2 - 12.D
     &    0*b*vdw - 2.D0*b*pi**2 )
      gcoeff1 = gcoeff1 + t2**(-2)*ro * ( 12.D0*b*vlwm**2 + 12.D0*b*vdw
     &     + 2.D0*b*pi**2 )
      gcoeff1 = gcoeff1 + t2**(-1)*ro * ( 12.D0*b*vlwm**2 + 12.D0*b*vdw
     &     + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 24.D0*vlpm*vlwm + 24.D0*
     &    vdmp + 24.D0*vdmb - 6.D0*pi**2 )
      gcoeff1 = gcoeff1 + t2**(-1)*UBAR**(-1)*xnsq**2 * ( 12.D0*b*vlwm
     &     )
      gcoeff1 = gcoeff1 + t2**(-1)*UBAR**(-1) * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t2**(-1)*INVG*xnsq**2 * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t2**(-1)*INVG * ( 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + t2**(-1)*xnsq * (  - 32.D0*XLF*b*TR*rmuom2*xn
     &     + 96.D0*b*vltm*vlwm - 240.D0*b*vlwm - 24.D0*b*pi**2 - 88.D0*
     &    b*rmuom2 + 48.D0*b )
      gcoeff1 = gcoeff1 + t2**(-1)*xnsq**2 * (  - 48.D0*b*vlsm*vlwm + 
     &    144.D0*b*vlwm + 12.D0*b*pi**2 + 88.D0*b*rmuom2 )
      gcoeff1 = gcoeff1 + t2**(-1) * ( 32.D0*XLF*b*TR*rmuom2*xn + 96.D0
     &    *b*vlwm - 48.D0*b + 24.D0*vlpm*vlsm - 48.D0*vlpm*vlwm - 48.D0
     &    *vdmp - 48.D0*vdmb + 12.D0*pi**2 )
      gcoeff1 = gcoeff1 + ro*TBAR**(-1)*xnsq * ( 3.D0*b*vltm )
      gcoeff1 = gcoeff1 + ro*TBAR**(-1)*xnsq**2 * (  - 3.D0*b*vltm )
      gcoeff1 = gcoeff1 + ro*xnsq * (  - 12.D0*b*vlpm**2*TR*xn - 12.D0*
     &    b*vltm**2 + 12.D0*b*vlwm**2 + 12.D0*b*vdw - 12.D0*b*vdt + 12.D
     &    0*b*TR*pi**2*xn - 96.D0*b*TR*xn - 12.D0*vlpm*vlsm + 24.D0*
     &    vlpm*vlwm + 48.D0*vlpm*TR*xn + 60.D0*vlpm + 6.D0*vlpm**2 + 24.
     &    D0*vdmp + 2.D0*pi**2 )
      gcoeff1 = gcoeff1 + ro**2*xnsq * (  - 48.D0*vlpm*TR*xn )
      gcoeff1 = gcoeff1 + TBAR**(-1)*xnsq * (  - 24.D0*b*vltm )
      gcoeff1 = gcoeff1 + TBAR**(-1)*xnsq**2 * ( 12.D0*b*vltm )
      gcoeff1 = gcoeff1 + TBAR**(-1) * ( 12.D0*b*vltm )
      gcoeff1 = gcoeff1 + UBAR**(-2)*xnsq * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + UBAR**(-2) * ( 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + UBAR**(-1)*INVG*xnsq**2 * (  - 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + UBAR**(-1)*INVG * ( 12.D0*b*vlwm )
      gcoeff1 = gcoeff1 + UBAR**(-1)*xnsq * ( 132.D0*b*vlwm + 12.D0*b )
      gcoeff1 = gcoeff1 + UBAR**(-1) * (  - 84.D0*b*vlwm - 12.D0*b )
      gcoeff1 = gcoeff1 + INVG*xnsq * ( 18.D0*b*vlpm**2 + 24.D0*b*vltm
     &     + 12.D0*b*vlwm + 72.D0*b*vlwm**2 + 72.D0*b*vdw - 6.D0*b*
     &    pi**2 + 48.D0*vlpm*vlsm - 96.D0*vlpm*vlwm - 12.D0*vlpm - 24.D0
     &    *vlpm**2 - 96.D0*vdmp - 8.D0*pi**2 )
      gcoeff1 = gcoeff1 + INVG*xnsq**2 * (  - 12.D0*b*vltm )
      gcoeff1 = gcoeff1 + INVG * ( 18.D0*b*vlpm**2 - 12.D0*b*vltm + 12.D
     &    0*b*vlwm + 72.D0*b*vlwm**2 + 72.D0*b*vdw - 6.D0*b*pi**2 + 48.D
     &    0*vlpm*vlsm - 96.D0*vlpm*vlwm - 12.D0*vlpm - 24.D0*vlpm**2 - 
     &    96.D0*vdmp - 8.D0*pi**2 )
      gcoeff1 = gcoeff1 + INVG**2*xnsq * ( 3.D0*b*vlpm**2 + 12.D0*b*
     &    vlwm**2 + 12.D0*b*vdw - b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vlwm - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff1 = gcoeff1 + INVG**2 * ( 3.D0*b*vlpm**2 + 12.D0*b*vlwm**2
     &     + 12.D0*b*vdw - b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm*vlwm
     &     - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff1 = gcoeff1 + xnsq * (  - 16.D0*XLF*b*TR*xn + 30.D0*b*
     &    vlpm**2 - 96.D0*b*vltm*vlwm + 60.D0*b*vltm + 24.D0*b*vltm**2
     &     + 96.D0*b*vlwm + 144.D0*b*vlwm**2 + 96.D0*b*vdw - 24.D0*b*
     &    vdt - 16.D0*b*TR*xn + 30.D0*b*pi**2 + 120.D0*vlpm*vlsm - 240.D
     &    0*vlpm*vlwm - 108.D0*vlpm - 60.D0*vlpm**2 - 240.D0*vdmp - 20.D
     &    0*pi**2 )
      gcoeff1 = gcoeff1 + xnsq**2 * ( 48.D0*b*vlsm*vlwm - 12.D0*b*vlsm*
     &    omro**(-1) + 48.D0*b*vlsm - 12.D0*b*vlsm**2 - 12.D0*b*vltm - 
     &    144.D0*b*vlwm + 48.D0*b*vdt - 4.D0*b*pi**2 - 16.D0*b + 3.D0*
     &    vlpm**2*omro**(-1) - 9.D0*vlpm**2 + 12.D0*vdmp*omro**(-1) - 
     &    36.D0*vdmp + omro**(-1)*pi**2 - 3.D0*pi**2 )
      gcoeff1 = gcoeff1 + 24.D0*b*vlpm**2 - 48.D0*b*vltm - 24.D0*b*
     &    vltm**2 + 48.D0*b*vlwm + 120.D0*b*vlwm**2 + 120.D0*b*vdw - 24.
     &    D0*b*vdt - 8.D0*b*pi**2 + 96.D0*vlpm*vlsm - 192.D0*vlpm*vlwm
     &     - 48.D0*vlpm - 48.D0*vlpm**2 - 192.D0*vdmp - 16.D0*pi**2

      gcoeff2 =  + t1**(-2)*ro*xnsq * (  - 6.D0*b*vltm**2 - 6.D0*b*vdt
     &     - b*pi**2 )
      gcoeff2 = gcoeff2 + t1**(-2)*ro * ( 6.D0*b*vltm**2 + 6.D0*b*vdt
     &     + b*pi**2 )
      gcoeff2 = gcoeff2 + t1**(-1)*ro*xnsq * (  - 6.D0*b*vltm**2 - 6.D0
     &    *b*vdt - b*pi**2 )
      gcoeff2 = gcoeff2 + t1**(-1)*ro * ( 12.D0*b*vltm**2 + 12.D0*b*vdt
     &     + 2.D0*b*pi**2 )
      gcoeff2 = gcoeff2 + t1**(-1)*TBAR**(-1)*xnsq * ( 24.D0*b*vltm )
      gcoeff2 = gcoeff2 + t1**(-1)*TBAR**(-1)*xnsq**2 * (  - 12.D0*b*
     &    vltm )
      gcoeff2 = gcoeff2 + t1**(-1)*TBAR**(-1) * (  - 12.D0*b*vltm )
      gcoeff2 = gcoeff2 + t1**(-1)*UBAR**(-1)*INVG*xnsq * ( 24.D0*b*
     &    vlwm )
      gcoeff2 = gcoeff2 + t1**(-1)*UBAR**(-1)*INVG*xnsq**2 * (  - 12.D0
     &    *b*vlwm )
      gcoeff2 = gcoeff2 + t1**(-1)*UBAR**(-1)*INVG * (  - 12.D0*b*vlwm
     &     )
      gcoeff2 = gcoeff2 + t1**(-1)*UBAR**(-1)*xnsq * ( 24.D0*b*vlwm )
      gcoeff2 = gcoeff2 + t1**(-1)*UBAR**(-1)*xnsq**2 * (  - 12.D0*b*
     &    vlwm )
      gcoeff2 = gcoeff2 + t1**(-1)*UBAR**(-1) * (  - 12.D0*b*vlwm )
      gcoeff2 = gcoeff2 + t1**(-1)*INVG*xnsq * (  - 24.D0*b*vltm - 24.D0
     &    *b*vlwm )
      gcoeff2 = gcoeff2 + t1**(-1)*INVG*xnsq**2 * ( 3.D0/2.D0*b*vlsm*
     &    omro**(-1) - 3.D0/2.D0*b*vlsm + 12.D0*b*vltm + 12.D0*b*vlwm
     &     - 3.D0/8.D0*vlpm**2*omro**(-1) + 3.D0/8.D0*vlpm**2 - 3.D0/2.D
     &    0*vdmp*omro**(-1) + 3.D0/2.D0*vdmp - 1.D0/8.D0*omro**(-1)*
     &    pi**2 + 1.D0/8.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**(-1)*INVG * ( 12.D0*b*vltm + 12.D0*b*vlwm
     &     )
      gcoeff2 = gcoeff2 + t1**(-1)*xnsq * (  - 24.D0*b*vltm + 24.D0*b )
      gcoeff2 = gcoeff2 + t1**(-1)*xnsq**2 * ( 6.D0*b*vlsm*omro**(-1)
     &     - 3.D0/2.D0*vlpm**2*omro**(-1) - 6.D0*vdmp*omro**(-1) - 1.D0/
     &    2.D0*omro**(-1)*pi**2 )
      gcoeff2 = gcoeff2 + t1**(-1) * ( 24.D0*b*vltm - 24.D0*b )
      gcoeff2 = gcoeff2 + t1*ro*xnsq * (  - 24.D0*b*vlpm**2*TR*xn + 24.D
     &    0*b*TR*pi**2*xn - 192.D0*b*TR*xn + 96.D0*vlpm*TR*xn + 72.D0*
     &    vlpm )
      gcoeff2 = gcoeff2 + t1*ro**2*xnsq * (  - 96.D0*vlpm*TR*xn )
      gcoeff2 = gcoeff2 + t1*INVG*xnsq * ( 186.D0*b*vlpm**2 - 72.D0*b*
     &    vltm*vlwm + 156.D0*b*vltm - 36.D0*b*vltm**2 + 60.D0*b*vlwm + 
     &    708.D0*b*vlwm**2 + 672.D0*b*vdw - 72.D0*b*vdt - 50.D0*b*pi**2
     &     + 72.D0*b + 756.D0*vlpm*vlsm - 1512.D0*vlpm*vlwm - 324.D0*
     &    vlpm - 378.D0*vlpm**2 - 1512.D0*vdmp - 126.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG*xnsq**2 * ( 24.D0*b*vlsm*vltm - 24.D0
     &    *b*vlsm*vlwm + 6.D0*b*vlsm - 60.D0*b*vltm + 36.D0*b*vlwm - 24.
     &    D0*b*vdw + 24.D0*b*vdt - 24.D0*b + 51.D0/2.D0*vlpm**2 + 102.D0
     &    *vdmp + 17.D0/2.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG * ( 168.D0*b*vlpm**2 - 120.D0*b*vltm
     &     + 48.D0*b*vltm**2 + 168.D0*b*vlwm + 624.D0*b*vlwm**2 + 624.D0
     &    *b*vdw + 48.D0*b*vdt - 56.D0*b*pi**2 - 24.D0*b + 648.D0*vlpm*
     &    vlsm + 48.D0*vlpm*vltm - 1344.D0*vlpm*vlwm - 240.D0*vlpm - 
     &    324.D0*vlpm**2 - 1296.D0*vdmp - 108.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG**2*xnsq * ( 84.D0*b*vlpm**2 + 48.D0*b
     &    *vltm + 24.D0*b*vlwm + 336.D0*b*vlwm**2 + 336.D0*b*vdw - 28.D0
     &    *b*pi**2 + 216.D0*vlpm*vlsm - 432.D0*vlpm*vlwm - 24.D0*vlpm
     &     - 108.D0*vlpm**2 - 432.D0*vdmp - 36.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG**2*xnsq**2 * (  - 24.D0*b*vlsm*vlwm
     &     + 6.D0*b*vlsm**2 - 24.D0*b*vltm - 24.D0*b*vdw + 2.D0*b*pi**2
     &     + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG**2 * ( 78.D0*b*vlpm**2 - 24.D0*b*vltm
     &     + 24.D0*b*vlwm + 312.D0*b*vlwm**2 + 312.D0*b*vdw - 26.D0*b*
     &    pi**2 + 204.D0*vlpm*vlsm - 408.D0*vlpm*vlwm - 24.D0*vlpm - 
     &    102.D0*vlpm**2 - 408.D0*vdmp - 34.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG**3*xnsq * ( 6.D0*b*vlpm**2 + 24.D0*b*
     &    vlwm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 + 12.D0*vlpm*vlsm - 24.D0
     &    *vlpm*vlwm - 6.D0*vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*INVG**3 * ( 6.D0*b*vlpm**2 + 24.D0*b*
     &    vlwm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 + 12.D0*vlpm*vlsm - 24.D0
     &    *vlpm*vlwm - 6.D0*vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*xnsq * (  - 32.D0*XLF*b*TR*xn + 48.D0*b*
     &    vlpm**2 + 48.D0*b*vltm + 48.D0*b*vlwm + 192.D0*b*vlwm**2 + 
     &    192.D0*b*vdw - 32.D0*b*TR*xn - 16.D0*b*pi**2 + 48.D0*b + 528.D
     &    0*vlpm*vlsm - 1056.D0*vlpm*vlwm - 600.D0*vlpm - 264.D0*
     &    vlpm**2 - 1056.D0*vdmp - 88.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1*xnsq**2 * ( 16.D0*b + 6.D0*vlpm**2 + 24.D0
     &    *vdmp + 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1 * ( 48.D0*b*vlpm**2 + 96.D0*b*vltm**2 + 96.
     &    D0*b*vlwm**2 + 96.D0*b*vdw + 96.D0*b*vdt - 16.D0*b*pi**2 + 
     &    288.D0*vlpm*vlsm + 192.D0*vlpm*vltm - 768.D0*vlpm*vlwm - 192.D
     &    0*vlpm - 144.D0*vlpm**2 - 576.D0*vdmp - 48.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG*xnsq * (  - 378.D0*b*vlpm**2 + 216.
     &    D0*b*vltm*vlwm - 264.D0*b*vltm + 156.D0*b*vltm**2 - 168.D0*b*
     &    vlwm - 1452.D0*b*vlwm**2 - 1344.D0*b*vdw + 264.D0*b*vdt + 90.D
     &    0*b*pi**2 - 72.D0*b - 2436.D0*vlpm*vlsm + 48.D0*vlpm*vltm + 
     &    4824.D0*vlpm*vlwm + 1356.D0*vlpm + 1218.D0*vlpm**2 + 4872.D0*
     &    vdmp + 406.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG*xnsq**2 * (  - 24.D0*b*vlsm*vltm
     &     + 24.D0*b*vlsm*vlwm + 48.D0*b*vltm - 48.D0*b*vlwm + 24.D0*b*
     &    vdw - 24.D0*b*vdt - 81.D0*vlpm**2 - 324.D0*vdmp - 27.D0*pi**2
     &     )
      gcoeff2 = gcoeff2 + t1**2*INVG * (  - 324.D0*b*vlpm**2 + 144.D0*b
     &    *vltm - 360.D0*b*vltm**2 - 144.D0*b*vlwm - 936.D0*b*vlwm**2
     &     - 936.D0*b*vdw - 360.D0*b*vdt + 108.D0*b*pi**2 - 1440.D0*
     &    vlpm*vlsm - 600.D0*vlpm*vltm + 3480.D0*vlpm*vlwm + 576.D0*
     &    vlpm + 720.D0*vlpm**2 + 2880.D0*vdmp + 240.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG**2*xnsq * (  - 390.D0*b*vlpm**2 + 
     &    120.D0*b*vltm*vlwm - 216.D0*b*vltm + 60.D0*b*vltm**2 - 96.D0*
     &    b*vlwm - 1500.D0*b*vlwm**2 - 1440.D0*b*vdw + 120.D0*b*vdt + 
     &    110.D0*b*pi**2 - 1308.D0*vlpm*vlsm + 2616.D0*vlpm*vlwm + 216.D
     &    0*vlpm + 654.D0*vlpm**2 + 2616.D0*vdmp + 218.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG**2*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm + 120.D0*b*vlsm*vlwm - 18.D0*b*vlsm**2 + 96.D0*b*vltm - 
     &    24.D0*b*vlwm + 120.D0*b*vdw - 48.D0*b*vdt - 6.D0*b*pi**2 - 66.
     &    D0*vlpm**2 - 264.D0*vdmp - 22.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG**2 * (  - 330.D0*b*vlpm**2 + 120.D0
     &    *b*vltm - 72.D0*b*vltm**2 - 120.D0*b*vlwm - 1248.D0*b*vlwm**2
     &     - 1248.D0*b*vdw - 72.D0*b*vdt + 110.D0*b*pi**2 - 1068.D0*
     &    vlpm*vlsm - 72.D0*vlpm*vltm + 2208.D0*vlpm*vlwm + 168.D0*vlpm
     &     + 534.D0*vlpm**2 + 2136.D0*vdmp + 178.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG**3*xnsq * (  - 60.D0*b*vlpm**2 - 
     &    240.D0*b*vlwm**2 - 240.D0*b*vdw + 20.D0*b*pi**2 - 144.D0*vlpm
     &    *vlsm + 288.D0*vlpm*vlwm + 72.D0*vlpm**2 + 288.D0*vdmp + 24.D0
     &    *pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG**3*xnsq**2 * ( 24.D0*b*vlsm*vlwm
     &     - 6.D0*b*vlsm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 - 6.D0*vlpm**2
     &     - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2*INVG**3 * (  - 54.D0*b*vlpm**2 - 216.D0
     &    *b*vlwm**2 - 216.D0*b*vdw + 18.D0*b*pi**2 - 132.D0*vlpm*vlsm
     &     + 264.D0*vlpm*vlwm + 66.D0*vlpm**2 + 264.D0*vdmp + 22.D0*
     &    pi**2 )
      gcoeff2 = gcoeff2 + t1**2*xnsq * ( 96.D0*b*vltm**2 - 96.D0*b*
     &    vlwm**2 - 96.D0*b*vdw + 96.D0*b*vdt - 720.D0*vlpm*vlsm + 192.D
     &    0*vlpm*vltm + 1248.D0*vlpm*vlwm + 1008.D0*vlpm + 360.D0*
     &    vlpm**2 + 1440.D0*vdmp + 120.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**2 * (  - 480.D0*vlpm*vltm + 480.D0*vlpm*
     &    vlwm )
      gcoeff2 = gcoeff2 + t1**3*INVG*xnsq * ( 252.D0*b*vlpm**2 - 144.D0
     &    *b*vltm*vlwm + 144.D0*b*vltm - 432.D0*b*vltm**2 + 144.D0*b*
     &    vlwm + 1296.D0*b*vlwm**2 + 1224.D0*b*vdw - 504.D0*b*vdt - 60.D
     &    0*b*pi**2 + 48.D0*b + 3864.D0*vlpm*vlsm - 600.D0*vlpm*vltm - 
     &    7128.D0*vlpm*vlwm - 2664.D0*vlpm - 1932.D0*vlpm**2 - 7728.D0*
     &    vdmp - 644.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG*xnsq**2 * ( 54.D0*vlpm**2 + 216.D0
     &    *vdmp + 18.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG * ( 216.D0*b*vlpm**2 + 432.D0*b*
     &    vltm**2 + 432.D0*b*vlwm**2 + 432.D0*b*vdw + 432.D0*b*vdt - 72.
     &    D0*b*pi**2 + 960.D0*vlpm*vlsm + 1728.D0*vlpm*vltm - 3648.D0*
     &    vlpm*vlwm - 384.D0*vlpm - 480.D0*vlpm**2 - 1920.D0*vdmp - 160.
     &    D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG**2*xnsq * ( 840.D0*b*vlpm**2 - 480.
     &    D0*b*vltm*vlwm + 360.D0*b*vltm - 312.D0*b*vltm**2 + 168.D0*b*
     &    vlwm + 3192.D0*b*vlwm**2 + 2952.D0*b*vdw - 552.D0*b*vdt - 200.
     &    D0*b*pi**2 + 3936.D0*vlpm*vlsm - 72.D0*vlpm*vltm - 7800.D0*
     &    vlpm*vlwm - 816.D0*vlpm - 1968.D0*vlpm**2 - 7872.D0*vdmp - 
     &    656.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG**2*xnsq**2 * ( 120.D0*b*vlsm*vltm
     &     - 168.D0*b*vlsm*vlwm + 12.D0*b*vlsm**2 - 120.D0*b*vltm + 72.D
     &    0*b*vlwm - 168.D0*b*vdw + 120.D0*b*vdt + 4.D0*b*pi**2 + 204.D0
     &    *vlpm**2 + 816.D0*vdmp + 68.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG**2 * ( 660.D0*b*vlpm**2 - 192.D0*b
     &    *vltm + 456.D0*b*vltm**2 + 192.D0*b*vlwm + 2184.D0*b*vlwm**2
     &     + 2184.D0*b*vdw + 456.D0*b*vdt - 220.D0*b*pi**2 + 2472.D0*
     &    vlpm*vlsm + 648.D0*vlpm*vltm - 5592.D0*vlpm*vlwm - 432.D0*
     &    vlpm - 1236.D0*vlpm**2 - 4944.D0*vdmp - 412.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG**3*xnsq * ( 252.D0*b*vlpm**2 - 48.D
     &    0*b*vltm*vlwm - 24.D0*b*vltm**2 + 984.D0*b*vlwm**2 + 960.D0*b
     &    *vdw - 48.D0*b*vdt - 76.D0*b*pi**2 + 744.D0*vlpm*vlsm - 1488.D
     &    0*vlpm*vlwm - 372.D0*vlpm**2 - 1488.D0*vdmp - 124.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG**3*xnsq**2 * ( 24.D0*b*vlsm*vltm
     &     - 120.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 120.D0*b*vdw + 24.D0
     &    *b*vdt + 8.D0*b*pi**2 + 48.D0*vlpm**2 + 192.D0*vdmp + 16.D0*
     &    pi**2 )
      gcoeff2 = gcoeff2 + t1**3*INVG**3 * ( 204.D0*b*vlpm**2 + 24.D0*b*
     &    vltm**2 + 792.D0*b*vlwm**2 + 792.D0*b*vdw + 24.D0*b*vdt - 68.D
     &    0*b*pi**2 + 600.D0*vlpm*vlsm + 24.D0*vlpm*vltm - 1224.D0*vlpm
     &    *vlwm - 300.D0*vlpm**2 - 1200.D0*vdmp - 100.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**3*xnsq * ( 480.D0*vlpm*vlsm - 480.D0*vlpm
     &    *vltm - 480.D0*vlpm*vlwm - 672.D0*vlpm - 240.D0*vlpm**2 - 960.
     &    D0*vdmp - 80.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG*xnsq * ( 432.D0*b*vltm**2 - 432.D0
     &    *b*vlwm**2 - 432.D0*b*vdw + 432.D0*b*vdt - 3360.D0*vlpm*vlsm
     &     + 1728.D0*vlpm*vltm + 4992.D0*vlpm*vlwm + 2640.D0*vlpm + 
     &    1680.D0*vlpm**2 + 6720.D0*vdmp + 560.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG * (  - 1344.D0*vlpm*vltm + 1344.D0
     &    *vlpm*vlwm )
      gcoeff2 = gcoeff2 + t1**4*INVG**2*xnsq * (  - 870.D0*b*vlpm**2 + 
     &    600.D0*b*vltm*vlwm - 288.D0*b*vltm + 756.D0*b*vltm**2 - 192.D0
     &    *b*vlwm - 3636.D0*b*vlwm**2 - 3336.D0*b*vdw + 1056.D0*b*vdt
     &     + 190.D0*b*pi**2 - 6780.D0*vlpm*vlsm + 648.D0*vlpm*vltm + 
     &    12912.D0*vlpm*vlwm + 1680.D0*vlpm + 3390.D0*vlpm**2 + 13560.D0
     &    *vdmp + 1130.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG**2*xnsq**2 * (  - 72.D0*b*vlsm*
     &    vltm + 72.D0*b*vlsm*vlwm + 48.D0*b*vltm - 48.D0*b*vlwm + 72.D0
     &    *b*vdw - 72.D0*b*vdt - 240.D0*vlpm**2 - 960.D0*vdmp - 80.D0*
     &    pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG**2 * (  - 660.D0*b*vlpm**2 + 96.D0
     &    *b*vltm - 888.D0*b*vltm**2 - 96.D0*b*vlwm - 1752.D0*b*vlwm**2
     &     - 1752.D0*b*vdw - 888.D0*b*vdt + 220.D0*b*pi**2 - 2640.D0*
     &    vlpm*vlsm - 2040.D0*vlpm*vltm + 7320.D0*vlpm*vlwm + 480.D0*
     &    vlpm + 1320.D0*vlpm**2 + 5280.D0*vdmp + 440.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG**3*xnsq * (  - 570.D0*b*vlpm**2 + 
     &    240.D0*b*vltm*vlwm + 144.D0*b*vltm**2 - 2184.D0*b*vlwm**2 - 
     &    2064.D0*b*vdw + 264.D0*b*vdt + 150.D0*b*pi**2 - 2172.D0*vlpm*
     &    vlsm + 24.D0*vlpm*vltm + 4320.D0*vlpm*vlwm + 1086.D0*vlpm**2
     &     + 4344.D0*vdmp + 362.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG**3*xnsq**2 * (  - 96.D0*b*vlsm*
     &    vltm + 216.D0*b*vlsm*vlwm - 30.D0*b*vlsm**2 + 216.D0*b*vdw - 
     &    96.D0*b*vdt - 10.D0*b*pi**2 - 150.D0*vlpm**2 - 600.D0*vdmp - 
     &    50.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**4*INVG**3 * (  - 420.D0*b*vlpm**2 - 168.D0
     &    *b*vltm**2 - 1512.D0*b*vlwm**2 - 1512.D0*b*vdw - 168.D0*b*vdt
     &     + 140.D0*b*pi**2 - 1440.D0*vlpm*vlsm - 216.D0*vlpm*vltm + 
     &    3096.D0*vlpm*vlwm + 720.D0*vlpm**2 + 2880.D0*vdmp + 240.D0*
     &    pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG*xnsq * ( 1344.D0*vlpm*vlsm - 1344.D
     &    0*vlpm*vltm - 1344.D0*vlpm*vlwm - 1056.D0*vlpm - 672.D0*
     &    vlpm**2 - 2688.D0*vdmp - 224.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG**2*xnsq * ( 348.D0*b*vlpm**2 - 240.
     &    D0*b*vltm*vlwm + 96.D0*b*vltm - 1008.D0*b*vltm**2 + 96.D0*b*
     &    vlwm + 2160.D0*b*vlwm**2 + 2040.D0*b*vdw - 1128.D0*b*vdt - 76.
     &    D0*b*pi**2 + 7080.D0*vlpm*vlsm - 2040.D0*vlpm*vltm - 12120.D0
     &    *vlpm*vlwm - 2016.D0*vlpm - 3540.D0*vlpm**2 - 14160.D0*vdmp
     &     - 1180.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG**2*xnsq**2 * ( 96.D0*vlpm**2 + 384.
     &    D0*vdmp + 32.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG**2 * ( 264.D0*b*vlpm**2 + 528.D0*b
     &    *vltm**2 + 528.D0*b*vlwm**2 + 528.D0*b*vdw + 528.D0*b*vdt - 
     &    88.D0*b*pi**2 + 1056.D0*vlpm*vlsm + 2688.D0*vlpm*vltm - 4800.D
     &    0*vlpm*vlwm - 192.D0*vlpm - 528.D0*vlpm**2 - 2112.D0*vdmp - 
     &    176.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG**3*xnsq * ( 732.D0*b*vlpm**2 - 432.
     &    D0*b*vltm*vlwm - 384.D0*b*vltm**2 + 2880.D0*b*vlwm**2 + 2664.D
     &    0*b*vdw - 600.D0*b*vdt - 172.D0*b*pi**2 + 3960.D0*vlpm*vlsm
     &     - 216.D0*vlpm*vltm - 7704.D0*vlpm*vlwm - 1980.D0*vlpm**2 - 
     &    7920.D0*vdmp - 660.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG**3*xnsq**2 * ( 120.D0*b*vlsm*vltm
     &     - 168.D0*b*vlsm*vlwm + 12.D0*b*vlsm**2 - 168.D0*b*vdw + 120.D
     &    0*b*vdt + 4.D0*b*pi**2 + 228.D0*vlpm**2 + 912.D0*vdmp + 76.D0
     &    *pi**2 )
      gcoeff2 = gcoeff2 + t1**5*INVG**3 * ( 504.D0*b*vlpm**2 + 432.D0*b
     &    *vltm**2 + 1584.D0*b*vlwm**2 + 1584.D0*b*vdw + 432.D0*b*vdt
     &     - 168.D0*b*pi**2 + 1920.D0*vlpm*vlsm + 768.D0*vlpm*vltm - 
     &    4608.D0*vlpm*vlwm - 960.D0*vlpm**2 - 3840.D0*vdmp - 320.D0*
     &    pi**2 )
      gcoeff2 = gcoeff2 + t1**6*INVG**2*xnsq * ( 528.D0*b*vltm**2 - 528.
     &    D0*b*vlwm**2 - 528.D0*b*vdw + 528.D0*b*vdt - 4368.D0*vlpm*
     &    vlsm + 2688.D0*vlpm*vltm + 6048.D0*vlpm*vlwm + 1344.D0*vlpm
     &     + 2184.D0*vlpm**2 + 8736.D0*vdmp + 728.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**6*INVG**2 * (  - 1248.D0*vlpm*vltm + 1248.
     &    D0*vlpm*vlwm )
      gcoeff2 = gcoeff2 + t1**6*INVG**3*xnsq * (  - 504.D0*b*vlpm**2 + 
     &    336.D0*b*vltm*vlwm + 600.D0*b*vltm**2 - 2280.D0*b*vlwm**2 - 
     &    2112.D0*b*vdw + 768.D0*b*vdt + 112.D0*b*pi**2 - 4704.D0*vlpm*
     &    vlsm + 768.D0*vlpm*vltm + 8640.D0*vlpm*vlwm + 2352.D0*vlpm**2
     &     + 9408.D0*vdmp + 784.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**6*INVG**3*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm + 48.D0*b*vlsm*vlwm + 48.D0*b*vdw - 48.D0*b*vdt - 168.D0
     &    *vlpm**2 - 672.D0*vdmp - 56.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**6*INVG**3 * (  - 336.D0*b*vlpm**2 - 480.D0
     &    *b*vltm**2 - 864.D0*b*vlwm**2 - 864.D0*b*vdw - 480.D0*b*vdt
     &     + 112.D0*b*pi**2 - 1344.D0*vlpm*vlsm - 1344.D0*vlpm*vltm + 
     &    4032.D0*vlpm*vlwm + 672.D0*vlpm**2 + 2688.D0*vdmp + 224.D0*
     &    pi**2 )
      gcoeff2 = gcoeff2 + t1**7*INVG**2*xnsq * ( 1248.D0*vlpm*vlsm - 
     &    1248.D0*vlpm*vltm - 1248.D0*vlpm*vlwm - 384.D0*vlpm - 624.D0*
     &    vlpm**2 - 2496.D0*vdmp - 208.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**7*INVG**3*xnsq * ( 144.D0*b*vlpm**2 - 96.D
     &    0*b*vltm*vlwm - 528.D0*b*vltm**2 + 1008.D0*b*vlwm**2 + 960.D0
     &    *b*vdw - 576.D0*b*vdt - 32.D0*b*pi**2 + 3648.D0*vlpm*vlsm - 
     &    1344.D0*vlpm*vltm - 5952.D0*vlpm*vlwm - 1824.D0*vlpm**2 - 
     &    7296.D0*vdmp - 608.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**7*INVG**3*xnsq**2 * ( 48.D0*vlpm**2 + 192.
     &    D0*vdmp + 16.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**7*INVG**3 * ( 96.D0*b*vlpm**2 + 192.D0*b*
     &    vltm**2 + 192.D0*b*vlwm**2 + 192.D0*b*vdw + 192.D0*b*vdt - 32.
     &    D0*b*pi**2 + 384.D0*vlpm*vlsm + 1152.D0*vlpm*vltm - 1920.D0*
     &    vlpm*vlwm - 192.D0*vlpm**2 - 768.D0*vdmp - 64.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**8*INVG**3*xnsq * ( 192.D0*b*vltm**2 - 192.
     &    D0*b*vlwm**2 - 192.D0*b*vdw + 192.D0*b*vdt - 1728.D0*vlpm*
     &    vlsm + 1152.D0*vlpm*vltm + 2304.D0*vlpm*vlwm + 864.D0*vlpm**2
     &     + 3456.D0*vdmp + 288.D0*pi**2 )
      gcoeff2 = gcoeff2 + t1**8*INVG**3 * (  - 384.D0*vlpm*vltm + 384.D0
     &    *vlpm*vlwm )
      gcoeff2 = gcoeff2 + t1**9*INVG**3*xnsq * ( 384.D0*vlpm*vlsm - 384.
     &    D0*vlpm*vltm - 384.D0*vlpm*vlwm - 192.D0*vlpm**2 - 768.D0*
     &    vdmp - 64.D0*pi**2 )
      gcoeff2 = gcoeff2 + t2**(-2)*ro*xnsq * ( 6.D0*b*vlwm**2 + 6.D0*b*
     &    vdw + b*pi**2 )
      gcoeff2 = gcoeff2 + t2**(-2)*ro * (  - 6.D0*b*vlwm**2 - 6.D0*b*
     &    vdw - b*pi**2 )
      gcoeff2 = gcoeff2 + t2**(-1)*ro*xnsq * ( 6.D0*b*vlwm**2 + 6.D0*b*
     &    vdw + b*pi**2 )
      gcoeff2 = gcoeff2 + t2**(-1)*ro * (  - 12.D0*b*vlwm**2 - 12.D0*b*
     &    vdw - 2.D0*b*pi**2 )
      gcoeff2 = gcoeff2 + t2**(-1)*xnsq * ( 24.D0*b*vlwm - 24.D0*b )
      gcoeff2 = gcoeff2 + t2**(-1) * (  - 24.D0*b*vlwm + 24.D0*b )
      gcoeff2 = gcoeff2 + ro*xnsq * ( 12.D0*b*vlpm**2*TR*xn + 12.D0*b*
     &    vltm**2 - 12.D0*b*vlwm**2 - 12.D0*b*vdw + 12.D0*b*vdt - 12.D0
     &    *b*TR*pi**2*xn + 96.D0*b*TR*xn - 48.D0*vlpm*TR*xn - 36.D0*
     &    vlpm )
      gcoeff2 = gcoeff2 + ro**2*xnsq * ( 48.D0*vlpm*TR*xn )
      gcoeff2 = gcoeff2 + TBAR**(-1)*xnsq * (  - 12.D0*b*vltm )
      gcoeff2 = gcoeff2 + TBAR**(-1)*xnsq**2 * ( 12.D0*b*vltm )
      gcoeff2 = gcoeff2 + UBAR**(-1)*INVG*xnsq * (  - 24.D0*b*vlwm )
      gcoeff2 = gcoeff2 + UBAR**(-1)*INVG*xnsq**2 * ( 12.D0*b*vlwm )
      gcoeff2 = gcoeff2 + UBAR**(-1)*INVG * ( 12.D0*b*vlwm )
      gcoeff2 = gcoeff2 + UBAR**(-1)*xnsq * ( 12.D0*b*vlwm )
      gcoeff2 = gcoeff2 + UBAR**(-1)*xnsq**2 * (  - 12.D0*b*vlwm )
      gcoeff2 = gcoeff2 + INVG*xnsq * (  - 30.D0*b*vlpm**2 + 12.D0*b*
     &    vltm - 48.D0*b*vlwm - 120.D0*b*vlwm**2 - 120.D0*b*vdw + 10.D0
     &    *b*pi**2 - 24.D0*b - 84.D0*vlpm*vlsm + 168.D0*vlpm*vlwm + 24.D
     &    0*vlpm + 42.D0*vlpm**2 + 168.D0*vdmp + 14.D0*pi**2 )
      gcoeff2 = gcoeff2 + INVG*xnsq**2 * (  - 3.D0*b*vlsm*omro**(-1) - 
     &    3.D0*b*vlsm + 12.D0*b*vlwm + 12.D0*b + 3.D0/4.D0*vlpm**2*
     &    omro**(-1) + 3.D0/4.D0*vlpm**2 + 3.D0*vdmp*omro**(-1) + 3.D0*
     &    vdmp + 1.D0/4.D0*omro**(-1)*pi**2 + 1.D0/4.D0*pi**2 )
      gcoeff2 = gcoeff2 + INVG * (  - 30.D0*b*vlpm**2 - 12.D0*b*vltm - 
     &    12.D0*b*vlwm - 120.D0*b*vlwm**2 - 120.D0*b*vdw + 10.D0*b*
     &    pi**2 + 12.D0*b - 84.D0*vlpm*vlsm + 168.D0*vlpm*vlwm + 24.D0*
     &    vlpm + 42.D0*vlpm**2 + 168.D0*vdmp + 14.D0*pi**2 )
      gcoeff2 = gcoeff2 + INVG**2*xnsq * (  - 6.D0*b*vlpm**2 - 24.D0*b*
     &    vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 24.D0
     &    *vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + INVG**2 * (  - 6.D0*b*vlpm**2 - 24.D0*b*
     &    vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 24.D0
     &    *vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff2 = gcoeff2 + xnsq * ( 16.D0*XLF*b*TR*xn - 24.D0*b*vlpm**2
     &     - 24.D0*b*vltm - 24.D0*b*vlwm - 96.D0*b*vlwm**2 - 96.D0*b*
     &    vdw + 16.D0*b*TR*xn + 8.D0*b*pi**2 - 24.D0*b - 144.D0*vlpm*
     &    vlsm + 288.D0*vlpm*vlwm + 132.D0*vlpm + 72.D0*vlpm**2 + 288.D0
     &    *vdmp + 24.D0*pi**2 )
      gcoeff2 = gcoeff2 + xnsq**2 * (  - 8.D0*b - 3.D0*vlpm**2 - 12.D0*
     &    vdmp - pi**2 )
      gcoeff2 = gcoeff2 - 24.D0*b*vlpm**2 + 48.D0*b*vltm - 48.D0*b*vlwm
     &     - 96.D0*b*vlwm**2 - 96.D0*b*vdw + 8.D0*b*pi**2 - 144.D0*vlpm
     &    *vlsm + 288.D0*vlpm*vlwm + 96.D0*vlpm + 72.D0*vlpm**2 + 288.D0
     &    *vdmp + 24.D0*pi**2

      gcoeff3 =  + t1**(-1)*ro*xnsq * (  - 12.D0*b*vltm**2 - 12.D0*b*
     &    vdt - 2.D0*b*pi**2 )
      gcoeff3 = gcoeff3 + t1**(-1)*ro * ( 12.D0*b*vltm**2 + 12.D0*b*vdt
     &     + 2.D0*b*pi**2 )
      gcoeff3 = gcoeff3 + t1**(-1)*UBAR**(-1)*INVG*xnsq * ( 12.D0*b*
     &    vlwm - 24.D0*b )
      gcoeff3 = gcoeff3 + t1**(-1)*UBAR**(-1)*INVG*xnsq**2 * ( 12.D0*b
     &     )
      gcoeff3 = gcoeff3 + t1**(-1)*UBAR**(-1)*INVG * (  - 12.D0*b*vlwm
     &     + 12.D0*b )
      gcoeff3 = gcoeff3 + t1**(-1)*UBAR**(-1)*xnsq * ( 12.D0*b*vlwm - 
     &    24.D0*b )
      gcoeff3 = gcoeff3 + t1**(-1)*UBAR**(-1)*xnsq**2 * ( 12.D0*b )
      gcoeff3 = gcoeff3 + t1**(-1)*UBAR**(-1) * (  - 12.D0*b*vlwm + 12.D
     &    0*b )
      gcoeff3 = gcoeff3 + t1**(-1)*INVG*xnsq * (  - 12.D0*b*vlwm + 24.D0
     &    *b )
      gcoeff3 = gcoeff3 + t1**(-1)*INVG*xnsq**2 * ( 9.D0/2.D0*b*vlsm*
     &    omro**(-1) - 9.D0/2.D0*b*vlsm - 12.D0*b - 9.D0/8.D0*vlpm**2*
     &    omro**(-1) + 9.D0/8.D0*vlpm**2 - 9.D0/2.D0*vdmp*omro**(-1) + 
     &    9.D0/2.D0*vdmp - 3.D0/8.D0*omro**(-1)*pi**2 + 3.D0/8.D0*pi**2
     &     )
      gcoeff3 = gcoeff3 + t1**(-1)*INVG * ( 12.D0*b*vlwm - 12.D0*b )
      gcoeff3 = gcoeff3 + t1**(-1)*xnsq**2 * ( 18.D0*b*vlsm*omro**(-1)
     &     - 9.D0/2.D0*vlpm**2*omro**(-1) - 18.D0*vdmp*omro**(-1) - 3.D0
     &    /2.D0*omro**(-1)*pi**2 )
      gcoeff3 = gcoeff3 + t1*ro*xnsq * (  - 24.D0*b*vlpm**2*TR*xn + 24.D
     &    0*b*TR*pi**2*xn - 192.D0*b*TR*xn + 96.D0*vlpm*TR*xn + 72.D0*
     &    vlpm )
      gcoeff3 = gcoeff3 + t1*ro**2*xnsq * (  - 96.D0*vlpm*TR*xn )
      gcoeff3 = gcoeff3 + t1*UBAR**(-2)*xnsq * (  - 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + t1*UBAR**(-2)*xnsq**2 * ( 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + t1*INVG*xnsq * ( 96.D0*b*vlpm**2 + 48.D0*b*
     &    vltm*vlwm + 24.D0*b*vltm + 24.D0*b*vltm**2 + 96.D0*b*vlwm + 
     &    408.D0*b*vlwm**2 + 432.D0*b*vdw + 48.D0*b*vdt - 40.D0*b*pi**2
     &     + 48.D0*b + 372.D0*vlpm*vlsm - 744.D0*vlpm*vlwm - 120.D0*
     &    vlpm - 186.D0*vlpm**2 - 744.D0*vdmp - 62.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1*INVG*xnsq**2 * (  - 24.D0*b*vlsm*vltm - 24.
     &    D0*b*vlsm*vlwm + 6.D0*b*vlsm + 12.D0*b*vlsm**2 - 24.D0*b*vlwm
     &     - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0*b*pi**2 - 24.D0*b + 9.D0/
     &    2.D0*vlpm**2 + 18.D0*vdmp + 3.D0/2.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1*INVG * ( 72.D0*b*vlpm**2 - 24.D0*b*vltm - 
     &    24.D0*b*vltm**2 + 72.D0*b*vlwm + 312.D0*b*vlwm**2 + 312.D0*b*
     &    vdw - 24.D0*b*vdt - 24.D0*b*pi**2 - 24.D0*b + 336.D0*vlpm*
     &    vlsm - 24.D0*vlpm*vltm - 648.D0*vlpm*vlwm - 96.D0*vlpm - 168.D
     &    0*vlpm**2 - 672.D0*vdmp - 56.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1*INVG**2*xnsq * ( 33.D0*b*vlpm**2 + 132.D0*
     &    b*vlwm**2 + 132.D0*b*vdw - 11.D0*b*pi**2 + 78.D0*vlpm*vlsm - 
     &    156.D0*vlpm*vlwm - 39.D0*vlpm**2 - 156.D0*vdmp - 13.D0*pi**2
     &     )
      gcoeff3 = gcoeff3 + t1*INVG**2*xnsq**2 * (  - 12.D0*b*vlsm*vlwm
     &     + 3.D0*b*vlsm**2 - 12.D0*b*vdw + b*pi**2 + 3.D0*vlpm**2 + 12.
     &    D0*vdmp + pi**2 )
      gcoeff3 = gcoeff3 + t1*INVG**2 * ( 30.D0*b*vlpm**2 + 120.D0*b*
     &    vlwm**2 + 120.D0*b*vdw - 10.D0*b*pi**2 + 72.D0*vlpm*vlsm - 
     &    144.D0*vlpm*vlwm - 36.D0*vlpm**2 - 144.D0*vdmp - 12.D0*pi**2
     &     )
      gcoeff3 = gcoeff3 + t1*xnsq * (  - 32.D0*XLF*b*TR*xn + 48.D0*b*
     &    vlpm**2 + 48.D0*b*vltm + 72.D0*b*vltm**2 + 48.D0*b*vlwm + 120.
     &    D0*b*vlwm**2 + 120.D0*b*vdw + 72.D0*b*vdt - 32.D0*b*TR*xn - 
     &    16.D0*b*pi**2 + 48.D0*b + 288.D0*vlpm*vlsm - 576.D0*vlpm*vlwm
     &     - 360.D0*vlpm - 144.D0*vlpm**2 - 576.D0*vdmp - 48.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1*xnsq**2 * ( 16.D0*b + 6.D0*vlpm**2 + 24.D0
     &    *vdmp + 2.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1 * ( 48.D0*b*vlpm**2 + 96.D0*b*vltm**2 + 96.
     &    D0*b*vlwm**2 + 96.D0*b*vdw + 96.D0*b*vdt - 16.D0*b*pi**2 + 
     &    288.D0*vlpm*vlsm - 144.D0*vlpm*vltm - 432.D0*vlpm*vlwm - 192.D
     &    0*vlpm - 144.D0*vlpm**2 - 576.D0*vdmp - 48.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG*xnsq * (  - 240.D0*b*vlpm**2 - 168.
     &    D0*b*vltm - 24.D0*b*vltm**2 - 72.D0*b*vlwm - 936.D0*b*vlwm**2
     &     - 936.D0*b*vdw - 24.D0*b*vdt + 80.D0*b*pi**2 - 24.D0*b - 
     &    1260.D0*vlpm*vlsm - 24.D0*vlpm*vltm + 2544.D0*vlpm*vlwm + 564.
     &    D0*vlpm + 630.D0*vlpm**2 + 2520.D0*vdmp + 210.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG*xnsq**2 * ( 24.D0*b*vlsm*vltm + 72.
     &    D0*b*vlsm*vlwm + 48.D0*b*vlsm - 24.D0*b*vlsm**2 - 96.D0*b*
     &    vlwm + 72.D0*b*vdw + 24.D0*b*vdt - 8.D0*b*pi**2 - 45.D0*
     &    vlpm**2 - 180.D0*vdmp - 15.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG * (  - 168.D0*b*vlpm**2 + 144.D0*b
     &    *vltm - 48.D0*b*vltm**2 - 144.D0*b*vlwm - 624.D0*b*vlwm**2 - 
     &    624.D0*b*vdw - 48.D0*b*vdt + 56.D0*b*pi**2 - 936.D0*vlpm*vlsm
     &     - 24.D0*vlpm*vltm + 1896.D0*vlpm*vlwm + 336.D0*vlpm + 468.D0
     &    *vlpm**2 + 1872.D0*vdmp + 156.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG**2*xnsq * (  - 165.D0*b*vlpm**2 + 
     &    24.D0*b*vltm*vlwm - 48.D0*b*vltm + 12.D0*b*vltm**2 - 24.D0*b*
     &    vlwm - 648.D0*b*vlwm**2 - 636.D0*b*vdw + 24.D0*b*vdt + 51.D0*
     &    b*pi**2 - 486.D0*vlpm*vlsm + 972.D0*vlpm*vlwm + 24.D0*vlpm + 
     &    243.D0*vlpm**2 + 972.D0*vdmp + 81.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG**2*xnsq**2 * (  - 12.D0*b*vlsm*
     &    vltm + 72.D0*b*vlsm*vlwm - 15.D0*b*vlsm**2 + 24.D0*b*vltm + 
     &    72.D0*b*vdw - 12.D0*b*vdt - 5.D0*b*pi**2 - 27.D0*vlpm**2 - 
     &    108.D0*vdmp - 9.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG**2 * (  - 138.D0*b*vlpm**2 + 24.D0
     &    *b*vltm - 12.D0*b*vltm**2 - 24.D0*b*vlwm - 540.D0*b*vlwm**2
     &     - 540.D0*b*vdw - 12.D0*b*vdt + 46.D0*b*pi**2 - 408.D0*vlpm*
     &    vlsm - 12.D0*vlpm*vltm + 828.D0*vlpm*vlwm + 24.D0*vlpm + 204.D
     &    0*vlpm**2 + 816.D0*vdmp + 68.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG**3*xnsq * (  - 6.D0*b*vlpm**2 - 24.
     &    D0*b*vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm
     &     + 24.D0*vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*INVG**3 * (  - 6.D0*b*vlpm**2 - 24.D0*b
     &    *vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 24.D
     &    0*vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2*xnsq * ( 96.D0*b*vltm**2 - 96.D0*b*
     &    vlwm**2 - 96.D0*b*vdw + 96.D0*b*vdt - 384.D0*vlpm*vlsm - 144.D
     &    0*vlpm*vltm + 912.D0*vlpm*vlwm + 528.D0*vlpm + 192.D0*vlpm**2
     &     + 768.D0*vdmp + 64.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**2 * (  - 480.D0*vlpm*vltm + 480.D0*vlpm*
     &    vlwm )
      gcoeff3 = gcoeff3 + t1**3*INVG*xnsq * ( 252.D0*b*vlpm**2 - 144.D0
     &    *b*vltm*vlwm + 144.D0*b*vltm - 120.D0*b*vltm**2 + 144.D0*b*
     &    vlwm + 984.D0*b*vlwm**2 + 912.D0*b*vdw - 192.D0*b*vdt - 60.D0
     &    *b*pi**2 + 48.D0*b + 2280.D0*vlpm*vlsm - 24.D0*vlpm*vltm - 
     &    4536.D0*vlpm*vlwm - 1368.D0*vlpm - 1140.D0*vlpm**2 - 4560.D0*
     &    vdmp - 380.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG*xnsq**2 * ( 54.D0*vlpm**2 + 216.D0
     &    *vdmp + 18.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG * ( 216.D0*b*vlpm**2 + 432.D0*b*
     &    vltm**2 + 432.D0*b*vlwm**2 + 432.D0*b*vdw + 432.D0*b*vdt - 72.
     &    D0*b*pi**2 + 960.D0*vlpm*vlsm + 720.D0*vlpm*vltm - 2640.D0*
     &    vlpm*vlwm - 384.D0*vlpm - 480.D0*vlpm**2 - 1920.D0*vdmp - 160.
     &    D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG**2*xnsq * ( 438.D0*b*vlpm**2 - 144.
     &    D0*b*vltm*vlwm + 168.D0*b*vltm - 84.D0*b*vltm**2 + 72.D0*b*
     &    vlwm + 1692.D0*b*vlwm**2 + 1620.D0*b*vdw - 156.D0*b*vdt - 122.
     &    D0*b*pi**2 + 1716.D0*vlpm*vlsm - 12.D0*vlpm*vltm - 3420.D0*
     &    vlpm*vlwm - 192.D0*vlpm - 858.D0*vlpm**2 - 3432.D0*vdmp - 286.
     &    D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG**2*xnsq**2 * ( 36.D0*b*vlsm*vltm
     &     - 180.D0*b*vlsm*vlwm + 36.D0*b*vlsm**2 - 72.D0*b*vltm + 24.D0
     &    *b*vlwm - 180.D0*b*vdw + 36.D0*b*vdt + 12.D0*b*pi**2 + 108.D0
     &    *vlpm**2 + 432.D0*vdmp + 36.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG**2 * ( 336.D0*b*vlpm**2 - 96.D0*b*
     &    vltm + 120.D0*b*vltm**2 + 96.D0*b*vlwm + 1224.D0*b*vlwm**2 + 
     &    1224.D0*b*vdw + 120.D0*b*vdt - 112.D0*b*pi**2 + 1224.D0*vlpm*
     &    vlsm + 144.D0*vlpm*vltm - 2592.D0*vlpm*vlwm - 144.D0*vlpm - 
     &    612.D0*vlpm**2 - 2448.D0*vdmp - 204.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG**3*xnsq * ( 54.D0*b*vlpm**2 + 216.D
     &    0*b*vlwm**2 + 216.D0*b*vdw - 18.D0*b*pi**2 + 132.D0*vlpm*vlsm
     &     - 264.D0*vlpm*vlwm - 66.D0*vlpm**2 - 264.D0*vdmp - 22.D0*
     &    pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG**3*xnsq**2 * (  - 24.D0*b*vlsm*
     &    vlwm + 6.D0*b*vlsm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 + 6.D0*
     &    vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**3*INVG**3 * ( 48.D0*b*vlpm**2 + 192.D0*b*
     &    vlwm**2 + 192.D0*b*vdw - 16.D0*b*pi**2 + 120.D0*vlpm*vlsm - 
     &    240.D0*vlpm*vlwm - 60.D0*vlpm**2 - 240.D0*vdmp - 20.D0*pi**2
     &     )
      gcoeff3 = gcoeff3 + t1**3*xnsq * ( 480.D0*vlpm*vlsm - 480.D0*vlpm
     &    *vltm - 480.D0*vlpm*vlwm - 672.D0*vlpm - 240.D0*vlpm**2 - 960.
     &    D0*vdmp - 80.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG*xnsq * ( 432.D0*b*vltm**2 - 432.D0
     &    *b*vlwm**2 - 432.D0*b*vdw + 432.D0*b*vdt - 2352.D0*vlpm*vlsm
     &     + 720.D0*vlpm*vltm + 3984.D0*vlpm*vlwm + 1776.D0*vlpm + 1176.
     &    D0*vlpm**2 + 4704.D0*vdmp + 392.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG * (  - 1344.D0*vlpm*vltm + 1344.D0
     &    *vlpm*vlwm )
      gcoeff3 = gcoeff3 + t1**4*INVG**2*xnsq * (  - 606.D0*b*vlpm**2 + 
     &    360.D0*b*vltm*vlwm - 192.D0*b*vltm + 300.D0*b*vltm**2 - 96.D0
     &    *b*vlwm - 2364.D0*b*vlwm**2 - 2184.D0*b*vdw + 480.D0*b*vdt + 
     &    142.D0*b*pi**2 - 3564.D0*vlpm*vlsm + 144.D0*vlpm*vltm + 6984.D
     &    0*vlpm*vlwm + 624.D0*vlpm + 1782.D0*vlpm**2 + 7128.D0*vdmp + 
     &    594.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG**2*xnsq**2 * (  - 24.D0*b*vlsm*
     &    vltm + 120.D0*b*vlsm*vlwm - 24.D0*b*vlsm**2 + 48.D0*b*vltm - 
     &    48.D0*b*vlwm + 120.D0*b*vdw - 24.D0*b*vdt - 8.D0*b*pi**2 - 
     &    180.D0*vlpm**2 - 720.D0*vdmp - 60.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG**2 * (  - 444.D0*b*vlpm**2 + 96.D0
     &    *b*vltm - 456.D0*b*vltm**2 - 96.D0*b*vlwm - 1320.D0*b*vlwm**2
     &     - 1320.D0*b*vdw - 456.D0*b*vdt + 148.D0*b*pi**2 - 1824.D0*
     &    vlpm*vlsm - 744.D0*vlpm*vltm + 4392.D0*vlpm*vlwm + 288.D0*
     &    vlpm + 912.D0*vlpm**2 + 3648.D0*vdmp + 304.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG**3*xnsq * (  - 198.D0*b*vlpm**2 + 
     &    48.D0*b*vltm*vlwm + 24.D0*b*vltm**2 - 768.D0*b*vlwm**2 - 744.D
     &    0*b*vdw + 48.D0*b*vdt + 58.D0*b*pi**2 - 612.D0*vlpm*vlsm + 
     &    1224.D0*vlpm*vlwm + 306.D0*vlpm**2 + 1224.D0*vdmp + 102.D0*
     &    pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG**3*xnsq**2 * (  - 24.D0*b*vlsm*
     &    vltm + 96.D0*b*vlsm*vlwm - 18.D0*b*vlsm**2 + 96.D0*b*vdw - 24.
     &    D0*b*vdt - 6.D0*b*pi**2 - 42.D0*vlpm**2 - 168.D0*vdmp - 14.D0
     &    *pi**2 )
      gcoeff3 = gcoeff3 + t1**4*INVG**3 * (  - 156.D0*b*vlpm**2 - 24.D0
     &    *b*vltm**2 - 600.D0*b*vlwm**2 - 600.D0*b*vdw - 24.D0*b*vdt + 
     &    52.D0*b*pi**2 - 480.D0*vlpm*vlsm - 24.D0*vlpm*vltm + 984.D0*
     &    vlpm*vlwm + 240.D0*vlpm**2 + 960.D0*vdmp + 80.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG*xnsq * ( 1344.D0*vlpm*vlsm - 1344.D
     &    0*vlpm*vltm - 1344.D0*vlpm*vlwm - 1056.D0*vlpm - 672.D0*
     &    vlpm**2 - 2688.D0*vdmp - 224.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG**2*xnsq * ( 348.D0*b*vlpm**2 - 240.
     &    D0*b*vltm*vlwm + 96.D0*b*vltm - 576.D0*b*vltm**2 + 96.D0*b*
     &    vlwm + 1728.D0*b*vlwm**2 + 1608.D0*b*vdw - 696.D0*b*vdt - 76.D
     &    0*b*pi**2 + 4440.D0*vlpm*vlsm - 744.D0*vlpm*vltm - 8136.D0*
     &    vlpm*vlwm - 1056.D0*vlpm - 2220.D0*vlpm**2 - 8880.D0*vdmp - 
     &    740.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG**2*xnsq**2 * ( 96.D0*vlpm**2 + 384.
     &    D0*vdmp + 32.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG**2 * ( 264.D0*b*vlpm**2 + 528.D0*b
     &    *vltm**2 + 528.D0*b*vlwm**2 + 528.D0*b*vdw + 528.D0*b*vdt - 
     &    88.D0*b*pi**2 + 1056.D0*vlpm*vlsm + 1632.D0*vlpm*vltm - 3744.D
     &    0*vlpm*vlwm - 192.D0*vlpm - 528.D0*vlpm**2 - 2112.D0*vdmp - 
     &    176.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG**3*xnsq * ( 372.D0*b*vlpm**2 - 192.
     &    D0*b*vltm*vlwm - 120.D0*b*vltm**2 + 1416.D0*b*vlwm**2 + 1320.D
     &    0*b*vdw - 216.D0*b*vdt - 92.D0*b*pi**2 + 1560.D0*vlpm*vlsm - 
     &    24.D0*vlpm*vltm - 3096.D0*vlpm*vlwm - 780.D0*vlpm**2 - 3120.D0
     &    *vdmp - 260.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG**3*xnsq**2 * ( 72.D0*b*vlsm*vltm
     &     - 120.D0*b*vlsm*vlwm + 12.D0*b*vlsm**2 - 120.D0*b*vdw + 72.D0
     &    *b*vdt + 4.D0*b*pi**2 + 108.D0*vlpm**2 + 432.D0*vdmp + 36.D0*
     &    pi**2 )
      gcoeff3 = gcoeff3 + t1**5*INVG**3 * ( 264.D0*b*vlpm**2 + 144.D0*b
     &    *vltm**2 + 912.D0*b*vlwm**2 + 912.D0*b*vdw + 144.D0*b*vdt - 
     &    88.D0*b*pi**2 + 960.D0*vlpm*vlsm + 192.D0*vlpm*vltm - 2112.D0
     &    *vlpm*vlwm - 480.D0*vlpm**2 - 1920.D0*vdmp - 160.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**6*INVG**2*xnsq * ( 528.D0*b*vltm**2 - 528.
     &    D0*b*vlwm**2 - 528.D0*b*vdw + 528.D0*b*vdt - 3312.D0*vlpm*
     &    vlsm + 1632.D0*vlpm*vltm + 4992.D0*vlpm*vlwm + 960.D0*vlpm + 
     &    1656.D0*vlpm**2 + 6624.D0*vdmp + 552.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**6*INVG**2 * (  - 1248.D0*vlpm*vltm + 1248.
     &    D0*vlpm*vlwm )
      gcoeff3 = gcoeff3 + t1**6*INVG**3*xnsq * (  - 360.D0*b*vlpm**2 + 
     &    240.D0*b*vltm*vlwm + 264.D0*b*vltm**2 - 1464.D0*b*vlwm**2 - 
     &    1344.D0*b*vdw + 384.D0*b*vdt + 80.D0*b*pi**2 - 2400.D0*vlpm*
     &    vlsm + 192.D0*vlpm*vltm + 4608.D0*vlpm*vlwm + 1200.D0*vlpm**2
     &     + 4800.D0*vdmp + 400.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**6*INVG**3*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm + 48.D0*b*vlsm*vlwm + 48.D0*b*vdw - 48.D0*b*vdt - 120.D0
     &    *vlpm**2 - 480.D0*vdmp - 40.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**6*INVG**3 * (  - 240.D0*b*vlpm**2 - 288.D0
     &    *b*vltm**2 - 672.D0*b*vlwm**2 - 672.D0*b*vdw - 288.D0*b*vdt
     &     + 80.D0*b*pi**2 - 960.D0*vlpm*vlsm - 576.D0*vlpm*vltm + 2496.
     &    D0*vlpm*vlwm + 480.D0*vlpm**2 + 1920.D0*vdmp + 160.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**7*INVG**2*xnsq * ( 1248.D0*vlpm*vlsm - 
     &    1248.D0*vlpm*vltm - 1248.D0*vlpm*vlwm - 384.D0*vlpm - 624.D0*
     &    vlpm**2 - 2496.D0*vdmp - 208.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**7*INVG**3*xnsq * ( 144.D0*b*vlpm**2 - 96.D
     &    0*b*vltm*vlwm - 336.D0*b*vltm**2 + 816.D0*b*vlwm**2 + 768.D0*
     &    b*vdw - 384.D0*b*vdt - 32.D0*b*pi**2 + 2304.D0*vlpm*vlsm - 
     &    576.D0*vlpm*vltm - 4032.D0*vlpm*vlwm - 1152.D0*vlpm**2 - 4608.
     &    D0*vdmp - 384.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**7*INVG**3*xnsq**2 * ( 48.D0*vlpm**2 + 192.
     &    D0*vdmp + 16.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**7*INVG**3 * ( 96.D0*b*vlpm**2 + 192.D0*b*
     &    vltm**2 + 192.D0*b*vlwm**2 + 192.D0*b*vdw + 192.D0*b*vdt - 32.
     &    D0*b*pi**2 + 384.D0*vlpm*vlsm + 768.D0*vlpm*vltm - 1536.D0*
     &    vlpm*vlwm - 192.D0*vlpm**2 - 768.D0*vdmp - 64.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**8*INVG**3*xnsq * ( 192.D0*b*vltm**2 - 192.
     &    D0*b*vlwm**2 - 192.D0*b*vdw + 192.D0*b*vdt - 1344.D0*vlpm*
     &    vlsm + 768.D0*vlpm*vltm + 1920.D0*vlpm*vlwm + 672.D0*vlpm**2
     &     + 2688.D0*vdmp + 224.D0*pi**2 )
      gcoeff3 = gcoeff3 + t1**8*INVG**3 * (  - 384.D0*vlpm*vltm + 384.D0
     &    *vlpm*vlwm )
      gcoeff3 = gcoeff3 + t1**9*INVG**3*xnsq * ( 384.D0*vlpm*vlsm - 384.
     &    D0*vlpm*vltm - 384.D0*vlpm*vlwm - 192.D0*vlpm**2 - 768.D0*
     &    vdmp - 64.D0*pi**2 )
      gcoeff3 = gcoeff3 + t2**(-2)*ro*xnsq * ( 12.D0*b*vlwm**2 + 12.D0*
     &    b*vdw + 2.D0*b*pi**2 )
      gcoeff3 = gcoeff3 + t2**(-2)*ro * (  - 12.D0*b*vlwm**2 - 12.D0*b*
     &    vdw - 2.D0*b*pi**2 )
      gcoeff3 = gcoeff3 + t2**(-1)*ro * (  - 12.D0*b*vlwm**2 - 12.D0*b*
     &    vdw - 2.D0*b*pi**2 )
      gcoeff3 = gcoeff3 + t2**(-1)*UBAR**(-1)*xnsq**2 * (  - 12.D0*b*
     &    vlwm )
      gcoeff3 = gcoeff3 + t2**(-1)*UBAR**(-1) * ( 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + t2**(-1)*INVG*xnsq**2 * ( 6.D0*b*vlsm*
     &    omro**(-1) - 6.D0*b*vlsm + 12.D0*b*vlwm - 3.D0/2.D0*vlpm**2*
     &    omro**(-1) + 3.D0/2.D0*vlpm**2 - 6.D0*vdmp*omro**(-1) + 6.D0*
     &    vdmp - 1.D0/2.D0*omro**(-1)*pi**2 + 1.D0/2.D0*pi**2 )
      gcoeff3 = gcoeff3 + t2**(-1)*INVG * (  - 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + t2**(-1)*xnsq * ( 48.D0*b*vlwm - 48.D0*b )
      gcoeff3 = gcoeff3 + t2**(-1)*xnsq**2 * ( 24.D0*b*vlsm*omro**(-1)
     &     - 6.D0*vlpm**2*omro**(-1) - 24.D0*vdmp*omro**(-1) - 2.D0*
     &    omro**(-1)*pi**2 )
      gcoeff3 = gcoeff3 + t2**(-1) * (  - 48.D0*b*vlwm + 48.D0*b )
      gcoeff3 = gcoeff3 + ro*xnsq * ( 12.D0*b*vlpm**2*TR*xn + 12.D0*b*
     &    vltm**2 - 12.D0*b*vlwm**2 - 12.D0*b*vdw + 12.D0*b*vdt - 12.D0
     &    *b*TR*pi**2*xn + 96.D0*b*TR*xn - 48.D0*vlpm*TR*xn - 36.D0*
     &    vlpm )
      gcoeff3 = gcoeff3 + ro**2*xnsq * ( 48.D0*vlpm*TR*xn )
      gcoeff3 = gcoeff3 + UBAR**(-2)*xnsq * (  - 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + UBAR**(-2) * ( 12.D0*b*vlwm )
      gcoeff3 = gcoeff3 + UBAR**(-1)*INVG*xnsq * (  - 36.D0*b*vlwm + 24.
     &    D0*b )
      gcoeff3 = gcoeff3 + UBAR**(-1)*INVG*xnsq**2 * ( 12.D0*b*vlwm - 12.
     &    D0*b )
      gcoeff3 = gcoeff3 + UBAR**(-1)*INVG * ( 24.D0*b*vlwm - 12.D0*b )
      gcoeff3 = gcoeff3 + UBAR**(-1)*xnsq * ( 24.D0*b*vlwm - 12.D0*b )
      gcoeff3 = gcoeff3 + UBAR**(-1)*xnsq**2 * ( 24.D0*b*vlwm + 12.D0*b
     &     )
      gcoeff3 = gcoeff3 + INVG*xnsq * (  - 18.D0*b*vlpm**2 - 72.D0*b*
     &    vlwm**2 - 72.D0*b*vdw + 6.D0*b*pi**2 + 24.D0*b - 48.D0*vlpm*
     &    vlsm + 96.D0*vlpm*vlwm + 12.D0*vlpm + 24.D0*vlpm**2 + 96.D0*
     &    vdmp + 8.D0*pi**2 )
      gcoeff3 = gcoeff3 + INVG*xnsq**2 * (  - 21.D0*b*vlsm*omro**(-1)
     &     + 3.D0*b*vlsm - 24.D0*b*vlwm - 12.D0*b + 21.D0/4.D0*vlpm**2*
     &    omro**(-1) - 3.D0/4.D0*vlpm**2 + 21.D0*vdmp*omro**(-1) - 3.D0
     &    *vdmp + 7.D0/4.D0*omro**(-1)*pi**2 - 1.D0/4.D0*pi**2 )
      gcoeff3 = gcoeff3 + INVG * (  - 18.D0*b*vlpm**2 - 72.D0*b*vlwm**2
     &     - 72.D0*b*vdw + 6.D0*b*pi**2 - 12.D0*b - 48.D0*vlpm*vlsm + 
     &    96.D0*vlpm*vlwm + 12.D0*vlpm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.
     &    D0*pi**2 )
      gcoeff3 = gcoeff3 + INVG**2*xnsq * (  - 3.D0*b*vlpm**2 - 12.D0*b*
     &    vlwm**2 - 12.D0*b*vdw + b*pi**2 - 6.D0*vlpm*vlsm + 12.D0*vlpm
     &    *vlwm + 3.D0*vlpm**2 + 12.D0*vdmp + pi**2 )
      gcoeff3 = gcoeff3 + INVG**2 * (  - 3.D0*b*vlpm**2 - 12.D0*b*
     &    vlwm**2 - 12.D0*b*vdw + b*pi**2 - 6.D0*vlpm*vlsm + 12.D0*vlpm
     &    *vlwm + 3.D0*vlpm**2 + 12.D0*vdmp + pi**2 )
      gcoeff3 = gcoeff3 + xnsq * ( 16.D0*XLF*b*TR*xn - 6.D0*b*vlpm**2
     &     - 48.D0*b*vltm - 72.D0*b*vltm**2 - 48.D0*b*vlwm - 96.D0*b*
     &    vlwm**2 - 96.D0*b*vdw - 72.D0*b*vdt + 16.D0*b*TR*xn - 22.D0*b
     &    *pi**2 + 24.D0*b - 72.D0*vlpm*vlsm + 144.D0*vlpm*vlwm + 108.D0
     &    *vlpm + 36.D0*vlpm**2 + 144.D0*vdmp + 12.D0*pi**2 )
      gcoeff3 = gcoeff3 + xnsq**2 * (  - 48.D0*b*vlsm*omro**(-1) - 8.D0
     &    *b + 12.D0*vlpm**2*omro**(-1) + 9.D0*vlpm**2 + 48.D0*vdmp*
     &    omro**(-1) + 36.D0*vdmp + 4.D0*omro**(-1)*pi**2 + 3.D0*pi**2
     &     )
      gcoeff3 = gcoeff3 + 12.D0*b*vlpm**2 + 48.D0*b*vltm + 72.D0*b*
     &    vltm**2 - 48.D0*b*vlwm - 24.D0*b*vlwm**2 - 24.D0*b*vdw + 72.D0
     &    *b*vdt - 4.D0*b*pi**2 - 72.D0*vlpm*vlsm + 144.D0*vlpm*vlwm + 
     &    48.D0*vlpm + 36.D0*vlpm**2 + 144.D0*vdmp + 12.D0*pi**2

      gcoeff4 =  + t1**(-1)*ro*xnsq * (  - 6.D0*b*vltm**2 - 6.D0*b*vdt
     &     - b*pi**2 )
      gcoeff4 = gcoeff4 + t1**(-1)*ro * ( 6.D0*b*vltm**2 + 6.D0*b*vdt
     &     + b*pi**2 )
      gcoeff4 = gcoeff4 + t1**(-1)*INVG*xnsq**2 * (  - 3.D0/2.D0*b*vlsm
     &    *omro**(-1) + 3.D0/2.D0*b*vlsm + 3.D0/8.D0*vlpm**2*omro**(-1)
     &     - 3.D0/8.D0*vlpm**2 + 3.D0/2.D0*vdmp*omro**(-1) - 3.D0/2.D0*
     &    vdmp + 1.D0/8.D0*omro**(-1)*pi**2 - 1.D0/8.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**(-1)*xnsq**2 * (  - 6.D0*b*vlsm*
     &    omro**(-1) + 3.D0/2.D0*vlpm**2*omro**(-1) + 6.D0*vdmp*
     &    omro**(-1) + 1.D0/2.D0*omro**(-1)*pi**2 )
      gcoeff4 = gcoeff4 + t1*ro*xnsq * (  - 48.D0*vlpm - 12.D0*vlpm**2
     &     + 48.D0*vdmb - 16.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1*UBAR**(-2)*xnsq**2 * ( 24.D0*b*vlwm )
      gcoeff4 = gcoeff4 + t1*UBAR**(-2) * (  - 24.D0*b*vlwm )
      gcoeff4 = gcoeff4 + t1*UBAR**(-1)*xnsq * ( 72.D0*b*vlwm + 24.D0*b
     &     )
      gcoeff4 = gcoeff4 + t1*UBAR**(-1)*xnsq**2 * (  - 168.D0*b*vlwm - 
     &    24.D0*b )
      gcoeff4 = gcoeff4 + t1*INVG*xnsq * (  - 12.D0*b*vlpm**2 + 48.D0*b
     &    *vltm*vlwm - 48.D0*b*vltm + 24.D0*b*vltm**2 + 12.D0*b*vlwm - 
     &    24.D0*b*vlwm**2 + 48.D0*b*vdt - 4.D0*b*pi**2 - 12.D0*vlpm*
     &    vlsm + 24.D0*vlpm*vlwm - 12.D0*vlpm + 6.D0*vlpm**2 + 24.D0*
     &    vdmp + 2.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1*INVG*xnsq**2 * (  - 24.D0*b*vlsm*vltm + 6.D
     &    0*b*vlsm**2 + 36.D0*b*vltm - 24.D0*b*vdt + 2.D0*b*pi**2 - 6.D0
     &    *vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1*INVG * (  - 18.D0*b*vlpm**2 + 12.D0*b*vltm
     &     - 24.D0*b*vltm**2 + 12.D0*b*vlwm - 48.D0*b*vlwm**2 - 48.D0*b
     &    *vdw - 24.D0*b*vdt + 6.D0*b*pi**2 - 24.D0*vlpm*vltm + 24.D0*
     &    vlpm*vlwm - 12.D0*vlpm )
      gcoeff4 = gcoeff4 + t1*INVG**2*xnsq * ( 3.D0*b*vlpm**2 + 12.D0*b*
     &    vlwm**2 + 12.D0*b*vdw - b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vlwm - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff4 = gcoeff4 + t1*INVG**2 * ( 3.D0*b*vlpm**2 + 12.D0*b*
     &    vlwm**2 + 12.D0*b*vdw - b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vlwm - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff4 = gcoeff4 + t1*xnsq * ( 64.D0*XLF*b*TR*rmuom2*xn - 18.D0*
     &    b*vlpm**2 + 24.D0*b*vltm + 120.D0*b*vltm**2 + 24.D0*b*vlwm - 
     &    48.D0*b*vlwm**2 - 48.D0*b*vdw + 120.D0*b*vdt + 30.D0*b*pi**2
     &     - 96.D0*b - 120.D0*vlpm*vlsm + 240.D0*vlpm*vlwm - 24.D0*vlpm
     &     + 84.D0*vlpm**2 + 240.D0*vdmp - 96.D0*vdmb + 52.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1*xnsq**2 * ( 72.D0*b*vlsm*omro**(-1) - 96.D0
     &    *b*vlsm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.D0*b*vdt - 16.D0
     &    *b*pi**2 - 176.D0*b*rmuom2 + 48.D0*b - 18.D0*vlpm**2*
     &    omro**(-1) - 12.D0*vlpm**2 - 72.D0*vdmp*omro**(-1) - 48.D0*
     &    vdmp - 6.D0*omro**(-1)*pi**2 - 4.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1 * (  - 12.D0*b*vlpm**2 - 24.D0*b*vltm**2
     &     - 24.D0*b*vlwm**2 - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0*b*pi**2
     &     - 24.D0*vlpm*vlsm - 144.D0*vlpm*vltm + 192.D0*vlpm*vlwm - 48.
     &    D0*vlpm + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**2*UBAR**(-2)*xnsq * ( 24.D0*b*vlwm )
      gcoeff4 = gcoeff4 + t1**2*UBAR**(-2)*xnsq**2 * (  - 24.D0*b*vlwm
     &     )
      gcoeff4 = gcoeff4 + t1**2*INVG*xnsq * ( 30.D0*b*vlpm**2 - 216.D0*
     &    b*vltm*vlwm + 60.D0*b*vltm - 132.D0*b*vltm**2 - 60.D0*b*vlwm
     &     + 36.D0*b*vlwm**2 - 72.D0*b*vdw - 240.D0*b*vdt + 26.D0*b*
     &    pi**2 + 96.D0*vlpm*vlsm - 24.D0*vlpm*vltm - 168.D0*vlpm*vlwm
     &     + 48.D0*vlpm - 48.D0*vlpm**2 - 192.D0*vdmp - 16.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**2*INVG*xnsq**2 * ( 120.D0*b*vlsm*vltm + 
     &    120.D0*b*vlsm*vlwm + 48.D0*b*vlsm - 60.D0*b*vlsm**2 - 84.D0*b
     &    *vltm - 12.D0*b*vlwm + 120.D0*b*vdw + 120.D0*b*vdt - 20.D0*b*
     &    pi**2 + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**2*INVG * ( 48.D0*b*vlpm**2 + 96.D0*b*
     &    vltm**2 + 96.D0*b*vlwm**2 + 96.D0*b*vdw + 96.D0*b*vdt - 16.D0
     &    *b*pi**2 + 144.D0*vlpm*vltm - 144.D0*vlpm*vlwm + 48.D0*vlpm )
      gcoeff4 = gcoeff4 + t1**2*INVG**2*xnsq * (  - 15.D0*b*vlpm**2 - 
     &    60.D0*b*vlwm**2 - 60.D0*b*vdw + 5.D0*b*pi**2 - 42.D0*vlpm*
     &    vlsm + 84.D0*vlpm*vlwm + 21.D0*vlpm**2 + 84.D0*vdmp + 7.D0*
     &    pi**2 )
      gcoeff4 = gcoeff4 + t1**2*INVG**2*xnsq**2 * ( 12.D0*b*vlsm*vlwm
     &     - 3.D0*b*vlsm**2 + 12.D0*b*vdw - b*pi**2 - 3.D0*vlpm**2 - 12.
     &    D0*vdmp - pi**2 )
      gcoeff4 = gcoeff4 + t1**2*INVG**2 * (  - 12.D0*b*vlpm**2 - 48.D0*
     &    b*vlwm**2 - 48.D0*b*vdw + 4.D0*b*pi**2 - 36.D0*vlpm*vlsm + 72.
     &    D0*vlpm*vlwm + 18.D0*vlpm**2 + 72.D0*vdmp + 6.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**2*xnsq * (  - 24.D0*b*vltm**2 + 24.D0*b*
     &    vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt + 192.D0*vlpm*vlsm - 144.D
     &    0*vlpm*vltm - 240.D0*vlpm*vlwm - 96.D0*vlpm**2 - 384.D0*vdmp
     &     - 32.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**2 * ( 48.D0*vlpm*vltm - 48.D0*vlpm*vlwm )
      gcoeff4 = gcoeff4 + t1**3*INVG*xnsq * (  - 6.D0*b*vlpm**2 + 120.D0
     &    *b*vltm*vlwm + 156.D0*b*vltm**2 - 60.D0*b*vlwm**2 + 216.D0*b*
     &    vdt - 18.D0*b*pi**2 - 240.D0*vlpm*vlsm + 144.D0*vlpm*vltm + 
     &    336.D0*vlpm*vlwm - 48.D0*vlpm + 120.D0*vlpm**2 + 480.D0*vdmp
     &     + 40.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**3*INVG*xnsq**2 * (  - 96.D0*b*vlsm*vltm
     &     - 96.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 48.D0*b*vlsm**2 + 48.D0
     &    *b*vltm + 48.D0*b*vlwm - 96.D0*b*vdw - 96.D0*b*vdt + 16.D0*b*
     &    pi**2 + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**3*INVG * (  - 12.D0*b*vlpm**2 - 24.D0*b*
     &    vltm**2 - 24.D0*b*vlwm**2 - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0*
     &    b*pi**2 + 24.D0*vlpm*vlsm - 240.D0*vlpm*vltm + 192.D0*vlpm*
     &    vlwm - 48.D0*vlpm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**3*INVG**2*xnsq * ( 27.D0*b*vlpm**2 + 48.D0
     &    *b*vltm*vlwm + 24.D0*b*vltm**2 + 132.D0*b*vlwm**2 + 156.D0*b*
     &    vdw + 48.D0*b*vdt - 17.D0*b*pi**2 + 114.D0*vlpm*vlsm - 228.D0
     &    *vlpm*vlwm - 57.D0*vlpm**2 - 228.D0*vdmp - 19.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**3*INVG**2*xnsq**2 * (  - 36.D0*b*vlsm*
     &    vltm - 72.D0*b*vlsm*vlwm + 27.D0*b*vlsm**2 - 72.D0*b*vdw - 36.
     &    D0*b*vdt + 9.D0*b*pi**2 + 15.D0*vlpm**2 + 60.D0*vdmp + 5.D0*
     &    pi**2 )
      gcoeff4 = gcoeff4 + t1**3*INVG**2 * ( 12.D0*b*vlpm**2 - 12.D0*b*
     &    vltm**2 + 60.D0*b*vlwm**2 + 60.D0*b*vdw - 12.D0*b*vdt - 4.D0*
     &    b*pi**2 + 84.D0*vlpm*vlsm - 12.D0*vlpm*vltm - 156.D0*vlpm*
     &    vlwm - 42.D0*vlpm**2 - 168.D0*vdmp - 14.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**3*xnsq * (  - 48.D0*vlpm*vlsm + 48.D0*
     &    vlpm*vltm + 48.D0*vlpm*vlwm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.D
     &    0*pi**2 )
      gcoeff4 = gcoeff4 + t1**4*INVG*xnsq * (  - 24.D0*b*vltm**2 + 24.D0
     &    *b*vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt + 240.D0*vlpm*vlsm - 
     &    240.D0*vlpm*vltm - 240.D0*vlpm*vlwm - 120.D0*vlpm**2 - 480.D0
     &    *vdmp - 40.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**4*INVG * ( 48.D0*vlpm*vltm - 48.D0*vlpm*
     &    vlwm )
      gcoeff4 = gcoeff4 + t1**4*INVG**2*xnsq * (  - 24.D0*b*vlpm**2 - 
     &    96.D0*b*vltm*vlwm - 60.D0*b*vltm**2 - 132.D0*b*vlwm**2 - 180.D
     &    0*b*vdw - 108.D0*b*vdt + 24.D0*b*pi**2 - 144.D0*vlpm*vlsm - 
     &    12.D0*vlpm*vltm + 300.D0*vlpm*vlwm + 72.D0*vlpm**2 + 288.D0*
     &    vdmp + 24.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**4*INVG**2*xnsq**2 * ( 84.D0*b*vlsm*vltm
     &     + 108.D0*b*vlsm*vlwm - 48.D0*b*vlsm**2 + 108.D0*b*vdw + 84.D0
     &    *b*vdt - 16.D0*b*pi**2 - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D0*
     &    pi**2 )
      gcoeff4 = gcoeff4 + t1**4*INVG**2 * ( 24.D0*b*vltm**2 - 24.D0*b*
     &    vlwm**2 - 24.D0*b*vdw + 24.D0*b*vdt - 96.D0*vlpm*vlsm + 48.D0
     &    *vlpm*vltm + 144.D0*vlpm*vlwm + 48.D0*vlpm**2 + 192.D0*vdmp
     &     + 16.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**5*INVG*xnsq * (  - 48.D0*vlpm*vlsm + 48.D0
     &    *vlpm*vltm + 48.D0*vlpm*vlwm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.
     &    D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**5*INVG**2*xnsq * ( 12.D0*b*vlpm**2 + 48.D0
     &    *b*vltm*vlwm + 48.D0*b*vltm**2 + 48.D0*b*vlwm**2 + 72.D0*b*
     &    vdw + 72.D0*b*vdt - 12.D0*b*pi**2 + 72.D0*vlpm*vlsm + 48.D0*
     &    vlpm*vltm - 192.D0*vlpm*vlwm - 36.D0*vlpm**2 - 144.D0*vdmp - 
     &    12.D0*pi**2 )
      gcoeff4 = gcoeff4 + t1**5*INVG**2*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm - 48.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.
     &    D0*b*vdt + 8.D0*b*pi**2 + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*
     &    pi**2 )
      gcoeff4 = gcoeff4 + t1**5*INVG**2 * ( 48.D0*vlpm*vlsm - 48.D0*
     &    vlpm*vltm - 48.D0*vlpm*vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D
     &    0*pi**2 )
      gcoeff4 = gcoeff4 + t1**6*INVG**2*xnsq * (  - 48.D0*vlpm*vltm + 
     &    48.D0*vlpm*vlwm )
      gcoeff4 = gcoeff4 + t2**(-2)*ro*xnsq * (  - 6.D0*b*vlwm**2 - 6.D0
     &    *b*vdw - b*pi**2 )
      gcoeff4 = gcoeff4 + t2**(-2)*ro * ( 6.D0*b*vlwm**2 + 6.D0*b*vdw
     &     + b*pi**2 )
      gcoeff4 = gcoeff4 + t2**(-1)*ro*xnsq * ( 12.D0*b*vlwm**2 + 12.D0*
     &    b*vdw + 2.D0*b*pi**2 )
      gcoeff4 = gcoeff4 + t2**(-1)*ro * (  - 6.D0*b*vlwm**2 - 6.D0*b*
     &    vdw - b*pi**2 - 12.D0*vlpm*vlsm + 24.D0*vlpm*vlwm + 24.D0*
     &    vdmp + 24.D0*vdmb - 6.D0*pi**2 )
      gcoeff4 = gcoeff4 + t2**(-1)*INVG*xnsq**2 * (  - 3.D0/2.D0*b*vlsm
     &    *omro**(-1) + 3.D0/2.D0*b*vlsm + 3.D0/8.D0*vlpm**2*omro**(-1)
     &     - 3.D0/8.D0*vlpm**2 + 3.D0/2.D0*vdmp*omro**(-1) - 3.D0/2.D0*
     &    vdmp + 1.D0/8.D0*omro**(-1)*pi**2 - 1.D0/8.D0*pi**2 )
      gcoeff4 = gcoeff4 + t2**(-1)*xnsq * (  - 32.D0*XLF*b*TR*rmuom2*xn
     &     + 96.D0*b*vltm*vlwm - 216.D0*b*vlwm - 24.D0*b*pi**2 - 88.D0*
     &    b*rmuom2 + 24.D0*b )
      gcoeff4 = gcoeff4 + t2**(-1)*xnsq**2 * (  - 48.D0*b*vlsm*vlwm - 6.
     &    D0*b*vlsm*omro**(-1) + 144.D0*b*vlwm + 12.D0*b*pi**2 + 88.D0*
     &    b*rmuom2 + 3.D0/2.D0*vlpm**2*omro**(-1) + 6.D0*vdmp*
     &    omro**(-1) + 1.D0/2.D0*omro**(-1)*pi**2 )
      gcoeff4 = gcoeff4 + t2**(-1) * ( 32.D0*XLF*b*TR*rmuom2*xn + 72.D0
     &    *b*vlwm - 24.D0*b + 24.D0*vlpm*vlsm - 48.D0*vlpm*vlwm - 48.D0
     &    *vdmp - 48.D0*vdmb + 12.D0*pi**2 )
      gcoeff4 = gcoeff4 + ro*xnsq * ( 6.D0*b*vltm**2 - 6.D0*b*vlwm**2
     &     - 6.D0*b*vdw + 6.D0*b*vdt - 12.D0*vlpm*vlsm + 24.D0*vlpm*
     &    vlwm + 24.D0*vlpm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff4 = gcoeff4 + UBAR**(-2)*xnsq * (  - 24.D0*b*vlwm )
      gcoeff4 = gcoeff4 + UBAR**(-2) * ( 24.D0*b*vlwm )
      gcoeff4 = gcoeff4 + UBAR**(-1)*xnsq * ( 168.D0*b*vlwm + 24.D0*b )
      gcoeff4 = gcoeff4 + UBAR**(-1) * (  - 72.D0*b*vlwm - 24.D0*b )
      gcoeff4 = gcoeff4 + INVG*xnsq**2 * ( 6.D0*b*vlsm*omro**(-1) - 3.D0
     &    /2.D0*vlpm**2*omro**(-1) - 6.D0*vdmp*omro**(-1) - 1.D0/2.D0*
     &    omro**(-1)*pi**2 )
      gcoeff4 = gcoeff4 + xnsq * ( 24.D0*b*vlpm**2 - 96.D0*b*vltm*vlwm
     &     - 24.D0*b*vltm - 48.D0*b*vltm**2 + 72.D0*b*vlwm + 48.D0*b*
     &    vlwm**2 - 96.D0*b*vdt + 8.D0*b*pi**2 + 24.D0*b + 48.D0*vlpm*
     &    vlsm - 96.D0*vlpm*vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D0*
     &    pi**2 )
      gcoeff4 = gcoeff4 + xnsq**2 * ( 48.D0*b*vlsm*vlwm - 48.D0*b*vlsm*
     &    omro**(-1) + 48.D0*b*vlsm - 12.D0*b*vlsm**2 - 144.D0*b*vlwm
     &     + 48.D0*b*vdt - 4.D0*b*pi**2 - 24.D0*b + 12.D0*vlpm**2*
     &    omro**(-1) + 48.D0*vdmp*omro**(-1) + 4.D0*omro**(-1)*pi**2 )
      gcoeff4 = gcoeff4 + 36.D0*b*vlpm**2 + 24.D0*b*vltm + 48.D0*b*
     &    vltm**2 - 24.D0*b*vlwm + 96.D0*b*vlwm**2 + 96.D0*b*vdw + 48.D0
     &    *b*vdt - 12.D0*b*pi**2 + 24.D0*vlpm*vlsm - 48.D0*vlpm*vlwm - 
     &    12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2

      gcoeff5 =  + t1**(-2)*ro*xnsq * (  - 6.D0*b*vltm**2 - 6.D0*b*vdt
     &     - b*pi**2 )
      gcoeff5 = gcoeff5 + t1**(-2)*ro * ( 6.D0*b*vltm**2 + 6.D0*b*vdt
     &     + b*pi**2 )
      gcoeff5 = gcoeff5 + t1**(-1)*ro * ( 6.D0*b*vltm**2 + 6.D0*b*vdt
     &     + b*pi**2 )
      gcoeff5 = gcoeff5 + t1**(-1)*UBAR**(-1)*INVG*xnsq**2 * (  - 24.D0
     &    *b*vlwm )
      gcoeff5 = gcoeff5 + t1**(-1)*UBAR**(-1)*INVG * ( 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t1**(-1)*UBAR**(-1)*xnsq**2 * (  - 24.D0*b*
     &    vlwm )
      gcoeff5 = gcoeff5 + t1**(-1)*UBAR**(-1) * ( 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t1**(-1)*INVG*xnsq**2 * ( 3.D0/2.D0*b*vlsm*
     &    omro**(-1) - 3.D0/2.D0*b*vlsm + 24.D0*b*vlwm - 3.D0/8.D0*
     &    vlpm**2*omro**(-1) + 3.D0/8.D0*vlpm**2 - 3.D0/2.D0*vdmp*
     &    omro**(-1) + 3.D0/2.D0*vdmp - 1.D0/8.D0*omro**(-1)*pi**2 + 1.D
     &    0/8.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**(-1)*INVG * (  - 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t1**(-1)*xnsq * (  - 24.D0*b*vltm + 24.D0*b )
      gcoeff5 = gcoeff5 + t1**(-1)*xnsq**2 * ( 6.D0*b*vlsm*omro**(-1)
     &     - 3.D0/2.D0*vlpm**2*omro**(-1) - 6.D0*vdmp*omro**(-1) - 1.D0/
     &    2.D0*omro**(-1)*pi**2 )
      gcoeff5 = gcoeff5 + t1**(-1) * ( 24.D0*b*vltm - 24.D0*b )
      gcoeff5 = gcoeff5 + t1*ro*xnsq * (  - 48.D0*vlpm - 12.D0*vlpm**2
     &     + 48.D0*vdmb - 16.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1*UBAR**(-2)*xnsq**2 * ( 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t1*UBAR**(-2) * (  - 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t1*UBAR**(-1)*xnsq * ( 72.D0*b*vlwm + 12.D0*b
     &     )
      gcoeff5 = gcoeff5 + t1*UBAR**(-1)*xnsq**2 * (  - 120.D0*b*vlwm - 
     &    12.D0*b )
      gcoeff5 = gcoeff5 + t1*INVG*xnsq * ( 30.D0*b*vlpm**2 + 24.D0*b*
     &    vltm*vlwm - 108.D0*b*vltm + 12.D0*b*vltm**2 + 72.D0*b*vlwm + 
     &    132.D0*b*vlwm**2 + 144.D0*b*vdw + 24.D0*b*vdt - 14.D0*b*pi**2
     &     + 156.D0*vlpm*vlsm - 312.D0*vlpm*vlwm - 60.D0*vlpm - 78.D0*
     &    vlpm**2 - 312.D0*vdmp - 26.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1*INVG*xnsq**2 * (  - 24.D0*b*vlsm*vltm - 48.
     &    D0*b*vlsm + 6.D0*b*vlsm**2 + 120.D0*b*vltm + 12.D0*b*vlwm - 
     &    24.D0*b*vdt + 2.D0*b*pi**2 )
      gcoeff5 = gcoeff5 + t1*INVG * ( 18.D0*b*vlpm**2 + 12.D0*b*vltm + 
     &    12.D0*b*vlwm + 72.D0*b*vlwm**2 + 72.D0*b*vdw - 6.D0*b*pi**2
     &     + 120.D0*vlpm*vlsm - 240.D0*vlpm*vlwm - 60.D0*vlpm - 60.D0*
     &    vlpm**2 - 240.D0*vdmp - 20.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1*INVG**2*xnsq * ( 18.D0*b*vlpm**2 + 72.D0*b
     &    *vlwm**2 + 72.D0*b*vdw - 6.D0*b*pi**2 + 48.D0*vlpm*vlsm - 96.D
     &    0*vlpm*vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1*INVG**2*xnsq**2 * (  - 12.D0*b*vlsm*vlwm
     &     + 3.D0*b*vlsm**2 - 12.D0*b*vdw + b*pi**2 + 3.D0*vlpm**2 + 12.
     &    D0*vdmp + pi**2 )
      gcoeff5 = gcoeff5 + t1*INVG**2 * ( 15.D0*b*vlpm**2 + 60.D0*b*
     &    vlwm**2 + 60.D0*b*vdw - 5.D0*b*pi**2 + 42.D0*vlpm*vlsm - 84.D0
     &    *vlpm*vlwm - 21.D0*vlpm**2 - 84.D0*vdmp - 7.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1*xnsq * ( 64.D0*XLF*b*TR*rmuom2*xn - 18.D0*
     &    b*vlpm**2 + 24.D0*b*vltm + 48.D0*b*vltm**2 + 24.D0*b*vlwm + 
     &    24.D0*b*vlwm**2 + 24.D0*b*vdw + 48.D0*b*vdt + 30.D0*b*pi**2
     &     - 96.D0*b + 24.D0*vlpm*vlsm - 48.D0*vlpm*vlwm - 24.D0*vlpm
     &     + 12.D0*vlpm**2 - 48.D0*vdmp - 96.D0*vdmb + 28.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1*xnsq**2 * ( 72.D0*b*vlsm*omro**(-1) - 96.D0
     &    *b*vlsm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.D0*b*vdt - 16.D0
     &    *b*pi**2 - 176.D0*b*rmuom2 + 48.D0*b - 18.D0*vlpm**2*
     &    omro**(-1) - 12.D0*vlpm**2 - 72.D0*vdmp*omro**(-1) - 48.D0*
     &    vdmp - 6.D0*omro**(-1)*pi**2 - 4.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1 * (  - 12.D0*b*vlpm**2 - 24.D0*b*vltm**2
     &     - 24.D0*b*vlwm**2 - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0*b*pi**2
     &     - 24.D0*vlpm*vlsm + 48.D0*vlpm*vlwm - 48.D0*vlpm + 12.D0*
     &    vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**2*UBAR**(-2)*xnsq * ( 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t1**2*UBAR**(-2)*xnsq**2 * (  - 12.D0*b*vlwm
     &     )
      gcoeff5 = gcoeff5 + t1**2*INVG*xnsq * (  - 12.D0*b*vlpm**2 - 192.D
     &    0*b*vltm*vlwm + 60.D0*b*vltm - 96.D0*b*vltm**2 - 60.D0*b*vlwm
     &     - 144.D0*b*vlwm**2 - 240.D0*b*vdw - 192.D0*b*vdt + 36.D0*b*
     &    pi**2 - 216.D0*vlpm*vlsm + 432.D0*vlpm*vlwm + 96.D0*vlpm + 
     &    108.D0*vlpm**2 + 432.D0*vdmp + 36.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**2*INVG*xnsq**2 * ( 120.D0*b*vlsm*vltm + 
     &    120.D0*b*vlsm*vlwm + 96.D0*b*vlsm - 60.D0*b*vlsm**2 - 132.D0*
     &    b*vltm - 60.D0*b*vlwm + 120.D0*b*vdw + 120.D0*b*vdt - 20.D0*b
     &    *pi**2 )
      gcoeff5 = gcoeff5 + t1**2*INVG * ( 12.D0*b*vlpm**2 + 24.D0*b*
     &    vltm**2 + 24.D0*b*vlwm**2 + 24.D0*b*vdw + 24.D0*b*vdt - 4.D0*
     &    b*pi**2 - 120.D0*vlpm*vlsm + 48.D0*vlpm*vltm + 192.D0*vlpm*
     &    vlwm + 96.D0*vlpm + 60.D0*vlpm**2 + 240.D0*vdmp + 20.D0*pi**2
     &     )
      gcoeff5 = gcoeff5 + t1**2*INVG**2*xnsq * (  - 42.D0*b*vlpm**2 - 
     &    48.D0*b*vltm*vlwm - 24.D0*b*vltm**2 - 192.D0*b*vlwm**2 - 216.D
     &    0*b*vdw - 48.D0*b*vdt + 22.D0*b*pi**2 - 156.D0*vlpm*vlsm + 
     &    312.D0*vlpm*vlwm + 78.D0*vlpm**2 + 312.D0*vdmp + 26.D0*pi**2
     &     )
      gcoeff5 = gcoeff5 + t1**2*INVG**2*xnsq**2 * ( 36.D0*b*vlsm*vltm
     &     + 84.D0*b*vlsm*vlwm - 30.D0*b*vlsm**2 + 84.D0*b*vdw + 36.D0*
     &    b*vdt - 10.D0*b*pi**2 - 18.D0*vlpm**2 - 72.D0*vdmp - 6.D0*
     &    pi**2 )
      gcoeff5 = gcoeff5 + t1**2*INVG**2 * (  - 24.D0*b*vlpm**2 + 12.D0*
     &    b*vltm**2 - 108.D0*b*vlwm**2 - 108.D0*b*vdw + 12.D0*b*vdt + 8.
     &    D0*b*pi**2 - 120.D0*vlpm*vlsm + 12.D0*vlpm*vltm + 228.D0*vlpm
     &    *vlwm + 60.D0*vlpm**2 + 240.D0*vdmp + 20.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**2*xnsq * (  - 24.D0*b*vltm**2 + 24.D0*b*
     &    vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt + 48.D0*vlpm*vlsm - 96.D0
     &    *vlpm*vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**2 * ( 48.D0*vlpm*vltm - 48.D0*vlpm*vlwm )
      gcoeff5 = gcoeff5 + t1**3*INVG*xnsq * (  - 6.D0*b*vlpm**2 + 120.D0
     &    *b*vltm*vlwm + 84.D0*b*vltm**2 + 12.D0*b*vlwm**2 + 72.D0*b*
     &    vdw + 144.D0*b*vdt - 18.D0*b*pi**2 + 48.D0*vlpm*vlsm + 48.D0*
     &    vlpm*vltm - 144.D0*vlpm*vlwm - 48.D0*vlpm - 24.D0*vlpm**2 - 
     &    96.D0*vdmp - 8.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**3*INVG*xnsq**2 * (  - 96.D0*b*vlsm*vltm
     &     - 96.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 48.D0*b*vlsm**2 + 48.D0
     &    *b*vltm + 48.D0*b*vlwm - 96.D0*b*vdw - 96.D0*b*vdt + 16.D0*b*
     &    pi**2 + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**3*INVG * (  - 12.D0*b*vlpm**2 - 24.D0*b*
     &    vltm**2 - 24.D0*b*vlwm**2 - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0*
     &    b*pi**2 + 24.D0*vlpm*vlsm - 96.D0*vlpm*vltm + 48.D0*vlpm*vlwm
     &     - 48.D0*vlpm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**3*INVG**2*xnsq * ( 51.D0*b*vlpm**2 + 144.D
     &    0*b*vltm*vlwm + 84.D0*b*vltm**2 + 264.D0*b*vlwm**2 + 336.D0*b
     &    *vdw + 156.D0*b*vdt - 41.D0*b*pi**2 + 258.D0*vlpm*vlsm + 12.D0
     &    *vlpm*vltm - 528.D0*vlpm*vlwm - 129.D0*vlpm**2 - 516.D0*vdmp
     &     - 43.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**3*INVG**2*xnsq**2 * (  - 120.D0*b*vlsm*
     &    vltm - 180.D0*b*vlsm*vlwm + 75.D0*b*vlsm**2 - 180.D0*b*vdw - 
     &    120.D0*b*vdt + 25.D0*b*pi**2 + 39.D0*vlpm**2 + 156.D0*vdmp + 
     &    13.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**3*INVG**2 * ( 12.D0*b*vlpm**2 - 36.D0*b*
     &    vltm**2 + 84.D0*b*vlwm**2 + 84.D0*b*vdw - 36.D0*b*vdt - 4.D0*
     &    b*pi**2 + 180.D0*vlpm*vlsm - 60.D0*vlpm*vltm - 300.D0*vlpm*
     &    vlwm - 90.D0*vlpm**2 - 360.D0*vdmp - 30.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**3*xnsq * (  - 48.D0*vlpm*vlsm + 48.D0*
     &    vlpm*vltm + 48.D0*vlpm*vlwm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.D
     &    0*pi**2 )
      gcoeff5 = gcoeff5 + t1**4*INVG*xnsq * (  - 24.D0*b*vltm**2 + 24.D0
     &    *b*vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt + 96.D0*vlpm*vlsm - 96.
     &    D0*vlpm*vltm - 96.D0*vlpm*vlwm - 48.D0*vlpm**2 - 192.D0*vdmp
     &     - 16.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**4*INVG * ( 48.D0*vlpm*vltm - 48.D0*vlpm*
     &    vlwm )
      gcoeff5 = gcoeff5 + t1**4*INVG**2*xnsq * (  - 36.D0*b*vlpm**2 - 
     &    144.D0*b*vltm*vlwm - 108.D0*b*vltm**2 - 180.D0*b*vlwm**2 - 
     &    252.D0*b*vdw - 180.D0*b*vdt + 36.D0*b*pi**2 - 216.D0*vlpm*
     &    vlsm - 60.D0*vlpm*vltm + 492.D0*vlpm*vlwm + 108.D0*vlpm**2 + 
     &    432.D0*vdmp + 36.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**4*INVG**2*xnsq**2 * ( 132.D0*b*vlsm*vltm
     &     + 156.D0*b*vlsm*vlwm - 72.D0*b*vlsm**2 + 156.D0*b*vdw + 132.D
     &    0*b*vdt - 24.D0*b*pi**2 - 36.D0*vlpm**2 - 144.D0*vdmp - 12.D0
     &    *pi**2 )
      gcoeff5 = gcoeff5 + t1**4*INVG**2 * ( 24.D0*b*vltm**2 - 24.D0*b*
     &    vlwm**2 - 24.D0*b*vdw + 24.D0*b*vdt - 144.D0*vlpm*vlsm + 96.D0
     &    *vlpm*vltm + 192.D0*vlpm*vlwm + 72.D0*vlpm**2 + 288.D0*vdmp
     &     + 24.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**5*INVG*xnsq * (  - 48.D0*vlpm*vlsm + 48.D0
     &    *vlpm*vltm + 48.D0*vlpm*vlwm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.
     &    D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**5*INVG**2*xnsq * ( 12.D0*b*vlpm**2 + 48.D0
     &    *b*vltm*vlwm + 48.D0*b*vltm**2 + 48.D0*b*vlwm**2 + 72.D0*b*
     &    vdw + 72.D0*b*vdt - 12.D0*b*pi**2 + 72.D0*vlpm*vlsm + 96.D0*
     &    vlpm*vltm - 240.D0*vlpm*vlwm - 36.D0*vlpm**2 - 144.D0*vdmp - 
     &    12.D0*pi**2 )
      gcoeff5 = gcoeff5 + t1**5*INVG**2*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm - 48.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.
     &    D0*b*vdt + 8.D0*b*pi**2 + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*
     &    pi**2 )
      gcoeff5 = gcoeff5 + t1**5*INVG**2 * ( 48.D0*vlpm*vlsm - 48.D0*
     &    vlpm*vltm - 48.D0*vlpm*vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D
     &    0*pi**2 )
      gcoeff5 = gcoeff5 + t1**6*INVG**2*xnsq * (  - 48.D0*vlpm*vltm + 
     &    48.D0*vlpm*vlwm )
      gcoeff5 = gcoeff5 + t2**(-2)*ro*xnsq * (  - 12.D0*b*vlwm**2 - 12.D
     &    0*b*vdw - 2.D0*b*pi**2 )
      gcoeff5 = gcoeff5 + t2**(-2)*ro * ( 12.D0*b*vlwm**2 + 12.D0*b*vdw
     &     + 2.D0*b*pi**2 )
      gcoeff5 = gcoeff5 + t2**(-1)*ro*xnsq * ( 18.D0*b*vlwm**2 + 18.D0*
     &    b*vdw + 3.D0*b*pi**2 )
      gcoeff5 = gcoeff5 + t2**(-1)*ro * (  - 6.D0*b*vlwm**2 - 6.D0*b*
     &    vdw - b*pi**2 - 12.D0*vlpm*vlsm + 24.D0*vlpm*vlwm + 24.D0*
     &    vdmp + 24.D0*vdmb - 6.D0*pi**2 )
      gcoeff5 = gcoeff5 + t2**(-1)*UBAR**(-1)*xnsq**2 * (  - 24.D0*b*
     &    vlwm )
      gcoeff5 = gcoeff5 + t2**(-1)*UBAR**(-1) * ( 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t2**(-1)*INVG*xnsq**2 * ( 3.D0/2.D0*b*vlsm*
     &    omro**(-1) - 3.D0/2.D0*b*vlsm + 24.D0*b*vlwm - 3.D0/8.D0*
     &    vlpm**2*omro**(-1) + 3.D0/8.D0*vlpm**2 - 3.D0/2.D0*vdmp*
     &    omro**(-1) + 3.D0/2.D0*vdmp - 1.D0/8.D0*omro**(-1)*pi**2 + 1.D
     &    0/8.D0*pi**2 )
      gcoeff5 = gcoeff5 + t2**(-1)*INVG * (  - 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + t2**(-1)*xnsq * (  - 32.D0*XLF*b*TR*rmuom2*xn
     &     + 96.D0*b*vltm*vlwm - 240.D0*b*vlwm - 24.D0*b*pi**2 - 88.D0*
     &    b*rmuom2 + 48.D0*b )
      gcoeff5 = gcoeff5 + t2**(-1)*xnsq**2 * (  - 48.D0*b*vlsm*vlwm + 6.
     &    D0*b*vlsm*omro**(-1) + 144.D0*b*vlwm + 12.D0*b*pi**2 + 88.D0*
     &    b*rmuom2 - 3.D0/2.D0*vlpm**2*omro**(-1) - 6.D0*vdmp*
     &    omro**(-1) - 1.D0/2.D0*omro**(-1)*pi**2 )
      gcoeff5 = gcoeff5 + t2**(-1) * ( 32.D0*XLF*b*TR*rmuom2*xn + 96.D0
     &    *b*vlwm - 48.D0*b + 24.D0*vlpm*vlsm - 48.D0*vlpm*vlwm - 48.D0
     &    *vdmp - 48.D0*vdmb + 12.D0*pi**2 )
      gcoeff5 = gcoeff5 + ro*TBAR**(-2)*xnsq * ( 6.D0*b*vltm )
      gcoeff5 = gcoeff5 + ro*TBAR**(-2)*xnsq**2 * (  - 3.D0*b*vltm )
      gcoeff5 = gcoeff5 + ro*TBAR**(-2) * (  - 3.D0*b*vltm )
      gcoeff5 = gcoeff5 + ro*TBAR**(-1)*xnsq * (  - 6.D0*b*vltm + 3.D0*
     &    b )
      gcoeff5 = gcoeff5 + ro*TBAR**(-1)*xnsq**2 * (  - 6.D0*b*vltm - 3.D
     &    0*b )
      gcoeff5 = gcoeff5 + ro*xnsq * ( 6.D0*b*vltm**2 - 6.D0*b*vlwm**2
     &     - 6.D0*b*vdw + 6.D0*b*vdt - 12.D0*vlpm*vlsm + 24.D0*vlpm*
     &    vlwm + 24.D0*vlpm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff5 = gcoeff5 + ro**2*TBAR**(-2)*xnsq * (  - 3.D0/4.D0*b*vltm
     &     )
      gcoeff5 = gcoeff5 + ro**2*TBAR**(-2)*xnsq**2 * ( 3.D0/4.D0*b*vltm
     &     )
      gcoeff5 = gcoeff5 + TBAR**(-1)*xnsq * (  - 24.D0*b*vltm - 24.D0*b
     &     )
      gcoeff5 = gcoeff5 + TBAR**(-1)*xnsq**2 * ( 36.D0*b*vltm + 12.D0*b
     &     )
      gcoeff5 = gcoeff5 + TBAR**(-1) * (  - 12.D0*b*vltm + 12.D0*b )
      gcoeff5 = gcoeff5 + UBAR**(-2)*xnsq * (  - 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + UBAR**(-2) * ( 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + UBAR**(-1)*INVG*xnsq**2 * ( 24.D0*b*vlwm )
      gcoeff5 = gcoeff5 + UBAR**(-1)*INVG * (  - 12.D0*b*vlwm )
      gcoeff5 = gcoeff5 + UBAR**(-1)*xnsq * ( 120.D0*b*vlwm + 12.D0*b )
      gcoeff5 = gcoeff5 + UBAR**(-1) * (  - 72.D0*b*vlwm - 12.D0*b )
      gcoeff5 = gcoeff5 + INVG*xnsq * (  - 12.D0*b*vlpm**2 + 48.D0*b*
     &    vltm - 12.D0*b*vlwm - 48.D0*b*vlwm**2 - 48.D0*b*vdw + 4.D0*b*
     &    pi**2 - 36.D0*vlpm*vlsm + 72.D0*vlpm*vlwm + 12.D0*vlpm + 18.D0
     &    *vlpm**2 + 72.D0*vdmp + 6.D0*pi**2 )
      gcoeff5 = gcoeff5 + INVG*xnsq**2 * (  - 6.D0*b*vlsm*omro**(-1) - 
     &    36.D0*b*vltm + 3.D0/2.D0*vlpm**2*omro**(-1) + 6.D0*vdmp*
     &    omro**(-1) + 1.D0/2.D0*omro**(-1)*pi**2 )
      gcoeff5 = gcoeff5 + INVG * (  - 12.D0*b*vlpm**2 - 12.D0*b*vltm - 
     &    12.D0*b*vlwm - 48.D0*b*vlwm**2 - 48.D0*b*vdw + 4.D0*b*pi**2
     &     - 36.D0*vlpm*vlsm + 72.D0*vlpm*vlwm + 12.D0*vlpm + 18.D0*
     &    vlpm**2 + 72.D0*vdmp + 6.D0*pi**2 )
      gcoeff5 = gcoeff5 + INVG**2*xnsq * (  - 3.D0*b*vlpm**2 - 12.D0*b*
     &    vlwm**2 - 12.D0*b*vdw + b*pi**2 - 6.D0*vlpm*vlsm + 12.D0*vlpm
     &    *vlwm + 3.D0*vlpm**2 + 12.D0*vdmp + pi**2 )
      gcoeff5 = gcoeff5 + INVG**2 * (  - 3.D0*b*vlpm**2 - 12.D0*b*
     &    vlwm**2 - 12.D0*b*vdw + b*pi**2 - 6.D0*vlpm*vlsm + 12.D0*vlpm
     &    *vlwm + 3.D0*vlpm**2 + 12.D0*vdmp + pi**2 )
      gcoeff5 = gcoeff5 + xnsq * ( 6.D0*b*vlpm**2 - 96.D0*b*vltm*vlwm
     &     - 12.D0*b*vltm + 24.D0*b*vltm**2 + 96.D0*b*vlwm + 48.D0*b*
     &    vlwm**2 - 24.D0*b*vdt + 38.D0*b*pi**2 - 12.D0*b - 24.D0*vlpm*
     &    vlsm + 48.D0*vlpm*vlwm + 24.D0*vlpm + 12.D0*vlpm**2 + 48.D0*
     &    vdmp + 4.D0*pi**2 )
      gcoeff5 = gcoeff5 + xnsq**2 * ( 48.D0*b*vlsm*vlwm - 24.D0*b*vlsm*
     &    omro**(-1) + 48.D0*b*vlsm - 12.D0*b*vlsm**2 - 36.D0*b*vltm - 
     &    144.D0*b*vlwm + 48.D0*b*vdt - 4.D0*b*pi**2 - 36.D0*b + 6.D0*
     &    vlpm**2*omro**(-1) - 12.D0*vlpm**2 + 24.D0*vdmp*omro**(-1) - 
     &    48.D0*vdmp + 2.D0*omro**(-1)*pi**2 - 4.D0*pi**2 )
      gcoeff5 = gcoeff5 + 24.D0*b*vltm - 24.D0*b*vltm**2 - 24.D0*b*vlwm
     &     + 24.D0*b*vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt - 48.D0*vlpm*
     &    vlsm + 96.D0*vlpm*vlwm + 48.D0*vlpm + 24.D0*vlpm**2 + 96.D0*
     &    vdmp + 8.D0*pi**2

      gcoeff6 =  + t1**(-2)*ro*xnsq * ( 6.D0*b*vltm**2 + 6.D0*b*vdt + b
     &    *pi**2 )
      gcoeff6 = gcoeff6 + t1**(-2)*ro * (  - 6.D0*b*vltm**2 - 6.D0*b*
     &    vdt - b*pi**2 )
      gcoeff6 = gcoeff6 + t1**(-1)*ro * (  - 6.D0*b*vltm**2 - 6.D0*b*
     &    vdt - b*pi**2 )
      gcoeff6 = gcoeff6 + t1**(-1)*UBAR**(-1)*INVG*xnsq**2 * (  - 30.D0
     &    *b*vlwm )
      gcoeff6 = gcoeff6 + t1**(-1)*UBAR**(-1)*INVG * (  - 36.D0*b*vlwm
     &     )
      gcoeff6 = gcoeff6 + t1**(-1)*UBAR**(-1)*xnsq**2 * (  - 30.D0*b*
     &    vlwm )
      gcoeff6 = gcoeff6 + t1**(-1)*UBAR**(-1) * (  - 36.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t1**(-1)*INVG*xnsq**2 * (  - 3.D0/2.D0*b*vlsm
     &    *omro**(-1) + 3.D0/2.D0*b*vlsm + 30.D0*b*vlwm + 3.D0/8.D0*
     &    vlpm**2*omro**(-1) - 3.D0/8.D0*vlpm**2 + 3.D0/2.D0*vdmp*
     &    omro**(-1) - 3.D0/2.D0*vdmp + 1.D0/8.D0*omro**(-1)*pi**2 - 1.D
     &    0/8.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**(-1)*INVG * ( 36.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t1**(-1)*xnsq * ( 72.D0*b*vltm + 24.D0*b )
      gcoeff6 = gcoeff6 + t1**(-1)*xnsq**2 * (  - 6.D0*b*vlsm*
     &    omro**(-1) - 48.D0*b*vltm - 24.D0*b + 3.D0/2.D0*vlpm**2*
     &    omro**(-1) + 6.D0*vdmp*omro**(-1) + 1.D0/2.D0*omro**(-1)*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1**(-1) * (  - 24.D0*b*vltm )
      gcoeff6 = gcoeff6 + t1*ro*xnsq * ( 24.D0*vlpm*vlsm - 24.D0*vlpm*
     &    vltm - 24.D0*vlpm*vlwm - 48.D0*vdmp - 48.D0*vdmb + 12.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1*UBAR**(-2)*xnsq**2 * (  - 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t1*UBAR**(-2) * ( 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t1*UBAR**(-1)*xnsq * (  - 72.D0*b*vlwm - 12.D0
     &    *b )
      gcoeff6 = gcoeff6 + t1*UBAR**(-1)*xnsq**2 * ( 120.D0*b*vlwm + 12.D
     &    0*b )
      gcoeff6 = gcoeff6 + t1*INVG*xnsq * (  - 102.D0*b*vlpm**2 + 168.D0
     &    *b*vltm*vlwm - 228.D0*b*vltm + 84.D0*b*vltm**2 + 24.D0*b*vlwm
     &     - 324.D0*b*vlwm**2 - 240.D0*b*vdw + 168.D0*b*vdt + 6.D0*b*
     &    pi**2 - 396.D0*vlpm*vlsm + 792.D0*vlpm*vlwm + 108.D0*vlpm + 
     &    198.D0*vlpm**2 + 792.D0*vdmp + 66.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1*INVG*xnsq**2 * (  - 120.D0*b*vlsm*vltm - 
     &    96.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 54.D0*b*vlsm**2 + 168.D0*b
     &    *vltm - 12.D0*b*vlwm - 96.D0*b*vdw - 120.D0*b*vdt + 18.D0*b*
     &    pi**2 - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1*INVG * (  - 90.D0*b*vlpm**2 + 84.D0*b*vltm
     &     - 48.D0*b*vltm**2 - 60.D0*b*vlwm - 312.D0*b*vlwm**2 - 312.D0
     &    *b*vdw - 48.D0*b*vdt + 30.D0*b*pi**2 - 312.D0*vlpm*vlsm - 48.D
     &    0*vlpm*vltm + 672.D0*vlpm*vlwm + 60.D0*vlpm + 156.D0*vlpm**2
     &     + 624.D0*vdmp + 52.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1*INVG**2*xnsq * (  - 30.D0*b*vlpm**2 - 120.D
     &    0*b*vlwm**2 - 120.D0*b*vdw + 10.D0*b*pi**2 - 72.D0*vlpm*vlsm
     &     + 144.D0*vlpm*vlwm + 36.D0*vlpm**2 + 144.D0*vdmp + 12.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1*INVG**2*xnsq**2 * ( 12.D0*b*vlsm*vlwm - 3.D
     &    0*b*vlsm**2 + 12.D0*b*vdw - b*pi**2 - 3.D0*vlpm**2 - 12.D0*
     &    vdmp - pi**2 )
      gcoeff6 = gcoeff6 + t1*INVG**2 * (  - 27.D0*b*vlpm**2 - 108.D0*b*
     &    vlwm**2 - 108.D0*b*vdw + 9.D0*b*pi**2 - 66.D0*vlpm*vlsm + 132.
     &    D0*vlpm*vlwm + 33.D0*vlpm**2 + 132.D0*vdmp + 11.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1*xnsq * (  - 64.D0*XLF*b*TR*rmuom2*xn - 18.D
     &    0*b*vlpm**2 + 144.D0*b*vltm*vlwm - 72.D0*b*vltm + 24.D0*b*
     &    vltm**2 - 72.D0*b*vlwm - 96.D0*b*vlwm**2 - 24.D0*b*vdw + 96.D0
     &    *b*vdt - 42.D0*b*pi**2 + 96.D0*b - 288.D0*vlpm*vlsm + 576.D0*
     &    vlpm*vlwm + 264.D0*vlpm + 120.D0*vlpm**2 + 576.D0*vdmp + 96.D0
     &    *vdmb + 16.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1*xnsq**2 * (  - 48.D0*b*vlsm*vltm - 48.D0*b
     &    *vlsm*vlwm - 72.D0*b*vlsm*omro**(-1) + 96.D0*b*vltm + 96.D0*b
     &    *vlwm + 24.D0*b*pi**2 + 176.D0*b*rmuom2 - 48.D0*b + 18.D0*
     &    vlpm**2*omro**(-1) - 12.D0*vlpm**2 + 72.D0*vdmp*omro**(-1) - 
     &    48.D0*vdmp + 6.D0*omro**(-1)*pi**2 - 4.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1 * (  - 36.D0*b*vlpm**2 - 72.D0*b*vltm**2
     &     - 72.D0*b*vlwm**2 - 72.D0*b*vdw - 72.D0*b*vdt + 12.D0*b*
     &    pi**2 - 72.D0*vlpm*vlsm - 192.D0*vlpm*vltm + 336.D0*vlpm*vlwm
     &     + 48.D0*vlpm + 36.D0*vlpm**2 + 144.D0*vdmp + 12.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**2*UBAR**(-2)*xnsq * (  - 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t1**2*UBAR**(-2)*xnsq**2 * ( 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t1**2*INVG*xnsq * ( 156.D0*b*vlpm**2 - 288.D0
     &    *b*vltm*vlwm + 180.D0*b*vltm - 192.D0*b*vltm**2 + 12.D0*b*
     &    vlwm + 528.D0*b*vlwm**2 + 384.D0*b*vdw - 336.D0*b*vdt - 4.D0*
     &    b*pi**2 + 1032.D0*vlpm*vlsm - 48.D0*vlpm*vltm - 2016.D0*vlpm*
     &    vlwm - 384.D0*vlpm - 516.D0*vlpm**2 - 2064.D0*vdmp - 172.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1**2*INVG*xnsq**2 * ( 168.D0*b*vlsm*vltm + 
     &    120.D0*b*vlsm*vlwm + 96.D0*b*vlsm - 72.D0*b*vlsm**2 - 156.D0*
     &    b*vltm - 36.D0*b*vlwm + 120.D0*b*vdw + 168.D0*b*vdt - 24.D0*b
     &    *pi**2 + 36.D0*vlpm**2 + 144.D0*vdmp + 12.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**2*INVG * ( 144.D0*b*vlpm**2 - 48.D0*b*
     &    vltm + 216.D0*b*vltm**2 + 48.D0*b*vlwm + 360.D0*b*vlwm**2 + 
     &    360.D0*b*vdw + 216.D0*b*vdt - 48.D0*b*pi**2 + 480.D0*vlpm*
     &    vlsm + 384.D0*vlpm*vltm - 1344.D0*vlpm*vlwm - 96.D0*vlpm - 
     &    240.D0*vlpm**2 - 960.D0*vdmp - 80.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**2*INVG**2*xnsq * ( 114.D0*b*vlpm**2 - 96.D
     &    0*b*vltm*vlwm - 48.D0*b*vltm**2 + 408.D0*b*vlwm**2 + 360.D0*b
     &    *vdw - 96.D0*b*vdt - 22.D0*b*pi**2 + 348.D0*vlpm*vlsm - 696.D0
     &    *vlpm*vlwm - 174.D0*vlpm**2 - 696.D0*vdmp - 58.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**2*INVG**2*xnsq**2 * ( 60.D0*b*vlsm*vltm
     &     + 12.D0*b*vlsm*vlwm - 18.D0*b*vlsm**2 + 12.D0*b*vdw + 60.D0*
     &    b*vdt - 6.D0*b*pi**2 + 18.D0*vlpm**2 + 72.D0*vdmp + 6.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1**2*INVG**2 * ( 96.D0*b*vlpm**2 + 36.D0*b*
     &    vltm**2 + 348.D0*b*vlwm**2 + 348.D0*b*vdw + 36.D0*b*vdt - 32.D
     &    0*b*pi**2 + 264.D0*vlpm*vlsm + 36.D0*vlpm*vltm - 564.D0*vlpm*
     &    vlwm - 132.D0*vlpm**2 - 528.D0*vdmp - 44.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**2*xnsq * (  - 72.D0*b*vltm**2 + 72.D0*b*
     &    vlwm**2 + 72.D0*b*vdw - 72.D0*b*vdt + 384.D0*vlpm*vlsm - 192.D
     &    0*vlpm*vltm - 576.D0*vlpm*vlwm - 384.D0*vlpm - 192.D0*vlpm**2
     &     - 768.D0*vdmp - 64.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**2 * ( 240.D0*vlpm*vltm - 240.D0*vlpm*vlwm
     &     )
      gcoeff6 = gcoeff6 + t1**3*INVG*xnsq * (  - 78.D0*b*vlpm**2 + 168.D
     &    0*b*vltm*vlwm - 48.D0*b*vltm + 300.D0*b*vltm**2 - 48.D0*b*
     &    vlwm - 444.D0*b*vlwm**2 - 360.D0*b*vdw + 384.D0*b*vdt - 2.D0*
     &    b*pi**2 - 1464.D0*vlpm*vlsm + 384.D0*vlpm*vltm + 2544.D0*vlpm
     &    *vlwm + 672.D0*vlpm + 732.D0*vlpm**2 + 2928.D0*vdmp + 244.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1**3*INVG*xnsq**2 * (  - 48.D0*b*vlsm*vltm
     &     - 48.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 24.D0*b*vlsm**2 + 48.D0
     &    *b*vltm + 48.D0*b*vlwm - 48.D0*b*vdw - 48.D0*b*vdt + 8.D0*b*
     &    pi**2 - 30.D0*vlpm**2 - 120.D0*vdmp - 10.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**3*INVG * (  - 84.D0*b*vlpm**2 - 168.D0*b*
     &    vltm**2 - 168.D0*b*vlwm**2 - 168.D0*b*vdw - 168.D0*b*vdt + 28.
     &    D0*b*pi**2 - 216.D0*vlpm*vlsm - 768.D0*vlpm*vltm + 1200.D0*
     &    vlpm*vlwm + 48.D0*vlpm + 108.D0*vlpm**2 + 432.D0*vdmp + 36.D0
     &    *pi**2 )
      gcoeff6 = gcoeff6 + t1**3*INVG**2*xnsq * (  - 207.D0*b*vlpm**2 + 
     &    288.D0*b*vltm*vlwm + 180.D0*b*vltm**2 - 720.D0*b*vlwm**2 - 
     &    576.D0*b*vdw + 324.D0*b*vdt + 21.D0*b*pi**2 - 906.D0*vlpm*
     &    vlsm + 36.D0*vlpm*vltm + 1776.D0*vlpm*vlwm + 453.D0*vlpm**2
     &     + 1812.D0*vdmp + 151.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**3*INVG**2*xnsq**2 * (  - 168.D0*b*vlsm*
     &    vltm - 108.D0*b*vlsm*vlwm + 69.D0*b*vlsm**2 - 108.D0*b*vdw - 
     &    168.D0*b*vdt + 23.D0*b*pi**2 - 39.D0*vlpm**2 - 156.D0*vdmp - 
     &    13.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**3*INVG**2 * (  - 168.D0*b*vlpm**2 - 156.D0
     &    *b*vltm**2 - 516.D0*b*vlwm**2 - 516.D0*b*vdw - 156.D0*b*vdt
     &     + 56.D0*b*pi**2 - 492.D0*vlpm*vlsm - 228.D0*vlpm*vltm + 1212.
     &    D0*vlpm*vlwm + 246.D0*vlpm**2 + 984.D0*vdmp + 82.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**3*xnsq * (  - 240.D0*vlpm*vlsm + 240.D0*
     &    vlpm*vltm + 240.D0*vlpm*vlwm + 192.D0*vlpm + 120.D0*vlpm**2
     &     + 480.D0*vdmp + 40.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**4*INVG*xnsq * (  - 168.D0*b*vltm**2 + 168.
     &    D0*b*vlwm**2 + 168.D0*b*vdw - 168.D0*b*vdt + 1200.D0*vlpm*
     &    vlsm - 768.D0*vlpm*vltm - 1632.D0*vlpm*vlwm - 576.D0*vlpm - 
     &    600.D0*vlpm**2 - 2400.D0*vdmp - 200.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**4*INVG * ( 432.D0*vlpm*vltm - 432.D0*vlpm
     &    *vlwm )
      gcoeff6 = gcoeff6 + t1**4*INVG**2*xnsq * ( 180.D0*b*vlpm**2 - 288.
     &    D0*b*vltm*vlwm - 300.D0*b*vltm**2 + 732.D0*b*vlwm**2 + 588.D0
     &    *b*vdw - 444.D0*b*vdt - 12.D0*b*pi**2 + 1416.D0*vlpm*vlsm - 
     &    228.D0*vlpm*vltm - 2604.D0*vlpm*vlwm - 708.D0*vlpm**2 - 2832.D
     &    0*vdmp - 236.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**4*INVG**2*xnsq**2 * ( 156.D0*b*vlsm*vltm
     &     + 132.D0*b*vlsm*vlwm - 72.D0*b*vlsm**2 + 132.D0*b*vdw + 156.D
     &    0*b*vdt - 24.D0*b*pi**2 + 36.D0*vlpm**2 + 144.D0*vdmp + 12.D0
     &    *pi**2 )
      gcoeff6 = gcoeff6 + t1**4*INVG**2 * ( 144.D0*b*vlpm**2 + 216.D0*b
     &    *vltm**2 + 360.D0*b*vlwm**2 + 360.D0*b*vdw + 216.D0*b*vdt - 
     &    48.D0*b*pi**2 + 432.D0*vlpm*vlsm + 528.D0*vlpm*vltm - 1392.D0
     &    *vlpm*vlwm - 216.D0*vlpm**2 - 864.D0*vdmp - 72.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**5*INVG*xnsq * (  - 432.D0*vlpm*vlsm + 432.
     &    D0*vlpm*vltm + 432.D0*vlpm*vlwm + 192.D0*vlpm + 216.D0*
     &    vlpm**2 + 864.D0*vdmp + 72.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**5*INVG**2*xnsq * (  - 60.D0*b*vlpm**2 + 
     &    96.D0*b*vltm*vlwm + 264.D0*b*vltm**2 - 408.D0*b*vlwm**2 - 360.
     &    D0*b*vdw + 312.D0*b*vdt + 4.D0*b*pi**2 - 1368.D0*vlpm*vlsm + 
     &    528.D0*vlpm*vltm + 2208.D0*vlpm*vlwm + 684.D0*vlpm**2 + 2736.D
     &    0*vdmp + 228.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**5*INVG**2*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm - 48.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.
     &    D0*b*vdt + 8.D0*b*pi**2 - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*
     &    pi**2 )
      gcoeff6 = gcoeff6 + t1**5*INVG**2 * (  - 48.D0*b*vlpm**2 - 96.D0*
     &    b*vltm**2 - 96.D0*b*vlwm**2 - 96.D0*b*vdw - 96.D0*b*vdt + 16.D
     &    0*b*pi**2 - 144.D0*vlpm*vlsm - 528.D0*vlpm*vltm + 816.D0*vlpm
     &    *vlwm + 72.D0*vlpm**2 + 288.D0*vdmp + 24.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**6*INVG**2*xnsq * (  - 96.D0*b*vltm**2 + 
     &    96.D0*b*vlwm**2 + 96.D0*b*vdw - 96.D0*b*vdt + 768.D0*vlpm*
     &    vlsm - 528.D0*vlpm*vltm - 1008.D0*vlpm*vlwm - 384.D0*vlpm**2
     &     - 1536.D0*vdmp - 128.D0*pi**2 )
      gcoeff6 = gcoeff6 + t1**6*INVG**2 * ( 192.D0*vlpm*vltm - 192.D0*
     &    vlpm*vlwm )
      gcoeff6 = gcoeff6 + t1**7*INVG**2*xnsq * (  - 192.D0*vlpm*vlsm + 
     &    192.D0*vlpm*vltm + 192.D0*vlpm*vlwm + 96.D0*vlpm**2 + 384.D0*
     &    vdmp + 32.D0*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-2)*ro*xnsq * ( 12.D0*b*vlwm**2 + 12.D0*
     &    b*vdw + 2.D0*b*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-2)*ro * (  - 12.D0*b*vlwm**2 - 12.D0*b*
     &    vdw - 2.D0*b*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-1)*ro*xnsq * (  - 18.D0*b*vlwm**2 - 18.D
     &    0*b*vdw - 3.D0*b*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-1)*ro * ( 6.D0*b*vlwm**2 + 6.D0*b*vdw
     &     + b*pi**2 + 12.D0*vlpm*vlsm - 24.D0*vlpm*vlwm - 24.D0*vdmp
     &     - 24.D0*vdmb + 6.D0*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-1)*UBAR**(-1)*xnsq**2 * (  - 30.D0*b*
     &    vlwm )
      gcoeff6 = gcoeff6 + t2**(-1)*UBAR**(-1) * (  - 36.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t2**(-1)*INVG*xnsq**2 * (  - 3.D0/2.D0*b*vlsm
     &    *omro**(-1) + 3.D0/2.D0*b*vlsm + 30.D0*b*vlwm + 3.D0/8.D0*
     &    vlpm**2*omro**(-1) - 3.D0/8.D0*vlpm**2 + 3.D0/2.D0*vdmp*
     &    omro**(-1) - 3.D0/2.D0*vdmp + 1.D0/8.D0*omro**(-1)*pi**2 - 1.D
     &    0/8.D0*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-1)*INVG * ( 36.D0*b*vlwm )
      gcoeff6 = gcoeff6 + t2**(-1)*xnsq * ( 32.D0*XLF*b*TR*rmuom2*xn - 
     &    96.D0*b*vltm*vlwm + 240.D0*b*vlwm + 24.D0*b*pi**2 + 88.D0*b*
     &    rmuom2 - 48.D0*b )
      gcoeff6 = gcoeff6 + t2**(-1)*xnsq**2 * ( 48.D0*b*vlsm*vlwm - 6.D0
     &    *b*vlsm*omro**(-1) - 144.D0*b*vlwm - 12.D0*b*pi**2 - 88.D0*b*
     &    rmuom2 + 3.D0/2.D0*vlpm**2*omro**(-1) + 6.D0*vdmp*omro**(-1)
     &     + 1.D0/2.D0*omro**(-1)*pi**2 )
      gcoeff6 = gcoeff6 + t2**(-1) * (  - 32.D0*XLF*b*TR*rmuom2*xn - 96.
     &    D0*b*vlwm + 48.D0*b - 24.D0*vlpm*vlsm + 48.D0*vlpm*vlwm + 48.D
     &    0*vdmp + 48.D0*vdmb - 12.D0*pi**2 )
      gcoeff6 = gcoeff6 + ro*TBAR**(-2)*xnsq * ( 6.D0*b*vltm )
      gcoeff6 = gcoeff6 + ro*TBAR**(-2)*xnsq**2 * (  - 3.D0*b*vltm )
      gcoeff6 = gcoeff6 + ro*TBAR**(-2) * (  - 3.D0*b*vltm )
      gcoeff6 = gcoeff6 + ro*TBAR**(-1)*xnsq * ( 12.D0*b*vltm + 3.D0*b
     &     )
      gcoeff6 = gcoeff6 + ro*TBAR**(-1)*xnsq**2 * (  - 24.D0*b*vltm - 3.
     &    D0*b )
      gcoeff6 = gcoeff6 + ro*xnsq * (  - 6.D0*b*vltm**2 + 6.D0*b*
     &    vlwm**2 + 6.D0*b*vdw - 6.D0*b*vdt )
      gcoeff6 = gcoeff6 + ro * (  - 24.D0*vlpm*vltm + 24.D0*vlpm*vlwm )
      gcoeff6 = gcoeff6 + ro**2*TBAR**(-2)*xnsq * (  - 3.D0/4.D0*b*vltm
     &     )
      gcoeff6 = gcoeff6 + ro**2*TBAR**(-2)*xnsq**2 * ( 3.D0/4.D0*b*vltm
     &     )
      gcoeff6 = gcoeff6 + TBAR**(-1)*xnsq * (  - 168.D0*b*vltm - 24.D0*
     &    b )
      gcoeff6 = gcoeff6 + TBAR**(-1)*xnsq**2 * ( 108.D0*b*vltm + 12.D0*
     &    b )
      gcoeff6 = gcoeff6 + TBAR**(-1) * ( 60.D0*b*vltm + 12.D0*b )
      gcoeff6 = gcoeff6 + UBAR**(-2)*xnsq * ( 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + UBAR**(-2) * (  - 12.D0*b*vlwm )
      gcoeff6 = gcoeff6 + UBAR**(-1)*INVG*xnsq**2 * ( 30.D0*b*vlwm )
      gcoeff6 = gcoeff6 + UBAR**(-1)*INVG * ( 36.D0*b*vlwm )
      gcoeff6 = gcoeff6 + UBAR**(-1)*xnsq * (  - 120.D0*b*vlwm - 12.D0*
     &    b )
      gcoeff6 = gcoeff6 + UBAR**(-1) * ( 72.D0*b*vlwm + 12.D0*b )
      gcoeff6 = gcoeff6 + INVG*xnsq * ( 24.D0*b*vlpm**2 + 96.D0*b*vltm
     &     + 12.D0*b*vlwm + 96.D0*b*vlwm**2 + 96.D0*b*vdw - 8.D0*b*
     &    pi**2 + 60.D0*vlpm*vlsm - 120.D0*vlpm*vlwm - 12.D0*vlpm - 30.D
     &    0*vlpm**2 - 120.D0*vdmp - 10.D0*pi**2 )
      gcoeff6 = gcoeff6 + INVG*xnsq**2 * ( 6.D0*b*vlsm*omro**(-1) - 60.D
     &    0*b*vltm - 3.D0/2.D0*vlpm**2*omro**(-1) - 6.D0*vdmp*
     &    omro**(-1) - 1.D0/2.D0*omro**(-1)*pi**2 )
      gcoeff6 = gcoeff6 + INVG * ( 24.D0*b*vlpm**2 - 36.D0*b*vltm + 12.D
     &    0*b*vlwm + 96.D0*b*vlwm**2 + 96.D0*b*vdw - 8.D0*b*pi**2 + 60.D
     &    0*vlpm*vlsm - 120.D0*vlpm*vlwm - 12.D0*vlpm - 30.D0*vlpm**2
     &     - 120.D0*vdmp - 10.D0*pi**2 )
      gcoeff6 = gcoeff6 + INVG**2*xnsq * ( 3.D0*b*vlpm**2 + 12.D0*b*
     &    vlwm**2 + 12.D0*b*vdw - b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vlwm - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff6 = gcoeff6 + INVG**2 * ( 3.D0*b*vlpm**2 + 12.D0*b*vlwm**2
     &     + 12.D0*b*vdw - b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm*vlwm
     &     - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff6 = gcoeff6 + xnsq * ( 6.D0*b*vlpm**2 + 48.D0*b*vltm*vlwm
     &     + 204.D0*b*vltm - 48.D0*b*vltm**2 - 192.D0*b*vlwm - 24.D0*b*
     &    vlwm**2 - 24.D0*b*vdt - 34.D0*b*pi**2 + 12.D0*b + 96.D0*vlpm*
     &    vlsm - 192.D0*vlpm*vlwm - 72.D0*vlpm - 48.D0*vlpm**2 - 192.D0
     &    *vdmp - 16.D0*pi**2 )
      gcoeff6 = gcoeff6 + xnsq**2 * ( 48.D0*b*vlsm*vltm - 48.D0*b*vlsm*
     &    vlwm + 72.D0*b*vlsm*omro**(-1) - 156.D0*b*vltm + 144.D0*b*
     &    vlwm + 36.D0*b - 18.D0*vlpm**2*omro**(-1) + 36.D0*vlpm**2 - 
     &    72.D0*vdmp*omro**(-1) + 144.D0*vdmp - 6.D0*omro**(-1)*pi**2
     &     + 12.D0*pi**2 )
      gcoeff6 = gcoeff6 + 12.D0*b*vlpm**2 - 72.D0*b*vltm + 24.D0*b*
     &    vltm**2 + 72.D0*b*vlwm + 24.D0*b*vlwm**2 + 24.D0*b*vdw + 24.D0
     &    *b*vdt - 4.D0*b*pi**2 + 120.D0*vlpm*vlsm - 240.D0*vlpm*vlwm
     &     - 48.D0*vlpm - 60.D0*vlpm**2 - 240.D0*vdmp - 20.D0*pi**2

      gcoeff7 =  + t1**(-1)*INVG*xnsq**2 * (  - 3.D0*b*vlsm*omro**(-1)
     &     + 3.D0*b*vlsm + 3.D0/4.D0*vlpm**2*omro**(-1) - 3.D0/4.D0*
     &    vlpm**2 + 3.D0*vdmp*omro**(-1) - 3.D0*vdmp + 1.D0/4.D0*
     &    omro**(-1)*pi**2 - 1.D0/4.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**(-1)*xnsq**2 * (  - 12.D0*b*vlsm*
     &    omro**(-1) + 3.D0*vlpm**2*omro**(-1) + 12.D0*vdmp*omro**(-1)
     &     + omro**(-1)*pi**2 )
      gcoeff7 = gcoeff7 + t1*ro*xnsq * (  - 12.D0*vlpm*vlsm + 12.D0*
     &    vlpm*vltm + 12.D0*vlpm*vlwm + 24.D0*vlpm + 6.D0*vlpm**2 + 24.D
     &    0*vdmp + 2.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1*UBAR**(-2)*xnsq**2 * (  - 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + t1*UBAR**(-2) * ( 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + t1*UBAR**(-1)*xnsq * (  - 36.D0*b*vlwm - 12.D0
     &    *b )
      gcoeff7 = gcoeff7 + t1*UBAR**(-1)*xnsq**2 * ( 84.D0*b*vlwm + 12.D0
     &    *b )
      gcoeff7 = gcoeff7 + t1*INVG*xnsq * ( 18.D0*b*vlpm**2 - 48.D0*b*
     &    vltm*vlwm + 72.D0*b*vltm - 24.D0*b*vltm**2 + 48.D0*b*vlwm**2
     &     + 24.D0*b*vdw - 48.D0*b*vdt + 2.D0*b*pi**2 + 36.D0*vlpm*vlsm
     &     - 72.D0*vlpm*vlwm - 18.D0*vlpm**2 - 72.D0*vdmp - 6.D0*pi**2
     &     )
      gcoeff7 = gcoeff7 + t1*INVG*xnsq**2 * ( 24.D0*b*vlsm*vltm - 6.D0*
     &    b*vlsm**2 - 48.D0*b*vltm + 24.D0*b*vdt - 2.D0*b*pi**2 + 6.D0*
     &    vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1*INVG * ( 24.D0*b*vlpm**2 - 24.D0*b*vltm + 
     &    24.D0*b*vltm**2 + 72.D0*b*vlwm**2 + 72.D0*b*vdw + 24.D0*b*vdt
     &     - 8.D0*b*pi**2 + 24.D0*vlpm*vlsm + 24.D0*vlpm*vltm - 72.D0*
     &    vlpm*vlwm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1*xnsq * ( 18.D0*b*vlpm**2 - 72.D0*b*vltm*
     &    vlwm + 24.D0*b*vltm - 60.D0*b*vltm**2 + 24.D0*b*vlwm + 60.D0*
     &    b*vlwm**2 + 24.D0*b*vdw - 96.D0*b*vdt + 6.D0*b*pi**2 + 132.D0
     &    *vlpm*vlsm - 24.D0*vlpm*vltm - 240.D0*vlpm*vlwm - 24.D0*vlpm
     &     - 66.D0*vlpm**2 - 264.D0*vdmp - 22.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1*xnsq**2 * ( 24.D0*b*vlsm*vltm + 24.D0*b*
     &    vlsm*vlwm + 48.D0*b*vlsm - 12.D0*b*vlsm**2 - 48.D0*b*vltm - 
     &    48.D0*b*vlwm + 24.D0*b*vdw + 24.D0*b*vdt - 4.D0*b*pi**2 + 12.D
     &    0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1 * ( 24.D0*b*vlpm**2 + 48.D0*b*vltm**2 + 48.
     &    D0*b*vlwm**2 + 48.D0*b*vdw + 48.D0*b*vdt - 8.D0*b*pi**2 + 48.D
     &    0*vlpm*vlsm + 96.D0*vlpm*vltm - 192.D0*vlpm*vlwm - 24.D0*
     &    vlpm**2 - 96.D0*vdmp - 8.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**2*UBAR**(-2)*xnsq * (  - 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + t1**2*UBAR**(-2)*xnsq**2 * ( 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + t1**2*INVG*xnsq * (  - 54.D0*b*vlpm**2 + 192.D
     &    0*b*vltm*vlwm - 96.D0*b*vltm + 120.D0*b*vltm**2 + 48.D0*b*
     &    vlwm - 144.D0*b*vlwm**2 - 48.D0*b*vdw + 216.D0*b*vdt - 14.D0*
     &    b*pi**2 - 228.D0*vlpm*vlsm + 24.D0*vlpm*vltm + 432.D0*vlpm*
     &    vlwm + 24.D0*vlpm + 114.D0*vlpm**2 + 456.D0*vdmp + 38.D0*
     &    pi**2 )
      gcoeff7 = gcoeff7 + t1**2*INVG*xnsq**2 * (  - 96.D0*b*vlsm*vltm
     &     - 72.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 42.D0*b*vlsm**2 + 96.D0
     &    *b*vltm - 72.D0*b*vdw - 96.D0*b*vdt + 14.D0*b*pi**2 - 18.D0*
     &    vlpm**2 - 72.D0*vdmp - 6.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**2*INVG * (  - 66.D0*b*vlpm**2 + 24.D0*b*
     &    vltm - 96.D0*b*vltm**2 - 24.D0*b*vlwm - 168.D0*b*vlwm**2 - 
     &    168.D0*b*vdw - 96.D0*b*vdt + 22.D0*b*pi**2 - 108.D0*vlpm*vlsm
     &     - 144.D0*vlpm*vltm + 360.D0*vlpm*vlwm + 54.D0*vlpm**2 + 216.D
     &    0*vdmp + 18.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**2*INVG**2*xnsq * (  - 6.D0*b*vlpm**2 - 24.
     &    D0*b*vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm
     &     + 24.D0*vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**2*INVG**2 * (  - 6.D0*b*vlpm**2 - 24.D0*b
     &    *vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 24.D
     &    0*vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**2*xnsq * ( 48.D0*b*vltm**2 - 48.D0*b*
     &    vlwm**2 - 48.D0*b*vdw + 48.D0*b*vdt - 216.D0*vlpm*vlsm + 96.D0
     &    *vlpm*vltm + 336.D0*vlpm*vlwm + 96.D0*vlpm + 108.D0*vlpm**2
     &     + 432.D0*vdmp + 36.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**2 * (  - 144.D0*vlpm*vltm + 144.D0*vlpm*
     &    vlwm )
      gcoeff7 = gcoeff7 + t1**3*INVG*xnsq * ( 42.D0*b*vlpm**2 - 144.D0*
     &    b*vltm*vlwm + 24.D0*b*vltm - 168.D0*b*vltm**2 + 24.D0*b*vlwm
     &     + 192.D0*b*vlwm**2 + 120.D0*b*vdw - 240.D0*b*vdt + 10.D0*b*
     &    pi**2 + 516.D0*vlpm*vlsm - 144.D0*vlpm*vltm - 888.D0*vlpm*
     &    vlwm - 120.D0*vlpm - 258.D0*vlpm**2 - 1032.D0*vdmp - 86.D0*
     &    pi**2 )
      gcoeff7 = gcoeff7 + t1**3*INVG*xnsq**2 * ( 72.D0*b*vlsm*vltm + 72.
     &    D0*b*vlsm*vlwm + 48.D0*b*vlsm - 36.D0*b*vlsm**2 - 48.D0*b*
     &    vltm - 48.D0*b*vlwm + 72.D0*b*vdw + 72.D0*b*vdt - 12.D0*b*
     &    pi**2 + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**3*INVG * ( 48.D0*b*vlpm**2 + 96.D0*b*
     &    vltm**2 + 96.D0*b*vlwm**2 + 96.D0*b*vdw + 96.D0*b*vdt - 16.D0
     &    *b*pi**2 + 96.D0*vlpm*vlsm + 336.D0*vlpm*vltm - 528.D0*vlpm*
     &    vlwm - 48.D0*vlpm**2 - 192.D0*vdmp - 16.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**3*INVG**2*xnsq * ( 30.D0*b*vlpm**2 - 72.D0
     &    *b*vltm*vlwm - 36.D0*b*vltm**2 + 84.D0*b*vlwm**2 + 48.D0*b*
     &    vdw - 72.D0*b*vdt + 2.D0*b*pi**2 + 84.D0*vlpm*vlsm - 168.D0*
     &    vlpm*vlwm - 42.D0*vlpm**2 - 168.D0*vdmp - 14.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**3*INVG**2*xnsq**2 * ( 48.D0*b*vlsm*vltm
     &     + 48.D0*b*vlsm*vlwm - 24.D0*b*vlsm**2 + 48.D0*b*vdw + 48.D0*
     &    b*vdt - 8.D0*b*pi**2 )
      gcoeff7 = gcoeff7 + t1**3*INVG**2 * ( 30.D0*b*vlpm**2 + 24.D0*b*
     &    vltm**2 + 96.D0*b*vlwm**2 + 96.D0*b*vdw + 24.D0*b*vdt - 10.D0
     &    *b*pi**2 + 60.D0*vlpm*vlsm + 24.D0*vlpm*vltm - 144.D0*vlpm*
     &    vlwm - 30.D0*vlpm**2 - 120.D0*vdmp - 10.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**3*xnsq * ( 144.D0*vlpm*vlsm - 144.D0*vlpm
     &    *vltm - 144.D0*vlpm*vlwm - 96.D0*vlpm - 72.D0*vlpm**2 - 288.D0
     &    *vdmp - 24.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**4*INVG*xnsq * ( 96.D0*b*vltm**2 - 96.D0*b
     &    *vlwm**2 - 96.D0*b*vdw + 96.D0*b*vdt - 552.D0*vlpm*vlsm + 336.
     &    D0*vlpm*vltm + 768.D0*vlpm*vlwm + 192.D0*vlpm + 276.D0*
     &    vlpm**2 + 1104.D0*vdmp + 92.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**4*INVG * (  - 240.D0*vlpm*vltm + 240.D0*
     &    vlpm*vlwm )
      gcoeff7 = gcoeff7 + t1**4*INVG**2*xnsq * (  - 48.D0*b*vlpm**2 + 
     &    144.D0*b*vltm*vlwm + 96.D0*b*vltm**2 - 144.D0*b*vlwm**2 - 72.D
     &    0*b*vdw + 168.D0*b*vdt - 8.D0*b*pi**2 - 240.D0*vlpm*vlsm + 24.
     &    D0*vlpm*vltm + 456.D0*vlpm*vlwm + 120.D0*vlpm**2 + 480.D0*
     &    vdmp + 40.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**4*INVG**2*xnsq**2 * (  - 96.D0*b*vlsm*
     &    vltm - 96.D0*b*vlsm*vlwm + 48.D0*b*vlsm**2 - 96.D0*b*vdw - 96.
     &    D0*b*vdt + 16.D0*b*pi**2 )
      gcoeff7 = gcoeff7 + t1**4*INVG**2 * (  - 48.D0*b*vlpm**2 - 72.D0*
     &    b*vltm**2 - 120.D0*b*vlwm**2 - 120.D0*b*vdw - 72.D0*b*vdt + 
     &    16.D0*b*pi**2 - 96.D0*vlpm*vlsm - 120.D0*vlpm*vltm + 312.D0*
     &    vlpm*vlwm + 48.D0*vlpm**2 + 192.D0*vdmp + 16.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**5*INVG*xnsq * ( 240.D0*vlpm*vlsm - 240.D0
     &    *vlpm*vltm - 240.D0*vlpm*vlwm - 96.D0*vlpm - 120.D0*vlpm**2
     &     - 480.D0*vdmp - 40.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**5*INVG**2*xnsq * ( 24.D0*b*vlpm**2 - 72.D0
     &    *b*vltm*vlwm - 108.D0*b*vltm**2 + 132.D0*b*vlwm**2 + 96.D0*b*
     &    vdw - 144.D0*b*vdt + 4.D0*b*pi**2 + 360.D0*vlpm*vlsm - 120.D0
     &    *vlpm*vltm - 600.D0*vlpm*vlwm - 180.D0*vlpm**2 - 720.D0*vdmp
     &     - 60.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**5*INVG**2*xnsq**2 * ( 48.D0*b*vlsm*vltm
     &     + 48.D0*b*vlsm*vlwm - 24.D0*b*vlsm**2 + 48.D0*b*vdw + 48.D0*
     &    b*vdt - 8.D0*b*pi**2 )
      gcoeff7 = gcoeff7 + t1**5*INVG**2 * ( 24.D0*b*vlpm**2 + 48.D0*b*
     &    vltm**2 + 48.D0*b*vlwm**2 + 48.D0*b*vdw + 48.D0*b*vdt - 8.D0*
     &    b*pi**2 + 48.D0*vlpm*vlsm + 192.D0*vlpm*vltm - 288.D0*vlpm*
     &    vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**6*INVG**2*xnsq * ( 48.D0*b*vltm**2 - 48.D0
     &    *b*vlwm**2 - 48.D0*b*vdw + 48.D0*b*vdt - 288.D0*vlpm*vlsm + 
     &    192.D0*vlpm*vltm + 384.D0*vlpm*vlwm + 144.D0*vlpm**2 + 576.D0
     &    *vdmp + 48.D0*pi**2 )
      gcoeff7 = gcoeff7 + t1**6*INVG**2 * (  - 96.D0*vlpm*vltm + 96.D0*
     &    vlpm*vlwm )
      gcoeff7 = gcoeff7 + t1**7*INVG**2*xnsq * ( 96.D0*vlpm*vlsm - 96.D0
     &    *vlpm*vltm - 96.D0*vlpm*vlwm - 48.D0*vlpm**2 - 192.D0*vdmp - 
     &    16.D0*pi**2 )
      gcoeff7 = gcoeff7 + t2**(-1)*INVG*xnsq**2 * (  - 3.D0*b*vlsm*
     &    omro**(-1) + 3.D0*b*vlsm + 3.D0/4.D0*vlpm**2*omro**(-1) - 3.D0
     &    /4.D0*vlpm**2 + 3.D0*vdmp*omro**(-1) - 3.D0*vdmp + 1.D0/4.D0*
     &    omro**(-1)*pi**2 - 1.D0/4.D0*pi**2 )
      gcoeff7 = gcoeff7 + t2**(-1)*xnsq * ( 24.D0*b*vlwm + 24.D0*b )
      gcoeff7 = gcoeff7 + t2**(-1)*xnsq**2 * (  - 12.D0*b*vlsm*
     &    omro**(-1) - 24.D0*b*vlwm - 12.D0*b + 3.D0*vlpm**2*omro**(-1)
     &     + 12.D0*vdmp*omro**(-1) + omro**(-1)*pi**2 )
      gcoeff7 = gcoeff7 + t2**(-1) * (  - 12.D0*b )
      gcoeff7 = gcoeff7 + ro*xnsq * ( 6.D0*vlpm*vlsm - 12.D0*vlpm*vlwm
     &     - 12.D0*vlpm - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff7 = gcoeff7 + ro * ( 12.D0*vlpm*vltm - 12.D0*vlpm*vlwm )
      gcoeff7 = gcoeff7 + UBAR**(-2)*xnsq * ( 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + UBAR**(-2) * (  - 12.D0*b*vlwm )
      gcoeff7 = gcoeff7 + UBAR**(-1)*xnsq * (  - 84.D0*b*vlwm - 12.D0*b
     &     )
      gcoeff7 = gcoeff7 + UBAR**(-1) * ( 36.D0*b*vlwm + 12.D0*b )
      gcoeff7 = gcoeff7 + INVG*xnsq**2 * ( 12.D0*b*vlsm*omro**(-1) - 3.D
     &    0*vlpm**2*omro**(-1) - 12.D0*vdmp*omro**(-1) - omro**(-1)*
     &    pi**2 )
      gcoeff7 = gcoeff7 + xnsq * (  - 12.D0*b*vlpm**2 + 48.D0*b*vltm*
     &    vlwm - 72.D0*b*vltm + 24.D0*b*vltm**2 + 48.D0*b*vlwm - 24.D0*
     &    b*vlwm**2 + 48.D0*b*vdt - 4.D0*b*pi**2 - 12.D0*b - 24.D0*vlpm
     &    *vlsm + 48.D0*vlpm*vlwm + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*
     &    pi**2 )
      gcoeff7 = gcoeff7 + xnsq**2 * (  - 24.D0*b*vlsm*vltm + 24.D0*b*
     &    vlsm*omro**(-1) - 24.D0*b*vlsm + 6.D0*b*vlsm**2 + 48.D0*b*
     &    vltm + 24.D0*b*vlwm - 24.D0*b*vdt + 2.D0*b*pi**2 + 12.D0*b - 
     &    6.D0*vlpm**2*omro**(-1) - 24.D0*vdmp*omro**(-1) - 2.D0*
     &    omro**(-1)*pi**2 )
      gcoeff7 = gcoeff7 - 18.D0*b*vlpm**2 + 24.D0*b*vltm - 24.D0*b*
     &    vltm**2 - 24.D0*b*vlwm - 48.D0*b*vlwm**2 - 48.D0*b*vdw - 24.D0
     &    *b*vdt + 6.D0*b*pi**2 - 12.D0*vlpm*vlsm - 24.D0*vlpm*vltm + 
     &    48.D0*vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2

      gcoeff8 =  + t1**(-1)*INVG*xnsq**2 * (  - 9.D0/4.D0*b*vlsm*
     &    omro**(-2) + 3.D0/4.D0*b*vlsm*omro**(-1) - 3.D0/2.D0*b*vlsm
     &     + 3.D0/2.D0*b*omro**(-1) - 3.D0/2.D0*b + 9.D0/16.D0*vlpm**2*
     &    omro**(-2) - 3.D0/8.D0*vlpm**2*omro**(-1) + 9.D0/16.D0*
     &    vlpm**2 + 9.D0/4.D0*vdmp*omro**(-2) - 3.D0/2.D0*vdmp*
     &    omro**(-1) + 9.D0/4.D0*vdmp + 3.D0/16.D0*omro**(-2)*pi**2 - 1.
     &    D0/8.D0*omro**(-1)*pi**2 + 3.D0/16.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**(-1)*INVG**2*xnsq**2 * (  - 3.D0/4.D0*b*
     &    vlsm*omro**(-1) + 3.D0/4.D0*b*vlsm + 3.D0/16.D0*vlpm**2*
     &    omro**(-1) - 3.D0/16.D0*vlpm**2 + 3.D0/4.D0*vdmp*omro**(-1)
     &     - 3.D0/4.D0*vdmp + 1.D0/16.D0*omro**(-1)*pi**2 - 1.D0/16.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**(-1)*xnsq**2 * (  - 9.D0*b*vlsm*
     &    omro**(-2) + 6.D0*b*vlsm*omro**(-1) + 6.D0*b*omro**(-1) + 9.D0
     &    /4.D0*vlpm**2*omro**(-2) - 9.D0/4.D0*vlpm**2*omro**(-1) + 9.D0
     &    *vdmp*omro**(-2) - 9.D0*vdmp*omro**(-1) + 3.D0/4.D0*
     &    omro**(-2)*pi**2 - 3.D0/4.D0*omro**(-1)*pi**2 )
      gcoeff8 = gcoeff8 + t1*ro*xnsq * ( 36.D0*b*vltm**2 - 36.D0*b*
     &    vlwm**2 - 36.D0*b*vdw + 36.D0*b*vdt - 72.D0*vlpm*vlsm + 24.D0
     &    *vlpm*vltm + 120.D0*vlpm*vlwm - 24.D0*vlpm + 36.D0*vlpm**2 + 
     &    144.D0*vdmp + 12.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*ro * (  - 72.D0*vlpm*vltm + 72.D0*vlpm*
     &    vlwm )
      gcoeff8 = gcoeff8 + t1*UBAR**(-1)*xnsq * ( 24.D0*b*vlwm )
      gcoeff8 = gcoeff8 + t1*UBAR**(-1)*xnsq**2 * (  - 24.D0*b*vlwm )
      gcoeff8 = gcoeff8 + t1*INVG*xnsq * (  - 54.D0*b*vlpm**2 - 96.D0*b
     &    *vltm*vlwm + 72.D0*b*vltm - 48.D0*b*vltm**2 - 72.D0*b*vlwm - 
     &    264.D0*b*vlwm**2 - 312.D0*b*vdw - 96.D0*b*vdt + 34.D0*b*pi**2
     &     - 48.D0*b - 204.D0*vlpm*vlsm + 408.D0*vlpm*vlwm + 72.D0*vlpm
     &     + 102.D0*vlpm**2 + 408.D0*vdmp + 34.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*INVG*xnsq**2 * ( 48.D0*b*vlsm*vltm + 24.D0
     &    *b*vlsm*vlwm + 6.D0*b*vlsm - 18.D0*b*vlsm**2 - 48.D0*b*vltm
     &     + 24.D0*b*vdw + 48.D0*b*vdt - 6.D0*b*pi**2 + 24.D0*b + 9.D0/
     &    2.D0*vlpm**2 + 18.D0*vdmp + 3.D0/2.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*INVG * (  - 36.D0*b*vlpm**2 - 24.D0*b*vltm
     &     + 48.D0*b*vltm**2 - 72.D0*b*vlwm - 192.D0*b*vlwm**2 - 192.D0
     &    *b*vdw + 48.D0*b*vdt + 12.D0*b*pi**2 + 24.D0*b - 216.D0*vlpm*
     &    vlsm + 48.D0*vlpm*vltm + 384.D0*vlpm*vlwm + 72.D0*vlpm + 108.D
     &    0*vlpm**2 + 432.D0*vdmp + 36.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*INVG**2*xnsq * (  - 18.D0*b*vlpm**2 - 72.D0
     &    *b*vlwm**2 - 72.D0*b*vdw + 6.D0*b*pi**2 - 36.D0*vlpm*vlsm + 
     &    72.D0*vlpm*vlwm + 18.D0*vlpm**2 + 72.D0*vdmp + 6.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*INVG**2 * (  - 18.D0*b*vlpm**2 - 72.D0*b*
     &    vlwm**2 - 72.D0*b*vdw + 6.D0*b*pi**2 - 36.D0*vlpm*vlsm + 72.D0
     &    *vlpm*vlwm + 18.D0*vlpm**2 + 72.D0*vdmp + 6.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*xnsq * ( 24.D0*b*vlpm**2 - 192.D0*b*vltm*
     &    vlwm + 96.D0*b*vltm - 144.D0*b*vltm**2 - 96.D0*b*vlwm + 48.D0
     &    *b*vlwm**2 - 48.D0*b*vdw - 240.D0*b*vdt + 24.D0*b*pi**2 - 48.D
     &    0*b - 24.D0*vlpm*vlsm - 48.D0*vlpm*vltm + 96.D0*vlpm*vlwm + 
     &    264.D0*vlpm + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1*xnsq**2 * ( 96.D0*b*vlsm*vltm + 48.D0*b*
     &    vlsm - 24.D0*b*vlsm**2 - 96.D0*b*vltm + 96.D0*b*vdt - 8.D0*b*
     &    pi**2 + 42.D0*vlpm**2 + 168.D0*vdmp + 14.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1 * ( 48.D0*b*vlpm**2 + 96.D0*b*vltm**2 + 96.
     &    D0*b*vlwm**2 + 96.D0*b*vdw + 96.D0*b*vdt - 16.D0*b*pi**2 - 
     &    144.D0*vlpm*vlsm + 192.D0*vlpm*vltm + 96.D0*vlpm*vlwm + 192.D0
     &    *vlpm + 72.D0*vlpm**2 + 288.D0*vdmp + 24.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*ro*xnsq * ( 72.D0*vlpm*vlsm - 72.D0*
     &    vlpm*vltm - 72.D0*vlpm*vlwm + 24.D0*vlpm - 36.D0*vlpm**2 - 
     &    144.D0*vdmp - 12.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG*xnsq * ( 186.D0*b*vlpm**2 + 240.D0
     &    *b*vltm*vlwm + 24.D0*b*vltm + 168.D0*b*vltm**2 + 168.D0*b*
     &    vlwm + 816.D0*b*vlwm**2 + 936.D0*b*vdw + 288.D0*b*vdt - 102.D0
     &    *b*pi**2 + 96.D0*b + 1020.D0*vlpm*vlsm + 48.D0*vlpm*vltm - 
     &    2088.D0*vlpm*vlwm - 624.D0*vlpm - 510.D0*vlpm**2 - 2040.D0*
     &    vdmp - 170.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG*xnsq**2 * (  - 168.D0*b*vlsm*vltm
     &     - 96.D0*b*vlsm*vlwm - 54.D0*b*vlsm + 66.D0*b*vlsm**2 + 72.D0
     &    *b*vltm + 72.D0*b*vlwm - 96.D0*b*vdw - 168.D0*b*vdt + 22.D0*b
     &    *pi**2 - 24.D0*b + 3.D0/2.D0*vlpm**2 + 6.D0*vdmp + 1.D0/2.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG * ( 108.D0*b*vlpm**2 - 96.D0*b*
     &    vltm - 72.D0*b*vltm**2 + 192.D0*b*vlwm + 504.D0*b*vlwm**2 + 
     &    504.D0*b*vdw - 72.D0*b*vdt - 36.D0*b*pi**2 - 24.D0*b + 936.D0
     &    *vlpm*vlsm - 168.D0*vlpm*vltm - 1704.D0*vlpm*vlwm - 456.D0*
     &    vlpm - 468.D0*vlpm**2 - 1872.D0*vdmp - 156.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG**2*xnsq * ( 168.D0*b*vlpm**2 + 96.D
     &    0*b*vltm + 48.D0*b*vlwm + 672.D0*b*vlwm**2 + 672.D0*b*vdw - 
     &    56.D0*b*pi**2 + 456.D0*vlpm*vlsm - 912.D0*vlpm*vlwm - 48.D0*
     &    vlpm - 228.D0*vlpm**2 - 912.D0*vdmp - 76.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG**2*xnsq**2 * (  - 72.D0*b*vlsm*
     &    vlwm + 18.D0*b*vlsm**2 - 48.D0*b*vltm - 72.D0*b*vdw + 6.D0*b*
     &    pi**2 + 18.D0*vlpm**2 + 72.D0*vdmp + 6.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG**2 * ( 150.D0*b*vlpm**2 - 48.D0*b*
     &    vltm + 48.D0*b*vlwm + 600.D0*b*vlwm**2 + 600.D0*b*vdw - 50.D0
     &    *b*pi**2 + 420.D0*vlpm*vlsm - 840.D0*vlpm*vlwm - 48.D0*vlpm
     &     - 210.D0*vlpm**2 - 840.D0*vdmp - 70.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG**3*xnsq * ( 12.D0*b*vlpm**2 + 48.D0
     &    *b*vlwm**2 + 48.D0*b*vdw - 4.D0*b*pi**2 + 24.D0*vlpm*vlsm - 
     &    48.D0*vlpm*vlwm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*INVG**3 * ( 12.D0*b*vlpm**2 + 48.D0*b*
     &    vlwm**2 + 48.D0*b*vdw - 4.D0*b*pi**2 + 24.D0*vlpm*vlsm - 48.D0
     &    *vlpm*vlwm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*xnsq * (  - 24.D0*b*vlpm**2 + 192.D0*b*
     &    vltm*vlwm + 192.D0*b*vltm**2 - 96.D0*b*vlwm**2 + 288.D0*b*vdt
     &     - 24.D0*b*pi**2 + 48.D0*b + 24.D0*vlpm*vlsm + 192.D0*vlpm*
     &    vltm - 240.D0*vlpm*vlwm - 744.D0*vlpm - 12.D0*vlpm**2 - 48.D0
     &    *vdmp - 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2*xnsq**2 * (  - 48.D0*b*vlsm*vltm - 48.D0
     &    *b*vlsm*vlwm - 48.D0*b*vlsm + 24.D0*b*vlsm**2 + 48.D0*b*vltm
     &     + 48.D0*b*vlwm - 48.D0*b*vdw - 48.D0*b*vdt + 8.D0*b*pi**2 - 
     &    42.D0*vlpm**2 - 168.D0*vdmp - 14.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**2 * (  - 48.D0*b*vlpm**2 - 96.D0*b*
     &    vltm**2 - 96.D0*b*vlwm**2 - 96.D0*b*vdw - 96.D0*b*vdt + 16.D0
     &    *b*pi**2 + 144.D0*vlpm*vlsm - 144.D0*vlpm*vltm - 144.D0*vlpm*
     &    vlwm - 192.D0*vlpm - 72.D0*vlpm**2 - 288.D0*vdmp - 24.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG*xnsq * (  - 264.D0*b*vlpm**2 - 288.
     &    D0*b*vltm*vlwm - 144.D0*b*vltm - 216.D0*b*vltm**2 - 240.D0*b*
     &    vlwm - 1128.D0*b*vlwm**2 - 1272.D0*b*vdw - 360.D0*b*vdt + 136.
     &    D0*b*pi**2 - 96.D0*b - 2304.D0*vlpm*vlsm - 168.D0*vlpm*vltm
     &     + 4776.D0*vlpm*vlwm + 1968.D0*vlpm + 1152.D0*vlpm**2 + 4608.D
     &    0*vdmp + 384.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG*xnsq**2 * ( 216.D0*b*vlsm*vltm + 
     &    168.D0*b*vlsm*vlwm + 96.D0*b*vlsm - 96.D0*b*vlsm**2 - 96.D0*b
     &    *vltm - 96.D0*b*vlwm + 168.D0*b*vdw + 216.D0*b*vdt - 32.D0*b*
     &    pi**2 - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG * (  - 144.D0*b*vlpm**2 + 96.D0*b*
     &    vltm - 96.D0*b*vltm**2 - 96.D0*b*vlwm - 480.D0*b*vlwm**2 - 
     &    480.D0*b*vdw - 96.D0*b*vdt + 48.D0*b*pi**2 - 1440.D0*vlpm*
     &    vlsm - 192.D0*vlpm*vltm + 3072.D0*vlpm*vlwm + 768.D0*vlpm + 
     &    720.D0*vlpm**2 + 2880.D0*vdmp + 240.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG**2*xnsq * (  - 600.D0*b*vlpm**2 + 
     &    144.D0*b*vltm*vlwm - 336.D0*b*vltm + 72.D0*b*vltm**2 - 144.D0
     &    *b*vlwm - 2328.D0*b*vlwm**2 - 2256.D0*b*vdw + 144.D0*b*vdt + 
     &    176.D0*b*pi**2 - 2208.D0*vlpm*vlsm + 4416.D0*vlpm*vlwm + 384.D
     &    0*vlpm + 1104.D0*vlpm**2 + 4416.D0*vdmp + 368.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG**2*xnsq**2 * (  - 24.D0*b*vlsm*
     &    vltm + 264.D0*b*vlsm*vlwm - 60.D0*b*vlsm**2 + 144.D0*b*vltm
     &     - 48.D0*b*vlwm + 264.D0*b*vdw - 24.D0*b*vdt - 20.D0*b*pi**2
     &     - 132.D0*vlpm**2 - 528.D0*vdmp - 44.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG**2 * (  - 480.D0*b*vlpm**2 + 192.D0
     &    *b*vltm - 120.D0*b*vltm**2 - 192.D0*b*vlwm - 1800.D0*b*
     &    vlwm**2 - 1800.D0*b*vdw - 120.D0*b*vdt + 160.D0*b*pi**2 - 
     &    1728.D0*vlpm*vlsm - 120.D0*vlpm*vltm + 3576.D0*vlpm*vlwm + 
     &    288.D0*vlpm + 864.D0*vlpm**2 + 3456.D0*vdmp + 288.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG**3*xnsq * (  - 108.D0*b*vlpm**2 - 
     &    432.D0*b*vlwm**2 - 432.D0*b*vdw + 36.D0*b*pi**2 - 264.D0*vlpm
     &    *vlsm + 528.D0*vlpm*vlwm + 132.D0*vlpm**2 + 528.D0*vdmp + 44.D
     &    0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG**3*xnsq**2 * ( 48.D0*b*vlsm*vlwm
     &     - 12.D0*b*vlsm**2 + 48.D0*b*vdw - 4.D0*b*pi**2 - 12.D0*
     &    vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**3*INVG**3 * (  - 96.D0*b*vlpm**2 - 384.D0
     &    *b*vlwm**2 - 384.D0*b*vdw + 32.D0*b*pi**2 - 240.D0*vlpm*vlsm
     &     + 480.D0*vlpm*vlwm + 120.D0*vlpm**2 + 480.D0*vdmp + 40.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**3*xnsq * (  - 96.D0*b*vltm**2 + 96.D0*b*
     &    vlwm**2 + 96.D0*b*vdw - 96.D0*b*vdt - 144.D0*vlpm*vltm + 144.D
     &    0*vlpm*vlwm + 960.D0*vlpm )
      gcoeff8 = gcoeff8 + t1**4*INVG*xnsq * ( 132.D0*b*vlpm**2 + 144.D0
     &    *b*vltm*vlwm + 96.D0*b*vltm - 24.D0*b*vltm**2 + 96.D0*b*vlwm
     &     + 696.D0*b*vlwm**2 + 768.D0*b*vdw + 48.D0*b*vdt - 68.D0*b*
     &    pi**2 + 48.D0*b + 2832.D0*vlpm*vlsm - 192.D0*vlpm*vltm - 5472.
     &    D0*vlpm*vlwm - 3144.D0*vlpm - 1416.D0*vlpm**2 - 5664.D0*vdmp
     &     - 472.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG*xnsq**2 * (  - 96.D0*b*vlsm*vltm
     &     - 96.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 48.D0*b*vlsm**2 + 48.D0
     &    *b*vltm + 48.D0*b*vlwm - 96.D0*b*vdw - 96.D0*b*vdt + 16.D0*b*
     &    pi**2 + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG * ( 72.D0*b*vlpm**2 + 144.D0*b*
     &    vltm**2 + 144.D0*b*vlwm**2 + 144.D0*b*vdw + 144.D0*b*vdt - 24.
     &    D0*b*pi**2 + 720.D0*vlpm*vlsm + 960.D0*vlpm*vltm - 2400.D0*
     &    vlpm*vlwm - 384.D0*vlpm - 360.D0*vlpm**2 - 1440.D0*vdmp - 120.
     &    D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG**2*xnsq * ( 1050.D0*b*vlpm**2 - 
     &    432.D0*b*vltm*vlwm + 480.D0*b*vltm - 336.D0*b*vltm**2 + 240.D0
     &    *b*vlwm + 4104.D0*b*vlwm**2 + 3888.D0*b*vdw - 552.D0*b*vdt - 
     &    278.D0*b*pi**2 + 5580.D0*vlpm*vlsm - 120.D0*vlpm*vltm - 11040.
     &    D0*vlpm*vlwm - 1296.D0*vlpm - 2790.D0*vlpm**2 - 11160.D0*vdmp
     &     - 930.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG**2*xnsq**2 * (  - 360.D0*b*vlsm*
     &    vlwm + 90.D0*b*vlsm**2 - 144.D0*b*vltm + 96.D0*b*vlwm - 360.D0
     &    *b*vdw + 30.D0*b*pi**2 + 306.D0*vlpm**2 + 1224.D0*vdmp + 102.D
     &    0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG**2 * ( 780.D0*b*vlpm**2 - 240.D0*b
     &    *vltm + 600.D0*b*vltm**2 + 240.D0*b*vlwm + 2520.D0*b*vlwm**2
     &     + 2520.D0*b*vdw + 600.D0*b*vdt - 260.D0*b*pi**2 + 3264.D0*
     &    vlpm*vlsm + 936.D0*vlpm*vltm - 7464.D0*vlpm*vlwm - 624.D0*
     &    vlpm - 1632.D0*vlpm**2 - 6528.D0*vdmp - 544.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG**3*xnsq * ( 408.D0*b*vlpm**2 - 96.D
     &    0*b*vltm*vlwm - 48.D0*b*vltm**2 + 1584.D0*b*vlwm**2 + 1536.D0
     &    *b*vdw - 96.D0*b*vdt - 120.D0*b*pi**2 + 1248.D0*vlpm*vlsm - 
     &    2496.D0*vlpm*vlwm - 624.D0*vlpm**2 - 2496.D0*vdmp - 208.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG**3*xnsq**2 * ( 48.D0*b*vlsm*vltm
     &     - 192.D0*b*vlsm*vlwm + 36.D0*b*vlsm**2 - 192.D0*b*vdw + 48.D0
     &    *b*vdt + 12.D0*b*pi**2 + 84.D0*vlpm**2 + 336.D0*vdmp + 28.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**4*INVG**3 * ( 324.D0*b*vlpm**2 + 48.D0*b*
     &    vltm**2 + 1248.D0*b*vlwm**2 + 1248.D0*b*vdw + 48.D0*b*vdt - 
     &    108.D0*b*pi**2 + 984.D0*vlpm*vlsm + 48.D0*vlpm*vltm - 2016.D0
     &    *vlpm*vlwm - 492.D0*vlpm**2 - 1968.D0*vdmp - 164.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**4*xnsq * (  - 480.D0*vlpm )
      gcoeff8 = gcoeff8 + t1**5*INVG*xnsq * ( 144.D0*b*vltm**2 - 144.D0
     &    *b*vlwm**2 - 144.D0*b*vdw + 144.D0*b*vdt - 2016.D0*vlpm*vlsm
     &     + 960.D0*vlpm*vltm + 3072.D0*vlpm*vlwm + 2592.D0*vlpm + 1008.
     &    D0*vlpm**2 + 4032.D0*vdmp + 336.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**5*INVG * (  - 672.D0*vlpm*vltm + 672.D0*
     &    vlpm*vlwm )
      gcoeff8 = gcoeff8 + t1**5*INVG**2*xnsq * (  - 900.D0*b*vlpm**2 + 
     &    432.D0*b*vltm*vlwm - 336.D0*b*vltm + 816.D0*b*vltm**2 - 240.D0
     &    *b*vlwm - 3984.D0*b*vlwm**2 - 3768.D0*b*vdw + 1032.D0*b*vdt
     &     + 228.D0*b*pi**2 - 8328.D0*vlpm*vlsm + 936.D0*vlpm*vltm + 
     &    15720.D0*vlpm*vlwm + 2400.D0*vlpm + 4164.D0*vlpm**2 + 16656.D0
     &    *vdmp + 1388.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**5*INVG**2*xnsq**2 * ( 72.D0*b*vlsm*vltm
     &     + 216.D0*b*vlsm*vlwm - 72.D0*b*vlsm**2 + 48.D0*b*vltm - 48.D0
     &    *b*vlwm + 216.D0*b*vdw + 72.D0*b*vdt - 24.D0*b*pi**2 - 288.D0
     &    *vlpm**2 - 1152.D0*vdmp - 96.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**5*INVG**2 * (  - 648.D0*b*vlpm**2 + 96.D0
     &    *b*vltm - 912.D0*b*vltm**2 - 96.D0*b*vlwm - 1680.D0*b*vlwm**2
     &     - 1680.D0*b*vdw - 912.D0*b*vdt + 216.D0*b*pi**2 - 2880.D0*
     &    vlpm*vlsm - 2496.D0*vlpm*vltm + 8256.D0*vlpm*vlwm + 576.D0*
     &    vlpm + 1440.D0*vlpm**2 + 5760.D0*vdmp + 480.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**5*INVG**3*xnsq * (  - 828.D0*b*vlpm**2 + 
     &    384.D0*b*vltm*vlwm + 240.D0*b*vltm**2 - 3168.D0*b*vlwm**2 - 
     &    2976.D0*b*vdw + 432.D0*b*vdt + 212.D0*b*pi**2 - 3336.D0*vlpm*
     &    vlsm + 48.D0*vlpm*vltm + 6624.D0*vlpm*vlwm + 1668.D0*vlpm**2
     &     + 6672.D0*vdmp + 556.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**5*INVG**3*xnsq**2 * (  - 144.D0*b*vlsm*
     &    vltm + 288.D0*b*vlsm*vlwm - 36.D0*b*vlsm**2 + 288.D0*b*vdw - 
     &    144.D0*b*vdt - 12.D0*b*pi**2 - 228.D0*vlpm**2 - 912.D0*vdmp
     &     - 76.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**5*INVG**3 * (  - 600.D0*b*vlpm**2 - 288.D0
     &    *b*vltm**2 - 2112.D0*b*vlwm**2 - 2112.D0*b*vdw - 288.D0*b*vdt
     &     + 200.D0*b*pi**2 - 2112.D0*vlpm*vlsm - 384.D0*vlpm*vltm + 
     &    4608.D0*vlpm*vlwm + 1056.D0*vlpm**2 + 4224.D0*vdmp + 352.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG*xnsq * ( 672.D0*vlpm*vlsm - 672.D0
     &    *vlpm*vltm - 672.D0*vlpm*vlwm - 864.D0*vlpm - 336.D0*vlpm**2
     &     - 1344.D0*vdmp - 112.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG**2*xnsq * ( 300.D0*b*vlpm**2 - 144.
     &    D0*b*vltm*vlwm + 96.D0*b*vltm - 984.D0*b*vltm**2 + 96.D0*b*
     &    vlwm + 2040.D0*b*vlwm**2 + 1968.D0*b*vdw - 1056.D0*b*vdt - 76.
     &    D0*b*pi**2 + 7704.D0*vlpm*vlsm - 2496.D0*vlpm*vltm - 12912.D0
     &    *vlpm*vlwm - 2592.D0*vlpm - 3852.D0*vlpm**2 - 15408.D0*vdmp
     &     - 1284.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG**2*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm - 48.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.
     &    D0*b*vdt + 8.D0*b*pi**2 + 96.D0*vlpm**2 + 384.D0*vdmp + 32.D0
     &    *pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG**2 * ( 216.D0*b*vlpm**2 + 432.D0*b
     &    *vltm**2 + 432.D0*b*vlwm**2 + 432.D0*b*vdw + 432.D0*b*vdt - 
     &    72.D0*b*pi**2 + 960.D0*vlpm*vlsm + 2736.D0*vlpm*vltm - 4656.D0
     &    *vlpm*vlwm - 192.D0*vlpm - 480.D0*vlpm**2 - 1920.D0*vdmp - 
     &    160.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG**3*xnsq * ( 948.D0*b*vlpm**2 - 576.
     &    D0*b*vltm*vlwm - 576.D0*b*vltm**2 + 3792.D0*b*vlwm**2 + 3504.D
     &    0*b*vdw - 864.D0*b*vdt - 220.D0*b*pi**2 + 5592.D0*vlpm*vlsm
     &     - 384.D0*vlpm*vltm - 10800.D0*vlpm*vlwm - 2796.D0*vlpm**2 - 
     &    11184.D0*vdmp - 932.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG**3*xnsq**2 * ( 144.D0*b*vlsm*vltm
     &     - 192.D0*b*vlsm*vlwm + 12.D0*b*vlsm**2 - 192.D0*b*vdw + 144.D
     &    0*b*vdt + 4.D0*b*pi**2 + 300.D0*vlpm**2 + 1200.D0*vdmp + 100.D
     &    0*pi**2 )
      gcoeff8 = gcoeff8 + t1**6*INVG**3 * ( 648.D0*b*vlpm**2 + 624.D0*b
     &    *vltm**2 + 1968.D0*b*vlwm**2 + 1968.D0*b*vdw + 624.D0*b*vdt
     &     - 216.D0*b*pi**2 + 2496.D0*vlpm*vlsm + 1200.D0*vlpm*vltm - 
     &    6192.D0*vlpm*vlwm - 1248.D0*vlpm**2 - 4992.D0*vdmp - 416.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**7*INVG**2*xnsq * ( 432.D0*b*vltm**2 - 432.
     &    D0*b*vlwm**2 - 432.D0*b*vdw + 432.D0*b*vdt - 4224.D0*vlpm*
     &    vlsm + 2736.D0*vlpm*vltm + 5712.D0*vlpm*vlwm + 1536.D0*vlpm
     &     + 2112.D0*vlpm**2 + 8448.D0*vdmp + 704.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**7*INVG**2 * (  - 1056.D0*vlpm*vltm + 1056.
     &    D0*vlpm*vlwm )
      gcoeff8 = gcoeff8 + t1**7*INVG**3*xnsq * (  - 576.D0*b*vlpm**2 + 
     &    384.D0*b*vltm*vlwm + 816.D0*b*vltm**2 - 2736.D0*b*vlwm**2 - 
     &    2544.D0*b*vdw + 1008.D0*b*vdt + 128.D0*b*pi**2 - 6144.D0*vlpm
     &    *vlsm + 1200.D0*vlpm*vltm + 11088.D0*vlpm*vlwm + 3072.D0*
     &    vlpm**2 + 12288.D0*vdmp + 1024.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**7*INVG**3*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm + 48.D0*b*vlsm*vlwm + 48.D0*b*vdw - 48.D0*b*vdt - 192.D0
     &    *vlpm**2 - 768.D0*vdmp - 64.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**7*INVG**3 * (  - 384.D0*b*vlpm**2 - 576.D0
     &    *b*vltm**2 - 960.D0*b*vlwm**2 - 960.D0*b*vdw - 576.D0*b*vdt
     &     + 128.D0*b*pi**2 - 1536.D0*vlpm*vlsm - 1824.D0*vlpm*vltm + 
     &    4896.D0*vlpm*vlwm + 768.D0*vlpm**2 + 3072.D0*vdmp + 256.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t1**8*INVG**2*xnsq * ( 1056.D0*vlpm*vlsm - 
     &    1056.D0*vlpm*vltm - 1056.D0*vlpm*vlwm - 384.D0*vlpm - 528.D0*
     &    vlpm**2 - 2112.D0*vdmp - 176.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**8*INVG**3*xnsq * ( 144.D0*b*vlpm**2 - 96.D
     &    0*b*vltm*vlwm - 624.D0*b*vltm**2 + 1104.D0*b*vlwm**2 + 1056.D0
     &    *b*vdw - 672.D0*b*vdt - 32.D0*b*pi**2 + 4416.D0*vlpm*vlsm - 
     &    1824.D0*vlpm*vltm - 7008.D0*vlpm*vlwm - 2208.D0*vlpm**2 - 
     &    8832.D0*vdmp - 736.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**8*INVG**3*xnsq**2 * ( 48.D0*vlpm**2 + 192.
     &    D0*vdmp + 16.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**8*INVG**3 * ( 96.D0*b*vlpm**2 + 192.D0*b*
     &    vltm**2 + 192.D0*b*vlwm**2 + 192.D0*b*vdw + 192.D0*b*vdt - 32.
     &    D0*b*pi**2 + 384.D0*vlpm*vlsm + 1344.D0*vlpm*vltm - 2112.D0*
     &    vlpm*vlwm - 192.D0*vlpm**2 - 768.D0*vdmp - 64.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**9*INVG**3*xnsq * ( 192.D0*b*vltm**2 - 192.
     &    D0*b*vlwm**2 - 192.D0*b*vdw + 192.D0*b*vdt - 1920.D0*vlpm*
     &    vlsm + 1344.D0*vlpm*vltm + 2496.D0*vlpm*vlwm + 960.D0*vlpm**2
     &     + 3840.D0*vdmp + 320.D0*pi**2 )
      gcoeff8 = gcoeff8 + t1**9*INVG**3 * (  - 384.D0*vlpm*vltm + 384.D0
     &    *vlpm*vlwm )
      gcoeff8 = gcoeff8 + t1**10*INVG**3*xnsq * ( 384.D0*vlpm*vlsm - 
     &    384.D0*vlpm*vltm - 384.D0*vlpm*vlwm - 192.D0*vlpm**2 - 768.D0
     &    *vdmp - 64.D0*pi**2 )
      gcoeff8 = gcoeff8 + t2**(-1)*INVG*xnsq**2 * (  - 9.D0/4.D0*b*vlsm
     &    *omro**(-2) + 3.D0/4.D0*b*vlsm*omro**(-1) - 3.D0/2.D0*b*vlsm
     &     + 3.D0/2.D0*b*omro**(-1) - 3.D0/2.D0*b + 9.D0/16.D0*vlpm**2*
     &    omro**(-2) - 3.D0/8.D0*vlpm**2*omro**(-1) + 9.D0/16.D0*
     &    vlpm**2 + 9.D0/4.D0*vdmp*omro**(-2) - 3.D0/2.D0*vdmp*
     &    omro**(-1) + 9.D0/4.D0*vdmp + 3.D0/16.D0*omro**(-2)*pi**2 - 1.
     &    D0/8.D0*omro**(-1)*pi**2 + 3.D0/16.D0*pi**2 )
      gcoeff8 = gcoeff8 + t2**(-1)*INVG**2*xnsq**2 * (  - 3.D0/4.D0*b*
     &    vlsm*omro**(-1) + 3.D0/4.D0*b*vlsm + 3.D0/16.D0*vlpm**2*
     &    omro**(-1) - 3.D0/16.D0*vlpm**2 + 3.D0/4.D0*vdmp*omro**(-1)
     &     - 3.D0/4.D0*vdmp + 1.D0/16.D0*omro**(-1)*pi**2 - 1.D0/16.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + t2**(-1)*xnsq**2 * (  - 9.D0*b*vlsm*
     &    omro**(-2) + 6.D0*b*vlsm*omro**(-1) + 6.D0*b*omro**(-1) + 9.D0
     &    /4.D0*vlpm**2*omro**(-2) - 9.D0/4.D0*vlpm**2*omro**(-1) + 9.D0
     &    *vdmp*omro**(-2) - 9.D0*vdmp*omro**(-1) + 3.D0/4.D0*
     &    omro**(-2)*pi**2 - 3.D0/4.D0*omro**(-1)*pi**2 )
      gcoeff8 = gcoeff8 + ro*TBAR**(-1)*xnsq * (  - 6.D0*b*vltm )
      gcoeff8 = gcoeff8 + ro*TBAR**(-1)*xnsq**2 * ( 6.D0*b*vltm )
      gcoeff8 = gcoeff8 + ro*xnsq * ( 9.D0*b*vlpm**2 - 36.D0*b*vltm**2
     &     - 36.D0*b*vdt - 15.D0*b*pi**2 + 24.D0*vlpm*vlsm - 48.D0*vlpm
     &    *vlwm + 6.D0*vlpm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff8 = gcoeff8 + ro*xnsq**2 * ( 3.D0*vlpm**2 + 12.D0*vdmp + 
     &    pi**2 )
      gcoeff8 = gcoeff8 + ro * ( 18.D0*b*vlpm**2 + 36.D0*b*vltm**2 + 36.
     &    D0*b*vlwm**2 + 36.D0*b*vdw + 36.D0*b*vdt - 6.D0*b*pi**2 + 12.D
     &    0*vlpm*vlsm + 24.D0*vlpm*vltm - 48.D0*vlpm*vlwm - 6.D0*
     &    vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff8 = gcoeff8 + TBAR**(-1)*xnsq * ( 48.D0*b*vltm )
      gcoeff8 = gcoeff8 + TBAR**(-1)*xnsq**2 * (  - 24.D0*b*vltm )
      gcoeff8 = gcoeff8 + TBAR**(-1) * (  - 24.D0*b*vltm )
      gcoeff8 = gcoeff8 + UBAR**(-1)*xnsq * ( 24.D0*b*vlwm )
      gcoeff8 = gcoeff8 + UBAR**(-1) * (  - 24.D0*b*vlwm )
      gcoeff8 = gcoeff8 + INVG*xnsq * ( 6.D0*b*vlpm**2 - 48.D0*b*vltm
     &     + 24.D0*b*vlwm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 + 12.D0*vlpm*
     &    vlsm - 24.D0*vlpm*vlwm - 6.D0*vlpm**2 - 24.D0*vdmp - 2.D0*
     &    pi**2 )
      gcoeff8 = gcoeff8 + INVG*xnsq**2 * ( 9.D0*b*vlsm*omro**(-2) - 15.D
     &    0/2.D0*b*vlsm*omro**(-1) + 3.D0/2.D0*b*vlsm + 24.D0*b*vltm - 
     &    6.D0*b*omro**(-1) - 9.D0/4.D0*vlpm**2*omro**(-2) + 21.D0/8.D0
     &    *vlpm**2*omro**(-1) - 3.D0/8.D0*vlpm**2 - 9.D0*vdmp*
     &    omro**(-2) + 21.D0/2.D0*vdmp*omro**(-1) - 3.D0/2.D0*vdmp - 3.D
     &    0/4.D0*omro**(-2)*pi**2 + 7.D0/8.D0*omro**(-1)*pi**2 - 1.D0/8.
     &    D0*pi**2 )
      gcoeff8 = gcoeff8 + INVG * ( 6.D0*b*vlpm**2 + 24.D0*b*vltm + 24.D0
     &    *b*vlwm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 + 12.D0*vlpm*vlsm - 
     &    24.D0*vlpm*vlwm - 6.D0*vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff8 = gcoeff8 + INVG**2*xnsq**2 * ( 3.D0*b*vlsm*omro**(-1) - 
     &    3.D0/4.D0*vlpm**2*omro**(-1) - 3.D0*vdmp*omro**(-1) - 1.D0/4.D
     &    0*omro**(-1)*pi**2 )
      gcoeff8 = gcoeff8 + xnsq * (  - 12.D0*b*vlpm**2 + 96.D0*b*vltm*
     &    vlwm - 120.D0*b*vltm + 48.D0*b*vltm**2 + 48.D0*b*vdw + 96.D0*
     &    b*vdt - 12.D0*b*pi**2 + 48.D0*b - 24.D0*vlpm )
      gcoeff8 = gcoeff8 + xnsq**2 * (  - 48.D0*b*vlsm*vltm + 18.D0*b*
     &    vlsm*omro**(-2) - 12.D0*b*vlsm*omro**(-1) - 12.D0*b*vlsm + 12.
     &    D0*b*vlsm**2 + 72.D0*b*vltm - 48.D0*b*vdt - 12.D0*b*
     &    omro**(-1) + 4.D0*b*pi**2 - 12.D0*b - 9.D0/2.D0*vlpm**2*
     &    omro**(-2) + 9.D0/2.D0*vlpm**2*omro**(-1) - 21.D0/2.D0*
     &    vlpm**2 - 18.D0*vdmp*omro**(-2) + 18.D0*vdmp*omro**(-1) - 42.D
     &    0*vdmp - 3.D0/2.D0*omro**(-2)*pi**2 + 3.D0/2.D0*omro**(-1)*
     &    pi**2 - 7.D0/2.D0*pi**2 )
      gcoeff8 = gcoeff8 - 24.D0*b*vlpm**2 + 48.D0*b*vltm - 48.D0*b*
     &    vltm**2 + 48.D0*b*vlwm - 48.D0*b*vlwm**2 - 48.D0*b*vdw - 48.D0
     &    *b*vdt + 8.D0*b*pi**2 - 24.D0*b + 24.D0*vlpm*vlsm - 48.D0*
     &    vlpm*vltm - 24.D0*vlpm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*
     &    pi**2

      gcoeff9 =  + t1**(-1)*UBAR**(-1)*INVG*xnsq * ( 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*UBAR**(-1)*INVG*xnsq**2 * (  - 24.D0
     &    *b*vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*UBAR**(-1)*INVG * (  - 24.D0*b*vlwm
     &     )
      gcoeff9 = gcoeff9 + t1**(-1)*UBAR**(-1)*xnsq * ( 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*UBAR**(-1)*xnsq**2 * (  - 24.D0*b*
     &    vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*UBAR**(-1) * (  - 24.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*INVG*xnsq * (  - 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*INVG*xnsq**2 * ( 9.D0/4.D0*b*vlsm*
     &    omro**(-2) - 15.D0/4.D0*b*vlsm*omro**(-1) + 9.D0/2.D0*b*vlsm
     &     + 24.D0*b*vlwm - 3.D0/2.D0*b*omro**(-1) + 3.D0/2.D0*b - 9.D0/
     &    16.D0*vlpm**2*omro**(-2) + 9.D0/8.D0*vlpm**2*omro**(-1) - 21.D
     &    0/16.D0*vlpm**2 - 9.D0/4.D0*vdmp*omro**(-2) + 9.D0/2.D0*vdmp*
     &    omro**(-1) - 21.D0/4.D0*vdmp - 3.D0/16.D0*omro**(-2)*pi**2 + 
     &    3.D0/8.D0*omro**(-1)*pi**2 - 7.D0/16.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**(-1)*INVG * ( 24.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1**(-1)*INVG**2*xnsq**2 * ( 3.D0/4.D0*b*vlsm
     &    *omro**(-1) - 3.D0/4.D0*b*vlsm - 3.D0/16.D0*vlpm**2*
     &    omro**(-1) + 3.D0/16.D0*vlpm**2 - 3.D0/4.D0*vdmp*omro**(-1)
     &     + 3.D0/4.D0*vdmp - 1.D0/16.D0*omro**(-1)*pi**2 + 1.D0/16.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**(-1)*xnsq**2 * ( 9.D0*b*vlsm*omro**(-2)
     &     - 18.D0*b*vlsm*omro**(-1) - 6.D0*b*omro**(-1) - 9.D0/4.D0*
     &    vlpm**2*omro**(-2) + 21.D0/4.D0*vlpm**2*omro**(-1) - 9.D0*
     &    vdmp*omro**(-2) + 21.D0*vdmp*omro**(-1) - 3.D0/4.D0*
     &    omro**(-2)*pi**2 + 7.D0/4.D0*omro**(-1)*pi**2 )
      gcoeff9 = gcoeff9 + t1*ro*xnsq * ( 36.D0*b*vltm**2 - 36.D0*b*
     &    vlwm**2 - 36.D0*b*vdw + 36.D0*b*vdt - 72.D0*vlpm*vlsm + 24.D0
     &    *vlpm*vltm + 120.D0*vlpm*vlwm - 24.D0*vlpm + 36.D0*vlpm**2 + 
     &    144.D0*vdmp + 12.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*ro * (  - 72.D0*vlpm*vltm + 72.D0*vlpm*
     &    vlwm )
      gcoeff9 = gcoeff9 + t1*UBAR**(-1)*xnsq * ( 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1*UBAR**(-1)*xnsq**2 * (  - 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + t1*INVG*xnsq * (  - 18.D0*b*vlpm**2 - 96.D0*b
     &    *vltm*vlwm + 144.D0*b*vltm - 48.D0*b*vltm**2 - 72.D0*b*vlwm
     &     - 120.D0*b*vlwm**2 - 168.D0*b*vdw - 96.D0*b*vdt + 22.D0*b*
     &    pi**2 - 84.D0*vlpm*vlsm + 168.D0*vlpm*vlwm + 24.D0*vlpm + 42.D
     &    0*vlpm**2 + 168.D0*vdmp + 14.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*INVG*xnsq**2 * ( 48.D0*b*vlsm*vltm + 24.D0
     &    *b*vlsm*vlwm - 18.D0*b*vlsm**2 - 96.D0*b*vltm + 24.D0*b*vlwm
     &     + 24.D0*b*vdw + 48.D0*b*vdt - 6.D0*b*pi**2 + 6.D0*vlpm**2 + 
     &    24.D0*vdmp + 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*INVG * (  - 48.D0*b*vltm + 48.D0*b*vltm**2
     &     - 48.D0*b*vlwm**2 - 48.D0*b*vdw + 48.D0*b*vdt - 96.D0*vlpm*
     &    vlsm + 48.D0*vlpm*vltm + 144.D0*vlpm*vlwm + 24.D0*vlpm + 48.D0
     &    *vlpm**2 + 192.D0*vdmp + 16.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*INVG**2*xnsq * (  - 6.D0*b*vlpm**2 - 24.D0
     &    *b*vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 
     &    24.D0*vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*INVG**2 * (  - 6.D0*b*vlpm**2 - 24.D0*b*
     &    vlwm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 - 12.D0*vlpm*vlsm + 24.D0
     &    *vlpm*vlwm + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*xnsq * ( 36.D0*b*vlpm**2 - 240.D0*b*vltm*
     &    vlwm + 96.D0*b*vltm - 168.D0*b*vltm**2 - 96.D0*b*vlwm + 72.D0
     &    *b*vlwm**2 - 48.D0*b*vdw - 288.D0*b*vdt + 28.D0*b*pi**2 + 96.D
     &    0*vlpm*vlsm - 48.D0*vlpm*vltm - 144.D0*vlpm*vlwm + 96.D0*vlpm
     &     - 48.D0*vlpm**2 - 192.D0*vdmp - 16.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1*xnsq**2 * ( 96.D0*b*vlsm*vltm - 72.D0*b*
     &    vlsm*omro**(-2) + 120.D0*b*vlsm*omro**(-1) + 48.D0*b*vlsm - 
     &    24.D0*b*vlsm**2 - 144.D0*b*vltm - 48.D0*b*vlwm + 96.D0*b*vdt
     &     + 48.D0*b*omro**(-1) - 8.D0*b*pi**2 - 48.D0*b + 18.D0*
     &    vlpm**2*omro**(-2) - 36.D0*vlpm**2*omro**(-1) + 54.D0*vlpm**2
     &     + 72.D0*vdmp*omro**(-2) - 144.D0*vdmp*omro**(-1) + 216.D0*
     &    vdmp + 6.D0*omro**(-2)*pi**2 - 12.D0*omro**(-1)*pi**2 + 18.D0
     &    *pi**2 )
      gcoeff9 = gcoeff9 + t1 * ( 72.D0*b*vlpm**2 + 144.D0*b*vltm**2 + 
     &    144.D0*b*vlwm**2 + 144.D0*b*vdw + 144.D0*b*vdt - 24.D0*b*
     &    pi**2 - 48.D0*vlpm*vlsm + 240.D0*vlpm*vltm - 144.D0*vlpm*vlwm
     &     + 96.D0*vlpm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*ro*xnsq * ( 72.D0*vlpm*vlsm - 72.D0*
     &    vlpm*vltm - 72.D0*vlpm*vlwm + 24.D0*vlpm - 36.D0*vlpm**2 - 
     &    144.D0*vdmp - 12.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*INVG*xnsq * ( 42.D0*b*vlpm**2 + 432.D0*
     &    b*vltm*vlwm - 144.D0*b*vltm + 264.D0*b*vltm**2 + 144.D0*b*
     &    vlwm + 336.D0*b*vlwm**2 + 552.D0*b*vdw + 480.D0*b*vdt - 86.D0
     &    *b*pi**2 + 48.D0*b + 252.D0*vlpm*vlsm + 48.D0*vlpm*vltm - 552.
     &    D0*vlpm*vlwm - 168.D0*vlpm - 126.D0*vlpm**2 - 504.D0*vdmp - 
     &    42.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*INVG*xnsq**2 * (  - 216.D0*b*vlsm*vltm
     &     - 144.D0*b*vlsm*vlwm - 102.D0*b*vlsm + 90.D0*b*vlsm**2 + 192.
     &    D0*b*vltm + 48.D0*b*vlwm - 144.D0*b*vdw - 216.D0*b*vdt + 30.D0
     &    *b*pi**2 - 24.D0*b - 33.D0/2.D0*vlpm**2 - 66.D0*vdmp - 11.D0/
     &    2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*INVG * (  - 36.D0*b*vlpm**2 - 216.D0*b*
     &    vltm**2 + 96.D0*b*vlwm + 72.D0*b*vlwm**2 + 72.D0*b*vdw - 216.D
     &    0*b*vdt + 12.D0*b*pi**2 - 24.D0*b + 360.D0*vlpm*vlsm - 312.D0
     &    *vlpm*vltm - 408.D0*vlpm*vlwm - 168.D0*vlpm - 180.D0*vlpm**2
     &     - 720.D0*vdmp - 60.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*INVG**2*xnsq * ( 48.D0*b*vlpm**2 + 192.D
     &    0*b*vlwm**2 + 192.D0*b*vdw - 16.D0*b*pi**2 + 120.D0*vlpm*vlsm
     &     - 240.D0*vlpm*vlwm - 60.D0*vlpm**2 - 240.D0*vdmp - 20.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**2*INVG**2*xnsq**2 * (  - 24.D0*b*vlsm*
     &    vlwm + 6.D0*b*vlsm**2 - 24.D0*b*vdw + 2.D0*b*pi**2 + 6.D0*
     &    vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*INVG**2 * ( 42.D0*b*vlpm**2 + 168.D0*b*
     &    vlwm**2 + 168.D0*b*vdw - 14.D0*b*pi**2 + 108.D0*vlpm*vlsm - 
     &    216.D0*vlpm*vlwm - 54.D0*vlpm**2 - 216.D0*vdmp - 18.D0*pi**2
     &     )
      gcoeff9 = gcoeff9 + t1**2*xnsq * (  - 24.D0*b*vlpm**2 + 192.D0*b*
     &    vltm*vlwm + 240.D0*b*vltm**2 - 144.D0*b*vlwm**2 - 48.D0*b*vdw
     &     + 336.D0*b*vdt - 24.D0*b*pi**2 + 48.D0*b - 264.D0*vlpm*vlsm
     &     + 240.D0*vlpm*vltm + 288.D0*vlpm*vlwm - 264.D0*vlpm + 132.D0
     &    *vlpm**2 + 528.D0*vdmp + 44.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2*xnsq**2 * (  - 48.D0*b*vlsm*vltm - 48.D0
     &    *b*vlsm*vlwm - 48.D0*b*vlsm + 24.D0*b*vlsm**2 + 48.D0*b*vltm
     &     + 48.D0*b*vlwm - 48.D0*b*vdw - 48.D0*b*vdt + 8.D0*b*pi**2 - 
     &    42.D0*vlpm**2 - 168.D0*vdmp - 14.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**2 * (  - 48.D0*b*vlpm**2 - 96.D0*b*
     &    vltm**2 - 96.D0*b*vlwm**2 - 96.D0*b*vdw - 96.D0*b*vdt + 16.D0
     &    *b*pi**2 + 144.D0*vlpm*vlsm - 432.D0*vlpm*vltm + 144.D0*vlpm*
     &    vlwm - 192.D0*vlpm - 72.D0*vlpm**2 - 288.D0*vdmp - 24.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG*xnsq * (  - 120.D0*b*vlpm**2 - 480.
     &    D0*b*vltm*vlwm - 48.D0*b*vltm - 456.D0*b*vltm**2 - 144.D0*b*
     &    vlwm - 504.D0*b*vlwm**2 - 744.D0*b*vdw - 696.D0*b*vdt + 120.D0
     &    *b*pi**2 - 48.D0*b - 576.D0*vlpm*vlsm - 312.D0*vlpm*vltm + 
     &    1464.D0*vlpm*vlwm + 648.D0*vlpm + 288.D0*vlpm**2 + 1152.D0*
     &    vdmp + 96.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG*xnsq**2 * ( 264.D0*b*vlsm*vltm + 
     &    216.D0*b*vlsm*vlwm + 144.D0*b*vlsm - 120.D0*b*vlsm**2 - 144.D0
     &    *b*vltm - 144.D0*b*vlwm + 216.D0*b*vdw + 264.D0*b*vdt - 40.D0
     &    *b*pi**2 + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG * ( 96.D0*b*vltm + 192.D0*b*
     &    vltm**2 - 96.D0*b*vlwm - 192.D0*b*vlwm**2 - 192.D0*b*vdw + 
     &    192.D0*b*vdt - 864.D0*vlpm*vlsm + 672.D0*vlpm*vltm + 1056.D0*
     &    vlpm*vlwm + 480.D0*vlpm + 432.D0*vlpm**2 + 1728.D0*vdmp + 144.
     &    D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG**2*xnsq * (  - 204.D0*b*vlpm**2 - 
     &    96.D0*b*vltm*vlwm - 96.D0*b*vltm - 48.D0*b*vltm**2 - 48.D0*b*
     &    vlwm - 864.D0*b*vlwm**2 - 912.D0*b*vdw - 96.D0*b*vdt + 84.D0*
     &    b*pi**2 - 648.D0*vlpm*vlsm + 1296.D0*vlpm*vlwm + 48.D0*vlpm
     &     + 324.D0*vlpm**2 + 1296.D0*vdmp + 108.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG**2*xnsq**2 * ( 72.D0*b*vlsm*vltm
     &     + 216.D0*b*vlsm*vlwm - 72.D0*b*vlsm**2 + 48.D0*b*vltm + 216.D
     &    0*b*vdw + 72.D0*b*vdt - 24.D0*b*pi**2 - 48.D0*vlpm**2 - 192.D0
     &    *vdmp - 16.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG**2 * (  - 156.D0*b*vlpm**2 + 48.D0
     &    *b*vltm + 24.D0*b*vltm**2 - 48.D0*b*vlwm - 648.D0*b*vlwm**2
     &     - 648.D0*b*vdw + 24.D0*b*vdt + 52.D0*b*pi**2 - 552.D0*vlpm*
     &    vlsm + 24.D0*vlpm*vltm + 1080.D0*vlpm*vlwm + 48.D0*vlpm + 276.
     &    D0*vlpm**2 + 1104.D0*vdmp + 92.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG**3*xnsq * (  - 12.D0*b*vlpm**2 - 
     &    48.D0*b*vlwm**2 - 48.D0*b*vdw + 4.D0*b*pi**2 - 24.D0*vlpm*
     &    vlsm + 48.D0*vlpm*vlwm + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**3*INVG**3 * (  - 12.D0*b*vlpm**2 - 48.D0*
     &    b*vlwm**2 - 48.D0*b*vdw + 4.D0*b*pi**2 - 24.D0*vlpm*vlsm + 48.
     &    D0*vlpm*vlwm + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**3*xnsq * (  - 96.D0*b*vltm**2 + 96.D0*b*
     &    vlwm**2 + 96.D0*b*vdw - 96.D0*b*vdt + 288.D0*vlpm*vlsm - 432.D
     &    0*vlpm*vltm - 144.D0*vlpm*vlwm + 480.D0*vlpm - 144.D0*vlpm**2
     &     - 576.D0*vdmp - 48.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG*xnsq * ( 132.D0*b*vlpm**2 + 144.D0
     &    *b*vltm*vlwm + 96.D0*b*vltm + 264.D0*b*vltm**2 + 96.D0*b*vlwm
     &     + 408.D0*b*vlwm**2 + 480.D0*b*vdw + 336.D0*b*vdt - 68.D0*b*
     &    pi**2 + 48.D0*b + 912.D0*vlpm*vlsm + 672.D0*vlpm*vltm - 2496.D
     &    0*vlpm*vlwm - 1416.D0*vlpm - 456.D0*vlpm**2 - 1824.D0*vdmp - 
     &    152.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG*xnsq**2 * (  - 96.D0*b*vlsm*vltm
     &     - 96.D0*b*vlsm*vlwm - 48.D0*b*vlsm + 48.D0*b*vlsm**2 + 48.D0
     &    *b*vltm + 48.D0*b*vlwm - 96.D0*b*vdw - 96.D0*b*vdt + 16.D0*b*
     &    pi**2 + 6.D0*vlpm**2 + 24.D0*vdmp + 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG * ( 72.D0*b*vlpm**2 + 144.D0*b*
     &    vltm**2 + 144.D0*b*vlwm**2 + 144.D0*b*vdw + 144.D0*b*vdt - 24.
     &    D0*b*pi**2 + 720.D0*vlpm*vlsm - 1440.D0*vlpm*vlwm - 384.D0*
     &    vlpm - 360.D0*vlpm**2 - 1440.D0*vdmp - 120.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG**2*xnsq * ( 498.D0*b*vlpm**2 + 48.D
     &    0*b*vltm*vlwm + 240.D0*b*vltm + 48.D0*b*vltm**2 + 96.D0*b*
     &    vlwm + 1992.D0*b*vlwm**2 + 2016.D0*b*vdw + 72.D0*b*vdt - 174.D
     &    0*b*pi**2 + 2076.D0*vlpm*vlsm + 24.D0*vlpm*vltm - 4176.D0*
     &    vlpm*vlwm - 336.D0*vlpm - 1038.D0*vlpm**2 - 4152.D0*vdmp - 
     &    346.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG**2*xnsq**2 * (  - 144.D0*b*vlsm*
     &    vltm - 408.D0*b*vlsm*vlwm + 138.D0*b*vlsm**2 - 96.D0*b*vltm
     &     + 48.D0*b*vlwm - 408.D0*b*vdw - 144.D0*b*vdt + 46.D0*b*pi**2
     &     + 162.D0*vlpm**2 + 648.D0*vdmp + 54.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG**2 * ( 348.D0*b*vlpm**2 - 144.D0*b
     &    *vltm + 72.D0*b*vltm**2 + 144.D0*b*vlwm + 1320.D0*b*vlwm**2
     &     + 1320.D0*b*vdw + 72.D0*b*vdt - 116.D0*b*pi**2 + 1536.D0*
     &    vlpm*vlsm + 24.D0*vlpm*vltm - 3096.D0*vlpm*vlwm - 240.D0*vlpm
     &     - 768.D0*vlpm**2 - 3072.D0*vdmp - 256.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG**3*xnsq * ( 96.D0*b*vlpm**2 + 384.D
     &    0*b*vlwm**2 + 384.D0*b*vdw - 32.D0*b*pi**2 + 240.D0*vlpm*vlsm
     &     - 480.D0*vlpm*vlwm - 120.D0*vlpm**2 - 480.D0*vdmp - 40.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG**3*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vlwm + 12.D0*b*vlsm**2 - 48.D0*b*vdw + 4.D0*b*pi**2 + 12.D0*
     &    vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**4*INVG**3 * ( 84.D0*b*vlpm**2 + 336.D0*b*
     &    vlwm**2 + 336.D0*b*vdw - 28.D0*b*pi**2 + 216.D0*vlpm*vlsm - 
     &    432.D0*vlpm*vlwm - 108.D0*vlpm**2 - 432.D0*vdmp - 36.D0*pi**2
     &     )
      gcoeff9 = gcoeff9 + t1**4*xnsq * (  - 480.D0*vlpm )
      gcoeff9 = gcoeff9 + t1**5*INVG*xnsq * ( 144.D0*b*vltm**2 - 144.D0
     &    *b*vlwm**2 - 144.D0*b*vdw + 144.D0*b*vdt - 1056.D0*vlpm*vlsm
     &     + 2112.D0*vlpm*vlwm + 1728.D0*vlpm + 528.D0*vlpm**2 + 2112.D0
     &    *vdmp + 176.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**5*INVG * (  - 672.D0*vlpm*vltm + 672.D0*
     &    vlpm*vlwm )
      gcoeff9 = gcoeff9 + t1**5*INVG**2*xnsq * (  - 624.D0*b*vlpm**2 + 
     &    192.D0*b*vltm*vlwm - 240.D0*b*vltm + 168.D0*b*vltm**2 - 144.D0
     &    *b*vlwm - 2472.D0*b*vlwm**2 - 2376.D0*b*vdw + 264.D0*b*vdt + 
     &    176.D0*b*pi**2 - 3936.D0*vlpm*vlsm + 24.D0*vlpm*vltm + 7848.D0
     &    *vlpm*vlwm + 960.D0*vlpm + 1968.D0*vlpm**2 + 7872.D0*vdmp + 
     &    656.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**5*INVG**2*xnsq**2 * ( 120.D0*b*vlsm*vltm
     &     + 264.D0*b*vlsm*vlwm - 96.D0*b*vlsm**2 + 48.D0*b*vltm - 48.D0
     &    *b*vlwm + 264.D0*b*vdw + 120.D0*b*vdt - 32.D0*b*pi**2 - 216.D0
     &    *vlpm**2 - 864.D0*vdmp - 72.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**5*INVG**2 * (  - 432.D0*b*vlpm**2 + 96.D0
     &    *b*vltm - 480.D0*b*vltm**2 - 96.D0*b*vlwm - 1248.D0*b*vlwm**2
     &     - 1248.D0*b*vdw - 480.D0*b*vdt + 144.D0*b*pi**2 - 2016.D0*
     &    vlpm*vlsm - 720.D0*vlpm*vltm + 4752.D0*vlpm*vlwm + 384.D0*
     &    vlpm + 1008.D0*vlpm**2 + 4032.D0*vdmp + 336.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**5*INVG**3*xnsq * (  - 312.D0*b*vlpm**2 + 
     &    96.D0*b*vltm*vlwm + 48.D0*b*vltm**2 - 1200.D0*b*vlwm**2 - 
     &    1152.D0*b*vdw + 96.D0*b*vdt + 88.D0*b*pi**2 - 1008.D0*vlpm*
     &    vlsm + 2016.D0*vlpm*vlwm + 504.D0*vlpm**2 + 2016.D0*vdmp + 
     &    168.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**5*INVG**3*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm + 144.D0*b*vlsm*vlwm - 24.D0*b*vlsm**2 + 144.D0*b*vdw - 
     &    48.D0*b*vdt - 8.D0*b*pi**2 - 72.D0*vlpm**2 - 288.D0*vdmp - 24.
     &    D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**5*INVG**3 * (  - 240.D0*b*vlpm**2 - 48.D0
     &    *b*vltm**2 - 912.D0*b*vlwm**2 - 912.D0*b*vdw - 48.D0*b*vdt + 
     &    80.D0*b*pi**2 - 768.D0*vlpm*vlsm - 48.D0*vlpm*vltm + 1584.D0*
     &    vlpm*vlwm + 384.D0*vlpm**2 + 1536.D0*vdmp + 128.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG*xnsq * ( 672.D0*vlpm*vlsm - 672.D0
     &    *vlpm*vltm - 672.D0*vlpm*vlwm - 864.D0*vlpm - 336.D0*vlpm**2
     &     - 1344.D0*vdmp - 112.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG**2*xnsq * ( 300.D0*b*vlpm**2 - 144.
     &    D0*b*vltm*vlwm + 96.D0*b*vltm - 552.D0*b*vltm**2 + 96.D0*b*
     &    vlwm + 1608.D0*b*vlwm**2 + 1536.D0*b*vdw - 624.D0*b*vdt - 76.D
     &    0*b*pi**2 + 4536.D0*vlpm*vlsm - 720.D0*vlpm*vltm - 8352.D0*
     &    vlpm*vlwm - 1440.D0*vlpm - 2268.D0*vlpm**2 - 9072.D0*vdmp - 
     &    756.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG**2*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm - 48.D0*b*vlsm*vlwm + 24.D0*b*vlsm**2 - 48.D0*b*vdw - 48.
     &    D0*b*vdt + 8.D0*b*pi**2 + 96.D0*vlpm**2 + 384.D0*vdmp + 32.D0
     &    *pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG**2 * ( 216.D0*b*vlpm**2 + 432.D0*b
     &    *vltm**2 + 432.D0*b*vlwm**2 + 432.D0*b*vdw + 432.D0*b*vdt - 
     &    72.D0*b*pi**2 + 960.D0*vlpm*vlsm + 1680.D0*vlpm*vltm - 3600.D0
     &    *vlpm*vlwm - 192.D0*vlpm - 480.D0*vlpm**2 - 1920.D0*vdmp - 
     &    160.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG**3*xnsq * ( 516.D0*b*vlpm**2 - 288.
     &    D0*b*vltm*vlwm - 192.D0*b*vltm**2 + 1968.D0*b*vlwm**2 + 1824.D
     &    0*b*vdw - 336.D0*b*vdt - 124.D0*b*pi**2 + 2328.D0*vlpm*vlsm
     &     - 48.D0*vlpm*vltm - 4608.D0*vlpm*vlwm - 1164.D0*vlpm**2 - 
     &    4656.D0*vdmp - 388.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG**3*xnsq**2 * ( 96.D0*b*vlsm*vltm
     &     - 144.D0*b*vlsm*vlwm + 12.D0*b*vlsm**2 - 144.D0*b*vdw + 96.D0
     &    *b*vdt + 4.D0*b*pi**2 + 156.D0*vlpm**2 + 624.D0*vdmp + 52.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**6*INVG**3 * ( 360.D0*b*vlpm**2 + 240.D0*b
     &    *vltm**2 + 1200.D0*b*vlwm**2 + 1200.D0*b*vdw + 240.D0*b*vdt
     &     - 120.D0*b*pi**2 + 1344.D0*vlpm*vlsm + 336.D0*vlpm*vltm - 
     &    3024.D0*vlpm*vlwm - 672.D0*vlpm**2 - 2688.D0*vdmp - 224.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**7*INVG**2*xnsq * ( 432.D0*b*vltm**2 - 432.
     &    D0*b*vlwm**2 - 432.D0*b*vdw + 432.D0*b*vdt - 3168.D0*vlpm*
     &    vlsm + 1680.D0*vlpm*vltm + 4656.D0*vlpm*vlwm + 1152.D0*vlpm
     &     + 1584.D0*vlpm**2 + 6336.D0*vdmp + 528.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**7*INVG**2 * (  - 1056.D0*vlpm*vltm + 1056.
     &    D0*vlpm*vlwm )
      gcoeff9 = gcoeff9 + t1**7*INVG**3*xnsq * (  - 432.D0*b*vlpm**2 + 
     &    288.D0*b*vltm*vlwm + 384.D0*b*vltm**2 - 1824.D0*b*vlwm**2 - 
     &    1680.D0*b*vdw + 528.D0*b*vdt + 96.D0*b*pi**2 - 3264.D0*vlpm*
     &    vlsm + 336.D0*vlpm*vltm + 6192.D0*vlpm*vlwm + 1632.D0*vlpm**2
     &     + 6528.D0*vdmp + 544.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**7*INVG**3*xnsq**2 * (  - 48.D0*b*vlsm*
     &    vltm + 48.D0*b*vlsm*vlwm + 48.D0*b*vdw - 48.D0*b*vdt - 144.D0
     &    *vlpm**2 - 576.D0*vdmp - 48.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**7*INVG**3 * (  - 288.D0*b*vlpm**2 - 384.D0
     &    *b*vltm**2 - 768.D0*b*vlwm**2 - 768.D0*b*vdw - 384.D0*b*vdt
     &     + 96.D0*b*pi**2 - 1152.D0*vlpm*vlsm - 864.D0*vlpm*vltm + 
     &    3168.D0*vlpm*vlwm + 576.D0*vlpm**2 + 2304.D0*vdmp + 192.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t1**8*INVG**2*xnsq * ( 1056.D0*vlpm*vlsm - 
     &    1056.D0*vlpm*vltm - 1056.D0*vlpm*vlwm - 384.D0*vlpm - 528.D0*
     &    vlpm**2 - 2112.D0*vdmp - 176.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**8*INVG**3*xnsq * ( 144.D0*b*vlpm**2 - 96.D
     &    0*b*vltm*vlwm - 432.D0*b*vltm**2 + 912.D0*b*vlwm**2 + 864.D0*
     &    b*vdw - 480.D0*b*vdt - 32.D0*b*pi**2 + 2880.D0*vlpm*vlsm - 
     &    864.D0*vlpm*vltm - 4896.D0*vlpm*vlwm - 1440.D0*vlpm**2 - 5760.
     &    D0*vdmp - 480.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**8*INVG**3*xnsq**2 * ( 48.D0*vlpm**2 + 192.
     &    D0*vdmp + 16.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**8*INVG**3 * ( 96.D0*b*vlpm**2 + 192.D0*b*
     &    vltm**2 + 192.D0*b*vlwm**2 + 192.D0*b*vdw + 192.D0*b*vdt - 32.
     &    D0*b*pi**2 + 384.D0*vlpm*vlsm + 960.D0*vlpm*vltm - 1728.D0*
     &    vlpm*vlwm - 192.D0*vlpm**2 - 768.D0*vdmp - 64.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**9*INVG**3*xnsq * ( 192.D0*b*vltm**2 - 192.
     &    D0*b*vlwm**2 - 192.D0*b*vdw + 192.D0*b*vdt - 1536.D0*vlpm*
     &    vlsm + 960.D0*vlpm*vltm + 2112.D0*vlpm*vlwm + 768.D0*vlpm**2
     &     + 3072.D0*vdmp + 256.D0*pi**2 )
      gcoeff9 = gcoeff9 + t1**9*INVG**3 * (  - 384.D0*vlpm*vltm + 384.D0
     &    *vlpm*vlwm )
      gcoeff9 = gcoeff9 + t1**10*INVG**3*xnsq * ( 384.D0*vlpm*vlsm - 
     &    384.D0*vlpm*vltm - 384.D0*vlpm*vlwm - 192.D0*vlpm**2 - 768.D0
     &    *vdmp - 64.D0*pi**2 )
      gcoeff9 = gcoeff9 + t2**(-1)*INVG*xnsq**2 * ( 9.D0/4.D0*b*vlsm*
     &    omro**(-2) - 9.D0/4.D0*b*vlsm*omro**(-1) + 3.D0*b*vlsm - 3.D0/
     &    2.D0*b*omro**(-1) + 3.D0/2.D0*b - 9.D0/16.D0*vlpm**2*
     &    omro**(-2) + 3.D0/4.D0*vlpm**2*omro**(-1) - 15.D0/16.D0*
     &    vlpm**2 - 9.D0/4.D0*vdmp*omro**(-2) + 3.D0*vdmp*omro**(-1) - 
     &    15.D0/4.D0*vdmp - 3.D0/16.D0*omro**(-2)*pi**2 + 1.D0/4.D0*
     &    omro**(-1)*pi**2 - 5.D0/16.D0*pi**2 )
      gcoeff9 = gcoeff9 + t2**(-1)*INVG**2*xnsq**2 * ( 3.D0/4.D0*b*vlsm
     &    *omro**(-1) - 3.D0/4.D0*b*vlsm - 3.D0/16.D0*vlpm**2*
     &    omro**(-1) + 3.D0/16.D0*vlpm**2 - 3.D0/4.D0*vdmp*omro**(-1)
     &     + 3.D0/4.D0*vdmp - 1.D0/16.D0*omro**(-1)*pi**2 + 1.D0/16.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + t2**(-1)*xnsq * (  - 48.D0*b*vlwm - 48.D0*b )
      gcoeff9 = gcoeff9 + t2**(-1)*xnsq**2 * ( 9.D0*b*vlsm*omro**(-2)
     &     - 12.D0*b*vlsm*omro**(-1) + 48.D0*b*vlwm - 6.D0*b*omro**(-1)
     &     + 24.D0*b - 9.D0/4.D0*vlpm**2*omro**(-2) + 15.D0/4.D0*
     &    vlpm**2*omro**(-1) - 9.D0*vdmp*omro**(-2) + 15.D0*vdmp*
     &    omro**(-1) - 3.D0/4.D0*omro**(-2)*pi**2 + 5.D0/4.D0*
     &    omro**(-1)*pi**2 )
      gcoeff9 = gcoeff9 + t2**(-1) * ( 24.D0*b )
      gcoeff9 = gcoeff9 + ro*xnsq * ( 9.D0*b*vlpm**2 - 36.D0*b*vltm**2
     &     - 36.D0*b*vdt - 15.D0*b*pi**2 + 24.D0*vlpm*vlsm - 48.D0*vlpm
     &    *vlwm + 6.D0*vlpm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff9 = gcoeff9 + ro*xnsq**2 * ( 3.D0*vlpm**2 + 12.D0*vdmp + 
     &    pi**2 )
      gcoeff9 = gcoeff9 + ro * ( 18.D0*b*vlpm**2 + 36.D0*b*vltm**2 + 36.
     &    D0*b*vlwm**2 + 36.D0*b*vdw + 36.D0*b*vdt - 6.D0*b*pi**2 + 12.D
     &    0*vlpm*vlsm + 24.D0*vlpm*vltm - 48.D0*vlpm*vlwm - 6.D0*
     &    vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + UBAR**(-1)*INVG*xnsq * (  - 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + UBAR**(-1)*INVG*xnsq**2 * ( 24.D0*b*vlwm )
      gcoeff9 = gcoeff9 + UBAR**(-1)*INVG * ( 24.D0*b*vlwm )
      gcoeff9 = gcoeff9 + UBAR**(-1)*xnsq * ( 72.D0*b*vlwm )
      gcoeff9 = gcoeff9 + UBAR**(-1)*xnsq**2 * (  - 24.D0*b*vlwm )
      gcoeff9 = gcoeff9 + UBAR**(-1) * (  - 48.D0*b*vlwm )
      gcoeff9 = gcoeff9 + INVG*xnsq * ( 6.D0*b*vlpm**2 - 48.D0*b*vlwm
     &     + 24.D0*b*vlwm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 + 12.D0*vlpm*
     &    vlsm - 24.D0*vlpm*vlwm - 6.D0*vlpm**2 - 24.D0*vdmp - 2.D0*
     &    pi**2 )
      gcoeff9 = gcoeff9 + INVG*xnsq**2 * (  - 9.D0*b*vlsm*omro**(-2) + 
     &    33.D0/2.D0*b*vlsm*omro**(-1) + 3.D0/2.D0*b*vlsm + 24.D0*b*
     &    vlwm + 6.D0*b*omro**(-1) + 9.D0/4.D0*vlpm**2*omro**(-2) - 39.D
     &    0/8.D0*vlpm**2*omro**(-1) - 3.D0/8.D0*vlpm**2 + 9.D0*vdmp*
     &    omro**(-2) - 39.D0/2.D0*vdmp*omro**(-1) - 3.D0/2.D0*vdmp + 3.D
     &    0/4.D0*omro**(-2)*pi**2 - 13.D0/8.D0*omro**(-1)*pi**2 - 1.D0/
     &    8.D0*pi**2 )
      gcoeff9 = gcoeff9 + INVG * ( 6.D0*b*vlpm**2 + 24.D0*b*vlwm + 24.D0
     &    *b*vlwm**2 + 24.D0*b*vdw - 2.D0*b*pi**2 + 12.D0*vlpm*vlsm - 
     &    24.D0*vlpm*vlwm - 6.D0*vlpm**2 - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff9 = gcoeff9 + INVG**2*xnsq**2 * (  - 3.D0*b*vlsm*omro**(-1)
     &     + 3.D0/4.D0*vlpm**2*omro**(-1) + 3.D0*vdmp*omro**(-1) + 1.D0/
     &    4.D0*omro**(-1)*pi**2 )
      gcoeff9 = gcoeff9 + xnsq * (  - 12.D0*b*vlpm**2 + 96.D0*b*vltm*
     &    vlwm - 144.D0*b*vltm + 48.D0*b*vltm**2 + 48.D0*b*vlwm + 48.D0
     &    *b*vdw + 96.D0*b*vdt - 12.D0*b*pi**2 + 48.D0*b - 24.D0*vlpm )
      gcoeff9 = gcoeff9 + xnsq**2 * (  - 48.D0*b*vlsm*vltm + 18.D0*b*
     &    vlsm*omro**(-2) - 12.D0*b*vlsm*omro**(-1) - 12.D0*b*vlsm + 12.
     &    D0*b*vlsm**2 + 96.D0*b*vltm - 48.D0*b*vlwm - 48.D0*b*vdt - 12.
     &    D0*b*omro**(-1) + 4.D0*b*pi**2 - 12.D0*b - 9.D0/2.D0*vlpm**2*
     &    omro**(-2) + 9.D0/2.D0*vlpm**2*omro**(-1) - 21.D0/2.D0*
     &    vlpm**2 - 18.D0*vdmp*omro**(-2) + 18.D0*vdmp*omro**(-1) - 42.D
     &    0*vdmp - 3.D0/2.D0*omro**(-2)*pi**2 + 3.D0/2.D0*omro**(-1)*
     &    pi**2 - 7.D0/2.D0*pi**2 )
      gcoeff9 = gcoeff9 - 24.D0*b*vlpm**2 + 48.D0*b*vltm - 48.D0*b*
     &    vltm**2 + 48.D0*b*vlwm - 48.D0*b*vlwm**2 - 48.D0*b*vdw - 48.D0
     &    *b*vdt + 8.D0*b*pi**2 - 24.D0*b + 24.D0*vlpm*vlsm - 48.D0*
     &    vlpm*vltm - 24.D0*vlpm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*
     &    pi**2

      gcoeff10 =  + t1**(-2)*ro*xnsq * ( 6.D0*b*vltm**2 + 6.D0*b*vdt + 
     &    b*pi**2 )
      gcoeff10 = gcoeff10 + t1**(-2)*ro * (  - 6.D0*b*vltm**2 - 6.D0*b*
     &    vdt - b*pi**2 )
      gcoeff10 = gcoeff10 + t1**(-1)*ro*xnsq * (  - 6.D0*b*vltm**2 - 6.D
     &    0*b*vdt - b*pi**2 )
      gcoeff10 = gcoeff10 + t1**(-1)*ro * ( 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vltm - 12.D0*vdmp - 12.D0*vdmb + 3.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1**(-1)*xnsq * ( 16.D0*XLF*b*TR*rmuom2*xn
     &     - 48.D0*b*vltm*vlwm + 120.D0*b*vltm + 12.D0*b*pi**2 + 44.D0*
     &    b*rmuom2 - 24.D0*b )
      gcoeff10 = gcoeff10 + t1**(-1)*xnsq**2 * ( 24.D0*b*vlsm*vltm - 72.
     &    D0*b*vltm - 6.D0*b*pi**2 - 44.D0*b*rmuom2 )
      gcoeff10 = gcoeff10 + t1**(-1) * (  - 16.D0*XLF*b*TR*rmuom2*xn - 
     &    48.D0*b*vltm + 24.D0*b - 12.D0*vlpm*vlsm + 24.D0*vlpm*vltm + 
     &    24.D0*vdmp + 24.D0*vdmb - 6.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1*UBAR**(-2)*xnsq**2 * (  - 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + t1*UBAR**(-2) * ( 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + t1*UBAR**(-1)*xnsq * (  - 36.D0*b*vlwm - 6.D
     &    0*b )
      gcoeff10 = gcoeff10 + t1*UBAR**(-1)*xnsq**2 * ( 60.D0*b*vlwm + 6.D
     &    0*b )
      gcoeff10 = gcoeff10 + t1*INVG*xnsq * (  - 24.D0*b*vltm*vlwm - 24.D
     &    0*b*vltm**2 - 12.D0*b*vdw - 36.D0*b*vdt + 4.D0*b*pi**2 + 12.D0
     &    *vlpm*vlsm - 12.D0*vlpm*vltm - 12.D0*vlpm*vlwm - 6.D0*vlpm**2
     &     - 24.D0*vdmp - 2.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1*INVG*xnsq**2 * (  - 12.D0*b*vlsm*vltm - 
     &    36.D0*b*vlsm*vlwm + 12.D0*b*vlsm**2 - 36.D0*b*vdw - 12.D0*b*
     &    vdt + 4.D0*b*pi**2 + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2
     &     )
      gcoeff10 = gcoeff10 + t1*INVG * ( 12.D0*b*vlpm**2 + 36.D0*b*
     &    vltm**2 + 12.D0*b*vlwm**2 + 12.D0*b*vdw + 36.D0*b*vdt - 4.D0*
     &    b*pi**2 - 24.D0*vlpm*vlsm + 60.D0*vlpm*vltm - 12.D0*vlpm*vlwm
     &     + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1*xnsq * (  - 24.D0*b*vltm**2 + 24.D0*b*
     &    vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt + 48.D0*vlpm*vlsm - 48.D0
     &    *vlpm*vltm - 48.D0*vlpm*vlwm - 24.D0*vlpm**2 - 96.D0*vdmp - 8.
     &    D0*pi**2 )
      gcoeff10 = gcoeff10 + t1 * ( 48.D0*vlpm*vltm - 48.D0*vlpm*vlwm )
      gcoeff10 = gcoeff10 + t1**2*UBAR**(-2)*xnsq * (  - 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + t1**2*UBAR**(-2)*xnsq**2 * ( 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + t1**2*INVG*xnsq * ( 24.D0*b*vltm*vlwm + 48.D
     &    0*b*vltm**2 - 24.D0*b*vlwm**2 - 12.D0*b*vdw + 60.D0*b*vdt - 4.
     &    D0*b*pi**2 - 60.D0*vlpm*vlsm + 60.D0*vlpm*vltm + 60.D0*vlpm*
     &    vlwm + 30.D0*vlpm**2 + 120.D0*vdmp + 10.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1**2*INVG*xnsq**2 * ( 24.D0*b*vlsm*vltm + 
     &    24.D0*b*vlsm*vlwm - 12.D0*b*vlsm**2 + 24.D0*b*vdw + 24.D0*b*
     &    vdt - 4.D0*b*pi**2 - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2
     &     )
      gcoeff10 = gcoeff10 + t1**2*INVG * (  - 12.D0*b*vlpm**2 - 24.D0*b
     &    *vltm**2 - 24.D0*b*vlwm**2 - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0
     &    *b*pi**2 + 24.D0*vlpm*vlsm - 96.D0*vlpm*vltm + 48.D0*vlpm*
     &    vlwm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1**2*xnsq * (  - 48.D0*vlpm*vlsm + 48.D0*
     &    vlpm*vltm + 48.D0*vlpm*vlwm + 24.D0*vlpm**2 + 96.D0*vdmp + 8.D
     &    0*pi**2 )
      gcoeff10 = gcoeff10 + t1**3*INVG*xnsq * (  - 24.D0*b*vltm**2 + 24.
     &    D0*b*vlwm**2 + 24.D0*b*vdw - 24.D0*b*vdt + 96.D0*vlpm*vlsm - 
     &    96.D0*vlpm*vltm - 96.D0*vlpm*vlwm - 48.D0*vlpm**2 - 192.D0*
     &    vdmp - 16.D0*pi**2 )
      gcoeff10 = gcoeff10 + t1**3*INVG * ( 48.D0*vlpm*vltm - 48.D0*vlpm
     &    *vlwm )
      gcoeff10 = gcoeff10 + t1**4*INVG*xnsq * (  - 48.D0*vlpm*vlsm + 48.
     &    D0*vlpm*vltm + 48.D0*vlpm*vlwm + 24.D0*vlpm**2 + 96.D0*vdmp
     &     + 8.D0*pi**2 )
      gcoeff10 = gcoeff10 + t2**(-2)*ro*xnsq * ( 6.D0*b*vlwm**2 + 6.D0*
     &    b*vdw + b*pi**2 )
      gcoeff10 = gcoeff10 + t2**(-2)*ro * (  - 6.D0*b*vlwm**2 - 6.D0*b*
     &    vdw - b*pi**2 )
      gcoeff10 = gcoeff10 + t2**(-1)*ro*xnsq * (  - 6.D0*b*vlwm**2 - 6.D
     &    0*b*vdw - b*pi**2 )
      gcoeff10 = gcoeff10 + t2**(-1)*ro * ( 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vlwm - 12.D0*vdmp - 12.D0*vdmb + 3.D0*pi**2 )
      gcoeff10 = gcoeff10 + t2**(-1)*xnsq * ( 16.D0*XLF*b*TR*rmuom2*xn
     &     - 48.D0*b*vltm*vlwm + 120.D0*b*vlwm + 12.D0*b*pi**2 + 44.D0*
     &    b*rmuom2 - 24.D0*b )
      gcoeff10 = gcoeff10 + t2**(-1)*xnsq**2 * ( 24.D0*b*vlsm*vlwm - 72.
     &    D0*b*vlwm - 6.D0*b*pi**2 - 44.D0*b*rmuom2 )
      gcoeff10 = gcoeff10 + t2**(-1) * (  - 16.D0*XLF*b*TR*rmuom2*xn - 
     &    48.D0*b*vlwm + 24.D0*b - 12.D0*vlpm*vlsm + 24.D0*vlpm*vlwm + 
     &    24.D0*vdmp + 24.D0*vdmb - 6.D0*pi**2 )
      gcoeff10 = gcoeff10 + ro*TBAR**(-2)*xnsq * ( 3.D0*b*vltm )
      gcoeff10 = gcoeff10 + ro*TBAR**(-2)*xnsq**2 * (  - 3.D0/2.D0*b*
     &    vltm )
      gcoeff10 = gcoeff10 + ro*TBAR**(-2) * (  - 3.D0/2.D0*b*vltm )
      gcoeff10 = gcoeff10 + ro*TBAR**(-1)*xnsq * ( 6.D0*b*vltm + 3.D0/2.
     &    D0*b )
      gcoeff10 = gcoeff10 + ro*TBAR**(-1)*xnsq**2 * (  - 12.D0*b*vltm
     &     - 3.D0/2.D0*b )
      gcoeff10 = gcoeff10 + ro*xnsq * ( 12.D0*vlpm*vlsm - 12.D0*vlpm*
     &    vltm - 12.D0*vlpm*vlwm - 24.D0*vdmp - 24.D0*vdmb + 6.D0*pi**2
     &     )
      gcoeff10 = gcoeff10 + ro**2*TBAR**(-2)*xnsq * (  - 3.D0/8.D0*b*
     &    vltm )
      gcoeff10 = gcoeff10 + ro**2*TBAR**(-2)*xnsq**2 * ( 3.D0/8.D0*b*
     &    vltm )
      gcoeff10 = gcoeff10 + TBAR**(-1)*xnsq * (  - 84.D0*b*vltm - 12.D0
     &    *b )
      gcoeff10 = gcoeff10 + TBAR**(-1)*xnsq**2 * ( 54.D0*b*vltm + 6.D0*
     &    b )
      gcoeff10 = gcoeff10 + TBAR**(-1) * ( 30.D0*b*vltm + 6.D0*b )
      gcoeff10 = gcoeff10 + UBAR**(-2)*xnsq * ( 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + UBAR**(-2) * (  - 6.D0*b*vlwm )
      gcoeff10 = gcoeff10 + UBAR**(-1)*xnsq * (  - 60.D0*b*vlwm - 6.D0*
     &    b )
      gcoeff10 = gcoeff10 + UBAR**(-1) * ( 36.D0*b*vlwm + 6.D0*b )
      gcoeff10 = gcoeff10 + INVG*xnsq * ( 24.D0*b*vltm*vlwm + 12.D0*b*
     &    vltm**2 + 12.D0*b*vlwm**2 + 24.D0*b*vdw + 24.D0*b*vdt - 4.D0*
     &    b*pi**2 )
      gcoeff10 = gcoeff10 + INVG*xnsq**2 * (  - 12.D0*b*vlsm*vltm + 3.D0
     &    *b*vlsm**2 - 12.D0*b*vdt + b*pi**2 - 3.D0*vlpm**2 - 12.D0*
     &    vdmp - pi**2 )
      gcoeff10 = gcoeff10 + INVG * (  - 3.D0*b*vlpm**2 - 12.D0*b*
     &    vltm**2 - 12.D0*b*vdt + b*pi**2 + 6.D0*vlpm*vlsm - 12.D0*vlpm
     &    *vltm - 3.D0*vlpm**2 - 12.D0*vdmp - pi**2 )
      gcoeff10 = gcoeff10 + xnsq * (  - 32.D0*XLF*b*TR*rmuom2*xn - 6.D0
     &    *b*vlpm**2 + 96.D0*b*vltm*vlwm - 18.D0*b*vltm - 24.D0*b*
     &    vltm**2 - 48.D0*b*vlwm - 48.D0*b*vlwm**2 + 24.D0*b*vdt - 38.D0
     &    *b*pi**2 + 54.D0*b - 24.D0*vlpm*vlsm + 24.D0*vlpm*vltm + 24.D0
     &    *vlpm*vlwm + 48.D0*vdmp + 48.D0*vdmb - 12.D0*pi**2 )
      gcoeff10 = gcoeff10 + xnsq**2 * (  - 24.D0*b*vlsm*vltm - 24.D0*b*
     &    vlsm*vlwm + 18.D0*b*vltm + 72.D0*b*vlwm + 12.D0*b*pi**2 + 88.D
     &    0*b*rmuom2 - 6.D0*b + 12.D0*vlpm**2 + 48.D0*vdmp + 4.D0*pi**2
     &     )
      gcoeff10 = gcoeff10 - 12.D0*b*vlpm**2 - 24.D0*b*vltm**2 - 24.D0*b
     &    *vlwm**2 - 24.D0*b*vdw - 24.D0*b*vdt + 4.D0*b*pi**2 + 24.D0*
     &    vlpm*vlsm - 48.D0*vlpm*vltm - 12.D0*vlpm**2 - 48.D0*vdmp - 4.D
     &    0*pi**2

      return
      end
