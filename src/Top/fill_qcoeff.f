      subroutine fill_qcoeff(
     . qcoeffe1,qcoeffe2,qcoeff0,qcoeff1,qcoeff2,
     . qcoeff3,qcoeff4,qcoeff5,qcoeff6)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      double precision
     . qcoeffe1,qcoeffe2,qcoeff0,qcoeff1,qcoeff2,
     . qcoeff3,qcoeff4,qcoeff5,qcoeff6

      qcoeffe1 =  + b**(-1)*vlpm*xn**(-1) * (  - 16.D0 + 8.D0*ro )
      qcoeffe1 = qcoeffe1 + vlsm*xn**(-1) * (  - 16.D0 )
      qcoeffe1 = qcoeffe1 + vltm*xn**(-1) * (  - 64.D0 )
      qcoeffe1 = qcoeffe1 + vltm*xn * ( 32.D0 )
      qcoeffe1 = qcoeffe1 + vlwm*xn**(-1) * ( 64.D0 )
      qcoeffe1 = qcoeffe1 + xn**(-1) * ( 40.D0 )
      qcoeffe1 = qcoeffe1 + xn * (  - 40.D0 )

      qcoeffe2 =  + xn**(-1) * ( 16.D0 )
      qcoeffe2 = qcoeffe2 + xn * (  - 16.D0 )

      qcoeff0 =  + b**(-2)*vlsm*xn**(-1) * (  - 8.D0*t1**(-1)*INVG - 32.
     &    D0*t1**(-1) - 32.D0*t1*INVG + 32.D0*INVG )
      qcoeff0 = qcoeff0 + b**(-1)*xn**(-1) * (  - 8.D0*pisq )
      qcoeff0 = qcoeff0 + b*vlpm*TR * ( 32.D0 )
      qcoeff0 = qcoeff0 + b*vlpm*xn**(-1) * (  - 24.D0 )
      qcoeff0 = qcoeff0 + b*xn**(-1)*bracks * ( 8.D0*t1**(-1)*INVG + 5.D
     &    0/2.D0*t1**(-1)*INVG**2 - 8.D0*t1**(-1) + 52.D0*t1*INVG**2 - 
     &    32.D0*t1**2*INVG**2 + 8.D0*t1**3*INVG**2 - 24.D0*INVG**2 )
      qcoeff0 = qcoeff0 + b*xn**(-1) * (  - 8.D0*pisq )
      qcoeff0 = qcoeff0 + b*xn*bracks * (  - 4.D0*t1**(-1)*INVG - 5.D0/
     &    4.D0*t1**(-1)*INVG**2 + 4.D0*t1**(-1) + 4.D0*t1*INVG - 9.D0*
     &    t1*INVG**2 + 4.D0*t1**2*INVG**2 - 8.D0*INVG + 6.D0*INVG**2 )
      qcoeff0 = qcoeff0 + b**3*vlpm*TR * (  - 32.D0/3.D0 )
      qcoeff0 = qcoeff0 + b**3*xn**(-1)*bracks * ( 16.D0*t1**(-1)*INVG
     &     + t1**(-1)*INVG**2 - 4.D0*t1*INVG**2 - 8.D0*INVG**2 )
      qcoeff0 = qcoeff0 + b**3*xn*bracks * (  - 4.D0*t1**(-1)*INVG + 1.D
     &    0/2.D0*t1**(-1)*INVG**2 + t1*INVG**2 + 2.D0*INVG**2 )
      qcoeff0 = qcoeff0 + b**5*xn**(-1)*bracks * (  - 7.D0/2.D0*
     &    t1**(-1)*INVG**2 )
      qcoeff0 = qcoeff0 + b**5*xn*bracks * ( 3.D0/4.D0*t1**(-1)*INVG**2
     &     )
      qcoeff0 = qcoeff0 + vlsm*TR * ( 64.D0/3.D0*XLF )
      qcoeff0 = qcoeff0 + vlsm*xn**(-1) * (  - 24.D0 + 16.D0*t1**(-1)*
     &    ro*INVG + 10.D0*t1**(-1)*ro*INVG**2 + 32.D0*t1**(-1)*ro - 8.D0
     &    *t1**(-1)*ro**2*INVG + 11.D0/2.D0*t1**(-1)*ro**2*INVG**2 - 4.D
     &    0*t1**(-1)*ro**3*INVG**2 + 48.D0*t1**(-1)*INVG - 24.D0*
     &    t1**(-1) - 28.D0*t1*ro*INVG**2 - 64.D0*t1*INVG + 80.D0*t1*
     &    INVG**2 - 8.D0*t1**3*INVG**2 - 32.D0*ro*INVG - 40.D0*ro*
     &    INVG**2 + 24.D0*ro**2*INVG**2 + 32.D0*INVG - 32.D0*INVG**2 )
      qcoeff0 = qcoeff0 + vlsm*xn * (  - 32.D0/3.D0 - 12.D0*t1**(-1)*ro
     &    *INVG + 2.D0*t1**(-1)*ro*INVG**2 - 8.D0*t1**(-1)*ro + 2.D0*
     &    t1**(-1)*ro**2*INVG - 15.D0/4.D0*t1**(-1)*ro**2*INVG**2 + 
     &    t1**(-1)*ro**3*INVG**2 + 8.D0*t1**(-1)*INVG + 12.D0*t1**(-1)
     &     + 5.D0*t1*ro*INVG**2 + 12.D0*t1*INVG - 8.D0*t1*INVG**2 - 4.D0
     &    *t1**2*INVG**2 + 18.D0*ro*INVG**2 - 6.D0*ro**2*INVG**2 + 8.D0
     &    *INVG - 8.D0*INVG**2 )
      qcoeff0 = qcoeff0 + vlsm**2*xn**(-1) * ( 8.D0 )
      qcoeff0 = qcoeff0 + vltm*TBAR**(-1)*xn**(-1) * (  - 32.D0 - 8.D0*
     &    t1**(-1)*ro*INVG + 16.D0*t1**(-1)*ro + 4.D0*t1**(-1)*ro**2*
     &    INVG - 32.D0*t1**(-1) + 8.D0*t1*ro*INVG - 64.D0*t1 - 24.D0*ro
     &    *INVG + 16.D0*ro + 32.D0*INVG )
      qcoeff0 = qcoeff0 + vltm*TBAR**(-1)*xn * ( 16.D0 + 4.D0*t1**(-1)*
     &    ro*INVG - 8.D0*t1**(-1)*ro - 2.D0*t1**(-1)*ro**2*INVG + 16.D0
     &    *t1**(-1) - 4.D0*t1*ro*INVG + 32.D0*t1 + 12.D0*ro*INVG - 8.D0
     &    *ro - 16.D0*INVG )
      qcoeff0 = qcoeff0 + vlwm*UBAR**(-1)*xn**(-1) * ( 32.D0 - 16.D0*
     &    t1**(-1)*ro - 4.D0*t1**(-1)*ro**2*INVG - 16.D0*t1*ro*INVG - 
     &    64.D0*t1 + 8.D0*t1**2*ro*INVG + 16.D0*ro*INVG - 8.D0*ro + 2.D0
     &    *ro**2*INVG )
      qcoeff0 = qcoeff0 + vdw*xn**(-1) * ( 64.D0 )
      qcoeff0 = qcoeff0 + vdt*xn**(-1) * (  - 64.D0 )
      qcoeff0 = qcoeff0 + vdt*xn * ( 32.D0 )
      qcoeff0 = qcoeff0 + TR*rmuom2 * (  - 64.D0/3.D0*XLF )
      qcoeff0 = qcoeff0 + TR * (  - 320.D0/9.D0 - 320.D0/9.D0*XLF - 64.D
     &    0/3.D0*ro )
      qcoeff0 = qcoeff0 + rmuom2*xn * ( 176.D0/3.D0 )
      qcoeff0 = qcoeff0 + xn**(-1)*brackt * ( 64.D0 - 40.D0*t1**(-1)*ro
     &    *INVG + 16.D0*t1**(-1)*ro*INVG**2 + 16.D0*t1**(-1)*ro + 6.D0*
     &    t1**(-1)*ro**2*INVG - 10.D0*t1**(-1)*ro**2*INVG**2 + 1.D0/2.D0
     &    *t1**(-1)*ro**3*INVG**2 + 64.D0*t1**(-1)*INVG + 16.D0*t1*ro*
     &    INVG - 40.D0*t1*ro*INVG**2 + 2.D0*t1*ro**2*INVG**2 - 32.D0*t1
     &    *INVG + 32.D0*t1*INVG**2 + 8.D0*t1**2*ro*INVG**2 + 48.D0*ro*
     &    INVG**2 + 32.D0*INVG - 64.D0*INVG**2 )
      qcoeff0 = qcoeff0 + xn**(-1)*brackw * (  - 64.D0 + 64.D0*t1**(-1)
     &    *ro*INVG - 16.D0*t1**(-1)*ro*INVG**2 - 8.D0*t1**(-1)*ro - 2.D0
     &    *t1**(-1)*ro**2*INVG + 14.D0*t1**(-1)*ro**2*INVG**2 - 64.D0*
     &    t1**(-1)*INVG + 32.D0*t1**(-1) + 80.D0*t1*ro*INVG**2 + 2.D0*
     &    t1*ro**2*INVG**2 + 96.D0*t1*INVG - 96.D0*t1*INVG**2 - 40.D0*
     &    t1**2*ro*INVG**2 + 32.D0*t1**2*INVG**2 + 8.D0*t1**3*ro*
     &    INVG**2 - 32.D0*ro*INVG - 48.D0*ro*INVG**2 - 8.D0*ro**2*
     &    INVG**2 - 64.D0*INVG + 64.D0*INVG**2 )
      qcoeff0 = qcoeff0 + xn**(-1) * ( 96.D0 - 8.D0*pisq )
      qcoeff0 = qcoeff0 + xn*brackt * (  - 32.D0 + 20.D0*t1**(-1)*ro*
     &    INVG - 8.D0*t1**(-1)*ro*INVG**2 - 8.D0*t1**(-1)*ro - 3.D0*
     &    t1**(-1)*ro**2*INVG + 5.D0*t1**(-1)*ro**2*INVG**2 - 1.D0/4.D0
     &    *t1**(-1)*ro**3*INVG**2 - 32.D0*t1**(-1)*INVG - 8.D0*t1*ro*
     &    INVG + 20.D0*t1*ro*INVG**2 - t1*ro**2*INVG**2 + 16.D0*t1*INVG
     &     - 16.D0*t1*INVG**2 - 4.D0*t1**2*ro*INVG**2 - 24.D0*ro*
     &    INVG**2 - 16.D0*INVG + 32.D0*INVG**2 )
      qcoeff0 = qcoeff0 + xn * ( 496.D0/9.D0 - 8.D0/3.D0*pisq )
      qcoeff0 = qcoeff0 + f1*xn**(-1) * (  - 16.D0 + 8.D0*ro )
      qcoeff0 = qcoeff0 + f2*xn**(-1) * (  - 40.D0*t1**(-1)*ro*INVG + 
     &    10.D0*t1**(-1)*ro*INVG**2 + 24.D0*t1**(-1)*ro - 32.D0*
     &    t1**(-1)*ro**2*INVG - 9.D0/2.D0*t1**(-1)*ro**2*INVG**2 - 32.D0
     &    *t1**(-1)*ro**2 + 8.D0*t1**(-1)*ro**3*INVG - 19.D0/2.D0*
     &    t1**(-1)*ro**3*INVG**2 + 4.D0*t1**(-1)*ro**4*INVG**2 + 40.D0*
     &    t1**(-1)*INVG - 88.D0*t1**(-1) + 32.D0*t1*ro*INVG - 108.D0*t1
     &    *ro*INVG**2 + 28.D0*t1*ro**2*INVG**2 - 128.D0*t1*INVG + 80.D0
     &    *t1*INVG**2 + 8.D0*t1**3*ro*INVG**2 - 8.D0*t1**3*INVG**2 - 32.
     &    D0*ro*INVG - 8.D0*ro*INVG**2 + 32.D0*ro**2*INVG + 64.D0*ro**2
     &    *INVG**2 - 24.D0*ro**3*INVG**2 + 96.D0*INVG - 32.D0*INVG**2 )
      qcoeff0 = qcoeff0 + f2*xn * ( 32.D0 - 16.D0*t1**(-1)*ro*INVG + 2.D
     &    0*t1**(-1)*ro*INVG**2 - 24.D0*t1**(-1)*ro + 13.D0*t1**(-1)*
     &    ro**2*INVG - 23.D0/4.D0*t1**(-1)*ro**2*INVG**2 + 8.D0*
     &    t1**(-1)*ro**2 - 2.D0*t1**(-1)*ro**3*INVG + 19.D0/4.D0*
     &    t1**(-1)*ro**3*INVG**2 - t1**(-1)*ro**4*INVG**2 + 8.D0*
     &    t1**(-1)*INVG + 28.D0*t1**(-1) - 8.D0*t1*ro*INVG + 13.D0*t1*
     &    ro*INVG**2 - 5.D0*t1*ro**2*INVG**2 + 12.D0*t1*INVG - 8.D0*t1*
     &    INVG**2 + 4.D0*t1**2*ro*INVG**2 - 4.D0*t1**2*INVG**2 + 26.D0*
     &    ro*INVG**2 - 8.D0*ro - 24.D0*ro**2*INVG**2 + 6.D0*ro**3*
     &    INVG**2 - 8.D0*INVG - 8.D0*INVG**2 )

      qcoeff1 =  + vlsm*xn**(-1) * ( 64.D0*t1**(-1)*ro*INVG + 16.D0*
     &    t1**(-1)*ro*INVG**2 + 16.D0*t1**(-1)*ro**2*INVG**2 + 64.D0*
     &    t1**(-1)*INVG + 64.D0*t1*INVG**2 + 64.D0*t1**2*INVG**2 + 64.D0
     &    *ro*INVG - 144.D0*ro*INVG**2 + 32.D0*ro**2*INVG**2 - 192.D0*
     &    INVG )
      qcoeff1 = qcoeff1 + vlsm*xn * (  - 64.D0*t1**(-1)*ro*INVG + 16.D0
     &    *t1**(-1)*ro*INVG**2 - 12.D0*t1**(-1)*ro**2*INVG**2 + 64.D0*
     &    t1**(-1)*INVG - 64.D0*t1**(-1) - 16.D0*ro*INVG + 64.D0*ro*
     &    INVG**2 - 8.D0*ro**2*INVG**2 + 128.D0*INVG - 64.D0*INVG**2 )
      qcoeff1 = qcoeff1 + vltm*TBAR**(-1)*xn**(-1) * (  - 128.D0 - 32.D0
     &    *t1**(-1)*ro*INVG + 32.D0*t1**(-1)*ro + 8.D0*t1**(-1)*ro**2*
     &    INVG - 128.D0*t1**(-1) - 64.D0*ro*INVG + 128.D0*INVG )
      qcoeff1 = qcoeff1 + vltm*TBAR**(-1)*xn * ( 64.D0 + 16.D0*t1**(-1)
     &    *ro*INVG - 16.D0*t1**(-1)*ro - 4.D0*t1**(-1)*ro**2*INVG + 64.D
     &    0*t1**(-1) + 32.D0*ro*INVG - 64.D0*INVG )
      qcoeff1 = qcoeff1 + vlwm*UBAR**(-1)*xn**(-1) * (  - 128.D0 - 32.D0
     &    *t1**(-1)*ro - 8.D0*t1**(-1)*ro**2*INVG - 64.D0*t1*INVG + 64.D
     &    0*INVG )
      qcoeff1 = qcoeff1 + xn**(-1)*brackt * (  - 32.D0*t1**(-1)*ro*INVG
     &     + 32.D0*t1**(-1)*ro*INVG**2 - 8.D0*t1**(-1)*ro**2*INVG**2 + 
     &    128.D0*t1**(-1)*INVG + 64.D0*ro*INVG**2 + 128.D0*INVG - 128.D0
     &    *INVG**2 )
      qcoeff1 = qcoeff1 + xn**(-1)*brackw * ( 32.D0*t1**(-1)*ro*INVG - 
     &    16.D0*t1**(-1)*ro*INVG**2 + 8.D0*t1**(-1)*ro**2*INVG**2 - 64.D
     &    0*t1**(-1)*INVG + 64.D0*t1*INVG**2 - 64.D0*t1**2*INVG**2 - 16.
     &    D0*ro*INVG**2 )
      qcoeff1 = qcoeff1 + xn*brackt * ( 16.D0*t1**(-1)*ro*INVG - 16.D0*
     &    t1**(-1)*ro*INVG**2 + 4.D0*t1**(-1)*ro**2*INVG**2 - 64.D0*
     &    t1**(-1)*INVG - 32.D0*ro*INVG**2 - 64.D0*INVG + 64.D0*INVG**2
     &     )
      qcoeff1 = qcoeff1 + f2*xn**(-1) * ( 16.D0*t1**(-1)*ro*INVG**2 - 
     &    64.D0*t1**(-1)*ro**2*INVG - 16.D0*t1**(-1)*ro**3*INVG**2 + 64.
     &    D0*t1**(-1)*INVG - 64.D0*t1*ro*INVG**2 + 64.D0*t1*INVG**2 - 
     &    64.D0*t1**2*ro*INVG**2 + 64.D0*t1**2*INVG**2 + 192.D0*ro*INVG
     &     - 144.D0*ro*INVG**2 - 64.D0*ro**2*INVG + 176.D0*ro**2*
     &    INVG**2 - 32.D0*ro**3*INVG**2 - 128.D0*INVG )
      qcoeff1 = qcoeff1 + f2*xn * (  - 112.D0*t1**(-1)*ro*INVG + 16.D0*
     &    t1**(-1)*ro*INVG**2 + 32.D0*t1**(-1)*ro + 56.D0*t1**(-1)*
     &    ro**2*INVG - 28.D0*t1**(-1)*ro**2*INVG**2 + 12.D0*t1**(-1)*
     &    ro**3*INVG**2 + 64.D0*t1**(-1)*INVG - 96.D0*ro*INVG + 128.D0*
     &    ro*INVG**2 + 16.D0*ro**2*INVG - 72.D0*ro**2*INVG**2 + 8.D0*
     &    ro**3*INVG**2 + 64.D0*INVG - 64.D0*INVG**2 )

      qcoeff2 =  + vlsm*xn**(-1) * ( 32.D0*t1**(-1)*ro*INVG - 20.D0*
     &    t1**(-1)*ro*INVG**2 + 8.D0*t1**(-1)*ro**2*INVG**2 - 80.D0*
     &    t1**(-1)*INVG - 32.D0*t1*INVG**2 - 16.D0*t1**2*INVG**2 - 20.D0
     &    *ro*INVG**2 - 16.D0*INVG + 64.D0*INVG**2 )
      qcoeff2 = qcoeff2 + vlsm*xn * (  - 8.D0*t1**(-1)*ro*INVG + 4.D0*
     &    t1**(-1)*ro*INVG**2 - 2.D0*t1**(-1)*ro**2*INVG**2 + 16.D0*
     &    t1**(-1)*INVG + 8.D0*t1*INVG**2 + 8.D0*ro*INVG**2 + 16.D0*
     &    INVG - 16.D0*INVG**2 )
      qcoeff2 = qcoeff2 + vltm*TBAR**(-1)*xn**(-1) * (  - 4.D0*t1**(-1)
     &    *ro*INVG - 16.D0*t1**(-1) - 8.D0*ro*INVG + 16.D0*INVG )
      qcoeff2 = qcoeff2 + vltm*TBAR**(-1)*xn * ( 2.D0*t1**(-1)*ro*INVG
     &     + 8.D0*t1**(-1) + 4.D0*ro*INVG - 8.D0*INVG )
      qcoeff2 = qcoeff2 + vlwm*UBAR**(-1)*xn**(-1) * ( 16.D0*t1*INVG + 
     &    8.D0*ro*INVG - 16.D0*INVG )
      qcoeff2 = qcoeff2 + xn**(-1)*brackt * (  - 16.D0*t1*INVG**2 + 8.D0
     &    *ro*INVG**2 + 16.D0*INVG )
      qcoeff2 = qcoeff2 + xn**(-1)*brackw * ( 4.D0*t1**(-1)*ro*INVG**2
     &     + 16.D0*t1**(-1)*INVG - 16.D0*t1*INVG**2 + 16.D0*t1**2*
     &    INVG**2 - 4.D0*ro*INVG**2 )
      qcoeff2 = qcoeff2 + xn*brackt * ( 8.D0*t1*INVG**2 - 4.D0*ro*
     &    INVG**2 - 8.D0*INVG )
      qcoeff2 = qcoeff2 + f2*xn**(-1) * ( 112.D0*t1**(-1)*ro*INVG - 20.D
     &    0*t1**(-1)*ro*INVG**2 - 32.D0*t1**(-1)*ro**2*INVG + 28.D0*
     &    t1**(-1)*ro**2*INVG**2 - 8.D0*t1**(-1)*ro**3*INVG**2 - 80.D0*
     &    t1**(-1)*INVG + 32.D0*t1*ro*INVG**2 - 32.D0*t1*INVG**2 + 16.D0
     &    *t1**2*ro*INVG**2 - 16.D0*t1**2*INVG**2 + 16.D0*ro*INVG - 84.D
     &    0*ro*INVG**2 + 20.D0*ro**2*INVG**2 - 16.D0*INVG + 64.D0*
     &    INVG**2 )
      qcoeff2 = qcoeff2 + f2*xn * (  - 24.D0*t1**(-1)*ro*INVG + 4.D0*
     &    t1**(-1)*ro*INVG**2 + 8.D0*t1**(-1)*ro**2*INVG - 6.D0*
     &    t1**(-1)*ro**2*INVG**2 + 2.D0*t1**(-1)*ro**3*INVG**2 + 16.D0*
     &    t1**(-1)*INVG - 8.D0*t1*ro*INVG**2 + 8.D0*t1*INVG**2 - 8.D0*
     &    ro*INVG + 24.D0*ro*INVG**2 - 8.D0*ro**2*INVG**2 + 8.D0*INVG
     &     - 16.D0*INVG**2 )

      qcoeff3 =  + vlsm*xn**(-1) * ( 8.D0*ro*INVG - 8.D0*INVG )
      qcoeff3 = qcoeff3 + vlsm*xn * (  - t1**(-1)*ro*INVG - 4.D0*
     &    t1**(-1) - 2.D0*ro*INVG + 4.D0*INVG )
      qcoeff3 = qcoeff3 + xn**(-1)*brackt * (  - 2.D0*t1**(-1)*ro*INVG
     &     - 8.D0*t1**(-1) + 8.D0*INVG )
      qcoeff3 = qcoeff3 + xn**(-1)*brackw * ( 2.D0*t1**(-1)*ro*INVG + 8.
     &    D0*t1**(-1) )
      qcoeff3 = qcoeff3 + xn*brackt * ( t1**(-1)*ro*INVG + 4.D0*
     &    t1**(-1) - 4.D0*INVG )
      qcoeff3 = qcoeff3 + f2*xn**(-1) * ( 16.D0*ro*INVG - 8.D0*ro**2*
     &    INVG - 8.D0*INVG )
      qcoeff3 = qcoeff3 + f2*xn * (  - t1**(-1)*ro*INVG + 4.D0*t1**(-1)
     &    *ro + t1**(-1)*ro**2*INVG - 4.D0*t1**(-1) - 6.D0*ro*INVG + 2.D
     &    0*ro**2*INVG + 4.D0*INVG )

      qcoeff4 =  + b**(-1)*vlpm*xn**(-1) * (  - 32.D0 )
      qcoeff4 = qcoeff4 + vlsm*xn**(-1) * (  - 96.D0*t1**(-1)*ro*INVG
     &     + 48.D0*t1**(-1)*ro*INVG**2 - 16.D0*t1**(-1)*ro**2*INVG**2
     &     + 192.D0*t1**(-1)*INVG - 128.D0*t1**(-1) + 64.D0*t1**2*
     &    INVG**2 + 48.D0*ro*INVG**2 + 128.D0*INVG - 128.D0*INVG**2 )
      qcoeff4 = qcoeff4 + vlsm*xn * ( 24.D0*t1**(-1)*ro*INVG - 8.D0*
     &    t1**(-1)*ro*INVG**2 + 4.D0*t1**(-1)*ro**2*INVG**2 - 32.D0*
     &    t1**(-1)*INVG + 32.D0*t1**(-1) - 24.D0*ro*INVG**2 - 64.D0*
     &    INVG + 32.D0*INVG**2 )
      qcoeff4 = qcoeff4 + vltm*TBAR**(-1)*xn**(-1) * ( 16.D0*t1**(-1)*
     &    ro*INVG + 64.D0*t1**(-1) + 16.D0*ro*INVG - 64.D0*INVG )
      qcoeff4 = qcoeff4 + vltm*TBAR**(-1)*xn * (  - 8.D0*t1**(-1)*ro*
     &    INVG - 32.D0*t1**(-1) - 8.D0*ro*INVG + 32.D0*INVG )
      qcoeff4 = qcoeff4 + vlwm*UBAR**(-1)*xn**(-1) * (  - 64.D0*t1*INVG
     &     - 16.D0*ro*INVG + 64.D0*INVG )
      qcoeff4 = qcoeff4 + xn**(-1)*brackt * (  - 16.D0*t1**(-1)*ro*
     &    INVG**2 - 64.D0*t1**(-1)*INVG - 16.D0*ro*INVG**2 - 64.D0*INVG
     &     + 64.D0*INVG**2 )
      qcoeff4 = qcoeff4 + xn**(-1)*brackw * ( 128.D0*t1*INVG**2 - 64.D0
     &    *t1**2*INVG**2 - 64.D0*INVG**2 )
      qcoeff4 = qcoeff4 + xn*brackt * ( 8.D0*t1**(-1)*ro*INVG**2 + 32.D0
     &    *t1**(-1)*INVG + 8.D0*ro*INVG**2 + 32.D0*INVG - 32.D0*INVG**2
     &     )
      qcoeff4 = qcoeff4 + f2*xn**(-1) * (  - 256.D0*t1**(-1)*ro*INVG + 
     &    48.D0*t1**(-1)*ro*INVG**2 + 128.D0*t1**(-1)*ro + 96.D0*
     &    t1**(-1)*ro**2*INVG - 64.D0*t1**(-1)*ro**2*INVG**2 + 16.D0*
     &    t1**(-1)*ro**3*INVG**2 + 192.D0*t1**(-1)*INVG - 64.D0*t1**2*
     &    ro*INVG**2 + 64.D0*t1**2*INVG**2 - 128.D0*ro*INVG + 176.D0*ro
     &    *INVG**2 - 48.D0*ro**2*INVG**2 + 64.D0*INVG - 128.D0*INVG**2
     &     )
      qcoeff4 = qcoeff4 + f2*xn * (  - 96.D0 + 48.D0*t1**(-1)*ro*INVG
     &     - 8.D0*t1**(-1)*ro*INVG**2 - 32.D0*t1**(-1)*ro - 24.D0*
     &    t1**(-1)*ro**2*INVG + 12.D0*t1**(-1)*ro**2*INVG**2 - 4.D0*
     &    t1**(-1)*ro**3*INVG**2 - 32.D0*t1**(-1)*INVG + 48.D0*ro*INVG
     &     - 56.D0*ro*INVG**2 + 24.D0*ro**2*INVG**2 - 32.D0*INVG + 32.D0
     &    *INVG**2 )
      qcoeff4 = qcoeff4 + f3*xn * ( 32.D0 )

      qcoeff5 =  + vlsm*xn**(-1) * (  - 8.D0*t1**(-1)*ro*INVG**2 - 32.D0
     &    *t1**(-1)*INVG - 32.D0*t1**2*INVG**2 + 24.D0*ro*INVG**2 + 32.D
     &    0*INVG )
      qcoeff5 = qcoeff5 + vlsm*xn * ( 16.D0*t1*INVG**2 - 8.D0*ro*
     &    INVG**2 - 16.D0*INVG )
      qcoeff5 = qcoeff5 + vltm*TBAR**(-1)*xn**(-1) * ( 8.D0*t1**(-1)*ro
     &    *INVG + 32.D0*t1**(-1) - 32.D0*INVG )
      qcoeff5 = qcoeff5 + vltm*TBAR**(-1)*xn * (  - 4.D0*t1**(-1)*ro*
     &    INVG - 16.D0*t1**(-1) + 16.D0*INVG )
      qcoeff5 = qcoeff5 + vlwm*UBAR**(-1)*xn**(-1) * ( 32.D0*t1*INVG - 
     &    32.D0*INVG )
      qcoeff5 = qcoeff5 + xn**(-1)*brackt * (  - 16.D0*t1**(-1)*ro*
     &    INVG**2 - 64.D0*t1**(-1)*INVG - 32.D0*t1*INVG**2 + 64.D0*
     &    INVG**2 )
      qcoeff5 = qcoeff5 + xn**(-1)*brackw * (  - 8.D0*t1**(-1)*ro*
     &    INVG**2 - 32.D0*t1**(-1)*INVG - 96.D0*t1*INVG**2 + 32.D0*
     &    t1**2*INVG**2 + 8.D0*ro*INVG**2 + 32.D0*INVG + 64.D0*INVG**2
     &     )
      qcoeff5 = qcoeff5 + xn*brackt * ( 8.D0*t1**(-1)*ro*INVG**2 + 32.D0
     &    *t1**(-1)*INVG + 16.D0*t1*INVG**2 - 32.D0*INVG**2 )
      qcoeff5 = qcoeff5 + f2*xn**(-1) * ( 32.D0*t1**(-1)*ro*INVG - 8.D0
     &    *t1**(-1)*ro*INVG**2 + 8.D0*t1**(-1)*ro**2*INVG**2 - 32.D0*
     &    t1**(-1)*INVG + 32.D0*t1**2*ro*INVG**2 - 32.D0*t1**2*INVG**2
     &     - 32.D0*ro*INVG + 24.D0*ro*INVG**2 - 24.D0*ro**2*INVG**2 + 
     &    32.D0*INVG )
      qcoeff5 = qcoeff5 + f2*xn * (  - 8.D0*t1**(-1)*ro*INVG - 32.D0*
     &    t1**(-1) - 16.D0*t1*ro*INVG**2 + 16.D0*t1*INVG**2 + 16.D0*ro*
     &    INVG - 8.D0*ro*INVG**2 + 8.D0*ro**2*INVG**2 )

      qcoeff6 =  + vlsm*xn**(-1) * ( 4.D0*t1**(-1)*ro*INVG**2 + 16.D0*
     &    t1**(-1)*INVG + 16.D0*t1**2*INVG**2 - 12.D0*ro*INVG**2 - 16.D0
     &    *INVG )
      qcoeff6 = qcoeff6 + vlsm*xn * (  - 8.D0*t1*INVG**2 + 4.D0*ro*
     &    INVG**2 + 8.D0*INVG )
      qcoeff6 = qcoeff6 + vltm*TBAR**(-1)*xn**(-1) * (  - 4.D0*t1**(-1)
     &    *ro*INVG - 16.D0*t1**(-1) + 16.D0*INVG )
      qcoeff6 = qcoeff6 + vltm*TBAR**(-1)*xn * ( 2.D0*t1**(-1)*ro*INVG
     &     + 8.D0*t1**(-1) - 8.D0*INVG )
      qcoeff6 = qcoeff6 + vlwm*UBAR**(-1)*xn**(-1) * (  - 16.D0*t1*INVG
     &     + 16.D0*INVG )
      qcoeff6 = qcoeff6 + xn**(-1)*brackt * ( 8.D0*t1**(-1)*ro*INVG**2
     &     + 32.D0*t1**(-1)*INVG + 16.D0*t1*INVG**2 - 32.D0*INVG**2 )
      qcoeff6 = qcoeff6 + xn**(-1)*brackw * ( 4.D0*t1**(-1)*ro*INVG**2
     &     + 16.D0*t1**(-1)*INVG + 48.D0*t1*INVG**2 - 16.D0*t1**2*
     &    INVG**2 - 4.D0*ro*INVG**2 - 16.D0*INVG - 32.D0*INVG**2 )
      qcoeff6 = qcoeff6 + xn*brackt * (  - 4.D0*t1**(-1)*ro*INVG**2 - 
     &    16.D0*t1**(-1)*INVG - 8.D0*t1*INVG**2 + 16.D0*INVG**2 )
      qcoeff6 = qcoeff6 + f2*xn**(-1) * (  - 16.D0*t1**(-1)*ro*INVG + 4.
     &    D0*t1**(-1)*ro*INVG**2 - 4.D0*t1**(-1)*ro**2*INVG**2 + 16.D0*
     &    t1**(-1)*INVG - 16.D0*t1**2*ro*INVG**2 + 16.D0*t1**2*INVG**2
     &     + 16.D0*ro*INVG - 12.D0*ro*INVG**2 + 12.D0*ro**2*INVG**2 - 
     &    16.D0*INVG )
      qcoeff6 = qcoeff6 + f2*xn * ( 4.D0*t1**(-1)*ro*INVG + 16.D0*
     &    t1**(-1) + 8.D0*t1*ro*INVG**2 - 8.D0*t1*INVG**2 - 8.D0*ro*
     &    INVG + 4.D0*ro*INVG**2 - 4.D0*ro**2*INVG**2 )

      return
      end
