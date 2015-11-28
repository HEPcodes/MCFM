      double precision function qcoeff0()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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

      return
      end
