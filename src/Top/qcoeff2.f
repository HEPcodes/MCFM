      double precision function qcoeff2()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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

      return
      end
