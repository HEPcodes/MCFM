      double precision function qcoeff5()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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

      return
      end
