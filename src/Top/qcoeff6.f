      double precision function qcoeff6()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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
