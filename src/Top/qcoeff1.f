      double precision function qcoeff1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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

      return
      end
