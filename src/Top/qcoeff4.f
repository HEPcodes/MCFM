      double precision function qcoeff4()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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

      return
      end
