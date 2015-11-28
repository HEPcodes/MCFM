      double precision function qcoeff3()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

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

      return
      end
