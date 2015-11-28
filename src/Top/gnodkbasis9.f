      double precision function gnodkbasis9()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis9 = 32.D0 - 2.D0*t1**(-2)*ro + 4.D0*t1**(-2)*ro**2 - 18.
     & D0*t1**(-1)*ro + 4.D0*t1**(-1)*ro**2 - 32.D0*t1 - 2.D0*t2**(-1)*
     & ro + 4.D0*t2**(-1)*ro**2

      return
      end
