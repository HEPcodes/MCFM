      double precision function gnodkbasis10()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis10 =  - 4.D0*t1**(-1)*ro + 8.D0*t1**(-1)*ro**2 + 64.D0*
     & t1 + 8.D0*t2**(-2)*ro**2 - 36.D0*t2**(-1)*ro + 8.D0*t2**(-1)*
     & ro**2 - 16.D0*t2**(-1)

      return
      end
