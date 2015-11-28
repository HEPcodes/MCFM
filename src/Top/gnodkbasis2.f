      double precision function gnodkbasis2()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis2 =  - 8.D0 - 2.D0*t1**(-2)*ro**2 + 8.D0*t1**(-1)*ro + 
     & 2.D0*t1**(-1)*ro**2 + 48.D0*t1 - 32.D0*t1**2 - 4.D0*t2**(-1)*ro
     &  + 2.D0*t2**(-1)*ro**2 - 16.D0*ro

      return
      end
