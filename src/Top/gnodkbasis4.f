      double precision function gnodkbasis4()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis4 = 2.D0*t1**(-2)*ro - 2.D0*t1**(-2)*ro**2 + 4.D0*
     & t1**(-1)*ro - 2.D0*t1**(-1)*ro**2 + 16.D0*t1 - 4.D0*t2**(-1)*ro
     &  - 2.D0*t2**(-1)*ro**2

      return
      end
