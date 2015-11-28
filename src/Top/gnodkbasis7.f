      double precision function gnodkbasis7()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis7 =  - 4.D0*t1**(-2)*ro - 4.D0*t1**(-1)*ro + 16.D0*
     & t1**(-1) - 4.D0*t2**(-1)*ro

      return
      end
