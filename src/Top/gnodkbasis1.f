      double precision function gnodkbasis1()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis1 =  - 4.D0*t1**(-2)*ro + 8.D0*t1**(-1)

      return
      end
