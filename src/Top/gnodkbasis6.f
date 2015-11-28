      double precision function gnodkbasis6()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis6 =  - 4.D0*t1**(-1)*ro**2 - 32.D0*t1 - 4.D0*t2**(-2)*
     & ro**2 + 16.D0*t2**(-1)*ro - 4.D0*t2**(-1)*ro**2 + 8.D0*t2**(-1)

      return
      end
