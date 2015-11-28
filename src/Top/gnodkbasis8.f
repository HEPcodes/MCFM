      double precision function gnodkbasis8()
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'virtexp.f'

      gnodkbasis8 =  - 32.D0 + 2.D0*t1**(-2)*ro - 4.D0*t1**(-2)*ro**2
     &  + 26.D0*t1**(-1)*ro - 4.D0*t1**(-1)*ro**2 - 8.D0*t1**(-1) + 32.D
     & 0*t1 + 10.D0*t2**(-1)*ro - 4.D0*t2**(-1)*ro**2 - 8.D0*t2**(-1)

      return
      end
