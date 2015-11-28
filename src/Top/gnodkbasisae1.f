      double precision function gnodkbasisae1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisae1=  + xnsq * ( 24.D0 - 8.D0*t1**(-1) - 16.D0*t1 + 16.D
     &    0*t1**2 - 8.D0*t2**(-1) )
      gnodkbasisae1 = gnodkbasisae1 - 8.D0 + 8.D0*t1**(-1) + 8.D0*
     &    t2**(-1)

      return
      end
