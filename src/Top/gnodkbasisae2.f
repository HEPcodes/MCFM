      double precision function gnodkbasisae2()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisae2=  + xnsq * (  - 8.D0 + 4.D0*t1**(-1) + 4.D0*
     &    t2**(-1) )
      gnodkbasisae2 = gnodkbasisae2 - 4.D0*t1**(-1) - 4.D0*t2**(-1)

      return
      end
