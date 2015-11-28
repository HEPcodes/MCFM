      double precision function gnodkbasisbe0()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisbe0=  + ro * ( 4.D0*t1**(-1) + 4.D0*t2**(-1) )
      gnodkbasisbe0 = gnodkbasisbe0 + ro**2 * (  - t1**(-2) - 2.D0*
     &    t1**(-1) - t2**(-2) - 2.D0*t2**(-1) )
      gnodkbasisbe0 = gnodkbasisbe0 - 8.D0 + 4.D0*t1**(-1) + 4.D0*
     &    t2**(-1)

      return
      end
