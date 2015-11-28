      double precision function gnodkbasisbe2()
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      gnodkbasisbe2=  + 4.D0*t1**(-1) + 4.D0*t2**(-1)

      return
      end
