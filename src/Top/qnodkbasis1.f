      double precision function qnodkbasis1(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qnodkbasis1=  + ro * ( 2.D0 - 4.D0*t1 )
      qnodkbasis1 = qnodkbasis1 - 8.D0*t1 + 24.D0*t1**2 - 16.D0*t1**3

      return
      end
