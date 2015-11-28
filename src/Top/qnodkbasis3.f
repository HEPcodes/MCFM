      double precision function qnodkbasis3(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qnodkbasis3=  + ro * ( 80.D0 )
      qnodkbasis3 = qnodkbasis3 + 64.D0 - 128.D0*t1 + 320.D0*t1**2

      return
      end
