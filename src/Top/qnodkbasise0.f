      double precision function qnodkbasise0(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qnodkbasise0=  + ro * ( 8.D0 )
      qnodkbasise0 = qnodkbasise0 + 16.D0 - 32.D0*t1 + 32.D0*t1**2

      return
      end
