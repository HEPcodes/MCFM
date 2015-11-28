      double precision function qnodkbasis6(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qnodkbasis6=  + ro * ( 8.D0 - 16.D0*t1**2 )
      qnodkbasis6 = qnodkbasis6 + ro**2 * (  - 4.D0 )

      return
      end
