      double precision function qnodkbasis4(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qnodkbasis4=  + ro * (  - 8.D0*t1 + 8.D0*t1**2 )
      qnodkbasis4 = qnodkbasis4 + ro**2 * ( 2.D0 )

      return
      end
