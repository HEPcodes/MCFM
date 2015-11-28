      double precision function qnodkbasis2(i1,i2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      qnodkbasis2=  + ro * (  - 8.D0*t1 )
      qnodkbasis2 = qnodkbasis2 - 32.D0*t1**3

      return
      end
