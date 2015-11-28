      subroutine fill_qnodkbasis(i1,i2,
     . qnodkbasise0,qnodkbasise1,qnodkbasis0,qnodkbasis1,
     . qnodkbasis2,qnodkbasis3,qnodkbasis4,qnodkbasis5,qnodkbasis6)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      include 'sprods_com.f'
      integer i1,i2
      double precision
     . qnodkbasise0,qnodkbasise1,qnodkbasis0,qnodkbasis1,
     . qnodkbasis2,qnodkbasis3,qnodkbasis4,qnodkbasis5,qnodkbasis6

      qnodkbasise0=  + ro * ( 8.D0 )
      qnodkbasise0 = qnodkbasise0 + 16.D0 - 32.D0*t1 + 32.D0*t1**2

      qnodkbasise1=  - 16.D0

      qnodkbasis0=  + ro * ( 8.D0 )
      qnodkbasis0 = qnodkbasis0 + 16.D0 - 32.D0*t1 + 32.D0*t1**2

      qnodkbasis1=  + ro * ( 2.D0 - 4.D0*t1 )
      qnodkbasis1 = qnodkbasis1 - 8.D0*t1 + 24.D0*t1**2 - 16.D0*t1**3

      qnodkbasis2=  + ro * (  - 8.D0*t1 )
      qnodkbasis2 = qnodkbasis2 - 32.D0*t1**3

      qnodkbasis3=  + ro * ( 80.D0 )
      qnodkbasis3 = qnodkbasis3 + 64.D0 - 128.D0*t1 + 320.D0*t1**2

      qnodkbasis4=  + ro * (  - 8.D0*t1 + 8.D0*t1**2 )
      qnodkbasis4 = qnodkbasis4 + ro**2 * ( 2.D0 )

      qnodkbasis5=  + ro * ( 4.D0 )

      qnodkbasis6=  + ro * ( 8.D0 - 16.D0*t1**2 )
      qnodkbasis6 = qnodkbasis6 + ro**2 * (  - 4.D0 )

      return
      end
