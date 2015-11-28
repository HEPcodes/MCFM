      double precision function qcoeffe2()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      qcoeffe2 =  + xn**(-1) * ( 16.D0 )
      qcoeffe2 = qcoeffe2 + xn * (  - 16.D0 )

      return
      end
