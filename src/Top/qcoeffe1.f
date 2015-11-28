      double precision function qcoeffe1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      qcoeffe1 =  + b**(-1)*vlpm*xn**(-1) * (  - 16.D0 + 8.D0*ro )
      qcoeffe1 = qcoeffe1 + vlsm*xn**(-1) * (  - 16.D0 )
      qcoeffe1 = qcoeffe1 + vltm*xn**(-1) * (  - 64.D0 )
      qcoeffe1 = qcoeffe1 + vltm*xn * ( 32.D0 )
      qcoeffe1 = qcoeffe1 + vlwm*xn**(-1) * ( 64.D0 )
      qcoeffe1 = qcoeffe1 + xn**(-1) * ( 40.D0 )
      qcoeffe1 = qcoeffe1 + xn * (  - 40.D0 )

      return
      end
