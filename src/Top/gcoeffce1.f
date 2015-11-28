      double precision function gcoeffce1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeffce1 = 96.D0*b*vltm*xnsq - 96.D0*b*vlwm*xnsq

      return
      end
