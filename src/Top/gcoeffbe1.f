      double precision function gcoeffbe1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeffbe1 = 96.D0*b*vlsm*xnsq - 96.D0*b*vltm*xnsq - 96.D0*b*vlwm*
     & xnsq - 48.D0*ro*vlpm*xnsq + 96.D0*vlpm*xnsq

      return
      end
