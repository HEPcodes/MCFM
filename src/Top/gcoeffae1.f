      double precision function gcoeffae1()
      implicit none
      include 'constants.f'
      include 'virtexp.f'

      gcoeffae1 =  - 96.D0*V*b + 128.D0*XLF*b*TR*xn + 96.D0*b*vlsm*xnsq
     &  + 96.D0*b*vltm*xnsq + 96.D0*b*vlwm*xnsq - 352.D0*b*xnsq + 48.D0
     & *ro*vlpm - 96.D0*vlpm

      return
      end
