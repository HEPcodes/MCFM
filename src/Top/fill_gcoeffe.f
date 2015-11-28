      subroutine fill_gcoeffe(
     . gcoeffae1,gcoeffbe1,gcoeffce1,
     . gcoeffae2,gcoeffbe2,gcoeffce2)
      implicit none
      include 'constants.f'
      include 'virtexp.f'
      double precision
     . gcoeffae1,gcoeffbe1,gcoeffce1,
     . gcoeffae2,gcoeffbe2,gcoeffce2

      gcoeffae1 =  - 96.D0*V*b + 128.D0*XLF*b*TR*xn + 96.D0*b*vlsm*xnsq
     &  + 96.D0*b*vltm*xnsq + 96.D0*b*vlwm*xnsq - 352.D0*b*xnsq + 48.D0
     & *ro*vlpm - 96.D0*vlpm

      gcoeffae2 =  - 192.D0*b*xnsq

      gcoeffbe1 = 96.D0*b*vlsm*xnsq - 96.D0*b*vltm*xnsq - 96.D0*b*vlwm*
     & xnsq - 48.D0*ro*vlpm*xnsq + 96.D0*vlpm*xnsq

      gcoeffbe2 =  0

      gcoeffce1 = 96.D0*b*vltm*xnsq - 96.D0*b*vlwm*xnsq

      gcoeffce2 =  0

      return
      end
