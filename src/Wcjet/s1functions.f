      subroutine Bfors1(taugs,msq,BB1,BB0)
      implicit none
      include 'epinv.f'
      include 'scale.f'
      double precision taugs,msq
      double complex BB1,BB0,lnrat
      BB0=  + lnrat( - taugs,msq) * (  - 1.D0 )
      BB0 = BB0 + lnrat(musq,msq) * ( 1.D0 )
      BB0 = BB0 + 2.D0 + epinv

      BB1=  + lnrat( - taugs,msq) * ( 1.D0/2.D0 )
      BB1 = BB1 + lnrat(musq,msq) * (  - 1.D0/2.D0 )
      BB1 = BB1 - 1.D0 - 1.D0/2.D0*epinv

      return
      end
