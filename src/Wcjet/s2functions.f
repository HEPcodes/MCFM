      subroutine Bfors2(taucg,msq,BB1,BB0)
      implicit none
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,msq
      double precision tcg
      double complex BB1,BB0,lnrat
      tcg=taucg+msq
      BB0=  + lnrat( - taucg,msq) * (  - taucg*tcg**(-1) )
      BB0 = BB0 + lnrat(musq,msq) * ( 1.D0 )
      BB0 = BB0 + 2.D0 + epinv

      BB1=  + lnrat( - taucg,msq) * ( 1.D0/2.D0*msq*taucg*tcg**(-2) + 1.
     &    D0/2.D0*taucg*tcg**(-1) )
      BB1 = BB1 + lnrat(musq,msq) * (  - 1.D0/2.D0 )
      BB1 = BB1 - 1.D0 - 1.D0/2.D0*epinv - 1.D0/2.D0*msq*tcg**(-1)

      return
      end
