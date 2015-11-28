      subroutine Cfort4(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
      include 'epinv.f'
      include 'scale.f'
C     C0(c,g,0,msq,msq)
      double precision taucg,taucs,taugs,msq
      double precision tcg
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat,C0fb2m
      tcg=taucg+msq
      C0=  + C0fb2m(tcg,msq)

      C1=  + lnrat( - taucg,msq) * ( tcg**(-1) )
      C1 = C1 - C0fb2m(tcg,msq)

      C2=  + lnrat( - taucg,msq) * (  - 2.D0*msq*taucg**(-1)*tcg**(-1)
     &     - tcg**(-1) )
      C2 = C2 + 2.D0*taucg**(-1) + 2.D0*C0fb2m(tcg,msq)*msq*taucg**(-1)

      C11=  + lnrat( - taucg,msq) * ( 1.D0/2.D0*msq*tcg**(-2) - 3.D0/2.D
     &    0*tcg**(-1) )
      C11 = C11 - 1.D0/2.D0*msq*taucg**(-1)*tcg**(-1) + 1.D0/2.D0*
     &    taucg**(-1) + C0fb2m(tcg,msq)

      C12=  + lnrat( - taucg,msq) * ( 4.D0*msq*taucg**(-1)*tcg**(-1) - 
     &    1.D0/2.D0*msq*tcg**(-2) - msq**2*taucg**(-1)*tcg**(-2) + 
     &    tcg**(-1) )
      C12 = C12 - msq*taucg**(-2) + 1.D0/2.D0*msq*taucg**(-1)*tcg**(-1)
     &     + msq**2*taucg**(-2)*tcg**(-1) - 5.D0/2.D0*taucg**(-1) - 3.D0
     &    *C0fb2m(tcg,msq)*msq*taucg**(-1)

      C22=  + lnrat( - taucg,msq) * (  - 2.D0*msq*taucg**(-1)*tcg**(-1)
     &     - 1.D0/2.D0*msq*tcg**(-2) - 6.D0*msq**2*taucg**(-2)*
     &    tcg**(-1) - msq**2*taucg**(-1)*tcg**(-2) + 1.D0/2.D0*
     &    tcg**(-1) )
      C22 = C22 + 5.D0*msq*taucg**(-2) + 1.D0/2.D0*msq*taucg**(-1)*
     &    tcg**(-1) + msq**2*taucg**(-2)*tcg**(-1) - taucg**(-1) + 6.D0
     &    *C0fb2m(tcg,msq)*msq**2*taucg**(-2)

      C00s=  + 1.D0/4.D0

      C00f=  + lnrat( - taucg,msq) * (  - 1.D0/2.D0*msq*tcg**(-1) - 1.D0
     &    /4.D0*taucg*tcg**(-1) )
      C00f = C00f + lnrat(musq,msq) * ( 1.D0/4.D0 )
      C00f = C00f + 3.D0/4.D0 + 1.D0/2.D0*C0fb2m(tcg,msq)*msq

      return
      end
