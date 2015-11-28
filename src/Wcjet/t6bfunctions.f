      subroutine Cfort6b(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
C     C0(g,c,0,0,msq)
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision tcg,ddilog
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat
      tcg=taucg+msq
      C0=  + 1.D0/2.D0*epinv**2*taucg**(-1) + 1.D0/2.D0*taucg**(-1)*
     &    pisqo6 + ddilog(msq**(-1)*tcg)*taucg**(-1) - lnrat( - taucg,
     &    msq)*epinv*taucg**(-1) + lnrat( - taucg,msq)**2*taucg**(-1)
     &     + 1.D0/2.D0*lnrat(musq,msq)*epinv*taucg**(-1) - lnrat(musq,
     &    msq)*lnrat( - taucg,msq)*taucg**(-1) + 1.D0/4.D0*lnrat(musq,
     &    msq)**2*taucg**(-1)

      C1=  + tcg**(-1) * ( lnrat( - taucg,msq) + 2.D0*lnrat( - taucg,
     &    msq)*msq*taucg**(-1) )
      C1 = C1 - epinv*taucg**(-1) - 1.D0/2.D0*epinv**2*taucg**(-1) - 1.D
     &    0/2.D0*taucg**(-1)*pisqo6 - 2.D0*taucg**(-1) - ddilog(
     &    msq**(-1)*tcg)*taucg**(-1) + lnrat( - taucg,msq)*epinv*
     &    taucg**(-1) - lnrat( - taucg,msq)**2*taucg**(-1) - 1.D0/2.D0*
     &    lnrat(musq,msq)*epinv*taucg**(-1) - lnrat(musq,msq)*
     &    taucg**(-1) + lnrat(musq,msq)*lnrat( - taucg,msq)*taucg**(-1)
     &     - 1.D0/4.D0*lnrat(musq,msq)**2*taucg**(-1)

      C2=  + tcg**(-1) * (  - lnrat( - taucg,msq) )

      C11=  + tcg**(-2) * ( 1.D0/2.D0*lnrat( - taucg,msq)*msq + lnrat(
     &     - taucg,msq)*msq**2*taucg**(-1) )
      C11 = C11 + tcg**(-1) * (  - 1.D0/2.D0*msq*taucg**(-1) - msq**2*
     &    taucg**(-2) - 3.D0/2.D0*lnrat( - taucg,msq) - 4.D0*lnrat( - 
     &    taucg,msq)*msq*taucg**(-1) )
      C11 = C11 + 3.D0/2.D0*epinv*taucg**(-1) + 1.D0/2.D0*epinv**2*
     &    taucg**(-1) + msq*taucg**(-2) + 1.D0/2.D0*taucg**(-1)*pisqo6
     &     + 3.D0*taucg**(-1) + ddilog(msq**(-1)*tcg)*taucg**(-1) - 
     &    lnrat( - taucg,msq)*epinv*taucg**(-1) + lnrat( - taucg,msq)**
     &    2*taucg**(-1) + 1.D0/2.D0*lnrat(musq,msq)*epinv*taucg**(-1)
     &     + 3.D0/2.D0*lnrat(musq,msq)*taucg**(-1) - lnrat(musq,msq)*
     &    lnrat( - taucg,msq)*taucg**(-1) + 1.D0/4.D0*lnrat(musq,msq)**
     &    2*taucg**(-1)

      C12=  + tcg**(-2) * (  - 1.D0/2.D0*lnrat( - taucg,msq)*msq )
      C12 = C12 + tcg**(-1) * ( 1.D0/2.D0*msq*taucg**(-1) + lnrat( - 
     &    taucg,msq) )
      C12 = C12 - 1.D0/2.D0*taucg**(-1)

      C22=  + tcg**(-2) * (  - 1.D0/2.D0*lnrat( - taucg,msq)*msq )
      C22 = C22 + tcg**(-1) * ( 1.D0/2.D0*msq*taucg**(-1) + 1.D0/2.D0*
     &    lnrat( - taucg,msq) )
      C22 = C22 - 1.D0/2.D0*taucg**(-1)

      C00s=  + 1.D0/4.D0

      C00f=  + tcg**(-1) * (  - 1.D0/4.D0*lnrat( - taucg,msq)*taucg )
      C00f = C00f + 3.D0/4.D0 + 1.D0/4.D0*lnrat(musq,msq)

      return
      end
