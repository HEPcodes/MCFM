      subroutine Cfort3(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
C     C0(s,g,0,0,0)
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat
      C0=  + lnrat( - taugs,msq)**2 * ( 1.D0/2.D0*taugs**(-1) )
      C0 = C0 + lnrat( - taugs,msq) * (  - epinv*taugs**(-1) )
      C0 = C0 + lnrat(musq,msq)*lnrat( - taugs,msq) * (  - taugs**(-1)
     &     )
      C0 = C0 + lnrat(musq,msq)**2 * ( 1.D0/2.D0*taugs**(-1) )
      C0 = C0 + lnrat(musq,msq) * ( epinv*taugs**(-1) )
      C0 = C0 + epinv**2*taugs**(-1)

      C1=  + lnrat( - taugs,msq)**2 * (  - 1.D0/2.D0*taugs**(-1) )
      C1 = C1 + lnrat( - taugs,msq) * ( epinv*taugs**(-1) + taugs**(-1)
     &     )
      C1 = C1 + lnrat(musq,msq)*lnrat( - taugs,msq) * ( taugs**(-1) )
      C1 = C1 + lnrat(musq,msq)**2 * (  - 1.D0/2.D0*taugs**(-1) )
      C1 = C1 + lnrat(musq,msq) * (  - epinv*taugs**(-1) - taugs**(-1)
     &     )
      C1 = C1 - epinv*taugs**(-1) - epinv**2*taugs**(-1) - 2.D0*
     &    taugs**(-1)

      C2=  + lnrat( - taugs,msq) * (  - taugs**(-1) )
      C2 = C2 + lnrat(musq,msq) * ( taugs**(-1) )
      C2 = C2 + epinv*taugs**(-1) + 2.D0*taugs**(-1)

      C11=  + lnrat( - taugs,msq)**2 * ( 1.D0/2.D0*taugs**(-1) )
      C11 = C11 + lnrat( - taugs,msq) * (  - epinv*taugs**(-1) - 3.D0/2.
     &    D0*taugs**(-1) )
      C11 = C11 + lnrat(musq,msq)*lnrat( - taugs,msq) * (  - 
     &    taugs**(-1) )
      C11 = C11 + lnrat(musq,msq)**2 * ( 1.D0/2.D0*taugs**(-1) )
      C11 = C11 + lnrat(musq,msq) * ( epinv*taugs**(-1) + 3.D0/2.D0*
     &    taugs**(-1) )
      C11 = C11 + 3.D0/2.D0*epinv*taugs**(-1) + epinv**2*taugs**(-1) + 
     &    3.D0*taugs**(-1)

      C12=  + lnrat( - taugs,msq) * ( taugs**(-1) )
      C12 = C12 + lnrat(musq,msq) * (  - taugs**(-1) )
      C12 = C12 - epinv*taugs**(-1) - 5.D0/2.D0*taugs**(-1)

      C22=  + lnrat( - taugs,msq) * ( 1.D0/2.D0*taugs**(-1) )
      C22 = C22 + lnrat(musq,msq) * (  - 1.D0/2.D0*taugs**(-1) )
      C22 = C22 - 1.D0/2.D0*epinv*taugs**(-1) - taugs**(-1)

      C00s=  + 1.D0/4.D0

      C00f=  + lnrat( - taugs,msq) * (  - 1.D0/4.D0 )
      C00f = C00f + lnrat(musq,msq) * ( 1.D0/4.D0 )
      C00f = C00f + 3.D0/4.D0

      return
      end
