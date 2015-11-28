      subroutine Cfort8b(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
C     C0(c,s,msq,0,0)
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision tcs,ddilog
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat
      tcs=taucs+msq
      C0=  + epinv**2 * ( 1.D0/2.D0*taucs**(-1) )
      C0 = C0 + lnrat( - taucs,msq)*epinv * (  - taucs**(-1) )
      C0 = C0 + lnrat( - taucs,msq)**2 * ( taucs**(-1) )
      C0 = C0 + lnrat(musq,msq)*epinv * ( 1.D0/2.D0*taucs**(-1) )
      C0 = C0 + lnrat(musq,msq)*lnrat( - taucs,msq) * (  - taucs**(-1)
     &     )
      C0 = C0 + lnrat(musq,msq)**2 * ( 1.D0/4.D0*taucs**(-1) )
      C0 = C0 + 1.D0/2.D0*taucs**(-1)*pisqo6 + ddilog(msq**(-1)*tcs)*
     &    taucs**(-1)

      C1=  + epinv**2 * (  - 1.D0/2.D0*taucs**(-1) )
      C1 = C1 + lnrat( - taucs,msq)*epinv * ( taucs**(-1) )
      C1 = C1 + lnrat( - taucs,msq)**2 * (  - taucs**(-1) )
      C1 = C1 + lnrat( - taucs,msq) * ( tcs**(-1) )
      C1 = C1 + lnrat(musq,msq)*epinv * (  - 1.D0/2.D0*taucs**(-1) )
      C1 = C1 + lnrat(musq,msq)*lnrat( - taucs,msq) * ( taucs**(-1) )
      C1 = C1 + lnrat(musq,msq)**2 * (  - 1.D0/4.D0*taucs**(-1) )
      C1 = C1 - 1.D0/2.D0*taucs**(-1)*pisqo6 - ddilog(msq**(-1)*tcs)*
     &    taucs**(-1)

      C2=  + epinv * ( taucs**(-1) )
      C2 = C2 + lnrat( - taucs,msq) * (  - 2.D0*msq*taucs**(-1)*
     &    tcs**(-1) - tcs**(-1) )
      C2 = C2 + lnrat(musq,msq) * ( taucs**(-1) )
      C2 = C2 + 2.D0*taucs**(-1)

      C11=  + epinv**2 * ( 1.D0/2.D0*taucs**(-1) )
      C11 = C11 + lnrat( - taucs,msq)*epinv * (  - taucs**(-1) )
      C11 = C11 + lnrat( - taucs,msq)**2 * ( taucs**(-1) )
      C11 = C11 + lnrat( - taucs,msq) * (  - 1.D0/2.D0*msq*tcs**(-2) - 
     &    3.D0/2.D0*tcs**(-1) )
      C11 = C11 + lnrat(musq,msq)*epinv * ( 1.D0/2.D0*taucs**(-1) )
      C11 = C11 + lnrat(musq,msq)*lnrat( - taucs,msq) * (  - 
     &    taucs**(-1) )
      C11 = C11 + lnrat(musq,msq)**2 * ( 1.D0/4.D0*taucs**(-1) )
      C11 = C11 + 1.D0/2.D0*msq*taucs**(-1)*tcs**(-1) + 1.D0/2.D0*
     &    taucs**(-1)*pisqo6 - 1.D0/2.D0*taucs**(-1) + ddilog(msq**(-1)
     &    *tcs)*taucs**(-1)

      C12=  + epinv * (  - taucs**(-1) )
      C12 = C12 + lnrat( - taucs,msq) * ( msq*taucs**(-1)*tcs**(-1) + 1.
     &    D0/2.D0*msq*tcs**(-2) + msq**2*taucs**(-1)*tcs**(-2) + 
     &    tcs**(-1) )
      C12 = C12 + lnrat(musq,msq) * (  - taucs**(-1) )
      C12 = C12 + msq*taucs**(-2) - 1.D0/2.D0*msq*taucs**(-1)*tcs**(-1)
     &     - msq**2*taucs**(-2)*tcs**(-1) - 5.D0/2.D0*taucs**(-1)

      C22=  + epinv * (  - 1.D0/2.D0*taucs**(-1) )
      C22 = C22 + lnrat( - taucs,msq) * ( 1.D0/2.D0*msq*tcs**(-2) + 
     &    msq**2*taucs**(-1)*tcs**(-2) + 1.D0/2.D0*tcs**(-1) )
      C22 = C22 + lnrat(musq,msq) * (  - 1.D0/2.D0*taucs**(-1) )
      C22 = C22 + msq*taucs**(-2) - 1.D0/2.D0*msq*taucs**(-1)*tcs**(-1)
     &     - msq**2*taucs**(-2)*tcs**(-1) - taucs**(-1)

      C00s=  + 1.D0/4.D0

      C00f=  + lnrat( - taucs,msq) * (  - 1.D0/4.D0*taucs*tcs**(-1) )
      C00f = C00f + lnrat(musq,msq) * ( 1.D0/4.D0 )
      C00f = C00f + 3.D0/4.D0

      return
      end
