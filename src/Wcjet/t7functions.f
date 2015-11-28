      subroutine Cfort7(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
C     C0(g,cs,0,0,msq)
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision tcs,qsqtcs,qsq,qsqhat
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat,C0fa2m
      qsq=taucg+taucs+taugs+msq
      qsqhat=taucg+taucs+taugs
      tcs=taucs+msq
      qsqtcs=qsqhat-taucs
      C0=  + lnrat( - taucs,msq)*epinv * ( qsqtcs**(-1) )
      C0 = C0 + lnrat( - taucs,msq)**2 * (  - qsqtcs**(-1) )
      C0 = C0 + lnrat( - qsqhat,msq)*epinv * (  - qsqtcs**(-1) )
      C0 = C0 + lnrat( - qsqhat,msq)**2 * ( qsqtcs**(-1) )
      C0 = C0 + lnrat(musq,msq)*lnrat( - taucs,msq) * ( qsqtcs**(-1) )
      C0 = C0 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * (  - 
     &    qsqtcs**(-1) )
      C0 = C0 + C0fa2m(tcs,qsq,msq) * (  - 1.D0 )

      C1=  + epinv * (  - qsqtcs**(-1) )
      C1 = C1 + lnrat( - taucs,msq)*epinv * (  - taucs*qsqtcs**(-2) - 
     &    qsqtcs**(-1) )
      C1 = C1 + lnrat( - taucs,msq)**2 * ( taucs*qsqtcs**(-2) + 
     &    qsqtcs**(-1) )
      C1 = C1 + lnrat( - taucs,msq) * (  - 2.D0*msq*taucs*tcs**(-1)*
     &    qsqtcs**(-2) - 2.D0*taucs**2*tcs**(-1)*qsqtcs**(-2) )
      C1 = C1 + lnrat( - qsqhat,msq)*epinv * ( taucs*qsqtcs**(-2) + 
     &    qsqtcs**(-1) )
      C1 = C1 + lnrat( - qsqhat,msq)**2 * (  - taucs*qsqtcs**(-2) - 
     &    qsqtcs**(-1) )
      C1 = C1 + lnrat( - qsqhat,msq) * ( 2.D0*msq*taucs*qsq**(-1)*
     &    qsqtcs**(-2) + 2.D0*msq*qsq**(-1)*qsqtcs**(-1) + 3.D0*taucs*
     &    qsq**(-1)*qsqtcs**(-1) + 2.D0*taucs**2*qsq**(-1)*qsqtcs**(-2)
     &     + qsq**(-1) )
      C1 = C1 + lnrat(musq,msq)*lnrat( - taucs,msq) * (  - taucs*
     &    qsqtcs**(-2) - qsqtcs**(-1) )
      C1 = C1 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * ( taucs*
     &    qsqtcs**(-2) + qsqtcs**(-1) )
      C1 = C1 + lnrat(musq,msq) * (  - qsqtcs**(-1) )
      C1 = C1 + C0fa2m(tcs,qsq,msq) * ( 1.D0 + taucs*qsqtcs**(-1) )
      C1 = C1 - 2.D0*qsqtcs**(-1)

      C2=  + lnrat( - taucs,msq) * ( taucs*tcs**(-1)*qsqtcs**(-1) )
      C2 = C2 + lnrat( - qsqhat,msq) * (  - taucs*qsq**(-1)*
     &    qsqtcs**(-1) - qsq**(-1) )

      C11=  + epinv * ( taucs*qsqtcs**(-2) + 3.D0/2.D0*qsqtcs**(-1) )
      C11 = C11 + lnrat( - taucs,msq)*epinv * ( 2.D0*taucs*qsqtcs**(-2)
     &     + taucs**2*qsqtcs**(-3) + qsqtcs**(-1) )
      C11 = C11 + lnrat( - taucs,msq)**2 * (  - 2.D0*taucs*qsqtcs**(-2)
     &     - taucs**2*qsqtcs**(-3) - qsqtcs**(-1) )
      C11 = C11 + lnrat( - taucs,msq) * ( 4.D0*msq*taucs*tcs**(-1)*
     &    qsqtcs**(-2) + 3.D0*msq*taucs**2*tcs**(-1)*qsqtcs**(-3) + 4.D0
     &    *taucs**2*tcs**(-1)*qsqtcs**(-2) + 3.D0*taucs**3*tcs**(-1)*
     &    qsqtcs**(-3) )
      C11 = C11 + lnrat( - qsqhat,msq)*epinv * (  - 2.D0*taucs*
     &    qsqtcs**(-2) - taucs**2*qsqtcs**(-3) - qsqtcs**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq)**2 * ( 2.D0*taucs*qsqtcs**(-2)
     &     + taucs**2*qsqtcs**(-3) + qsqtcs**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq) * (  - 6.D0*msq*taucs*qsq**(-1)*
     &    qsqtcs**(-2) - 3.D0*msq*taucs**2*qsq**(-1)*qsqtcs**(-3) - 7.D0
     &    /2.D0*msq*qsq**(-1)*qsqtcs**(-1) - msq**2*taucs*qsq**(-2)*
     &    qsqtcs**(-2) - 1.D0/2.D0*msq**2*qsq**(-2)*qsqtcs**(-1) + 
     &    msq**2*qsq**(-1)*qsqtcs**(-2) - msq**3*qsq**(-2)*qsqtcs**(-2)
     &     - 13.D0/2.D0*taucs*qsq**(-1)*qsqtcs**(-1) - 8.D0*taucs**2*
     &    qsq**(-1)*qsqtcs**(-2) - 3.D0*taucs**3*qsq**(-1)*qsqtcs**(-3)
     &     - 3.D0/2.D0*qsq**(-1) )
      C11 = C11 + lnrat(musq,msq)*lnrat( - taucs,msq) * ( 2.D0*taucs*
     &    qsqtcs**(-2) + taucs**2*qsqtcs**(-3) + qsqtcs**(-1) )
      C11 = C11 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * (  - 2.D0*
     &    taucs*qsqtcs**(-2) - taucs**2*qsqtcs**(-3) - qsqtcs**(-1) )
      C11 = C11 + lnrat(musq,msq) * ( taucs*qsqtcs**(-2) + 3.D0/2.D0*
     &    qsqtcs**(-1) )
      C11 = C11 + C0fa2m(tcs,qsq,msq) * (  - 1.D0 - 2.D0*taucs*
     &    qsqtcs**(-1) - taucs**2*qsqtcs**(-2) )
      C11 = C11 - msq*taucs*qsq**(-1)*qsqtcs**(-2) - 1.D0/2.D0*msq*
     &    qsq**(-1)*qsqtcs**(-1) + msq*qsqtcs**(-2) - msq**2*qsq**(-1)*
     &    qsqtcs**(-2) + 3.D0*taucs*qsqtcs**(-2) + 3.D0*qsqtcs**(-1)

      C12=  + lnrat( - taucs,msq) * (  - taucs*tcs**(-1)*qsqtcs**(-1)
     &     - 1.D0/2.D0*taucs**2*tcs**(-1)*qsqtcs**(-2) )
      C12 = C12 + lnrat( - qsqhat,msq) * (  - 1.D0/2.D0*msq*qsq**(-1)*
     &    qsqtcs**(-1) + 1.D0/2.D0*msq**2*qsq**(-2)*qsqtcs**(-1) + 3.D0/
     &    2.D0*taucs*qsq**(-1)*qsqtcs**(-1) + 1.D0/2.D0*taucs**2*
     &    qsq**(-1)*qsqtcs**(-2) + qsq**(-1) )
      C12 = C12 + 1.D0/2.D0*msq*qsq**(-1)*qsqtcs**(-1) - 1.D0/2.D0*
     &    qsqtcs**(-1)

      C22=  + lnrat( - taucs,msq) * ( 1.D0/2.D0*msq*taucs*tcs**(-2)*
     &    qsqtcs**(-1) - 1.D0/2.D0*taucs*tcs**(-1)*qsqtcs**(-1) )
      C22 = C22 + lnrat( - qsqhat,msq) * (  - 1.D0/2.D0*msq*qsq**(-1)*
     &    qsqtcs**(-1) + 1.D0/2.D0*msq**2*qsq**(-2)*qsqtcs**(-1) + 1.D0/
     &    2.D0*taucs*qsq**(-1)*qsqtcs**(-1) + 1.D0/2.D0*qsq**(-1) )
      C22 = C22 - 1.D0/2.D0*msq*tcs**(-1)*qsqtcs**(-1) + 1.D0/2.D0*msq*
     &    qsq**(-1)*qsqtcs**(-1)

      C00s=  + 1.D0/4.D0

      C00f=  + lnrat( - taucs,msq) * ( 1.D0/4.D0*taucs**2*tcs**(-1)*
     &    qsqtcs**(-1) )
      C00f = C00f + lnrat( - qsqhat,msq) * (  - 1.D0/2.D0*taucs*
     &    qsq**(-1) - 1.D0/4.D0*taucs**2*qsq**(-1)*qsqtcs**(-1) - 1.D0/
     &    4.D0*taucg*qsq**(-1) - 1.D0/4.D0*taugs*qsq**(-1) )
      C00f = C00f + lnrat(musq,msq) * ( 1.D0/4.D0 )
      C00f = C00f + 3.D0/4.D0

      return
      end
