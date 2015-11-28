      subroutine Cfort1(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
      include 'epinv.f'
      include 'scale.f'
C     C0(c,gs,msq,0,0)
      double precision taucg,taucs,taugs,msq
      double precision qsq,qsqhat,gramdet
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat,I3me
      qsq=taucg+taucs+taugs+msq
      qsqhat=taucg+taucs+taugs
      gramdet=((taucg+taucs)**2-4d0*taugs*msq)
      C0=  + I3me(msq,taugs,qsq) * (  - 1.D0 )

      C1=  + I3me(msq,taugs,qsq)*gramdet**(-1) * ( taucs*taugs + taucg*
     &    taugs )
      C1 = C1 + I3me(msq,taugs,qsq) * ( 1.D0 )
      C1 = C1 + lnrat( - taugs,msq)*gramdet**(-1) * (  - 2.D0*taugs )
      C1 = C1 + lnrat( - qsqhat,msq)*gramdet**(-1) * ( 4.D0*msq*taugs*
     &    qsq**(-1) + 3.D0*taucs*taugs*qsq**(-1) + 3.D0*taucg*taugs*
     &    qsq**(-1) + 2.D0*taugs**2*qsq**(-1) )
      C1 = C1 + lnrat( - qsqhat,msq) * ( qsq**(-1) )

      C2=  + I3me(msq,taugs,qsq)*gramdet**(-1) * (  - 2.D0*msq*taugs )
      C2 = C2 + lnrat( - taugs,msq)*gramdet**(-1) * ( taucs + taucg )
      C2 = C2 + lnrat( - qsqhat,msq)*gramdet**(-1) * (  - 2.D0*msq*
     &    taucs*qsq**(-1) - 2.D0*msq*taucg*qsq**(-1) - 6.D0*msq*taugs*
     &    qsq**(-1) - taucs*taugs*qsq**(-1) - taucg*taugs*qsq**(-1) )
      C2 = C2 + lnrat( - qsqhat,msq) * (  - qsq**(-1) )

      C11=  + gramdet**(-1) * ( 1.D0/2.D0*msq*taucs*qsq**(-1) + 1.D0/2.D
     &    0*msq*taucg*qsq**(-1) + msq*taugs*qsq**(-1) - 1.D0/2.D0*taucs
     &     - 1.D0/2.D0*taucg + taugs )
      C11 = C11 + I3me(msq,taugs,qsq)*gramdet**(-2) * (  - 6.D0*msq*
     &    taugs**3 )
      C11 = C11 + I3me(msq,taugs,qsq)*gramdet**(-1) * (  - 2.D0*taucs*
     &    taugs - 2.D0*taucg*taugs - taugs**2 )
      C11 = C11 + I3me(msq,taugs,qsq) * (  - 1.D0 )
      C11 = C11 + lnrat( - taugs,msq)*gramdet**(-2) * ( 3.D0*taucs*
     &    taugs**2 + 3.D0*taucg*taugs**2 )
      C11 = C11 + lnrat( - taugs,msq)*gramdet**(-1) * ( 4.D0*taugs )
      C11 = C11 + lnrat( - qsqhat,msq)*gramdet**(-2) * (  - 6.D0*msq*
     &    taucs*taugs**2*qsq**(-1) - 6.D0*msq*taucg*taugs**2*qsq**(-1)
     &     - 18.D0*msq*taugs**3*qsq**(-1) - 3.D0*taucs*taugs**3*
     &    qsq**(-1) - 3.D0*taucg*taugs**3*qsq**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq)*gramdet**(-1) * (  - 1.D0/2.D0*
     &    msq*taucs*qsq**(-1) - 1.D0/2.D0*msq*taucg*qsq**(-1) - 7.D0*
     &    msq*taugs*qsq**(-1) + 1.D0/2.D0*msq**2*taucs*qsq**(-2) + 1.D0/
     &    2.D0*msq**2*taucg*qsq**(-2) + msq**2*taugs*qsq**(-2) - 13.D0/
     &    2.D0*taucs*taugs*qsq**(-1) - 13.D0/2.D0*taucg*taugs*qsq**(-1)
     &     - 8.D0*taugs**2*qsq**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq) * (  - 3.D0/2.D0*qsq**(-1) )

      C12=  + gramdet**(-1) * (  - 1.D0/2.D0*msq*taucs*qsq**(-1) - 1.D0/
     &    2.D0*msq*taucg*qsq**(-1) + msq - msq**2*qsq**(-1) - 1.D0/2.D0
     &    *taucs - 1.D0/2.D0*taucg )
      C12 = C12 + I3me(msq,taugs,qsq)*gramdet**(-2) * ( 3.D0*msq*taucs*
     &    taugs**2 + 3.D0*msq*taucg*taugs**2 )
      C12 = C12 + I3me(msq,taugs,qsq)*gramdet**(-1) * ( 2.D0*msq*taugs
     &     )
      C12 = C12 + lnrat( - taugs,msq)*gramdet**(-2) * (  - 6.D0*msq*
     &    taugs**2 )
      C12 = C12 + lnrat( - taugs,msq)*gramdet**(-1) * (  - taucs - 
     &    taucg - 1.D0/2.D0*taugs )
      C12 = C12 + lnrat( - qsqhat,msq)*gramdet**(-2) * ( 9.D0*msq*taucs
     &    *taugs**2*qsq**(-1) + 9.D0*msq*taucg*taugs**2*qsq**(-1) + 6.D0
     &    *msq*taugs**3*qsq**(-1) + 12.D0*msq**2*taugs**2*qsq**(-1) )
      C12 = C12 + lnrat( - qsqhat,msq)*gramdet**(-1) * ( 3.D0/2.D0*msq*
     &    taucs*qsq**(-1) + 3.D0/2.D0*msq*taucg*qsq**(-1) + 8.D0*msq*
     &    taugs*qsq**(-1) - 1.D0/2.D0*msq**2*taucs*qsq**(-2) - 1.D0/2.D0
     &    *msq**2*taucg*qsq**(-2) + msq**2*qsq**(-1) - msq**3*qsq**(-2)
     &     + 3.D0/2.D0*taucs*taugs*qsq**(-1) + 3.D0/2.D0*taucg*taugs*
     &    qsq**(-1) + 1.D0/2.D0*taugs**2*qsq**(-1) )
      C12 = C12 + lnrat( - qsqhat,msq) * ( qsq**(-1) )

      C22=  + gramdet**(-1) * (  - 1.D0/2.D0*msq*taucs*qsq**(-1) - 1.D0/
     &    2.D0*msq*taucg*qsq**(-1) + msq - msq**2*qsq**(-1) )
      C22 = C22 + I3me(msq,taugs,qsq)*gramdet**(-2) * (  - 6.D0*msq**2*
     &    taugs**2 )
      C22 = C22 + lnrat( - taugs,msq)*gramdet**(-2) * ( 3.D0*msq*taucs*
     &    taugs + 3.D0*msq*taucg*taugs )
      C22 = C22 + lnrat( - taugs,msq)*gramdet**(-1) * (  - 1.D0/2.D0*
     &    taucs - 1.D0/2.D0*taucg )
      C22 = C22 + lnrat( - qsqhat,msq)*gramdet**(-2) * (  - 3.D0*msq*
     &    taucs*taugs**2*qsq**(-1) - 3.D0*msq*taucg*taugs**2*qsq**(-1)
     &     - 6.D0*msq**2*taucs*taugs*qsq**(-1) - 6.D0*msq**2*taucg*
     &    taugs*qsq**(-1) - 18.D0*msq**2*taugs**2*qsq**(-1) )
      C22 = C22 + lnrat( - qsqhat,msq)*gramdet**(-1) * ( 1.D0/2.D0*msq*
     &    taucs*qsq**(-1) + 1.D0/2.D0*msq*taucg*qsq**(-1) - msq*taugs*
     &    qsq**(-1) - 1.D0/2.D0*msq**2*taucs*qsq**(-2) - 1.D0/2.D0*
     &    msq**2*taucg*qsq**(-2) + msq**2*qsq**(-1) - msq**3*qsq**(-2)
     &     + 1.D0/2.D0*taucs*taugs*qsq**(-1) + 1.D0/2.D0*taucg*taugs*
     &    qsq**(-1) )
      C22 = C22 + lnrat( - qsqhat,msq) * ( 1.D0/2.D0*qsq**(-1) )

      C00s=  + 1.D0/4.D0

      C00f=  + I3me(msq,taugs,qsq)*gramdet**(-1) * (  - 1.D0/2.D0*msq*
     &    taugs**2 )
      C00f = C00f + lnrat( - taugs,msq)*gramdet**(-1) * ( 1.D0/4.D0*
     &    taucs*taugs + 1.D0/4.D0*taucg*taugs )
      C00f = C00f + lnrat( - qsqhat,msq)*gramdet**(-1) * (  - 1.D0/2.D0
     &    *msq*taucs*taugs*qsq**(-1) - 1.D0/2.D0*msq*taucg*taugs*
     &    qsq**(-1) - 3.D0/2.D0*msq*taugs**2*qsq**(-1) - 1.D0/4.D0*
     &    taucs*taugs**2*qsq**(-1) - 1.D0/4.D0*taucg*taugs**2*qsq**(-1)
     &     )
      C00f = C00f + lnrat( - qsqhat,msq) * (  - 1.D0/4.D0*taucs*
     &    qsq**(-1) - 1.D0/4.D0*taucg*qsq**(-1) - 1.D0/2.D0*taugs*
     &    qsq**(-1) )
      C00f = C00f + lnrat(musq,msq) * ( 1.D0/4.D0 )
      C00f = C00f + 3.D0/4.D0

      return
      end
