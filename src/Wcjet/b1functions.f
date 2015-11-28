      subroutine Dforb1(taucg,taucs,taugs,msq,DD3,DD2,DD1,DD0)
      implicit none
C      D0(g,s,c,0,0,0,msq)
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision qsq,qsqhat,tcs,qsqtcs,gramdet,ddilog
      double complex DD1,DD2,DD3,DD0,lnrat,Lsm1x2m,Lsm1_2m,C0fa2m,I3me
      qsq=taucg+taucs+taugs+msq
      qsqhat=taucg+taucs+taugs
      gramdet=taucs-msq*taugs/taucg
      qsqtcs=qsqhat-taucs
      tcs=taucs+msq
      Lsm1x2m=Lsm1_2m(taugs,taucs,qsq,msq)
      DD0=  + epinv**2 * ( 3.D0/2.D0*taucs**(-1)*taugs**(-1) )
      DD0 = DD0 + lnrat( - taucs,msq)*epinv*qsqtcs**(-1) * (  - 
     &    taucs**(-1)*taucg*taugs**(-1) - taucs**(-1) )
      DD0 = DD0 + lnrat( - taucs,msq)*epinv * (  - taucs**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat( - taucs,msq)**2*qsqtcs**(-1) * ( taucs**(-1)*
     &    taucg*taugs**(-1) + taucs**(-1) )
      DD0 = DD0 + lnrat( - taucs,msq)**2 * ( taucs**(-1)*taugs**(-1) )
      DD0 = DD0 + lnrat( - taugs,msq)*epinv * (  - taucs**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat( - taugs,msq)**2 * ( 1.D0/2.D0*taucs**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat( - qsqhat,msq)*epinv*qsqtcs**(-1) * ( 
     &    taucs**(-1)*taucg*taugs**(-1) + taucs**(-1) )
      DD0 = DD0 + lnrat( - qsqhat,msq)**2*qsqtcs**(-1) * (  - 
     &    taucs**(-1)*taucg*taugs**(-1) - taucs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*epinv * ( 3.D0/2.D0*taucs**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taucs,msq)*qsqtcs**(-1) * ( 
     &     - taucs**(-1)*taucg*taugs**(-1) - taucs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taucs,msq) * (  - 
     &    taucs**(-1)*taugs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taugs,msq) * (  - 
     &    taucs**(-1)*taugs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - qsqhat,msq)*qsqtcs**(-1)
     &  * ( taucs**(-1)*taucg*taugs**(-1) + taucs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)**2 * ( 3.D0/4.D0*taucs**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + Lsm1x2m * ( taucs**(-1)*taugs**(-1) )
      DD0 = DD0 + 1.D0/2.D0*taucs**(-1)*taugs**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcs)*taucs**(-1)*taugs**(-1) + C0fa2m(tcs,qsq,msq)*
     &    taucs**(-1)*taucg*taugs**(-1) + C0fa2m(tcs,qsq,msq)*
     &    taucs**(-1)

      DD1=  + epinv**2*gramdet**(-1) * ( 3.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-1) - 3.D0/2.D0*taugs**(-1) )
      DD1 = DD1 + gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-1)*pisqo6 - 1.D0/2.D0*taugs**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcs)*msq*taucs**(-1)*taucg**(-1) - ddilog(msq**(-1)
     &    *tcs)*taugs**(-1) - C0fa2m(tcs,qsq,msq) + C0fa2m(tcs,qsq,msq)
     &    *msq*taucs**(-1)*taucg**(-1)*taugs + C0fa2m(tcs,qsq,msq)*msq*
     &    taucs**(-1) + C0fa2m(tcs,qsq,msq)*msq*taucg**(-1) - C0fa2m(
     &    tcs,qsq,msq)*taucs*taugs**(-1) - C0fa2m(tcs,qsq,msq)*taucg*
     &    taugs**(-1) )
      DD1 = DD1 + I3me(msq,taugs,qsq)*gramdet**(-1) * (  - msq*
     &    taucg**(-1) + 1.D0/2.D0*taucs*taugs**(-1) + 1.D0/2.D0*
     &    taucs**2*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + lnrat( - taucs,msq)*epinv*gramdet**(-1)*qsqtcs**(-1)
     &  * ( 1.D0 - msq*taucs**(-1)*taucg**(-1)*taugs - msq*taucs**(-1)
     &     - msq*taucg**(-1) + taucs*taugs**(-1) + taucg*taugs**(-1) )
      DD1 = DD1 + lnrat( - taucs,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucs**(-1)*taucg**(-1) + taugs**(-1) )
      DD1 = DD1 + lnrat( - taucs,msq)**2*gramdet**(-1)*qsqtcs**(-1)
     &  * (  - 1.D0 + msq*taucs**(-1)*taucg**(-1)*taugs + msq*
     &    taucs**(-1) + msq*taucg**(-1) - taucs*taugs**(-1) - taucg*
     &    taugs**(-1) )
      DD1 = DD1 + lnrat( - taucs,msq)**2*gramdet**(-1) * ( msq*
     &    taucs**(-1)*taucg**(-1) - taugs**(-1) )
      DD1 = DD1 + lnrat( - taugs,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucs**(-1)*taucg**(-1) + taugs**(-1) )
      DD1 = DD1 + lnrat( - taugs,msq)**2*gramdet**(-1) * ( 1.D0/2.D0*
     &    msq*taucs**(-1)*taucg**(-1) - 1.D0/2.D0*taugs**(-1) )
      DD1 = DD1 + lnrat( - qsqhat,msq)*epinv*gramdet**(-1)*qsqtcs**(-1)
     &  * (  - 1.D0 + msq*taucs**(-1)*taucg**(-1)*taugs + msq*
     &    taucs**(-1) + msq*taucg**(-1) - taucs*taugs**(-1) - taucg*
     &    taugs**(-1) )
      DD1 = DD1 + lnrat( - qsqhat,msq)**2*gramdet**(-1)*qsqtcs**(-1)
     &  * ( 1.D0 - msq*taucs**(-1)*taucg**(-1)*taugs - msq*taucs**(-1)
     &     - msq*taucg**(-1) + taucs*taugs**(-1) + taucg*taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*epinv*gramdet**(-1) * ( 3.D0/2.D0*msq
     &    *taucs**(-1)*taucg**(-1) - 3.D0/2.D0*taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taucs,msq)*gramdet**(-1)*
     & qsqtcs**(-1) * ( 1.D0 - msq*taucs**(-1)*taucg**(-1)*taugs - msq*
     &    taucs**(-1) - msq*taucg**(-1) + taucs*taugs**(-1) + taucg*
     &    taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taucs,msq)*gramdet**(-1)
     &  * (  - msq*taucs**(-1)*taucg**(-1) + taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taugs,msq)*gramdet**(-1)
     &  * (  - msq*taucs**(-1)*taucg**(-1) + taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - qsqhat,msq)*gramdet**(-1)*
     & qsqtcs**(-1) * (  - 1.D0 + msq*taucs**(-1)*taucg**(-1)*taugs + 
     &    msq*taucs**(-1) + msq*taucg**(-1) - taucs*taugs**(-1) - taucg
     &    *taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)**2*gramdet**(-1) * ( 3.D0/4.D0*msq*
     &    taucs**(-1)*taucg**(-1) - 3.D0/4.D0*taugs**(-1) )
      DD1 = DD1 + Lsm1x2m*gramdet**(-1) * ( msq*taucs**(-1)*taucg**(-1)
     &     - 1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1) - taugs**(-1) )

      DD2=  + epinv**2*gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-1) - 1.D0/2.D0*taugs**(-1) )
      DD2 = DD2 + gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-1)*pisqo6 - 1.D0/2.D0*taugs**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcs)*msq*taucs**(-1)*taucg**(-1) - ddilog(msq**(-1)
     &    *tcs)*taugs**(-1) )
      DD2 = DD2 + I3me(msq,taugs,qsq)*gramdet**(-1) * ( msq*taucg**(-1)
     &     - 1.D0/2.D0*taucs*taugs**(-1) - 1.D0/2.D0*taucg*taugs**(-1)
     &     )
      DD2 = DD2 + lnrat( - taucs,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucs**(-1)*taucg**(-1) + taugs**(-1) )
      DD2 = DD2 + lnrat( - taucs,msq)**2*gramdet**(-1) * ( msq*
     &    taucs**(-1)*taucg**(-1) - taugs**(-1) )
      DD2 = DD2 + lnrat(musq,msq)*epinv*gramdet**(-1) * ( 1.D0/2.D0*msq
     &    *taucs**(-1)*taucg**(-1) - 1.D0/2.D0*taugs**(-1) )
      DD2 = DD2 + lnrat(musq,msq)*lnrat( - taucs,msq)*gramdet**(-1)
     &  * (  - msq*taucs**(-1)*taucg**(-1) + taugs**(-1) )
      DD2 = DD2 + lnrat(musq,msq)**2*gramdet**(-1) * ( 1.D0/4.D0*msq*
     &    taucs**(-1)*taucg**(-1) - 1.D0/4.D0*taugs**(-1) )
      DD2 = DD2 + Lsm1x2m*gramdet**(-1) * (  - 1.D0/2.D0*taugs**(-1) )

      DD3=  + I3me(msq,taugs,qsq)*gramdet**(-1) * ( 1.D0/2.D0 - 1.D0/2.D
     &    0*taucs*taucg**(-1) )
      DD3 = DD3 + Lsm1x2m*gramdet**(-1) * ( 1.D0/2.D0*taucg**(-1) )

      return
      end
