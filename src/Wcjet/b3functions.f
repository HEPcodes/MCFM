      subroutine Dforb3(taucg,taucs,taugs,msq,DD3,DD2,DD1,DD0)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision qsq,qsqhat,tcg,qsqtcg,gramdet,ddilog
      double complex DD1,DD2,DD3,DD0,lnrat,Lsm1x2m,Lsm1_2m,C0fa2m,I3me
      qsq=taucg+taucs+taugs+msq
      qsqhat=taucg+taucs+taugs
      gramdet=taucs-msq*taugs/taucg
      tcg=taucg+msq
      qsqtcg=qsqhat-taucg
      Lsm1x2m=Lsm1_2m(taugs,taucg,qsq,msq)
      DD0=  + epinv**2 * ( 3.D0/2.D0*taucg**(-1)*taugs**(-1) )
      DD0 = DD0 + lnrat( - taucg,msq)*epinv * (  - taucs*taucg**(-1)*
     &    taugs**(-1)*qsqtcg**(-1) - taucg**(-1)*taugs**(-1) - 
     &    taucg**(-1)*qsqtcg**(-1) )
      DD0 = DD0 + lnrat( - taucg,msq)**2 * ( taucs*taucg**(-1)*
     &    taugs**(-1)*qsqtcg**(-1) + taucg**(-1)*taugs**(-1) + 
     &    taucg**(-1)*qsqtcg**(-1) )
      DD0 = DD0 + lnrat( - taugs,msq)*epinv * (  - taucg**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat( - taugs,msq)**2 * ( 1.D0/2.D0*taucg**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat( - qsqhat,msq)*epinv * ( taucs*taucg**(-1)*
     &    taugs**(-1)*qsqtcg**(-1) + taucg**(-1)*qsqtcg**(-1) )
      DD0 = DD0 + lnrat( - qsqhat,msq)**2 * (  - taucs*taucg**(-1)*
     &    taugs**(-1)*qsqtcg**(-1) - taucg**(-1)*qsqtcg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*epinv * ( 3.D0/2.D0*taucg**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taucg,msq) * (  - taucs*
     &    taucg**(-1)*taugs**(-1)*qsqtcg**(-1) - taucg**(-1)*
     &    taugs**(-1) - taucg**(-1)*qsqtcg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taugs,msq) * (  - 
     &    taucg**(-1)*taugs**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * ( taucs*
     &    taucg**(-1)*taugs**(-1)*qsqtcg**(-1) + taucg**(-1)*
     &    qsqtcg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)**2 * ( 3.D0/4.D0*taucg**(-1)*
     &    taugs**(-1) )
      DD0 = DD0 + Lsm1x2m * ( taucg**(-1)*taugs**(-1) )
      DD0 = DD0 + C0fa2m(tcg,qsq,msq) * ( taucs*taucg**(-1)*taugs**(-1)
     &     + taucg**(-1) )
      DD0 = DD0 + 1.D0/2.D0*taucg**(-1)*taugs**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcg)*taucg**(-1)*taugs**(-1)

      DD1=  + epinv**2*gramdet**(-1) * ( 3.D0/2.D0*msq*taucg**(-2) - 3.D
     &    0/2.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + gramdet**(-1) * ( 1.D0/2.D0*msq*taucg**(-2)*pisqo6 - 
     &    1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcg)*msq*taucg**(-2) - ddilog(msq**(-1)*tcg)*taucs*
     &    taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + I3me(msq,taugs,qsq)*gramdet**(-1) * (  - msq*
     &    taucg**(-1) + 1.D0/2.D0*taucs*taugs**(-1) + 1.D0/2.D0*taucg*
     &    taugs**(-1) )
      DD1 = DD1 + lnrat( - taucg,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucs*taucg**(-2)*qsqtcg**(-1) - msq*taucg**(-2)*taugs*
     &    qsqtcg**(-1) - msq*taucg**(-2) - msq*taucg**(-1)*qsqtcg**(-1)
     &     + taucs*taucg**(-1)*taugs**(-1) + taucs*taucg**(-1)*
     &    qsqtcg**(-1) + taucs*taugs**(-1)*qsqtcg**(-1) + taucs**2*
     &    taucg**(-1)*taugs**(-1)*qsqtcg**(-1) )
      DD1 = DD1 + lnrat( - taucg,msq)**2*gramdet**(-1) * ( msq*taucs*
     &    taucg**(-2)*qsqtcg**(-1) + msq*taucg**(-2)*taugs*qsqtcg**(-1)
     &     + msq*taucg**(-2) + msq*taucg**(-1)*qsqtcg**(-1) - taucs*
     &    taucg**(-1)*taugs**(-1) - taucs*taucg**(-1)*qsqtcg**(-1) - 
     &    taucs*taugs**(-1)*qsqtcg**(-1) - taucs**2*taucg**(-1)*
     &    taugs**(-1)*qsqtcg**(-1) )
      DD1 = DD1 + lnrat( - taugs,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucg**(-2) + taucs*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + lnrat( - taugs,msq)**2*gramdet**(-1) * ( 1.D0/2.D0*
     &    msq*taucg**(-2) - 1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + lnrat( - qsqhat,msq)*epinv*gramdet**(-1) * ( msq*
     &    taucs*taucg**(-2)*qsqtcg**(-1) + msq*taucg**(-2)*taugs*
     &    qsqtcg**(-1) + msq*taucg**(-1)*qsqtcg**(-1) - taucs*
     &    taucg**(-1)*qsqtcg**(-1) - taucs*taugs**(-1)*qsqtcg**(-1) - 
     &    taucs**2*taucg**(-1)*taugs**(-1)*qsqtcg**(-1) )
      DD1 = DD1 + lnrat( - qsqhat,msq)**2*gramdet**(-1) * (  - msq*
     &    taucs*taucg**(-2)*qsqtcg**(-1) - msq*taucg**(-2)*taugs*
     &    qsqtcg**(-1) - msq*taucg**(-1)*qsqtcg**(-1) + taucs*
     &    taucg**(-1)*qsqtcg**(-1) + taucs*taugs**(-1)*qsqtcg**(-1) + 
     &    taucs**2*taucg**(-1)*taugs**(-1)*qsqtcg**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*epinv*gramdet**(-1) * ( 3.D0/2.D0*msq
     &    *taucg**(-2) - 3.D0/2.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taucg,msq)*gramdet**(-1)
     &  * (  - msq*taucs*taucg**(-2)*qsqtcg**(-1) - msq*taucg**(-2)*
     &    taugs*qsqtcg**(-1) - msq*taucg**(-2) - msq*taucg**(-1)*
     &    qsqtcg**(-1) + taucs*taucg**(-1)*taugs**(-1) + taucs*
     &    taucg**(-1)*qsqtcg**(-1) + taucs*taugs**(-1)*qsqtcg**(-1) + 
     &    taucs**2*taucg**(-1)*taugs**(-1)*qsqtcg**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taugs,msq)*gramdet**(-1)
     &  * (  - msq*taucg**(-2) + taucs*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - qsqhat,msq)*gramdet**(-1)
     &  * ( msq*taucs*taucg**(-2)*qsqtcg**(-1) + msq*taucg**(-2)*taugs*
     &    qsqtcg**(-1) + msq*taucg**(-1)*qsqtcg**(-1) - taucs*
     &    taucg**(-1)*qsqtcg**(-1) - taucs*taugs**(-1)*qsqtcg**(-1) - 
     &    taucs**2*taucg**(-1)*taugs**(-1)*qsqtcg**(-1) )
      DD1 = DD1 + lnrat(musq,msq)**2*gramdet**(-1) * ( 3.D0/4.D0*msq*
     &    taucg**(-2) - 3.D0/4.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD1 = DD1 + Lsm1x2m*gramdet**(-1) * ( msq*taucg**(-2) - taucs*
     &    taucg**(-1)*taugs**(-1) - 1.D0/2.D0*taugs**(-1) )
      DD1 = DD1 + C0fa2m(tcg,qsq,msq)*gramdet**(-1) * ( msq*taucs*
     &    taucg**(-2) + msq*taucg**(-2)*taugs + msq*taucg**(-1) - taucs
     &    *taucg**(-1) - taucs*taugs**(-1) - taucs**2*taucg**(-1)*
     &    taugs**(-1) )

      DD2=  + epinv**2*gramdet**(-1) * ( 1.D0/2.D0*msq*taucg**(-2) - 1.D
     &    0/2.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + gramdet**(-1) * ( 1.D0/2.D0*msq*taucg**(-2)*pisqo6 - 
     &    1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcg)*msq*taucg**(-2) - ddilog(msq**(-1)*tcg)*taucs*
     &    taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + I3me(msq,taugs,qsq)*gramdet**(-1) * ( msq*taucg**(-1)
     &     - 1.D0/2.D0*taucs*taugs**(-1) - 1.D0/2.D0*taucs**2*
     &    taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + lnrat( - taucg,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucg**(-2) + taucs*taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + lnrat( - taucg,msq)**2*gramdet**(-1) * ( msq*
     &    taucg**(-2) - taucs*taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + lnrat(musq,msq)*epinv*gramdet**(-1) * ( 1.D0/2.D0*msq
     &    *taucg**(-2) - 1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + lnrat(musq,msq)*lnrat( - taucg,msq)*gramdet**(-1)
     &  * (  - msq*taucg**(-2) + taucs*taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + lnrat(musq,msq)**2*gramdet**(-1) * ( 1.D0/4.D0*msq*
     &    taucg**(-2) - 1.D0/4.D0*taucs*taucg**(-1)*taugs**(-1) )
      DD2 = DD2 + Lsm1x2m*gramdet**(-1) * (  - 1.D0/2.D0*taucs*
     &    taucg**(-1)*taugs**(-1) )

      DD3=  + I3me(msq,taugs,qsq)*gramdet**(-1) * (  - 1.D0/2.D0 + 1.D0/
     &    2.D0*taucs*taucg**(-1) )
      DD3 = DD3 + Lsm1x2m*gramdet**(-1) * ( 1.D0/2.D0*taucg**(-1) )

      return
      end
