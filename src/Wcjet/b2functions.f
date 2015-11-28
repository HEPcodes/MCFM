      subroutine Dforb2(taucg,taucs,taugs,msq,DD3,DD2,DD1,DD0)
      implicit none
C     D0(g,c,s,msq,msq,0,0)
      include 'constants.f'
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision qsq,qsqhat,tcg,tcs,gramdet,ddilog
      double complex DD1,DD2,DD3,DD0,lnrat,Lsm2_2m,Lsm2x2m,C0fa2m,C0fb2m
      qsq=taucg+taucs+taugs+msq
      qsqhat=taucg+taucs+taugs
      gramdet=taucs-msq*taugs/taucg
      tcg=taucg+msq
      tcs=taucs+msq
      Lsm2x2m=Lsm2_2m(taucs,taucg,qsq,msq)
      DD0=  + epinv**2 * ( 1.D0/2.D0*taucs**(-1)*taucg**(-1) )
      DD0 = DD0 + lnrat( - taucs,msq)*epinv * (  - taucs**(-1)*
     &    taucg**(-1) )
      DD0 = DD0 + lnrat( - taucs,msq)**2 * ( taucs**(-1)*taucg**(-1) )
      DD0 = DD0 + lnrat( - taucg,msq)*epinv * (  - 1/( - taucg + qsqhat
     &    )*taucs**(-1)*taucg**(-1)*taugs - 1/( - taucg + qsqhat)*
     &    taucg**(-1) )
      DD0 = DD0 + lnrat( - taucg,msq)**2 * ( 1/( - taucg + qsqhat)*
     &    taucs**(-1)*taucg**(-1)*taugs + 1/( - taucg + qsqhat)*
     &    taucg**(-1) )
      DD0 = DD0 + lnrat( - qsqhat,msq)*epinv * ( 1/( - taucg + qsqhat)*
     &    taucs**(-1)*taucg**(-1)*taugs + 1/( - taucg + qsqhat)*
     &    taucg**(-1) )
      DD0 = DD0 + lnrat( - qsqhat,msq)**2 * (  - 1/( - taucg + qsqhat)*
     &    taucs**(-1)*taucg**(-1)*taugs - 1/( - taucg + qsqhat)*
     &    taucg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*epinv * ( 1.D0/2.D0*taucs**(-1)*
     &    taucg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taucs,msq) * (  - 
     &    taucs**(-1)*taucg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - taucg,msq) * (  - 1/( - 
     &    taucg + qsqhat)*taucs**(-1)*taucg**(-1)*taugs - 1/( - taucg
     &     + qsqhat)*taucg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * ( 1/( - taucg
     &     + qsqhat)*taucs**(-1)*taucg**(-1)*taugs + 1/( - taucg + 
     &    qsqhat)*taucg**(-1) )
      DD0 = DD0 + lnrat(musq,msq)**2 * ( 1.D0/4.D0*taucs**(-1)*
     &    taucg**(-1) )
      DD0 = DD0 + Lsm2x2m * ( taucs**(-1)*taucg**(-1) )
      DD0 = DD0 + C0fa2m(tcg,qsq,msq) * ( taucs**(-1)*taucg**(-1)*taugs
     &     + taucg**(-1) )
      DD0 = DD0 + 1.D0/2.D0*taucs**(-1)*taucg**(-1)*pisqo6 + ddilog(
     &    msq**(-1)*tcs)*taucs**(-1)*taucg**(-1)

      DD1=  + epinv**2*gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-2)*taugs - 1.D0/2.D0*taucg**(-1) )
      DD1 = DD1 + gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-2)*taugs*pisqo6 - 1.D0/2.D0*taucg**(-1)*pisqo6 + 
     &    ddilog(msq**(-1)*tcs)*msq*taucs**(-1)*taucg**(-2)*taugs - 
     &    ddilog(msq**(-1)*tcs)*taucg**(-1) )
      DD1 = DD1 + lnrat( - taucs,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucs**(-1)*taucg**(-2)*taugs + taucg**(-1) )
      DD1 = DD1 + lnrat( - taucs,msq)**2*gramdet**(-1) * ( msq*
     &    taucs**(-1)*taucg**(-2)*taugs - taucg**(-1) )
      DD1 = DD1 + lnrat( - taucg,msq)*epinv*gramdet**(-1) * (  - 1/( - 
     &    taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 - 1/( - 
     &    taucg + qsqhat)*msq*taucg**(-2)*taugs + 1/( - taucg + qsqhat)
     &    *taucs*taucg**(-1) + 1/( - taucg + qsqhat)*taucg**(-1)*taugs
     &     )
      DD1 = DD1 + lnrat( - taucg,msq)**2*gramdet**(-1) * ( 1/( - taucg
     &     + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 + 1/( - taucg
     &     + qsqhat)*msq*taucg**(-2)*taugs - 1/( - taucg + qsqhat)*
     &    taucs*taucg**(-1) - 1/( - taucg + qsqhat)*taucg**(-1)*taugs )
      DD1 = DD1 + lnrat( - qsqhat,msq)*epinv*gramdet**(-1) * ( 1/( - 
     &    taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 + 1/( - 
     &    taucg + qsqhat)*msq*taucg**(-2)*taugs - 1/( - taucg + qsqhat)
     &    *taucs*taucg**(-1) - 1/( - taucg + qsqhat)*taucg**(-1)*taugs
     &     )
      DD1 = DD1 + lnrat( - qsqhat,msq)**2*gramdet**(-1) * (  - 1/( - 
     &    taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 - 1/( - 
     &    taucg + qsqhat)*msq*taucg**(-2)*taugs + 1/( - taucg + qsqhat)
     &    *taucs*taucg**(-1) + 1/( - taucg + qsqhat)*taucg**(-1)*taugs
     &     )
      DD1 = DD1 + lnrat(musq,msq)*epinv*gramdet**(-1) * ( 1.D0/2.D0*msq
     &    *taucs**(-1)*taucg**(-2)*taugs - 1.D0/2.D0*taucg**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taucs,msq)*gramdet**(-1)
     &  * (  - msq*taucs**(-1)*taucg**(-2)*taugs + taucg**(-1) )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - taucg,msq)*gramdet**(-1)
     &  * (  - 1/( - taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*
     &    taugs**2 - 1/( - taucg + qsqhat)*msq*taucg**(-2)*taugs + 1/(
     &     - taucg + qsqhat)*taucs*taucg**(-1) + 1/( - taucg + qsqhat)*
     &    taucg**(-1)*taugs )
      DD1 = DD1 + lnrat(musq,msq)*lnrat( - qsqhat,msq)*gramdet**(-1)
     &  * ( 1/( - taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2
     &     + 1/( - taucg + qsqhat)*msq*taucg**(-2)*taugs - 1/( - taucg
     &     + qsqhat)*taucs*taucg**(-1) - 1/( - taucg + qsqhat)*
     &    taucg**(-1)*taugs )
      DD1 = DD1 + lnrat(musq,msq)**2*gramdet**(-1) * ( 1.D0/4.D0*msq*
     &    taucs**(-1)*taucg**(-2)*taugs - 1.D0/4.D0*taucg**(-1) )
      DD1 = DD1 + Lsm2x2m*gramdet**(-1) * ( msq*taucs**(-1)*taucg**(-2)
     &    *taugs - 1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1) - 
     &    taucg**(-1) )
      DD1 = DD1 + C0fb2m(tcg,msq)*gramdet**(-1) * (  - msq*taucg**(-1)
     &     + 1.D0/2.D0*taucs*taugs**(-1) )
      DD1 = DD1 + C0fa2m(tcs,qsq,msq)*gramdet**(-1) * ( msq*taucg**(-1)
     &     + 1.D0/2.D0*taucs*taucg**(-1) - 1.D0/2.D0*taucs*taugs**(-1)
     &     )
      DD1 = DD1 + C0fa2m(tcg,qsq,msq)*gramdet**(-1) * ( msq*taucs**(-1)
     &    *taucg**(-2)*taugs**2 + msq*taucg**(-2)*taugs - taucs*
     &    taucg**(-1) - taucg**(-1)*taugs )

      DD2=  + epinv**2*gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-2)*taugs - 1.D0/2.D0*taucg**(-1) )
      DD2 = DD2 + gramdet**(-1) * ( 1.D0/2.D0*msq*taucs**(-1)*
     &    taucg**(-2)*taugs*pisqo6 - 1.D0/2.D0*taucg**(-1)*pisqo6 + 
     &    ddilog(msq**(-1)*tcs)*msq*taucs**(-1)*taucg**(-2)*taugs - 
     &    ddilog(msq**(-1)*tcs)*taucg**(-1) )
      DD2 = DD2 + lnrat( - taucs,msq)*epinv*gramdet**(-1) * (  - msq*
     &    taucs**(-1)*taucg**(-2)*taugs + taucg**(-1) )
      DD2 = DD2 + lnrat( - taucs,msq)**2*gramdet**(-1) * ( msq*
     &    taucs**(-1)*taucg**(-2)*taugs - taucg**(-1) )
      DD2 = DD2 + lnrat( - taucg,msq)*epinv*gramdet**(-1) * (  - 1/( - 
     &    taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 - 1/( - 
     &    taucg + qsqhat)*msq*taucg**(-2)*taugs + 1/( - taucg + qsqhat)
     &    *taucs*taucg**(-1) + 1/( - taucg + qsqhat)*taucg**(-1)*taugs
     &     )
      DD2 = DD2 + lnrat( - taucg,msq)**2*gramdet**(-1) * ( 1/( - taucg
     &     + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 + 1/( - taucg
     &     + qsqhat)*msq*taucg**(-2)*taugs - 1/( - taucg + qsqhat)*
     &    taucs*taucg**(-1) - 1/( - taucg + qsqhat)*taucg**(-1)*taugs )
      DD2 = DD2 + lnrat( - qsqhat,msq)*epinv*gramdet**(-1) * ( 1/( - 
     &    taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 + 1/( - 
     &    taucg + qsqhat)*msq*taucg**(-2)*taugs - 1/( - taucg + qsqhat)
     &    *taucs*taucg**(-1) - 1/( - taucg + qsqhat)*taucg**(-1)*taugs
     &     )
      DD2 = DD2 + lnrat( - qsqhat,msq)**2*gramdet**(-1) * (  - 1/( - 
     &    taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2 - 1/( - 
     &    taucg + qsqhat)*msq*taucg**(-2)*taugs + 1/( - taucg + qsqhat)
     &    *taucs*taucg**(-1) + 1/( - taucg + qsqhat)*taucg**(-1)*taugs
     &     )
      DD2 = DD2 + lnrat(musq,msq)*epinv*gramdet**(-1) * ( 1.D0/2.D0*msq
     &    *taucs**(-1)*taucg**(-2)*taugs - 1.D0/2.D0*taucg**(-1) )
      DD2 = DD2 + lnrat(musq,msq)*lnrat( - taucs,msq)*gramdet**(-1)
     &  * (  - msq*taucs**(-1)*taucg**(-2)*taugs + taucg**(-1) )
      DD2 = DD2 + lnrat(musq,msq)*lnrat( - taucg,msq)*gramdet**(-1)
     &  * (  - 1/( - taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*
     &    taugs**2 - 1/( - taucg + qsqhat)*msq*taucg**(-2)*taugs + 1/(
     &     - taucg + qsqhat)*taucs*taucg**(-1) + 1/( - taucg + qsqhat)*
     &    taucg**(-1)*taugs )
      DD2 = DD2 + lnrat(musq,msq)*lnrat( - qsqhat,msq)*gramdet**(-1)
     &  * ( 1/( - taucg + qsqhat)*msq*taucs**(-1)*taucg**(-2)*taugs**2
     &     + 1/( - taucg + qsqhat)*msq*taucg**(-2)*taugs - 1/( - taucg
     &     + qsqhat)*taucs*taucg**(-1) - 1/( - taucg + qsqhat)*
     &    taucg**(-1)*taugs )
      DD2 = DD2 + lnrat(musq,msq)**2*gramdet**(-1) * ( 1.D0/4.D0*msq*
     &    taucs**(-1)*taucg**(-2)*taugs - 1.D0/4.D0*taucg**(-1) )
      DD2 = DD2 + Lsm2x2m*gramdet**(-1) * ( msq*taucs**(-1)*taucg**(-2)
     &    *taugs - 1.D0/2.D0*taucg**(-1) )
      DD2 = DD2 + C0fb2m(tcg,msq)*gramdet**(-1) * ( 1.D0/2.D0 )
      DD2 = DD2 + C0fa2m(tcs,qsq,msq)*gramdet**(-1) * (  - 1.D0/2.D0 - 
     &    1.D0/2.D0*taucg**(-1)*taugs )
      DD2 = DD2 + C0fa2m(tcg,qsq,msq)*gramdet**(-1) * ( msq*taucs**(-1)
     &    *taucg**(-2)*taugs**2 + msq*taucg**(-2)*taugs - taucs*
     &    taucg**(-1) - taucg**(-1)*taugs )

      DD3=  + lnrat( - taucg,msq)*epinv*gramdet**(-1) * (  - 1/( - 
     &    taucg + qsqhat) + 1/( - taucg + qsqhat)*msq*taucs**(-1)*
     &    taucg**(-1)*taugs )
      DD3 = DD3 + lnrat( - taucg,msq)**2*gramdet**(-1) * ( 1/( - taucg
     &     + qsqhat) - 1/( - taucg + qsqhat)*msq*taucs**(-1)*
     &    taucg**(-1)*taugs )
      DD3 = DD3 + lnrat( - qsqhat,msq)*epinv*gramdet**(-1) * ( 1/( - 
     &    taucg + qsqhat) - 1/( - taucg + qsqhat)*msq*taucs**(-1)*
     &    taucg**(-1)*taugs )
      DD3 = DD3 + lnrat( - qsqhat,msq)**2*gramdet**(-1) * (  - 1/( - 
     &    taucg + qsqhat) + 1/( - taucg + qsqhat)*msq*taucs**(-1)*
     &    taucg**(-1)*taugs )
      DD3 = DD3 + lnrat(musq,msq)*lnrat( - taucg,msq)*gramdet**(-1)
     &  * (  - 1/( - taucg + qsqhat) + 1/( - taucg + qsqhat)*msq*
     &    taucs**(-1)*taucg**(-1)*taugs )
      DD3 = DD3 + lnrat(musq,msq)*lnrat( - qsqhat,msq)*gramdet**(-1)
     &  * ( 1/( - taucg + qsqhat) - 1/( - taucg + qsqhat)*msq*
     &    taucs**(-1)*taucg**(-1)*taugs )
      DD3 = DD3 + Lsm2x2m*gramdet**(-1) * (  - msq*taucs**(-1)*
     &    taucg**(-1) + 1.D0/2.D0*taugs**(-1) )
      DD3 = DD3 + C0fb2m(tcg,msq)*gramdet**(-1) * (  - 1.D0/2.D0*taucg*
     &    taugs**(-1) )
      DD3 = DD3 + C0fa2m(tcs,qsq,msq)*gramdet**(-1) * ( 1.D0/2.D0 + 1.D0
     &    /2.D0*taucg*taugs**(-1) )
      DD3 = DD3 + C0fa2m(tcg,qsq,msq)*gramdet**(-1) * ( 1.D0 - msq*
     &    taucs**(-1)*taucg**(-1)*taugs )

      return
      end
