      subroutine Cfort2(taucg,taucs,taugs,msq,
     & C00s,C00f,C11,C12,C22,C1,C2,C0)
      implicit none
C     C0(cg,s,msq,0,0)
      include 'epinv.f'
      include 'scale.f'
      double precision taucg,taucs,taugs,msq
      double precision qsq,qsqhat,tcg,qsqtcg
      double complex C00s,C00f,C11,C12,C22,C1,C2,C0,lnrat,C0fa2m
      qsq=taucg+taucs+taugs+msq
      qsqhat=taucg+taucs+taugs
      tcg=taucg+msq
      qsqtcg=qsqhat-taucg
      C0=  + lnrat( - taucg,msq)*epinv * ( qsqtcg**(-1) )
      C0 = C0 + lnrat( - taucg,msq)**2 * (  - qsqtcg**(-1) )
      C0 = C0 + lnrat( - qsqhat,msq)*epinv * (  - qsqtcg**(-1) )
      C0 = C0 + lnrat( - qsqhat,msq)**2 * ( qsqtcg**(-1) )
      C0 = C0 + lnrat(musq,msq)*lnrat( - taucg,msq) * ( qsqtcg**(-1) )
      C0 = C0 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * (  - 
     &    qsqtcg**(-1) )
      C0 = C0 - C0fa2m(tcg,qsq,msq)

      C1=  + lnrat( - taucg,msq)*epinv * (  - qsqtcg**(-1) )
      C1 = C1 + lnrat( - taucg,msq)**2 * ( qsqtcg**(-1) )
      C1 = C1 + lnrat( - taucg,msq) * (  - taucg*tcg**(-1)*qsqtcg**(-1)
     &     )
      C1 = C1 + lnrat( - qsqhat,msq)*epinv * ( qsqtcg**(-1) )
      C1 = C1 + lnrat( - qsqhat,msq)**2 * (  - qsqtcg**(-1) )
      C1 = C1 + lnrat( - qsqhat,msq) * ( taucg*qsq**(-1)*qsqtcg**(-1)
     &     + qsq**(-1) )
      C1 = C1 + lnrat(musq,msq)*lnrat( - taucg,msq) * (  - qsqtcg**(-1)
     &     )
      C1 = C1 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * ( qsqtcg**(-1) )
      C1 = C1 + C0fa2m(tcg,qsq,msq)

      C2=  + epinv * ( qsqtcg**(-1) )
      C2 = C2 + lnrat( - taucg,msq)*epinv * ( taucg*qsqtcg**(-2) )
      C2 = C2 + lnrat( - taucg,msq)**2 * (  - taucg*qsqtcg**(-2) )
      C2 = C2 + lnrat( - taucg,msq) * ( 2.D0*msq*taucg*tcg**(-1)*
     &    qsqtcg**(-2) + 2.D0*taucg**2*tcg**(-1)*qsqtcg**(-2) )
      C2 = C2 + lnrat( - qsqhat,msq)*epinv * (  - taucg*qsqtcg**(-2) )
      C2 = C2 + lnrat( - qsqhat,msq)**2 * ( taucg*qsqtcg**(-2) )
      C2 = C2 + lnrat( - qsqhat,msq) * (  - 2.D0*msq*taucg*qsq**(-1)*
     &    qsqtcg**(-2) - 2.D0*msq*qsq**(-1)*qsqtcg**(-1) - 3.D0*taucg*
     &    qsq**(-1)*qsqtcg**(-1) - 2.D0*taucg**2*qsq**(-1)*qsqtcg**(-2)
     &     - qsq**(-1) )
      C2 = C2 + lnrat(musq,msq)*lnrat( - taucg,msq) * ( taucg*
     &    qsqtcg**(-2) )
      C2 = C2 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * (  - taucg*
     &    qsqtcg**(-2) )
      C2 = C2 + lnrat(musq,msq) * ( qsqtcg**(-1) )
      C2 = C2 + 2.D0*qsqtcg**(-1) - C0fa2m(tcg,qsq,msq)*taucg*
     &    qsqtcg**(-1)

      C11=  + lnrat( - taucg,msq)*epinv * ( qsqtcg**(-1) )
      C11 = C11 + lnrat( - taucg,msq)**2 * (  - qsqtcg**(-1) )
      C11 = C11 + lnrat( - taucg,msq) * ( 1.D0/2.D0*msq*taucg*tcg**(-2)
     &    *qsqtcg**(-1) + 3.D0/2.D0*taucg*tcg**(-1)*qsqtcg**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq)*epinv * (  - qsqtcg**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq)**2 * ( qsqtcg**(-1) )
      C11 = C11 + lnrat( - qsqhat,msq) * (  - 1.D0/2.D0*msq*qsq**(-1)*
     &    qsqtcg**(-1) + 1.D0/2.D0*msq**2*qsq**(-2)*qsqtcg**(-1) - 3.D0/
     &    2.D0*taucg*qsq**(-1)*qsqtcg**(-1) - 3.D0/2.D0*qsq**(-1) )
      C11 = C11 + lnrat(musq,msq)*lnrat( - taucg,msq) * ( qsqtcg**(-1)
     &     )
      C11 = C11 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * (  - 
     &    qsqtcg**(-1) )
      C11 = C11 - 1.D0/2.D0*msq*tcg**(-1)*qsqtcg**(-1) + 1.D0/2.D0*msq*
     &    qsq**(-1)*qsqtcg**(-1) - C0fa2m(tcg,qsq,msq)

      C12=  + epinv * (  - qsqtcg**(-1) )
      C12 = C12 + lnrat( - taucg,msq)*epinv * (  - taucg*qsqtcg**(-2) )
      C12 = C12 + lnrat( - taucg,msq)**2 * ( taucg*qsqtcg**(-2) )
      C12 = C12 + lnrat( - taucg,msq) * (  - msq*taucg*tcg**(-1)*
     &    qsqtcg**(-2) - msq*taucg**2*tcg**(-2)*qsqtcg**(-2) - msq**2*
     &    taucg*tcg**(-2)*qsqtcg**(-2) - 5.D0/2.D0*taucg**2*tcg**(-1)*
     &    qsqtcg**(-2) )
      C12 = C12 + lnrat( - qsqhat,msq)*epinv * ( taucg*qsqtcg**(-2) )
      C12 = C12 + lnrat( - qsqhat,msq)**2 * (  - taucg*qsqtcg**(-2) )
      C12 = C12 + lnrat( - qsqhat,msq) * ( 2.D0*msq*taucg*qsq**(-1)*
     &    qsqtcg**(-2) + 3.D0/2.D0*msq*qsq**(-1)*qsqtcg**(-1) - msq**2*
     &    taucg*qsq**(-2)*qsqtcg**(-2) - 1.D0/2.D0*msq**2*qsq**(-2)*
     &    qsqtcg**(-1) + msq**2*qsq**(-1)*qsqtcg**(-2) - msq**3*
     &    qsq**(-2)*qsqtcg**(-2) + 7.D0/2.D0*taucg*qsq**(-1)*
     &    qsqtcg**(-1) + 5.D0/2.D0*taucg**2*qsq**(-1)*qsqtcg**(-2) + 
     &    qsq**(-1) )
      C12 = C12 + lnrat(musq,msq)*lnrat( - taucg,msq) * (  - taucg*
     &    qsqtcg**(-2) )
      C12 = C12 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * ( taucg*
     &    qsqtcg**(-2) )
      C12 = C12 + lnrat(musq,msq) * (  - qsqtcg**(-1) )
      C12 = C12 + msq*taucg*tcg**(-1)*qsqtcg**(-2) - msq*taucg*
     &    qsq**(-1)*qsqtcg**(-2) - 1.D0/2.D0*msq*qsq**(-1)*qsqtcg**(-1)
     &     + msq**2*tcg**(-1)*qsqtcg**(-2) - msq**2*qsq**(-1)*
     &    qsqtcg**(-2) - 5.D0/2.D0*qsqtcg**(-1) + C0fa2m(tcg,qsq,msq)*
     &    taucg*qsqtcg**(-1)

      C22=  + epinv * ( taucg*qsqtcg**(-2) - 1.D0/2.D0*qsqtcg**(-1) )
      C22 = C22 + lnrat( - taucg,msq)*epinv * ( taucg**2*qsqtcg**(-3) )
      C22 = C22 + lnrat( - taucg,msq)**2 * (  - taucg**2*qsqtcg**(-3) )
      C22 = C22 + lnrat( - taucg,msq) * ( 3.D0*msq*taucg**2*tcg**(-1)*
     &    qsqtcg**(-3) + 3.D0*taucg**3*tcg**(-1)*qsqtcg**(-3) )
      C22 = C22 + lnrat( - qsqhat,msq)*epinv * (  - taucg**2*
     &    qsqtcg**(-3) )
      C22 = C22 + lnrat( - qsqhat,msq)**2 * ( taucg**2*qsqtcg**(-3) )
      C22 = C22 + lnrat( - qsqhat,msq) * (  - 2.D0*msq*taucg*qsq**(-1)*
     &    qsqtcg**(-2) - 3.D0*msq*taucg**2*qsq**(-1)*qsqtcg**(-3) + 1.D0
     &    /2.D0*msq*qsq**(-1)*qsqtcg**(-1) - msq**2*taucg*qsq**(-2)*
     &    qsqtcg**(-2) - 1.D0/2.D0*msq**2*qsq**(-2)*qsqtcg**(-1) + 
     &    msq**2*qsq**(-1)*qsqtcg**(-2) - msq**3*qsq**(-2)*qsqtcg**(-2)
     &     - 1.D0/2.D0*taucg*qsq**(-1)*qsqtcg**(-1) - 4.D0*taucg**2*
     &    qsq**(-1)*qsqtcg**(-2) - 3.D0*taucg**3*qsq**(-1)*qsqtcg**(-3)
     &     + 1.D0/2.D0*qsq**(-1) )
      C22 = C22 + lnrat(musq,msq)*lnrat( - taucg,msq) * ( taucg**2*
     &    qsqtcg**(-3) )
      C22 = C22 + lnrat(musq,msq)*lnrat( - qsqhat,msq) * (  - taucg**2*
     &    qsqtcg**(-3) )
      C22 = C22 + lnrat(musq,msq) * ( taucg*qsqtcg**(-2) - 1.D0/2.D0*
     &    qsqtcg**(-1) )
      C22 = C22 - msq*taucg*qsq**(-1)*qsqtcg**(-2) - 1.D0/2.D0*msq*
     &    qsq**(-1)*qsqtcg**(-1) + msq*qsqtcg**(-2) - msq**2*qsq**(-1)*
     &    qsqtcg**(-2) + 3.D0*taucg*qsqtcg**(-2) - qsqtcg**(-1) - 
     &    C0fa2m(tcg,qsq,msq)*taucg**2*qsqtcg**(-2)

      C00s=  + 1.D0/4.D0

      C00f=  + lnrat( - taucg,msq) * ( 1.D0/4.D0*taucg**2*tcg**(-1)*
     &    qsqtcg**(-1) )
      C00f = C00f + lnrat( - qsqhat,msq) * (  - 1.D0/4.D0*taucs*
     &    qsq**(-1) - 1.D0/2.D0*taucg*qsq**(-1) - 1.D0/4.D0*taucg**2*
     &    qsq**(-1)*qsqtcg**(-1) - 1.D0/4.D0*taugs*qsq**(-1) )
      C00f = C00f + lnrat(musq,msq) * ( 1.D0/4.D0 )
      C00f = C00f + 3.D0/4.D0

      return
      end
