      double complex function virt_pp(mQ,ig,is,ie,in,ic,p)
      implicit none
************************************************************************
*     Author: F. Tramontano                                            *
*     Dicember, 2008                                                   *
************************************************************************
c---- One-loop ++ helicity amplitude for W+c production
c---- hg=+1  with s(mu) as gauge vector
c---- hQ=+1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5) 
c----
      include 'constants.f'
      include 'epinv.f'
      include 'b0.f'
      include 'scheme.f'
      include 'scale.f'
      include 'zprods_decl.f'
      integer is,ig,ic,ie,in,nu,jpart
      double precision p(mxpart,4),q(mxpart,4)
      double precision dot,taucg,taugs,taucs,msq,mQ
      double precision qsq,tcg,tcs,qsqhat,qcg,qcs,cgs1
      double precision ddilog,epin,epin2,xlog
      double complex spp,lopp
      double complex app1,app2,app3,app4
      double complex ffcg1,ffcg3,ffcs2,fun1,fun2
      double complex fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      double complex lnrat,C0fa2m,C0fb2m,I3me
c      double precision mass2,width2,mass3,width3
c      common/breit/n2,n3,mass2,width2,mass3,width3
      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv**2+epinv*xlog+half*xlog**2
      taugs=+2d0*dot(p,is,ig)
      taucs=+2d0*dot(p,ic,is)
      taucg=+2d0*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      ffcs2=fun2(tcs,qsq,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart.eq.ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      app1=za(ig,ie)*za(is,ic)*zb(is,in)/za(ig,is)/za(ig,ic)**2
      app2=za(ig,ie)*zb(ig,in)/za(ig,is)/za(ig,ic)/za(is,ic)/
     .     zb(is,ic)
      app3=za(ig,ie)*zb(is,in)/za(ig,is)/za(ig,ic)**2/zb(is,ic)
      app4=za(is,ie)*zb(is,in)/za(ig,is)/za(ig,ic)/za(is,ic)/
     .     zb(is,ic)
      spp=  + taucs*taucg*taugs**(-1)*mQ*xn**(-1)*app3*fL6m2 + 1.D0/2.D0
     &    *taucs*mQ*xn**(-1)*app3*fL6m1 + 1.D0/2.D0*taucs*mQ*xn*app3*
     &    fL6m3 + 2.D0*taucg**(-2)*taugs**2*qsq*mQ*cf*app4*ffcg3 + 
     &    taucg**(-1)*taugs*tcs*mQ*xn**(-1)*app2*ffcs2 + 2.D0*
     &    taucg**(-1)*taugs*qsq*mQ*cf*app3*ffcg3 - taucg**(-1)*taugs*
     &    qsq*mQ*xn**(-1)*app4*ffcg3 + taucg**(-1)*taugs**2*mQ*cf*app4*
     &    ffcg1 - 3.D0*taucg**(-1)*taugs**2*mQ*cf*app4*ffcg3 - 
     &    taucg**(-1)*mQ*qsqhat*cf*cgs1*app2 - taucg**(-1)*mQ*cf*qcs**2
     &    *app4 + taucg**(-1)*mQ*cf*qcg*cgs1*app2 - 1.D0/2.D0*
     &    taucg**(-1)*mQ**3*xn**(-1)*qcs*app4*fL6m2 - 1.D0/2.D0*
     &    taucg**(-1)*mQ**3*xn**(-1)*qcg*app2*fL6m2 + taucg*tcg**(-1)*
     &    qsq*mQ*cf*app4 - taucg*tcg**(-1)*mQ*cf*qcg*app3 - taugs**(-1)
     &    *mQ*xn**(-1)*qcg*cgs1*app2*fL6m2 + taugs*mQ*cf*app3*ffcg3 + 
     &    taugs*mQ*cf*app3 - 2.D0*taugs*mQ*cf*app4*ffcg3 - 1.D0/2.D0*
     &    tcs*mQ*xn**(-1)*app3*fL6m2 + xlog*mQ*b0*cgs1*app2 - xlog*mQ*
     &    b0*app1
      spp = spp - 3.D0/2.D0*mQ*epin*cf*cgs1*app2 + 3.D0/2.D0*mQ*epin*cf
     &    *app1 - mQ*epin*b0*cgs1*app2 + mQ*epin*b0*app1 + 1.D0/2.D0*mQ
     &    *epin*xn**(-1)*cgs1*app2 - 1.D0/2.D0*mQ*epin*xn**(-1)*app1 - 
     &    1.D0/2.D0*mQ*epin*xn*cgs1*app2 + 1.D0/2.D0*mQ*epin*xn*app1 + 
     &    1.D0/2.D0*mQ*epin2*xn**(-1)*cgs1*app2 - 1.D0/2.D0*mQ*epin2*
     &    xn**(-1)*app1 - 3.D0/2.D0*mQ*epin2*xn*cgs1*app2 + 3.D0/2.D0*
     &    mQ*epin2*xn*app1 - 1.D0/2.D0*mQ*cf*cgs1*app2 + 5.D0/2.D0*mQ*
     &    cf*app1 + 1.D0/2.D0*mQ*xn**(-1)*cgs1*pisqo6*app2 + 1.D0/2.D0*
     &    mQ*xn**(-1)*cgs1*app2*fL6m1 + 1.D0/2.D0*mQ*xn**(-1)*cgs1*app2
     &    *fL6m2 + 3.D0/2.D0*mQ*xn**(-1)*cgs1*app2 + mQ*xn**(-1)*cgs1*
     &    app3*fL6m2 + 1.D0/2.D0*mQ*xn**(-1)*cgs1*app3 - 1.D0/2.D0*mQ*
     &    xn**(-1)*pisqo6*app1 - 1.D0/2.D0*mQ*xn*cgs1*pisqo6*app2 - 1.D0
     &    /2.D0*mQ*xn*cgs1*app2*fL6m3 - mQ*xn*cgs1*app3*fL6m3 + 1.D0/2.D
     &    0*mQ*xn*pisqo6*app1 + ddilog(msq**(-1)*tcs)*mQ*xn**(-1)*cgs1*
     &    app2 - ddilog(msq**(-1)*tcs)*mQ*xn**(-1)*app1 - ddilog(
     &    msq**(-1)*tcg)*mQ*xn*cgs1*app2
      spp = spp + ddilog(msq**(-1)*tcg)*mQ*xn*app1 - I3me(msq,taugs,qsq
     &    )*taucs**(-1)*taugs*mQ**3*xn**(-1)*cgs1*app2 + 4.D0*I3me(msq,
     &    taugs,qsq)*taucg**(-1)*taugs*mQ**3*cf*cgs1*app3 + 2.D0*I3me(
     &    msq,taugs,qsq)*taucg**(-1)*taugs*mQ**3*xn**(-1)*cgs1*app3 + 
     &    I3me(msq,taugs,qsq)*taucg**(-1)*taugs*mQ**3*xn*cgs1*app2 - 
     &    lnrat( - taucs,msq)*taucs**2*taucg**(-1)*mQ*xn**(-1)*app2 - 
     &    lnrat( - taucs,msq)*taucg*tcs**(-1)*mQ*xn**(-1)*cgs1*app3 + 
     &    lnrat( - taucs,msq)*taugs*tcs**(-1)*mQ*xn**(-1)*cgs1*app2 + 
     &    lnrat( - taucs,msq)*tcs**(-1)*qsq*mQ*xn**(-1)*cgs1*app4 - 
     &    lnrat( - taucs,msq)*mQ*epin*xn**(-1)*cgs1*app2 + lnrat( - 
     &    taucs,msq)*mQ*epin*xn**(-1)*app1 + lnrat( - taucs,msq)**2*mQ*
     &    xn**(-1)*cgs1*app2 - lnrat( - taucs,msq)**2*mQ*xn**(-1)*app1
     &     - 2.D0*lnrat( - taucg,msq)*taucg**(-1)*tcg**(-1)*qsq*mQ*cf*
     &    qcs**2*app4 - 2.D0*lnrat( - taucg,msq)*taucg**(-1)*mQ*cf*
     &    qcs**2*app4 - lnrat( - taucg,msq)*taucg*tcg**(-1)*qsq*mQ*cf*
     &    app2
      spp = spp - 2.D0*lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*cf*qcs*
     & app3 + lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*cf*cgs1*app2 - 
     &    lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*xn**(-1)*cgs1*app3 + 
     &    lnrat( - taucg,msq)*taucg**2*tcg**(-2)*qsq*mQ*cf*app3 - 
     &    lnrat( - taucg,msq)*taucg**2*tcg**(-2)*qsq*mQ*cf*app4 - 2.D0*
     &    lnrat( - taucg,msq)*tcg**(-1)*qsq*mQ*cf*qcs*app3 + 4.D0*
     &    lnrat( - taucg,msq)*tcg**(-1)*qsq*mQ*cf*qcs*app4 + lnrat( - 
     &    taucg,msq)*tcg**(-1)*qsq*mQ*xn**(-1)*qcs*app4 + 2.D0*lnrat(
     &     - taucg,msq)*tcg**(-1)*mQ*cf*qcs**2*app4 + lnrat( - taucg,
     &    msq)*mQ*epin*xn*cgs1*app2 - lnrat( - taucg,msq)*mQ*epin*xn*
     &    app1 + lnrat( - taucg,msq)*mQ*cf*qcs*app2 + 2.D0*lnrat( - 
     &    taucg,msq)*mQ*cf*cgs1*app2 - lnrat( - taucg,msq)*mQ*xn**(-1)*
     &    cgs1*app2 - lnrat( - taucg,msq)**2*mQ*xn*cgs1*app2 + lnrat(
     &     - taucg,msq)**2*mQ*xn*app1 + lnrat( - taugs,msq)*mQ*epin*xn*
     &    cgs1*app2 - lnrat( - taugs,msq)*mQ*epin*xn*app1 - 1.D0/2.D0*
     &    lnrat(
     &  - taugs,msq)**2*mQ*xn*cgs1*app2 + 1.D0/2.D0*lnrat( - taugs,msq)
     &    **2*mQ*xn*app1 + lnrat( - qsqhat,msq)*taucs*taucg**(-1)*tcs*
     &    qsq**(-1)*mQ*qsqhat*xn**(-1)*app2 + lnrat( - qsqhat,msq)*
     &    taucg**(-2)*qsq**(-1)*mQ**3*cf*qcs**2*qcg*app4 + lnrat( - 
     &    qsqhat,msq)*taucg**(-2)*mQ*cf*qcs**2*qcg*app4 + lnrat( - 
     &    qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ*qsqhat*cf*qcs*qcg*app3
     &     + lnrat( - qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ**3*cf*qcs**2
     &    *app4 + lnrat( - qsqhat,msq)*taucg**(-1)*mQ*qsqhat*cf*qcs*
     &    app3 - 2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*mQ*qsqhat*cf*qcs
     &    *app4 - lnrat( - qsqhat,msq)*taucg**(-1)*mQ*qsqhat*xn**(-1)*
     &    qcs*app4 + lnrat( - qsqhat,msq)*taucg**(-1)*mQ*cf*qcs**2*app4
     &     + lnrat( - qsqhat,msq)*qsq**(-1)*mQ*qsqhat*cf*qcs*app3 + 
     &    lnrat( - qsqhat,msq)*qsq**(-1)*mQ*qsqhat*xn**(-1)*cgs1*app2
     &     + lnrat( - qsqhat,msq)*mQ*qsqhat*cf*app3 - C0fb2m(tcg,msq)*
     &    taucs**(-1)*taugs*mQ**3*xn**(-1)*cgs1*app3 - C0fb2m(tcg,msq)*
     &    taucs**(-1)*mQ**3*xn**(-1)*qcs*cgs1*app3
      spp = spp - C0fb2m(tcg,msq)*taucs**(-1)*mQ**3*xn**(-1)*qcs*cgs1*
     & app4 + 2.D0*C0fb2m(tcg,msq)*mQ**3*xn**(-1)*cgs1*app2 + 2.D0*
     &    C0fa2m(tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*taugs*mQ**3*
     &    xn**(-1)*qcs*cgs1*app3 + C0fa2m(tcs,qsq,msq)*taucs**(-1)*
     &    taucg**(-1)*mQ**3*xn**(-1)*qcs**2*cgs1*app4 + C0fa2m(tcs,qsq,
     &    msq)*taucs**(-1)*mQ**3*xn**(-1)*qcs*cgs1*app3 - C0fa2m(tcs,
     &    qsq,msq)*taucg**(-1)*taugs*mQ**3*xn**(-1)*cgs1*app2 - 2.D0*
     &    C0fa2m(tcs,qsq,msq)*mQ**3*xn**(-1)*cgs1*app2

      
      lopp=-mQ*za(ig,ie)/za(ig,is)/za(ig,ic)**2*
     . (za(is,ic)*zb(is,in)+za(ig,ic)*zb(ig,in))
      
c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      spp=spp + lopp * xn/6d0
      
c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - cf * alphas(MSbar) / (2*pi))
      spp=spp - lopp * cf
      
      virt_pp=spp
      return
      end
