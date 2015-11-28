      double precision function virtsqwcg(is,ig,ie,in,ic,p)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2004.                                                     *
************************************************************************
c---- One-loop matrix element for W+c production, including c mass
c---- averaged over initial colours and spins
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + cbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5) 
c----
      include 'constants.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'scale.f'
      include 'masses.f'
      include 'nflav.f'
      include 'b0.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer is,ig,ic,ie,in,nu,jpart,i,j,js
      integer t1,t2,t3,t4,t6,t7,t8,t9
      integer t1b,t2b,t3b,t4b,t6b,t8b
      integer b1,b2,b3,s1,s2
      parameter(b1=1,b2=2,b3=3,s1=1,s2=2)
      parameter(t1=1, t2=2,  t3=3,  t4=4,  t6=5,  t8=6)
      parameter(t1b=7,t2b=8,t3b=9,t4b=10,t6b=11,t8b=12)
      parameter(t7=13,t9=14)
      double precision p(mxpart,4),q(mxpart,4),fac,prop
      double precision dot,taucg,taugs,taucs,msq
      double precision epin2,epin,poles_qg,poles_gq,w1cjet,temp
      double complex lnrat,lomm,lopp,lopm,lomp
      double complex smm,spp,spm,smp
      double complex a1mm,a1pp,a1pm,a1mp
      double complex a2mm,a2pp,a2pm,a2mp
      double complex a3mm,a3pp,a3pm,a3mp
      double complex a4mm,a4pp,a4pm,a4mp
      double complex D0(3),D1(3),D2(3),D3(3),
     & C0(14),C1(14),C2(14),C11(14),C12(14),C22(14),C00s(14),C00f(14),
     & BB0(2),BB1(2)
c      double precision mass2,width2,mass3,width3
c      common/breit/n2,n3,mass2,width2,mass3,width3
      scheme='dred'
      call dotem(5,p,s)
      temp=-w1cjet(is,ic,ie,in,ig)
      taugs=+2d0*dot(p,is,ig)
      taucs=+2d0*dot(p,ic,is)
      taucg=+2d0*dot(p,ic,ig)
      msq=mc**2
      call Dforb1(taucg,taucs,taugs,msq,D3(b1),D2(b1),D1(b1),D0(b1))
      call Dforb2(taucg,taucs,taugs,msq,D3(b2),D2(b2),D1(b2),D0(b2))
      call Dforb3(taucg,taucs,taugs,msq,D3(b3),D2(b3),D1(b3),D0(b3))
      call Cfort1(taucg,taucs,taugs,msq,C00s(t1),C00f(t1),
     . C11(t1),C12(t1),C22(t1),C1(t1),C2(t1),C0(t1))
      call Cfort2(taucg,taucs,taugs,msq,C00s(t2),C00f(t2),
     . C11(t2),C12(t2),C22(t2),C1(t2),C2(t2),C0(t2))
      call Cfort3(taucg,taucs,taugs,msq,C00s(t3),C00f(t3),
     . C11(t3),C12(t3),C22(t3),C1(t3),C2(t3),C0(t3))
      call Cfort4(taucg,taucs,taugs,msq,C00s(t4),C00f(t4),
     . C11(t4),C12(t4),C22(t4),C1(t4),C2(t4),C0(t4))
      call Cfort6(taucg,taucs,taugs,msq,C00s(t6),C00f(t6),
     . C11(t6),C12(t6),C22(t6),C1(t6),C2(t6),C0(t6))
      call Cfort8(taucg,taucs,taugs,msq,C00s(t8),C00f(t8),
     . C11(t8),C12(t8),C22(t8),C1(t8),C2(t8),C0(t8))
      call Cfort7(taucg,taucs,taugs,msq,C00s(t7),C00f(t7),
     . C11(t7),C12(t7),C22(t7),C1(t7),C2(t7),C0(t7))
      call Cfort9(taucg,taucs,taugs,msq,C00s(t9),C00f(t9),
     . C11(t9),C12(t9),C22(t9),C1(t9),C2(t9),C0(t9))
      do j=1,6
      js=j+6
      C11(js) =C22(j)+2d0*C2(j)+C0(j)
      C22(js) =C11(j)+2d0*C1(j)+C0(j)
      C12(js) =C12(j)+C1(j)+C2(j)+C0(j)
      C1(js)  =-C2(j)-C0(j)
      C2(js)  =-C1(j)-C0(j)
      C0(js)  =C0(j)
      C00s(js)=C00s(j)
      C00f(js)=C00f(j)
c      write(6,*) js,C11(js),C22(js),C12(js),C1(js),C2(js),C0(js),
c     & C00s(js),C00f(js)
      enddo
c      call Cfort1b(taucg,taucs,taugs,msq,C00s(t1b),C00f(t1b),
c     . C11(t1b),C12(t1b),C22(t1b),C1(t1b),C2(t1b),C0(t1b))
c      call Cfort2b(taucg,taucs,taugs,msq,C00s(t2b),C00f(t2b),
c     . C11(t2b),C12(t2b),C22(t2b),C1(t2b),C2(t2b),C0(t2b))
c      call Cfort3b(taucg,taucs,taugs,msq,C00s(t3b),C00f(t3b),
c     . C11(t3b),C12(t3b),C22(t3b),C1(t3b),C2(t3b),C0(t3b))
c      call Cfort4b(taucg,taucs,taugs,msq,C00s(t4b),C00f(t4b),
c     . C11(t4b),C12(t4b),C22(t4b),C1(t4b),C2(t4b),C0(t4b))
c      call Cfort6b(taucg,taucs,taugs,msq,C00s(t6b),C00f(t6b),
c     . C11(t6b),C12(t6b),C22(t6b),C1(t6b),C2(t6b),C0(t6b))
c      call Cfort8b(taucg,taucs,taugs,msq,C00s(t8b),C00f(t8b),
c     . C11(t8b),C12(t8b),C22(t8b),C1(t8b),C2(t8b),C0(t8b))
c      do j=1,6
c      js=j+6
c      write(6,*) js,C11(js),C22(js),C12(js),C1(js),C2(js),C0(js),
c     & C00s(js),C00f(js)
c      enddo
c      pause
      call Bfors1(taugs,msq,BB1(s1),BB0(s1))
      call Bfors2(taucg,msq,BB1(s2),BB0(s2))
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart.eq.ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      a1mm=za(ig,ie)*zb(ig,in)
     . /za(ig,is)/za(is,ic)/zb(ig,is)/zb(ig,ic)**2
      a2mm=za(is,ie)*zb(is,in)
     . /za(ig,is)/za(is,ic)/zb(ig,is)/zb(ig,ic)**2
      a3mm=za(ie,ic)*zb(is,in)/za(is,ic)/zb(ig,is)/zb(ig,ic)
      a4mm=za(ie,ic)*zb(in,ic)
     . /za(ig,is)/za(is,ic)/zb(ig,is)/zb(ig,ic)**2
      a1pp=za(ig,ie)*zb(ig,in)/za(ig,is)/za(ig,ic)/za(is,ic)/zb(is,ic)
      a2pp=za(ig,ie)*zb(ig,ic)*zb(is,in)/za(ig,is)/za(ig,ic)/zb(is,ic)
      a3pp=za(is,ie)*zb(is,in)/za(ig,is)/za(ig,ic)/za(is,ic)/zb(is,ic)
      a4pp=za(ie,ic)*zb(in,ic)/za(ig,is)/za(ig,ic)/za(is,ic)/zb(is,ic)
      a1pm=za(ig,ie)*zb(ig,in)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      a2pm=za(ig,ie)*zb(is,in)/za(is,ic)/zb(ig,is)
      a3pm=za(is,ie)*zb(is,in)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      a4pm=za(ie,ic)*zb(in,ic)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      a1mp =za(ig,ie)*zb(ig,in)/za(ig,is)**2/zb(is,ic)
      a2mp=za(is,ie)*zb(is,in)/za(ig,is)**2/zb(is,ic)
      a3mp=za(ie,ic)*zb(ig,ic)*zb(is,in)/za(ig,is)/zb(is,ic)
      a4mp=za(ie,ic)*zb(in,ic)/za(ig,is)**2/zb(is,ic)
      spp=  + mc*xn**(-1) * ( 3.D0/4.D0*epinv*taucs*taucg**(-1)*a2pp + 
     &    3.D0/4.D0*epinv*taucs*a1pp + 3.D0/2.D0*epinv*taucs*a3pp + 5.D0
     &    /2.D0*taucs*taucg**(-1)*a2pp + 5.D0/2.D0*taucs*a3pp + 2.D0*
     &    C00s(t1)*epinv*taucs*taucg**(-1)*a2pp - 2.D0*C00s(t1)*epinv*
     &    taucs*a1pp + 2.D0*C00s(t1)*taucs*taucg**(-1)*taugs*a1pp + 2.D0
     &    *C00s(t1)*taucs*taucg**(-1)*taugs*a3pp + 2.D0*C00s(t1)*taucs*
     &    taucg**(-1)*a2pp + 2.D0*C00s(t1)*taucs*a1pp + 2.D0*C00s(t1)*
     &    taucs*a3pp + 2.D0*C00s(t1)*taucs**2*taucg**(-1)*a1pp - 2.D0*
     &    C00s(t1b)*epinv*taucs*taucg**(-1)*a2pp - 2.D0*C00s(t1b)*epinv
     &    *taucs**2*taucg**(-1)*a1pp + 2.D0*C00s(t2)*epinv*taucs*
     &    taucg**(-1)*taugs*a1pp + 2.D0*C00s(t2)*epinv*taucs*
     &    taucg**(-1)*a2pp + 2.D0*C00s(t2)*epinv*taucs**2*taucg**(-1)*
     &    a1pp + 2.D0*C00s(t3)*epinv*taucs*taucg**(-1)*a2pp - 2.D0*
     &    C00s(t3)*epinv*taucs*a1pp - 2.D0*C00s(t3)*taucs*taucg**(-1)*
     &    a2pp + 2.D0*C00s(t3b)*epinv*taucs*taucg**(-1)*taugs*a1pp - 2.D
     &    0*C00s(t7)*epinv*taucs*taucg**(-1)*taugs*a1pp )
      spp = spp + mc*xn**(-1) * ( 2.D0*C00s(t8)*epinv*taucs*taucg**(-1)
     &    *a2pp + 2.D0*C00s(t8)*epinv*taucs**2*taucg**(-1)*a1pp - 2.D0*
     &    C00s(t8b)*epinv*taucs*taucg**(-1)*a2pp - 2.D0*C00s(t8b)*epinv
     &    *taucs**2*taucg**(-1)*a1pp - 2.D0*C00s(t9)*epinv*taucs*
     &    taucg**(-1)*taugs*a1pp + 2.D0*C00f(t1)*taucs*taucg**(-1)*a2pp
     &     - 2.D0*C00f(t1)*taucs*a1pp - 2.D0*C00f(t1b)*taucs*
     &    taucg**(-1)*a2pp - 2.D0*C00f(t1b)*taucs**2*taucg**(-1)*a1pp
     &     + 2.D0*C00f(t2)*taucs*taucg**(-1)*taugs*a1pp + 2.D0*C00f(t2)
     &    *taucs*taucg**(-1)*a2pp + 2.D0*C00f(t2)*taucs**2*taucg**(-1)*
     &    a1pp + 2.D0*C00f(t3)*taucs*taucg**(-1)*a2pp - 2.D0*C00f(t3)*
     &    taucs*a1pp + 2.D0*C00f(t3b)*taucs*taucg**(-1)*taugs*a1pp - 2.D
     &    0*C00f(t7)*taucs*taucg**(-1)*taugs*a1pp + 2.D0*C00f(t8)*taucs
     &    *taucg**(-1)*a2pp + 2.D0*C00f(t8)*taucs**2*taucg**(-1)*a1pp
     &     - 2.D0*C00f(t8b)*taucs*taucg**(-1)*a2pp - 2.D0*C00f(t8b)*
     &    taucs**2*taucg**(-1)*a1pp - 2.D0*C00f(t9)*taucs*taucg**(-1)*
     &    taugs*a1pp )
      spp = spp + mc*xn**(-1) * ( 3.D0/2.D0*lnrat(musq,msq)*taucs*
     &    taucg**(-1)*a2pp + 3.D0/2.D0*lnrat(musq,msq)*taucs*a3pp )
      spp = spp + mc*xn * (  - 3.D0/4.D0*epinv*taucs*taucg**(-1)*a2pp
     &     - 3.D0/4.D0*epinv*taucs*a1pp - 3.D0/2.D0*epinv*taucs*a3pp - 
     &    8.D0/3.D0*taucs*taucg**(-1)*a2pp + 1.D0/6.D0*taucs*a1pp - 5.D0
     &    /2.D0*taucs*a3pp - 2.D0*C00s(t1)*epinv*taucs*taucg**(-1)*a2pp
     &     + 2.D0*C00s(t1)*epinv*taucs*a1pp - 2.D0*C00s(t1)*taucs*
     &    taucg**(-1)*taugs*a1pp - 2.D0*C00s(t1)*taucs*taucg**(-1)*
     &    taugs*a3pp - 2.D0*C00s(t1)*taucs*taucg**(-1)*a2pp - 2.D0*
     &    C00s(t1)*taucs*a1pp - 2.D0*C00s(t1)*taucs*a3pp - 2.D0*C00s(t1
     &    )*taucs**2*taucg**(-1)*a1pp + 2.D0*C00s(t1b)*epinv*taucs*
     &    taucg**(-1)*taugs*a1pp + 2.D0*C00s(t1b)*epinv*taucs*
     &    taucg**(-1)*taugs*a3pp + 4.D0*C00s(t1b)*epinv*taucs*
     &    taucg**(-1)*a2pp + 4.D0*C00s(t1b)*epinv*taucs**2*taucg**(-1)*
     &    a1pp - 4.D0*C00s(t2b)*epinv*taucs*taucg**(-1)*taugs*a1pp - 2.D
     &    0*C00s(t2b)*epinv*taucs*taucg**(-1)*taugs*a3pp - 4.D0*C00s(
     &    t2b)*epinv*taucs*taucg**(-1)*a2pp - 4.D0*C00s(t2b)*epinv*
     &    taucs**2*taucg**(-1)*a1pp )
      spp = spp + mc*xn * ( 2.D0*C00s(t3)*epinv*taucs*taucg**(-1)*taugs
     &    *a1pp - 10.D0*C00s(t3)*epinv*taucs*taucg**(-1)*a2pp + 10.D0*
     &    C00s(t3)*epinv*taucs*a1pp + 4.D0*C00s(t3)*taucs*taucg**(-1)*
     &    a2pp - 6.D0*C00s(t3)*taucs*a1pp - 2.D0*C00f(t1)*taucs*
     &    taucg**(-1)*a2pp + 2.D0*C00f(t1)*taucs*a1pp + 2.D0*C00f(t1b)*
     &    taucs*taucg**(-1)*taugs*a1pp + 2.D0*C00f(t1b)*taucs*
     &    taucg**(-1)*taugs*a3pp + 4.D0*C00f(t1b)*taucs*taucg**(-1)*
     &    a2pp + 4.D0*C00f(t1b)*taucs**2*taucg**(-1)*a1pp - 4.D0*C00f(
     &    t2b)*taucs*taucg**(-1)*taugs*a1pp - 2.D0*C00f(t2b)*taucs*
     &    taucg**(-1)*taugs*a3pp - 4.D0*C00f(t2b)*taucs*taucg**(-1)*
     &    a2pp - 4.D0*C00f(t2b)*taucs**2*taucg**(-1)*a1pp + 2.D0*C00f(
     &    t3)*taucs*taucg**(-1)*taugs*a1pp - 10.D0*C00f(t3)*taucs*
     &    taucg**(-1)*a2pp + 10.D0*C00f(t3)*taucs*a1pp - 3.D0/2.D0*
     &    lnrat(musq,msq)*taucs*taucg**(-1)*a2pp - 3.D0/2.D0*lnrat(musq
     &    ,msq)*taucs*a3pp )
      spp = spp + mc * ( epinv*taucs*taucg**(-1)*b0*a2pp - epinv*taucs*
     &    b0*a1pp )
      spp = spp + mc**3*xn**(-1) * (  - 3.D0/4.D0*epinv*taucg**(-2)*
     &    taugs*a2pp - 3.D0/4.D0*epinv*taucg**(-1)*taugs*a1pp - 3.D0/2.D
     &    0*epinv*taucg**(-1)*taugs*a3pp - 5.D0/2.D0*taucg**(-2)*taugs*
     &    a2pp - 5.D0/2.D0*taucg**(-1)*taugs*a3pp - 2.D0*C00s(t1)*epinv
     &    *taucg**(-2)*taugs*a2pp + 2.D0*C00s(t1)*epinv*taucg**(-1)*
     &    taugs*a1pp - 2.D0*C00s(t1)*taucs*taucg**(-2)*taugs*a1pp - 2.D0
     &    *C00s(t1)*taucg**(-2)*taugs*a2pp - 2.D0*C00s(t1)*taucg**(-2)*
     &    taugs**2*a1pp - 2.D0*C00s(t1)*taucg**(-2)*taugs**2*a3pp - 2.D0
     &    *C00s(t1)*taucg**(-1)*taugs*a1pp - 2.D0*C00s(t1)*taucg**(-1)*
     &    taugs*a3pp + 2.D0*C00s(t1b)*epinv*taucs*taucg**(-2)*taugs*
     &    a1pp + 2.D0*C00s(t1b)*epinv*taucg**(-2)*taugs*a2pp - 2.D0*
     &    C00s(t2)*epinv*taucs*taucg**(-2)*taugs*a1pp - 2.D0*C00s(t2)*
     &    epinv*taucg**(-2)*taugs*a2pp - 2.D0*C00s(t2)*epinv*
     &    taucg**(-2)*taugs**2*a1pp - 2.D0*C00s(t3)*epinv*taucg**(-2)*
     &    taugs*a2pp + 2.D0*C00s(t3)*epinv*taucg**(-1)*taugs*a1pp + 2.D0
     &    *C00s(t3)*taucg**(-2)*taugs*a2pp )
      spp = spp + mc**3*xn**(-1) * (  - 2.D0*C00s(t3b)*epinv*
     &    taucg**(-2)*taugs**2*a1pp + 2.D0*C00s(t7)*epinv*taucg**(-2)*
     &    taugs**2*a1pp - 2.D0*C00s(t8)*epinv*taucs*taucg**(-2)*taugs*
     &    a1pp - 2.D0*C00s(t8)*epinv*taucg**(-2)*taugs*a2pp + 2.D0*
     &    C00s(t8b)*epinv*taucs*taucg**(-2)*taugs*a1pp + 2.D0*C00s(t8b)
     &    *epinv*taucg**(-2)*taugs*a2pp + 2.D0*C00s(t9)*epinv*
     &    taucg**(-2)*taugs**2*a1pp - 2.D0*C00f(t1)*taucg**(-2)*taugs*
     &    a2pp + 2.D0*C00f(t1)*taucg**(-1)*taugs*a1pp + 2.D0*C00f(t1b)*
     &    taucs*taucg**(-2)*taugs*a1pp + 2.D0*C00f(t1b)*taucg**(-2)*
     &    taugs*a2pp - 2.D0*C00f(t2)*taucs*taucg**(-2)*taugs*a1pp - 2.D0
     &    *C00f(t2)*taucg**(-2)*taugs*a2pp - 2.D0*C00f(t2)*taucg**(-2)*
     &    taugs**2*a1pp - 2.D0*C00f(t3)*taucg**(-2)*taugs*a2pp + 2.D0*
     &    C00f(t3)*taucg**(-1)*taugs*a1pp - 2.D0*C00f(t3b)*taucg**(-2)*
     &    taugs**2*a1pp + 2.D0*C00f(t7)*taucg**(-2)*taugs**2*a1pp - 2.D0
     &    *C00f(t8)*taucs*taucg**(-2)*taugs*a1pp - 2.D0*C00f(t8)*
     &    taucg**(-2)*taugs*a2pp )
      spp = spp + mc**3*xn**(-1) * ( 2.D0*C00f(t8b)*taucs*taucg**(-2)*
     &    taugs*a1pp + 2.D0*C00f(t8b)*taucg**(-2)*taugs*a2pp + 2.D0*
     &    C00f(t9)*taucg**(-2)*taugs**2*a1pp - 3.D0/2.D0*lnrat(musq,msq
     &    )*taucg**(-2)*taugs*a2pp - 3.D0/2.D0*lnrat(musq,msq)*
     &    taucg**(-1)*taugs*a3pp )
      spp = spp + mc**3*xn * ( 3.D0/4.D0*epinv*taucg**(-2)*taugs*a2pp
     &     + 3.D0/4.D0*epinv*taucg**(-1)*taugs*a1pp + 3.D0/2.D0*epinv*
     &    taucg**(-1)*taugs*a3pp + 8.D0/3.D0*taucg**(-2)*taugs*a2pp - 1.
     &    D0/6.D0*taucg**(-1)*taugs*a1pp + 5.D0/2.D0*taucg**(-1)*taugs*
     &    a3pp + 2.D0*C00s(t1)*epinv*taucg**(-2)*taugs*a2pp - 2.D0*
     &    C00s(t1)*epinv*taucg**(-1)*taugs*a1pp + 2.D0*C00s(t1)*taucs*
     &    taucg**(-2)*taugs*a1pp + 2.D0*C00s(t1)*taucg**(-2)*taugs*a2pp
     &     + 2.D0*C00s(t1)*taucg**(-2)*taugs**2*a1pp + 2.D0*C00s(t1)*
     &    taucg**(-2)*taugs**2*a3pp + 2.D0*C00s(t1)*taucg**(-1)*taugs*
     &    a1pp + 2.D0*C00s(t1)*taucg**(-1)*taugs*a3pp - 4.D0*C00s(t1b)*
     &    epinv*taucs*taucg**(-2)*taugs*a1pp - 4.D0*C00s(t1b)*epinv*
     &    taucg**(-2)*taugs*a2pp - 2.D0*C00s(t1b)*epinv*taucg**(-2)*
     &    taugs**2*a1pp - 2.D0*C00s(t1b)*epinv*taucg**(-2)*taugs**2*
     &    a3pp + 4.D0*C00s(t2b)*epinv*taucs*taucg**(-2)*taugs*a1pp + 4.D
     &    0*C00s(t2b)*epinv*taucg**(-2)*taugs*a2pp + 4.D0*C00s(t2b)*
     &    epinv*taucg**(-2)*taugs**2*a1pp )
      spp = spp + mc**3*xn * ( 2.D0*C00s(t2b)*epinv*taucg**(-2)*
     &    taugs**2*a3pp + 10.D0*C00s(t3)*epinv*taucg**(-2)*taugs*a2pp
     &     - 2.D0*C00s(t3)*epinv*taucg**(-2)*taugs**2*a1pp - 10.D0*
     &    C00s(t3)*epinv*taucg**(-1)*taugs*a1pp - 4.D0*C00s(t3)*
     &    taucg**(-2)*taugs*a2pp + 6.D0*C00s(t3)*taucg**(-1)*taugs*a1pp
     &     + 2.D0*C00f(t1)*taucg**(-2)*taugs*a2pp - 2.D0*C00f(t1)*
     &    taucg**(-1)*taugs*a1pp - 4.D0*C00f(t1b)*taucs*taucg**(-2)*
     &    taugs*a1pp - 4.D0*C00f(t1b)*taucg**(-2)*taugs*a2pp - 2.D0*
     &    C00f(t1b)*taucg**(-2)*taugs**2*a1pp - 2.D0*C00f(t1b)*
     &    taucg**(-2)*taugs**2*a3pp + 4.D0*C00f(t2b)*taucs*taucg**(-2)*
     &    taugs*a1pp + 4.D0*C00f(t2b)*taucg**(-2)*taugs*a2pp + 4.D0*
     &    C00f(t2b)*taucg**(-2)*taugs**2*a1pp + 2.D0*C00f(t2b)*
     &    taucg**(-2)*taugs**2*a3pp + 10.D0*C00f(t3)*taucg**(-2)*taugs*
     &    a2pp - 2.D0*C00f(t3)*taucg**(-2)*taugs**2*a1pp - 10.D0*C00f(
     &    t3)*taucg**(-1)*taugs*a1pp + 3.D0/2.D0*lnrat(musq,msq)*
     &    taucg**(-2)*taugs*a2pp )
      spp = spp + mc**3*xn * ( 3.D0/2.D0*lnrat(musq,msq)*taucg**(-1)*
     &    taugs*a3pp )
      spp = spp + mc**3 * (  - epinv*taucg**(-2)*taugs*b0*a2pp + epinv*
     &    taucg**(-1)*taugs*b0*a1pp )
      spp = spp + BB1(s1)*mc*xn**(-1) * ( taucs*taucg**(-1)*a2pp - 
     &    taucs*a1pp )
      spp = spp + BB1(s1)*mc*xn * (  - taucs*taucg**(-1)*a2pp + taucs*
     &    a1pp )
      spp = spp + BB1(s1)*mc**3*xn**(-1) * (  - taucg**(-2)*taugs*a2pp
     &     + taucg**(-1)*taugs*a1pp )
      spp = spp + BB1(s1)*mc**3*xn * ( taucg**(-2)*taugs*a2pp - 
     &    taucg**(-1)*taugs*a1pp )
      spp = spp + BB1(s2)*mc*xn**(-1) * ( taucs*taucg**(-1)*a2pp + 
     &    taucs*a3pp )
      spp = spp + BB1(s2)*mc*xn * (  - taucs*taucg**(-1)*a2pp - taucs*
     &    a3pp )
      spp = spp + BB1(s2)*mc**3*xn**(-1) * (  - taucg**(-2)*taugs*a2pp
     &     - taucg**(-1)*taugs*a3pp )
      spp = spp + BB1(s2)*mc**3*xn * ( taucg**(-2)*taugs*a2pp + 
     &    taucg**(-1)*taugs*a3pp )
      spp = spp + BB0(s2)*mc*xn**(-1) * (  - taucs*taucg**(-1)*a2pp - 
     &    taucs*a3pp )
      spp = spp + BB0(s2)*mc*xn * ( taucs*taucg**(-1)*a2pp + taucs*a3pp
     &     )
      spp = spp + BB0(s2)*mc**3*xn**(-1) * ( taucg**(-2)*taugs*a2pp + 
     &    taucg**(-1)*taugs*a3pp )
      spp = spp + BB0(s2)*mc**3*xn * (  - taucg**(-2)*taugs*a2pp - 
     &    taucg**(-1)*taugs*a3pp )
      spp = spp + C11(t1)*mc*xn**(-1) * (  - taucs*taucg*a1pp - taucs*
     &    taucg*a3pp - taucs**2*a1pp - taucs**2*a3pp )
      spp = spp + C11(t1)*mc*xn * ( taucs*taucg*a1pp + taucs*taucg*a3pp
     &     + taucs**2*a1pp + taucs**2*a3pp )
      spp = spp + C11(t1)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pp + 2.D0*taucs*taucg**(-1)*taugs*a3pp + taucs*
     &    taucg**(-1)*a2pp - taucs*a1pp + taugs*a1pp + taugs*a3pp )
      spp = spp + C11(t1)*mc**3*xn * (  - 2.D0*taucs*taucg**(-1)*taugs*
     &    a1pp - 2.D0*taucs*taucg**(-1)*taugs*a3pp - taucs*taucg**(-1)*
     &    a2pp + taucs*a1pp - taugs*a1pp - taugs*a3pp )
      spp = spp + C11(t1)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs*a2pp
     &     - taucg**(-2)*taugs**2*a1pp - taucg**(-2)*taugs**2*a3pp + 
     &    taucg**(-1)*taugs*a1pp )
      spp = spp + C11(t1)*mc**5*xn * ( taucg**(-2)*taugs*a2pp + 
     &    taucg**(-2)*taugs**2*a1pp + taucg**(-2)*taugs**2*a3pp - 
     &    taucg**(-1)*taugs*a1pp )
      spp = spp + C11(t1b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pp + taucs*taucg**(-1)*taugs**2*a3pp )
      spp = spp + C11(t1b)*mc*xn * ( 2.D0*taucs*taucg**(-1)*taugs*a2pp
     &     + taucs*taucg**(-1)*taugs**2*a1pp + taucs*taucg**(-1)*
     &    taugs**2*a3pp + 2.D0*taucs**2*taucg**(-1)*taugs*a1pp )
      spp = spp + C11(t1b)*mc**3*xn**(-1) * (  - taucg**(-2)*taugs**3*
     &    a1pp - taucg**(-2)*taugs**3*a3pp )
      spp = spp + C11(t1b)*mc**3*xn * (  - 2.D0*taucs*taucg**(-2)*
     &    taugs**2*a1pp - 2.D0*taucg**(-2)*taugs**2*a2pp - taucg**(-2)*
     &    taugs**3*a1pp - taucg**(-2)*taugs**3*a3pp )
      spp = spp + C11(t2)*mc*xn**(-1) * (  - taucs*taucg*a3pp - taucs*
     &    taugs*a3pp - taucs**2*a3pp )
      spp = spp + C11(t2)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3pp - taucs*a3pp + taucg**(-1)*taugs**2*a3pp + taugs*a3pp )
      spp = spp + C11(t2)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a3pp )
      spp = spp + C11(t8b)*mc*xn**(-1) * ( taucs**2*a1pp + taucs**2*
     &    a3pp )
      spp = spp + C11(t8b)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*
     &    taugs*a1pp - taucs*taucg**(-1)*taugs*a3pp + taucs*a3pp )
      spp = spp + C11(t8b)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs*a3pp
     &     )
      spp = spp + C22(t1)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs*a2pp
     &     + taucs*taugs*a3pp )
      spp = spp + C22(t1)*mc*xn * (  - taucs*taucg**(-1)*taugs*a2pp - 
     &    taucs*taugs*a3pp )
      spp = spp + C22(t1)*mc**3*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a2pp - taucg**(-1)*taugs**2*a3pp )
      spp = spp + C22(t1)*mc**3*xn * ( taucg**(-2)*taugs**2*a2pp + 
     &    taucg**(-1)*taugs**2*a3pp )
      spp = spp + C22(t1b)*mc*xn**(-1) * ( taucs*taucg*a1pp + taucs*
     &    taucg*a3pp + taucs**2*a1pp + taucs**2*a3pp )
      spp = spp + C22(t1b)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*
     &    taugs*a1pp - taucs*taucg**(-1)*taugs*a3pp + taucs*a1pp + 
     &    taucs*a3pp - taugs*a1pp - taugs*a3pp )
      spp = spp + C22(t1b)*mc**3*xn * (  - taucs*a3pp + taucs**2*
     &    taucg**(-1)*a1pp )
      spp = spp + C22(t1b)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs*a1pp
     &     - taucg**(-1)*taugs*a3pp )
      spp = spp + C22(t1b)*mc**5*xn * (  - taucs*taucg**(-2)*taugs*a1pp
     &     + taucg**(-1)*taugs*a3pp )
      spp = spp + C22(t2b)*mc*xn * ( taucs*taucg*a3pp - taucs*taugs*
     &    a1pp - taucs**2*a1pp )
      spp = spp + C22(t2b)*mc**3*xn * ( taucs*a3pp - taucs**2*
     &    taucg**(-1)*a1pp + taucg**(-1)*taugs**2*a1pp - taugs*a3pp )
      spp = spp + C22(t2b)*mc**5*xn * ( taucs*taucg**(-2)*taugs*a1pp + 
     &    taucg**(-2)*taugs**2*a1pp - taucg**(-1)*taugs*a3pp )
      spp = spp + C22(t4b)*mc*xn**(-1) * ( taucs*taucg*a1pp + taucs*
     &    taucg*a3pp )
      spp = spp + C22(t4b)*mc**3*xn**(-1) * ( taucs*a1pp - taugs*a1pp
     &     - taugs*a3pp )
      spp = spp + C22(t4b)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs*a1pp
     &     )
      spp = spp + C22(t7)*mc*xn**(-1) * (  - taucs*taucg*a1pp - taucs*
     &    taugs*a1pp - taucs**2*taucg**(-1)*taugs*a1pp - taucs**2*a1pp
     &     )
      spp = spp + C22(t7)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs**2
     &    *a1pp - taucs*a1pp + taucg**(-1)*taugs**2*a1pp + taugs*a1pp )
      spp = spp + C22(t7)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a1pp
     &     + taucg**(-1)*taugs*a1pp )
      spp = spp + C22(t8)*mc*xn**(-1) * (  - taucs**2*a1pp - taucs**2*
     &    a3pp )
      spp = spp + C22(t8)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1pp + taucs*taucg**(-1)*taugs*a3pp - taucs*a3pp )
      spp = spp + C22(t8)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a3pp )
      spp = spp + C22(t9)*mc*xn**(-1) * (  - taucs*taucg*a1pp - taucs*
     &    taugs*a1pp - taucs**2*taucg**(-1)*taugs*a1pp - taucs**2*a1pp
     &     )
      spp = spp + C22(t9)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs**2
     &    *a1pp - taucs*a1pp + taucg**(-1)*taugs**2*a1pp + taugs*a1pp )
      spp = spp + C22(t9)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a1pp
     &     + taucg**(-1)*taugs*a1pp )
      spp = spp + C12(t1)*mc*xn**(-1) * ( taucs*taucg*a3pp - taucs*
     &    taugs*a1pp - taucs*taugs*a3pp + taucs*a2pp + taucs**2*
     &    taucg**(-1)*a2pp + taucs**2*a3pp )
      spp = spp + C12(t1)*mc*xn * (  - taucs*taucg*a3pp + taucs*taugs*
     &    a1pp + taucs*taugs*a3pp - taucs*a2pp - taucs**2*taucg**(-1)*
     &    a2pp - taucs**2*a3pp )
      spp = spp + C12(t1)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a2pp - taucs*taucg**(-1)*taugs*a1pp - 2.D0*taucs*taucg**(-1)
     &    *taugs*a3pp - taucg**(-1)*taugs*a2pp + taucg**(-1)*taugs**2*
     &    a1pp + taucg**(-1)*taugs**2*a3pp - taugs*a3pp )
      spp = spp + C12(t1)*mc**3*xn * ( taucs*taucg**(-2)*taugs*a2pp + 
     &    taucs*taucg**(-1)*taugs*a1pp + 2.D0*taucs*taucg**(-1)*taugs*
     &    a3pp + taucg**(-1)*taugs*a2pp - taucg**(-1)*taugs**2*a1pp - 
     &    taucg**(-1)*taugs**2*a3pp + taugs*a3pp )
      spp = spp + C12(t1)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a1pp
     &     + taucg**(-2)*taugs**2*a3pp )
      spp = spp + C12(t1)*mc**5*xn * (  - taucg**(-2)*taugs**2*a1pp - 
     &    taucg**(-2)*taugs**2*a3pp )
      spp = spp + C12(t1b)*mc*xn**(-1) * ( 2.D0*taucs*taugs*a1pp + 2.D0
     &    *taucs*taugs*a3pp + taucs**2*taucg**(-1)*taugs*a1pp + 
     &    taucs**2*taucg**(-1)*taugs*a3pp )
      spp = spp + C12(t1b)*mc*xn * ( taucs*taugs*a1pp + taucs*taugs*
     &    a3pp + taucs*a2pp + taucs**2*taucg**(-1)*a2pp + taucs**2*a1pp
     &     + taucs**3*taucg**(-1)*a1pp )
      spp = spp + C12(t1b)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1pp - taucs*taucg**(-2)*taugs**2*a3pp + taucs*
     &    taucg**(-1)*taugs*a1pp + taucs*taucg**(-1)*taugs*a3pp - 2.D0*
     &    taucg**(-1)*taugs**2*a1pp - 2.D0*taucg**(-1)*taugs**2*a3pp )
      spp = spp + C12(t1b)*mc**3*xn * (  - taucs*taucg**(-2)*taugs*a2pp
     &     - taucs*taucg**(-1)*taugs*a3pp - taucs**2*taucg**(-2)*taugs*
     &    a1pp - taucg**(-1)*taugs*a2pp - taucg**(-1)*taugs**2*a1pp - 
     &    taucg**(-1)*taugs**2*a3pp )
      spp = spp + C12(t1b)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a1pp - taucg**(-2)*taugs**2*a3pp )
      spp = spp + C12(t1b)*mc**5*xn * (  - taucg**(-2)*taugs**2*a1pp + 
     &    taucg**(-2)*taugs**2*a3pp )
      spp = spp + C12(t2)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2
     &    *a3pp - taucs*taugs*a3pp - taucs**2*taucg**(-1)*taugs*a3pp )
      spp = spp + C12(t2)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs**2
     &    *a3pp - taucs*taucg**(-1)*taugs*a3pp + taucg**(-2)*taugs**3*
     &    a3pp + taucg**(-1)*taugs**2*a3pp )
      spp = spp + C12(t2)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a3pp
     &     )
      spp = spp + C12(t2b)*mc*xn * (  - taucs*taucg**(-1)*taugs*a2pp - 
     &    taucs*taucg**(-1)*taugs**2*a1pp + taucs*taugs*a3pp - 2.D0*
     &    taucs**2*taucg**(-1)*taugs*a1pp - taucs**2*taucg**(-1)*a2pp
     &     - taucs**3*taucg**(-1)*a1pp )
      spp = spp + C12(t2b)*mc**3*xn * ( taucs*taucg**(-2)*taugs*a2pp + 
     &    2.D0*taucs*taucg**(-2)*taugs**2*a1pp + taucs*taucg**(-1)*
     &    taugs*a3pp + taucs**2*taucg**(-2)*taugs*a1pp + taucg**(-2)*
     &    taugs**2*a2pp + taucg**(-2)*taugs**3*a1pp - taucg**(-1)*
     &    taugs**2*a3pp )
      spp = spp + C12(t2b)*mc**5*xn * (  - taucg**(-2)*taugs**2*a3pp )
      spp = spp + C12(t3)*mc*xn**(-1) * (  - taucs*taugs*a1pp )
      spp = spp + C12(t3)*mc*xn * (  - 3.D0*taucs*taucg**(-1)*taugs*
     &    a2pp + taucs*taucg**(-1)*taugs**2*a1pp + 2.D0*taucs*taugs*
     &    a1pp )
      spp = spp + C12(t3)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pp
     &     )
      spp = spp + C12(t3)*mc**3*xn * ( 3.D0*taucg**(-2)*taugs**2*a2pp
     &     - taucg**(-2)*taugs**3*a1pp - 2.D0*taucg**(-1)*taugs**2*a1pp
     &     )
      spp = spp + C12(t3b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pp )
      spp = spp + C12(t3b)*mc**3*xn**(-1) * (  - taucg**(-2)*taugs**3*
     &    a1pp )
      spp = spp + C12(t7)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2
     &    *a1pp - taucs*taugs*a1pp )
      spp = spp + C12(t7)*mc**3*xn**(-1) * ( taucg**(-2)*taugs**3*a1pp
     &     + taucg**(-1)*taugs**2*a1pp )
      spp = spp + C12(t8)*mc*xn**(-1) * (  - taucs**2*taucg**(-1)*taugs
     &    *a1pp - taucs**2*taucg**(-1)*taugs*a3pp )
      spp = spp + C12(t8)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs**2
     &    *a1pp + taucs*taucg**(-2)*taugs**2*a3pp - taucs*taucg**(-1)*
     &    taugs*a3pp )
      spp = spp + C12(t8)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a3pp
     &     )
      spp = spp + C12(t8b)*mc*xn**(-1) * ( taucs**2*taucg**(-1)*taugs*
     &    a1pp + taucs**2*taucg**(-1)*taugs*a3pp )
      spp = spp + C12(t8b)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1pp - taucs*taucg**(-2)*taugs**2*a3pp + taucs*
     &    taucg**(-1)*taugs*a3pp )
      spp = spp + C12(t8b)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a3pp )
      spp = spp + C12(t9)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2
     &    *a1pp - taucs*taugs*a1pp )
      spp = spp + C12(t9)*mc**3*xn**(-1) * ( taucg**(-2)*taugs**3*a1pp
     &     + taucg**(-1)*taugs**2*a1pp )
      spp = spp + C0(t2)*mc**3*xn**(-1) * (  - 3.D0*taucs*taucg**(-1)*
     &    taugs*a1pp - taucs*taucg**(-1)*taugs*a3pp - taucs*taucg**(-1)
     &    *a2pp + taucs*a3pp - taucs**2*taucg**(-1)*a1pp )
      spp = spp + C0(t2)*mc**5*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a1pp + taucg**(-2)*taugs*a2pp + 3.D0*taucg**(-2)*taugs**2*
     &    a1pp + taucg**(-2)*taugs**2*a3pp - taucg**(-1)*taugs*a3pp )
      spp = spp + C0(t3)*mc*xn**(-1) * (  - taucs*taugs*a1pp )
      spp = spp + C0(t3)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pp )
      spp = spp + C0(t4)*mc*xn**(-1) * (  - taucs*taucg*a3pp - taucs*
     &    a2pp )
      spp = spp + C0(t4)*mc**3*xn**(-1) * ( taucg**(-1)*taugs*a2pp + 
     &    taugs*a3pp )
      spp = spp + C0(t4b)*mc*xn**(-1) * (  - taucs*taucg*a3pp - taucs*
     &    taugs*a1pp - taucs*taugs*a3pp - taucs*a2pp )
      spp = spp + C0(t4b)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs
     &    *a1pp + taucg**(-1)*taugs*a2pp + taucg**(-1)*taugs**2*a1pp + 
     &    taucg**(-1)*taugs**2*a3pp + taugs*a3pp )
      spp = spp + C0(t4b)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a1pp
     &     )
      spp = spp + C0(t6)*mc*xn * ( taucs*taucg*a3pp + taucs*a2pp )
      spp = spp + C0(t6)*mc**3*xn * (  - taucg**(-1)*taugs*a2pp - taugs
     &    *a3pp )
      spp = spp + C0(t6b)*mc*xn * ( taucs*taucg**(-1)*taugs**2*a1pp + 
     &    taucs*taucg**(-1)*taugs**2*a3pp + taucs**2*taucg**(-1)*taugs*
     &    a1pp )
      spp = spp + C0(t6b)*mc**3*xn * (  - taucs*taucg**(-2)*taugs**2*
     &    a1pp - taucg**(-2)*taugs**3*a1pp - taucg**(-2)*taugs**3*a3pp
     &     )
      spp = spp + C0(t8)*mc*xn**(-1) * (  - taucs*a2pp - taucs**2*
     &    taucg**(-1)*taugs*a1pp )
      spp = spp + C0(t8)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs**2*
     &    a1pp + taucg**(-1)*taugs*a2pp )
      spp = spp + C0(t8b)*mc**3*xn**(-1) * ( 3.D0*taucs*taucg**(-1)*
     &    taugs*a1pp + taucs*taucg**(-1)*taugs*a3pp + taucs*taucg**(-1)
     &    *a2pp - taucs*a3pp + taucs**2*taucg**(-1)*a1pp )
      spp = spp + C0(t8b)*mc**5*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a1pp - taucg**(-2)*taugs*a2pp - 3.D0*taucg**(-2)*taugs**2*
     &    a1pp - taucg**(-2)*taugs**2*a3pp + taucg**(-1)*taugs*a3pp )
      spp = spp + C0(t9)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pp )
      spp = spp + C0(t9)*mc**5*xn**(-1) * (  - 2.D0*taucg**(-2)*
     &    taugs**2*a1pp )
      spp = spp + C1(t1)*mc*xn**(-1) * (  - 2.D0*taucs*taucg*a1pp - 
     &    taucs*taucg*a3pp + taucs*a2pp + taucs**2*taucg**(-1)*a2pp - 2.
     &    D0*taucs**2*a1pp - taucs**2*a3pp )
      spp = spp + C1(t1)*mc*xn * ( 2.D0*taucs*taucg*a1pp + taucs*taucg*
     &    a3pp - taucs*a2pp - taucs**2*taucg**(-1)*a2pp + 2.D0*taucs**2
     &    *a1pp + taucs**2*a3pp )
      spp = spp + C1(t1)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a2pp + 4.D0*taucs*taucg**(-1)*taugs*a1pp + 2.D0*taucs*
     &    taucg**(-1)*taugs*a3pp + taucs*taucg**(-1)*a2pp - taucs*a1pp
     &     - taucg**(-1)*taugs*a2pp + 2.D0*taugs*a1pp + taugs*a3pp )
      spp = spp + C1(t1)*mc**3*xn * ( taucs*taucg**(-2)*taugs*a2pp - 4.D
     &    0*taucs*taucg**(-1)*taugs*a1pp - 2.D0*taucs*taucg**(-1)*taugs
     &    *a3pp - taucs*taucg**(-1)*a2pp + taucs*a1pp + taucg**(-1)*
     &    taugs*a2pp - 2.D0*taugs*a1pp - taugs*a3pp )
      spp = spp + C1(t1)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs*a2pp
     &     - 2.D0*taucg**(-2)*taugs**2*a1pp - taucg**(-2)*taugs**2*a3pp
     &     + taucg**(-1)*taugs*a1pp )
      spp = spp + C1(t1)*mc**5*xn * ( taucg**(-2)*taugs*a2pp + 2.D0*
     &    taucg**(-2)*taugs**2*a1pp + taucg**(-2)*taugs**2*a3pp - 
     &    taucg**(-1)*taugs*a1pp )
      spp = spp + C1(t1b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pp + taucs*taucg**(-1)*taugs**2*a3pp + taucs*taucg*a1pp - 
     &    taucs*a2pp - taucs**2*taucg**(-1)*a2pp + taucs**2*a1pp )
      spp = spp + C1(t1b)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pp + 2.D0
     &    *taucs*taucg**(-1)*taugs**2*a1pp + 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a3pp + taucs*taucg*a3pp + taucs*taugs*a1pp + taucs*
     &    taugs*a3pp + taucs*a2pp + 2.D0*taucs**2*taucg**(-1)*taugs*
     &    a1pp + taucs**2*taucg**(-1)*taugs*a3pp + taucs**2*taucg**(-1)
     &    *a2pp + taucs**2*a3pp )
      spp = spp + C1(t1b)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a2pp - 2.D0*taucs*taucg**(-1)*taugs*a1pp - taucg**(-2)*
     &    taugs**3*a1pp - taucg**(-2)*taugs**3*a3pp + taucg**(-1)*taugs
     &    *a2pp - taugs*a1pp )
      spp = spp + C1(t1b)*mc**3*xn * (  - taucs*taucg**(-2)*taugs*a2pp
     &     - 2.D0*taucs*taucg**(-2)*taugs**2*a1pp - taucs*taucg**(-2)*
     &    taugs**2*a3pp + taucs*taucg**(-1)*taugs*a1pp - taucs*
     &    taucg**(-1)*taugs*a3pp - taucg**(-2)*taugs**2*a2pp - 2.D0*
     &    taucg**(-2)*taugs**3*a1pp - 2.D0*taucg**(-2)*taugs**3*a3pp - 
     &    taucg**(-1)*taugs*a2pp - taucg**(-1)*taugs**2*a1pp - 
     &    taucg**(-1)*taugs**2*a3pp - taugs*a3pp )
      spp = spp + C1(t1b)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a1pp
     &     )
      spp = spp + C1(t1b)*mc**5*xn * (  - taucg**(-2)*taugs**2*a1pp )
      spp = spp + C1(t2)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs*a2pp
     &     - 2.D0*taucs*taucg*a3pp - taucs*taugs*a1pp - 2.D0*taucs*
     &    taugs*a3pp + taucs**2*taucg**(-1)*a2pp - taucs**2*a1pp - 2.D0
     &    *taucs**2*a3pp )
      spp = spp + C1(t2)*mc*xn * ( taucs*taugs*a3pp )
      spp = spp + C1(t2)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a2pp + taucs*taucg**(-1)*taugs*a3pp - taucs*taucg**(-1)*a2pp
     &     - 3.D0*taucs*a3pp - taucs**2*taucg**(-1)*a1pp - taucg**(-2)*
     &    taugs**2*a2pp + taucg**(-1)*taugs**2*a1pp + 2.D0*taucg**(-1)*
     &    taugs**2*a3pp + 2.D0*taugs*a3pp )
      spp = spp + C1(t2)*mc**3*xn * (  - taucg**(-1)*taugs**2*a3pp )
      spp = spp + C1(t2)*mc**5*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a1pp + taucg**(-2)*taugs*a2pp + taucg**(-2)*taugs**2*a1pp + 
     &    taucg**(-2)*taugs**2*a3pp + 3.D0*taucg**(-1)*taugs*a3pp )
      spp = spp + C1(t2b)*mc*xn * (  - taucs*taucg**(-1)*taugs*a2pp - 
     &    taucs*taucg**(-1)*taugs**2*a3pp - 2.D0*taucs*taugs*a3pp - 
     &    taucs**2*taucg**(-1)*taugs*a3pp - taucs**2*taucg**(-1)*a2pp
     &     - taucs**2*a3pp )
      spp = spp + C1(t2b)*mc**3*xn * ( taucs*taucg**(-2)*taugs*a2pp + 
     &    taucs*taucg**(-2)*taugs**2*a3pp + taucs*taucg**(-1)*taugs*
     &    a3pp + taucg**(-2)*taugs**2*a2pp + taucg**(-2)*taugs**3*a3pp
     &     + 2.D0*taucg**(-1)*taugs**2*a3pp )
      spp = spp + C1(t3)*mc*xn**(-1) * (  - taucs*taugs*a1pp )
      spp = spp + C1(t3)*mc*xn * ( taucs*taucg**(-1)*taugs**2*a3pp + 
     &    taucs*taugs*a1pp + taucs*taugs*a3pp )
      spp = spp + C1(t3)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pp )
      spp = spp + C1(t3)*mc**3*xn * (  - taucg**(-2)*taugs**3*a3pp - 
     &    taucg**(-1)*taugs**2*a1pp - taucg**(-1)*taugs**2*a3pp )
      spp = spp + C1(t3b)*mc*xn**(-1) * ( taucs*taugs*a1pp )
      spp = spp + C1(t3b)*mc**3*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1pp )
      spp = spp + C1(t4)*mc*xn**(-1) * (  - taucs*taucg*a3pp - taucs*
     &    a2pp )
      spp = spp + C1(t4)*mc**3*xn**(-1) * ( taucg**(-1)*taugs*a2pp + 
     &    taugs*a3pp )
      spp = spp + C1(t6)*mc*xn * ( taucs*taucg*a3pp + taucs*a2pp )
      spp = spp + C1(t6)*mc**3*xn * (  - taucg**(-1)*taugs*a2pp - taugs
     &    *a3pp )
      spp = spp + C1(t6b)*mc*xn * (  - taucs*taugs*a1pp )
      spp = spp + C1(t6b)*mc**3*xn * ( taucg**(-1)*taugs**2*a1pp )
      spp = spp + C1(t7)*mc*xn**(-1) * (  - taucs*taucg*a1pp - taucs*
     &    taugs*a1pp )
      spp = spp + C1(t7)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pp
     &     + taugs*a1pp )
      spp = spp + C1(t8)*mc*xn**(-1) * (  - taucs**2*taucg**(-1)*taugs*
     &    a1pp + taucs**2*taucg**(-1)*a2pp )
      spp = spp + C1(t8)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a2pp + taucs*taucg**(-2)*taugs**2*a1pp )
      spp = spp + C1(t8b)*mc*xn**(-1) * (  - taucs**2*taucg**(-1)*a2pp
     &     + 3.D0*taucs**2*a1pp + 2.D0*taucs**2*a3pp )
      spp = spp + C1(t8b)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a2pp - taucs*taucg**(-1)*taugs*a1pp - taucs*taucg**(-1)*taugs
     &    *a3pp + taucs*taucg**(-1)*a2pp + 3.D0*taucs*a3pp + taucs**2*
     &    taucg**(-1)*a1pp )
      spp = spp + C1(t8b)*mc**5*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a1pp - taucg**(-2)*taugs*a2pp - 2.D0*taucg**(-2)*taugs**2*
     &    a1pp - taucg**(-2)*taugs**2*a3pp - 3.D0*taucg**(-1)*taugs*
     &    a3pp )
      spp = spp + C2(t1)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs*a2pp
     &     + taucs*taucg*a1pp + taucs*taucg*a3pp - taucs*taugs*a1pp + 
     &    taucs**2*a1pp + taucs**2*a3pp )
      spp = spp + C2(t1)*mc*xn * (  - taucs*taucg**(-1)*taugs*a2pp - 
     &    taucs*taucg*a1pp - taucs*taucg*a3pp + taucs*taugs*a1pp - 
     &    taucs**2*a1pp - taucs**2*a3pp )
      spp = spp + C2(t1)*mc**3*xn**(-1) * (  - 3.D0*taucs*taucg**(-1)*
     &    taugs*a1pp - 2.D0*taucs*taucg**(-1)*taugs*a3pp - taucg**(-2)*
     &    taugs**2*a2pp + taucg**(-1)*taugs**2*a1pp - taugs*a1pp - 
     &    taugs*a3pp )
      spp = spp + C2(t1)*mc**3*xn * ( 3.D0*taucs*taucg**(-1)*taugs*a1pp
     &     + 2.D0*taucs*taucg**(-1)*taugs*a3pp + taucg**(-2)*taugs**2*
     &    a2pp - taucg**(-1)*taugs**2*a1pp + taugs*a1pp + taugs*a3pp )
      spp = spp + C2(t1)*mc**5*xn**(-1) * ( 2.D0*taucg**(-2)*taugs**2*
     &    a1pp + taucg**(-2)*taugs**2*a3pp )
      spp = spp + C2(t1)*mc**5*xn * (  - 2.D0*taucg**(-2)*taugs**2*a1pp
     &     - taucg**(-2)*taugs**2*a3pp )
      spp = spp + C2(t1b)*mc*xn**(-1) * ( taucs*taucg*a3pp + 2.D0*taucs
     &    *taugs*a1pp + 2.D0*taucs*taugs*a3pp + taucs*a2pp + taucs**2*
     &    taucg**(-1)*taugs*a1pp + taucs**2*taucg**(-1)*taugs*a3pp + 
     &    taucs**2*taucg**(-1)*a2pp + taucs**2*a3pp )
      spp = spp + C2(t1b)*mc*xn * (  - taucs*taucg*a3pp - taucs*a2pp - 
     &    taucs**2*taucg**(-1)*a2pp - taucs**2*a3pp )
      spp = spp + C2(t1b)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a2pp - taucs*taucg**(-2)*taugs**2*a1pp - taucs*taucg**(-2)*
     &    taugs**2*a3pp + 2.D0*taucs*taucg**(-1)*taugs*a1pp - taucs*
     &    taucg**(-1)*taugs*a3pp - taucs*taucg**(-1)*a2pp + taucs*a1pp
     &     - taucg**(-1)*taugs*a2pp - 2.D0*taucg**(-1)*taugs**2*a1pp - 
     &    2.D0*taucg**(-1)*taugs**2*a3pp - taugs*a3pp )
      spp = spp + C2(t1b)*mc**3*xn * ( taucs*taucg**(-2)*taugs*a2pp + 
     &    taucs*taucg**(-1)*taugs*a1pp + 3.D0*taucs*taucg**(-1)*taugs*
     &    a3pp + taucs*taucg**(-1)*a2pp + taucs*a1pp + 2.D0*taucs*a3pp
     &     + taucg**(-1)*taugs*a2pp + taugs*a3pp )
      spp = spp + C2(t1b)*mc**5*xn**(-1) * ( taucg**(-2)*taugs*a2pp - 2.
     &    D0*taucg**(-2)*taugs**2*a1pp - taucg**(-1)*taugs*a1pp )
      spp = spp + C2(t1b)*mc**5*xn * (  - taucg**(-2)*taugs*a2pp - 
     &    taucg**(-2)*taugs**2*a1pp - 2.D0*taucg**(-2)*taugs**2*a3pp - 
     &    taucg**(-1)*taugs*a1pp - 2.D0*taucg**(-1)*taugs*a3pp )
      spp = spp + C2(t2)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pp - taucs*taucg**(-1)*taugs**2*a1pp - taucs*taucg**(-1)*
     &    taugs**2*a3pp - taucs*taugs*a3pp - 2.D0*taucs**2*taucg**(-1)*
     &    taugs*a1pp - taucs**2*taucg**(-1)*taugs*a3pp - taucs**2*
     &    taucg**(-1)*a2pp - taucs**3*taucg**(-1)*a1pp )
      spp = spp + C2(t2)*mc*xn * (  - taucs*taugs*a3pp )
      spp = spp + C2(t2)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a2pp + 2.D0*taucs*taucg**(-2)*taugs**2*a1pp + taucs*
     &    taucg**(-2)*taugs**2*a3pp - taucs*taucg**(-1)*taugs*a3pp + 
     &    taucs**2*taucg**(-2)*taugs*a1pp + taucg**(-2)*taugs**2*a2pp
     &     + taucg**(-2)*taugs**3*a1pp + taucg**(-2)*taugs**3*a3pp + 
     &    taucg**(-1)*taugs**2*a3pp )
      spp = spp + C2(t2)*mc**3*xn * ( taucg**(-1)*taugs**2*a3pp )
      spp = spp + C2(t2)*mc**5*xn**(-1) * ( taucg**(-2)*taugs**2*a3pp )
      spp = spp + C2(t2b)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pp - 
     &    taucs*taucg*a3pp + taucs*taugs*a1pp + taucs*taugs*a3pp + 
     &    taucs**2*taucg**(-1)*a2pp + taucs**2*a1pp + taucs**2*a3pp )
      spp = spp + C2(t2b)*mc**3*xn * (  - taucs*taucg**(-2)*taugs*a2pp
     &     - taucs*taucg**(-1)*taugs*a1pp - 3.D0*taucs*taucg**(-1)*
     &    taugs*a3pp - taucs*taucg**(-1)*a2pp - 2.D0*taucs*a3pp - 
     &    taucg**(-2)*taugs**2*a2pp - taucg**(-1)*taugs**2*a1pp - 
     &    taucg**(-1)*taugs**2*a3pp + taugs*a3pp )
      spp = spp + C2(t2b)*mc**5*xn * ( taucg**(-2)*taugs*a2pp + 2.D0*
     &    taucg**(-2)*taugs**2*a3pp + 2.D0*taucg**(-1)*taugs*a3pp )
      spp = spp + C2(t3)*mc*xn**(-1) * (  - taucs*taugs*a1pp )
      spp = spp + C2(t3)*mc*xn * (  - taucs*taucg**(-1)*taugs*a2pp )
      spp = spp + C2(t3)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pp )
      spp = spp + C2(t3)*mc**3*xn * ( taucg**(-2)*taugs**2*a2pp )
      spp = spp + C2(t4b)*mc*xn**(-1) * ( 2.D0*taucs*taucg*a1pp + taucs
     &    *taucg*a3pp - taucs*a2pp )
      spp = spp + C2(t4b)*mc**3*xn**(-1) * ( 2.D0*taucs*a1pp + 
     &    taucg**(-1)*taugs*a2pp - 2.D0*taugs*a1pp - taugs*a3pp )
      spp = spp + C2(t4b)*mc**5*xn**(-1) * (  - 2.D0*taucg**(-1)*taugs*
     &    a1pp )
      spp = spp + C2(t6b)*mc*xn * ( taucs*a2pp + taucs**2*a1pp )
      spp = spp + C2(t6b)*mc**3*xn * (  - 3.D0*taucs*taucg**(-1)*taugs*
     &    a1pp - taucs*a1pp - taucg**(-1)*taugs*a2pp )
      spp = spp + C2(t6b)*mc**5*xn * ( 2.D0*taucg**(-2)*taugs**2*a1pp
     &     + taucg**(-1)*taugs*a1pp )
      spp = spp + C2(t7)*mc*xn**(-1) * (  - taucs**2*a1pp )
      spp = spp + C2(t7)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1pp - taucs*a1pp )
      spp = spp + C2(t7)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a1pp )
      spp = spp + C2(t8)*mc*xn**(-1) * (  - taucs*taucg*a3pp - taucs*
     &    taugs*a1pp - taucs*taugs*a3pp - taucs*a2pp - taucs**2*
     &    taucg**(-1)*taugs*a1pp - taucs**2*taucg**(-1)*taugs*a3pp - 
     &    taucs**2*taucg**(-1)*a2pp - taucs**2*a3pp )
      spp = spp + C2(t8)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a2pp + taucs*taucg**(-2)*taugs**2*a1pp + taucs*taucg**(-2)*
     &    taugs**2*a3pp - 2.D0*taucs*taucg**(-1)*taugs*a1pp + taucs*
     &    taucg**(-1)*taugs*a3pp + taucs*taucg**(-1)*a2pp + taucg**(-1)
     &    *taugs*a2pp + taucg**(-1)*taugs**2*a1pp + taucg**(-1)*
     &    taugs**2*a3pp + taugs*a3pp )
      spp = spp + C2(t8)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs*a2pp
     &     + 2.D0*taucg**(-2)*taugs**2*a1pp )
      spp = spp + C2(t8b)*mc*xn**(-1) * ( 2.D0*taucs**2*taucg**(-1)*
     &    taugs*a1pp + taucs**2*taucg**(-1)*taugs*a3pp + taucs**2*
     &    taucg**(-1)*a2pp + taucs**3*taucg**(-1)*a1pp )
      spp = spp + C2(t8b)*mc**3*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a2pp - 2.D0*taucs*taucg**(-2)*taugs**2*a1pp - taucs*
     &    taucg**(-2)*taugs**2*a3pp + taucs*taucg**(-1)*taugs*a3pp - 
     &    taucs**2*taucg**(-2)*taugs*a1pp )
      spp = spp + C2(t8b)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a3pp )
      spp = spp + C2(t9)*mc*xn**(-1) * (  - 2.D0*taucs*taucg*a1pp - 2.D0
     &    *taucs*taugs*a1pp - 2.D0*taucs**2*a1pp )
      spp = spp + C2(t9)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pp - 2.D0*taucs*a1pp + 2.D0*taucg**(-1)*taugs**2*a1pp
     &     + 2.D0*taugs*a1pp )
      spp = spp + C2(t9)*mc**5*xn**(-1) * ( 2.D0*taucg**(-1)*taugs*a1pp
     &     )
      spp = spp + D0(b2)*mc*xn**(-1) * ( taucs*taucg*taugs*a1pp + 2.D0*
     &    taucs*taucg*taugs*a3pp + taucs*taucg*a2pp + taucs*taucg**2*
     &    a3pp + taucs*taugs*a2pp + taucs*taugs**2*a1pp + taucs*
     &    taugs**2*a3pp + taucs**2*taucg*a3pp + taucs**2*taugs*a1pp + 
     &    taucs**2*taugs*a3pp + taucs**2*a2pp )
      spp = spp + D0(b2)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pp - taucs*taucg**(-1)*taugs**2*a3pp + taucs*taucg*a3pp - 
     &    taucs*a2pp + taucs**2*taucg**(-1)*taugs*a1pp - 2.D0*taucs**2*
     &    a1pp - taucg**(-1)*taugs**2*a2pp - taucg**(-1)*taugs**3*a1pp
     &     - taucg**(-1)*taugs**3*a3pp - taucg*taugs*a3pp - taugs*a2pp
     &     - taugs**2*a1pp - 2.D0*taugs**2*a3pp )
      spp = spp + D0(b2)*mc**5*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1pp - taucg**(-2)*taugs**3*a1pp + taucg**(-1)*taugs
     &    *a2pp - taucg**(-1)*taugs**2*a3pp - taugs*a3pp )
      spp = spp + D0(b2)*mc**7*xn**(-1) * ( 2.D0*taucg**(-2)*taugs**2*
     &    a1pp )
      spp = spp + D1(b1)*mc*xn**(-1) * (  - taucs**2*taugs*a1pp )
      spp = spp + D1(b1)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pp )
      spp = spp + D1(b2)*mc*xn**(-1) * ( taucs*taucg*taugs*a1pp + 2.D0*
     &    taucs*taucg*taugs*a3pp + taucs*taucg*a2pp + taucs*taucg**2*
     &    a3pp + taucs*taugs*a2pp + taucs*taugs**2*a1pp + taucs*
     &    taugs**2*a3pp )
      spp = spp + D1(b2)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pp - taucg**(-1)*taugs**2*a2pp - taucg**(-1)*taugs**3*a1pp
     &     - taucg**(-1)*taugs**3*a3pp - taucg*taugs*a3pp - taugs*a2pp
     &     - taugs**2*a1pp - 2.D0*taugs**2*a3pp )
      spp = spp + D1(b2)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs**3*
     &    a1pp )
      spp = spp + D1(b3)*mc*xn * (  - taucs*taucg*taugs*a3pp - taucs*
     &    taugs*a2pp - taucs*taugs**2*a3pp )
      spp = spp + D1(b3)*mc**3*xn * ( taucg**(-1)*taugs**2*a2pp + 
     &    taucg**(-1)*taugs**3*a3pp + taugs**2*a3pp )
      spp = spp + D2(b1)*mc*xn**(-1) * (  - taucs*taugs*a2pp )
      spp = spp + D2(b1)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a2pp )
      spp = spp + D2(b2)*mc*xn**(-1) * (  - 2.D0*taucs**2*taucg*a1pp + 
     &    taucs**2*taugs*a1pp + taucs**2*taugs*a3pp + 2.D0*taucs**2*
     &    a2pp )
      spp = spp + D2(b2)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a2pp - taucs*taucg**(-1)*taugs**2*a1pp - taucs*
     &    taucg**(-1)*taugs**2*a3pp - 2.D0*taucs*taucg*a3pp + 3.D0*
     &    taucs*taugs*a1pp + taucs*taugs*a3pp - taucs*a2pp + taucs**2*
     &    taucg**(-1)*taugs*a1pp - 2.D0*taucs**2*a1pp )
      spp = spp + D2(b2)*mc**5*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1pp + taucg**(-1)*taugs*a2pp - taucg**(-1)*taugs**2
     &    *a1pp - taucg**(-1)*taugs**2*a3pp + 2.D0*taugs*a3pp )
      spp = spp + D2(b2)*mc**7*xn**(-1) * ( 2.D0*taucg**(-2)*taugs**2*
     &    a1pp )
      spp = spp + D2(b3)*mc*xn * ( taucs*taucg**(-1)*taugs**3*a1pp + 
     &    taucs*taucg**(-1)*taugs**3*a3pp + taucs*taucg*taugs*a1pp + 
     &    taucs*taucg*taugs*a3pp + taucs*taugs**2*a1pp + taucs*taugs**2
     &    *a3pp + taucs**2*taucg**(-1)*taugs**2*a1pp )
      spp = spp + D2(b3)*mc**3*xn * (  - taucs*taucg**(-2)*taugs**3*
     &    a1pp - taucg**(-2)*taugs**4*a1pp - taucg**(-2)*taugs**4*a3pp
     &     - taucg**(-1)*taugs**3*a1pp - taucg**(-1)*taugs**3*a3pp - 
     &    taugs**2*a1pp - taugs**2*a3pp )
      spp = spp + D3(b1)*mc*xn**(-1) * (  - taucs*taucg*a2pp - taucs**2
     &    *taucg*a1pp - taucs**2*a2pp - taucs**3*a1pp )
      spp = spp + D3(b1)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a2pp + 2.D0*taucs*taugs*a1pp + taucs**2*taucg**(-1)*taugs*
     &    a1pp + taugs*a2pp )
      spp = spp + D3(b1)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1pp )
      spp = spp + D3(b2)*mc*xn**(-1) * ( taucs**2*taucg*a3pp - 2.D0*
     &    taucs**2*taugs*a1pp - taucs**2*a2pp - 2.D0*taucs**3*a1pp )
      spp = spp + D3(b2)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a2pp + 3.D0*taucs*taucg**(-1)*taugs**2*a1pp - 4.D0*taucs*
     &    taugs*a3pp + taucs**2*taucg**(-1)*taugs*a1pp )
      spp = spp + D3(b2)*mc**5*xn**(-1) * ( taucs*taucg**(-2)*taugs**2*
     &    a1pp - taucg**(-2)*taugs**3*a1pp + 3.D0*taucg**(-1)*taugs**2*
     &    a3pp )
      spp = spp + D3(b3)*mc*xn * (  - taucs*taucg*taugs*a3pp - taucs*
     &    taucg*a2pp - taucs*taucg**2*a3pp + taucs**2*taucg**(-1)*
     &    taugs**2*a1pp + taucs**2*taucg**(-1)*taugs**2*a3pp - taucs**2
     &    *taucg*a3pp + taucs**2*taugs*a1pp - taucs**2*taugs*a3pp - 
     &    taucs**2*a2pp + taucs**3*taucg**(-1)*taugs*a1pp )
      spp = spp + D3(b3)*mc**3*xn * (  - taucs*taucg**(-2)*taugs**3*
     &    a1pp - taucs*taucg**(-2)*taugs**3*a3pp + taucs*taucg**(-1)*
     &    taugs*a2pp - taucs*taucg**(-1)*taugs**2*a1pp + 3.D0*taucs*
     &    taucg**(-1)*taugs**2*a3pp + taucs*taugs*a1pp + 3.D0*taucs*
     &    taugs*a3pp - taucs**2*taucg**(-2)*taugs**2*a1pp + taucg*taugs
     &    *a3pp + taugs*a2pp + taugs**2*a3pp )
      spp = spp + D3(b3)*mc**5*xn * (  - 2.D0*taucg**(-2)*taugs**3*a3pp
     &     - taucg**(-1)*taugs**2*a1pp - 2.D0*taucg**(-1)*taugs**2*a3pp
     &     )

      smm=  + mc**2*xn**(-1) * ( 3.D0*epinv*taucs*taucg**(-1)*a3mm + 3.D
     &    0/4.D0*epinv*taucg**(-1)*taugs*a3mm - 3.D0/2.D0*epinv*taugs*
     &    a2mm + 5.D0*taucs*taucg**(-1)*a3mm + 5.D0/2.D0*taucg**(-1)*
     &    taugs*a3mm - 5.D0/2.D0*taugs*a2mm + 2.D0*C00s(t1b)*epinv*
     &    taugs*a1mm + 2.D0*C00s(t2)*epinv*taucg**(-1)*taugs*a3mm + 2.D0
     &    *C00s(t2)*taucs*taucg**(-1)*taugs*a1mm + 2.D0*C00s(t2)*
     &    taucg**(-1)*taugs*a3mm + 2.D0*C00s(t2)*taucg**(-1)*taugs**2*
     &    a1mm + 2.D0*C00s(t2)*taucg**(-1)*taugs**2*a2mm - 2.D0*C00s(t2
     &    )*taugs*a2mm + 2.D0*C00s(t3)*taugs*a1mm + 2.D0*C00s(t4)*epinv
     &    *taucg**(-1)*taugs*a3mm - 2.D0*C00s(t4)*taucg**(-1)*taugs*
     &    a3mm - 2.D0*C00s(t7)*epinv*taugs*a1mm + 2.D0*C00f(t1b)*taugs*
     &    a1mm + 2.D0*C00f(t2)*taucg**(-1)*taugs*a3mm + 2.D0*C00f(t4)*
     &    taucg**(-1)*taugs*a3mm - 2.D0*C00f(t7)*taugs*a1mm + 3.D0*
     &    lnrat(musq,msq)*taucs*taucg**(-1)*a3mm + 3.D0/2.D0*lnrat(musq
     &    ,msq)*taucg**(-1)*taugs*a3mm - 3.D0/2.D0*lnrat(musq,msq)*
     &    taugs*a2mm )
      smm = smm + mc**2*xn * (  - 3.D0*epinv*taucs*taucg**(-1)*a3mm - 3.
     &    D0/4.D0*epinv*taucg**(-1)*taugs*a3mm + 3.D0/2.D0*epinv*taugs*
     &    a2mm - 5.D0*taucs*taucg**(-1)*a3mm - 8.D0/3.D0*taucg**(-1)*
     &    taugs*a3mm + 5.D0/2.D0*taugs*a2mm + 2.D0*C00s(t1b)*epinv*
     &    taugs*a2mm - 2.D0*C00s(t2)*epinv*taucg**(-1)*taugs*a3mm - 2.D0
     &    *C00s(t2)*taucs*taucg**(-1)*taugs*a1mm - 2.D0*C00s(t2)*
     &    taucg**(-1)*taugs*a3mm - 2.D0*C00s(t2)*taucg**(-1)*taugs**2*
     &    a1mm - 2.D0*C00s(t2)*taucg**(-1)*taugs**2*a2mm + 2.D0*C00s(t2
     &    )*taugs*a2mm - 2.D0*C00s(t2b)*epinv*taugs*a2mm + 2.D0*C00s(t3
     &    )*taugs*a1mm - 10.D0*C00s(t6)*epinv*taucg**(-1)*taugs*a3mm + 
     &    4.D0*C00s(t6)*taucg**(-1)*taugs*a3mm + 2.D0*C00f(t1b)*taugs*
     &    a2mm - 2.D0*C00f(t2)*taucg**(-1)*taugs*a3mm - 2.D0*C00f(t2b)*
     &    taugs*a2mm - 10.D0*C00f(t6)*taucg**(-1)*taugs*a3mm - 3.D0*
     &    lnrat(musq,msq)*taucs*taucg**(-1)*a3mm - 3.D0/2.D0*lnrat(musq
     &    ,msq)*taucg**(-1)*taugs*a3mm + 3.D0/2.D0*lnrat(musq,msq)*
     &    taugs*a2mm )
      smm = smm + mc**2 * ( epinv*taucg**(-1)*taugs*b0*a3mm )
      smm = smm + mc**4*xn**(-1) * (  - 3.D0*epinv*taucg**(-2)*taugs*
     &    a3mm - 5.D0*taucg**(-2)*taugs*a3mm - 3.D0*lnrat(musq,msq)*
     &    taucg**(-2)*taugs*a3mm )
      smm = smm + mc**4*xn * ( 3.D0*epinv*taucg**(-2)*taugs*a3mm + 5.D0
     &    *taucg**(-2)*taugs*a3mm + 3.D0*lnrat(musq,msq)*taucg**(-2)*
     &    taugs*a3mm )
      smm = smm + xn**(-1) * ( 3.D0/4.D0*epinv*taucs*a3mm - 3.D0/4.D0*
     &    epinv*taucg*taugs*a2mm + 3.D0/4.D0*epinv*taugs*a3mm + 2.D0*
     &    C00s(t1b)*epinv*taucs*taucg*a2mm - 2.D0*C00s(t1b)*epinv*taucs
     &    *a3mm - 2.D0*C00s(t1b)*epinv*taucs**2*a1mm - 2.D0*C00s(t2)*
     &    epinv*taucs*taucg*a2mm + 4.D0*C00s(t2)*epinv*taucs*taugs*a1mm
     &     + 2.D0*C00s(t2)*epinv*taucs**2*a1mm + 2.D0*C00s(t2)*epinv*
     &    taugs*a3mm + 2.D0*C00s(t2)*epinv*taugs**2*a1mm + 2.D0*C00s(t2
     &    )*epinv*taugs**2*a2mm + 2.D0*C00s(t2)*taucs*taucg*a2mm - 4.D0
     &    *C00s(t2)*taucs*taugs*a1mm - 2.D0*C00s(t2)*taucs*taugs*a2mm
     &     - 2.D0*C00s(t2)*taucs*a3mm - 2.D0*C00s(t2)*taucs**2*a1mm - 2.
     &    D0*C00s(t2)*taugs*a3mm - 2.D0*C00s(t2)*taugs**2*a1mm - 2.D0*
     &    C00s(t2)*taugs**2*a2mm + 2.D0*C00s(t3)*taucg*taugs*a1mm + 2.D0
     &    *C00s(t3)*taucg*taugs*a2mm + 2.D0*C00s(t3)*taucg*a3mm + 6.D0*
     &    C00s(t3b)*epinv*taucs*taugs*a1mm - 2.D0*C00s(t3b)*epinv*taucg
     &    *taugs*a2mm + 6.D0*C00s(t3b)*epinv*taugs*a3mm + 4.D0*C00s(t3b
     &    )*epinv*taugs**2*a1mm )
      smm = smm + xn**(-1) * ( 4.D0*C00s(t3b)*epinv*taugs**2*a2mm - 2.D0
     &    *C00s(t4)*epinv*taucs*a3mm + 2.D0*C00s(t4)*epinv*taucg*taugs*
     &    a2mm - 2.D0*C00s(t4)*epinv*taugs*a3mm + 2.D0*C00s(t4)*taucs*
     &    a3mm - 6.D0*C00s(t7)*epinv*taucs*taugs*a1mm + 2.D0*C00s(t7)*
     &    epinv*taucg*taugs*a2mm - 6.D0*C00s(t7)*epinv*taugs*a3mm - 4.D0
     &    *C00s(t7)*epinv*taugs**2*a1mm - 4.D0*C00s(t7)*epinv*taugs**2*
     &    a2mm - 2.D0*C00s(t8)*epinv*taucs*taucg*a2mm + 2.D0*C00s(t8)*
     &    epinv*taucs*a3mm + 2.D0*C00s(t8)*epinv*taucs**2*a1mm + 2.D0*
     &    C00s(t8b)*epinv*taucs*taucg*a2mm - 2.D0*C00s(t8b)*epinv*taucs
     &    *a3mm - 2.D0*C00s(t8b)*epinv*taucs**2*a1mm - 4.D0*C00s(t9)*
     &    epinv*taucs*taugs*a1mm + 2.D0*C00s(t9)*epinv*taucg*taugs*a2mm
     &     - 4.D0*C00s(t9)*epinv*taugs*a3mm - 2.D0*C00s(t9)*epinv*
     &    taugs**2*a1mm - 2.D0*C00s(t9)*epinv*taugs**2*a2mm + 2.D0*
     &    C00f(t1b)*taucs*taucg*a2mm - 2.D0*C00f(t1b)*taucs*a3mm - 2.D0
     &    *C00f(t1b)*taucs**2*a1mm - 2.D0*C00f(t2)*taucs*taucg*a2mm + 4.
     &    D0*C00f(t2)*taucs*taugs*a1mm )
      smm = smm + xn**(-1) * ( 2.D0*C00f(t2)*taucs**2*a1mm + 2.D0*C00f(
     &    t2)*taugs*a3mm + 2.D0*C00f(t2)*taugs**2*a1mm + 2.D0*C00f(t2)*
     &    taugs**2*a2mm + 6.D0*C00f(t3b)*taucs*taugs*a1mm - 2.D0*C00f(
     &    t3b)*taucg*taugs*a2mm + 6.D0*C00f(t3b)*taugs*a3mm + 4.D0*
     &    C00f(t3b)*taugs**2*a1mm + 4.D0*C00f(t3b)*taugs**2*a2mm - 2.D0
     &    *C00f(t4)*taucs*a3mm + 2.D0*C00f(t4)*taucg*taugs*a2mm - 2.D0*
     &    C00f(t4)*taugs*a3mm - 6.D0*C00f(t7)*taucs*taugs*a1mm + 2.D0*
     &    C00f(t7)*taucg*taugs*a2mm - 6.D0*C00f(t7)*taugs*a3mm - 4.D0*
     &    C00f(t7)*taugs**2*a1mm - 4.D0*C00f(t7)*taugs**2*a2mm - 2.D0*
     &    C00f(t8)*taucs*taucg*a2mm + 2.D0*C00f(t8)*taucs*a3mm + 2.D0*
     &    C00f(t8)*taucs**2*a1mm + 2.D0*C00f(t8b)*taucs*taucg*a2mm - 2.D
     &    0*C00f(t8b)*taucs*a3mm - 2.D0*C00f(t8b)*taucs**2*a1mm - 4.D0*
     &    C00f(t9)*taucs*taugs*a1mm + 2.D0*C00f(t9)*taucg*taugs*a2mm - 
     &    4.D0*C00f(t9)*taugs*a3mm - 2.D0*C00f(t9)*taugs**2*a1mm - 2.D0
     &    *C00f(t9)*taugs**2*a2mm )
      smm = smm + xn * (  - 3.D0/4.D0*epinv*taucs*a3mm + 3.D0/4.D0*
     &    epinv*taucg*taugs*a2mm - 3.D0/4.D0*epinv*taugs*a3mm + 1.D0/6.D
     &    0*taucs*a3mm - 1.D0/6.D0*taucg*taugs*a2mm + 1.D0/6.D0*taugs*
     &    a3mm - 4.D0*C00s(t1b)*epinv*taucs*taucg*a2mm + 2.D0*C00s(t1b)
     &    *epinv*taucs*taugs*a1mm + 2.D0*C00s(t1b)*epinv*taucs*taugs*
     &    a2mm + 4.D0*C00s(t1b)*epinv*taucs*a3mm + 4.D0*C00s(t1b)*epinv
     &    *taucs**2*a1mm + 2.D0*C00s(t2)*epinv*taucs*a3mm - 2.D0*C00s(
     &    t2)*epinv*taucg*taugs*a2mm + 2.D0*C00s(t2)*epinv*taugs*a3mm
     &     - 2.D0*C00s(t2)*taucs*taucg*a2mm + 4.D0*C00s(t2)*taucs*taugs
     &    *a1mm + 2.D0*C00s(t2)*taucs*taugs*a2mm + 2.D0*C00s(t2)*taucs*
     &    a3mm + 2.D0*C00s(t2)*taucs**2*a1mm + 2.D0*C00s(t2)*taugs*a3mm
     &     + 2.D0*C00s(t2)*taugs**2*a1mm + 2.D0*C00s(t2)*taugs**2*a2mm
     &     + 4.D0*C00s(t2b)*epinv*taucs*taucg*a2mm - 8.D0*C00s(t2b)*
     &    epinv*taucs*taugs*a1mm - 2.D0*C00s(t2b)*epinv*taucs*taugs*
     &    a2mm - 4.D0*C00s(t2b)*epinv*taucs*a3mm - 4.D0*C00s(t2b)*epinv
     &    *taucs**2*a1mm )
      smm = smm + xn * ( 2.D0*C00s(t2b)*epinv*taucg*taugs*a2mm - 6.D0*
     &    C00s(t2b)*epinv*taugs*a3mm - 4.D0*C00s(t2b)*epinv*taugs**2*
     &    a1mm - 4.D0*C00s(t2b)*epinv*taugs**2*a2mm + 6.D0*C00s(t3)*
     &    epinv*taucs*taugs*a1mm - 2.D0*C00s(t3)*epinv*taucg*taugs*a2mm
     &     + 6.D0*C00s(t3)*epinv*taugs*a3mm + 4.D0*C00s(t3)*epinv*
     &    taugs**2*a1mm + 4.D0*C00s(t3)*epinv*taugs**2*a2mm + 2.D0*
     &    C00s(t3)*taucg*taugs*a1mm + 2.D0*C00s(t3)*taucg*taugs*a2mm + 
     &    2.D0*C00s(t3)*taucg*a3mm + 10.D0*C00s(t6)*epinv*taucs*a3mm - 
     &    10.D0*C00s(t6)*epinv*taucg*taugs*a2mm + 10.D0*C00s(t6)*epinv*
     &    taugs*a3mm - 4.D0*C00s(t6)*taucs*a3mm + 6.D0*C00s(t6)*taucg*
     &    taugs*a2mm - 6.D0*C00s(t6)*taugs*a3mm - 4.D0*C00f(t1b)*taucs*
     &    taucg*a2mm + 2.D0*C00f(t1b)*taucs*taugs*a1mm + 2.D0*C00f(t1b)
     &    *taucs*taugs*a2mm + 4.D0*C00f(t1b)*taucs*a3mm + 4.D0*C00f(t1b
     &    )*taucs**2*a1mm + 2.D0*C00f(t2)*taucs*a3mm - 2.D0*C00f(t2)*
     &    taucg*taugs*a2mm + 2.D0*C00f(t2)*taugs*a3mm + 4.D0*C00f(t2b)*
     &    taucs*taucg*a2mm )
      smm = smm + xn * (  - 8.D0*C00f(t2b)*taucs*taugs*a1mm - 2.D0*
     &    C00f(t2b)*taucs*taugs*a2mm - 4.D0*C00f(t2b)*taucs*a3mm - 4.D0
     &    *C00f(t2b)*taucs**2*a1mm + 2.D0*C00f(t2b)*taucg*taugs*a2mm - 
     &    6.D0*C00f(t2b)*taugs*a3mm - 4.D0*C00f(t2b)*taugs**2*a1mm - 4.D
     &    0*C00f(t2b)*taugs**2*a2mm + 6.D0*C00f(t3)*taucs*taugs*a1mm - 
     &    2.D0*C00f(t3)*taucg*taugs*a2mm + 6.D0*C00f(t3)*taugs*a3mm + 4.
     &    D0*C00f(t3)*taugs**2*a1mm + 4.D0*C00f(t3)*taugs**2*a2mm + 10.D
     &    0*C00f(t6)*taucs*a3mm - 10.D0*C00f(t6)*taucg*taugs*a2mm + 10.D
     &    0*C00f(t6)*taugs*a3mm )
      smm = smm + BB1(s2)*mc**2*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    a3mm - taugs*a2mm )
      smm = smm + BB1(s2)*mc**2*xn * (  - 2.D0*taucs*taucg**(-1)*a3mm
     &     + taugs*a2mm )
      smm = smm + BB1(s2)*mc**4*xn**(-1) * (  - 2.D0*taucg**(-2)*taugs*
     &    a3mm )
      smm = smm + BB1(s2)*mc**4*xn * ( 2.D0*taucg**(-2)*taugs*a3mm )
      smm = smm + BB1(s2)*xn**(-1) * ( taucs*a3mm - taucg*taugs*a2mm + 
     &    taugs*a3mm )
      smm = smm + BB1(s2)*xn * (  - taucs*a3mm + taucg*taugs*a2mm - 
     &    taugs*a3mm )
      smm = smm + BB0(s2)*mc**2*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    a3mm - 2.D0*taucg**(-1)*taugs*a3mm + taugs*a2mm )
      smm = smm + BB0(s2)*mc**2*xn * ( 2.D0*taucs*taucg**(-1)*a3mm + 2.D
     &    0*taucg**(-1)*taugs*a3mm - taugs*a2mm )
      smm = smm + BB0(s2)*mc**4*xn**(-1) * ( 2.D0*taucg**(-2)*taugs*
     &    a3mm )
      smm = smm + BB0(s2)*mc**4*xn * (  - 2.D0*taucg**(-2)*taugs*a3mm )
      smm = smm + BB0(s2)*xn**(-1) * ( taucs*a3mm - taucg*taugs*a2mm + 
     &    taugs*a3mm )
      smm = smm + BB0(s2)*xn * (  - taucs*a3mm + taucg*taugs*a2mm - 
     &    taugs*a3mm )
      smm = smm + C11(t1b)*mc**2*xn**(-1) * ( taugs**2*a1mm + taugs**2*
     &    a2mm )
      smm = smm + C11(t1b)*mc**2*xn * ( taugs**2*a1mm + taugs**2*a2mm )
      smm = smm + C11(t1b)*xn**(-1) * ( taucs*taugs**2*a1mm + taucs*
     &    taugs**2*a2mm )
      smm = smm + C11(t1b)*xn * (  - 2.D0*taucs*taucg*taugs*a2mm + 2.D0
     &    *taucs*taugs*a3mm + taucs*taugs**2*a1mm + taucs*taugs**2*a2mm
     &     + 2.D0*taucs**2*taugs*a1mm )
      smm = smm + C11(t2)*mc**2*xn**(-1) * (  - taucs*taucg*a2mm - 
     &    taucs*taugs*a2mm - taucs*a3mm + 2.D0*taugs*a3mm + 2.D0*
     &    taugs**2*a2mm )
      smm = smm + C11(t2)*mc**2*xn * ( taucs*taugs*a2mm + taucs*a3mm - 
     &    taucg*taugs*a2mm )
      smm = smm + C11(t2)*mc**4*xn**(-1) * ( taucg**(-1)*taugs*a3mm + 
     &    taucg**(-1)*taugs**2*a2mm )
      smm = smm + C11(t2)*mc**4*xn * (  - taucg**(-1)*taugs*a3mm - 
     &    taucg**(-1)*taugs**2*a2mm )
      smm = smm + C11(t2)*xn**(-1) * (  - 2.D0*taucs*taucg*taugs*a2mm
     &     - taucs*taucg*a3mm - taucs*taucg**2*a2mm - taucs**2*taucg*
     &    a2mm + taucg*taugs*a3mm + taucg*taugs**2*a2mm )
      smm = smm + C11(t2)*xn * ( taucs*taucg*a3mm + taucg*taugs*a3mm - 
     &    taucg**2*taugs*a2mm )
      smm = smm + C11(t4)*mc**2*xn**(-1) * ( taucs*a3mm + taucg*taugs*
     &    a2mm - taugs*a3mm )
      smm = smm + C11(t4)*mc**4*xn**(-1) * (  - taucg**(-1)*taugs*a3mm
     &     )
      smm = smm + C11(t4b)*xn**(-1) * (  - taucg**2*taugs*a1mm )
      smm = smm + C11(t6)*mc**2*xn * ( 4.D0*taucs*a3mm - 2.D0*taucg*
     &    taugs*a2mm + 2.D0*taugs*a3mm )
      smm = smm + C11(t6)*mc**4*xn * (  - 4.D0*taucg**(-1)*taugs*a3mm )
      smm = smm + C11(t8b)*mc**2*xn**(-1) * ( taucs*taucg*a2mm )
      smm = smm + C11(t8b)*xn**(-1) * ( taucs**2*taucg*a1mm + taucs**2*
     &    taucg*a2mm )
      smm = smm + C11(t9)*xn**(-1) * ( taucg*taugs**2*a1mm + taucg**2*
     &    taugs*a1mm )
      smm = smm + C22(t1b)*mc**2*xn**(-1) * ( taucs*taucg*a1mm + taucs*
     &    taucg*a2mm - taucg*taugs*a1mm - taucg*taugs*a2mm )
      smm = smm + C22(t1b)*mc**2*xn * (  - taucs*taucg*a2mm - taucs*
     &    taugs*a1mm - taucs*taugs*a2mm + taucs**2*a1mm )
      smm = smm + C22(t1b)*xn**(-1) * ( taucs*taucg**2*a1mm + taucs*
     &    taucg**2*a2mm + taucs**2*taucg*a1mm + taucs**2*taucg*a2mm )
      smm = smm + C22(t2b)*mc**2*xn * ( taucs*taucg*a2mm - 2.D0*taucs*
     &    taugs*a1mm + taucs*taugs*a2mm - taucs**2*a1mm + taucg*taugs*
     &    a2mm - 2.D0*taugs*a3mm - taugs**2*a1mm - taugs**2*a2mm )
      smm = smm + C22(t2b)*xn * (  - 2.D0*taucs*taucg*taugs*a1mm + 
     &    taucs*taucg*taugs*a2mm + taucs*taucg**2*a2mm - taucs**2*taucg
     &    *a1mm - 2.D0*taucg*taugs*a3mm - taucg*taugs**2*a1mm - taucg*
     &    taugs**2*a2mm + taucg**2*taugs*a2mm )
      smm = smm + C22(t4b)*mc**2*xn**(-1) * ( taucs*taucg*a1mm - taucg*
     &    taugs*a1mm - taucg*taugs*a2mm )
      smm = smm + C22(t4b)*mc**4*xn**(-1) * (  - 2.D0*taugs*a1mm )
      smm = smm + C22(t4b)*xn**(-1) * ( taucs*taucg**2*a1mm + taucs*
     &    taucg**2*a2mm )
      smm = smm + C22(t7)*mc**2*xn**(-1) * (  - taucs*taucg*a1mm - 2.D0
     &    *taucs*taugs*a1mm + taucg*taugs*a1mm + taucg*taugs*a2mm - 2.D0
     &    *taugs*a3mm - taugs**2*a1mm - taugs**2*a2mm )
      smm = smm + C22(t7)*xn**(-1) * ( taucs*taucg*taugs*a2mm - taucs*
     &    taucg**2*a1mm - 2.D0*taucs*taugs*a3mm - taucs*taugs**2*a1mm
     &     - taucs*taugs**2*a2mm - taucs**2*taucg*a1mm - 2.D0*taucs**2*
     &    taugs*a1mm )
      smm = smm + C22(t8)*mc**2*xn**(-1) * (  - taucs*taucg*a2mm )
      smm = smm + C22(t8)*xn**(-1) * (  - taucs**2*taucg*a1mm - 
     &    taucs**2*taucg*a2mm )
      smm = smm + C22(t9)*mc**2*xn**(-1) * (  - taucs*taucg*a1mm + 2.D0
     &    *taucs*taugs*a1mm + taucg*taugs*a1mm + taucg*taugs*a2mm + 
     &    taugs**2*a1mm + taugs**2*a2mm )
      smm = smm + C22(t9)*mc**4*xn**(-1) * ( 2.D0*taugs*a1mm )
      smm = smm + C22(t9)*xn**(-1) * ( taucs*taucg*taugs*a2mm - taucs*
     &    taucg**2*a1mm + taucs*taugs**2*a1mm + taucs*taugs**2*a2mm - 
     &    taucs**2*taucg*a1mm )
      smm = smm + C12(t1b)*mc**2*xn**(-1) * ( taucs*taugs*a1mm + taucs*
     &    taugs*a2mm + taucg*taugs*a1mm + taucg*taugs*a2mm - taugs**2*
     &    a1mm - taugs**2*a2mm )
      smm = smm + C12(t1b)*mc**2*xn * ( 2.D0*taucs*taugs*a1mm - 
     &    taugs**2*a1mm - taugs**2*a2mm )
      smm = smm + C12(t1b)*xn**(-1) * ( 2.D0*taucs*taucg*taugs*a1mm + 2.
     &    D0*taucs*taucg*taugs*a2mm + taucs**2*taugs*a1mm + taucs**2*
     &    taugs*a2mm )
      smm = smm + C12(t1b)*xn * ( taucs*taucg*taugs*a1mm + taucs*taucg*
     &    taugs*a2mm + taucs*taucg*a3mm - taucs*taucg**2*a2mm + 
     &    taucs**2*taucg*a1mm - taucs**2*taucg*a2mm + taucs**2*a3mm + 
     &    taucs**3*a1mm )
      smm = smm + C12(t2)*mc**2*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3mm + taucg**(-1)*taugs**2*a3mm - taugs**2*a2mm )
      smm = smm + C12(t2)*mc**2*xn * (  - taucs*taucg**(-1)*taugs*a3mm
     &     - taucs*taugs*a2mm - taucg**(-1)*taugs**2*a3mm )
      smm = smm + C12(t2)*mc**4*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a2mm )
      smm = smm + C12(t2)*mc**4*xn * ( taucg**(-1)*taugs**2*a2mm )
      smm = smm + C12(t2)*xn**(-1) * (  - taucs*taugs*a3mm - taucs*
     &    taugs**2*a2mm - taucs**2*taugs*a2mm - taucs**2*a3mm )
      smm = smm + C12(t2)*xn * (  - taucs*taucg*taugs*a2mm + 2.D0*taucs
     &    *taugs*a3mm + taucs**2*a3mm - taucg*taugs**2*a2mm + taugs**2*
     &    a3mm )
      smm = smm + C12(t2b)*xn * ( 2.D0*taucs*taucg*taugs*a2mm - 3.D0*
     &    taucs*taugs*a3mm - 3.D0*taucs*taugs**2*a1mm - taucs*taugs**2*
     &    a2mm + taucs**2*taucg*a2mm - 3.D0*taucs**2*taugs*a1mm - 
     &    taucs**2*a3mm - taucs**3*a1mm + taucg*taugs**2*a2mm - 2.D0*
     &    taugs**2*a3mm - taugs**3*a1mm - taugs**3*a2mm )
      smm = smm + C12(t3)*xn * ( 2.D0*taucs*taugs**2*a1mm - taucg*
     &    taugs**2*a2mm + 2.D0*taugs**2*a3mm + taugs**3*a1mm + taugs**3
     &    *a2mm )
      smm = smm + C12(t3b)*xn**(-1) * ( 2.D0*taucs*taugs**2*a1mm - 
     &    taucg*taugs**2*a2mm + 2.D0*taugs**2*a3mm + taugs**3*a1mm + 
     &    taugs**3*a2mm )
      smm = smm + C12(t4)*xn**(-1) * (  - taucg*taugs*a3mm + taucg**2*
     &    taugs*a2mm )
      smm = smm + C12(t4b)*mc**2*xn**(-1) * (  - 2.D0*taucg*taugs*a1mm
     &     )
      smm = smm + C12(t6)*mc**2*xn * (  - 3.D0*taugs*a3mm )
      smm = smm + C12(t6)*xn * ( 3.D0*taucs*taucg*a3mm + 2.D0*taucg*
     &    taugs*a3mm - 2.D0*taucg**2*taugs*a2mm )
      smm = smm + C12(t7)*mc**2*xn**(-1) * (  - taucg*taugs*a1mm )
      smm = smm + C12(t7)*xn**(-1) * (  - 3.D0*taucs*taucg*taugs*a1mm
     &     - 2.D0*taucs*taugs**2*a1mm - 2.D0*taucg*taugs*a3mm - taucg*
     &    taugs**2*a1mm + taucg**2*taugs*a2mm - 2.D0*taugs**2*a3mm - 
     &    taugs**3*a1mm - taugs**3*a2mm )
      smm = smm + C12(t8)*mc**2*xn**(-1) * (  - taucs*taugs*a2mm )
      smm = smm + C12(t8)*xn**(-1) * (  - taucs**2*taugs*a1mm - 
     &    taucs**2*taugs*a2mm )
      smm = smm + C12(t8b)*mc**2*xn**(-1) * ( taucs*taugs*a2mm )
      smm = smm + C12(t8b)*xn**(-1) * ( taucs**2*taugs*a1mm + taucs**2*
     &    taugs*a2mm )
      smm = smm + C12(t9)*mc**2*xn**(-1) * ( 2.D0*taucg*taugs*a1mm + 
     &    taugs**2*a1mm )
      smm = smm + C12(t9)*xn**(-1) * (  - taucs*taugs**2*a1mm - taucg*
     &    taugs*a3mm + taucg*taugs**2*a2mm + taucg**2*taugs*a2mm - 
     &    taugs**2*a3mm )
      smm = smm + C0(t2)*mc**2*xn**(-1) * ( taucs*taucg*a2mm - 2.D0*
     &    taucs*taugs*a1mm - taucs*taugs*a2mm - taucs*a3mm - taucs**2*
     &    a1mm - taugs*a3mm - taugs**2*a1mm - taugs**2*a2mm )
      smm = smm + C0(t2)*mc**4*xn**(-1) * (  - taugs*a2mm )
      smm = smm + C0(t4)*mc**2*xn**(-1) * (  - 2.D0*taucs*a3mm + taucg*
     &    taugs*a2mm - taugs*a3mm )
      smm = smm + C0(t4)*mc**4*xn**(-1) * ( 2.D0*taucg**(-1)*taugs*a3mm
     &     )
      smm = smm + C0(t4)*xn**(-1) * (  - taucg*taugs*a3mm + taucg**2*
     &    taugs*a2mm )
      smm = smm + C0(t4b)*mc**2*xn**(-1) * ( taucg*taugs*a1mm + taugs*
     &    a3mm )
      smm = smm + C0(t4b)*mc**4*xn**(-1) * (  - taugs*a1mm )
      smm = smm + C0(t6)*mc**2*xn * ( 2.D0*taucs*a3mm - taucg*taugs*
     &    a2mm + taugs*a3mm )
      smm = smm + C0(t6)*mc**4*xn * (  - 2.D0*taucg**(-1)*taugs*a3mm )
      smm = smm + C0(t6b)*mc**2*xn * ( 2.D0*taucs*taugs*a1mm - taugs*
     &    a3mm + taugs**2*a1mm )
      smm = smm + C0(t6b)*xn * (  - taucs*taucg*taugs*a1mm + taucs*
     &    taucg*taugs*a2mm + taucs*taucg*a3mm - taucg*taugs*a3mm - 
     &    taucg*taugs**2*a1mm - taucg*taugs**2*a2mm )
      smm = smm + C0(t8)*mc**2*xn**(-1) * ( taucg*taugs*a2mm )
      smm = smm + C0(t8)*xn**(-1) * ( taucs*taucg*taugs*a2mm - taucs*
     &    taugs*a3mm - 2.D0*taucs**2*taugs*a1mm )
      smm = smm + C0(t8b)*mc**2*xn**(-1) * (  - taucs*taucg*a2mm + 4.D0
     &    *taucs*taugs*a1mm + taucs*taugs*a2mm + taucs*a3mm + taucs**2*
     &    a1mm - taucg*taugs*a2mm + 3.D0*taugs*a3mm + 2.D0*taugs**2*
     &    a1mm + taugs**2*a2mm )
      smm = smm + C0(t8b)*mc**4*xn**(-1) * ( taugs*a2mm )
      smm = smm + C0(t8b)*xn**(-1) * ( taucs*taugs*a3mm )
      smm = smm + C0(t9)*mc**2*xn**(-1) * ( taucs*taugs*a1mm - taucg*
     &    taugs*a1mm - taucg*taugs*a2mm + taugs*a3mm + taugs**2*a1mm + 
     &    taugs**2*a2mm )
      smm = smm + C0(t9)*mc**4*xn**(-1) * ( taugs*a1mm )
      smm = smm + C1(t1b)*mc**2*xn**(-1) * (  - taucs*taugs*a1mm + 
     &    taucs*taugs*a2mm + taucg*taugs*a2mm + taugs*a3mm + taugs**2*
     &    a2mm )
      smm = smm + C1(t1b)*mc**2*xn * ( 2.D0*taucs*taugs*a1mm + taucg*
     &    taugs*a1mm - taugs*a3mm + 2.D0*taugs**2*a1mm + taugs**2*a2mm
     &     )
      smm = smm + C1(t1b)*xn**(-1) * (  - taucs*taucg*a3mm + 2.D0*taucs
     &    *taugs**2*a1mm + 2.D0*taucs*taugs**2*a2mm - taucs**2*a3mm )
      smm = smm + C1(t1b)*xn * (  - taucs*taucg*taugs*a2mm + taucs*
     &    taucg*a3mm + taucs*taugs*a3mm + taucs**2*taugs*a1mm + 
     &    taucs**2*a3mm )
      smm = smm + C1(t2)*mc**2*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3mm - 2.D0*taucs*taugs*a1mm - 3.D0*taucs*taugs*a2mm - 2.D0*
     &    taucs*a3mm - taucs**2*a1mm + taucg**(-1)*taugs**2*a3mm - 
     &    taucg*taugs*a2mm + taugs*a3mm - taugs**2*a1mm )
      smm = smm + C1(t2)*mc**2*xn * (  - taucs*taucg**(-1)*taugs*a3mm
     &     + taucs*taugs*a2mm + taucs*a3mm - taucg**(-1)*taugs**2*a3mm
     &     - taucg*taugs*a2mm )
      smm = smm + C1(t2)*mc**4*xn**(-1) * ( taucg**(-1)*taugs*a3mm + 
     &    taucg**(-1)*taugs**2*a2mm - taugs*a2mm )
      smm = smm + C1(t2)*mc**4*xn * (  - taucg**(-1)*taugs*a3mm - 
     &    taucg**(-1)*taugs**2*a2mm )
      smm = smm + C1(t2)*xn**(-1) * (  - 2.D0*taucs*taucg*taugs*a1mm - 
     &    2.D0*taucs*taucg*taugs*a2mm - taucs*taucg*a3mm - taucs**2*
     &    taucg*a1mm - taucs**2*taucg*a2mm - taucg*taugs**2*a1mm )
      smm = smm + C1(t2)*xn * (  - taucs*taucg*taugs*a2mm + taucs*taucg
     &    *a3mm + 2.D0*taucs*taugs*a3mm + taucs**2*a3mm + taucg*taugs*
     &    a3mm - taucg*taugs**2*a2mm - taucg**2*taugs*a2mm + taugs**2*
     &    a3mm )
      smm = smm + C1(t2b)*xn * ( taucs*taucg*taugs*a2mm - 2.D0*taucs*
     &    taugs*a3mm - taucs**2*a3mm + taucg*taugs**2*a2mm - taugs**2*
     &    a3mm )
      smm = smm + C1(t3)*xn * ( taucs*taugs*a3mm - taucg*taugs**2*a2mm
     &     + taugs**2*a3mm )
      smm = smm + C1(t3b)*xn**(-1) * ( taucs*taugs**2*a1mm )
      smm = smm + C1(t4)*mc**2*xn**(-1) * (  - 2.D0*taucs*a3mm + 3.D0*
     &    taucg*taugs*a2mm - 3.D0*taugs*a3mm )
      smm = smm + C1(t4)*mc**4*xn**(-1) * ( 2.D0*taucg**(-1)*taugs*a3mm
     &     )
      smm = smm + C1(t4)*xn**(-1) * (  - taucg*taugs*a3mm + taucg**2*
     &    taugs*a2mm )
      smm = smm + C1(t4b)*mc**2*xn**(-1) * (  - 2.D0*taucg*taugs*a1mm )
      smm = smm + C1(t4b)*xn**(-1) * ( taucg*taugs*a3mm - taucg**2*
     &    taugs*a1mm - taucg**2*taugs*a2mm )
      smm = smm + C1(t6)*mc**2*xn * ( 6.D0*taucs*a3mm - 3.D0*taucg*
     &    taugs*a2mm + 2.D0*taugs*a3mm )
      smm = smm + C1(t6)*mc**4*xn * (  - 6.D0*taucg**(-1)*taugs*a3mm )
      smm = smm + C1(t6)*xn * ( taucs*taucg*a3mm + taucg*taugs*a3mm - 
     &    taucg**2*taugs*a2mm )
      smm = smm + C1(t6b)*mc**2*xn * (  - taucg*taugs*a1mm )
      smm = smm + C1(t6b)*xn * (  - taucs*taucg*taugs*a1mm - taucg*
     &    taugs*a3mm - taucg*taugs**2*a1mm - taucg*taugs**2*a2mm )
      smm = smm + C1(t7)*xn**(-1) * (  - taucs*taucg*taugs*a1mm - taucs
     &    *taugs**2*a1mm )
      smm = smm + C1(t8)*mc**2*xn**(-1) * (  - taucs*taugs*a2mm )
      smm = smm + C1(t8)*xn**(-1) * (  - taucs*taugs*a3mm - taucs*
     &    taugs**2*a1mm - taucs*taugs**2*a2mm - 2.D0*taucs**2*taugs*
     &    a1mm + taucs**2*a3mm )
      smm = smm + C1(t8b)*mc**2*xn**(-1) * ( 4.D0*taucs*taugs*a1mm + 2.D
     &    0*taucs*taugs*a2mm + taucs*a3mm + taucs**2*a1mm - taucg*taugs
     &    *a2mm + 3.D0*taugs*a3mm + 2.D0*taugs**2*a1mm + 2.D0*taugs**2*
     &    a2mm )
      smm = smm + C1(t8b)*mc**4*xn**(-1) * ( taugs*a2mm )
      smm = smm + C1(t8b)*xn**(-1) * ( taucs**2*taucg*a1mm + taucs**2*
     &    taucg*a2mm - taucs**2*a3mm )
      smm = smm + C1(t9)*mc**2*xn**(-1) * ( 2.D0*taucg*taugs*a1mm )
      smm = smm + C1(t9)*xn**(-1) * ( taucs*taucg*taugs*a1mm - taucg*
     &    taugs*a3mm + taucg*taugs**2*a1mm + taucg*taugs**2*a2mm + 
     &    taucg**2*taugs*a1mm + taucg**2*taugs*a2mm - taugs**2*a3mm )
      smm = smm + C2(t1b)*mc**2*xn**(-1) * ( 2.D0*taucs*taugs*a1mm + 
     &    taucs*taugs*a2mm - taucs*a3mm - taucg*taugs*a1mm - taugs*a3mm
     &     - taugs**2*a1mm - taugs**2*a2mm )
      smm = smm + C2(t1b)*mc**2*xn * ( taucs*a3mm + taucg*taugs*a1mm + 
     &    taugs*a3mm - taugs**2*a1mm - taugs**2*a2mm )
      smm = smm + C2(t1b)*mc**4*xn**(-1) * (  - taugs*a1mm + 2.D0*taugs
     &    *a2mm )
      smm = smm + C2(t1b)*mc**4*xn * ( 3.D0*taugs*a1mm )
      smm = smm + C2(t1b)*xn**(-1) * ( 2.D0*taucs*taucg*taugs*a1mm + 2.D
     &    0*taucs*taucg*taugs*a2mm + taucs*taucg*a3mm + 2.D0*taucs**2*
     &    taugs*a1mm + 2.D0*taucs**2*taugs*a2mm + taucs**2*a3mm )
      smm = smm + C2(t1b)*xn * (  - taucs*taucg*a3mm - taucs**2*taugs*
     &    a1mm - taucs**2*taugs*a2mm - taucs**2*a3mm )
      smm = smm + C2(t2)*mc**2*xn**(-1) * (  - taugs**2*a2mm )
      smm = smm + C2(t2)*mc**2*xn * (  - taucs*taugs*a2mm )
      smm = smm + C2(t2)*mc**4*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a2mm )
      smm = smm + C2(t2)*mc**4*xn * ( taucg**(-1)*taugs**2*a2mm )
      smm = smm + C2(t2)*xn**(-1) * ( taucs*taucg*taugs*a2mm - 2.D0*
     &    taucs*taugs*a3mm - 3.D0*taucs*taugs**2*a1mm - 2.D0*taucs*
     &    taugs**2*a2mm + taucs**2*taucg*a2mm - 3.D0*taucs**2*taugs*
     &    a1mm - taucs**2*taugs*a2mm - taucs**2*a3mm - taucs**3*a1mm - 
     &    taugs**2*a3mm - taugs**3*a1mm - taugs**3*a2mm )
      smm = smm + C2(t2b)*mc**2*xn * ( taucs*taugs*a1mm - taucs*a3mm + 
     &    taucg*taugs*a2mm - taugs*a3mm + taugs**2*a1mm )
      smm = smm + C2(t2b)*xn * ( taucs*taucg*taugs*a2mm + taucs*taugs*
     &    a3mm + taucs*taugs**2*a2mm + taucs**2*taugs*a2mm + taucs**2*
     &    a3mm - taucg*taugs*a3mm + taucg**2*taugs*a2mm )
      smm = smm + C2(t3)*mc**2*xn * (  - taugs**2*a1mm )
      smm = smm + C2(t3)*xn * ( taucs*taugs**2*a1mm )
      smm = smm + C2(t3b)*mc**2*xn**(-1) * (  - taugs**2*a2mm )
      smm = smm + C2(t3b)*xn**(-1) * ( taucs*taugs*a3mm - taucg*
     &    taugs**2*a2mm + taugs**2*a3mm )
      smm = smm + C2(t4)*xn**(-1) * (  - taucg*taugs*a3mm + taucg**2*
     &    taugs*a2mm )
      smm = smm + C2(t4b)*mc**2*xn**(-1) * (  - taucg*taugs*a1mm - 
     &    taucg*taugs*a2mm + taugs*a3mm )
      smm = smm + C2(t4b)*mc**4*xn**(-1) * (  - 2.D0*taugs*a1mm )
      smm = smm + C2(t6)*mc**2*xn * (  - taugs*a3mm )
      smm = smm + C2(t6)*xn * ( taucs*taucg*a3mm + taucg*taugs*a3mm - 
     &    taucg**2*taugs*a2mm )
      smm = smm + C2(t6b)*mc**2*xn * (  - 2.D0*taucs*taugs*a1mm - taucg
     &    *taugs*a1mm - 3.D0*taugs*a3mm - 2.D0*taugs**2*a1mm - 2.D0*
     &    taugs**2*a2mm )
      smm = smm + C2(t6b)*mc**4*xn * (  - 3.D0*taugs*a1mm )
      smm = smm + C2(t6b)*xn * ( taucs*taucg*taugs*a1mm + taucs*taucg*
     &    taugs*a2mm + 2.D0*taucs*taucg*a3mm - taucs*taucg**2*a2mm + 
     &    taucs**2*taucg*a1mm )
      smm = smm + C2(t7)*mc**2*xn**(-1) * ( taucg*taugs*a1mm + 2.D0*
     &    taucg*taugs*a2mm - taugs*a3mm + taugs**2*a1mm + taugs**2*a2mm
     &     )
      smm = smm + C2(t7)*mc**4*xn**(-1) * ( taugs*a1mm )
      smm = smm + C2(t7)*xn**(-1) * ( taucs*taucg*taugs*a2mm - taucs*
     &    taugs*a3mm - taucs**2*taugs*a1mm - taucg*taugs*a3mm + taucg*
     &    taugs**2*a2mm + taucg**2*taugs*a2mm - taugs**2*a3mm )
      smm = smm + C2(t8)*mc**2*xn**(-1) * (  - 3.D0*taucs*taugs*a1mm - 
     &    taucs*taugs*a2mm + taucs*a3mm - 2.D0*taugs*a3mm - 2.D0*
     &    taugs**2*a1mm - 2.D0*taugs**2*a2mm )
      smm = smm + C2(t8)*mc**4*xn**(-1) * (  - 2.D0*taugs*a2mm )
      smm = smm + C2(t8)*xn**(-1) * (  - 2.D0*taucs**2*taugs*a1mm - 2.D0
     &    *taucs**2*taugs*a2mm - taucs**2*a3mm )
      smm = smm + C2(t8b)*mc**2*xn**(-1) * ( taucs*taugs*a2mm )
      smm = smm + C2(t8b)*xn**(-1) * (  - taucs*taucg*taugs*a2mm + 2.D0
     &    *taucs*taugs*a3mm + taucs*taugs**2*a1mm + taucs*taugs**2*a2mm
     &     - taucs**2*taucg*a2mm + 3.D0*taucs**2*taugs*a1mm + taucs**2*
     &    taugs*a2mm + taucs**2*a3mm + taucs**3*a1mm )
      smm = smm + C2(t9)*mc**2*xn**(-1) * ( 3.D0*taucs*taugs*a1mm + 
     &    taucg*taugs*a1mm + taucg*taugs*a2mm + taugs*a3mm + 3.D0*
     &    taugs**2*a1mm + 2.D0*taugs**2*a2mm )
      smm = smm + C2(t9)*mc**4*xn**(-1) * ( 2.D0*taugs*a1mm )
      smm = smm + C2(t9)*xn**(-1) * ( taucs*taucg*taugs*a1mm + taucs*
     &    taucg*taugs*a2mm + 2.D0*taucs*taugs**2*a1mm + taucs*taugs**2*
     &    a2mm + taucs**2*taugs*a1mm + taucg*taugs*a3mm + taucg*
     &    taugs**2*a1mm + taucg*taugs**2*a2mm + taugs**2*a3mm + 
     &    taugs**3*a1mm + taugs**3*a2mm )
      smm = smm + D0(b2)*mc**2*xn**(-1) * ( taucs*taucg*taugs*a2mm - 
     &    taucs*taugs*a3mm - taucg*taugs**2*a1mm - taucg*taugs**2*a2mm
     &     )
      smm = smm + D0(b2)*mc**4*xn**(-1) * ( taucs*taugs*a1mm + taucg*
     &    taugs*a2mm - taugs*a3mm - taugs**2*a2mm )
      smm = smm + D1(b1)*xn**(-1) * (  - taucs**2*taugs**2*a1mm )
      smm = smm + D1(b2)*mc**2*xn**(-1) * ( taucs*taucg*taugs*a1mm - 
     &    taucg*taugs**2*a1mm - 2.D0*taucg*taugs**2*a2mm + taugs**2*
     &    a3mm )
      smm = smm + D1(b2)*xn**(-1) * (  - taucs*taucg*taugs*a3mm + taucs
     &    *taucg**2*taugs*a2mm )
      smm = smm + D1(b3)*mc**2*xn * (  - taugs**2*a3mm )
      smm = smm + D2(b1)*mc**2*xn**(-1) * ( taucs*taugs**2*a2mm + taucg
     &    *taugs**2*a2mm )
      smm = smm + D2(b1)*xn**(-1) * ( taucs*taucg*taugs**2*a2mm + taucs
     &    *taugs**3*a1mm + taucs*taugs**3*a2mm - taucs**2*taugs*a3mm + 
     &    taucs**2*taugs**2*a1mm )
      smm = smm + D2(b2)*mc**2*xn**(-1) * (  - taucs*taucg*taugs*a1mm
     &     - taugs**2*a3mm )
      smm = smm + D2(b2)*mc**4*xn**(-1) * ( taucs*taugs*a1mm + taucg*
     &    taugs*a2mm - taugs*a3mm - taugs**2*a2mm )
      smm = smm + D2(b2)*xn**(-1) * ( taucs*taucg*taugs*a3mm )
      smm = smm + D2(b3)*mc**2*xn * ( 2.D0*taucs*taugs**2*a1mm + taucg*
     &    taugs**2*a1mm + taugs**3*a1mm )
      smm = smm + D2(b3)*xn * ( taucs*taucg*taugs*a3mm + taucs*taucg*
     &    taugs**2*a2mm )
      smm = smm + D3(b1)*mc**2*xn**(-1) * ( taucs*taucg*taugs*a2mm + 2.D
     &    0*taucs*taugs**2*a1mm + taucs**2*taugs*a1mm + taucg**2*taugs*
     &    a2mm )
      smm = smm + D3(b1)*xn**(-1) * ( taucs*taucg*taugs**2*a1mm + taucs
     &    *taucg*taugs**2*a2mm + taucs*taucg**2*taugs*a2mm + taucs**2*
     &    taucg*taugs*a2mm + taucs**2*taugs*a3mm + taucs**2*taugs**2*
     &    a1mm + taucs**2*taugs**2*a2mm )
      smm = smm + D3(b2)*mc**2*xn**(-1) * ( taucs*taucg*taugs*a2mm - 
     &    taucs*taugs*a3mm - taucs*taugs**2*a2mm - taugs**3*a2mm )
      smm = smm + D3(b2)*xn**(-1) * ( taucs*taucg*taugs**2*a2mm + 
     &    taucs**2*taucg*taugs*a2mm )
      smm = smm + D3(b3)*mc**2*xn * ( taucs*taucg*taugs*a1mm + taucs*
     &    taugs*a3mm + taucs*taugs**2*a1mm + 2.D0*taucs**2*taugs*a1mm
     &     + taucg*taugs**2*a1mm + taugs**2*a3mm )
      smm = smm + D3(b3)*mc**4*xn * ( taugs**2*a1mm )
      smm = smm + D3(b3)*xn * (  - taucs*taucg*taugs*a3mm + taucs*
     &    taucg**2*taugs*a2mm + taucs**2*taucg*taugs*a2mm )
      smm = smm - epinv*taucs*b0*a3mm + epinv*taucg*taugs*b0*a2mm - 
     &    epinv*taugs*b0*a3mm

      smp=  + mc**2*xn**(-1) * (  - 3.D0/4.D0*epinv*taucg**(-2)*taugs*
     &    a3mp + 3.D0/4.D0*epinv*taucg**(-1)*taugs*a1mp + 2.D0*C00s(t1)
     &    *epinv*taucg**(-2)*taugs*a3mp - 2.D0*C00s(t1)*epinv*
     &    taucg**(-1)*taugs*a1mp + 2.D0*C00s(t1)*taucs*taucg**(-2)*
     &    taugs*a1mp + 2.D0*C00s(t1)*taucg**(-2)*taugs*a3mp + 2.D0*
     &    C00s(t1)*taucg**(-2)*taugs**2*a1mp + 2.D0*C00s(t1)*
     &    taucg**(-2)*taugs**2*a2mp + 2.D0*C00s(t1)*taucg**(-1)*taugs*
     &    a1mp - 2.D0*C00s(t1)*taucg**(-1)*taugs*a2mp - 2.D0*C00s(t1b)*
     &    epinv*taucs*taucg**(-2)*taugs*a1mp - 2.D0*C00s(t1b)*epinv*
     &    taucg**(-2)*taugs*a3mp - 2.D0*C00s(t1b)*epinv*taucg**(-1)*
     &    taugs*a1mp + 2.D0*C00s(t1b)*epinv*taucg**(-1)*taugs*a2mp + 2.D
     &    0*C00s(t2)*epinv*taucs*taucg**(-2)*taugs*a1mp + 2.D0*C00s(t2)
     &    *epinv*taucg**(-2)*taugs*a3mp + 2.D0*C00s(t2)*epinv*
     &    taucg**(-2)*taugs**2*a1mp - 2.D0*C00s(t2)*epinv*taucg**(-1)*
     &    taugs*a2mp + 2.D0*C00s(t3)*epinv*taucg**(-2)*taugs*a3mp - 2.D0
     &    *C00s(t3)*epinv*taucg**(-1)*taugs*a1mp )
      smp = smp + mc**2*xn**(-1) * (  - 2.D0*C00s(t3)*taucg**(-2)*taugs
     &    *a3mp + 2.D0*C00s(t3b)*epinv*taucg**(-2)*taugs**2*a1mp - 2.D0
     &    *C00s(t4b)*epinv*taucg**(-1)*taugs*a1mp - 2.D0*C00s(t7)*epinv
     &    *taucg**(-2)*taugs**2*a1mp + 2.D0*C00s(t7)*epinv*taucg**(-1)*
     &    taugs*a1mp + 2.D0*C00s(t8)*epinv*taucs*taucg**(-2)*taugs*a1mp
     &     + 2.D0*C00s(t8)*epinv*taucg**(-2)*taugs*a3mp - 2.D0*C00s(t8)
     &    *epinv*taucg**(-1)*taugs*a2mp - 2.D0*C00s(t8b)*epinv*taucs*
     &    taucg**(-2)*taugs*a1mp - 2.D0*C00s(t8b)*epinv*taucg**(-2)*
     &    taugs*a3mp + 2.D0*C00s(t8b)*epinv*taucg**(-1)*taugs*a2mp - 2.D
     &    0*C00s(t9)*epinv*taucg**(-2)*taugs**2*a1mp + 2.D0*C00s(t9)*
     &    epinv*taucg**(-1)*taugs*a1mp + 2.D0*C00f(t1)*taucg**(-2)*
     &    taugs*a3mp - 2.D0*C00f(t1)*taucg**(-1)*taugs*a1mp - 2.D0*
     &    C00f(t1b)*taucs*taucg**(-2)*taugs*a1mp - 2.D0*C00f(t1b)*
     &    taucg**(-2)*taugs*a3mp - 2.D0*C00f(t1b)*taucg**(-1)*taugs*
     &    a1mp + 2.D0*C00f(t1b)*taucg**(-1)*taugs*a2mp + 2.D0*C00f(t2)*
     &    taucs*taucg**(-2)*taugs*a1mp )
      smp = smp + mc**2*xn**(-1) * ( 2.D0*C00f(t2)*taucg**(-2)*taugs*
     &    a3mp + 2.D0*C00f(t2)*taucg**(-2)*taugs**2*a1mp - 2.D0*C00f(t2
     &    )*taucg**(-1)*taugs*a2mp + 2.D0*C00f(t3)*taucg**(-2)*taugs*
     &    a3mp - 2.D0*C00f(t3)*taucg**(-1)*taugs*a1mp + 2.D0*C00f(t3b)*
     &    taucg**(-2)*taugs**2*a1mp - 2.D0*C00f(t4b)*taucg**(-1)*taugs*
     &    a1mp - 2.D0*C00f(t7)*taucg**(-2)*taugs**2*a1mp + 2.D0*C00f(t7
     &    )*taucg**(-1)*taugs*a1mp + 2.D0*C00f(t8)*taucs*taucg**(-2)*
     &    taugs*a1mp + 2.D0*C00f(t8)*taucg**(-2)*taugs*a3mp - 2.D0*
     &    C00f(t8)*taucg**(-1)*taugs*a2mp - 2.D0*C00f(t8b)*taucs*
     &    taucg**(-2)*taugs*a1mp - 2.D0*C00f(t8b)*taucg**(-2)*taugs*
     &    a3mp + 2.D0*C00f(t8b)*taucg**(-1)*taugs*a2mp - 2.D0*C00f(t9)*
     &    taucg**(-2)*taugs**2*a1mp + 2.D0*C00f(t9)*taucg**(-1)*taugs*
     &    a1mp )
      smp = smp + mc**2*xn * ( 3.D0/4.D0*epinv*taucg**(-2)*taugs*a3mp
     &     - 3.D0/4.D0*epinv*taucg**(-1)*taugs*a1mp - 1.D0/6.D0*
     &    taucg**(-2)*taugs*a3mp + 1.D0/6.D0*taucg**(-1)*taugs*a1mp - 2.
     &    D0*C00s(t1)*epinv*taucg**(-2)*taugs*a3mp + 2.D0*C00s(t1)*
     &    epinv*taucg**(-1)*taugs*a1mp - 2.D0*C00s(t1)*taucs*
     &    taucg**(-2)*taugs*a1mp - 2.D0*C00s(t1)*taucg**(-2)*taugs*a3mp
     &     - 2.D0*C00s(t1)*taucg**(-2)*taugs**2*a1mp - 2.D0*C00s(t1)*
     &    taucg**(-2)*taugs**2*a2mp - 2.D0*C00s(t1)*taucg**(-1)*taugs*
     &    a1mp + 2.D0*C00s(t1)*taucg**(-1)*taugs*a2mp + 4.D0*C00s(t1b)*
     &    epinv*taucs*taucg**(-2)*taugs*a1mp + 4.D0*C00s(t1b)*epinv*
     &    taucg**(-2)*taugs*a3mp + 2.D0*C00s(t1b)*epinv*taucg**(-2)*
     &    taugs**2*a1mp + 2.D0*C00s(t1b)*epinv*taucg**(-2)*taugs**2*
     &    a2mp - 6.D0*C00s(t1b)*epinv*taucg**(-1)*taugs*a2mp - 4.D0*
     &    C00s(t2b)*epinv*taucs*taucg**(-2)*taugs*a1mp - 4.D0*C00s(t2b)
     &    *epinv*taucg**(-2)*taugs*a3mp - 4.D0*C00s(t2b)*epinv*
     &    taucg**(-2)*taugs**2*a1mp )
      smp = smp + mc**2*xn * (  - 2.D0*C00s(t2b)*epinv*taucg**(-2)*
     &    taugs**2*a2mp + 6.D0*C00s(t2b)*epinv*taucg**(-1)*taugs*a2mp
     &     - 10.D0*C00s(t3)*epinv*taucg**(-2)*taugs*a3mp + 2.D0*C00s(t3
     &    )*epinv*taucg**(-2)*taugs**2*a1mp + 10.D0*C00s(t3)*epinv*
     &    taucg**(-1)*taugs*a1mp + 4.D0*C00s(t3)*taucg**(-2)*taugs*a3mp
     &     - 6.D0*C00s(t3)*taucg**(-1)*taugs*a1mp - 2.D0*C00f(t1)*
     &    taucg**(-2)*taugs*a3mp + 2.D0*C00f(t1)*taucg**(-1)*taugs*a1mp
     &     + 4.D0*C00f(t1b)*taucs*taucg**(-2)*taugs*a1mp + 4.D0*C00f(
     &    t1b)*taucg**(-2)*taugs*a3mp + 2.D0*C00f(t1b)*taucg**(-2)*
     &    taugs**2*a1mp + 2.D0*C00f(t1b)*taucg**(-2)*taugs**2*a2mp - 6.D
     &    0*C00f(t1b)*taucg**(-1)*taugs*a2mp - 4.D0*C00f(t2b)*taucs*
     &    taucg**(-2)*taugs*a1mp - 4.D0*C00f(t2b)*taucg**(-2)*taugs*
     &    a3mp - 4.D0*C00f(t2b)*taucg**(-2)*taugs**2*a1mp - 2.D0*C00f(
     &    t2b)*taucg**(-2)*taugs**2*a2mp + 6.D0*C00f(t2b)*taucg**(-1)*
     &    taugs*a2mp - 10.D0*C00f(t3)*taucg**(-2)*taugs*a3mp + 2.D0*
     &    C00f(t3)*taucg**(-2)*taugs**2*a1mp )
      smp = smp + mc**2*xn * ( 10.D0*C00f(t3)*taucg**(-1)*taugs*a1mp )
      smp = smp + mc**2 * ( epinv*taucg**(-2)*taugs*b0*a3mp - epinv*
     &    taucg**(-1)*taugs*b0*a1mp )
      smp = smp + xn**(-1) * ( 3.D0/4.D0*epinv*taucs*taucg**(-1)*a3mp
     &     + 3.D0/4.D0*epinv*taugs*a1mp + 3.D0/4.D0*epinv*taugs*a2mp + 
     &    3.D0/4.D0*epinv*a3mp - 2.D0*C00s(t1)*epinv*taucs*taucg**(-1)*
     &    a3mp - 2.D0*C00s(t1)*epinv*taugs*a1mp - 2.D0*C00s(t1)*epinv*
     &    taugs*a2mp - 2.D0*C00s(t1)*epinv*a3mp - 2.D0*C00s(t1)*taucs*
     &    taucg**(-1)*taugs*a1mp - 2.D0*C00s(t1)*taucs*taucg**(-1)*
     &    taugs*a2mp - 2.D0*C00s(t1)*taucs*taucg**(-1)*a3mp - 2.D0*
     &    C00s(t1)*taucs*a1mp + 2.D0*C00s(t1)*taucs*a2mp - 2.D0*C00s(t1
     &    )*taucs**2*taucg**(-1)*a1mp + 2.D0*C00s(t1)*taucg*a2mp - 2.D0
     &    *C00s(t1)*a3mp + 2.D0*C00s(t1b)*epinv*taucs*taucg**(-1)*a3mp
     &     + 4.D0*C00s(t1b)*epinv*taucs*a1mp - 2.D0*C00s(t1b)*epinv*
     &    taucs*a2mp + 2.D0*C00s(t1b)*epinv*taucs**2*taucg**(-1)*a1mp
     &     - 4.D0*C00s(t1b)*epinv*taucg*a2mp + 2.D0*C00s(t1b)*epinv*
     &    taugs*a1mp + 2.D0*C00s(t1b)*epinv*taugs*a2mp + 4.D0*C00s(t1b)
     &    *epinv*a3mp - 2.D0*C00s(t2)*epinv*taucs*taucg**(-1)*taugs*
     &    a1mp )
      smp = smp + xn**(-1) * (  - 2.D0*C00s(t2)*epinv*taucs*taucg**(-1)
     &    *a3mp + 2.D0*C00s(t2)*epinv*taucs*a2mp - 2.D0*C00s(t2)*epinv*
     &    taucs**2*taucg**(-1)*a1mp - 2.D0*C00s(t3)*epinv*taucs*
     &    taucg**(-1)*a3mp - 2.D0*C00s(t3)*epinv*taugs*a1mp - 2.D0*
     &    C00s(t3)*epinv*taugs*a2mp - 2.D0*C00s(t3)*epinv*a3mp + 2.D0*
     &    C00s(t3)*taucs*taucg**(-1)*a3mp - 2.D0*C00s(t3b)*epinv*taucs*
     &    taucg**(-1)*taugs*a1mp - 2.D0*C00s(t3b)*epinv*taugs*a2mp + 2.D
     &    0*C00s(t4)*taucg**(-1)*taugs*a3mp - 2.D0*C00s(t4)*taugs*a2mp
     &     + 4.D0*C00s(t4b)*epinv*taucs*a1mp - 4.D0*C00s(t4b)*epinv*
     &    taucg*a2mp + 2.D0*C00s(t4b)*epinv*taugs*a1mp + 2.D0*C00s(t4b)
     &    *epinv*taugs*a2mp + 6.D0*C00s(t4b)*epinv*a3mp + 2.D0*C00s(t7)
     &    *epinv*taucs*taucg**(-1)*taugs*a1mp - 4.D0*C00s(t7)*epinv*
     &    taucs*a1mp + 4.D0*C00s(t7)*epinv*taucg*a2mp - 2.D0*C00s(t7)*
     &    epinv*taugs*a1mp - 4.D0*C00s(t7)*epinv*a3mp - 2.D0*C00s(t8)*
     &    epinv*taucs*taucg**(-1)*a3mp + 2.D0*C00s(t8)*epinv*taucs*a2mp
     &     - 2.D0*C00s(t8)*epinv*taucs**2*taucg**(-1)*a1mp )
      smp = smp + xn**(-1) * ( 2.D0*C00s(t8b)*epinv*taucs*taucg**(-1)*
     &    a3mp - 2.D0*C00s(t8b)*epinv*taucs*a2mp + 2.D0*C00s(t8b)*epinv
     &    *taucs**2*taucg**(-1)*a1mp + 2.D0*C00s(t9)*epinv*taucs*
     &    taucg**(-1)*taugs*a1mp - 4.D0*C00s(t9)*epinv*taucs*a1mp + 4.D0
     &    *C00s(t9)*epinv*taucg*a2mp - 2.D0*C00s(t9)*epinv*taugs*a1mp
     &     - 2.D0*C00s(t9)*epinv*taugs*a2mp - 6.D0*C00s(t9)*epinv*a3mp
     &     - 2.D0*C00f(t1)*taucs*taucg**(-1)*a3mp - 2.D0*C00f(t1)*taugs
     &    *a1mp - 2.D0*C00f(t1)*taugs*a2mp - 2.D0*C00f(t1)*a3mp + 2.D0*
     &    C00f(t1b)*taucs*taucg**(-1)*a3mp + 4.D0*C00f(t1b)*taucs*a1mp
     &     - 2.D0*C00f(t1b)*taucs*a2mp + 2.D0*C00f(t1b)*taucs**2*
     &    taucg**(-1)*a1mp - 4.D0*C00f(t1b)*taucg*a2mp + 2.D0*C00f(t1b)
     &    *taugs*a1mp + 2.D0*C00f(t1b)*taugs*a2mp + 4.D0*C00f(t1b)*a3mp
     &     - 2.D0*C00f(t2)*taucs*taucg**(-1)*taugs*a1mp - 2.D0*C00f(t2)
     &    *taucs*taucg**(-1)*a3mp + 2.D0*C00f(t2)*taucs*a2mp - 2.D0*
     &    C00f(t2)*taucs**2*taucg**(-1)*a1mp - 2.D0*C00f(t3)*taucs*
     &    taucg**(-1)*a3mp )
      smp = smp + xn**(-1) * (  - 2.D0*C00f(t3)*taugs*a1mp - 2.D0*C00f(
     &    t3)*taugs*a2mp - 2.D0*C00f(t3)*a3mp - 2.D0*C00f(t3b)*taucs*
     &    taucg**(-1)*taugs*a1mp - 2.D0*C00f(t3b)*taugs*a2mp + 4.D0*
     &    C00f(t4b)*taucs*a1mp - 4.D0*C00f(t4b)*taucg*a2mp + 2.D0*C00f(
     &    t4b)*taugs*a1mp + 2.D0*C00f(t4b)*taugs*a2mp + 6.D0*C00f(t4b)*
     &    a3mp + 2.D0*C00f(t7)*taucs*taucg**(-1)*taugs*a1mp - 4.D0*
     &    C00f(t7)*taucs*a1mp + 4.D0*C00f(t7)*taucg*a2mp - 2.D0*C00f(t7
     &    )*taugs*a1mp - 4.D0*C00f(t7)*a3mp - 2.D0*C00f(t8)*taucs*
     &    taucg**(-1)*a3mp + 2.D0*C00f(t8)*taucs*a2mp - 2.D0*C00f(t8)*
     &    taucs**2*taucg**(-1)*a1mp + 2.D0*C00f(t8b)*taucs*taucg**(-1)*
     &    a3mp - 2.D0*C00f(t8b)*taucs*a2mp + 2.D0*C00f(t8b)*taucs**2*
     &    taucg**(-1)*a1mp + 2.D0*C00f(t9)*taucs*taucg**(-1)*taugs*a1mp
     &     - 4.D0*C00f(t9)*taucs*a1mp + 4.D0*C00f(t9)*taucg*a2mp - 2.D0
     &    *C00f(t9)*taugs*a1mp - 2.D0*C00f(t9)*taugs*a2mp - 6.D0*C00f(
     &    t9)*a3mp )
      smp = smp + xn * (  - 3.D0/4.D0*epinv*taucs*taucg**(-1)*a3mp - 3.D
     &    0/4.D0*epinv*taugs*a1mp - 3.D0/4.D0*epinv*taugs*a2mp - 3.D0/4.
     &    D0*epinv*a3mp + 1.D0/6.D0*taucs*taucg**(-1)*a3mp + 1.D0/6.D0*
     &    taugs*a1mp + 1.D0/6.D0*taugs*a2mp + 1.D0/6.D0*a3mp + 2.D0*
     &    C00s(t1)*epinv*taucs*taucg**(-1)*a3mp + 2.D0*C00s(t1)*epinv*
     &    taugs*a1mp + 2.D0*C00s(t1)*epinv*taugs*a2mp + 2.D0*C00s(t1)*
     &    epinv*a3mp + 2.D0*C00s(t1)*taucs*taucg**(-1)*taugs*a1mp + 2.D0
     &    *C00s(t1)*taucs*taucg**(-1)*taugs*a2mp + 2.D0*C00s(t1)*taucs*
     &    taucg**(-1)*a3mp + 2.D0*C00s(t1)*taucs*a1mp - 2.D0*C00s(t1)*
     &    taucs*a2mp + 2.D0*C00s(t1)*taucs**2*taucg**(-1)*a1mp - 2.D0*
     &    C00s(t1)*taucg*a2mp + 2.D0*C00s(t1)*a3mp - 2.D0*C00s(t1b)*
     &    epinv*taucs*taucg**(-1)*taugs*a1mp - 2.D0*C00s(t1b)*epinv*
     &    taucs*taucg**(-1)*taugs*a2mp - 4.D0*C00s(t1b)*epinv*taucs*
     &    taucg**(-1)*a3mp - 2.D0*C00s(t1b)*epinv*taucs*a1mp + 4.D0*
     &    C00s(t1b)*epinv*taucs*a2mp - 4.D0*C00s(t1b)*epinv*taucs**2*
     &    taucg**(-1)*a1mp )
      smp = smp + xn * ( 2.D0*C00s(t1b)*epinv*taucg*a2mp - 2.D0*C00s(
     &    t1b)*epinv*taugs*a1mp - 2.D0*C00s(t1b)*epinv*taugs*a2mp - 2.D0
     &    *C00s(t1b)*epinv*a3mp + 4.D0*C00s(t2b)*epinv*taucs*
     &    taucg**(-1)*taugs*a1mp + 2.D0*C00s(t2b)*epinv*taucs*
     &    taucg**(-1)*taugs*a2mp + 4.D0*C00s(t2b)*epinv*taucs*
     &    taucg**(-1)*a3mp - 4.D0*C00s(t2b)*epinv*taucs*a2mp + 4.D0*
     &    C00s(t2b)*epinv*taucs**2*taucg**(-1)*a1mp + 2.D0*C00s(t2b)*
     &    epinv*taugs*a2mp - 2.D0*C00s(t3)*epinv*taucs*taucg**(-1)*
     &    taugs*a1mp + 10.D0*C00s(t3)*epinv*taucs*taucg**(-1)*a3mp + 10.
     &    D0*C00s(t3)*epinv*taugs*a1mp + 8.D0*C00s(t3)*epinv*taugs*a2mp
     &     + 10.D0*C00s(t3)*epinv*a3mp - 4.D0*C00s(t3)*taucs*
     &    taucg**(-1)*a3mp - 6.D0*C00s(t3)*taugs*a1mp - 6.D0*C00s(t3)*
     &    taugs*a2mp - 6.D0*C00s(t3)*a3mp + 2.D0*C00s(t6)*taucg**(-1)*
     &    taugs*a3mp - 2.D0*C00s(t6)*taugs*a2mp + 2.D0*C00s(t6b)*epinv*
     &    taucs*a1mp - 2.D0*C00s(t6b)*epinv*taucg*a2mp + 2.D0*C00s(t6b)
     &    *epinv*taugs*a1mp )
      smp = smp + xn * ( 2.D0*C00s(t6b)*epinv*taugs*a2mp + 2.D0*C00s(
     &    t6b)*epinv*a3mp + 2.D0*C00f(t1)*taucs*taucg**(-1)*a3mp + 2.D0
     &    *C00f(t1)*taugs*a1mp + 2.D0*C00f(t1)*taugs*a2mp + 2.D0*C00f(
     &    t1)*a3mp - 2.D0*C00f(t1b)*taucs*taucg**(-1)*taugs*a1mp - 2.D0
     &    *C00f(t1b)*taucs*taucg**(-1)*taugs*a2mp - 4.D0*C00f(t1b)*
     &    taucs*taucg**(-1)*a3mp - 2.D0*C00f(t1b)*taucs*a1mp + 4.D0*
     &    C00f(t1b)*taucs*a2mp - 4.D0*C00f(t1b)*taucs**2*taucg**(-1)*
     &    a1mp + 2.D0*C00f(t1b)*taucg*a2mp - 2.D0*C00f(t1b)*taugs*a1mp
     &     - 2.D0*C00f(t1b)*taugs*a2mp - 2.D0*C00f(t1b)*a3mp + 4.D0*
     &    C00f(t2b)*taucs*taucg**(-1)*taugs*a1mp + 2.D0*C00f(t2b)*taucs
     &    *taucg**(-1)*taugs*a2mp + 4.D0*C00f(t2b)*taucs*taucg**(-1)*
     &    a3mp - 4.D0*C00f(t2b)*taucs*a2mp + 4.D0*C00f(t2b)*taucs**2*
     &    taucg**(-1)*a1mp + 2.D0*C00f(t2b)*taugs*a2mp - 2.D0*C00f(t3)*
     &    taucs*taucg**(-1)*taugs*a1mp + 10.D0*C00f(t3)*taucs*
     &    taucg**(-1)*a3mp + 10.D0*C00f(t3)*taugs*a1mp + 8.D0*C00f(t3)*
     &    taugs*a2mp )
      smp = smp + xn * ( 10.D0*C00f(t3)*a3mp + 2.D0*C00f(t6b)*taucs*
     &    a1mp - 2.D0*C00f(t6b)*taucg*a2mp + 2.D0*C00f(t6b)*taugs*a1mp
     &     + 2.D0*C00f(t6b)*taugs*a2mp + 2.D0*C00f(t6b)*a3mp )
      smp = smp + BB1(s1)*mc**2*xn**(-1) * ( taucg**(-2)*taugs*a3mp - 
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + BB1(s1)*mc**2*xn * (  - taucg**(-2)*taugs*a3mp + 
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + BB1(s1)*xn**(-1) * (  - taucs*taucg**(-1)*a3mp - 
     &    taugs*a1mp - taugs*a2mp - a3mp )
      smp = smp + BB1(s1)*xn * ( taucs*taucg**(-1)*a3mp + taugs*a1mp + 
     &    taugs*a2mp + a3mp )
      smp = smp + C11(t1)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs
     &    *a1mp - taucs*taucg**(-1)*taugs*a2mp - taucs*taucg**(-1)*a3mp
     &     - taugs*a1mp - taugs*a2mp - a3mp )
      smp = smp + C11(t1)*mc**2*xn * ( taucs*taucg**(-1)*taugs*a1mp + 
     &    taucs*taucg**(-1)*taugs*a2mp + taucs*taucg**(-1)*a3mp + taugs
     &    *a1mp + taugs*a2mp + a3mp )
      smp = smp + C11(t1)*mc**4*xn**(-1) * ( taucg**(-2)*taugs*a3mp + 
     &    taucg**(-2)*taugs**2*a1mp + taucg**(-2)*taugs**2*a2mp - 
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + C11(t1)*mc**4*xn * (  - taucg**(-2)*taugs*a3mp - 
     &    taucg**(-2)*taugs**2*a1mp - taucg**(-2)*taugs**2*a2mp + 
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + C11(t1b)*mc**2*xn**(-1) * ( taucg**(-2)*taugs**3*a1mp
     &     + taucg**(-2)*taugs**3*a2mp - taucg**(-1)*taugs**2*a1mp - 
     &    taucg**(-1)*taugs**2*a2mp )
      smp = smp + C11(t1b)*mc**2*xn * ( 2.D0*taucs*taucg**(-2)*taugs**2
     &    *a1mp + 2.D0*taucg**(-2)*taugs**2*a3mp + taucg**(-2)*taugs**3
     &    *a1mp + taucg**(-2)*taugs**3*a2mp - taucg**(-1)*taugs**2*a1mp
     &     - 3.D0*taucg**(-1)*taugs**2*a2mp )
      smp = smp + C11(t1b)*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2*
     &    a1mp - taucs*taucg**(-1)*taugs**2*a2mp + 2.D0*taucs*taugs*
     &    a1mp - 2.D0*taucg*taugs*a2mp + 2.D0*taugs*a3mp + taugs**2*
     &    a1mp + taugs**2*a2mp )
      smp = smp + C11(t1b)*xn * (  - 2.D0*taucs*taucg**(-1)*taugs*a3mp
     &     - taucs*taucg**(-1)*taugs**2*a1mp - taucs*taucg**(-1)*
     &    taugs**2*a2mp + 2.D0*taucs*taugs*a2mp - 2.D0*taucs**2*
     &    taucg**(-1)*taugs*a1mp - taugs**2*a1mp - taugs**2*a2mp )
      smp = smp + C11(t2)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs
     &    *a2mp + taucs*a2mp - taucg**(-1)*taugs**2*a2mp - taugs*a2mp )
      smp = smp + C11(t2)*mc**4*xn**(-1) * (  - taucg**(-1)*taugs*a2mp
     &     )
      smp = smp + C11(t2)*xn**(-1) * ( taucs*taucg*a2mp + taucs*taugs*
     &    a2mp + taucs**2*a2mp )
      smp = smp + C11(t8b)*mc**2*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1mp + taucs*taucg**(-1)*taugs*a2mp - taucs*a2mp )
      smp = smp + C11(t8b)*mc**4*xn**(-1) * ( taucg**(-1)*taugs*a2mp )
      smp = smp + C11(t8b)*xn**(-1) * (  - taucs**2*a1mp - taucs**2*
     &    a2mp )
      smp = smp + C22(t1)*mc**2*xn**(-1) * ( taucg**(-2)*taugs**2*a3mp
     &     - taucg**(-1)*taugs**2*a1mp )
      smp = smp + C22(t1)*mc**2*xn * (  - taucg**(-2)*taugs**2*a3mp + 
     &    taucg**(-1)*taugs**2*a1mp )
      smp = smp + C22(t1)*xn**(-1) * (  - taucs*taucg**(-1)*taugs*a3mp
     &     - taugs*a3mp - taugs**2*a1mp - taugs**2*a2mp )
      smp = smp + C22(t1)*xn * ( taucs*taucg**(-1)*taugs*a3mp + taugs*
     &    a3mp + taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + C22(t1b)*mc**2*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1mp + taucs*taucg**(-1)*taugs*a2mp - taucs*a2mp - taucg*a2mp
     &     + 2.D0*taugs*a1mp + 2.D0*taugs*a2mp )
      smp = smp + C22(t1b)*mc**2*xn * ( taucs*taucg**(-1)*taugs*a1mp + 
     &    taucs*taucg**(-1)*taugs*a2mp - taucs*a1mp - taucs**2*
     &    taucg**(-1)*a1mp )
      smp = smp + C22(t1b)*mc**4*xn**(-1) * ( taucg**(-1)*taugs*a1mp + 
     &    taucg**(-1)*taugs*a2mp )
      smp = smp + C22(t1b)*mc**4*xn * ( taucs*taucg**(-2)*taugs*a1mp - 
     &    taucg**(-1)*taugs*a2mp )
      smp = smp + C22(t1b)*xn**(-1) * (  - taucs*taucg*a1mp - taucs*
     &    taucg*a2mp - taucs**2*a1mp - taucs**2*a2mp )
      smp = smp + C22(t1b)*xn * (  - taucs*taucg*a1mp - taucs*taucg*
     &    a2mp - taucs**2*a1mp - taucs**2*a2mp )
      smp = smp + C22(t2b)*mc**2*xn * (  - taucs*taucg**(-1)*taugs*a2mp
     &     + taucs**2*taucg**(-1)*a1mp - taucg**(-1)*taugs**2*a1mp - 
     &    taucg**(-1)*taugs**2*a2mp + 2.D0*taugs*a2mp )
      smp = smp + C22(t2b)*mc**4*xn * (  - taucs*taucg**(-2)*taugs*a1mp
     &     - taucg**(-2)*taugs**2*a1mp + taucg**(-1)*taugs*a2mp )
      smp = smp + C22(t2b)*xn * ( taucs*taugs*a1mp + taucs*taugs*a2mp
     &     + taucs**2*a1mp + taucs**2*a2mp + taucg*taugs*a2mp )
      smp = smp + C22(t4b)*mc**2*xn**(-1) * (  - taucg*a2mp + 2.D0*
     &    taugs*a1mp + 2.D0*taugs*a2mp + 2.D0*a3mp )
      smp = smp + C22(t4b)*mc**4*xn**(-1) * ( taucg**(-1)*taugs*a1mp )
      smp = smp + C22(t4b)*xn**(-1) * (  - taucs*taucg*a1mp - taucs*
     &    taucg*a2mp )
      smp = smp + C22(t6b)*mc**2*xn * ( taucs*a1mp )
      smp = smp + C22(t6b)*xn * ( taucs*taucg*a1mp + taucs*taucg*a2mp )
      smp = smp + C22(t7)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1mp - 2.D0*taucg**(-1)*taugs**2*a1mp + taucg*a2mp
     &     - 2.D0*taugs*a1mp + taugs*a2mp )
      smp = smp + C22(t7)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a1mp - taucg**(-1)*taugs*a1mp )
      smp = smp + C22(t7)*xn**(-1) * ( taucs*taucg*a1mp + taucs*taucg*
     &    a2mp + taucs*taugs*a1mp + taucs*taugs*a2mp + taucs**2*
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + C22(t8)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs
     &    *a1mp - taucs*taucg**(-1)*taugs*a2mp + taucs*a2mp )
      smp = smp + C22(t8)*mc**4*xn**(-1) * (  - taucg**(-1)*taugs*a2mp
     &     )
      smp = smp + C22(t8)*xn**(-1) * ( taucs**2*a1mp + taucs**2*a2mp )
      smp = smp + C22(t9)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1mp - 2.D0*taucg**(-1)*taugs**2*a1mp + taucg*a2mp
     &     - 2.D0*taugs*a1mp - taugs*a2mp - 2.D0*a3mp )
      smp = smp + C22(t9)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a1mp - taucg**(-1)*taugs*a1mp )
      smp = smp + C22(t9)*xn**(-1) * ( taucs*taucg*a1mp + taucs*taucg*
     &    a2mp + taucs*taugs*a1mp - taucs*taugs*a2mp - 2.D0*taucs*a3mp
     &     + taucs**2*taucg**(-1)*taugs*a1mp )
      smp = smp + C12(t1)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a3mp + taucs*taucg**(-1)*taugs*a2mp + taucg**(-1)*taugs*a3mp
     &     - taugs*a1mp )
      smp = smp + C12(t1)*mc**2*xn * (  - taucs*taucg**(-2)*taugs*a3mp
     &     - taucs*taucg**(-1)*taugs*a2mp - taucg**(-1)*taugs*a3mp + 
     &    taugs*a1mp )
      smp = smp + C12(t1)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a1mp - taucg**(-2)*taugs**2*a2mp )
      smp = smp + C12(t1)*mc**4*xn * ( taucg**(-2)*taugs**2*a1mp + 
     &    taucg**(-2)*taugs**2*a2mp )
      smp = smp + C12(t1)*xn**(-1) * (  - taucs*taugs*a1mp - taucs*
     &    taugs*a2mp - 2.D0*taucs*a3mp - taucs**2*taucg**(-1)*a3mp - 
     &    taucg*taugs*a1mp - taucg*taugs*a2mp - taucg*a3mp )
      smp = smp + C12(t1)*xn * ( taucs*taugs*a1mp + taucs*taugs*a2mp + 
     &    2.D0*taucs*a3mp + taucs**2*taucg**(-1)*a3mp + taucg*taugs*
     &    a1mp + taucg*taugs*a2mp + taucg*a3mp )
      smp = smp + C12(t1b)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*
     &    taugs**2*a1mp + taucs*taucg**(-2)*taugs**2*a2mp - taucs*
     &    taucg**(-1)*taugs*a1mp - taucs*taucg**(-1)*taugs*a2mp + 3.D0*
     &    taucg**(-1)*taugs**2*a1mp + 3.D0*taucg**(-1)*taugs**2*a2mp - 
     &    2.D0*taugs*a2mp )
      smp = smp + C12(t1b)*mc**2*xn * ( taucs*taucg**(-2)*taugs*a3mp - 
     &    taucs*taucg**(-1)*taugs*a1mp - taucs*taucg**(-1)*taugs*a2mp
     &     + taucs**2*taucg**(-2)*taugs*a1mp + taucg**(-1)*taugs*a3mp
     &     + 2.D0*taucg**(-1)*taugs**2*a1mp + 2.D0*taucg**(-1)*taugs**2
     &    *a2mp - taugs*a1mp - 2.D0*taugs*a2mp )
      smp = smp + C12(t1b)*mc**4*xn**(-1) * ( taucg**(-2)*taugs**2*a1mp
     &     + taucg**(-2)*taugs**2*a2mp )
      smp = smp + C12(t1b)*mc**4*xn * ( taucg**(-2)*taugs**2*a1mp - 
     &    taucg**(-2)*taugs**2*a2mp )
      smp = smp + C12(t1b)*xn**(-1) * ( taucs*taucg*a1mp - taucs*taucg*
     &    a2mp - 2.D0*taucs*taugs*a1mp - 2.D0*taucs*taugs*a2mp + taucs*
     &    a3mp - taucs**2*taucg**(-1)*taugs*a1mp - taucs**2*taucg**(-1)
     &    *taugs*a2mp + taucs**2*a1mp + taucg*taugs*a1mp + taucg*taugs*
     &    a2mp + taucg*a3mp - taucg**2*a2mp )
      smp = smp + C12(t1b)*xn * ( taucs*taucg*a2mp - 3.D0*taucs*taugs*
     &    a1mp - 3.D0*taucs*taugs*a2mp - taucs*a3mp - taucs**2*
     &    taucg**(-1)*a3mp - taucs**2*a1mp + taucs**2*a2mp - taucs**3*
     &    taucg**(-1)*a1mp - taucg*taugs*a1mp - taucg*taugs*a2mp )
      smp = smp + C12(t2)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a2mp + taucs*taucg**(-1)*taugs*a2mp - taucg**(-2)*
     &    taugs**3*a2mp - taucg**(-1)*taugs**2*a2mp )
      smp = smp + C12(t2)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a2mp )
      smp = smp + C12(t2)*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*a2mp
     &     + taucs*taugs*a2mp + taucs**2*taucg**(-1)*taugs*a2mp )
      smp = smp + C12(t2b)*mc**2*xn * (  - taucs*taucg**(-2)*taugs*a3mp
     &     - 2.D0*taucs*taucg**(-2)*taugs**2*a1mp + taucs*taucg**(-1)*
     &    taugs*a2mp - taucs**2*taucg**(-2)*taugs*a1mp - taucg**(-2)*
     &    taugs**2*a3mp - taucg**(-2)*taugs**3*a1mp + 3.D0*taucg**(-1)*
     &    taugs**2*a2mp )
      smp = smp + C12(t2b)*mc**4*xn * ( taucg**(-2)*taugs**2*a2mp )
      smp = smp + C12(t2b)*xn * ( taucs*taucg**(-1)*taugs*a3mp + taucs*
     &    taucg**(-1)*taugs**2*a1mp - taucs*taugs*a2mp + 2.D0*taucs**2*
     &    taucg**(-1)*taugs*a1mp + taucs**2*taucg**(-1)*a3mp - taucs**2
     &    *a2mp + taucs**3*taucg**(-1)*a1mp + taugs**2*a2mp )
      smp = smp + C12(t3)*mc**2*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1mp )
      smp = smp + C12(t3)*mc**2*xn * (  - 3.D0*taucg**(-2)*taugs**2*
     &    a3mp + taucg**(-2)*taugs**3*a1mp + 2.D0*taucg**(-1)*taugs**2*
     &    a1mp )
      smp = smp + C12(t3)*xn**(-1) * (  - taugs*a3mp - taugs**2*a1mp - 
     &    taugs**2*a2mp )
      smp = smp + C12(t3)*xn * ( 3.D0*taucs*taucg**(-1)*taugs*a3mp - 
     &    taucs*taucg**(-1)*taugs**2*a1mp + 2.D0*taugs*a3mp + 2.D0*
     &    taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + C12(t3b)*mc**2*xn**(-1) * ( taucg**(-2)*taugs**3*a1mp
     &     )
      smp = smp + C12(t3b)*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2*
     &    a1mp - taugs**2*a2mp )
      smp = smp + C12(t4b)*xn**(-1) * ( taucs*taucg*a1mp + taucg*taugs*
     &    a1mp + taucg*taugs*a2mp + 2.D0*taucg*a3mp - taucg**2*a2mp )
      smp = smp + C12(t6b)*mc**2*xn * ( taugs*a1mp )
      smp = smp + C12(t6b)*xn * ( taucg*taugs*a1mp + taucg*taugs*a2mp )
      smp = smp + C12(t7)*mc**2*xn**(-1) * (  - taucg**(-2)*taugs**3*
     &    a1mp )
      smp = smp + C12(t7)*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*a1mp
     &     - taucs*taucg*a1mp - taucs*taugs*a1mp - taucg*taugs*a1mp + 
     &    taucg*taugs*a2mp - taucg*a3mp + taucg**2*a2mp - taugs*a3mp - 
     &    taugs**2*a1mp )
      smp = smp + C12(t8)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1mp - taucs*taucg**(-2)*taugs**2*a2mp + taucs*
     &    taucg**(-1)*taugs*a2mp )
      smp = smp + C12(t8)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a2mp )
      smp = smp + C12(t8)*xn**(-1) * ( taucs**2*taucg**(-1)*taugs*a1mp
     &     + taucs**2*taucg**(-1)*taugs*a2mp )
      smp = smp + C12(t8b)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*
     &    taugs**2*a1mp + taucs*taucg**(-2)*taugs**2*a2mp - taucs*
     &    taucg**(-1)*taugs*a2mp )
      smp = smp + C12(t8b)*mc**4*xn**(-1) * ( taucg**(-2)*taugs**2*a2mp
     &     )
      smp = smp + C12(t8b)*xn**(-1) * (  - taucs**2*taucg**(-1)*taugs*
     &    a1mp - taucs**2*taucg**(-1)*taugs*a2mp )
      smp = smp + C12(t9)*mc**2*xn**(-1) * (  - taucg**(-2)*taugs**3*
     &    a1mp )
      smp = smp + C12(t9)*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*a1mp
     &     - taucs*taucg*a1mp - taucs*taugs*a1mp - taucg*taugs*a1mp - 2.
     &    D0*taucg*a3mp + taucg**2*a2mp - 2.D0*taugs*a3mp - taugs**2*
     &    a1mp - taugs**2*a2mp )
      smp = smp + C0(t2)*mc**2*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a1mp + taucs*taucg**(-1)*taugs*a2mp + taucs*taucg**(-1)
     &    *a3mp - taucs*a2mp + taucs**2*taucg**(-1)*a1mp )
      smp = smp + C0(t2)*mc**4*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a1mp - taucg**(-2)*taugs*a3mp - 3.D0*taucg**(-2)*taugs**2*
     &    a1mp - taucg**(-2)*taugs**2*a2mp + taucg**(-1)*taugs*a2mp )
      smp = smp + C0(t3)*mc**2*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1mp )
      smp = smp + C0(t3)*xn**(-1) * (  - taugs*a3mp - taugs**2*a1mp - 
     &    taugs**2*a2mp )
      smp = smp + C0(t4b)*mc**2*xn**(-1) * (  - taucs*a1mp - 
     &    taucg**(-1)*taugs*a3mp - taucg**(-1)*taugs**2*a1mp - 
     &    taucg**(-1)*taugs**2*a2mp + taucg*a2mp - a3mp )
      smp = smp + C0(t4b)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a1mp + taucg**(-1)*taugs*a1mp )
      smp = smp + C0(t6b)*mc**2*xn * ( taucs*taucg**(-2)*taugs**2*a1mp
     &     + taucg**(-2)*taugs**3*a1mp + taucg**(-2)*taugs**3*a2mp )
      smp = smp + C0(t6b)*xn * (  - taucs*taucg**(-1)*taugs**2*a1mp - 
     &    taucs*taucg**(-1)*taugs**2*a2mp - taucs*taugs*a1mp - taucs*
     &    taugs*a2mp - taucs*a3mp - taucs**2*taucg**(-1)*taugs*a1mp + 
     &    taucg*taugs*a2mp )
      smp = smp + C0(t8)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1mp - taucg**(-1)*taugs*a3mp )
      smp = smp + C0(t8)*xn**(-1) * ( taucs*a3mp + taucs**2*taucg**(-1)
     &    *taugs*a1mp - taucg*taugs*a2mp - taugs**2*a1mp - taugs**2*
     &    a2mp )
      smp = smp + C0(t8b)*mc**2*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a1mp - taucs*taucg**(-1)*taugs*a2mp - taucs*taucg**(-1)
     &    *a3mp + taucs*a2mp - taucs**2*taucg**(-1)*a1mp )
      smp = smp + C0(t8b)*mc**4*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a1mp + taucg**(-2)*taugs*a3mp + 3.D0*taucg**(-2)*taugs**2*
     &    a1mp + taucg**(-2)*taugs**2*a2mp - taucg**(-1)*taugs*a2mp )
      smp = smp + C0(t8b)*xn**(-1) * ( taucs*taucg*a2mp - taucs*a3mp - 
     &    taucs**2*a1mp )
      smp = smp + C0(t9)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a1mp + taucs*a1mp - taucg*a2mp + a3mp )
      smp = smp + C0(t9)*mc**4*xn**(-1) * ( 2.D0*taucg**(-2)*taugs**2*
     &    a1mp - taucg**(-1)*taugs*a1mp )
      smp = smp + C1(t1)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a3mp - 2.D0*taucs*taucg**(-1)*taugs*a1mp - taucs*taucg**(-1)*
     &    taugs*a2mp - taucs*taucg**(-1)*a3mp + 2.D0*taucg**(-1)*taugs*
     &    a3mp + taucg**(-1)*taugs**2*a1mp + taucg**(-1)*taugs**2*a2mp
     &     - 2.D0*taugs*a1mp - taugs*a2mp - a3mp )
      smp = smp + C1(t1)*mc**2*xn * (  - taucs*taucg**(-2)*taugs*a3mp
     &     + 2.D0*taucs*taucg**(-1)*taugs*a1mp + taucs*taucg**(-1)*
     &    taugs*a2mp + taucs*taucg**(-1)*a3mp - 2.D0*taucg**(-1)*taugs*
     &    a3mp - taucg**(-1)*taugs**2*a1mp - taucg**(-1)*taugs**2*a2mp
     &     + 2.D0*taugs*a1mp + taugs*a2mp + a3mp )
      smp = smp + C1(t1)*mc**4*xn**(-1) * ( taucg**(-2)*taugs*a3mp + 2.D
     &    0*taucg**(-2)*taugs**2*a1mp + taucg**(-2)*taugs**2*a2mp - 
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + C1(t1)*mc**4*xn * (  - taucg**(-2)*taugs*a3mp - 2.D0*
     &    taucg**(-2)*taugs**2*a1mp - taucg**(-2)*taugs**2*a2mp + 
     &    taucg**(-1)*taugs*a1mp )
      smp = smp + C1(t1)*xn**(-1) * (  - taucs*taugs*a1mp - taucs*taugs
     &    *a2mp - 2.D0*taucs*a3mp - taucs**2*taucg**(-1)*a3mp - taucg*
     &    taugs*a1mp - taucg*taugs*a2mp - taucg*a3mp )
      smp = smp + C1(t1)*xn * ( taucs*taugs*a1mp + taucs*taugs*a2mp + 2.
     &    D0*taucs*a3mp + taucs**2*taucg**(-1)*a3mp + taucg*taugs*a1mp
     &     + taucg*taugs*a2mp + taucg*a3mp )
      smp = smp + C1(t1b)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a3mp + taucs*taucg**(-1)*taugs*a1mp + taucg**(-2)*taugs**3*
     &    a1mp + taucg**(-2)*taugs**3*a2mp - 2.D0*taucg**(-1)*taugs*
     &    a3mp - taucg**(-1)*taugs**2*a1mp - 2.D0*taucg**(-1)*taugs**2*
     &    a2mp + taugs*a1mp )
      smp = smp + C1(t1b)*mc**2*xn * ( taucs*taucg**(-2)*taugs*a3mp + 2.
     &    D0*taucs*taucg**(-2)*taugs**2*a1mp + taucs*taucg**(-1)*taugs*
     &    a2mp + taucg**(-2)*taugs**2*a3mp + 2.D0*taucg**(-2)*taugs**3*
     &    a1mp + 2.D0*taucg**(-2)*taugs**3*a2mp + 2.D0*taucg**(-1)*
     &    taugs*a3mp + taucg**(-1)*taugs**2*a1mp + taugs*a2mp )
      smp = smp + C1(t1b)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a1mp )
      smp = smp + C1(t1b)*mc**4*xn * ( taucg**(-2)*taugs**2*a1mp )
      smp = smp + C1(t1b)*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2*
     &    a1mp - taucs*taucg**(-1)*taugs**2*a2mp + taucs*taugs*a1mp + 2.
     &    D0*taucs*a3mp + taucs**2*taucg**(-1)*a3mp - taucg*taugs*a2mp
     &     + taucg*a3mp + taugs*a3mp - taugs**2*a1mp - taugs**2*a2mp )
      smp = smp + C1(t1b)*xn * (  - taucs*taucg**(-1)*taugs*a3mp - 2.D0
     &    *taucs*taucg**(-1)*taugs**2*a1mp - 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a2mp - 2.D0*taucs*taugs*a1mp - 2.D0*taucs*a3mp - 2.D0
     &    *taucs**2*taucg**(-1)*taugs*a1mp - taucs**2*taucg**(-1)*a3mp
     &     - taucg*taugs*a1mp - taucg*taugs*a2mp - taucg*a3mp - 
     &    taugs**2*a1mp - taugs**2*a2mp )
      smp = smp + C1(t2)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a3mp - taucs*taucg**(-1)*taugs*a2mp + taucs*taucg**(-1)*a3mp
     &     + taucs*a2mp + taucs**2*taucg**(-1)*a1mp + taucg**(-2)*
     &    taugs**2*a3mp - taucg**(-1)*taugs**2*a1mp - 2.D0*taucg**(-1)*
     &    taugs**2*a2mp - taugs*a2mp )
      smp = smp + C1(t2)*mc**4*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a1mp - taucg**(-2)*taugs*a3mp - taucg**(-2)*taugs**2*a1mp - 
     &    taucg**(-2)*taugs**2*a2mp - taucg**(-1)*taugs*a2mp )
      smp = smp + C1(t2)*xn**(-1) * (  - taucs*taucg**(-1)*taugs*a3mp
     &     + taucs*taucg*a2mp + taucs*taugs*a1mp + 2.D0*taucs*taugs*
     &    a2mp - taucs**2*taucg**(-1)*a3mp + taucs**2*a1mp + 2.D0*
     &    taucs**2*a2mp )
      smp = smp + C1(t2b)*mc**2*xn * (  - taucs*taucg**(-2)*taugs*a3mp
     &     - taucs*taucg**(-1)*taugs*a2mp - taucg**(-2)*taugs**2*a3mp
     &     - taucg**(-1)*taugs**2*a2mp )
      smp = smp + C1(t2b)*xn * ( taucs*taucg**(-1)*taugs*a3mp + 
     &    taucs**2*taucg**(-1)*a3mp )
      smp = smp + C1(t3)*mc**2*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1mp )
      smp = smp + C1(t3)*mc**2*xn * ( taucg**(-1)*taugs**2*a1mp + 
     &    taucg**(-1)*taugs**2*a2mp )
      smp = smp + C1(t3)*xn**(-1) * (  - taugs*a3mp - taugs**2*a1mp - 
     &    taugs**2*a2mp )
      smp = smp + C1(t3)*xn * ( taugs*a3mp + taugs**2*a1mp + taugs**2*
     &    a2mp )
      smp = smp + C1(t3b)*mc**2*xn**(-1) * ( taucg**(-1)*taugs**2*a1mp
     &     )
      smp = smp + C1(t3b)*xn**(-1) * ( taugs*a3mp + taugs**2*a2mp )
      smp = smp + C1(t6b)*mc**2*xn * (  - taucg**(-1)*taugs**2*a1mp )
      smp = smp + C1(t6b)*xn * ( taucs*taugs*a1mp + taucg*taugs*a1mp + 
     &    2.D0*taucg*taugs*a2mp + taucg*a3mp )
      smp = smp + C1(t7)*mc**2*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1mp - taugs*a1mp )
      smp = smp + C1(t7)*xn**(-1) * (  - taucg*taugs*a2mp - taucg*a3mp
     &     - taugs*a3mp - taugs**2*a2mp )
      smp = smp + C1(t8)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a3mp - taucs*taucg**(-2)*taugs**2*a1mp )
      smp = smp + C1(t8)*xn**(-1) * ( taucs*taugs*a2mp + taucs**2*
     &    taucg**(-1)*taugs*a1mp - taucs**2*taucg**(-1)*a3mp )
      smp = smp + C1(t8b)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*taugs
     &    *a3mp + taucs*taucg**(-1)*taugs*a1mp + taucs*taucg**(-1)*
     &    taugs*a2mp - taucs*taucg**(-1)*a3mp - taucs*a2mp - taucs**2*
     &    taucg**(-1)*a1mp )
      smp = smp + C1(t8b)*mc**4*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a1mp + taucg**(-2)*taugs*a3mp + 2.D0*taucg**(-2)*taugs**2*
     &    a1mp + taucg**(-2)*taugs**2*a2mp + taucg**(-1)*taugs*a2mp )
      smp = smp + C1(t8b)*xn**(-1) * ( taucs*taucg*a2mp - taucs*a3mp + 
     &    taucs**2*taucg**(-1)*a3mp - 3.D0*taucs**2*a1mp - 2.D0*
     &    taucs**2*a2mp )
      smp = smp + C2(t1)*mc**2*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1mp + taucs*taucg**(-1)*taugs*a2mp + taucg**(-2)*taugs**2*
     &    a3mp - taucg**(-1)*taugs*a3mp - 2.D0*taucg**(-1)*taugs**2*
     &    a1mp - taucg**(-1)*taugs**2*a2mp )
      smp = smp + C2(t1)*mc**2*xn * (  - taucs*taucg**(-1)*taugs*a1mp
     &     - taucs*taucg**(-1)*taugs*a2mp - taucg**(-2)*taugs**2*a3mp
     &     + taucg**(-1)*taugs*a3mp + 2.D0*taucg**(-1)*taugs**2*a1mp + 
     &    taucg**(-1)*taugs**2*a2mp )
      smp = smp + C2(t1)*mc**4*xn**(-1) * (  - 2.D0*taucg**(-2)*
     &    taugs**2*a1mp - taucg**(-2)*taugs**2*a2mp )
      smp = smp + C2(t1)*mc**4*xn * ( 2.D0*taucg**(-2)*taugs**2*a1mp + 
     &    taucg**(-2)*taugs**2*a2mp )
      smp = smp + C2(t1)*xn**(-1) * (  - taucs*taucg**(-1)*taugs*a3mp
     &     - taugs*a3mp - taugs**2*a1mp - taugs**2*a2mp )
      smp = smp + C2(t1)*xn * ( taucs*taucg**(-1)*taugs*a3mp + taugs*
     &    a3mp + taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + C2(t1b)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a3mp + taucs*taucg**(-2)*taugs**2*a1mp + taucs*taucg**(-2)*
     &    taugs**2*a2mp - taucs*taucg**(-1)*taugs*a1mp + taucs*
     &    taucg**(-1)*a3mp + 2.D0*taucg**(-1)*taugs*a3mp + 3.D0*
     &    taucg**(-1)*taugs**2*a1mp + 3.D0*taucg**(-1)*taugs**2*a2mp - 
     &    taugs*a1mp - taugs*a2mp + a3mp )
      smp = smp + C2(t1b)*mc**2*xn * (  - taucs*taucg**(-2)*taugs*a3mp
     &     - 2.D0*taucs*taucg**(-1)*taugs*a1mp - taucs*taucg**(-1)*a3mp
     &     - 2.D0*taucg**(-1)*taugs*a3mp - taugs*a1mp - 2.D0*taugs*a2mp
     &     - a3mp )
      smp = smp + C2(t1b)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs*a3mp
     &     + 2.D0*taucg**(-2)*taugs**2*a1mp + taucg**(-1)*taugs*a1mp )
      smp = smp + C2(t1b)*mc**4*xn * ( taucg**(-2)*taugs*a3mp + 
     &    taucg**(-2)*taugs**2*a1mp + taucg**(-1)*taugs*a1mp + 2.D0*
     &    taucg**(-1)*taugs*a2mp )
      smp = smp + C2(t1b)*xn**(-1) * (  - 3.D0*taucs*taugs*a1mp - 3.D0*
     &    taucs*taugs*a2mp - 2.D0*taucs*a3mp - taucs**2*taucg**(-1)*
     &    taugs*a1mp - taucs**2*taucg**(-1)*taugs*a2mp - taucs**2*
     &    taucg**(-1)*a3mp - taucg*taugs*a1mp - taucg*taugs*a2mp - 
     &    taucg*a3mp )
      smp = smp + C2(t1b)*xn * (  - taucs*taugs*a1mp - taucs*taugs*a2mp
     &     + 2.D0*taucs*a3mp + taucs**2*taucg**(-1)*a3mp + taucg*a3mp )
      smp = smp + C2(t2)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a3mp - 2.D0*taucs*taucg**(-2)*taugs**2*a1mp - taucs*
     &    taucg**(-2)*taugs**2*a2mp + 2.D0*taucs*taucg**(-1)*taugs*a2mp
     &     - taucs**2*taucg**(-2)*taugs*a1mp - taucg**(-2)*taugs**2*
     &    a3mp - taucg**(-2)*taugs**3*a1mp - taucg**(-2)*taugs**3*a2mp
     &     )
      smp = smp + C2(t2)*mc**4*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a2mp )
      smp = smp + C2(t2)*xn**(-1) * ( taucs*taucg**(-1)*taugs*a3mp + 
     &    taucs*taucg**(-1)*taugs**2*a1mp + taucs*taucg**(-1)*taugs**2*
     &    a2mp + 2.D0*taucs**2*taucg**(-1)*taugs*a1mp + taucs**2*
     &    taucg**(-1)*taugs*a2mp + taucs**2*taucg**(-1)*a3mp - taucs**2
     &    *a2mp + taucs**3*taucg**(-1)*a1mp )
      smp = smp + C2(t2b)*mc**2*xn * ( taucs*taucg**(-2)*taugs*a3mp + 
     &    taucs*taucg**(-1)*taugs*a1mp + taucs*taucg**(-1)*a3mp + 
     &    taucg**(-2)*taugs**2*a3mp + taucg**(-1)*taugs**2*a1mp + 
     &    taucg**(-1)*taugs**2*a2mp - 2.D0*taugs*a2mp )
      smp = smp + C2(t2b)*mc**4*xn * (  - taucg**(-2)*taugs*a3mp - 2.D0
     &    *taucg**(-1)*taugs*a2mp )
      smp = smp + C2(t2b)*xn * (  - taucs*taucg**(-1)*taugs*a3mp - 
     &    taucs*taugs*a2mp - taucs**2*taucg**(-1)*a3mp )
      smp = smp + C2(t3)*mc**2*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a1mp )
      smp = smp + C2(t3)*mc**2*xn * (  - taucg**(-2)*taugs**2*a3mp )
      smp = smp + C2(t3)*xn**(-1) * (  - taugs*a3mp - taugs**2*a1mp - 
     &    taugs**2*a2mp )
      smp = smp + C2(t3)*xn * ( taucs*taucg**(-1)*taugs*a3mp + taugs*
     &    a3mp + taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + C2(t4b)*mc**2*xn**(-1) * (  - taucs*a1mp - 
     &    taucg**(-1)*taugs*a3mp + 3.D0*taugs*a1mp + 2.D0*taugs*a2mp + 
     &    a3mp )
      smp = smp + C2(t4b)*mc**4*xn**(-1) * ( 2.D0*taucg**(-1)*taugs*
     &    a1mp )
      smp = smp + C2(t4b)*xn**(-1) * (  - taucs*taucg*a1mp - taucs*
     &    taucg*a2mp + taucs*a3mp + taucg*taugs*a1mp + taucg*taugs*a2mp
     &     + taucg*a3mp )
      smp = smp + C2(t6b)*mc**2*xn * ( 3.D0*taucs*taucg**(-1)*taugs*
     &    a1mp + 2.D0*taucg**(-1)*taugs*a3mp + taugs*a1mp + taugs*a2mp
     &     + a3mp )
      smp = smp + C2(t6b)*mc**4*xn * (  - 2.D0*taucg**(-2)*taugs**2*
     &    a1mp - taucg**(-1)*taugs*a1mp )
      smp = smp + C2(t6b)*xn * ( taucs*taucg*a2mp - 2.D0*taucs*a3mp - 
     &    taucs**2*a1mp - taucg*a3mp )
      smp = smp + C2(t7)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a1mp + taugs*a1mp - taugs*a2mp - a3mp )
      smp = smp + C2(t7)*mc**4*xn**(-1) * (  - taucg**(-1)*taugs*a1mp )
      smp = smp + C2(t7)*xn**(-1) * ( taucs*taugs*a1mp + taucg*taugs*
     &    a1mp + taucg*taugs*a2mp + taucg*a3mp + taugs*a3mp + taugs**2*
     &    a1mp + taugs**2*a2mp )
      smp = smp + C2(t8)*mc**2*xn**(-1) * (  - taucs*taucg**(-2)*taugs*
     &    a3mp - taucs*taucg**(-2)*taugs**2*a1mp - taucs*taucg**(-2)*
     &    taugs**2*a2mp + taucs*taucg**(-1)*taugs*a1mp - taucs*
     &    taucg**(-1)*a3mp - taucg**(-1)*taugs*a3mp - taucg**(-1)*
     &    taugs**2*a1mp - taucg**(-1)*taugs**2*a2mp + 2.D0*taugs*a2mp )
      smp = smp + C2(t8)*mc**4*xn**(-1) * ( taucg**(-2)*taugs*a3mp - 2.D
     &    0*taucg**(-2)*taugs**2*a1mp )
      smp = smp + C2(t8)*xn**(-1) * (  - taucs*taucg*a2mp + 2.D0*taucs*
     &    taugs*a1mp + 2.D0*taucs*taugs*a2mp + 2.D0*taucs*a3mp + 
     &    taucs**2*taucg**(-1)*taugs*a1mp + taucs**2*taucg**(-1)*taugs*
     &    a2mp + taucs**2*taucg**(-1)*a3mp + taucs**2*a1mp )
      smp = smp + C2(t8b)*mc**2*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a3mp + 2.D0*taucs*taucg**(-2)*taugs**2*a1mp + taucs*
     &    taucg**(-2)*taugs**2*a2mp - 2.D0*taucs*taucg**(-1)*taugs*a2mp
     &     + taucs**2*taucg**(-2)*taugs*a1mp )
      smp = smp + C2(t8b)*mc**4*xn**(-1) * ( taucg**(-2)*taugs**2*a2mp
     &     )
      smp = smp + C2(t8b)*xn**(-1) * (  - 2.D0*taucs**2*taucg**(-1)*
     &    taugs*a1mp - taucs**2*taucg**(-1)*taugs*a2mp - taucs**2*
     &    taucg**(-1)*a3mp + taucs**2*a2mp - taucs**3*taucg**(-1)*a1mp
     &     )
      smp = smp + C2(t9)*mc**2*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a1mp + taucs*a1mp - 2.D0*taucg**(-1)*taugs**2*a1mp - 3.D
     &    0*taugs*a1mp - taugs*a2mp - a3mp )
      smp = smp + C2(t9)*mc**4*xn**(-1) * (  - 2.D0*taucg**(-1)*taugs*
     &    a1mp )
      smp = smp + C2(t9)*xn**(-1) * ( taucs*taucg*a1mp - taucs*taugs*
     &    a2mp - taucs*a3mp + taucs**2*a1mp - taucg*taugs*a1mp - taucg*
     &    taugs*a2mp - taucg*a3mp - taugs*a3mp - taugs**2*a1mp - 
     &    taugs**2*a2mp )
      smp = smp + D0(b2)*mc**2*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3mp + taucs*taucg**(-1)*taugs**2*a1mp + taucs*taucg**(-1)*
     &    taugs**2*a2mp - taucs*taucg*a2mp - taucs*taugs*a2mp + taucs*
     &    a3mp + taucs**2*a1mp + taucg**(-1)*taugs**2*a3mp + 
     &    taucg**(-1)*taugs**3*a1mp + taucg**(-1)*taugs**3*a2mp )
      smp = smp + D0(b2)*mc**4*xn**(-1) * ( taucs*taucg**(-2)*taugs**2*
     &    a1mp - taucs*taucg**(-1)*taugs*a1mp + taucg**(-2)*taugs**3*
     &    a1mp - taucg**(-1)*taugs*a3mp - taucg**(-1)*taugs**2*a1mp + 
     &    taucg**(-1)*taugs**2*a2mp + taugs*a2mp )
      smp = smp + D0(b2)*mc**6*xn**(-1) * (  - 2.D0*taucg**(-2)*
     &    taugs**2*a1mp )
      smp = smp + D1(b1)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*
     &    taugs**2*a1mp )
      smp = smp + D1(b1)*xn**(-1) * (  - taucs*taugs*a3mp - taucs*
     &    taugs**2*a2mp )
      smp = smp + D1(b2)*mc**2*xn**(-1) * ( taucg**(-1)*taugs**2*a3mp
     &     + taucg**(-1)*taugs**3*a1mp + taucg**(-1)*taugs**3*a2mp )
      smp = smp + D1(b2)*mc**4*xn**(-1) * ( taucg**(-2)*taugs**3*a1mp
     &     - taucg**(-1)*taugs**2*a1mp )
      smp = smp + D1(b3)*mc**2*xn * (  - taucg**(-1)*taugs**2*a3mp - 
     &    taugs**2*a2mp )
      smp = smp + D1(b3)*xn * ( taucs*taugs*a3mp )
      smp = smp + D2(b1)*mc**2*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a3mp )
      smp = smp + D2(b1)*xn**(-1) * ( taucs*taugs*a3mp - taucs*taugs**2
     &    *a1mp - taucs*taugs**2*a2mp - taucg*taugs**2*a2mp - taugs**3*
     &    a1mp - taugs**3*a2mp )
      smp = smp + D2(b2)*mc**2*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a3mp + taucs*taucg**(-1)*taugs**2*a1mp + taucs*
     &    taucg**(-1)*taugs**2*a2mp - 2.D0*taucs*taugs*a1mp - 2.D0*
     &    taucs*taugs*a2mp + taucs*a3mp + taucs**2*a1mp )
      smp = smp + D2(b2)*mc**4*xn**(-1) * ( taucs*taucg**(-2)*taugs**2*
     &    a1mp - taucs*taucg**(-1)*taugs*a1mp - taucg**(-1)*taugs*a3mp
     &     + taucg**(-1)*taugs**2*a1mp + taucg**(-1)*taugs**2*a2mp )
      smp = smp + D2(b2)*mc**6*xn**(-1) * (  - 2.D0*taucg**(-2)*
     &    taugs**2*a1mp )
      smp = smp + D2(b2)*xn**(-1) * ( taucs**2*taucg*a1mp + taucs**2*
     &    taucg*a2mp - taucs**2*a3mp )
      smp = smp + D2(b3)*mc**2*xn * ( taucs*taucg**(-2)*taugs**3*a1mp
     &     + taucg**(-2)*taugs**4*a1mp + taucg**(-2)*taugs**4*a2mp + 
     &    taucg**(-1)*taugs**3*a1mp + taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + D2(b3)*xn * (  - taucs*taucg**(-1)*taugs**3*a1mp - 
     &    taucs*taucg**(-1)*taugs**3*a2mp - taucs*taugs*a3mp - 2.D0*
     &    taucs*taugs**2*a1mp - taucs*taugs**2*a2mp - taucs**2*
     &    taucg**(-1)*taugs**2*a1mp )
      smp = smp + D3(b1)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a3mp - taucs*taugs*a1mp - taucs**2*taucg**(-1)*taugs*a1mp - 
     &    taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + D3(b1)*mc**4*xn**(-1) * ( taucg**(-1)*taugs**2*a1mp )
      smp = smp + D3(b1)*xn**(-1) * (  - 2.D0*taucs*taucg*taugs*a2mp - 
     &    taucs**2*taugs*a2mp - taucg*taugs**2*a1mp - taucg*taugs**2*
     &    a2mp - taucg**2*taugs*a2mp )
      smp = smp + D3(b2)*mc**2*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a3mp - 2.D0*taucs*taucg**(-1)*taugs**2*a1mp + 2.D0*taucs*
     &    taugs*a2mp - taucs**2*taucg**(-1)*taugs*a1mp )
      smp = smp + D3(b2)*mc**4*xn**(-1) * (  - taucs*taucg**(-2)*
     &    taugs**2*a1mp + taucg**(-2)*taugs**3*a1mp - taucg**(-1)*
     &    taugs**2*a2mp )
      smp = smp + D3(b2)*xn**(-1) * (  - taucs**2*taucg*a2mp + taucs**2
     &    *taugs*a1mp + taucs**2*a3mp + taucs**3*a1mp )
      smp = smp + D3(b3)*mc**2*xn * ( taucs*taucg**(-2)*taugs**3*a1mp
     &     + taucs*taucg**(-2)*taugs**3*a2mp - taucs*taucg**(-1)*taugs*
     &    a3mp + taucs*taucg**(-1)*taugs**2*a1mp - taucs*taugs*a2mp + 
     &    taucs**2*taucg**(-2)*taugs**2*a1mp - taucg*taugs*a2mp + 
     &    taugs**2*a1mp + taugs**2*a2mp )
      smp = smp + D3(b3)*mc**4*xn * ( taucg**(-1)*taugs**2*a1mp + 2.D0*
     &    taucg**(-1)*taugs**2*a2mp )
      smp = smp + D3(b3)*xn * (  - taucs*taucg*taugs*a1mp - taucs*taucg
     &    *taugs*a2mp - taucs**2*taucg**(-1)*taugs**2*a1mp - taucs**2*
     &    taucg**(-1)*taugs**2*a2mp - 2.D0*taucs**2*taugs*a1mp - 
     &    taucs**2*taugs*a2mp - taucs**3*taucg**(-1)*taugs*a1mp )
      smp = smp - epinv*taucs*taucg**(-1)*b0*a3mp - epinv*taugs*b0*a1mp
     &     - epinv*taugs*b0*a2mp - epinv*b0*a3mp

      spm=  + mc*xn**(-1) * ( 9.D0/4.D0*epinv*taucs*taucg**(-1)*a2pm + 
     &    5.D0/2.D0*taucs*taucg**(-1)*a2pm - 2.D0*C00s(t1b)*epinv*taucs
     &    *taucg**(-1)*a2pm + 4.D0*C00s(t1b)*epinv*taucs*a1pm - 2.D0*
     &    C00s(t1b)*epinv*taucs**2*taucg**(-1)*a1pm + 2.D0*C00s(t1b)*
     &    epinv*taucg*a3pm + 2.D0*C00s(t1b)*epinv*taugs*a1pm + 2.D0*
     &    C00s(t1b)*epinv*taugs*a3pm + 4.D0*C00s(t1b)*epinv*a2pm + 4.D0
     &    *C00s(t2)*epinv*taucs*taucg**(-1)*taugs*a1pm + 2.D0*C00s(t2)*
     &    epinv*taucs**2*taucg**(-1)*a1pm + 4.D0*C00s(t2)*epinv*
     &    taucg**(-1)*taugs*a2pm + 2.D0*C00s(t2)*epinv*taucg**(-1)*
     &    taugs**2*a1pm + 2.D0*C00s(t2)*epinv*taucg**(-1)*taugs**2*a3pm
     &     + 2.D0*C00s(t2)*epinv*taugs*a3pm - 2.D0*C00s(t2)*taucs*
     &    taucg**(-1)*taugs*a1pm - 2.D0*C00s(t2)*taucs*taucg**(-1)*
     &    taugs*a3pm - 2.D0*C00s(t2)*taucs*taucg**(-1)*a2pm - 2.D0*
     &    C00s(t2)*taucs*a3pm - 2.D0*C00s(t2)*taucs**2*taucg**(-1)*a1pm
     &     - 2.D0*C00s(t3)*taucs*a1pm + 6.D0*C00s(t3b)*epinv*taucs*
     &    taucg**(-1)*taugs*a1pm )
      spm = spm + mc*xn**(-1) * ( 6.D0*C00s(t3b)*epinv*taucg**(-1)*
     &    taugs*a2pm + 4.D0*C00s(t3b)*epinv*taucg**(-1)*taugs**2*a1pm
     &     + 4.D0*C00s(t3b)*epinv*taucg**(-1)*taugs**2*a3pm + 4.D0*
     &    C00s(t3b)*epinv*taugs*a3pm - 2.D0*C00s(t4)*epinv*taucs*
     &    taucg**(-1)*a2pm + 2.D0*C00s(t4)*taucs*taucg**(-1)*a2pm + 4.D0
     &    *C00s(t4b)*epinv*taucs*a1pm + 4.D0*C00s(t4b)*epinv*taucg*a3pm
     &     + 4.D0*C00s(t4b)*epinv*taugs*a1pm + 4.D0*C00s(t4b)*epinv*
     &    taugs*a3pm + 6.D0*C00s(t4b)*epinv*a2pm - 6.D0*C00s(t7)*epinv*
     &    taucs*taucg**(-1)*taugs*a1pm - 4.D0*C00s(t7)*epinv*taucs*a1pm
     &     - 6.D0*C00s(t7)*epinv*taucg**(-1)*taugs*a2pm - 4.D0*C00s(t7)
     &    *epinv*taucg**(-1)*taugs**2*a1pm - 4.D0*C00s(t7)*epinv*
     &    taucg**(-1)*taugs**2*a3pm - 2.D0*C00s(t7)*epinv*taucg*a3pm - 
     &    2.D0*C00s(t7)*epinv*taugs*a1pm - 6.D0*C00s(t7)*epinv*taugs*
     &    a3pm - 4.D0*C00s(t7)*epinv*a2pm + 2.D0*C00s(t8)*epinv*taucs*
     &    taucg**(-1)*a2pm + 2.D0*C00s(t8)*epinv*taucs**2*taucg**(-1)*
     &    a1pm )
      spm = spm + mc*xn**(-1) * (  - 2.D0*C00s(t8b)*epinv*taucs*
     &    taucg**(-1)*a2pm - 2.D0*C00s(t8b)*epinv*taucs**2*taucg**(-1)*
     &    a1pm - 4.D0*C00s(t9)*epinv*taucs*taucg**(-1)*taugs*a1pm - 4.D0
     &    *C00s(t9)*epinv*taucs*a1pm - 4.D0*C00s(t9)*epinv*taucg**(-1)*
     &    taugs*a2pm - 2.D0*C00s(t9)*epinv*taucg**(-1)*taugs**2*a1pm - 
     &    2.D0*C00s(t9)*epinv*taucg**(-1)*taugs**2*a3pm - 4.D0*C00s(t9)
     &    *epinv*taucg*a3pm - 4.D0*C00s(t9)*epinv*taugs*a1pm - 6.D0*
     &    C00s(t9)*epinv*taugs*a3pm - 6.D0*C00s(t9)*epinv*a2pm - 2.D0*
     &    C00f(t1b)*taucs*taucg**(-1)*a2pm + 4.D0*C00f(t1b)*taucs*a1pm
     &     - 2.D0*C00f(t1b)*taucs**2*taucg**(-1)*a1pm + 2.D0*C00f(t1b)*
     &    taucg*a3pm + 2.D0*C00f(t1b)*taugs*a1pm + 2.D0*C00f(t1b)*taugs
     &    *a3pm + 4.D0*C00f(t1b)*a2pm + 4.D0*C00f(t2)*taucs*taucg**(-1)
     &    *taugs*a1pm + 2.D0*C00f(t2)*taucs**2*taucg**(-1)*a1pm + 4.D0*
     &    C00f(t2)*taucg**(-1)*taugs*a2pm + 2.D0*C00f(t2)*taucg**(-1)*
     &    taugs**2*a1pm + 2.D0*C00f(t2)*taucg**(-1)*taugs**2*a3pm + 2.D0
     &    *C00f(t2)*taugs*a3pm )
      spm = spm + mc*xn**(-1) * ( 6.D0*C00f(t3b)*taucs*taucg**(-1)*
     &    taugs*a1pm + 6.D0*C00f(t3b)*taucg**(-1)*taugs*a2pm + 4.D0*
     &    C00f(t3b)*taucg**(-1)*taugs**2*a1pm + 4.D0*C00f(t3b)*
     &    taucg**(-1)*taugs**2*a3pm + 4.D0*C00f(t3b)*taugs*a3pm - 2.D0*
     &    C00f(t4)*taucs*taucg**(-1)*a2pm + 4.D0*C00f(t4b)*taucs*a1pm
     &     + 4.D0*C00f(t4b)*taucg*a3pm + 4.D0*C00f(t4b)*taugs*a1pm + 4.D
     &    0*C00f(t4b)*taugs*a3pm + 6.D0*C00f(t4b)*a2pm - 6.D0*C00f(t7)*
     &    taucs*taucg**(-1)*taugs*a1pm - 4.D0*C00f(t7)*taucs*a1pm - 6.D0
     &    *C00f(t7)*taucg**(-1)*taugs*a2pm - 4.D0*C00f(t7)*taucg**(-1)*
     &    taugs**2*a1pm - 4.D0*C00f(t7)*taucg**(-1)*taugs**2*a3pm - 2.D0
     &    *C00f(t7)*taucg*a3pm - 2.D0*C00f(t7)*taugs*a1pm - 6.D0*C00f(
     &    t7)*taugs*a3pm - 4.D0*C00f(t7)*a2pm + 2.D0*C00f(t8)*taucs*
     &    taucg**(-1)*a2pm + 2.D0*C00f(t8)*taucs**2*taucg**(-1)*a1pm - 
     &    2.D0*C00f(t8b)*taucs*taucg**(-1)*a2pm - 2.D0*C00f(t8b)*
     &    taucs**2*taucg**(-1)*a1pm - 4.D0*C00f(t9)*taucs*taucg**(-1)*
     &    taugs*a1pm )
      spm = spm + mc*xn**(-1) * (  - 4.D0*C00f(t9)*taucs*a1pm - 4.D0*
     &    C00f(t9)*taucg**(-1)*taugs*a2pm - 2.D0*C00f(t9)*taucg**(-1)*
     &    taugs**2*a1pm - 2.D0*C00f(t9)*taucg**(-1)*taugs**2*a3pm - 4.D0
     &    *C00f(t9)*taucg*a3pm - 4.D0*C00f(t9)*taugs*a1pm - 6.D0*C00f(
     &    t9)*taugs*a3pm - 6.D0*C00f(t9)*a2pm + 3.D0/2.D0*lnrat(musq,
     &    msq)*taucs*taucg**(-1)*a2pm )
      spm = spm + mc*xn * (  - 9.D0/4.D0*epinv*taucs*taucg**(-1)*a2pm
     &     - 7.D0/3.D0*taucs*taucg**(-1)*a2pm + 2.D0*C00s(t1b)*epinv*
     &    taucs*taucg**(-1)*taugs*a1pm + 2.D0*C00s(t1b)*epinv*taucs*
     &    taucg**(-1)*taugs*a3pm + 4.D0*C00s(t1b)*epinv*taucs*
     &    taucg**(-1)*a2pm - 2.D0*C00s(t1b)*epinv*taucs*a1pm + 4.D0*
     &    C00s(t1b)*epinv*taucs**2*taucg**(-1)*a1pm - 2.D0*C00s(t1b)*
     &    epinv*taucg*a3pm - 2.D0*C00s(t1b)*epinv*taugs*a1pm - 2.D0*
     &    C00s(t1b)*epinv*taugs*a3pm - 2.D0*C00s(t1b)*epinv*a2pm + 2.D0
     &    *C00s(t2)*epinv*taucs*taucg**(-1)*a2pm + 2.D0*C00s(t2)*taucs*
     &    taucg**(-1)*taugs*a1pm + 2.D0*C00s(t2)*taucs*taucg**(-1)*
     &    taugs*a3pm + 2.D0*C00s(t2)*taucs*taucg**(-1)*a2pm + 2.D0*
     &    C00s(t2)*taucs*a3pm + 2.D0*C00s(t2)*taucs**2*taucg**(-1)*a1pm
     &     - 8.D0*C00s(t2b)*epinv*taucs*taucg**(-1)*taugs*a1pm - 2.D0*
     &    C00s(t2b)*epinv*taucs*taucg**(-1)*taugs*a3pm - 4.D0*C00s(t2b)
     &    *epinv*taucs*taucg**(-1)*a2pm - 4.D0*C00s(t2b)*epinv*taucs**2
     &    *taucg**(-1)*a1pm )
      spm = spm + mc*xn * (  - 6.D0*C00s(t2b)*epinv*taucg**(-1)*taugs*
     &    a2pm - 4.D0*C00s(t2b)*epinv*taucg**(-1)*taugs**2*a1pm - 4.D0*
     &    C00s(t2b)*epinv*taucg**(-1)*taugs**2*a3pm - 4.D0*C00s(t2b)*
     &    epinv*taugs*a3pm + 6.D0*C00s(t3)*epinv*taucs*taucg**(-1)*
     &    taugs*a1pm + 6.D0*C00s(t3)*epinv*taucg**(-1)*taugs*a2pm + 4.D0
     &    *C00s(t3)*epinv*taucg**(-1)*taugs**2*a1pm + 4.D0*C00s(t3)*
     &    epinv*taucg**(-1)*taugs**2*a3pm + 4.D0*C00s(t3)*epinv*taugs*
     &    a3pm - 2.D0*C00s(t3)*taucs*a1pm + 10.D0*C00s(t6)*epinv*taucs*
     &    taucg**(-1)*a2pm - 4.D0*C00s(t6)*taucs*taucg**(-1)*a2pm + 2.D0
     &    *C00s(t6b)*epinv*taucs*a1pm + 2.D0*C00s(t6b)*epinv*taucg*a3pm
     &     + 2.D0*C00s(t6b)*epinv*taugs*a1pm + 2.D0*C00s(t6b)*epinv*
     &    taugs*a3pm + 2.D0*C00s(t6b)*epinv*a2pm + 2.D0*C00f(t1b)*taucs
     &    *taucg**(-1)*taugs*a1pm + 2.D0*C00f(t1b)*taucs*taucg**(-1)*
     &    taugs*a3pm + 4.D0*C00f(t1b)*taucs*taucg**(-1)*a2pm - 2.D0*
     &    C00f(t1b)*taucs*a1pm + 4.D0*C00f(t1b)*taucs**2*taucg**(-1)*
     &    a1pm )
      spm = spm + mc*xn * (  - 2.D0*C00f(t1b)*taucg*a3pm - 2.D0*C00f(
     &    t1b)*taugs*a1pm - 2.D0*C00f(t1b)*taugs*a3pm - 2.D0*C00f(t1b)*
     &    a2pm + 2.D0*C00f(t2)*taucs*taucg**(-1)*a2pm - 8.D0*C00f(t2b)*
     &    taucs*taucg**(-1)*taugs*a1pm - 2.D0*C00f(t2b)*taucs*
     &    taucg**(-1)*taugs*a3pm - 4.D0*C00f(t2b)*taucs*taucg**(-1)*
     &    a2pm - 4.D0*C00f(t2b)*taucs**2*taucg**(-1)*a1pm - 6.D0*C00f(
     &    t2b)*taucg**(-1)*taugs*a2pm - 4.D0*C00f(t2b)*taucg**(-1)*
     &    taugs**2*a1pm - 4.D0*C00f(t2b)*taucg**(-1)*taugs**2*a3pm - 4.D
     &    0*C00f(t2b)*taugs*a3pm + 6.D0*C00f(t3)*taucs*taucg**(-1)*
     &    taugs*a1pm + 6.D0*C00f(t3)*taucg**(-1)*taugs*a2pm + 4.D0*
     &    C00f(t3)*taucg**(-1)*taugs**2*a1pm + 4.D0*C00f(t3)*
     &    taucg**(-1)*taugs**2*a3pm + 4.D0*C00f(t3)*taugs*a3pm + 10.D0*
     &    C00f(t6)*taucs*taucg**(-1)*a2pm + 2.D0*C00f(t6b)*taucs*a1pm
     &     + 2.D0*C00f(t6b)*taucg*a3pm + 2.D0*C00f(t6b)*taugs*a1pm + 2.D
     &    0*C00f(t6b)*taugs*a3pm + 2.D0*C00f(t6b)*a2pm - 3.D0/2.D0*
     &    lnrat(musq,msq)*taucs*taucg**(-1)*a2pm )
      spm = spm + mc * (  - epinv*taucs*taucg**(-1)*b0*a2pm )
      spm = spm + mc**3*xn**(-1) * ( 3.D0*epinv*taucs*taucg**(-2)*a2pm
     &     - 9.D0/4.D0*epinv*taucg**(-2)*taugs*a2pm + 5.D0*taucs*
     &    taucg**(-2)*a2pm - 5.D0/2.D0*taucg**(-2)*taugs*a2pm + 2.D0*
     &    C00s(t1b)*epinv*taucg**(-1)*taugs*a1pm + 2.D0*C00s(t2)*epinv*
     &    taucg**(-2)*taugs*a2pm + 2.D0*C00s(t2)*taucs*taucg**(-2)*
     &    taugs*a1pm + 2.D0*C00s(t2)*taucg**(-2)*taugs*a2pm + 2.D0*
     &    C00s(t2)*taucg**(-2)*taugs**2*a1pm + 2.D0*C00s(t2)*
     &    taucg**(-2)*taugs**2*a3pm + 2.D0*C00s(t2)*taucg**(-1)*taugs*
     &    a3pm + 2.D0*C00s(t3)*taucg**(-1)*taugs*a1pm + 2.D0*C00s(t4)*
     &    epinv*taucg**(-2)*taugs*a2pm - 2.D0*C00s(t4)*taucg**(-2)*
     &    taugs*a2pm - 2.D0*C00s(t7)*epinv*taucg**(-1)*taugs*a1pm + 2.D0
     &    *C00f(t1b)*taucg**(-1)*taugs*a1pm + 2.D0*C00f(t2)*taucg**(-2)
     &    *taugs*a2pm + 2.D0*C00f(t4)*taucg**(-2)*taugs*a2pm - 2.D0*
     &    C00f(t7)*taucg**(-1)*taugs*a1pm + 3.D0*lnrat(musq,msq)*taucs*
     &    taucg**(-2)*a2pm - 3.D0/2.D0*lnrat(musq,msq)*taucg**(-2)*
     &    taugs*a2pm )
      spm = spm + mc**3*xn * (  - 3.D0*epinv*taucs*taucg**(-2)*a2pm + 9.
     &    D0/4.D0*epinv*taucg**(-2)*taugs*a2pm - 5.D0*taucs*taucg**(-2)
     &    *a2pm + 7.D0/3.D0*taucg**(-2)*taugs*a2pm + 2.D0*C00s(t1b)*
     &    epinv*taucg**(-1)*taugs*a3pm - 2.D0*C00s(t2)*epinv*
     &    taucg**(-2)*taugs*a2pm - 2.D0*C00s(t2)*taucs*taucg**(-2)*
     &    taugs*a1pm - 2.D0*C00s(t2)*taucg**(-2)*taugs*a2pm - 2.D0*
     &    C00s(t2)*taucg**(-2)*taugs**2*a1pm - 2.D0*C00s(t2)*
     &    taucg**(-2)*taugs**2*a3pm - 2.D0*C00s(t2)*taucg**(-1)*taugs*
     &    a3pm - 2.D0*C00s(t2b)*epinv*taucg**(-1)*taugs*a3pm + 2.D0*
     &    C00s(t3)*taucg**(-1)*taugs*a1pm - 10.D0*C00s(t6)*epinv*
     &    taucg**(-2)*taugs*a2pm + 4.D0*C00s(t6)*taucg**(-2)*taugs*a2pm
     &     + 2.D0*C00f(t1b)*taucg**(-1)*taugs*a3pm - 2.D0*C00f(t2)*
     &    taucg**(-2)*taugs*a2pm - 2.D0*C00f(t2b)*taucg**(-1)*taugs*
     &    a3pm - 10.D0*C00f(t6)*taucg**(-2)*taugs*a2pm - 3.D0*lnrat(
     &    musq,msq)*taucs*taucg**(-2)*a2pm + 3.D0/2.D0*lnrat(musq,msq)*
     &    taucg**(-2)*taugs*a2pm )
      spm = spm + mc**3 * ( epinv*taucg**(-2)*taugs*b0*a2pm )
      spm = spm + mc**5*xn**(-1) * (  - 3.D0*epinv*taucg**(-3)*taugs*
     &    a2pm - 5.D0*taucg**(-3)*taugs*a2pm - 3.D0*lnrat(musq,msq)*
     &    taucg**(-3)*taugs*a2pm )
      spm = spm + mc**5*xn * ( 3.D0*epinv*taucg**(-3)*taugs*a2pm + 5.D0
     &    *taucg**(-3)*taugs*a2pm + 3.D0*lnrat(musq,msq)*taucg**(-3)*
     &    taugs*a2pm )
      spm = spm + BB1(s2)*mc*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*a2pm )
      spm = spm + BB1(s2)*mc*xn * (  - 2.D0*taucs*taucg**(-1)*a2pm )
      spm = spm + BB1(s2)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-2)*
     &    a2pm - 2.D0*taucg**(-2)*taugs*a2pm )
      spm = spm + BB1(s2)*mc**3*xn * (  - 2.D0*taucs*taucg**(-2)*a2pm
     &     + 2.D0*taucg**(-2)*taugs*a2pm )
      spm = spm + BB1(s2)*mc**5*xn**(-1) * (  - 2.D0*taucg**(-3)*taugs*
     &    a2pm )
      spm = spm + BB1(s2)*mc**5*xn * ( 2.D0*taucg**(-3)*taugs*a2pm )
      spm = spm + BB0(s2)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-2)*
     &    a2pm )
      spm = spm + BB0(s2)*mc**3*xn * ( 2.D0*taucs*taucg**(-2)*a2pm )
      spm = spm + BB0(s2)*mc**5*xn**(-1) * ( 2.D0*taucg**(-3)*taugs*
     &    a2pm )
      spm = spm + BB0(s2)*mc**5*xn * (  - 2.D0*taucg**(-3)*taugs*a2pm )
      spm = spm + C11(t1b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pm + taucs*taucg**(-1)*taugs**2*a3pm + 2.D0*taucs*taugs*
     &    a1pm + taucg*taugs*a1pm + taucg*taugs*a3pm + 2.D0*taugs*a2pm
     &     + taugs**2*a1pm + taugs**2*a3pm )
      spm = spm + C11(t1b)*mc*xn * ( 2.D0*taucs*taucg**(-1)*taugs*a2pm
     &     + taucs*taucg**(-1)*taugs**2*a1pm + taucs*taucg**(-1)*
     &    taugs**2*a3pm + 2.D0*taucs**2*taucg**(-1)*taugs*a1pm - taucg*
     &    taugs*a1pm - taucg*taugs*a3pm - taugs**2*a1pm - taugs**2*a3pm
     &     )
      spm = spm + C11(t1b)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pm
     &     + taucg**(-1)*taugs**2*a3pm )
      spm = spm + C11(t1b)*mc**3*xn * ( taucg**(-1)*taugs**2*a1pm + 
     &    taucg**(-1)*taugs**2*a3pm )
      spm = spm + C11(t2)*mc*xn**(-1) * (  - taucs*taucg*a3pm - 2.D0*
     &    taucs*taugs*a3pm - taucs*a2pm + taucg*taugs*a3pm + 2.D0*taugs
     &    *a2pm + taugs**2*a3pm )
      spm = spm + C11(t2)*mc*xn * ( taucs*a2pm - taucs**2*a3pm )
      spm = spm + C11(t2)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a3pm - taucs*taucg**(-1)*a2pm - taucs*a3pm + 3.D0*
     &    taucg**(-1)*taugs*a2pm + 2.D0*taucg**(-1)*taugs**2*a3pm + 
     &    taugs*a3pm )
      spm = spm + C11(t2)*mc**3*xn * ( 2.D0*taucs*taucg**(-1)*taugs*
     &    a3pm + taucs*taucg**(-1)*a2pm - taucg**(-1)*taugs*a2pm )
      spm = spm + C11(t2)*mc**5*xn**(-1) * ( taucg**(-2)*taugs*a2pm + 
     &    taucg**(-2)*taugs**2*a3pm )
      spm = spm + C11(t2)*mc**5*xn * (  - taucg**(-2)*taugs*a2pm - 
     &    taucg**(-2)*taugs**2*a3pm )
      spm = spm + C11(t3b)*mc*xn**(-1) * ( taugs**2*a1pm )
      spm = spm + C11(t4)*mc*xn**(-1) * ( taucs*a2pm )
      spm = spm + C11(t4)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*a2pm - 
     &    taucg**(-1)*taugs*a2pm )
      spm = spm + C11(t4)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs*a2pm
     &     )
      spm = spm + C11(t4b)*mc*xn**(-1) * (  - taucg*taugs*a1pm )
      spm = spm + C11(t6)*mc*xn * ( taucs*a2pm )
      spm = spm + C11(t6)*mc**3*xn * ( 4.D0*taucs*taucg**(-1)*a2pm - 
     &    taucg**(-1)*taugs*a2pm )
      spm = spm + C11(t6)*mc**5*xn * (  - 4.D0*taucg**(-2)*taugs*a2pm )
      spm = spm + C11(t6b)*mc*xn * ( taucg*taugs*a1pm )
      spm = spm + C11(t7)*mc*xn**(-1) * (  - taucg*taugs*a1pm - 
     &    taugs**2*a1pm )
      spm = spm + C11(t8b)*mc*xn**(-1) * ( taucs**2*a1pm + taucs**2*
     &    a3pm )
      spm = spm + C11(t8b)*mc**3*xn**(-1) * ( taucs*a3pm )
      spm = spm + C11(t9)*mc*xn**(-1) * ( taucg*taugs*a1pm + taugs**2*
     &    a1pm )
      spm = spm + C22(t1b)*mc*xn**(-1) * ( taucs**2*a1pm + taucs**2*
     &    a3pm )
      spm = spm + C22(t1b)*mc*xn * (  - taucs**2*a1pm - taucs**2*a3pm )
      spm = spm + C22(t1b)*mc**3*xn**(-1) * ( 2.D0*taucs*a1pm + taucs*
     &    a3pm - taucg*a3pm - taugs*a1pm - taugs*a3pm )
      spm = spm + C22(t1b)*mc**3*xn * (  - taucs*taucg**(-1)*taugs*a1pm
     &     - taucs*taucg**(-1)*taugs*a3pm - taucs*a1pm - 2.D0*taucs*
     &    a3pm + taucs**2*taucg**(-1)*a1pm )
      spm = spm + C22(t2b)*mc*xn * ( 2.D0*taucs*taucg*a3pm - 2.D0*taucs
     &    *taugs*a1pm + 3.D0*taucs*taugs*a3pm - taucs**2*a1pm + 
     &    taucs**2*a3pm - 2.D0*taugs*a2pm - taugs**2*a1pm )
      spm = spm + C22(t2b)*mc**3*xn * (  - 2.D0*taucs*taucg**(-1)*taugs
     &    *a1pm + taucs*taucg**(-1)*taugs*a3pm + 2.D0*taucs*a3pm - 
     &    taucs**2*taucg**(-1)*a1pm - 2.D0*taucg**(-1)*taugs*a2pm - 
     &    taucg**(-1)*taugs**2*a1pm - taucg**(-1)*taugs**2*a3pm )
      spm = spm + C22(t3)*mc*xn * ( taugs**2*a1pm )
      spm = spm + C22(t4b)*mc**3*xn**(-1) * ( 2.D0*taucs*a1pm + taucg*
     &    a3pm + taugs*a1pm + taugs*a3pm + 2.D0*a2pm )
      spm = spm + C22(t4b)*mc**5*xn**(-1) * (  - 2.D0*taucg**(-1)*taugs
     &    *a1pm )
      spm = spm + C22(t6b)*mc**3*xn * ( taucs*a1pm )
      spm = spm + C22(t7)*mc*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a2pm - taucs*taucg**(-1)*taugs**2*a1pm - taucs*
     &    taucg**(-1)*taugs**2*a3pm + taucs*taucg*a3pm + taucs*taugs*
     &    a1pm - 2.D0*taucs**2*taucg**(-1)*taugs*a1pm - 2.D0*taucs**2*
     &    a1pm )
      spm = spm + C22(t7)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pm - 2.D0*taucs*a1pm - 2.D0*taucg**(-1)*taugs*a2pm - 
     &    taucg**(-1)*taugs**2*a1pm - taucg**(-1)*taugs**2*a3pm + taucg
     &    *a3pm + taugs*a1pm )
      spm = spm + C22(t8)*mc*xn**(-1) * (  - taucs**2*a1pm - taucs**2*
     &    a3pm )
      spm = spm + C22(t8)*mc**3*xn**(-1) * (  - taucs*a3pm )
      spm = spm + C22(t9)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pm + taucs*taucg**(-1)*taugs**2*a3pm - taucs*taucg*a3pm - 
     &    taucs*taugs*a1pm - 2.D0*taucs*a2pm - 2.D0*taucs**2*a1pm )
      spm = spm + C22(t9)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pm - 2.D0*taucs*a1pm + taucg**(-1)*taugs**2*a1pm + 
     &    taucg**(-1)*taugs**2*a3pm - taucg*a3pm - taugs*a1pm - 2.D0*
     &    a2pm )
      spm = spm + C22(t9)*mc**5*xn**(-1) * ( 2.D0*taucg**(-1)*taugs*
     &    a1pm )
      spm = spm + C12(t1b)*mc*xn**(-1) * ( 2.D0*taucs*taucg*a1pm + 
     &    taucs*taucg*a3pm + 2.D0*taucs*taugs*a1pm + 2.D0*taucs*taugs*
     &    a3pm + taucs*a2pm + taucs**2*taucg**(-1)*taugs*a1pm + 
     &    taucs**2*taucg**(-1)*taugs*a3pm + taucs**2*a1pm + taucg*a2pm
     &     )
      spm = spm + C12(t1b)*mc*xn * (  - taucs*taucg*a1pm - taucs*taucg*
     &    a3pm - taucs*taugs*a1pm - taucs*taugs*a3pm + taucs*a2pm + 
     &    taucs**2*taucg**(-1)*a2pm + taucs**2*a1pm + taucs**3*
     &    taucg**(-1)*a1pm )
      spm = spm + C12(t1b)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1pm + taucs*taucg**(-1)*taugs*a3pm - taucg**(-1)*taugs**2*
     &    a1pm - taucg**(-1)*taugs**2*a3pm + 2.D0*taugs*a1pm )
      spm = spm + C12(t1b)*mc**3*xn * ( 2.D0*taucs*taucg**(-1)*taugs*
     &    a1pm - taucg**(-1)*taugs**2*a1pm - taucg**(-1)*taugs**2*a3pm
     &     - taugs*a1pm - taugs*a3pm )
      spm = spm + C12(t2)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2
     &    *a3pm - taucs**2*taucg**(-1)*taugs*a3pm - taucs**2*
     &    taucg**(-1)*a2pm - taucs**2*a3pm + taucg**(-1)*taugs**2*a2pm
     &     )
      spm = spm + C12(t2)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pm + 
     &    taucs**2*taucg**(-1)*a2pm + taucs**2*a3pm )
      spm = spm + C12(t2)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a2pm + taucs*taucg**(-1)*taugs*a3pm + taucg**(-2)*taugs**2*
     &    a2pm - taucg**(-1)*taugs**2*a3pm )
      spm = spm + C12(t2)*mc**3*xn * (  - taucs*taucg**(-2)*taugs*a2pm
     &     - 2.D0*taucs*taucg**(-1)*taugs*a3pm - taucg**(-2)*taugs**2*
     &    a2pm )
      spm = spm + C12(t2)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a3pm )
      spm = spm + C12(t2)*mc**5*xn * ( taucg**(-2)*taugs**2*a3pm )
      spm = spm + C12(t2b)*mc*xn * (  - 3.D0*taucs*taucg**(-1)*taugs*
     &    a2pm - 3.D0*taucs*taucg**(-1)*taugs**2*a1pm - taucs*
     &    taucg**(-1)*taugs**2*a3pm - taucs*taugs*a3pm - 3.D0*taucs**2*
     &    taucg**(-1)*taugs*a1pm - taucs**2*taucg**(-1)*a2pm - taucs**3
     &    *taucg**(-1)*a1pm - 2.D0*taucg**(-1)*taugs**2*a2pm - 
     &    taucg**(-1)*taugs**3*a1pm - taucg**(-1)*taugs**3*a3pm - 
     &    taugs**2*a3pm )
      spm = spm + C12(t3)*mc*xn * ( 2.D0*taucs*taucg**(-1)*taugs**2*
     &    a1pm + 2.D0*taucg**(-1)*taugs**2*a2pm + taucg**(-1)*taugs**3*
     &    a1pm + taucg**(-1)*taugs**3*a3pm + taugs**2*a3pm )
      spm = spm + C12(t3b)*mc*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a1pm + 2.D0*taucg**(-1)*taugs**2*a2pm + taucg**(-1)*
     &    taugs**3*a1pm + taucg**(-1)*taugs**3*a3pm + taugs**2*a3pm )
      spm = spm + C12(t4b)*mc*xn**(-1) * ( 2.D0*taucs*taucg*a1pm + 
     &    taucg*taugs*a1pm + taucg*taugs*a3pm + 2.D0*taucg*a2pm + 
     &    taucg**2*a3pm )
      spm = spm + C12(t4b)*mc**3*xn**(-1) * (  - 3.D0*taugs*a1pm )
      spm = spm + C12(t6)*mc*xn * ( 3.D0*taucs*a2pm )
      spm = spm + C12(t6)*mc**3*xn * (  - 3.D0*taucg**(-1)*taugs*a2pm )
      spm = spm + C12(t6b)*mc*xn * ( taucs*taucg*a1pm )
      spm = spm + C12(t6b)*mc**3*xn * ( taugs*a1pm )
      spm = spm + C12(t7)*mc*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a1pm - 2.D0*taucs*taucg*a1pm - 6.D0*taucs*taugs*a1pm
     &     - 2.D0*taucg**(-1)*taugs**2*a2pm - taucg**(-1)*taugs**3*a1pm
     &     - taucg**(-1)*taugs**3*a3pm - taucg*taugs*a3pm - taucg*a2pm
     &     - 3.D0*taugs*a2pm - taugs**2*a1pm - 2.D0*taugs**2*a3pm )
      spm = spm + C12(t7)*mc**3*xn**(-1) * (  - 2.D0*taugs*a1pm )
      spm = spm + C12(t8)*mc*xn**(-1) * (  - taucs**2*taucg**(-1)*taugs
     &    *a1pm - taucs**2*taucg**(-1)*taugs*a3pm )
      spm = spm + C12(t8)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs
     &    *a3pm )
      spm = spm + C12(t8b)*mc*xn**(-1) * ( taucs**2*taucg**(-1)*taugs*
     &    a1pm + taucs**2*taucg**(-1)*taugs*a3pm )
      spm = spm + C12(t8b)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3pm )
      spm = spm + C12(t9)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2
     &    *a1pm - 2.D0*taucs*taucg*a1pm - taucs*taugs*a1pm - 
     &    taucg**(-1)*taugs**2*a2pm - taucg*taugs*a1pm - 2.D0*taucg*
     &    taugs*a3pm - 2.D0*taucg*a2pm - taucg**2*a3pm - 3.D0*taugs*
     &    a2pm - taugs**2*a1pm - taugs**2*a3pm )
      spm = spm + C12(t9)*mc**3*xn**(-1) * ( taucg**(-1)*taugs**2*a1pm
     &     + 3.D0*taugs*a1pm )
      spm = spm + C0(t2)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pm - taucs*taucg**(-1)*taugs*a3pm - taucs*taucg**(-1)
     &    *a2pm - taucs**2*taucg**(-1)*a1pm - taucg**(-1)*taugs*a2pm - 
     &    taucg**(-1)*taugs**2*a1pm - taucg**(-1)*taugs**2*a3pm - taugs
     &    *a3pm )
      spm = spm + C0(t2)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs*a3pm )
      spm = spm + C0(t4)*mc*xn**(-1) * (  - taucs*a2pm )
      spm = spm + C0(t4)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    a2pm + taucg**(-1)*taugs*a2pm )
      spm = spm + C0(t4)*mc**5*xn**(-1) * ( 2.D0*taucg**(-2)*taugs*a2pm
     &     )
      spm = spm + C0(t4b)*mc*xn**(-1) * (  - taucs*a2pm )
      spm = spm + C0(t4b)*mc**3*xn**(-1) * (  - 2.D0*taucs*a1pm + 
     &    taucg**(-1)*taugs*a2pm - taucg*a3pm - taugs*a1pm - taugs*a3pm
     &     - a2pm )
      spm = spm + C0(t4b)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs*a1pm
     &     )
      spm = spm + C0(t6)*mc*xn * ( taucs*a2pm )
      spm = spm + C0(t6)*mc**3*xn * ( 2.D0*taucs*taucg**(-1)*a2pm - 
     &    taucg**(-1)*taugs*a2pm )
      spm = spm + C0(t6)*mc**5*xn * (  - 2.D0*taucg**(-2)*taugs*a2pm )
      spm = spm + C0(t6b)*mc*xn * (  - 2.D0*taucs*taugs*a1pm - taucg*
     &    taugs*a3pm - 2.D0*taugs*a2pm - taugs**2*a1pm - taugs**2*a3pm
     &     )
      spm = spm + C0(t6b)*mc**3*xn * (  - taucg**(-1)*taugs*a2pm + 
     &    taucg**(-1)*taugs**2*a1pm )
      spm = spm + C0(t8)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pm - 2.D0*taucs**2*taucg**(-1)*taugs*a1pm )
      spm = spm + C0(t8)*mc**3*xn**(-1) * ( taugs*a3pm )
      spm = spm + C0(t8b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs*a2pm
     &     - taucs*taucg*a3pm - 2.D0*taucs*taugs*a1pm - taucs*taugs*
     &    a3pm - taucs*a2pm - taucs**2*a1pm )
      spm = spm + C0(t8b)*mc**3*xn**(-1) * ( 4.D0*taucs*taucg**(-1)*
     &    taugs*a1pm + taucs*taucg**(-1)*taugs*a3pm + taucs*taucg**(-1)
     &    *a2pm + taucs**2*taucg**(-1)*a1pm + 3.D0*taucg**(-1)*taugs*
     &    a2pm + 2.D0*taucg**(-1)*taugs**2*a1pm + taucg**(-1)*taugs**2*
     &    a3pm + taugs*a3pm )
      spm = spm + C0(t8b)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a3pm )
      spm = spm + C0(t9)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1pm + 2.D0*taucs*a1pm + taucg**(-1)*taugs*a2pm + taucg**(-1)
     &    *taugs**2*a1pm + taucg**(-1)*taugs**2*a3pm + taucg*a3pm + 
     &    taugs*a1pm + 2.D0*taugs*a3pm + a2pm )
      spm = spm + C0(t9)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a1pm )
      spm = spm + C1(t1b)*mc*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a1pm + 2.D0*taucs*taucg**(-1)*taugs**2*a3pm + 2.D0*
     &    taucs*taugs*a1pm + taucs*taugs*a3pm - taucs**2*taucg**(-1)*
     &    a2pm + taucs**2*a1pm + taucg*taugs*a1pm + taugs*a2pm )
      spm = spm + C1(t1b)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pm + 
     &    taucs*taucg*a1pm + taucs*taucg*a3pm - taucs*taugs*a1pm + 
     &    taucs**2*taucg**(-1)*taugs*a1pm + taucs**2*taucg**(-1)*taugs*
     &    a3pm + taucs**2*taucg**(-1)*a2pm + taucs**2*a3pm - 2.D0*taucg
     &    *taugs*a1pm - taucg*taugs*a3pm - taugs**2*a1pm - taugs**2*
     &    a3pm )
      spm = spm + C1(t1b)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs
     &    *a1pm + taucs*taucg**(-1)*taugs*a3pm + taucg**(-1)*taugs*a2pm
     &     + taucg**(-1)*taugs**2*a3pm + taugs*a3pm )
      spm = spm + C1(t1b)*mc**3*xn * (  - taucs*taucg**(-1)*taugs*a3pm
     &     - taucg**(-1)*taugs*a2pm + 2.D0*taucg**(-1)*taugs**2*a1pm + 
     &    taucg**(-1)*taugs**2*a3pm - taugs*a1pm - taugs*a3pm )
      spm = spm + C1(t2)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs*a2pm
     &     - taucs*taucg*a3pm - 2.D0*taucs*taugs*a1pm - 2.D0*taucs*
     &    taugs*a3pm - taucs*a2pm - taucs**2*a1pm + taucg**(-1)*
     &    taugs**2*a2pm + taugs*a2pm - taugs**2*a1pm )
      spm = spm + C1(t2)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pm + taucs
     &    *a2pm + taucs**2*taucg**(-1)*a2pm - taucs**2*a3pm )
      spm = spm + C1(t2)*mc**3*xn**(-1) * ( taucs*taucg**(-2)*taugs*
     &    a2pm - 2.D0*taucs*taucg**(-1)*taugs*a1pm - 4.D0*taucs*
     &    taucg**(-1)*taugs*a3pm - 2.D0*taucs*taucg**(-1)*a2pm - 2.D0*
     &    taucs*a3pm - taucs**2*taucg**(-1)*a1pm + taucg**(-2)*taugs**2
     &    *a2pm + 2.D0*taucg**(-1)*taugs*a2pm - taucg**(-1)*taugs**2*
     &    a1pm - taugs*a3pm )
      spm = spm + C1(t2)*mc**3*xn * (  - taucs*taucg**(-2)*taugs*a2pm
     &     + 2.D0*taucs*taucg**(-1)*taugs*a3pm + taucs*taucg**(-1)*a2pm
     &     - taucg**(-2)*taugs**2*a2pm - taucg**(-1)*taugs*a2pm )
      spm = spm + C1(t2)*mc**5*xn**(-1) * ( taucg**(-2)*taugs*a2pm + 
     &    taucg**(-2)*taugs**2*a3pm - taucg**(-1)*taugs*a3pm )
      spm = spm + C1(t2)*mc**5*xn * (  - taucg**(-2)*taugs*a2pm - 
     &    taucg**(-2)*taugs**2*a3pm )
      spm = spm + C1(t2b)*mc*xn * (  - 2.D0*taucs*taucg**(-1)*taugs*
     &    a2pm - taucs*taucg**(-1)*taugs**2*a3pm - taucs*taugs*a3pm - 
     &    taucs**2*taucg**(-1)*taugs*a3pm - taucs**2*taucg**(-1)*a2pm
     &     - taucs**2*a3pm - taucg**(-1)*taugs**2*a2pm )
      spm = spm + C1(t2b)*mc**3*xn * ( taucs*taucg**(-1)*taugs*a3pm + 
     &    taucg**(-1)*taugs**2*a3pm )
      spm = spm + C1(t3)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pm + taucs
     &    *taucg**(-1)*taugs**2*a3pm + taucs*taugs*a3pm + taucg**(-1)*
     &    taugs**2*a2pm )
      spm = spm + C1(t3)*mc**3*xn * (  - taucg**(-1)*taugs**2*a3pm )
      spm = spm + C1(t3b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a1pm + taugs**2*a1pm )
      spm = spm + C1(t4)*mc**3*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    a2pm )
      spm = spm + C1(t4)*mc**5*xn**(-1) * ( 2.D0*taucg**(-2)*taugs*a2pm
     &     )
      spm = spm + C1(t4b)*mc*xn**(-1) * ( taucs*taucg*a1pm - taucg*
     &    taugs*a1pm + taucg*a2pm + taugs*a2pm )
      spm = spm + C1(t4b)*mc**3*xn**(-1) * (  - 2.D0*taugs*a1pm )
      spm = spm + C1(t6)*mc*xn * ( 3.D0*taucs*a2pm )
      spm = spm + C1(t6)*mc**3*xn * ( 6.D0*taucs*taucg**(-1)*a2pm - 3.D0
     &    *taucg**(-1)*taugs*a2pm )
      spm = spm + C1(t6)*mc**5*xn * (  - 6.D0*taucg**(-2)*taugs*a2pm )
      spm = spm + C1(t6b)*mc*xn * (  - taucs*taucg*a1pm - taucs*taugs*
     &    a1pm + 2.D0*taucg*taugs*a1pm - taucg*taugs*a3pm - 2.D0*taugs*
     &    a2pm - taugs**2*a1pm - taugs**2*a3pm )
      spm = spm + C1(t6b)*mc**3*xn * ( taugs*a1pm )
      spm = spm + C1(t7)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs**2*
     &    a1pm - 2.D0*taucs*taugs*a1pm - taucg*taugs*a1pm - taugs**2*
     &    a1pm )
      spm = spm + C1(t8)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pm - taucs*taucg**(-1)*taugs**2*a1pm - taucs*taucg**(-1)*
     &    taugs**2*a3pm - taucs*taugs*a3pm - 2.D0*taucs**2*taucg**(-1)*
     &    taugs*a1pm + taucs**2*taucg**(-1)*a2pm )
      spm = spm + C1(t8)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a3pm )
      spm = spm + C1(t8b)*mc*xn**(-1) * (  - taucs*taucg*a3pm - taucs*
     &    taugs*a1pm - taucs*taugs*a3pm - taucs*a2pm - taucs**2*
     &    taucg**(-1)*a2pm + taucs**2*a1pm + taucs**2*a3pm )
      spm = spm + C1(t8b)*mc**3*xn**(-1) * ( 4.D0*taucs*taucg**(-1)*
     &    taugs*a1pm + 2.D0*taucs*taucg**(-1)*taugs*a3pm + taucs*
     &    taucg**(-1)*a2pm + 2.D0*taucs*a3pm + taucs**2*taucg**(-1)*
     &    a1pm + 3.D0*taucg**(-1)*taugs*a2pm + 2.D0*taucg**(-1)*
     &    taugs**2*a1pm + 2.D0*taucg**(-1)*taugs**2*a3pm + 2.D0*taugs*
     &    a3pm )
      spm = spm + C1(t8b)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a3pm )
      spm = spm + C1(t9)*mc*xn**(-1) * (  - taucs*taucg*a1pm - 
     &    taucg**(-1)*taugs**2*a2pm + taucg*taugs*a1pm - taucg*a2pm - 2.
     &    D0*taugs*a2pm + taugs**2*a1pm )
      spm = spm + C1(t9)*mc**3*xn**(-1) * ( 2.D0*taugs*a1pm )
      spm = spm + C2(t1b)*mc*xn**(-1) * ( taucs*taucg*a1pm + taucs*
     &    taugs*a1pm + taucs*taugs*a3pm + 2.D0*taucs**2*taucg**(-1)*
     &    taugs*a1pm + 2.D0*taucs**2*taucg**(-1)*taugs*a3pm + taucs**2*
     &    taucg**(-1)*a2pm + taucs**2*a3pm - taucg*taugs*a1pm - taucg*
     &    taugs*a3pm - taucg**2*a3pm )
      spm = spm + C2(t1b)*mc*xn * (  - taucs*taucg*a1pm - taucs*taugs*
     &    a1pm - taucs*taugs*a3pm - taucs**2*taucg**(-1)*taugs*a1pm - 
     &    taucs**2*taucg**(-1)*taugs*a3pm - taucs**2*taucg**(-1)*a2pm
     &     - taucs**2*a3pm + taucg*taugs*a1pm + taucg*taugs*a3pm + 
     &    taucg**2*a3pm )
      spm = spm + C2(t1b)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs*a1pm + taucs*taucg**(-1)*taugs*a3pm - taucs*taucg**(-1)
     &    *a2pm + taucs*a1pm - taucg**(-1)*taugs*a2pm - taucg**(-1)*
     &    taugs**2*a1pm - taucg**(-1)*taugs**2*a3pm - taugs*a1pm + a2pm
     &     )
      spm = spm + C2(t1b)*mc**3*xn * ( 2.D0*taucs*taucg**(-1)*taugs*
     &    a3pm + taucs*taucg**(-1)*a2pm + taucs*a1pm + 2.D0*taucs*a3pm
     &     + taucg**(-1)*taugs*a2pm - taucg**(-1)*taugs**2*a1pm - 
     &    taucg**(-1)*taugs**2*a3pm - 2.D0*taugs*a3pm - a2pm )
      spm = spm + C2(t1b)*mc**5*xn**(-1) * (  - taucg**(-1)*taugs*a1pm
     &     + 2.D0*taucg**(-1)*taugs*a3pm )
      spm = spm + C2(t1b)*mc**5*xn * (  - taucg**(-1)*taugs*a1pm - 2.D0
     &    *taucg**(-1)*taugs*a3pm )
      spm = spm + C2(t2)*mc*xn**(-1) * (  - 2.D0*taucs*taucg**(-1)*
     &    taugs*a2pm - 3.D0*taucs*taucg**(-1)*taugs**2*a1pm - 2.D0*
     &    taucs*taucg**(-1)*taugs**2*a3pm - taucs*taugs*a3pm - 3.D0*
     &    taucs**2*taucg**(-1)*taugs*a1pm - taucs**2*taucg**(-1)*taugs*
     &    a3pm - taucs**2*taucg**(-1)*a2pm - taucs**2*a3pm - taucs**3*
     &    taucg**(-1)*a1pm - taucg**(-1)*taugs**2*a2pm - taucg**(-1)*
     &    taugs**3*a1pm - taucg**(-1)*taugs**3*a3pm - taugs**2*a3pm )
      spm = spm + C2(t2)*mc*xn * ( taucs**2*a3pm )
      spm = spm + C2(t2)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3pm - taucg**(-1)*taugs**2*a3pm )
      spm = spm + C2(t2)*mc**3*xn * (  - 2.D0*taucs*taucg**(-1)*taugs*
     &    a3pm )
      spm = spm + C2(t2)*mc**5*xn**(-1) * (  - taucg**(-2)*taugs**2*
     &    a3pm )
      spm = spm + C2(t2)*mc**5*xn * ( taucg**(-2)*taugs**2*a3pm )
      spm = spm + C2(t2b)*mc*xn * ( taucs*taucg**(-1)*taugs*a2pm + 
     &    taucs*taucg**(-1)*taugs**2*a3pm - taucs*taucg*a3pm + taucs*
     &    taugs*a1pm + taucs*taugs*a3pm + taucs**2*taucg**(-1)*taugs*
     &    a3pm + taucs**2*taucg**(-1)*a2pm + taucs**2*a1pm + taucs**2*
     &    a3pm + taucg*taugs*a3pm - taugs*a2pm + taugs**2*a3pm )
      spm = spm + C2(t2b)*mc**3*xn * (  - taucs*taucg**(-1)*taugs*a1pm
     &     - 2.D0*taucs*taucg**(-1)*taugs*a3pm - taucs*taucg**(-1)*a2pm
     &     - 2.D0*taucs*a3pm - taucg**(-1)*taugs*a2pm - taucg**(-1)*
     &    taugs**2*a1pm + 3.D0*taugs*a3pm )
      spm = spm + C2(t2b)*mc**5*xn * ( 2.D0*taucg**(-1)*taugs*a3pm )
      spm = spm + C2(t3)*mc*xn * ( taucs*taucg**(-1)*taugs**2*a1pm - 
     &    taucs*taugs*a1pm + taugs**2*a1pm )
      spm = spm + C2(t3)*mc**3*xn * ( taucg**(-1)*taugs**2*a1pm )
      spm = spm + C2(t3b)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs*a2pm
     &     + taucg**(-1)*taugs**2*a2pm )
      spm = spm + C2(t3b)*mc**3*xn**(-1) * (  - taucg**(-1)*taugs**2*
     &    a3pm )
      spm = spm + C2(t4b)*mc**3*xn**(-1) * ( taucs*a1pm + taucg**(-1)*
     &    taugs*a2pm - taugs*a1pm + a2pm )
      spm = spm + C2(t4b)*mc**5*xn**(-1) * (  - 2.D0*taucg**(-1)*taugs*
     &    a1pm )
      spm = spm + C2(t6)*mc*xn * ( taucs*a2pm )
      spm = spm + C2(t6)*mc**3*xn * (  - taucg**(-1)*taugs*a2pm )
      spm = spm + C2(t6b)*mc*xn * ( taucs*taucg*a1pm + taucs*taugs*a1pm
     &     + taucs*taugs*a3pm + taucs*a2pm + taucs**2*a1pm - taucg*
     &    taugs*a1pm - taucg*taugs*a3pm - taucg**2*a3pm )
      spm = spm + C2(t6b)*mc**3*xn * (  - 2.D0*taucs*taucg**(-1)*taugs*
     &    a1pm - taucs*a1pm - 3.D0*taucg**(-1)*taugs*a2pm - 2.D0*
     &    taucg**(-1)*taugs**2*a1pm - 2.D0*taucg**(-1)*taugs**2*a3pm - 
     &    taugs*a3pm + a2pm )
      spm = spm + C2(t6b)*mc**5*xn * ( taucg**(-1)*taugs*a1pm )
      spm = spm + C2(t7)*mc*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pm - taucs*taucg*a1pm + taucs*taucg*a3pm + taucs*taugs*a3pm
     &     - taucs**2*taucg**(-1)*taugs*a1pm - taucs**2*a1pm - 
     &    taucg**(-1)*taugs**2*a2pm + taucg*taugs*a1pm + 2.D0*taucg*
     &    taugs*a3pm + taucg**2*a3pm - taugs*a2pm + taugs**2*a1pm + 
     &    taugs**2*a3pm )
      spm = spm + C2(t7)*mc**3*xn**(-1) * (  - taucs*a1pm - taucg**(-1)
     &    *taugs*a2pm + taucg**(-1)*taugs**2*a1pm + taucg**(-1)*
     &    taugs**2*a3pm + taugs*a1pm + taugs*a3pm - a2pm )
      spm = spm + C2(t7)*mc**5*xn**(-1) * ( taucg**(-1)*taugs*a1pm )
      spm = spm + C2(t8)*mc*xn**(-1) * ( taucs*taucg*a3pm + taucs*taugs
     &    *a1pm + taucs*taugs*a3pm + taucs*a2pm - 2.D0*taucs**2*
     &    taucg**(-1)*taugs*a1pm - 2.D0*taucs**2*taucg**(-1)*taugs*a3pm
     &     - taucs**2*taucg**(-1)*a2pm + taucs**2*a1pm - taucs**2*a3pm
     &     )
      spm = spm + C2(t8)*mc**3*xn**(-1) * (  - 3.D0*taucs*taucg**(-1)*
     &    taugs*a1pm - taucs*taucg**(-1)*taugs*a3pm + taucs*taucg**(-1)
     &    *a2pm - 2.D0*taucg**(-1)*taugs*a2pm - 2.D0*taucg**(-1)*
     &    taugs**2*a1pm - 2.D0*taucg**(-1)*taugs**2*a3pm - 2.D0*taugs*
     &    a3pm )
      spm = spm + C2(t8)*mc**5*xn**(-1) * (  - 2.D0*taucg**(-1)*taugs*
     &    a3pm )
      spm = spm + C2(t8b)*mc*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*taugs*
     &    a2pm + taucs*taucg**(-1)*taugs**2*a1pm + taucs*taucg**(-1)*
     &    taugs**2*a3pm + taucs*taugs*a3pm + 3.D0*taucs**2*taucg**(-1)*
     &    taugs*a1pm + taucs**2*taucg**(-1)*taugs*a3pm + taucs**2*
     &    taucg**(-1)*a2pm + taucs**3*taucg**(-1)*a1pm )
      spm = spm + C2(t8b)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a3pm )
      spm = spm + C2(t9)*mc*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a1pm + taucs*taucg**(-1)*taugs**2*a3pm + taucs*taugs
     &    *a1pm + taucs*taugs*a3pm - taucs*a2pm + taucs**2*taucg**(-1)*
     &    taugs*a1pm - taucs**2*a1pm + taucg**(-1)*taugs**2*a2pm + 
     &    taucg**(-1)*taugs**3*a1pm + taucg**(-1)*taugs**3*a3pm + taucg
     &    *taugs*a3pm + taugs*a2pm + taugs**2*a1pm + 2.D0*taugs**2*a3pm
     &     )
      spm = spm + C2(t9)*mc**3*xn**(-1) * ( 3.D0*taucs*taucg**(-1)*
     &    taugs*a1pm - taucs*a1pm + taucg**(-1)*taugs*a2pm + 3.D0*
     &    taucg**(-1)*taugs**2*a1pm + 2.D0*taucg**(-1)*taugs**2*a3pm + 
     &    taugs*a1pm + 2.D0*taugs*a3pm - a2pm )
      spm = spm + C2(t9)*mc**5*xn**(-1) * ( 2.D0*taucg**(-1)*taugs*a1pm
     &     )
      spm = spm + D0(b2)*mc*xn**(-1) * (  - taucs**2*taucg*a1pm + 
     &    taucs**2*a2pm )
      spm = spm + D0(b2)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pm + 2.D0*taucs*taugs*a1pm )
      spm = spm + D0(b2)*mc**5*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1pm - taucg**(-1)*taugs*a2pm - taucg**(-1)*taugs**2*a3pm - 
     &    taugs*a3pm )
      spm = spm + D1(b1)*mc*xn**(-1) * (  - taucs**2*taucg**(-1)*
     &    taugs**2*a1pm )
      spm = spm + D1(b2)*mc*xn**(-1) * (  - taucs**2*taucg*a1pm )
      spm = spm + D1(b2)*mc**3*xn**(-1) * ( 4.D0*taucs*taugs*a1pm + 
     &    taucg**(-1)*taugs**2*a2pm + taugs*a2pm )
      spm = spm + D1(b3)*mc*xn * (  - taucs*taucg*taugs*a3pm - taucs*
     &    taugs**2*a3pm - taugs**2*a2pm )
      spm = spm + D1(b3)*mc**3*xn * (  - taucg**(-1)*taugs**2*a2pm + 
     &    taugs**2*a3pm )
      spm = spm + D2(b1)*mc*xn**(-1) * ( taucs*taucg**(-1)*taugs**3*
     &    a1pm + taucs*taucg**(-1)*taugs**3*a3pm + taucs*taugs**2*a3pm
     &     - taucs**2*taucg**(-1)*taugs*a2pm + taucs**2*taucg**(-1)*
     &    taugs**2*a1pm )
      spm = spm + D2(b1)*mc**3*xn**(-1) * ( taucs*taucg**(-1)*taugs**2*
     &    a3pm + taugs**2*a3pm )
      spm = spm + D2(b2)*mc*xn**(-1) * (  - taucs*taucg*taugs*a1pm + 
     &    taucs*taugs*a2pm - taucs**2*taucg*a1pm + taucs**2*a2pm )
      spm = spm + D2(b2)*mc**3*xn**(-1) * (  - taucs*taucg*a3pm - 2.D0*
     &    taucs*taugs*a1pm - taucs*taugs*a3pm - taucg**(-1)*taugs**2*
     &    a2pm - taucg*taugs*a3pm - taugs*a2pm - taugs**2*a3pm )
      spm = spm + D2(b2)*mc**5*xn**(-1) * ( taucs*taucg**(-1)*taugs*
     &    a1pm - taucg**(-1)*taugs*a2pm - taucg**(-1)*taugs**2*a3pm - 
     &    taugs*a3pm )
      spm = spm + D2(b3)*mc*xn * ( taucs*taucg*taugs*a1pm + taucs*taucg
     &    *taugs*a3pm - taucs*taugs**2*a1pm + taucs*taugs**2*a3pm - 
     &    taucg*taugs**2*a1pm + taugs**2*a2pm )
      spm = spm + D2(b3)*mc**3*xn * ( taucg**(-1)*taugs**3*a1pm - 
     &    taugs**2*a1pm - taugs**2*a3pm )
      spm = spm + D3(b1)*mc*xn**(-1) * ( taucs*taucg*taugs*a3pm + taucs
     &    *taugs**2*a1pm + taucs*taugs**2*a3pm + taucs**2*taucg**(-1)*
     &    taugs*a2pm + taucs**2*taucg**(-1)*taugs**2*a1pm + taucs**2*
     &    taucg**(-1)*taugs**2*a3pm - taucs**2*taugs*a1pm + taucs**2*
     &    taugs*a3pm - taucs**2*a2pm - taucs**3*a1pm )
      spm = spm + D3(b1)*mc**3*xn**(-1) * ( 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a1pm + taucs*taugs*a3pm + taucs**2*taucg**(-1)*taugs
     &    *a1pm + taucg*taugs*a3pm )
      spm = spm + D3(b2)*mc*xn**(-1) * (  - taucs*taugs*a2pm - taucs*
     &    taugs**2*a1pm - 2.D0*taucs**2*taugs*a1pm - taucs**3*a1pm )
      spm = spm + D3(b2)*mc**3*xn**(-1) * (  - taucs*taucg**(-1)*taugs*
     &    a2pm - taucs*taucg**(-1)*taugs**2*a3pm - taucs*taugs*a3pm - 
     &    taucg**(-1)*taugs**3*a3pm - taugs**2*a3pm )
      spm = spm + D3(b3)*mc*xn * (  - taucs*taucg*taugs*a1pm - taucs*
     &    taucg**2*a3pm + taucs**2*taucg*a1pm - taucs**2*taucg*a3pm - 
     &    taucs**2*taugs*a1pm - taucs**2*taugs*a3pm - taucs**2*a2pm + 
     &    taucg*taugs**2*a1pm + taucg*taugs**2*a3pm + taucg**2*taugs*
     &    a3pm )
      spm = spm + D3(b3)*mc**3*xn * ( taucs*taucg**(-1)*taugs*a2pm + 
     &    taucs*taucg**(-1)*taugs**2*a1pm + 2.D0*taucs*taucg**(-1)*
     &    taugs**2*a3pm - 2.D0*taucs*taugs*a1pm + 3.D0*taucs*taugs*a3pm
     &     + taucg**(-1)*taugs**2*a2pm + taucg*taugs*a3pm - taugs*a2pm
     &     + taugs**2*a1pm - taugs**2*a3pm )
      spm = spm + D3(b3)*mc**5*xn * ( taucg**(-1)*taugs**2*a1pm - 2.D0*
     &    taucg**(-1)*taugs**2*a3pm )

      lopp=  + mc * (  - taucs*taucg**(-1)*a2pp + taucs*a1pp )
      lopp = lopp + mc**3 * ( taucg**(-2)*taugs*a2pp - taucg**(-1)*
     &    taugs*a1pp )

      lomm=  + mc**2 * (  - taucg**(-1)*taugs*a3mm )
      lomm = lomm + taucs*a3mm - taucg*taugs*a2mm + taugs*a3mm

      lomp=  + mc**2 * (  - taucg**(-2)*taugs*a3mp + taucg**(-1)*taugs*
     &    a1mp )
      lomp = lomp + taucs*taucg**(-1)*a3mp + taugs*a1mp + taugs*a2mp + 
     &    a3mp

      lopm=  + mc * ( taucs*taucg**(-1)*a2pm )
      lopm = lopm + mc**3 * (  - taucg**(-2)*taugs*a2pm )

c      lomp=za(ie,ic)
c     . /za(ig,is)/za(ig,ic)*(za(is,ic)*zb(is,in)+za(ig,ic)*zb(ig,in)) 
c      lopp=-mc/za(ig,ic)*za(ig,ie)
c     . /za(ig,is)/za(ig,ic)*(za(is,ic)*zb(is,in)+za(ig,ic)*zb(ig,in))
c      lomm=-(za(ig,ie)*zb(ig,is)+za(ie,ic)*zb(is,ic))*zb(is,in)
c     . /zb(ig,is)/zb(ig,ic)
c      lopm =-mc/taucg*za(ig,ie)*zb(is,in)*zb(is,ic)/zb(ig,is)
c      write(6,*) 'lomm,lopp,lomp,lopm',lomm,lopp,lomp,lopm
      prop=(dble(za(ie,in)*zb(in,ie))-wmass**2)**2+(wmass*wwidth)**2 
c      virtsqwcg=+abs(lomm)**2+abs(lopp)**2+abs(lomp)**2+abs(lopm)**2
c      write(6,*) virtsqwcg,temp,virtsqwcg/temp
      virtsqwcg=dble(smm*Dconjg(lomm)+spp*Dconjg(lopp)
     .              +smp*Dconjg(lomp)+spm*Dconjg(lopm))
      virtsqwcg=virtsqwcg/prop
      return
      end
