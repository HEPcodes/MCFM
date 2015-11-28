      subroutine qqb_ww_v(p,msqv)
      implicit none
C----Author R.K.Ellis December 1998
C----modified by JMC to include supplementary diagrams February 1999
c----Matrix element for WW production
c----in the notation of DKS
C----averaged over initial colours and spins
C----massless b quarks
c     q(-p1)+qbar(-p2)-->q'(p5)+bar{q'}(p6)+n(p3)+ebar(p4)
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'dprodx.f'
      include 'sprodx.f'
      include 'epinv.f'
      include 'zerowidth.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'scale.f'
      logical msbar
      common/msbar/msbar
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),mp(nf),
     . p(mxpart,4),qdks(mxpart,4),facnlo,sub,ave,xl12,virt
      double complex AWWM,AWWP,BWWM,BWWP
      double complex prop12,prop34,prop56
      double complex a6treea,a6treeb,a6loopa,a6loopb
      double complex propwp,propwm,propzg,cprop,Vpole,Vpole12,suppl
      double complex Fa123456,Fa213456,Fb123456,Fb213456
      double complex Fa126543,Fa216543,Fb126543,Fb216543
      double complex La123456,La213456,Lb123456,Lb213456
      double complex La126543,La216543,Lb126543,Lb216543
      double complex Fa341256,Fa653421,Fa346521,Fa651243
      double complex Fa342156,Fa653412,Fa346512,Fa652143
      double complex cl(2),cr(2),clgamz(2),crgamz(2),clz(2),crz(2)
      double precision FAC
      integer j,k,jk
      parameter(ave=0.25d0/xn)
      data mp/-1d0,+1d0,-1d0,+1d0,-1d0/

      FAC=gw**8*ave
      facnlo=ason2pi*cf

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msqv(j,k)=0d0
      enddo
      enddo

c---calculate the lowest order matrix element
      call qqb_ww(p,msq)

c---add result of integrating subtraction terms
      xl12=log(s(1,2)/musq)
      sub = epinv**2+epinv*(1.5d0-xl12)
     . +half*xl12**2-pisqo6-0.5d0
      if (msbar) sub=sub+1d0
c---note extra colour factor
      sub=sub*ason2pi*two*cf


C----Change the momenta to DKS notation 
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)

      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(5,j)
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(3,j)
      enddo


c--   s returned from sprod (common block) is 2*dot product
      call spinoru(6,qdks,za,zb)

c--   calculate propagators

      if     (zerowidth  .eqv. .true.) then
      prop12=s(1,2)/(s(1,2)-zmass**2+im*zmass*zwidth)
      prop34=s(3,4)/(s(3,4)-wmass**2+im*wmass*wwidth)
      prop56=s(5,6)/(s(5,6)-wmass**2+im*wmass*wwidth)
      cprop=dcmplx(1d0)
      elseif (zerowidth .neqv. .true.) then
      prop12=dcmplx(s(1,2)/(s(1,2)-zmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-wmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-wmass**2))
      propwm=(s(3,4)-wmass**2)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
      propwp=(s(5,6)-wmass**2)/(s(5,6)-wmass**2+im*wmass*wwidth)
      propzg=(s(1,2)-zmass**2)/(s(1,2)-zmass**2+im*zmass*zwidth)
      cprop=propwp*propwm*propzg
      endif
      
c-- couplings with or without photon pole
      do j=1,2
      cl(j)=mp(j)*(two*Q(j)*xw+l(j)*sin2w*prop12)
      cr(j)=mp(j)*two*Q(j)*xw*(1d0-prop12)
      if (zerowidth .neqv. .true.) then
      clgamz(j)=two*xw*(-Q(j)+le*L(j)*prop12)
      crgamz(j)=two*xw*(-Q(j)+le*R(j)*prop12)
      clz(j)=two*xw*ln*L(j)*prop12
      crz(j)=two*xw*ln*R(j)*prop12
      endif
      enddo

c---case dbar-d and d-dbar
      Fa126543=A6treea(1,2,6,5,4,3,za,zb)
      Fa216543=A6treea(2,1,6,5,4,3,za,zb)
      Fa123456=A6treea(1,2,3,4,5,6,za,zb)
      Fa213456=A6treea(2,1,3,4,5,6,za,zb)

      Fb123456=A6treeb(1,2,3,4,5,6,za,zb)
      Fb126543=-Fb123456
      Fb213456=A6treeb(2,1,3,4,5,6,za,zb)
      Fb216543=-Fb213456

      La126543=A6loopa(1,2,6,5,4,3,za,zb)
      La216543=A6loopa(2,1,6,5,4,3,za,zb)
      La123456=A6loopa(1,2,3,4,5,6,za,zb)
      La213456=A6loopa(2,1,3,4,5,6,za,zb)

      Lb123456=A6loopb(1,2,3,4,5,6,za,zb)
      Lb213456=A6loopb(2,1,3,4,5,6,za,zb)
      Lb126543=A6loopb(1,2,6,5,4,3,za,zb)
      Lb216543=A6loopb(2,1,6,5,4,3,za,zb)

      if (zerowidth .neqv. .true.) then
c---for supplementary diagrams.
      Fa341256=A6treea(3,4,1,2,5,6,za,zb)
      Fa653421=A6treea(6,5,3,4,2,1,za,zb)
      Fa346521=A6treea(3,4,6,5,2,1,za,zb)
      Fa651243=A6treea(6,5,1,2,4,3,za,zb)
      Fa342156=A6treea(3,4,2,1,5,6,za,zb)
      Fa653412=A6treea(6,5,3,4,1,2,za,zb)
      Fa346512=A6treea(3,4,6,5,1,2,za,zb)
      Fa652143=A6treea(6,5,2,1,4,3,za,zb)
c---loop diagrams just tree*Vpole since they're all triangle-type
      Vpole12=Vpole(s(1,2))
      endif

      do j=-nf,nf
      k=-j
      if (j.eq.0) go to 20
      jk=max(j,k)

c--assign values
      if (j .gt. 0) then
         if (tau(jk) .eq. +1d0) then
           AWWM=(Fa213456+cl(2)*Fb213456)*prop56*prop34
           AWWP=cr(2)*Fb123456*prop56*prop34
           BWWM=(La213456+cl(2)*Lb213456)*prop56*prop34
           BWWP=cr(2)*Lb123456*prop56*prop34
         elseif (tau(jk) .eq. -1d0) then
           AWWM=(Fa216543+cl(1)*Fb216543)*prop56*prop34
           AWWP=cr(1)*Fb126543*prop56*prop34
           BWWM=(La216543+cl(1)*Lb216543)*prop56*prop34
           BWWP=cr(1)*Lb126543*prop56*prop34
	 endif
      elseif (j .lt. 0) then
         if (tau(jk) .eq. +1d0) then
C---remember that this is the basic process, M refers to the helicity
c---of the outgoing quark
           AWWM=(Fa123456+cl(2)*Fb123456)*prop56*prop34
           AWWP=cr(2)*Fb213456*prop56*prop34
           BWWM=(La123456+cl(2)*Lb123456)*prop56*prop34
           BWWP=cr(2)*Lb213456*prop56*prop34
         elseif (tau(jk) .eq. -1d0) then
           AWWM=(Fa126543+cl(1)*Fb126543)*prop56*prop34
           AWWP=cr(1)*Fb216543*prop56*prop34
           BWWM=(La126543+cl(1)*Lb126543)*prop56*prop34
           BWWP=cr(1)*Lb216543*prop56*prop34
         endif
      endif

      if (zerowidth .neqv. .true.) then
c---we need supplementary diagrams for gauge invariance.
      if     (j .gt. 0) then
          if     (tau(jk) .eq. +1d0) then
                   suppl=clgamz(2)*(Fa342156*prop56+Fa653412*prop34)
     .                  +clz(2)*(Fa346512*prop56+Fa652143*prop34)
                   AWWM=AWWM+suppl
                   BWWM=BWWM+suppl*Vpole12
                   suppl=crgamz(2)*(Fa341256*prop56+Fa653421*prop34)
     .                  +crz(2)*(Fa346521*prop56+Fa651243*prop34)
                   AWWP=AWWP+suppl
                   BWWP=BWWP+suppl*Vpole12
          elseif (tau(jk) .eq. -1d0) then
                   suppl=clgamz(1)*(Fa342156*prop56+Fa653412*prop34)
     .                  +clz(1)*(Fa346512*prop56+Fa652143*prop34)
c                   suppl=clgamz(1)*(Fa652143*prop34+Fa346512*prop56)
c     .                  +clz(1)*(Fa653412*prop34+Fa342156*prop56)
                   AWWM=AWWM+suppl
                   BWWM=BWWM+suppl*Vpole12
                   suppl=crgamz(1)*(Fa341256*prop56+Fa653421*prop34)
     .                  +crz(1)*(Fa346521*prop56+Fa651243*prop34)
c                   suppl=crgamz(1)*(Fa651243*prop34+Fa346521*prop56)
c     .                  +crz(1)*(Fa653421*prop34+Fa341256*prop56)
                   AWWP=AWWP+suppl
                   BWWP=BWWP+suppl*Vpole12
         endif
      elseif (j .lt. 0) then
          if     (tau(jk) .eq. +1d0) then
                   suppl=clgamz(2)*(Fa341256*prop56+Fa653421*prop34)
     .                  +clz(2)*(Fa346521*prop56+Fa651243*prop34)
                   AWWM=AWWM+suppl
                   BWWM=BWWM+suppl*Vpole12
                   suppl=crgamz(2)*(Fa342156*prop56+Fa653412*prop34)
     .                  +crz(2)*(Fa346512*prop56+Fa652143*prop34)
                   AWWP=AWWP+suppl
                   BWWP=BWWP+suppl*Vpole12
          elseif (tau(jk) .eq. -1d0) then
                   suppl=clgamz(1)*(Fa341256*prop56+Fa653421*prop34)
     .                  +clz(1)*(Fa346521*prop56+Fa651243*prop34)
c                   suppl=clgamz(1)*(Fa651243*prop34+Fa346521*prop56)
c     .                  +clz(1)*(Fa653421*prop34+Fa341256*prop56)
                   AWWM=AWWM+suppl
                   BWWM=BWWM+suppl*Vpole12
                   suppl=crgamz(1)*(Fa342156*prop56+Fa653412*prop34)
     .                  +crz(1)*(Fa346512*prop56+Fa652143*prop34)
c                   suppl=crgamz(1)*(Fa652143*prop34+Fa346512*prop56)
c     .                  +crz(1)*(Fa653412*prop34+Fa342156*prop56)
                   AWWP=AWWP+suppl
                   BWWP=BWWP+suppl*Vpole12
          endif
      endif
      endif

C-- Inclusion of width for W's a la Baur and Zeppenfeld
c-- also include coupling in l1 to account for non-leptonic W decays
      AWWM=l1*cprop*AWWM
      AWWP=l1*cprop*AWWP
      BWWM=l1*cprop*BWWM
      BWWP=l1*cprop*BWWP
      
      virt=facnlo*fac*two*dble(dconjg(AWWM)*BWWM+dconjg(AWWP)*BWWP)

      msqv(j,k)=sub*msq(j,k)+virt

 20   continue
      enddo
	
      return
      end
