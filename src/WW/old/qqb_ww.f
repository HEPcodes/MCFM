      subroutine qqb_ww(p,msq)
      implicit none
C----Author R.K.Ellis December 1998
C----modified by JMC to include supplementary diagrams February 1999
c----Matrix element for WW production
c----in the notation of DKS
C----averaged over initial colours and spins
C----massless final state particles
c     q(-p1)+qbar(-p2)-->q'(p5)+bar{q'}(p6)+n(p3)+ebar(p4)
c--- note that non-leptonic W decays do not include scattering diagrams
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'

      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4),ave
      double complex prop12,prop34,prop56
      double complex AWW(2),a6treea,a6treeb
      double complex propwp,propwm,propzg,cprop
      double complex Fa123456,Fa213456,Fb123456,Fb213456
      double complex Fa126543,Fa216543,Fb126543,Fb216543
      double complex Fa341256,Fa653421,Fa346521,Fa651243
      double complex Fa342156,Fa653412,Fa346512,Fa652143
      double complex cs(2,2),cgamz(2,2),cz(2,2)
      double precision fac,mp(nf)
      integer j,k,jk,tjk,minus,mplus
      parameter(ave=0.25d0/xn)
      data minus,mplus/1,2/
      data mp/-1d0,+1d0,-1d0,+1d0,-1d0/
      fac=gw**8*ave

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      enddo
      enddo

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

      call spinoru(6,qdks,za,zb)
c--   s returned from sprod (common block) is 2*dot product
      
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
      cs(minus,j)=mp(j)*(2d0*Q(j)*xw+l(j)*sin2w*prop12)
      cs(mplus,j)=mp(j)*2d0*Q(j)*xw*(1d0-prop12)
      if (zerowidth .neqv. .true.) then
      cz(minus,j)=2d0*xw*ln*L(j)*prop12
      cz(mplus,j)=2d0*xw*ln*R(j)*prop12
      cgamz(minus,j)=2d0*xw*(-Q(j)+le*L(j)*prop12)
      cgamz(mplus,j)=2d0*xw*(-Q(j)+le*R(j)*prop12)
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
      endif

      do j=-nf,nf
      k=-j
c--Exclude gluon-gluon initial state
      if (j.eq.0) go to 20
      jk=max(j,k)

c--assign values
c---Remember that base process is ub-u so this has the natural 123456 order
      if (j .gt. 0) then
          if         (tau(jk) .eq. +1d0) then
                AWW(minus)=(Fa213456+cs(minus,2)*Fb213456)*prop56*prop34
                AWW(mplus)=cs(mplus,2)*Fb123456*prop56*prop34
          elseif     (tau(jk) .eq. -1d0) then
                AWW(minus)=(Fa216543+cs(minus,1)*Fb216543)*prop56*prop34
                AWW(mplus)=cs(mplus,1)*Fb126543*prop56*prop34
          endif
      elseif (j .lt. 0) then
          if     (tau(jk) .eq. +1d0) then
                AWW(minus)=(Fa123456+cs(minus,2)*Fb123456)*prop56*prop34
                AWW(mplus)=cs(mplus,2)*Fb213456*prop56*prop34
          elseif (tau(jk) .eq. -1d0) then
                AWW(minus)=(Fa126543+cs(minus,1)*Fb126543)*prop56*prop34
                AWW(mplus)=cs(mplus,1)*Fb216543*prop56*prop34
          endif
      endif
     
      if (zerowidth .neqv. .true.) then
c---we need supplementary diagrams for gauge invariance.
C---tjk is equal to 2 (u,c) or 1 (d,s,b)
      tjk=2-mod(abs(jk),2)     
      if (j .gt. 0) then
          AWW(minus)=AWW(minus)
     .              +cgamz(minus,tjk)*(Fa342156*prop56+Fa653412*prop34)
     .                 +cz(minus,tjk)*(Fa346512*prop56+Fa652143*prop34)
          AWW(mplus)=AWW(mplus)
     .              +cgamz(mplus,tjk)*(Fa341256*prop56+Fa653421*prop34)
     .                 +cz(mplus,tjk)*(Fa346521*prop56+Fa651243*prop34)
      elseif (j .lt. 0) then
          AWW(minus)=AWW(minus)
     .             +cgamz(minus,tjk)*(Fa341256*prop56+Fa653421*prop34)
     .                +cz(minus,tjk)*(Fa346521*prop56+Fa651243*prop34)
          AWW(mplus)=AWW(mplus)
     .             +cgamz(mplus,tjk)*(Fa342156*prop56+Fa653412*prop34)
     .                +cz(mplus,tjk)*(Fa346512*prop56+Fa652143*prop34)
      endif
      endif
C-- Inclusion of width for W's a la Baur and Zeppenfeld with cprop.
c-- also include coupling in l1 to account for non-leptonic W decays
      AWW(minus)=l1*cprop*AWW(minus)
      AWW(mplus)=l1*cprop*AWW(mplus)

      msq(j,k)=fac*(abs(AWW(minus))**2+abs(AWW(mplus))**2)
 20   continue

      enddo
      
      return
      end


