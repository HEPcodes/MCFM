      subroutine qqb_wgam(p,msq)
      implicit none
c----Matrix element for W gam production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5)
C For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'ckm.f'
      include 'prods.f'
      include 'nwz.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),qbq,qqb
      double precision ubdmsq,fac
      double complex agamtree

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(5,p,za,zb)
      fac=aveqq*2d0*xn*gwsq**2*esq
      
      if (nwz .eq. -1) then 
         qbq=fac*(abs(agamtree(1,2,3,4,5,za,zb,-1))**2
     .           +abs(agamtree(1,2,3,4,5,za,zb,+1))**2)
         qqb=fac*(abs(agamtree(2,1,3,4,5,za,zb,-1))**2
     .           +abs(agamtree(2,1,3,4,5,za,zb,+1))**2)
      elseif (nwz .eq. +1) then 
         qbq=fac*(abs(agamtree(2,1,4,3,5,zb,za,-1))**2
     .           +abs(agamtree(2,1,4,3,5,zb,za,+1))**2)
         qqb=fac*(abs(agamtree(1,2,4,3,5,zb,za,-1))**2
     .           +abs(agamtree(1,2,4,3,5,zb,za,+1))**2)
      endif

      do j=-nf,nf
      do k=-nf,nf
          if ((j .gt. 0) .and. (k .lt. 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo
      return
      end

      double complex function agamtree(p1,p2,p3,p4,p5,za,zb,hgamma)
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'masses.f'
      double complex prp12,prp34
      integer p1,p2,p3,p4,p5,hgamma
C  -i * Matrix element is what is implemented 
      prp12=dcmplx((s(p1,p2)-wmass**2),wmass*wwidth)
      prp34=dcmplx((s(p3,p4)-wmass**2),wmass*wwidth)
      if (hgamma .eq. -1) then
      agamtree=Qu*za(p2,p5)*za(p3,p4)/zb(p1,p5)/prp34
     .  +Qd*za(p1,p5)*za(p3,p4)/zb(p2,p5)/prp34
     .  +(Qu-Qd)*za(p1,p2)*za(p4,p5)/zb(p3,p5)/prp12
      agamtree=-zb(p2,p4)**2*agamtree/(s(p1,p2)-s(p3,p4))
      elseif (hgamma .eq. +1) then
      agamtree=Qu*zb(p2,p5)*zb(p3,p4)/za(p1,p5)/prp34
     .  +Qd*zb(p1,p5)*zb(p3,p4)/za(p2,p5)/prp34
     .  +(Qu-Qd)*zb(p1,p2)*zb(p4,p5)/za(p3,p5)/prp12
      agamtree=agamtree*za(p1,p3)**2/(s(p1,p2)-s(p3,p4))
      endif
      return
      end

      double precision function ubdmsq(i1,i2,i3,i4,i5)
      implicit none
C     Matrix element for 
C     ub(-p1)+d(-p2)=e-(p3)+nu~(p4)+gamma(p5)
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'ckm.f'
      include 'prods.f'
      integer i1,i2,i3,i4,i5
      double precision fac,propsq12,propsq34,propsq1234 
c     ans=
c      +2*xn*gw^4*ee^2*(s24^2+s13^2)*[2*p1.p5+2*p2.p5]^-2*( 
c       s12*[Qu-Qd]^2*s45/s35/[s12-mw^2]^2
c
c       +s34*(Qu/s15+Qd/s25)^2*s15*s25/[s34-mw^2]^2
c
C       +((s23*s45-s24*s35)*s15-(s13*s45-s14*s35)*s25)
C       *[Qu-Qd]*(Qu/s15+Qd/s25)/s35/[s12-mw^2]/[s34-mw^2])
      propsq12=(s(i1,i2)-wmass**2)**2+(wmass*wwidth)**2
      propsq34=(s(i3,i4)-wmass**2)**2+(wmass*wwidth)**2
      propsq1234=
     . (s(i1,i2)-wmass**2)*(s(i3,i4)-wmass**2)+(wmass*wwidth)**2
      ubdmsq=
     . +s(i1,i2)*(Qu-Qd)**2*s(i4,i5)/(s(i3,i5)*propsq12)
     .     +s(i3,i4)*(Qu/s(i1,i5)+Qd/s(i2,i5))**2
     .      *s(i1,i5)*s(i2,i5)/propsq34
     . +((s(i2,i3)*s(i4,i5)-s(i2,i4)*s(i3,i5))*s(i1,i5)
     .  -(s(i1,i3)*s(i4,i5)-s(i1,i4)*s(i3,i5))*s(i2,i5))
     . *(Qu-Qd)*(Qu/s(i1,i5)+Qd/s(i2,i5))
     . /(s(i3,i5)*propsq1234)
CCCC Epsilon piece missing! 

      ubdmsq=(s(i2,i4)**2+s(i1,i3)**2)/(s(i1,i5)+s(i2,i5))**2*ubdmsq
      fac=aveqq*2d0*xn*gwsq**2*esq
      ubdmsq=fac*ubdmsq
      return
      end
