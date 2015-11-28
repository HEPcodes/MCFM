      subroutine qqb_wgam_g(p,msq)
      implicit none
c----Matrix element for W gam production
C----averaged over initial colours and spins
C For nwz=+1
c     u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + g(p6)
C For nwz=-1
c     ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + g(p6)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'nwz.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision qbq,qqb,qg,gq,qbg,gqb
      double precision ubdgmsq


      call spinoru(6,p,za,zb)
      fac=V/2d0*gwsq**2*4d0*gsq*esq



      if (nwz .eq. -1) then 
      qbq=aveqq*fac*ubdgmsq(1,2,3,4,5,6,za,zb)
      qqb=aveqq*fac*ubdgmsq(2,1,3,4,5,6,za,zb)
      gq=aveqg*fac*ubdgmsq(6,2,3,4,5,1,za,zb)
      qg=aveqg*fac*ubdgmsq(6,1,3,4,5,2,za,zb)
      gqb=aveqg*fac*ubdgmsq(2,6,3,4,5,1,za,zb)
      qbg=aveqg*fac*ubdgmsq(1,6,3,4,5,2,za,zb)

      elseif (nwz .eq. +1) then 
      qbq=aveqq*fac*ubdgmsq(2,1,4,3,5,6,zb,za)
      qqb=aveqq*fac*ubdgmsq(1,2,4,3,5,6,zb,za)
      gq=aveqg*fac*ubdgmsq(2,6,4,3,5,1,zb,za)
      qg=aveqg*fac*ubdgmsq(1,6,4,3,5,2,zb,za)
      gqb=aveqg*fac*ubdgmsq(6,2,4,3,5,1,zb,za)
      qbg=aveqg*fac*ubdgmsq(6,1,4,3,5,2,zb,za)

      endif


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
          if ((j .gt. 0) .and. (k .lt. 0)) then
            msq(j,k)=Vsq(j,k)*qqb
          elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=Vsum(k)*gqb
          elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=Vsum(k)*gq
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            msq(j,k)=Vsq(j,k)*qbq
          elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Vsum(j)*qg
         elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=Vsum(j)*qbg
          endif
      enddo
      enddo
      return
      end

      double precision function ubdgmsq(p1,p2,p3,p4,p5,p6,za,zb)
      implicit none
C     Matrix element for 
C     ub(-p1)+d(-p2)=e-(p3)+nu~(p4)+gamma(p5)+g(p6)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      integer p1,p2,p3,p4,p5,p6
      double complex aLL,aRR,aRL,aLR,prp34,prp345
      double precision s156,s256

      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)
      s256=s(p2,p5)+s(p2,p6)+s(p5,p6)
      prp34=dcmplx(s(p3,p4)-wmass**2,wmass*wwidth)
      prp345=dcmplx(s(p3,p4)+s(p3,p5)+s(p4,p5)-wmass**2,wmass*wwidth)

      if (zerowidth) then
      prp345=dcmplx(s(p3,p5)+s(p4,p5),0d0)



C  Eqs. 4.9-4.12 Multiplied by i 
      aLR=+Qd/zb(p2,p5)*(
     .   +za(p1,p3)*zb(p6,p2)*za(p1,p6)/s256
     .  *(za(p5,p1)*zb(p1,p4)+za(p5,p3)*zb(p3,p4))

     . +(za(p1,p2)*zb(p2,p4)+za(p1,p6)*zb(p6,p4))/prp345
     . *((za(p3,p4)*za(p1,p5)*zb(p4,p2)+za(p3,p5)*za(p1,p6)*zb(p6,p2))
     . ))
 
     .    -Qu/zb(p1,p5)*(
     .  +za(p1,p5)*zb(p2,p4)*za(p6,p2)/s156
     .  *(za(p3,p1)*zb(p1,p6)+za(p3,p5)*zb(p5,p6))

     .  +(za(p1,p2)*zb(p2,p4)+za(p1,p6)*zb(p6,p4))
     . *((za(p5,p2)*zb(p2,p4)+za(p5,p6)*zb(p6,p4))*za(p4,p3)
     . +za(p5,p3)*za(p2,p6)*zb(p6,p2))/prp345)
      aLR=aLR/(prp34*za(p1,p6)*za(p6,p2))

      aRL=+Qu/za(p1,p5)*(
     .   +zb(p2,p4)*za(p6,p1)*zb(p2,p6)/s156
     .  *(zb(p5,p2)*za(p2,p3)+zb(p5,p4)*za(p4,p3))

     . +(zb(p2,p1)*za(p1,p3)+zb(p2,p6)*za(p6,p3))/prp345
     . *((zb(p4,p3)*zb(p2,p5)*za(p3,p1)+zb(p4,p5)*zb(p2,p6)*za(p6,p1))
     . ))
 
     .    -Qd/za(p2,p5)*(
     .  +zb(p2,p5)*za(p1,p3)*zb(p6,p1)/s256
     .  *(zb(p4,p2)*za(p2,p6)+zb(p4,p5)*za(p5,p6))

     .  +(zb(p2,p1)*za(p1,p3)+zb(p2,p6)*za(p6,p3))
     . *((zb(p5,p1)*za(p1,p3)+zb(p5,p6)*za(p6,p3))*zb(p3,p4)
     . +zb(p5,p4)*zb(p1,p6)*za(p6,p1))/prp345)
      aRL=aRL/(prp34*zb(p2,p6)*zb(p6,p1))


      aRR=+Qu/za(p1,p5)*(za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5))
     .    -Qd/za(p2,p5)*(za(p2,p3)*zb(p3,p5)+za(p2,p4)*zb(p4,p5))
      aRR=aRR*za(p1,p3)**2*zb(p4,p3)/(prp34*prp345*za(p1,p6)*za(p6,p2))

      aLL=+Qd/zb(p2,p5)*(zb(p2,p4)*za(p4,p5)+zb(p2,p3)*za(p3,p5))
     .    -Qu/zb(p1,p5)*(zb(p1,p4)*za(p4,p5)+zb(p1,p3)*za(p3,p5))
      aLL=aLL*zb(p2,p4)**2*za(p3,p4)/(prp34*prp345*zb(p1,p6)*zb(p6,p2))

      else
c      prp345=dcmplx(s(p3,p4)+s(p3,p5)+s(p4,p5)-wmass**2,wmass*wwidth)
      aLR=+Qd/zb(p2,p5)*(
     .   +za(p1,p3)*za(p1,p6)*zb(p2,p6)/prp34/s256
     .  *(za(p1,p5)*zb(p1,p4)+za(p3,p5)*zb(p3,p4))

     .  +(za(p1,p2)*zb(p2,p4)-za(p1,p6)*zb(p4,p6))/prp345/prp34
     .  *((-za(p1,p5)*za(p3,p4)*zb(p2,p4)-za(p1,p6)*za(p3,p5)*zb(p2,p6))
     .  -prp34/zb(p3,p5)
     .  *(za(p1,p3)*zb(p2,p3)+za(p1,p5)*zb(p2,p5))))
 
     .    -Qu/zb(p1,p5)*(
     .  +za(p1,p5)*za(p2,p6)*zb(p2,p4)/prp34/s156
     .  *(za(p1,p3)*zb(p1,p6)-za(p3,p5)*zb(p5,p6))
     .  +(za(p1,p2)*zb(p2,p4)-za(p1,p6)*zb(p4,p6))/prp345/prp34
     . *(za(p3,p5)*(-za(p2,p6)*zb(p6,p2)-za(p3,p4)*zb(p3,p4))
     .  -za(p1,p5)*za(p3,p4)*zb(p1,p4)
     .  -prp34/zb(p3,p5)
     .  *(za(p1,p5)*zb(p1,p5)+za(p1,p3)*zb(p1,p3)+za(p3,p5)*zb(p3,p5))))
      aLR=aLR/(za(p1,p6)*za(p6,p2))

 

      aRL=
     .      +Qu/za(p1,p5)*(
     .  +zb(p2,p4)*zb(p2,p6)*za(p1,p6)/prp34/s156
     .  *(zb(p2,p5)*za(p2,p3)+zb(p4,p5)*za(p4,p3))
     . +(zb(p2,p1)*za(p1,p3)-zb(p2,p6)*za(p3,p6))/prp345/prp34
     . *((zb(p2,p5)*zb(p3,p4)*za(p1,p3)-zb(p2,p6)*zb(p4,p5)*za(p1,p6))
     .    -prp34/za(p3,p5)*za(p1,p3)*zb(p2,p4)))

     .    -Qd/za(p2,p5)*(
     .   +zb(p2,p5)*zb(p1,p6)*za(p1,p3)/prp34/s256
     .  *(zb(p2,p4)*za(p2,p6)-zb(p4,p5)*za(p5,p6))
     .  +(zb(p2,p1)*za(p1,p3)-zb(p2,p6)*za(p3,p6))/prp345/prp34
     . *(zb(p4,p5)*(-zb(p1,p6)*za(p6,p1)-zb(p3,p4)*za(p3,p4))
     .     -zb(p2,p5)*zb(p4,p3)*za(p2,p3)
     .   -prp34/za(p3,p5)
     .  *(za(p2,p3)*zb(p2,p4)+za(p3,p5)*zb(p4,p5))))
      aRL=aRL/(zb(p1,p6)*zb(p6,p2))

      aRR=+Qu/za(p1,p5)*(
     .    +zb(p4,p3)*(za(p1,p3)*zb(p3,p5)+za(p1,p4)*zb(p4,p5))/prp34
     .    -(za(p1,p5)*zb(p4,p5)-za(p1,p3)*zb(p3,p4))/za(p3,p5))

     .        -Qd/za(p2,p5)*(
     .    +zb(p4,p3)*(za(p2,p3)*zb(p3,p5)+za(p2,p4)*zb(p4,p5))/prp34
     .     -(za(p2,p5)*zb(p4,p5)-za(p2,p3)*zb(p3,p4))/za(p3,p5))
      aRR=aRR*za(p1,p3)**2/(prp345*za(p1,p6)*za(p6,p2))

      aLL=+Qd/zb(p2,p5)*(
     .   +za(p3,p4)*(zb(p2,p4)*za(p4,p5)+zb(p2,p3)*za(p3,p5))/prp34
     .   +(zb(p2,p5)*za(p4,p5)-zb(p2,p3)*za(p3,p4))/zb(p3,p5))

     .        -Qu/zb(p1,p5)*(
     .   +za(p3,p4)*(zb(p1,p4)*za(p4,p5)+zb(p1,p3)*za(p3,p5))/prp34
     .     +(zb(p1,p5)*za(p4,p5)-zb(p1,p3)*za(p3,p4))/zb(p3,p5))
      aLL=aLL*zb(p2,p4)**2/(prp345*zb(p1,p6)*zb(p6,p2))
      endif
      ubdgmsq=cdabs(aLL)**2+cdabs(aRR)**2+cdabs(aRL)**2+cdabs(aLR)**2
      return
      end

