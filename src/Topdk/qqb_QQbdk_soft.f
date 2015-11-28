      subroutine qqb_QQbdk_soft(p,msq)
      implicit none
c----Matrix element for tt~g production
c----in the soft approximation
C----averaged over initial colours and spins
C***********************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2002.                                                      *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     q(-p1) +q(-p2)=t(nu(p3)+e^+(p4)+b(p5))                           *
*                   +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)                         *
*     t and tbar are assumed to be on shell                            * 
************************************************************************

      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'msq_cs.f'
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),ps(mxpart,4)
      double precision wtgg,eik(5,5),feikqqb,feikqbq,p3Dp3,p5Dp5

C---set all matrix elements equal to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      do j=1,2
      ps(j,nu)=p(j,nu)
      enddo
      do j=3,mxpart
      ps(j,nu)=zip
      enddo
      enddo
      
      do nu=1,4
c t momentum
      ps(3,nu)=(p(3,nu)+p(4,nu)+p(5,nu))
c positron
      ps(4,nu)=p(4,nu)
c t-bar momentum
      ps(5,nu)=(p(6,nu)+p(7,nu)+p(8,nu))
c electron
      ps(6,nu)=p(7,nu)
c gluon
      ps(7,nu)=p(9,nu)
      enddo      

      p3Dp3=ps(3,4)**2-ps(3,1)**2-ps(3,2)**2-ps(3,3)**2
      p5Dp5=ps(5,4)**2-ps(5,1)**2-ps(5,2)**2-ps(5,3)**2

      call dotem(7,ps,s)

      eik(1,2)=+2d0*gsq*s(1,2)/(s(1,7)*s(2,7))
      eik(1,3)=+2d0*gsq*s(1,3)/(s(1,7)*s(3,7))
      eik(1,5)=+2d0*gsq*s(1,5)/(s(1,7)*s(5,7))
      eik(2,3)=+2d0*gsq*s(2,3)/(s(2,7)*s(3,7))
      eik(2,5)=+2d0*gsq*s(2,5)/(s(2,7)*s(5,7))
      eik(3,5)=+2d0*gsq*s(3,5)/(s(3,7)*s(5,7))

      eik(1,1)=0d0
      eik(2,2)=0d0
      eik(3,3)=gsq*p3Dp3*(2d0/s(3,7))**2
      eik(5,5)=gsq*p5dp5*(2d0/s(5,7))**2


      call qqb_QQbdk(p,msq)

      wtgg=
     .  + msq_cs(1,0,0)*xn*(eik(1,2)+eik(1,5)+eik(2,3)+eik(3,5))
     .  + msq_cs(2,0,0)*xn*(eik(1,2)+eik(1,3)+eik(2,5)+eik(3,5))
     .  + msq(0,0)*xn*(-half*eik(3,3)-half*eik(5,5)-eik(3,5))
     .  + msq(0,0)/xn*(half*eik(3,3)+half*eik(5,5)-eik(3,5))
     .  + msq_cs(0,0,0)*xn*(eik(1,3)+eik(1,5)+eik(2,3)+eik(2,5)) 

      feikqqb=CF*(2d0*eik(1,3)+2d0*eik(2,5)-eik(3,3)-eik(5,5))
     . +(2d0*eik(1,5)+2d0*eik(2,3)
     . -eik(1,3)-eik(2,5)-eik(1,2)-eik(3,5))/xn
      feikqbq=CF*(2d0*eik(2,3)+2d0*eik(1,5)-eik(3,3)-eik(5,5))
     . +(2d0*eik(2,5)+2d0*eik(1,3)
     . -eik(2,3)-eik(1,5)-eik(1,2)-eik(3,5))/xn



      do j=-nf,nf
      k=-j      
      if (j .lt.0) then
      msq(j,k)=feikqbq*msq(-1,1)
      elseif (j.gt.0) then
      msq(j,k)=feikqqb*msq(1,-1)
      elseif (j.eq.0) then
      msq(0,0)=wtgg
      endif
      enddo
      return
      end
 
