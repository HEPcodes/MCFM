      subroutine qqb_tth(pin,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=tbar(bbar(p6)+e-(p7)+nubar(p8))
C                      +t(b(p5)+nu(p3)+e+(p4))
C                      +H(b(p9)+bbar(p10))
C  
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),pin(mxpart,4),p(mxpart,4)
      double precision pw1(4),pw2(4),p12(4),q(4),a(4),r(4),b(4),h(4),
     . x(4),y(4)
      double precision shat,sh,sw1,sw2,qDq,aDa,rDr,bDb,yDy,xDx,densq,
     . p3Dp5,p6Dp8
      double precision wtqqb,wtgg,hdecay
      double precision gamr1,gamr2,gamb1,gamb2,gamx1,gamx2,
     . gamy1,gamy2,gamq4,gama7,dot
      double precision fac,denr,denb
      double precision p1Dr,p2Dr,p1Db,p2Db,p1Dx,p2Dx,p1Dy,p2Dy,p4Dq,p7Da


      integer q4,a7,r1,r2,b1,b2
      parameter(q4=3,a7=5,r1=6,r2=8,b1=9,b2=10)

      do nu=1,4
      do j=1,mxpart
      p(j,nu)=pin(j,nu)
      enddo      
      h(nu)=p(9,nu)+p(10,nu)
      p12(nu)=p(1,nu)+p(2,nu)
      pw1(nu)=p(3,nu)+p(4,nu)
      pw2(nu)=p(7,nu)+p(8,nu)

      q(nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      r(nu)=q(nu)+h(nu)
      x(nu)=q(nu)+p(1,nu)

      a(nu)=-p(6,nu)-p(7,nu)-p(8,nu)
      b(nu)=a(nu)-h(nu)
      y(nu)=a(nu)-p(2,nu)
      enddo      


      shat=(p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2)
      sh=(h(4)**2-h(1)**2-h(2)**2-h(3)**2)
      sw1=(pw1(4)**2-pw1(1)**2-pw1(2)**2-pw1(3)**2)
      sw2=(pw2(4)**2-pw2(1)**2-pw2(2)**2-pw2(3)**2)
      qDq=(q(4)**2-q(1)**2-q(2)**2-q(3)**2)
      aDa=(a(4)**2-a(1)**2-a(2)**2-a(3)**2)
      rDr=(r(4)**2-r(1)**2-r(2)**2-r(3)**2)
      bDb=(b(4)**2-b(1)**2-b(2)**2-b(3)**2)
      xDx=(x(4)**2-x(1)**2-x(2)**2-x(3)**2)
      yDy=(y(4)**2-y(1)**2-y(2)**2-y(3)**2)



      p3Dp5=dot(p,3,5)
      p6Dp8=dot(p,6,8)

      densq=shat**2*((sw1-wmass**2)**2+wmass**2*wwidth**2)
      densq=densq*((sw2-wmass**2)**2+wmass**2*wwidth**2)
      densq=densq*((sh-hmass**2)**2+hmass**2*hwidth**2)
      densq=densq*((qDq-mt**2)**2+mt**2*twidth**2)
      densq=densq*((aDa-mt**2)**2+mt**2*twidth**2)

      denr=sqrt((rDr-mt**2)**2+mt**2*twidth**2)
      denb=sqrt((bDb-mt**2)**2+mt**2*twidth**2)

      hdecay=2d0*(sh-4d0*mb**2)*xn 
      fac=V/4d0/densq*p3Dp5*p6Dp8
     . *gwsq**6*gsq**2*mt**2*mbsq/wmass**4*hdecay
C (gw/rt2)^4*(gw/2)^2*16 from amplitude

      p4Dq=p(4,4)*q(4)-p(4,1)*q(1)-p(4,2)*q(2)-p(4,3)*q(3)
      p7Da=p(7,4)*a(4)-p(7,1)*a(1)-p(7,2)*a(2)-p(7,3)*a(3)

      p1Dr=p(1,4)*r(4)-p(1,1)*r(1)-p(1,2)*r(2)-p(1,3)*r(3) 
      p2Dr=p(2,4)*r(4)-p(2,1)*r(1)-p(2,2)*r(2)-p(2,3)*r(3)

      p1Db=p(1,4)*b(4)-p(1,1)*b(1)-p(1,2)*b(2)-p(1,3)*b(3)
      p2Db=p(2,4)*b(4)-p(2,1)*b(1)-p(2,2)*b(2)-p(2,3)*b(3)

      p1Dx=p(1,4)*x(4)-p(1,1)*x(1)-p(1,2)*x(2)-p(1,3)*x(3)
      p2Dx=p(2,4)*x(4)-p(2,1)*x(1)-p(2,2)*x(2)-p(2,3)*x(3)

      p1Dy=p(1,4)*y(4)-p(1,1)*y(1)-p(1,2)*y(2)-p(1,3)*y(3)
      p2Dy=p(2,4)*y(4)-p(2,1)*y(1)-p(2,2)*y(2)-p(2,3)*y(3)

      gamq4=qDq/(2d0*p4Dq)
      gama7=aDa/(2d0*p7Da)

      gamb1=bDb/(2d0*p1Db)
      gamb2=bDb/(2d0*p2Db)
      gamr1=rDr/(2d0*p1Dr)
      gamr2=rDr/(2d0*p2Dr)
      gamx1=xDx/(2d0*p1Dx)
      gamx2=xDx/(2d0*p2Dx)
      gamy1=yDy/(2d0*p1Dy)
      gamy2=yDy/(2d0*p2Dy)

C     now the momenta 3,5,6,8,9,10 are no longer needed
C     so set       
      do nu=1,4
      p(q4,nu)=q(nu)-gamq4*p(4,nu)
      p(a7,nu)=a(nu)-gama7*p(7,nu)
      p(r1,nu)=r(nu)-gamr1*p(1,nu)
      p(r2,nu)=r(nu)-gamr2*p(2,nu)
      p(b1,nu)=b(nu)-gamb1*p(1,nu)
      p(b2,nu)=b(nu)-gamb2*p(2,nu)
      enddo      

      call qqbtth(p,denr,denb,wtqqb)
      wtgg=0d0
C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=aveqq*fac*wtqqb
      elseif (j .eq. 0) then
          msq(j,j)=avegg*fac*wtgg
      elseif (j .gt. 0) then
          msq(j,-j)=aveqq*fac*wtqqb
      endif
      enddo
      return
      end
