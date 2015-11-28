      subroutine qqb_w_cjet_gvec(p,n,in,msq)
      implicit none
c----Matrix element for W production
C----averaged over initial colours and spins
c    contracted with the vector v(mu)
C For nwz=+1
c     c(-p1)+sbar(-p2)--> g(p5)+ W^+(n(p3)+e^+(p4))
C For nwz=-1
c     s(-p1)+cbar(-p2)--> g(p5)+ W^-(e^-(p3)+nbar(p4))
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'ckm.f'
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision wcjetn,p1p2(-1:1,-1:1),n(4)
     
      double precision FAC

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      fac=2d0*gsq*V*gwsq**2
      call dotem(5,p,s)

C----needs more work (not correct yet)
      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*wcjetn(5,2,1,p,n)
      p1p2(0,+1)=-aveqg*fac*wcjetn(2,5,1,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*wcjetn(1,5,2,p,n)
      p1p2(-1,0)=-aveqg*fac*wcjetn(5,1,2,p,n)
c      elseif (in .eq. 5) then     
c      p1p2(1,-1)=+aveqq*fac*wcjetn(1,2,5,p,n)
c      p1p2(-1,1)=+aveqq*fac*wcjetn(2,1,5,p,n)
      endif

      do j=-(nf-2),(nf-2)
      do k=-(nf-2),(nf-2)
c      if     ((j .gt. 0) .and. (k .lt. 0)) then
c          msq(j,k)=Vsq(j,k)*p1p2(1,-1)
c      elseif ((j .lt. 0) .and. (k .gt. 0)) then
c          msq(j,k)=Vsq(j,k)*p1p2(-1,1)
      if     ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    Vsq(j,-4)*p1p2(+1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    Vsq(j,+4)*p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    Vsq(-4,k)*p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    Vsq(+4,k)*p1p2(0,-1)
      endif

      enddo
      enddo
      
      return
      end
 
      double precision function wcjetn(p1,p2,p5,p,n)
      implicit none 
C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> W(l(p3)+a(p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5
      double precision p(mxpart,4),n(4),nDn,prop,
     .                 nDp1,nDp2,nDp3,nDp4,nDp5
      p3=3
      p4=4
      nDp1=n(4)*p(p1,4)-n(3)*p(p1,3)-n(2)*p(p1,2)-n(1)*p(p1,1)
      nDp2=n(4)*p(p2,4)-n(3)*p(p2,3)-n(2)*p(p2,2)-n(1)*p(p2,1)
      nDp3=n(4)*p(p3,4)-n(3)*p(p3,3)-n(2)*p(p3,2)-n(1)*p(p3,1)
      nDp4=n(4)*p(p4,4)-n(3)*p(p4,3)-n(2)*p(p4,2)-n(1)*p(p4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      nDp5=n(4)*p(p5,4)-n(3)*p(p5,3)-n(2)*p(p5,2)-n(1)*p(p5,1)

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp5).gt.1d-3*abs(p(p1,4))) then 
         write(*,*) 'Error for :',p1,p2,p3,p4,p5
         write(*,*) 'cutoff',1d-3*abs(p(p1,4))
         write(6,*) 'nDp5',nDp5
         call flush(6)
         stop
      endif

c---calculate the propagator
      prop=((s(p3,p4)-wmass**2)**2+(wmass*wwidth)**2)

      wcjetn=-nDp3*nDp4+0.25d0*s(p3,p4)*nDn
     . +s(p1,p4)*nDp2**2*(s(p3,p2)+s(p3,p5))/s(p2,p5)**2
     . +s(p3,p2)*nDp1**2*(s(p1,p4)+s(p4,p5))/s(p1,p5)**2
 
     . -(nDp2*(nDp3*s(p1,p4)-(s(p3,p2)+s(p3,p5))*nDp4)
     . +s(p3,p5)*0.25d0*nDn*(s(p4,p2)+s(p1,p4)))/s(p2,p5)
 
     . -(nDp1*(s(p3,p2)*nDp4-(s(p1,p4)+s(p4,p5))*nDp3)
     . +s(p4,p5)*0.25d0*nDn*(s(p3,p2)+s(p3,p1)))/s(p1,p5)

     . -(nDp1*nDp2*(s(p3,p2)*s(p4,p5)+s(p3,p5)*s(p1,p4)
     . +s(p3,p5)*s(p4,p5)+2d0*s(p3,p2)*s(p1,p4))
     . -s(p3,p5)*s(p4,p5)*s(p1,p2)*0.25d0*nDn)/s(p1,p5)/s(p2,p5)

c--- massless expression
c      wcjetn=(nDp1*s(p2,p3)/s(p1,p5)-nDp2*s(p1,p4)/s(p2,p5))**2
c     . +(s(p2,p3)*nDp1/s(p1,p5)-s(p1,p4)*nDp2/s(p2,p5))
c     . *(nDp2+nDp3-nDp4-nDp1)
c     . -(s(p1,p4)-s(p2,p3))**2*s(p3,p4)*nDn/4d0/s(p1,p5)/s(p2,p5)
c     . -(nDp1+nDp4)*(nDp2+nDp3)

      wcjetn=wcjetn/prop

      return
      end

