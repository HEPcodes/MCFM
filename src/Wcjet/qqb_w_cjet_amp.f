      subroutine qqb_w_cjet_amp(p,msq)
      implicit none
c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + cbar(p5)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5) 
c---
      include 'constants.f'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision qgWq,qbgWqb,gqbWqb,gqWq,w1cjetamp
     

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      fac=gwsq**2*gsq*V
      gqbWqb=aveqg*fac*w1cjetamp(2,5,4,3,1,p)
      qgWq=  aveqg*fac*w1cjetamp(1,5,3,4,2,p)
      qbgWqb=aveqg*fac*w1cjetamp(1,5,4,3,2,p)
      gqWq=  aveqg*fac*w1cjetamp(2,5,3,4,1,p)

      do j=-nf,nf
      do k=-nf,nf
      if ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=Vsq(j,-4)*qgWq
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=Vsq(j,+4)*qbgWqb
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=Vsq(-4,k)*gqWq
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=Vsq(+4,k)*gqbWqb
      endif

      enddo
      enddo
      return
      end
 
      double precision function w1cjetamp(js,jc,je,jn,jg,p)
C     Matrix element squared for s(1) cbar(2) -> e-(3) nu(4) g(5)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer js,jc,jg,je,jn,nu,j
      double precision masssq,mass,scg,p(mxpart,4),q(mxpart,4),prop
      double complex amp(2,2)

      masssq=p(jc,4)**2-p(jc,1)**2-p(jc,2)**2-p(jc,3)**2
      scg=2d0*(p(jc,4)*p(jg,4)
     .        -p(jc,1)*p(jg,1)
     .        -p(jc,2)*p(jg,2)
     .        -p(jc,3)*p(jg,3))
      mass=sqrt(masssq)

C     setup massless vector basis
      do j=1,5
        do nu=1,4
          if (j.eq.jc) then
             q(j,nu)=p(jc,nu)-masssq/scg*p(jg,nu)
          else
             q(j,nu)=p(j,nu)
          endif
       enddo
      enddo

      call spinoru(5,q,za,zb)
      prop=((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      
C     First argument is c helicity, second is gluon helicity,
C     (s helicity is always negative)

      amp(1,1)=zb(jn,js)/(zb(jg,jc)*zb(jg,js))
     . *(zb(js,jc)*za(jc,je)+zb(js,jg)*za(jg,je))

      amp(2,1)=mass*zb(jn,js)/(zb(jg,jc)**2*zb(jg,js))
     . *za(jg,je)*zb(jc,js)

      amp(1,2)=za(jc,je)/(za(jg,jc)*za(jg,js))
     . *(zb(jn,js)*za(js,jc)+zb(jn,jg)*za(jg,jc))

      amp(2,2)=mass*za(jg,je)/(za(jg,jc)**2*za(jg,js))
     . *(zb(jn,js)*za(js,jc)+zb(jn,jg)*za(jg,jc))

      w1cjetamp=+cdabs(amp(1,1))**2+cdabs(amp(2,1))**2
     .          +cdabs(amp(1,2))**2+cdabs(amp(2,2))**2
      w1cjetamp=w1cjetamp/prop

      return
      end

