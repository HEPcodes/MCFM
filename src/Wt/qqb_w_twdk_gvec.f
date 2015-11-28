      subroutine qqb_w_twdk_gvec(p,n,in,msq)
      implicit none
************************************************************************
*     Authors: John Campbell & Francesco Tramontano                    *
*     February, 2005.                                                  *
*    Matrix element squared and averaged over initial colours and spins*
*     with gluon index contracted with vector n                        *
*         ip emitter                                                   *
*         kp spectator                                                 *
*         in label of gluon which is contracted with n                 *
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),n(4),
     . p1p2(-1:1,-1:1),wtgvecn
     
      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

C---- Need to check the overall factors here
      if (in .eq. 1) then
c      p1p2(0,-1)=-aveqg*fac*wtgvecn(1,2,3,4,6,5,7,p,n)
      p1p2(0,+1)=aveqg*wtgvecn(1,2,3,4,5,6,7,p,n)
c      write(6,*) 'pg ',p(1,4),p(1,1),p(1,2),p(1,3)
c      write(6,*) 'vec',n(4),n(1),n(2),n(3)
c      write(6,*) 'wtgvecn(1,2,3,4,6,5,7,p,n)',wtgvecn(1,2,3,4,6,5,7,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=aveqg*wtgvecn(2,1,3,4,5,6,7,p,n)
c      p1p2(-1,0)=-aveqg*fac*wtgvecn(5,1,2,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=p1p2(+1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=p1p2(0,-1)
      endif

      enddo
      enddo
      
      return
      end
 
      double precision function wtgvecn(ig,is,ie,in,jn,je,jb,p,vec)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'zprods_com.f'
      integer ig,is,ie,in,je,jn,jb
      double precision p(mxpart,4),vec(4),prop,nDpg,fac
      double complex amp
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba

      nDpg=vec(4)*p(ig,4)-vec(1)*p(ig,1)-vec(2)*p(ig,2)-vec(3)*p(ig,3)

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDpg).gt.1d-3*abs(p(ig,4))) then
         write(*,*) 'Error for :',ig,is,ie,in,je,jn,jb
         write(*,*) 'cutoff',1d-3*abs(p(ig,4))
         write(6,*) 'nDpg',nDpg
         call flush(6)
         stop
      endif

      call spinoru(7,p,za,zb)
      call spinork(7,p,zab,zba,vec)
      call ampsn(p,ig,is,ie,in,jn,je,jb,amp)

      prop=(dble(za(ie,in)*zb(in,ie))-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((dble(za(jn,je)*zb(je,jn))-wmass**2)**2
     . +(wmass*wwidth)**2)
      prop=prop*(
     .(dble(za(jn,je)*zb(je,jn)+za(jn,jb)*zb(jb,jn)+za(je,jb)*zb(jb,je))
     . -mt**2)**2+(mt*twidth)**2)

      fac=xn*cf*gsq*gw**8
      
      wtgvecn=fac*abs(amp)**2/prop

      return
      end


      subroutine ampsn(p,ig,is,ie,in,jn,je,jb,amp)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex amp
      double precision p(mxpart,4)
      double precision dot,taugt,taugs
      integer ig,is,ie,in,je,jn,jb
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba
      taugs=two*dot(p,ig,is)
      taugt=two*(dot(p,ig,jn)+dot(p,ig,je)+dot(p,ig,jb))
      amp= za(ig,ie)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(jn,ig)*
     & taugt**(-1) + za(ig,ie)*za(jn,jb)*zb(is,in)*zb(je,jb)*zab(jb,ig)
     & *taugt**(-1) - za(ie,je)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(jn,je
     & )*taugt**(-1) - za(ie,je)*za(jn,jb)*zb(is,in)*zb(je,jb)*zab(jb,
     & je)*taugt**(-1) + za(ie,jn)*za(jn,jb)*zb(ig,in)*zb(je,jn)*zab(ig
     & ,is)*taugs**(-1) + za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(
     & is,is)*taugs**(-1) - za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jn)*
     & zab(jn,jn)*taugt**(-1) - za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jb)
     & *zab(jb,jn)*taugt**(-1) + za(ie,jb)*za(jn,jb)*zb(ig,in)*zb(je,jb
     & )*zab(ig,is)*taugs**(-1) - za(ie,jb)*za(jn,jb)*zb(is,in)*zb(je,
     & jn)*zab(jn,jb)*taugt**(-1) + za(ie,jb)*za(jn,jb)*zb(is,in)*zb(je
     & ,jb)*zab(is,is)*taugs**(-1) - za(ie,jb)*za(jn,jb)*zb(is,in)*zb(
     & je,jb)*zab(jb,jb)*taugt**(-1) + za(jn,jb)*zb(is,in)*zab(ie,je)*
     & mt**2*taugt**(-1)

      return
      end
