      subroutine qqb_z1jet_v(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John Campbell                            *
*     May, 2001.                                                       *
*     Matrix element for Z + jet production                            *
*     in order alpha_s^2                                               *
*     averaged over initial colours and spins                          *
*     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+g(p5)                        *
************************************************************************
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'prods.f'
      include 'cutoff.f'
      include 'epinv.f'
      include 'scale.f'
      include 'flags.f'
      integer j,k,
     . iqqbgLL(5),iqqbgLR(5),iqqbgRL(5),iqqbgRR(5),
     . iqgqLL(5),iqgqLR(5),iqgqRL(5),iqgqRR(5),
     . igqqLL(5),igqqLR(5),igqqRL(5),igqqRR(5)
      double precision msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     . p(mxpart,4),fac,sz,virt5
      double precision qqbZgLL,qqbZgRR,qqbZgLR,qqbZgRL
      double precision gqZqLL,gqZqRR,gqZqLR,gqZqRL
      double precision qgZqLL,qgZqRR,qgZqLR,qgZqRL
      double precision qbqZgLL,qbqZgRR,qbqZgLR,qbqZgRL
      double precision gqbZqbLL,gqbZqbRR,gqbZqbLR,gqbZqbRL
      double precision qbgZqbLL,qbgZqbRR,qbgZqbLR,qbgZqbRL,
     . xl12,xl15,xl25,subqqb,subqbq,subqg,subgq,subqbg,subgqb,subuv,
     . subqqb_QQb,subqbq_QBQ,subqg_qg,subgq_gq,subqbg_qbg,subgqb_gqb,
     . ii_qg,ii_gq,ii_gg,if_qg,if_gg,fi_qg,fi_gg,fi_gq
      double complex prop
      
      data iqqbgLL/1,2,3,4,5/,iqqbgRR/2,1,4,3,5/
      data iqqbgRL/2,1,3,4,5/,iqqbgLR/1,2,4,3,5/

      data iqgqLL/1,5,3,4,2/,iqgqRR/5,1,4,3,2/
      data iqgqRL/5,1,3,4,2/,iqgqLR/1,5,4,3,2/

      data igqqLL/2,5,3,4,1/,igqqRR/5,2,4,3,1/
      data igqqRL/5,2,3,4,1/,igqqLR/2,5,4,3,1/

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--calculate spinor and dot-products (using BDK type notation)

      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities

      if ((abs(s(1,5)).lt.cutoff) .or. (abs(s(2,5)).lt.cutoff)) return

c--- calculate lowest order
      call qqb_z_g(p,msq0)
      
c---add result of integrating subtraction terms
      xl12=dlog(+s(1,2)/musq)
      xl15=dlog(-s(1,5)/musq)
      xl25=dlog(-s(2,5)/musq)

      subqqb=0d0
      subqg=0d0
      subgq=0d0

      if (Gflag) then
      subqqb=subqqb+xn*(if_qg(one,xl15,1)+0.5d0*fi_gg(one,xl15,1)
     &                +if_qg(one,xl25,1)+0.5d0*fi_gg(one,xl25,1))
     &                -2d0/xn*(ii_qg(one,xl12,1))
     &                +2d0*tr*ii_gq(one,xl12,1)

      subqg=subqg+xn*(ii_qg(one,xl12,1)+ii_gg(one,xl12,1)
     &               +if_gg(one,xl25,1)+fi_qg(one,xl25,1))
     &     -1d0/xn*(if_qg(one,xl15,1)+fi_qg(one,xl15,1))
     &           +2d0*tr*ii_gq(one,xl12,1)


      subgq=subgq+xn*(ii_qg(one,xl12,1)+ii_gg(one,xl12,1)
     &               +if_gg(one,xl15,1)+fi_qg(one,xl15,1))
     &     -1d0/xn*(if_qg(one,xl25,1)+fi_qg(one,xl25,1))
     &           +2d0*tr*ii_gq(one,xl12,1)


      endif

      if (Qflag) then
c--- all the ii_gq terms are zero, hence commented out
c--- the fi_gq contribution is already included in the fi_gg
c--- piece and just provides the UV subtraction
c      subqqb=subqqb+2d0*tr*dfloat(nf)*fi_gq(one,xl25,1)
c      subqg=subqg+(xn-1d0/xn)*ii_gq(one,xl25,1)
c     .           +(xn-1d0/xn)*ii_gq(one,xl12,1)
c      subgq=subgq+(xn-1d0/xn)*ii_gq(one,xl15,1)
c     .           +(xn-1d0/xn)*ii_gq(one,xl12,1)
      endif

      subqqb=subqqb*ason2pi*half
      subqg=subqg*ason2pi*half
      subgq=subgq*ason2pi*half

c--- UV counter-term
      
      subuv=ason2pi*xn*(epinv-log(fourpi))*
     .              (11d0-2d0*dble(nf)/xn)/3d0
      subqqb=subqqb-0.5d0*subuv
      subqg =subqg-0.5d0*subuv
      subgq =subgq-0.5d0*subuv

      subqbq=subqqb
      subqbg=subqg
      subgqb=subgq

c--   calculate propagator
      sz=s(3,4)
      prop=sz/dcmplx((sz-zmass**2),zmass*zwidth)

c---??Needs to be modified?? - Think I fixed, XN -> XNSQ, JMC
      FAC=8d0*CF*XNSQ*esq**2*gsq

      qqbZgLL=aveqq*fac*virt5(iqqbgLL,za,zb)
      qqbZgLR=aveqq*fac*virt5(iqqbgLR,za,zb)
      qqbZgRL=aveqq*fac*virt5(iqqbgRL,za,zb)
      qqbZgRR=aveqq*fac*virt5(iqqbgRR,za,zb)

      qbqZgLL=qqbZgRL
      qbqZgLR=qqbZgRR
      qbqZgRL=qqbZgLL
      qbqZgRR=qqbZgLR

      gqZqLL=aveqg*fac*virt5(igqqLL,za,zb)
      gqZqLR=aveqg*fac*virt5(igqqLR,za,zb)
      gqZqRL=aveqg*fac*virt5(igqqRL,za,zb)
      gqZqRR=aveqg*fac*virt5(igqqRR,za,zb)

      gqbZqbRL=gqZqLL
      gqbZqbRR=gqZqLR
      gqbZqbLL=gqZqRL
      gqbZqbLR=gqZqRR

      qgZqLL=aveqg*fac*virt5(iqgqLL,za,zb)
      qgZqLR=aveqg*fac*virt5(iqgqLR,za,zb)
      qgZqRL=aveqg*fac*virt5(iqgqRL,za,zb)
      qgZqRR=aveqg*fac*virt5(iqgqRR,za,zb)

      qbgZqbRL=qgZqLL
      qbgZqbRR=qgZqLR
      qbgZqbLL=qgZqRL
      qbgZqbLR=qgZqRR
      


      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=0d0
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=+cdabs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgLL
     .             +cdabs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgRR
     .             +cdabs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgLR
     .             +cdabs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgRL
     .             +subqqb*msq0(j,k)
c          write(*,79) j,k,msq(j,k)-subqqb*msq0(j,k),
c     .                subqqb*msq0(j,k),msq(j,k)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=+cdabs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgLL
     .             +cdabs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgRR
     .             +cdabs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgLR
     .             +cdabs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgRL
     .             +subqbq*msq0(j,k)
c          write(*,79) j,k,msq(j,k)-subqbq*msq0(j,k),
c     .                subqbq*msq0(j,k),msq(j,k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+cdabs(Q(j)*q1+L(j)*l1*prop)**2*qgZqLL
     .             +cdabs(Q(j)*q1+R(j)*r1*prop)**2*qgZqRR
     .             +cdabs(Q(j)*q1+L(j)*r1*prop)**2*qgZqLR
     .             +cdabs(Q(j)*q1+R(j)*l1*prop)**2*qgZqRL
     .             +subqg*msq0(j,k)
c          write(*,79) j,k,msq(j,k)-subqg*msq0(j,k),
c     .                subqg*msq0(j,k),msq(j,k)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+cdabs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbLL
     .             +cdabs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbRR
     .             +cdabs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbLR
     .             +cdabs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbRL
     .             +subqbg*msq0(j,k)
c          write(*,79) j,k,msq(j,k)-subqbg*msq0(j,k),
c     .                subqbg*msq0(j,k),msq(j,k)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=+cdabs(Q(k)*q1+L(k)*l1*prop)**2*gqZqLL
     .             +cdabs(Q(k)*q1+R(k)*r1*prop)**2*gqZqRR
     .             +cdabs(Q(k)*q1+L(k)*r1*prop)**2*gqZqLR
     .             +cdabs(Q(k)*q1+R(k)*l1*prop)**2*gqZqRL
     .             +subgq*msq0(j,k)
c          write(*,79) j,k,msq(j,k)-subgq*msq0(j,k),
c     .                subgq*msq0(j,k),msq(j,k)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=+cdabs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbLL
     .             +cdabs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbRR
     .             +cdabs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbLR
     .             +cdabs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbRL
     .             +subgqb*msq0(j,k)
c          write(*,79) j,k,msq(j,k)-subgqb*msq0(j,k),
c     .                subgqb*msq0(j,k),msq(j,k)
      endif
      
   19 continue
      enddo
      enddo

c      pause
   79 format(2i3,3f14.9)   
   
      return
      end
