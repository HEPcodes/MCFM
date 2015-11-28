      subroutine qqb_z1jet_v(p,msq)
      implicit none
c----Matrix element for Z + jet production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+g(p5)
c---
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
      integer j,k,
     . iqqbgLL(5),iqqbgLR(5),iqqbgRL(5),iqqbgRR(5),
     . iqgqLL(5),iqgqLR(5),iqgqRL(5),iqgqRR(5),
     . igqqLL(5),igqqLR(5),igqqRL(5),igqqRR(5)
      double precision msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     . p(mxpart,4),fac,prop,sz,virt5
      double precision qqbZgLL,qqbZgRR,qqbZgLR,qqbZgRL
      double precision gqZqLL,gqZqRR,gqZqLR,gqZqRL
      double precision qgZqLL,qgZqRR,qgZqLR,qgZqRL
      double precision qbqZgLL,qbqZgRR,qbqZgLR,qbqZgRL
      double precision gqbZqbLL,gqbZqbRR,gqbZqbLR,gqbZqbRL
      double precision qbgZqbLL,qbgZqbRR,qbgZqbLR,qbgZqbRL,
     . xl12,xl15,xl25,subqqb,subqbq,subqg,subgq,subqbg,subgqb,
     . ii_qg,ii_gq,ii_gg,if_qg,if_gg,fi_qg,fi_gg

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

      subqqb=0.5d0*xn*(if_qg(one,xl15,1)+0.5d0*fi_gg(one,xl15,1)
     &                +if_qg(one,xl25,1)+0.5d0*fi_gg(one,xl25,1))
     &        -one/xn*(ii_qg(one,xl12,1))
     &            +tr*ii_gq(one,xl12,1)
      subqbq=subqqb

      subqg=0.5d0*xn*(ii_qg(one,xl12,1)+ii_gg(one,xl12,1)
     &               +if_gg(one,xl25,1)+fi_qg(one,xl25,1))
     &     -0.5d0/xn*(if_qg(one,xl15,1)+fi_qg(one,xl15,1))
     &           +tr*ii_gq(one,xl12,1)

      subqbg=subqg

      subgq=0.5d0*xn*(ii_qg(one,xl12,1)+ii_gg(one,xl12,1)
     &               +if_gg(one,xl15,1)+fi_qg(one,xl15,1))
     &     -0.5d0/xn*(if_qg(one,xl25,1)+fi_qg(one,xl25,1))
     &           +tr*ii_gq(one,xl12,1)

      subgqb=subgq

      subqqb=subqqb*ason2pi
      subqbq=subqbq*ason2pi
      subqg =subqg *ason2pi
      subqbg=subqbg*ason2pi
      subgq =subgq *ason2pi
      subgqb=subgqb*ason2pi


c--   calculate propagator
      sz=s(3,4)
      prop=sz/sqrt((sz-zmass**2)**2+(zmass*zwidth)**2)

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
          msq(j,k)=+(Q(j)*q1+L(j)*l1*prop)**2*qqbZgLL
     .             +(Q(j)*q1+R(j)*r1*prop)**2*qqbZgRR
     .             +(Q(j)*q1+L(j)*r1*prop)**2*qqbZgLR
     .             +(Q(j)*q1+R(j)*l1*prop)**2*qqbZgRL
     .             +subqqb*msq0(j,k)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(Q(k)*q1+L(k)*l1*prop)**2*qbqZgLL
     .             +(Q(k)*q1+R(k)*r1*prop)**2*qbqZgRR
     .             +(Q(k)*q1+L(k)*r1*prop)**2*qbqZgLR
     .             +(Q(k)*q1+R(k)*l1*prop)**2*qbqZgRL
     .             +subqbq*msq0(j,k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+(Q(j)*q1+L(j)*l1*prop)**2*qgZqLL
     .             +(Q(j)*q1+R(j)*r1*prop)**2*qgZqRR
     .             +(Q(j)*q1+L(j)*r1*prop)**2*qgZqLR
     .             +(Q(j)*q1+R(j)*l1*prop)**2*qgZqRL
     .             +subqg*msq0(j,k)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbLL
     .             +(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbRR
     .             +(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbLR
     .             +(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbRL
     .             +subqbg*msq0(j,k)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(Q(k)*q1+L(k)*l1*prop)**2*gqZqLL
     .             +(Q(k)*q1+R(k)*r1*prop)**2*gqZqRR
     .             +(Q(k)*q1+L(k)*r1*prop)**2*gqZqLR
     .             +(Q(k)*q1+R(k)*l1*prop)**2*gqZqRL
     .             +subgq*msq0(j,k)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=+(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbLL
     .             +(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbRR
     .             +(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbLR
     .             +(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbRL
     .             +subgqb*msq0(j,k)
      endif
      
   19 continue
      enddo
      enddo
      
      return
      end
