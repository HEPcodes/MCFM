      subroutine qqb_w1jet_v(p,msq)
      implicit none
c----Matrix element for W + jet production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->W^+(nu(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'prods.f'
      include 'cutoff.f'
      include 'epinv.f'
      include 'scale.f'
      integer j,k,iqqbg(5),iqbqg(5),iqgq(5),igqq(5),
     . igqbqb(5),iqbgqb(5)
      double precision msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     . p(mxpart,4),fac,sw,prop,
     . virt5,xl12,xl15,xl25,
     . subqqb,subqbq,subqg,subgq,subqbg,subgqb,
     . qqbWg,qbqWg,qgWq,gqWq,qbgWqb,gqbWqb,
     . ii_qg,ii_gq,ii_gg,if_qg,if_gg,fi_qg,fi_gg

      data iqqbg/1,2,3,4,5/
      data iqgq/1,5,3,4,2/
      data igqq/2,5,3,4,1/

      data iqbqg/2,1,3,4,5/
      data iqbgqb/5,1,3,4,2/
      data igqbqb/5,2,3,4,1/

c--calculate spinor and dot-products (using BDK type notation)

      call spinoru(5,p,za,zb)

c---protect from soft and collinear singularities
      if ((abs(s(1,5)) .lt. cutoff).or.(abs(s(2,5)) .lt. cutoff)) return

c--- calculate lowest order
      call qqb_w_g(p,msq0)
      
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

      sw=s(3,4)
      prop=sw**2/((sw-wmass**2)**2+(wmass*wwidth)**2)
      FAC=2d0*CF*xnsq*gwsq**2*gsq*prop

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo


      qqbWg=aveqq*fac*virt5(iqqbg,za,zb)
      qbqWg=aveqq*fac*virt5(iqbqg,za,zb)
      gqWq=aveqg*fac*virt5(igqq,za,zb)
      qgWq=aveqg*fac*virt5(iqgq,za,zb)
      gqbWqb=aveqg*fac*virt5(igqbqb,za,zb)
      qbgWqb=aveqg*fac*virt5(iqbgqb,za,zb)
      
      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=Vsq(j,k)*qqbWg+subqqb*msq0(j,k)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=Vsq(j,k)*qbqWg+subqbq*msq0(j,k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
     &     +subqg*msq0(j,k)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
     &     +subqbg*msq0(j,k)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
     &     +subgq*msq0(j,k)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
     &     +subgqb*msq0(j,k)
      endif

      enddo
      enddo
      
      return
      end
