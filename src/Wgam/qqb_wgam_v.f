      subroutine qqb_wgam_v(p,msqv)
      implicit none
C----Author R.K.Ellis August 2002
C====Virtual corrections to
c     q(-p1)+qbar(-p2)-->e^-(p3)+nu(p4)+gamma(p5)
      include 'constants.f'
      include 'scheme.f'
      include 'zerowidth.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'masses.f'
      include 'dprodx.f'
      include 'sprodx.f'
      include 'epinv.f'
      include 'anomcoup.f'
      include 'nwz.f'
      integer j,k
      double precision msqv(-nf:nf,-nf:nf),p(mxpart,4),qbq,qqb,
     . msq(-nf:nf,-nf:nf)
      double precision fac
      double complex agamtree,agamvirt

      call spinoru(5,p,za,zb)

      fac=ason2pi*2d0*cf*aveqq*2d0*xn*gwsq**2*esq
      scheme='dred'
      if (nwz .eq. -1) then
C ie ub-d
      qbq=+fac*dble(
     . Conjg(agamtree(1,2,3,4,5,za,zb,+1))*agamvirt(1,2,3,4,5,za,zb,+1))
     .    +fac*dble(
     . Conjg(agamtree(1,2,3,4,5,za,zb,-1))*agamvirt(1,2,3,4,5,za,zb,-1))
C ie d-ub
      qqb=+fac*dble(
     . Conjg(agamtree(2,1,3,4,5,za,zb,+1))*agamvirt(2,1,3,4,5,za,zb,+1))
     .    +fac*dble(
     . Conjg(agamtree(2,1,3,4,5,za,zb,-1))*agamvirt(2,1,3,4,5,za,zb,-1))

      elseif (nwz .eq. +1) then 
C ie db-u
      qbq=+fac*dble(
     . Conjg(agamtree(2,1,4,3,5,zb,za,+1))*agamvirt(2,1,4,3,5,zb,za,+1))
     .    +fac*dble(
     . Conjg(agamtree(2,1,4,3,5,zb,za,-1))*agamvirt(2,1,4,3,5,zb,za,-1))
C ie u-db
      qqb=+fac*dble(
     . Conjg(agamtree(1,2,4,3,5,zb,za,+1))*agamvirt(1,2,4,3,5,zb,za,+1))
     .    +fac*dble(
     . Conjg(agamtree(1,2,4,3,5,zb,za,-1))*agamvirt(1,2,4,3,5,zb,za,-1))
      endif
      do j=-nf,nf
      do k=-nf,nf
c--- set msqv=0 to initalize
      msqv(j,k)=0d0
          if ((j .gt. 0) .and. (k .lt. 0)) then
            msqv(j,k)=Vsq(j,k)*qqb
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
            msqv(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo

      return
      end


      

      double complex function agamvirt(p1,p2,p3,p4,p5,za,zb,hgamma)
      implicit none
C----Author R.K.Ellis August 2002
c     q(-p1)+qbar(-p2)-->e^-(p3)+nu(p4)+gamma(p5)
      include 'constants.f'
      include 'nwz.f'
      include 'dprodx.f'
      include 'sprodx.f'
      double complex agamtree,vpole,fagamma,fbgamma
      integer p1,p2,p3,p4,p5,hgamma

          if (hgamma .eq. +1) then
          agamvirt=vpole(s(p1,p2))*agamtree(p1,p2,p3,p4,p5,za,zb,+1)
     .    +Qd*fagamma(p1,p2,p3,p4,p5,za,zb)
     .    +Qu*fbgamma(p1,p2,p3,p4,p5,za,zb)
          elseif (hgamma .eq. -1) then
          agamvirt=vpole(s(p2,p1))*agamtree(p1,p2,p3,p4,p5,za,zb,-1)
     .    +Qu*fagamma(p2,p1,p4,p3,p5,zb,za)
     .    +Qd*fbgamma(p2,p1,p4,p3,p5,zb,za)
          endif
      return
      end
