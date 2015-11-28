      subroutine qqb_w_cjet_massless_gs(p,msq)
************************************************************************
*     Author: J.M. Campbell                                            *
*     October, 2003.                                                   *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W + c [or c~](p5) + g(p6)
c                           |
c                            -->l(p3)+a(p4)

      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'ckm.f'
      integer j,k,nd
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),stat
      double precision 
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq26_5(-nf:nf,-nf:nf),msq16_5(-nf:nf,-nf:nf),
     & msq56_1v(-nf:nf,-nf:nf),msq56_2v(-nf:nf,-nf:nf),
     & msq26_5v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & msq16_2v(-nf:nf,-nf:nf),msq16_5v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),sub16_5(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub26_1v,sub16_5v,sub16_2v
      external qqb_w_cjet_massless,qqb_w_cjet_massless_gvec

c--- Note that the subtractions here are very similar to the ones
c--- for W+1 jet, except some singularities are not present

      ndmax=6

c--- calculate the initial-initial dipoles
      call dips(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     . qqb_w_cjet_massless,qqb_w_cjet_massless_gvec)
      call dips(4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     . qqb_w_cjet_massless,qqb_w_cjet_massless_gvec)

c--- now the initial-final ones
      call dips(5,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     . qqb_w_cjet_massless,qqb_w_cjet_massless_gvec)
      call dips(6,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     . qqb_w_cjet_massless,qqb_w_cjet_massless_gvec)

c--- now the final-initial ones
      call dips(5,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     . qqb_w_cjet_massless,qqb_w_cjet_massless_gvec)
      call dips(6,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     . qqb_w_cjet_massless,qqb_w_cjet_massless_gvec)

      do j=-nf,nf
      do k=-nf,nf      
      do nd=1,6
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo

      do j=-(nf-2),(nf-2)
      do k=-(nf-2),(nf-2)

      if     ((k .eq. 0).and. (j .ne. 0)) then
c--- q-g and qb-g cases
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-msq16_5(j,k)*(sub16_5(qq)+sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     .              +msq26_5(j,k)*sub56_2(qq))

      elseif ((j .eq. 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     .              +msq16_5(j,k)*sub56_1(qq))
      msq(6,j,k)=-msq26_5(j,k)*(sub26_5(qq)+sub56_2(qq))/xn

      elseif ((j .eq. 0).and.(k .eq. 0)) then
c--- g-g case
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     .           +msq16_2(-2,k)+msq16_2(-1,k)
     .           +msq16_2(+5,k)+msq16_2(+4,k)+msq16_2(+3,k)
     .           +msq16_2(+2,k)+msq16_2(+1,k))*sub16_2(qg)*2d0*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
     .           +msq26_1(k,-2)+msq26_1(k,-1)
     .           +msq26_1(k,+5)+msq26_1(k,+4)+msq26_1(k,+3)
     .           +msq26_1(k,+2)+msq26_1(k,+1))*sub26_1(qg)*2d0*tr

      elseif ((j .ne. 0) .and. (k .ne. 0)) then
c--- subtraction terms for 4-quark matrix elements
      msq(3,j,k)=2d0*cf*(
     .            msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=2d0*cf*(
     .            msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)

      endif     
      
      enddo
      enddo

      return
      end
