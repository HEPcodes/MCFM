      subroutine qqb_hbbbar_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c----for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H --> (b(p3)+bbar(p4))+g(p5)
c---
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     . msq15_2v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     . sub15_2(4),sub25_1(4),sub15_2v,sub25_1v
      external qqb_hbbbar,qqb_hbbbar_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     . qqb_hbbbar,qqb_hbbbar_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     . qqb_hbbbar,qqb_hbbbar_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo


      if     ((j .ne. 0) .and. (k .eq. 0)) then
         msq(1,j,k)=two*cf
     .   *(msq15_2(0,0)*sub15_2(gq)+msq15_2v(0,0)*sub15_2v)
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(2,j,k)=two*cf
     .   *(msq25_1(0,0)*sub25_1(gq)+msq25_1v(0,0)*sub25_1v)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
         msq(1,j,k)=two*xn
     .   *(msq15_2(j,k)*sub15_2(gg)+msq15_2v(j,k)*sub15_2v)
         msq(2,j,k)=two*xn
     .   *(msq25_1(j,k)*sub25_1(gg)+msq25_1v(j,k)*sub25_1v)
      endif


      enddo
      enddo

      return      
      end


