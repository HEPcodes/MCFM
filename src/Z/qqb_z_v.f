      subroutine qqb_z_v(p,msqv)
      implicit none
      integer j,k
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      double precision p(mxpart,4),
     . msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     . xlog,virt
      common/xlog/xlog
c---  calculate lowest order matrix element
      call qqb_z(p,msq)
c---calculate the multiple of the lowest order
      virt=ason2pi*cf*(two*pisq/3d0-eight+3d0*xlog)
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=virt*msq(j,k)
      enddo
      enddo
      end

