      subroutine storeqq(mqq_dip,iqq_dip,mqqx_dip)
      implicit none
      include 'constants.f'
      include 'isq.f'
      integer i,j,k,l,m,n
      double precision mqq(0:2,-nf:nf,-nf:nf),iqq_dip(-nf:nf,-nf:nf)
      double precision mqq_dip(0:2,-nf:nf,-nf:nf)
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision mqqx_dip(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      common/mqq/mqq
      common/msqx/msqx
      do j=-nf,nf
      do k=-nf,nf
        do i=0,2
          mqq_dip(i,j,k)=mqq(i,j,k)
c          do l=-nf,nf
c          do m=-nf,nf
c            mqqx_dip(i,j,k,l,m)=msqx(i,j,k,l,m)
c          enddo
c          enddo
        enddo
        iqq_dip(j,k)=isq(j,k)
      enddo
      enddo
      return
      end
