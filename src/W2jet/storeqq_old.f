      subroutine storeqq(mqq_dip,iqq_dip)
      implicit none
      include 'constants.f'
      include 'isq.f'
      integer i,j,k,n
      double precision mqq(0:2,-nf:nf,-nf:nf),iqq_dip(-nf:nf,-nf:nf)
      double precision mqq_dip(0:2,-nf:nf,-nf:nf)
      common/mqq/mqq
      do j=-nf,nf
      do k=-nf,nf
        do i=0,2
          mqq_dip(i,j,k)=mqq(i,j,k)
        enddo
          iqq_dip(j,k)=isq(j,k)
      enddo
      enddo
      return
      end
