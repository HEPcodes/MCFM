      subroutine qqb_QQbdk_mad(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2002.                                                      *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'flags.f'
      integer i,j,k,nu,hq,Qh,hg,lh,n1,n2
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision mmsq_gg,mmsq_qqb,mmsq_qbq,mmsq_qg
      double precision mmsq_gq,mmsq_gqb,mmsq_qbg
      double precision fac,pswap(mxpart,4),LRq(2),LRb(2),lr1(2),Vfac
      real*8 P1(0:3,8),ans
      double complex prop,czq,czb

c--- implement the momentum exchange      

      do i=1,4
        if (i.lt.4) then
          j=i
        else
          j=0
        endif 
        p1(j,1)=-p(1,i)
        p1(j,2)=-p(2,i)
        p1(j,3)=p(3,i)
        p1(j,4)=p(4,i)
        p1(j,5)=p(5,i)
        p1(j,6)=p(6,i)
        p1(j,7)=p(7,i)
        p1(j,8)=p(8,i)
      enddo

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      call sgg_lo(p1,ans)
      do j=-nf,nf
      do k=-nf,nf
      
      if     ((j .eq. 0) .and. (k .eq. 0)) then
        msq(j,k)=ans
      endif
      
      enddo
      enddo

      return
      end
