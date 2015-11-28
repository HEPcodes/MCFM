      subroutine qqb_QQbg_mad(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2001.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  T(p3)+TB(p4)+_g(p5)
C -------Madgraph version
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'hardscale.f'
      include 'flags.f'
      integer i,j,k,nu,hq,Qh,hg,lh,n1,n2,nquark,j1,k2,
     . jj(-nf:nf),kk(-nf:nf)
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)
      double precision sgg_ttbg,suub_ttbg,sug_ttbu,subg_ttbub
      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data kk/-1,-2,-1,-2,-1,0,1,2,1,2,1/

c--- implement the momentum exchange      
      do i=1,4
        if (i.lt.4) then
          j=i
        else
          j=0
        endif 
        p1(j)=-p(1,i)
        p2(j)=-p(2,i)
        p3(j)=p(3,i)
        p4(j)=p(4,i)
        p5(j)=p(5,i)
      enddo

      call initialize

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf

      if ((j .eq. 0) .and. (k .eq. 0)) then
c-gg
      msq(j,k)=sgg_ttbg(p1,p2,p3,p4,p5)
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
c-qqb
      if (k.eq.-j) msq(j,k)=suub_ttbg(p1,p2,p3,p4,p5)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
c-qbq
      if (k.eq.-j) msq(j,k)=suub_ttbg(p2,p1,p3,p4,p5)


      elseif ((j .gt. 0) .and. (k .eq. 0)) then
C-qg
              msq(j,k)=sug_ttbu(p1,p2,p3,p4,p5)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
c-qbg
              msq(j,k)=subg_ttbub(p1,p2,p3,p4,p5)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
c-gq
              msq(j,k)=sug_ttbu(p2,p1,p3,p4,p5)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
c-gqb
              msq(j,k)=subg_ttbub(p2,p1,p3,p4,p5)
      endif

   
      enddo
      enddo

      return
      end
