      subroutine qqb_w_twdk(p,msq)
      implicit none
c----Matrix element for W + t production
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + t~(e^-(p5)+m(p6)+b~(p7))
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ t(m(p5)+e^+(p6)+b(p7))
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision qgWq,qbgWqb,gqbWqb,gqWq
      double precision ubtdg_h

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
C--setup z products through common block
      call spinoru(7,p,za,zb)

      fac=gwsq**4*gsq*V
      
C ubtdg_h(ju,jb,jn,je,jc,jd,jg,p)

      if     (nwz .eq. -1) then
      qgWq=aveqg*fac*ubtdg_h(4,1,5,6,7,3,2,p)
      gqWq=aveqg*fac*ubtdg_h(4,2,5,6,7,3,1,p)
      elseif (nwz .eq. +1) then
      gqbWqb=aveqg*fac*ubtdg_h(3,2,6,5,7,4,1,p)
      qbgWqb=aveqg*fac*ubtdg_h(3,1,6,5,7,4,2,p)
      else
        write(6,*) 'Problem with nwz in qqb_w_tdk.f: nwz=',nwz
        stop
      endif


      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j .eq. +5) .and. (k .eq. 0) .and. (nwz .eq. -1)) then
          msq(j,k)=qgWq
      elseif ((j .eq. -5) .and. (k .eq. 0) .and. (nwz .eq. +1)) then
          msq(j,k)=qbgWqb
      elseif ((j .eq. 0) .and. (k .eq. +5) .and. (nwz .eq. -1)) then
          msq(j,k)=gqWq
      elseif ((j .eq. 0) .and. (k .eq. -5) .and. (nwz .eq. +1)) then
          msq(j,k)=gqbWqb
      endif

      enddo
      enddo
      return
      end
 
