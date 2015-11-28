      subroutine qqb_wbb(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> nu(p3)+e^+(p4)+b(p5)+bb(p6)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ckm.f'
      include 'prods.f'
      include 'hardscale.f'
      integer j,k
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),msqwbb
      double precision qqb,qbq

C---Initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---Fill spinor products
      call spinoru(6,p,za,zb)

c ensure that we have a hard process
      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return

C--calculate matrix element squared
      qqb=msqwbb(1,2,5,6)
      qbq=msqwbb(2,1,5,6)

      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .lt. 0)) then
               msq(j,k)=Vsq(j,k)*qqb
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
               msq(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo
      return
      end







