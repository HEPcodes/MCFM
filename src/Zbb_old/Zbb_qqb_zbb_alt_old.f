      subroutine qqb_zbb_alt(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> b(p4)+bb(p5)+e^-(p6)+e^+(p7)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k,polq,polb,polz
      double precision p(mxpart,4),q(6,4),msq(-nf:nf,-nf:nf),faclo
      double precision v2(2),vQ(nf,2)
      double complex atreez,amp

      faclo=4d0*V*gsq**2*esq**2*aveqq

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c---transform momenta into BDKW notation
      do k=1,4
      q(1,k)=p(2,k)
      q(2,k)=p(5,k)
      q(3,k)=p(4,k)
      q(4,k)=p(1,k)
      q(5,k)=p(7,k)
      q(6,k)=p(6,k)
      enddo      

      call spinoru(6,q,za,zb)

c ensure that we have a hard process
      if (
     .      (s(2,3) .lt. four*mbsq) 
     . .or. (s(3,4)*s(1,3)/s(1,4) .lt. mbsq) 
     . .or. (s(2,4)*s(1,2)/s(1,4) .lt. mbsq) ) return

      v2(1)=l1
      v2(2)=r1

      do j=1,nf
        vQ(j,1)=L(j)
        vQ(j,2)=R(j)
      enddo

      do j=-nf,nf
      k=-j

      do polq=1,2
      do polb=1,2
      do polz=1,2
        if     ((j .eq. 0) .and. (k .eq. 0)) then
          amp=0d0
        elseif ((j .gt. 0) .and. (k .lt. 0)) then
          amp=atreez(polq,polb,polz,1,2,3,4,5,6,za,zb)
     .        *vQ(j,polq)*v2(polz)
        elseif ((j .lt. 0) .and. (k. gt. 0)) then
          amp=atreez(3-polq,polb,polz,1,2,3,4,5,6,za,zb)
     .        *vQ(k,polq)*v2(polz)
        endif
        msq(j,k)=msq(j,k)+cdabs(amp)**2
      enddo
      enddo
      enddo
      
      msq(j,k)=msq(j,k)*faclo
      
      enddo

      return
      end







