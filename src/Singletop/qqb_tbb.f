      subroutine qqb_tbb(p,msq)
      implicit none
c     Matrix element for t-bbar production
c      u(-p1)+dbar(-p2)-->n(p3)+e^+(p4)+b(p5)+bbar(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision s1t,s2t,s4t,s345,prop
     
      integer j,k
      double precision fac,qqb,qbq,d

      d(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
 

c + 0.5d0*spinave*gw^8*s16*s53*[s12-mw^2]^-2*[ttDtt-mt^2]^-2*[s34-mw^2]^-2
c  * (s24*ttDtt-s4t*s2t )

      s1t=d(1,3)+d(1,4)+d(1,5)
      s2t=d(2,3)+d(2,4)+d(2,5)
      s4t=d(3,4)+d(4,4)+d(4,5)
      s345=0.5d0*d(5,5)+d(3,5)+d(4,5)+d(3,4)
      prop=     ((d(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s345-mt**2)**2+(mt*twidth)**2)
      prop=prop*((d(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      fac=spinave*gw**8/prop*d(5,3)
      qqb=fac*d(1,6)*(s2t*s4t-d(2,4)*s345)
      qbq=fac*d(2,6)*(s1t*s4t-d(1,4)*s345)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      if     ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsq(j,k)*qqb 
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsq(j,k)*qbq 
      endif
      enddo
      enddo

      return
      end
