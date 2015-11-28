      subroutine boosta(p1,pin,pout)
c---Boost pin to centre of mass frame of p1 and call resultant 
c---four-vector pout
      implicit none
      include 'constants.f'
      integer k
      double precision p1(4),pin(4),pout(4),beta(3),gam,bdotp,m1
      m1=dsqrt(p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2)
      gam=p1(4)/m1
      bdotp=0d0
      do k=1,3
      beta(k)=p1(k)/p1(4)
      bdotp=bdotp+beta(k)*pin(k)
      enddo
      pout(4)=gam*(pin(4)-bdotp)

      do k=1,3
      pout(k)=pin(k)+gam*beta(k)*(gam/(gam+one)*bdotp-pin(4))
      enddo

      return
      end
 
