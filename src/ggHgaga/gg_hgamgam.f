      subroutine gg_hgamgam(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c    g(-p1)+g(-p2)-->H --> gamma(p3) + gamma(p4) 
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s12
      double precision decay,gg,Asq,msqgamgam
c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      s12=2d0*(p(1,4)*p(2,4)-p(1,1)*p(2,1)-p(1,2)*p(2,2)-p(1,3)*p(2,3))

      decay=msqgamgam(hmass)/((s12-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(3d0*pi))**2/vevsq
      gg=0.5d0*V*Asq*s12**2

c---calculate propagators
      msq(0,0)=avegg*gg*decay

      return
      end

      double precision function msqgamgam(mh)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      double complex Iw,Iq,Ftriangle,ctwo
      double precision prefac,mh
      double precision x_t,x_b,x_w,x
      parameter(ctwo=(2d0,0d0))
C---statement functions
      Iq(x)=dcmplx(4d0*x)*(ctwo+dcmplx(4d0*x-1d0)*Ftriangle(x))
      Iw(x)=-ctwo*(dcmplx(6d0*x+1d0)
     & +dcmplx(6d0*x*(2d0*x-1d0))*Ftriangle(x))
C---end statement functions

C---Total width multiplied by a factor of 16*pi*mh 
C---to get matrix element squared.
c---maybe it would be better to add esq at a higher scale.
      prefac=(esq/(4d0*pi))**2*Gf*mh**4/(8d0*rt2*pi**2)
      x_b=mbsq/mh**2
      x_t=(mt/mh)**2
      x_w=(wmass/mh)**2
      msqgamgam=prefac
     & *abs(xn*(Q(1)**2*Iq(x_b)+Q(2)**2*Iq(x_t))+Iw(x_w))**2
      return
      end

      double complex function Ftriangle(x)
      implicit none 
      include 'constants.f'
      double precision x,y
      double complex arg
      y=1d0-4d0*x
      if (y .gt. 0d0) then
      arg=dcmplx((1d0+sqrt(y))/(1d0-sqrt(y)))
      Ftriangle=+dcmplx(0.5d0)*(log(arg)-impi)**2
      elseif (y .le. 0d0) then
      Ftriangle=-dcmplx(2d0)*dcmplx(asin(half/sqrt(x)))**2
      endif
      return 
      end 


