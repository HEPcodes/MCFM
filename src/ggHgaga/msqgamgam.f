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

