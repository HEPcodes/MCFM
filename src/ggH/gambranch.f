      double precision function gamwidth()
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcharge.f'
      integer j
      double precision arg,fac,x(6)
      double complex temp,Isubq,Isubw
C     calculates the width of the SM Higgs into gamma-gamma      
      fac=aem**2*Gf*hmass**3/(128d0*sqrt(2d0)*pi**3)
      x(1)=(md/hmass)**2
      x(2)=(mu/hmass)**2
      x(3)=(ms/hmass)**2
      x(4)=(mc/hmass)**2
      x(5)=(mb/hmass)**2
      x(6)=(mt/hmass)**2
      arg=(wmass/hmass)**2
      temp=Isubw(arg)
      do j=1,2
      temp=temp+xn*Q(j)*(Isubq(x(j))+Isubq(x(j+2))+Isubq(x(j+4)))
      enddo
      gamwidth=fac*abs(temp)**2
      return
      end

      double complex function Isubq(x)
      implicit none
      double precision x
      double complex Fgaga
      Isubq=4d0*x*(2d0+(4d0*x-1d0)*Fgaga(x))
      return
      end

      double complex function Isubw(x)
      implicit none
      double precision x
      double complex Fgaga
      Isubw=-2d0*(6d0*x+1d0+6d0*x*(2d0*x-1d0)*Fgaga(x))
      return
      end

      double complex function Fgaga(x)
      implicit none
      include 'constants.f'
      double precision x,fourx,arg
      fourx=4d0*x
      if (fourx .ge. 1d0) then
          arg=sqrt(0.25d0/x)
          Fgaga=-2d0*asin(arg)
      elseif ((fourx .lt. 1d0) .and. (fourx .gt. 0d0)) then
          arg=sqrt(1d0-fourx)
          arg=(1d0+arg)/(1d0-arg)
          Fgaga=dcmplx(0.5d0*log(arg),pi)
      else 
          Fgaga=czip
      endif
      return
      end
