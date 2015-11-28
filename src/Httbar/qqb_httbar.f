      subroutine qqb_httbar(p,msq)
      implicit none
c----Lowest order matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H --> ttbar--->  W^+ (nu(p3)+e^+(p4))+b(p5)
c                                      ++bbar(p6)+W^- (e^-(p7)+nubar(p8))
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k,nu
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),s,s12
      double precision decay,aw,gg,rat,fac,
     . s4x,s7x,xDx,qDq,aDa,pq(4),pa(4),px(4),props,Ifunsq
C--statment function
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c---set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      pq(nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      pa(nu)=-p(6,nu)-p(7,nu)-p(8,nu)
      px(nu)=pq(nu)+pa(nu)
      enddo
      qDq=pq(4)**2-pq(1)**2-pq(2)**2-pq(3)**2
      aDa=pa(4)**2-pa(1)**2-pa(2)**2-pa(3)**2
      xDx=px(4)**2-px(1)**2-px(2)**2-px(3)**2
      s4x=2d0*(p(4,4)*px(4)-p(4,1)*px(1)-p(4,2)*px(2)-p(4,3)*px(3))
      s7x=2d0*(p(7,4)*px(4)-p(7,1)*px(1)-p(7,2)*px(2)-p(7,3)*px(3))
      props=
     . +((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
     . *((s(7,8)-wmass**2)**2+(wmass*wwidth)**2)
     . *((qDq-mt**2)**2+(mt*twidth)**2)
     . *((aDa-mt**2)**2+(mt*twidth)**2)
C---matrix element squared for ttbar decay
      decay=xn*gwsq**5*mt**4/(4d0*wmass**2)
     .  *s(3,5)*s(4,6)*(s4x*s7x-xDx*s(4,7))/props

C---production cross section in the heavy top limit
      aw=gwsq/(4d0*pi)
      s12=s(1,2)
      rat=(mt**2/s12)
      fac=Ifunsq(rat)
      gg=aw*as**2*V/(18d0*pi)*(s12/wmass)**2
     . /((s12-hmass**2)**2+(hmass*hwidth)**2)
c---calculate propagators
      msq(0,0)=avegg*gg*fac*decay
      return
      end

      double precision function Ifunsq(x)
      implicit none
      double precision x
      double complex ffun,Ifun
      Ifun=3d0*x*(2d0+(4d0*x-1d0)*ffun(x))
      Ifunsq=abs(Ifun)**2
      return
      end

      double complex function ffun(x)
      implicit none
      double precision x,pi,discr,arg
      pi=2d0*asin(1d0)
      discr=1d0-4d0*x
      if (discr .ge. 0d0) then
      arg=log((1d0+sqrt(discr))/(1d0-sqrt(discr)))
      ffun=0.5d0*(dcmplx(arg,-pi))**2
      elseif (discr .lt. 0d0) then
      arg=0.5d0/sqrt(x)
      ffun=-2d0*asin(arg)**2
      endif
      return
      end
