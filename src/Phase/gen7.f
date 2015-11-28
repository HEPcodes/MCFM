      subroutine gen7(r,q,wt7,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'debug.f'
      integer nu

      double precision r(mxdim)
      double precision wt7,q(mxpart,4)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     . p9(4),pswt,xjac,p1ext(4),p2ext(4)
      double precision xx(2),xmin,tau,taumin,x1mx2,surd
      double precision lntaum
      character*6 case
      common/process/case
      common/pext/p1ext,p2ext
      common/x1x2/xx
      common/taumin/taumin
      parameter(xmin=1d-5)
      data p1,p2,p3,p4,p5,p6,p7,p8,p9/36*0d0/

      wt7=0d0

      lntaum=log(taumin)
      tau=exp(lntaum*(one-r(18)))
      xjac=-lntaum*tau

c      tau=(one-taumin)*r(14)**2+taumin
c      xjac=2*r(13)*(one-taumin)

      x1mx2=two*r(19)-one
      surd=sqrt(x1mx2**2+four*tau) 
           
      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)
      xjac=xjac*two/surd

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 1 

      do nu=1,4
      p1(nu)=xx(1)*p1ext(nu)
      p2(nu)=xx(2)*p2ext(nu)
      enddo

      if     (case .eq. 'qq_ttg') then
        call  phase7(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      elseif (case .eq. 'hlljet') then
        call phase7m(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999)
      else
        write(*,*) 'Bad process name in gen7.f'
        stop
      endif
      
      do nu=1,4
      q(1,nu)=p1(nu)
      q(2,nu)=p2(nu)
      q(3,nu)=p3(nu)
      q(4,nu)=p4(nu)
      q(5,nu)=p5(nu)
      q(6,nu)=p6(nu)
      q(7,nu)=p7(nu)
      q(8,nu)=p8(nu)
      q(9,nu)=p9(nu)
      q(10,nu)=0d0

      enddo 
      wt7=xjac*pswt
      if (debug) write(6,*) 'wt7 in gen7',wt7
      return

 999  return 1
      end

