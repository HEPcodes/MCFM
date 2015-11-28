      subroutine gen6(r,q,wt6,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'phasemin.f'
      integer nu
      double precision r(mxdim)
      double precision wt6,q(mxpart,4)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      double precision pswt,xjac,p1ext(4),p2ext(4)
      double precision xx(2),tau,x1mx2,surd
      double precision lntaum
      common/pext/p1ext,p2ext
      common/x1x2/xx
      data p3/0d0,0d0,0d0,0d0/

      wt6=0d0

      lntaum=dlog(taumin)
      tau=dexp(lntaum*(one-r(9)))
      xjac=-lntaum*tau

c      tau=(one-taumin)*r(14)**2+taumin
c      xjac=2*r(13)*(one-taumin)

      x1mx2=two*r(10)-one
      surd=dsqrt(x1mx2**2+four*tau) 
           
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


      call phase6(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999) 
c      write(6,*) 'p1sq',p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
c      write(6,*) 'p2sq',p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
c      write(6,*) 'p4sq',p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
c      write(6,*) 'p4sq',p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
c      write(6,*) 'p5sq',p5(4)**2-p5(1)**2-p5(2)**2-p5(3)**2
c      write(6,*) 'p6sq',p6(4)**2-p6(1)**2-p6(2)**2-p6(3)**2
c      write(6,*) 'p7sq',p7(4)**2-p7(1)**2-p7(2)**2-p7(3)**2
c      write(6,*) 'p8sq',p8(4)**2-p8(1)**2-p8(2)**2-p8(3)**2

c      write(6,*) 'p34',2d0*(p3(4)*p4(4)-p3(3)*p4(3)
c     . -p3(2)*p4(2)-p3(1)*p4(1))

c      write(6,*) 'p78',2d0*(p8(4)*p7(4)-p8(3)*p7(3)
c     . -p8(2)*p7(2)-p8(1)*p7(1))
c      pause

      do nu=1,4
      q(1,nu)=p1(nu)
      q(2,nu)=p2(nu)
      q(3,nu)=p3(nu)
      q(4,nu)=p4(nu)
      q(5,nu)=p5(nu)
      q(6,nu)=p6(nu)
      q(7,nu)=p7(nu)
      q(8,nu)=p8(nu)

      enddo 
      wt6=xjac*pswt

      if (debug) write(6,*) 'wt6 in gen6',wt6
      return

 999  return 1
      end

