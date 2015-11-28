      subroutine gen7(r,q,wt7,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'process.f'
      include 'phasemin.f'
      include 'debug.f'
      integer nu
      double precision r(mxdim)
      double precision wt7,q(mxpart,4)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     . p9(4),pswt,xjac
      double precision xx(2),tau,sqrts,y
      common/x1x2/xx
      common/energy/sqrts
      data p1,p2,p3,p4,p5,p6,p7,p8,p9/36*0d0/

      wt7=0d0

      tau=dexp(dlog(taumin)*r(9))
      y=0.5d0*dlog(tau)*(1d0-2d0*r(10))
      xjac=dlog(taumin)*tau*dlog(tau)

      xx(1)=dsqrt(tau)*dexp(+y)
      xx(2)=dsqrt(tau)*dexp(-y)

c--- phase space volume only checked for x1=x2=1
      if ((case .eq. 'vlchwg') .or. (case .eq. 'vlchwh')) then
        xx(1)=1d0
        xx(2)=1d0
        xjac=1d0
      endif

c---if x's out of normal range alternative return
      if   ((xx(1) .gt. 1d0)
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      if  ((case .eq. 'qq_HWW') .or. (case .eq. 'qq_HZZ')
     ..or. (case .eq. 'HWW2jt') .or. (case .eq. 'HZZ2jt')
     ..or. (case .eq. 'WpWp3j')) then
        call  phase7a(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      elseif ((case .eq. 'WH__WW') .or. (case .eq. 'ZH__WW')) then
        call  phase7b(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      elseif ((case .eq. 'WH__ZZ') .or. (case .eq. 'ZH__ZZ')) then
        call  phase7b(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,pswt,*999) 
      else
        write(6,*) 'Unanticipated process in gen7.f!'
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
      enddo 
      
      wt7=xjac*pswt
      
      if (debug) write(6,*) 'wt7 in gen7',wt7
      
      return

 999  return 1
      end

