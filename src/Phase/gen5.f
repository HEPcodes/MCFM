      subroutine gen5(r,p,wt5,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      integer nu

      double precision r(mxdim),sqrts,wt5,
     . p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision pswt,xjac
      double precision xx(2),xmin,tau,y,taumin

      character*6 case
      common/process/case
      common/energy/sqrts
      common/x1x2/xx
      common/taumin/taumin
      parameter(xmin=1d-5)

      wt5=0d0
      tau=exp(log(taumin)*r(9))
      y=0.5d0*log(tau)*(1d0-2d0*r(10))
      xjac=log(taumin)*tau*log(tau)

      xx(1)=sqrt(tau)*exp(+y)
      xx(2)=sqrt(tau)*exp(-y)

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


      if (case .eq. 't_bbar') then
      write(6,*) 'higher order not implemented yet'
      stop
      elseif (case .eq. 'qg_tbb') then
      call phase51(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      elseif (case .eq. 'vlchk5')  then  
      call phase5(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      else
      call phase5(r,p1,p2,p3,p4,p5,p6,p7,pswt)
      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      enddo 
      wt5=xjac*pswt
      if(wt5 .eq. 0) return 1
      return
      end
