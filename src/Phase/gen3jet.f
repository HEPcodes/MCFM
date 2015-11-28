      subroutine gen3jet(r,p,wt3,*)
C---generate three particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4+p5
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      integer j,nu

      double precision r(mxdim),p(mxpart,4),xx(2),
     . sqrts,ymax,ymin,xjac,y3,y4,y5,phi,phi45,wt0,wt3,
     . pt3,pt4,pt5,xt3,xt4,xt5,xmin,rtson2,cphi,sphi,cphi45,sphi45,
     . deltay
      common/energy/sqrts
      parameter(xmin=1d-5,wt0=1d0/512d0/pi**3)
      common/x1x2/xx

      do j=1,mxpart     
      do nu=1,4     
      p(j,nu)=0d0
      enddo     
      enddo     

      phi=2d0*pi*r(1)
      cphi=dcos(phi)
      sphi=dsin(phi)
      phi45=2d0*pi*r(2)
      cphi45=dcos(phi45)
      sphi45=dsin(phi45)
      xjac=sqrts**2
      ymax=5d0
      ymin=-5d0
      Deltay=ymax-ymin	  
      y3=ymin+Deltay*r(3)
      y4=ymin+Deltay*r(4)
      y5=ymin+Deltay*r(5)
      
      xjac=xjac*Deltay**3

      xt4=r(6)
      xt5=r(7)
      xjac=xjac*xt4*xt5
      rtson2=0.5d0*sqrts
      pt4=rtson2*xt4
      pt5=rtson2*xt5

      p(4,1)=rtson2*xt4*sphi
      p(4,2)=rtson2*xt4*cphi

      p(5,1)=rtson2*xt5*(+cphi45*sphi+sphi45*cphi)
      p(5,2)=rtson2*xt5*(-sphi45*sphi+cphi45*cphi)


      p(3,1)=-p(4,1)-p(5,1)
      p(3,2)=-p(4,2)-p(5,2)
      pt3=dsqrt(p(3,1)**2+p(3,2)**2)
      xt3=pt3/rtson2

      xx(1)=half*(+xt3*exp(+y3)+xt4*exp(+y4)+xt5*exp(+y5))
      xx(2)=half*(+xt3*exp(-y3)+xt4*exp(-y4)+xt5*exp(-y5))

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
c      write(6,*) 'problems with xx(1),xx(2) in gen3',xx(1),xx(2)  
      return 1 
      endif
	  
      p(1,4)=-0.5d0*xx(1)*sqrts
      p(1,1)=0d0
      p(1,2)=0d0
      p(1,3)=-0.5d0*xx(1)*sqrts
      
      p(2,4)=-0.5d0*xx(2)*sqrts
      p(2,1)=0d0
      p(2,2)=0d0
      p(2,3)=+0.5d0*xx(2)*sqrts

      p(3,4)=+pt3*cosh(y3)
      p(3,3)=+pt3*sinh(y3)

      p(4,4)=+pt4*cosh(y4)
      p(4,3)=+pt4*sinh(y4)

      p(5,4)=+pt5*cosh(y5)
      p(5,3)=+pt5*sinh(y5)

      wt3=wt0*xjac
      return

      end
