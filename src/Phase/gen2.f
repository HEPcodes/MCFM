      subroutine gen2(r,p,wt2,*)
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4
      implicit none
      include 'constants.f'
      include 'limits.f'
      include 'mxdim.f'
      integer n2,n3,j,nu

      double precision r(mxdim),p(mxpart,4),xx(2)
      double precision sqrts,ymax,yave,ydif,xjac,y3,y4,phi,wt0,wt2,w3
      double precision pt,s34,xmin,rtshat,udif
      common/energy/sqrts
      double precision mass2,width2,mass3,width3 
      common/breit/n2,n3,mass2,width2,mass3,width3 
      parameter(xmin=1d-5,wt0=1d0/16d0/pi)
      common/x1x2/xx

      do j=1,mxpart     
      do nu=1,4     
      p(j,nu)=0d0
      enddo     
      enddo     

      wt2=0d0
      if (n3.eq.0) then
         w3=(wsqmax-wsqmin)
         s34=(wsqmax-wsqmin)*r(3)+wsqmin
      elseif (n3.eq.1) then 
         call breitw(r(3),wsqmin,wsqmax,mass3,width3,s34,w3)
      endif

      rtshat=dsqrt(s34)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(1)-1d0)
      
c----udif==tanh(ydif)
      udif=(two*r(2)-1d0)
      ydif=half*dlog((1d0+udif)/(1d0-udif))
      xjac=four*ymax
	  
      y3=yave+ydif
      y4=yave-ydif
	  
      xjac=xjac*w3
      phi=2d0*pi*r(4)

      pt=rtshat/(2d0*cosh(ydif))
      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
c      write(6,*) 'problems with xx(1),xx(2) in gen2',xx(1),xx(2)  
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

      p(3,4)=+pt*cosh(y3)
      p(3,1)=+pt*sin(phi)
      p(3,2)=+pt*cos(phi)
      p(3,3)=+pt*sinh(y3)

      p(4,4)=+pt*cosh(y4)
      p(4,1)=-pt*sin(phi)
      p(4,2)=-pt*cos(phi)
      p(4,3)=+pt*sinh(y4)

      wt2=wt0*xjac/sqrts**2
      return

      end
