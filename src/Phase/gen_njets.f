      subroutine gen_njets(r,njets,p,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'process.f'
      include 'limits.f'
      include 'cutoff.f'
c---- generate phase space for 2-->2+n process
c---- with (34) being a vector boson and 5,..,4+n the jets
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
      double precision r(mxdim)
      double precision p(mxpart,4),p3(4),p34(4),psumjet(4),pcm(4),Q(4)
      double precision wt,wt3456,wt34,wt56,wt0,xxmin
      double precision xmin,xmax,delx,x,sqrts,pt,ymax,ymin,xx(2)
      double precision y,sinhy,coshy,phi,mv2,wtbw,mjets
      double precision ybar,ptsumjet2,shat,acm,ycm,sumpst,q0st,rshat
      double precision costh,sinth,dely
      double precision ptjetmin,yjetmax
      double precision plstar,estar,plstarsq,y5starmax,y5starmin
      double precision mass2,width2,mass3,width3
      integer j,nu,njets,ijet,n2,n3
      logical first
      character*4 part
      common/part/part
      common/energy/sqrts
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/x1x2/xx
      parameter(wt0=1d0/twopi**2)
      parameter(xxmin=1d-5)
      data first/.true./
      save first,ptjetmin,yjetmax

      if (first) then
        first=.false.
        call read_jetcuts(ptjetmin,yjetmax)
        if (part .eq. 'real') then
          ptjetmin=dsqrt(cutoff)
          yjetmax=10d0
        endif
      endif        

      do nu=1,4
        do j=1,4+njets
          p(j,nu)=0d0
        enddo
        psumjet(nu)=0d0
        pcm(nu)=0d0
      enddo 

      wt=2d0*pi
            
      do ijet=1,njets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16d0/pi**3
        xmin=2d0/sqrts
        xmax=1d0/ptjetmin
        delx=xmax-xmin
        x=xmin+r(ijet)*delx
        pt=1d0/x
        ymax=sqrts/2d0/pt
        ymax=dlog(ymax+dsqrt(ymax**2-1d0))
        
        ymax=min(ymax,yjetmax)
        y=ymax*(2d0*r(njets+ijet)-1d0)
        wt=wt*2d0*ymax
        
        sinhy=dsinh(y)
        coshy=dsqrt(1+sinhy**2)
        
        p(4+ijet,4)=pt*coshy
        wt=wt*delx*pt**3
        
        phi=2d0*pi*r(2*njets+ijet)
        wt=wt*2d0*pi
        
        p(4+ijet,1)=pt*dcos(phi)
        p(4+ijet,2)=pt*dsin(phi)
        p(4+ijet,3)=pt*sinhy
        
        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(4+ijet,nu)
        enddo
      enddo
c--- now generate Breit-Wigner        
      call breitw(r(3*njets+1),wsqmin,wsqmax,mass3,width3,mv2,wtbw)
      wt=wt*wtbw/2d0/pi
c--- invariant mass of jets
      mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
      mjets=dsqrt(abs(mjets))
      
      ybar=0.5d0*dlog((psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3)))
      ptsumjet2=psumjet(1)**2+psumjet(2)**2
      plstarsq=((sqrts**2-mv2-mjets**2)**2
     . -4d0*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4d0*sqrts**2)
      if (plstarsq .le. 0d0) then
        wt=0d0
        return 1
      endif
      plstar=dsqrt(plstarsq)
      Estar=dsqrt(plstarsq+ptsumjet2+mjets**2)
      y5starmax=0.5d0*dlog((Estar+plstar)/(Estar-plstar))
      y5starmin=-y5starmax

      ymax=ybar-y5starmin
      ymin=ybar-y5starmax
      dely=ymax-ymin
      ycm=ymin+r(3*njets+2)*dely     
      sinhy=dsinh(ycm)
      coshy=dsqrt(1d0+sinhy**2)
      
c--- now make the initial state momenta
      sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
      q0st=dsqrt(mv2+sumpst)
      rshat=q0st+dsqrt(mjets**2+sumpst)
      pcm(4)=rshat*coshy
      pcm(3)=rshat*sinhy
            
      xx(1)=(pcm(4)+pcm(3))/sqrts
      xx(2)=(pcm(4)-pcm(3))/sqrts
      
      if   (xx(1)*xx(2) .gt. 1d0) then
      write(6,*) 'gen_njetsL:xx(1)*xx(2),xx(1),xx(2)',
     . xx(1)*xx(2),xx(1),xx(2)  
      endif
      if   ((xx(1) .gt. 1d0) .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin).or. (xx(2) .lt. xmin)) then
         wt=0d0
        return 1
      endif 
      
      wt=wt*dely
      do j=1,4
        Q(j)=pcm(j)-psumjet(j)
      enddo
      
      p(1,4)=-xx(1)*sqrts/2d0
      p(1,3)=p(1,4)
      p(2,4)=-xx(2)*sqrts/2d0
      p(2,3)=-p(2,4)
      
      wt=wt*rshat/(sqrts**2*q0st)
c--- decay boson into leptons, in boson rest frame
      costh=2d0*r(3*njets+3)-1d0
      sinth=dsqrt(1d0-costh**2)
      phi=2d0*pi*r(3*njets+4)
      p34(4)=dsqrt(mv2)/2d0
      p34(1)=p34(4)*sinth*dcos(phi)
      p34(2)=p34(4)*sinth*dsin(phi)
      p34(3)=p34(4)*costh
      
c--- boost into lab frame    
      call boost(dsqrt(mv2),Q,p34,p3)
      do j=1,4
      p(3,j)=p3(j)
      p(4,j)=Q(j)-p(3,j)
      enddo

      wt=wt/8d0/pi
            
      return
      end
      
      
      
      
      
      
      
      
      
      
      
