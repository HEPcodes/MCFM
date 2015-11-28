      subroutine phi1_2m(m2,x3,xth,xphi,s3min,p1,p2,p3,wt,*)
c     massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
c     with invariant mass of particle three s3 integrated over.
c     s3min is the minimum value of s3.
c     Vectors returned p2 and p3 are in the same frame as p1 is supplied.
c     Expression evaluated is 
c     ds2 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-s2) delta(p3^2-s3)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'debug.f'
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x3,xth,xphi,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w3
      double precision s3max,s3min,xx,xexp
      double precision m1,m2,m3,s1,s2,s3,lambda,xjac,rtxth,
     . mass2,width2,mass3,width3
c      double precision Eg
      integer j,n2,n3
      integer jbranch
c      logical first
      common/breit/n2,n3,mass2,width2,mass3,width3
      parameter(wt0=one/8.d0/pi)
      data jbranch/1/
c      data first/.true./
      save jbranch,xexp
c      if (first) then
c      first=.false.
c      write(6,*) 'Enter exponent for reweighting'
c      read(5,*) xexp
c      write(6,*) 'xexp',xexp
c      endif
      xexp=1d0
      wt=0d0
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) return 1
      m1=sqrt(s1)
      s2=m2**2
      s3max=(m2-m1)**2
      if (s3min .gt. s3max) return 1
      if (n3 .eq. 0) then
         w3=s3max-s3min
         s3=s3max*x3+s3min*(1d0-x3)
         xx=0
      elseif (n3 .eq. 1) then
        xx=1
        call breitw(x3,s3min,s3max,mass3,width3,s3,w3) 
      endif

      m3=sqrt(s3)


      if (jbranch .eq. 1) then
        jbranch=2
        rtxth=xth**xexp 
        xjac=1d0/(xexp*xth**(xexp-1d0))
      elseif (jbranch .eq. 2) then 
        jbranch=1
        rtxth=1d0-xth**xexp
        xjac=1d0/(xexp*xth**(xexp-1d0))
      endif



      costh=two*rtxth-one      
      phi=twopi*xphi
      sinth=sqrt(one-costh**2)
      cphi=cos(phi)
      sphi=sin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)
      if ((lambda .lt. 0d0) .or. debug) then
      write(6,*) 'lambda in phi1_2m',lambda
      write(6,*) 's1 in phi1_2m',s1
      write(6,*) 's2 in phi1_2m',s2
      write(6,*) 's3 in phi1_2m',s3
      write(6,*) 'm1 in phi1_2m',m1
      write(6,*) 'm2 in phi1_2m',m2
      write(6,*) 'm3 in phi1_2m',m3
      write(6,*) 'm1-m2-m3 in phi1_2m',m1-m2-m3
      write(6,*) 'xx in phi1_2m',xx
      write(6,*) 'x3 in phi1_2m',x3
      write(6,*) 'n3 in phi1_2m',n3
      write(6,*) 'mass3 in phi1_2m',mass3
      return 1
      endif
      lambda=dsqrt(lambda)
c      Eg=s1+s2-s3

      wt=wt0*w3*lambda/s1/xjac

c      write(6,*) 'xjac',xjac
c      write(6,*) 'wt in phi1_2m',wt
c      pause

      if(debug) write(6,*) 'wt in phi1_2m',wt
c      write(6,*) 'xjac',xjac
c      write(6,*) 'lambda',lambda
c      write(6,*) 'Eg',Eg
c      write(6,*) 'xjac',xjac
c      write(6,*) 'wt in phi1_2m',wt
c      write(6,*) jbranch
c      pause

      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh


      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo


c      write(6,*) '13',p1(4)*p3(4)-p1(3)*p3(3)-p1(2)*p3(2)-p1(2)*p3(2)
c      write(6,*) '23',p2(4)*p3(4)-p2(3)*p3(3)-p2(2)*p3(2)-p2(2)*p3(2)
c      write(6,*) 'costh',costh
c      write(6,*) 'xth',xth
c      pause
      if (  (p1(4) .lt. 0.d0) 
     & .or. (p2(4) .lt. 0.d0) 
     & .or. (p3(4) .lt. 0.d0)) then  
c      write(6,*) 'p1(4)',p1(4)
c      write(6,*) 'p2(4)',p2(4)
c      write(6,*) 'p3(4)',p3(4)
c      write(6,*) 'p1sq',p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
c      write(6,*) 'p2sq',p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,
c      write(6,*) 'p3sq',p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
c      write(6,*) 'in phi1_2m.f'
c      write(6,*) n2,n3
      return 1
      endif
      return
      end



