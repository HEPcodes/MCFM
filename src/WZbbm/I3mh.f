      double complex function I3Mh(k1sq,k2sq,k3sq,mmsq)
      implicit none
C     Triangle function with one internal mass and two spacelike
C     offshell lines and one timelike of mass k3sq
C 
c                         || k2
c                         ||
c                         /\
c                        /__\
C               k1  *****----$$$$$$ k3
C                         m
C
c
C   Adapted from
c   %\cite{Denner:kt}
c   \bibitem{Denner:kt}
c   A.~Denner,
c   %``Techniques For Calculation Of Electroweak Radiative 
C   Corrections At The One Loop Level And Results For W Physics At Lep-200,''
c   Fortsch.\ Phys.\  {\bf 41}, 307 (1993).
c   %%CITATION = FPYKA,41,307;%%

      include 'constants.f'
      integer i,j,k,jj(0:2),kk(0:2)
      double precision x,y,z,k1sq,k2sq,k3sq,mmsq,msq(0:2),htheta,epsilon
      double complex lambda,alpha,y0(0:2),xp(0:2),xm(0:2),polylog,
     . yip(0:2),yim(0:2),alphai(0:2),temp
      double precision psq(0:2,0:2),theta,eta
      double complex a,b
      data jj/1,2,0/
      data kk/2,0,1/
      data epsilon/1d-15/
      lambda(x,y,z)=sqrt(x**2+y**2+z**2-2d0*(x*y+y*z+z*x))
C--- define Heaviside htheta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+half*sign(one,x)

      eta(a,b)=dcmplx(0d0,twopi*(
     . +htheta(-dimag(a))*htheta(-dimag(b))*htheta(+dimag(a*b))
     . -htheta(+dimag(a))*htheta(+dimag(b))*htheta(-dimag(a*b))))
      
c      if ((k1sq .ge. 0d0).or.(k3sq .ge. 0d0)) then 
c      write(6,*) 'k1sq,k2sq,k3sq in I3mh',k1sq,k2sq,k3sq
c      stop
c      endif

      I3mh=0d0
      psq(0,1)=k1sq
      psq(1,2)=k2sq
      psq(2,0)=k3sq
      msq(0)=mmsq
      msq(1)=0d0
      msq(2)=0d0

      alpha=lambda(psq(0,1),psq(1,2),psq(2,0))
      do i=0,2
      j=jj(i)
      k=kk(i)
      alphai(i)=lambda(psq(j,k),msq(j),msq(k))
     . *(1+dcmplx(0d0,sign(epsilon,psq(j,k))))
      y0(i)=0.5d0/alpha/psq(j,k)
     . *(psq(j,k)*(psq(j,k)-psq(k,i)-psq(i,j)
     . +2d0*msq(i)-msq(j)-msq(k))
     . -(psq(k,i)-psq(i,j))*(msq(j)-msq(k))
     . +alpha*(psq(j,k)-msq(j)+msq(k)))
      xp(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)+alphai(i))
      xm(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)-alphai(i))
      yip(i)=y0(i)-xp(i)
      yim(i)=y0(i)-xm(i)
      I3mh=I3mh+polylog(y0(i),yip(i))+polylog(y0(i),yim(i))
     .+eta(1d0-xp(i),1d0/yip(i))*log((y0(i)-1d0)/yip(i))
     .-eta(-xp(i),1d0/yip(i))*log(y0(i)/yip(i))
     .+eta(1d0-xm(i),1d0/yim(i))*log((y0(i)-1d0)/yim(i))
     .-eta(-xm(i),1d0/yim(i))*log(y0(i)/yim(i))
      temp=
     .-(eta(-xp(i),-xm(i))-eta(yip(i),yim(i))
     . -twopi*Im*htheta(-psq(j,k))*htheta(-dimag(yip(i)*yim(i))))
      if (temp.eq.czip) goto 20
      I3mh=I3mh+temp*(log(-(1d0-y0(i))/y0(i)))
 20   enddo 
C    minus sign for (-1)^n
      I3mh=-I3mh/alpha
      return
      end

      double complex function polylog(y,z)
      implicit none
      include 'constants.f'
      double precision rarg1,rarg2
      double complex y,z,arg1,arg2,dilog1,dilog2,wgplg
      arg1=(y-1d0)/z
      arg2=y/z
      rarg1=dble(arg1)
      rarg2=dble(arg2)
      if (rarg1 .gt. 1d0) then
        dilog1=-wgplg(1,1,1d0-rarg1)+pisqo6-log(rarg1)*log(1d0-arg1)
      else 
        dilog1=wgplg(1,1,rarg1)
      endif
      if (rarg2 .gt. 1d0) then
        dilog2=-wgplg(1,1,1d0-rarg2)+pisqo6-log(rarg2)*log(1d0-arg2)
      else 
        dilog2=wgplg(1,1,rarg2)
      endif
      polylog=dilog1-dilog2
      return
      end

