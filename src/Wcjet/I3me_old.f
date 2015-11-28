      double precision function I3Me(k1sq,k2sq,k3sq)
      implicit none
C     Triangle function with one internal mass and two spacelike
C     offshell lines and one timelike of mass k1sq
C 
c                         | k2
c                         |
c                         /\
c                        /__\____
C               k3  *****-------- k1,(k1sq=msq)  
C
C   Mass of the one massive internal line is equal to k1sq 
c
C   Adapted from
c   %\cite{Denner:kt}
c   \bibitem{Denner:kt}
c   A.~Denner,
c   %``Techniques For Calculation Of Electroweak Radiative 
C   Corrections At The One Loop Level And Results For W Physics At Lep-200,''
c   Fortsch.\ Phys.\  {\bf 41}, 307 (1993).
c   %%CITATION = FPYKA,41,307;%%

      integer i,j,k,jj(0:2),kk(0:2)
      double precision lambda,alpha,y0(0:2),xp(0:2),xm(0:2),
     . yip(0:2),yim(0:2),alphai(0:2),x,y,z,k1sq,k2sq,k3sq
      double precision psq(0:2,0:2),msq(0:2),polylog,ddilog
      data jj/1,2,0/
      data kk/2,0,1/
      polylog(y,z)=ddilog((y-1d0)/z)-ddilog(y/z)
      lambda(x,y,z)=sqrt(x**2+y**2+z**2-2d0*(x*y+y*z+z*x))
      
      if ((k2sq .ge. 0d0).or.(k3sq .ge. 0d0)) then 
      write(6,*) 'k1sq,k2sq,k3sq in I3me',k1sq,k2sq,k3sq
      stop
      endif

      I3me=0d0
      psq(0,1)=k1sq
      psq(1,2)=k2sq
      psq(2,0)=k3sq
      msq(0)=k1sq
      msq(1)=0d0
      msq(2)=0d0

      alpha=lambda(psq(0,1),psq(1,2),psq(2,0))
      do i=0,2
      j=jj(i)
      k=kk(i)
      alphai(i)=lambda(psq(j,k),msq(j),msq(k))
      y0(i)=0.5/alpha/psq(j,k)
     . *(psq(j,k)*(psq(j,k)-psq(k,i)-psq(i,j)
     . +2d0*msq(i)-msq(j)-msq(k))
     . -(psq(k,i)-psq(i,j))*(msq(j)-msq(k))
     . +alpha*(psq(j,k)-msq(j)+msq(k)))
      xp(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)+alphai(i))
      xm(i)=0.5d0/psq(j,k)*(psq(j,k)-msq(j)+msq(k)-alphai(i))
      yip(i)=y0(i)-xp(i)
      yim(i)=y0(i)-xm(i)
      I3me=I3me+polylog(y0(i),yip(i))+polylog(y0(i),yim(i))
      enddo 
C---- Minus sign for (-1)^n
      I3me=-I3me/alpha
      return
      end

