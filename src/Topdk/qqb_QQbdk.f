      subroutine qqb_QQbdk(p,msq)
      implicit none
c----Matrix element for tt~ production
C----averaged over initial colours and spins
C***********************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2002.                                                      *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     g(-p1) +g(-p2)=t(nu(p3)+e^+(p4)+b(p5))                           *
*                   +t~(b~(p6)+e^-(p7)+nu(p8))                         *
*     t and tbar are assumed to be on shell                            * 
************************************************************************

      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'msq_cs.f'
      integer j,k,in,nu
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),ps(mxpart,4)
      double precision fac,wtgg,wtqqb,wtqbq,ro,t1,t2,s12,nab1,nab2,ab
      double precision test
     
      double precision s34,s78,s35,s68,propw

C---set all matrix elements equal to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      do j=1,2
      ps(j,nu)=p(j,nu)
      enddo
      do j=3,mxpart
      ps(j,nu)=zip
      enddo
      enddo
      
      s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1)) 
      s78=2d0*(p(7,4)*p(8,4)-p(7,3)*p(8,3)-p(7,2)*p(8,2)-p(7,1)*p(8,1)) 
      s35=2d0*(p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1))
      s68=2d0*(p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1)) 
      propw=((s34-wmass**2)**2+(wmass*wwidth)**2)
     .     *((s78-wmass**2)**2+(wmass*wwidth)**2)
      fac=gwsq**4*s35*s68/((mt*twidth)**4*propw)*gsq**2*V/16d0/xn

c      we will have no further need for p3 and p5 
c      we will have no further need for p6 and p8 

      do nu=1,4
c t momentum
      ps(3,nu)=p(3,nu)+p(4,nu)+p(5,nu)
c t-bar momentum
      ps(4,nu)=p(6,nu)+p(7,nu)+p(8,nu)
c positron
      ps(5,nu)=p(4,nu)
c electron
      ps(6,nu)=p(7,nu)
      enddo      
      

      call dotem(6,ps,s)
      s12=s(1,2)
      ro=4d0*mt**2/s12
      t1=-s(2,4)/s12
      t2=-s(2,3)/s12


cId,bit=g^4*V/N*(V/t1/t2-2*N^2)*(
c       + 4*p3Dp4*p7Dp6*(t1^2+t2^2)
c       + ro*(2*t2*p1Dp4*p2Dp7+2*t1*p1Dp7*p2Dp4
c         -p3Dp4*p3Dp7-p6Dp4*p6Dp7
c         -p3Dp4*p6Dp7*(-2+1/t1/t2))

c       + 1/2*ro^2/t1/t2 * (
c           -p1Dp4*p2Dp7*t2-p1Dp7*p2Dp4*t1
c            +p3Dp4*p3Dp7 
c            +p6Dp4*p6Dp7 )
c       -(1-1/2*ro/t1/t2+ro^2/t1^2/t2^2/8)*t1*t2*ro*psDps*p4Dp7 )

                 
      msq_cs(0,0,0)=16d0*s(3,5)*(-s(4,6))*(1d0/t1/t2-2d0)
     . +4d0*ro/t1/t2*((s(1,5)+s(2,5))*(s(1,6)+s(2,6))-s(3,6)*s(4,5))
     . -8d0*ro*(s(2,5)*s(1,6)/t2+s(1,5)*s(2,6)/t1-s12*s(5,6))
     . +4d0*ro/t1**2/t2**2*(3*t1*t2-1)*(-s(3,5)*s(4,6))
     . -2d0*(ro/(t1*t2))**2
     . *(s(2,5)*s(2,6)+s(1,5)*s(1,6)-s(4,5)*s(3,6)-s(3,5)*s(4,6)
     . +s(2,5)*s(1,6)*t2+s(1,5)*s(2,6)*t1+(2d0*t1*t2-ro/2d0)*s12*s(5,6))

      msq_cs(0,0,0)=fac*avegg*msq_cs(0,0,0)

      msq_cs(1,0,0)=-16d0*s(4,6)*s(3,5)*t1/t2*(2d0*t1*t2-1d0)
     . -4d0*ro*t1/t2*(s(2,5)*(s(2,6)+s(3,6))+s(1,5)*(s(1,6)+s(3,6))
     . +s(3,5)*s(3,6))
     . +4d0*ro/t2*(s(2,5)*s(1,6)-s(2,6)*s(1,5))*t1*(t1-t2)
     . -4d0*ro/t2**2*(1d0-3d0*t1*t2)*s(4,6)*s(3,5)
     . -8d0*ro*t1**2*s(5,6)*s12
     . +2d0*(ro/t2)**2*(s(2,5)*s(2,6)+s(1,5)*s(1,6)
     . +t2*s(2,5)*s(1,6)+t1*s(1,5)*s(2,6)
     . -s(4,6)*s(3,5)-s(4,5)*s(3,6))
     . +4d0*t1/t2*s(5,6)*ro**2*s12*(1d0-ro/4d0/t1/t2)

      msq_cs(1,0,0)=fac*xn**2*avegg*msq_cs(1,0,0)

      msq_cs(2,0,0)=-16d0*s(4,6)*s(3,5)*t2/t1*(2d0*t2*t1-1d0)
     . -4d0*ro*t2/t1*(s(1,5)*(s(1,6)+s(3,6))+s(2,5)*(s(2,6)+s(3,6))
     . +s(3,5)*s(3,6))
     . +4d0*ro/t1*(s(1,5)*s(2,6)-s(1,6)*s(2,5))*t2*(t2-t1)
     . +4d0*ro/t1**2*(1d0-3d0*t2*t1)*(s(1,6)+s(2,6)+s(3,6))*s(3,5)
     . -8d0*ro*t2**2*s(5,6)*s12
     . +2d0*(ro/t1)**2*(s(1,5)*s(1,6)+s(2,5)*s(2,6)
     . +t1*s(1,5)*s(2,6)+t2*s(2,5)*s(1,6)
     . -s(4,6)*s(3,5)-s(4,5)*s(3,6))
     . +4d0*t2/t1*s(5,6)*ro**2*s12*(1d0-ro/4d0/t2/t1)

      msq_cs(2,0,0)=fac*xn**2*avegg*msq_cs(2,0,0)
      wtgg=msq_cs(0,0,0)+msq_cs(1,0,0)+msq_cs(2,0,0)

c      wtgg=fac*avegg*(V/(t1*t2)-2d0*xn**2)*(
c     .  + 16d0*s(3,5)*s(4,6)*(t1**2+t2**2)
c     .  + 4d0*ro*(2d0*t2*s(1,5)*s(2,6)+2d0*t1*s(1,6)*s(2,5)
c     .    -s(3,5)*s(3,6)-s(4,5)*s(4,6)
c     .    -s(3,5)*s(4,6)*(1d0/(t1*t2)-2d0))
c     .  + 2d0*ro**2/(t1*t2)*(
c     .      -s(1,5)*s(2,6)*t2-s(1,6)*s(2,5)*t1
c     .      +s(3,5)*s(3,6)+s(4,5)*s(4,6))
c     . -4d0*(2d0-(ro/(t1*t2)-(ro/(2d0*t1*t2))**2))*t1*t2*ro*s12*s(5,6))

      wtqqb=fac*aveqq*xn*(16d0*s(3,5)*s(4,6)*(t1**2+t2**2+ro/2d0)
     . -4d0*ro*(s(3,5)*(s(3,6)+s(4,6))+s(4,6)*(s(3,5)+s(4,5)))
     . +8d0*ro*(t2*s(1,5)*s(2,6)+t1*s(1,6)*s(2,5)-t1*t2*s12*s(5,6))
     . +2d0*ro**2*s12*s(5,6))
      wtqbq=wtqqb

      do j=-nf,nf
      k=-j      
      if (j .lt.0) then
      msq(j,k)=wtqbq
      elseif (j.gt.0) then
      msq(j,k)=wtqqb
      elseif (j.eq.0) then
      msq(0,0)=wtgg
      endif
      enddo
      return
      end
 
