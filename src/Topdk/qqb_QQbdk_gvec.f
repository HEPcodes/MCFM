      subroutine qqb_QQbdk_gvec(p,n,in,msq)
      implicit none
c----Matrix element for tt~ production
C----averaged over initial colours and spins
c    contracted with the vector n(mu)
C***********************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2002.                                                      *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     g(-p1) +g(-p2)=t(nu(p3)+e^+(p4)+b(p5))                           *
*                   +t~(b~(p6)+e^-(p7)+nu(p8))                         *
*                                                                      * 
************************************************************************

c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'msqv_cs.f'
      integer j,k,in,nu
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),ps(mxpart,4)
      double precision n(4),fac,msqn(0:2)
     
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
      fac=gwsq**4*s35*s68/((mt*twidth)**4*propw)*gsq**2*V/4d0/xn
c      we will have no further need for p3 and p5 
c      we will have no further need for p6 and p8 
      do nu=1,4
c t momentum
      ps(3,nu)=p(3,nu)+p(4,nu)+p(5,nu)
c positron
      ps(4,nu)=p(4,nu)
c t-bar momentum
      ps(5,nu)=p(6,nu)+p(7,nu)+p(8,nu)
c vector n
      ps(6,nu)=n(nu)
c elctron
      ps(7,nu)=p(7,nu)
      enddo      
      
      if (in .eq. 1) then
      call ttbdkn(1,2,ps,msqn)
      elseif (in .eq. 2) then
      call ttbdkn(2,1,ps,msqn)
      endif
      do j=0,2
      msqv_cs(j,0,0)=avegg*fac*msqn(j)
      enddo 
      msq(0,0)=msqv_cs(0,0,0)+msqv_cs(1,0,0)+msqv_cs(2,0,0)
      return
      end
 
