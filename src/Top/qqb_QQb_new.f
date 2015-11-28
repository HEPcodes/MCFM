      subroutine qqb_QQb_new(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=t(nu(p3)+t~(p4)                                *
C                                                                      * 
C***********************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex m1(2,2,2,2),m2(2,2,2,2),mqqb(2,2,2),mqbq(2,2,2)
      integer nu,g3,g4,t1,t2,r1,r2,j,k,l1,l2,l3,l4
      parameter(g3=1,g4=2,t1=3,t2=5)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision q(mxpart,4),qqb,qbq
      double precision p1Dp2,p3Dq5,p4Dq6,mass,ggg,facqq,facgg
      mass=mt
      r1=1
      r2=1
      do nu=1,4
      q(1,nu)=p(1,nu)
      q(2,nu)=p(2,nu)
      q(5,nu)=p(r1,nu)
      q(6,nu)=p(r2,nu)
      enddo

      p3Dq5=p(3,4)*q(5,4)-p(3,1)*q(5,1)-p(3,2)*q(5,2)-p(3,3)*q(5,3)
      p4Dq6=p(4,4)*q(6,4)-p(4,1)*q(6,1)-p(4,2)*q(6,2)-p(4,3)*q(6,3)

      do nu=1,4
      q(3,nu)=p(3,nu)-0.5d0*mass**2/p3Dq5*q(5,nu)
      q(4,nu)=p(4,nu)-0.5d0*mass**2/p4Dq6*q(6,nu)
      enddo

      
c--- set up spinors
      call spinoru(6,q,za,zb)
      call qqamps(mqqb,mqbq)
      call hhhamp(m1,m2)
      facqq=aveqq*V/4d0*gsq**2
      facgg=avegg*V/4d0/xn*gsq**2
      qqb=0d0
      qbq=0d0
      ggg=0d0

      do l3=1,2
      do l2=1,2
      do l1=1,2
      qqb=qqb+facqq*abs(mqqb(l1,l2,l3))**2
      qbq=qbq+facqq*abs(mqbq(l1,l2,l3))**2
      do l4=1,2
      ggg=ggg+facgg*(
     . +xn**2*(cdabs(m1(l1,l2,l3,l4))**2+cdabs(m2(l1,l2,l3,l4))**2)
     . -cdabs(m1(l1,l2,l3,l4)+m2(l1,l2,l3,l4))**2)
      enddo
      enddo
      enddo
      enddo
     
      pause


      do j=-nf,nf
          k=-j
          if (j.gt.0) msq(j,k)=qqb
          if (j.eq.0) msq(j,k)=ggg
          if (j.lt.0) msq(j,k)=qbq
      enddo

      return
      end


      subroutine qqamps(mqqb,mqbq)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex mqqb(2,2,2),mqbq(2,2,2)
      double precision mass,p1Dp2
      integer l,r,swap(2),l1,l2,l3
      parameter(r=2,l=1)
      data swap/2,1/
      
      p1Dp2=0.5d0*za(1,2)*zb(2,1)
      mass=mt
      mqqb(r,r,r)=-zb(1,3)*zb(2,4)/P1DP2
      mqqb(l,l,l)=-za(2,4)*zb(1,3)*zb(1,4)/za(1,4)/P1DP2
      mqqb(r,l,l)=-zb(1,4)*zb(2,3)/P1DP2
      mqqb(l,r,r)=-za(2,3)*zb(1,3)*zb(1,4)/za(1,3)/P1DP2
      mqqb(r,r,l)=mass
     . *(zb(1,2)*zb(1,3)/za(1,4)+zb(1,2)*zb(1,4)/za(1,3))/P1DP2
      mqqb(l,l,r)=mass 
     . *(za(1,2)*zb(1,3)/za(1,4)+za(1,2)*zb(1,4)/za(1,3))/P1DP2
      mqqb(r,l,r) =czip
      mqqb(l,r,l) =czip

      do l1=1,2
      do l2=1,2
      do l3=1,2
      mqbq(l1,l2,l3)=mqqb(swap(l1),l2,l3)
      enddo
      enddo
      enddo

      return
      end 


      subroutine hhhamp(m1,m2)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex m1(2,2,2,2),m2(2,2,2,2)
      double precision XM,XM2,s12,s13,s35,s46,rt2o3t3,rt2o4t4
      integer g3,g4,t1,t2,r1,r2,j1,j2,l,r,loop,l3,l4
      parameter(l=1,r=2)
      XM=mt
      XM2=mt**2
      do loop=1,2
      if (loop.eq.1) then
      j1=2
      j2=1
      else
      j1=1
      j2=2
      endif

      s12=s(1,2)
      if (j1 .eq. 2) s13=s(1,4)
      if (j1 .eq. 1) s13=s(1,3)
      s35=s(3,5)
      s46=s(4,6)
      rt2o3t3=1d0/sqrt(abs(s35))
      rt2o4t4=1d0/sqrt(abs(s46))

      m1(r,r,r,r) =
     &  + s12**(-1)*rt2O3T3*rt2O4T4 * ( za(j1,6)*zb(j1,5)*zb(j1,j2)
     &    /za(j1,j2)*XM2 - za(j2,6)*zb(j1,j2)*zb(j2,5)/za(j1,j2)*XM2
     &     + za(4,6)*za(j1,3)*zb(3,5)*zb(j1,4)*zb(j1,j2)/za(j1,j2) - 
     &    za(4,6)*za(j2,3)*zb(3,5)*zb(j1,j2)*zb(j2,4)/za(j1,j2) )
     &  + s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,6)*
     &    za(j2,5)*zb(j1,5)*zb(j2,5)/za(j1,j2)**2*XM2**2 + 2*za(4,6)*
     &    za(j1,5)*za(j2,3)*zb(3,5)*zb(j1,5)*zb(j2,4)/za(j1,j2)**2*XM2
     &     )
     &  + s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,6)*za(j2,3)*zb(
     &    j1,5)*zb(j2,3)/za(j1,j2)**2*XM2 + 2*za(j1,6)*za(j2,3)*zb(3,5
     &    )*zb(j1,j2)/za(j1,j2)**2*XM2 + 2*za(j1,6)*zb(j1,5)*zb(j1,j2)
     &    /za(j1,j2)*XM2 + 2*za(4,6)*za(j1,3)*za(j2,3)*zb(3,5)*zb(j1,3
     &    )*zb(j2,4)/za(j1,j2)**2 - 2*za(4,6)*zb(j1,5)*zb(j2,4)/za(j1
     &    ,j2)*XM2 )

      m1(l,l,l,l) =
     &  + s12**(-1)*rt2O3T3*rt2O4T4 * (  - za(j1,j2)*za(j2,5)*zb(j2,6)
     &    /zb(j1,j2)*XM2 + za(j1,5)*za(j1,j2)*zb(j1,6)/zb(j1,j2)*XM2
     &     - za(3,5)*za(j1,j2)*za(j2,4)*zb(4,6)*zb(j2,3)/zb(j1,j2) + 
     &    za(3,5)*za(j1,4)*za(j1,j2)*zb(4,6)*zb(j1,3)/zb(j1,j2) )
     &  + s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)*
     &    za(j2,5)*zb(j1,6)*zb(j2,5)/zb(j1,j2)**2*XM2**2 + 2*za(3,5)*
     &    za(j1,5)*za(j2,4)*zb(4,6)*zb(j1,5)*zb(j2,3)/zb(j1,j2)**2*XM2
     &     )
     &  + s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)*za(j1,j2)*
     &    zb(j1,6)/zb(j1,j2)*XM2 + 2*za(j1,5)*za(j2,3)*zb(j1,6)*zb(j2,
     &    3)/zb(j1,j2)**2*XM2 - 2*za(j1,5)*za(j2,4)*zb(4,6)/zb(j1,j2)
     &    *XM2 + 2*za(3,5)*za(j1,j2)*zb(j1,6)*zb(j2,3)/zb(j1,j2)**2*
     &    XM2 + 2*za(3,5)*za(j1,3)*za(j2,4)*zb(4,6)*zb(j1,3)*zb(j2,3)/
     &    zb(j1,j2)**2 )

      m1(r,r,l,r) =
     &  + XM*s12**(-1)*rt2O3T3*rt2O4T4 * (  - za(3,5)*za(j1,6)*zb(j1,3)
     &    *zb(j1,j2)/za(j1,j2) + za(3,5)*za(j2,6)*zb(j1,j2)*zb(j2,3)
     &    /za(j1,j2) - za(4,6)*za(j1,5)*zb(j1,4)*zb(j1,j2)/za(j1,j2)
     &     + za(4,6)*za(j2,5)*zb(j1,j2)*zb(j2,4)/za(j1,j2) )
     &  + XM*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(3,
     &    5)*za(j1,6)*za(j2,5)*zb(j1,3)*zb(j2,5)/za(j1,j2)**2*XM2 - 2*
     &    za(4,6)*za(j1,5)*za(j2,5)*zb(j1,5)*zb(j2,4)/za(j1,j2)**2*XM2
     &     )
     &  + XM*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j1,6)*za(j2,
     &    5)*zb(j1,j2)/za(j1,j2)**2*XM2 - 2*za(3,5)*za(j1,6)*za(j2,3)*
     &    zb(j1,3)*zb(j2,3)/za(j1,j2)**2 - 2*za(3,5)*za(j1,6)*zb(j1,3)
     &    *zb(j1,j2)/za(j1,j2) + 2*za(3,5)*za(4,6)*zb(j1,3)*zb(j2,4)
     &    /za(j1,j2) - 2*za(4,6)*za(j1,3)*za(j2,5)*zb(j1,3)*zb(j2,4)
     &    /za(j1,j2)**2 )

      m1(l,l,r,l) =
     &  + XM*s12**(-1)*rt2O3T3*rt2O4T4 * ( za(j1,j2)*za(j2,3)*zb(3,5)*
     &    zb(j2,6)/zb(j1,j2) + za(j1,j2)*za(j2,4)*zb(4,6)*zb(j2,5)
     &    /zb(j1,j2) - za(j1,3)*za(j1,j2)*zb(3,5)*zb(j1,6)/zb(j1,j2)
     &     - za(j1,4)*za(j1,j2)*zb(4,6)*zb(j1,5)/zb(j1,j2) )
     &  + XM*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j1
     &    ,3)*za(j2,5)*zb(3,5)*zb(j1,6)*zb(j2,5)/zb(j1,j2)**2*XM2 - 2*
     &    za(j1,5)*za(j2,4)*zb(4,6)*zb(j1,5)*zb(j2,5)/zb(j1,j2)**2*XM2
     &     )
     &  + XM*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j1,j2)*zb(j1
     &    ,6)*zb(j2,5)/zb(j1,j2)**2*XM2 - 2*za(j1,3)*za(j1,j2)*zb(3,5)
     &    *zb(j1,6)/zb(j1,j2) - 2*za(j1,3)*za(j2,3)*zb(3,5)*zb(j1,6)*
     &    zb(j2,3)/zb(j1,j2)**2 + 2*za(j1,3)*za(j2,4)*zb(3,5)*zb(4,6)
     &    /zb(j1,j2) - 2*za(j1,3)*za(j2,4)*zb(4,6)*zb(j1,3)*zb(j2,5)
     &    /zb(j1,j2)**2 )

      m1(r,r,l,l) =
     &  + s12**(-1)*rt2O3T3*rt2O4T4 * ( za(j1,5)*zb(j1,6)*zb(j1,j2)
     &    /za(j1,j2)*XM2 - za(j2,5)*zb(j1,j2)*zb(j2,6)/za(j1,j2)*XM2
     &     + za(3,5)*za(j1,4)*zb(4,6)*zb(j1,3)*zb(j1,j2)/za(j1,j2) - 
     &    za(3,5)*za(j2,4)*zb(4,6)*zb(j1,j2)*zb(j2,3)/za(j1,j2) )
     &  + s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)*
     &    za(j2,5)*zb(j1,5)*zb(j2,6)/za(j1,j2)**2*XM2**2 + 2*za(3,5)*
     &    za(j1,4)*za(j2,5)*zb(4,6)*zb(j1,3)*zb(j2,5)/za(j1,j2)**2*XM2
     &     )
     &  + s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,3)*za(j2,5)*zb(
     &    j1,3)*zb(j2,6)/za(j1,j2)**2*XM2 + 2*za(j1,4)*za(j2,5)*zb(4,6
     &    )*zb(j1,j2)/za(j1,j2)**2*XM2 + 2*za(3,5)*za(j1,4)*za(j2,3)*
     &    zb(4,6)*zb(j1,3)*zb(j2,3)/za(j1,j2)**2 + 2*za(3,5)*za(j1,4)*
     &    zb(4,6)*zb(j1,3)*zb(j1,j2)/za(j1,j2) - 2*za(3,5)*zb(j1,3)*
     &    zb(j2,6)/za(j1,j2)*XM2 )

      m1(l,l,r,r) =
     &  + s12**(-1)*rt2O3T3*rt2O4T4 * (  - za(j1,j2)*za(j2,6)*zb(j2,5)
     &    /zb(j1,j2)*XM2 + za(j1,6)*za(j1,j2)*zb(j1,5)/zb(j1,j2)*XM2
     &     - za(4,6)*za(j1,j2)*za(j2,3)*zb(3,5)*zb(j2,4)/zb(j1,j2) + 
     &    za(4,6)*za(j1,3)*za(j1,j2)*zb(3,5)*zb(j1,4)/zb(j1,j2) )
     &  + s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)*
     &    za(j2,6)*zb(j1,5)*zb(j2,5)/zb(j1,j2)**2*XM2**2 + 2*za(4,6)*
     &    za(j1,3)*za(j2,5)*zb(3,5)*zb(j1,4)*zb(j2,5)/zb(j1,j2)**2*XM2
     &     )
     &  + s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,3)*za(j2,6)*zb(
     &    j1,3)*zb(j2,5)/zb(j1,j2)**2*XM2 - 2*za(j1,3)*za(j2,6)*zb(3,5
     &    )/zb(j1,j2)*XM2 + 2*za(4,6)*za(j1,j2)*zb(j1,4)*zb(j2,5)/zb(
     &    j1,j2)**2*XM2 + 2*za(4,6)*za(j1,3)*za(j1,j2)*zb(3,5)*zb(j1,4)
     &    /zb(j1,j2) + 2*za(4,6)*za(j1,3)*za(j2,3)*zb(3,5)*zb(j1,4)*
     &    zb(j2,3)/zb(j1,j2)**2 )

      m1(r,r,r,l) =
     &  + XM*s12**(-1)*rt2O3T3*rt2O4T4 * (  - za(j1,3)*zb(3,5)*zb(j1,6)
     &    *zb(j1,j2)/za(j1,j2) - za(j1,4)*zb(4,6)*zb(j1,5)*zb(j1,j2)
     &    /za(j1,j2) + za(j2,3)*zb(3,5)*zb(j1,j2)*zb(j2,6)/za(j1,j2)
     &     + za(j2,4)*zb(4,6)*zb(j1,j2)*zb(j2,5)/za(j1,j2) )
     & + XM*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j1
     &    ,4)*za(j2,5)*zb(4,6)*zb(j1,5)*zb(j2,5)/za(j1,j2)**2*XM2 - 2*
     &    za(j1,5)*za(j2,3)*zb(3,5)*zb(j1,5)*zb(j2,6)/za(j1,j2)**2*XM2
     &     )
     & + XM*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j1,3)*za(j2,
     &    3)*zb(3,5)*zb(j1,3)*zb(j2,6)/za(j1,j2)**2 - 2*za(j1,4)*za(j2
     &    ,3)*zb(3,5)*zb(4,6)*zb(j1,j2)/za(j1,j2)**2 - 2*za(j1,4)*za(
     &    j2,3)*zb(4,6)*zb(j1,5)*zb(j2,3)/za(j1,j2)**2 - 2*za(j1,4)*
     &    zb(4,6)*zb(j1,5)*zb(j1,j2)/za(j1,j2) + 2*zb(j1,5)*zb(j2,6)
     &    /za(j1,j2)*XM2 )

      m1(l,l,l,r) =
     &  + XM*s12**(-1)*rt2O3T3*rt2O4T4 * ( za(3,5)*za(j1,j2)*za(j2,6)*
     &    zb(j2,3)/zb(j1,j2) - za(3,5)*za(j1,6)*za(j1,j2)*zb(j1,3)
     &    /zb(j1,j2) + za(4,6)*za(j1,j2)*za(j2,5)*zb(j2,4)/zb(j1,j2)
     &     - za(4,6)*za(j1,5)*za(j1,j2)*zb(j1,4)/zb(j1,j2) )
     &  + XM*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(3,
     &    5)*za(j1,5)*za(j2,6)*zb(j1,5)*zb(j2,3)/zb(j1,j2)**2*XM2 - 2*
     &    za(4,6)*za(j1,5)*za(j2,5)*zb(j1,4)*zb(j2,5)/zb(j1,j2)**2*XM2
     &     )
     &  + XM*s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)*za(j2,6)
     &    /zb(j1,j2)*XM2 - 2*za(3,5)*za(j1,3)*za(j2,6)*zb(j1,3)*zb(j2,3
     &    )/zb(j1,j2)**2 - 2*za(3,5)*za(4,6)*za(j1,j2)*zb(j1,4)*zb(j2,
     &    3)/zb(j1,j2)**2 - 2*za(4,6)*za(j1,5)*za(j1,j2)*zb(j1,4)/zb(
     &    j1,j2) - 2*za(4,6)*za(j1,5)*za(j2,3)*zb(j1,4)*zb(j2,3)/zb(j1
     &    ,j2)**2 )

      m1(l,r,r,r) =
     &  + s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)*
     &    za(j1,6)*zb(j2,5)**2*XM2**2 + 2*za(4,6)*za(j1,3)*za(j1,5)*zb(
     &    3,5)*zb(j2,4)*zb(j2,5)*XM2 )
     &  + s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,3)*
     &    za(j1,6)*zb(j2,3)*zb(j2,5)*XM2 + 2*za(4,6)*za(j1,3)**2*zb(3,5
     &    )*zb(j2,3)*zb(j2,4) )

      m1(r,l,r,r) =
     &  + s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j2,5)*
     &    za(j2,6)*zb(j1,5)**2*XM2**2 + 2*za(4,6)*za(j2,3)*za(j2,5)*zb(
     &    3,5)*zb(j1,4)*zb(j1,5)*XM2 )
     &  + s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j2,3)*
     &    za(j2,6)*zb(j1,3)*zb(j1,5)*XM2 + 2*za(4,6)*za(j2,3)**2*zb(3,5
     &    )*zb(j1,3)*zb(j1,4) )

      m1(l,r,l,l) =
     &  + s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,5)
     &    **2*zb(j2,5)*zb(j2,6)*XM2**2 + 2*za(3,5)*za(j1,4)*za(j1,5)*
     &    zb(4,6)*zb(j2,3)*zb(j2,5)*XM2 )
     &  + s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j1,3)*
     &    za(j1,5)*zb(j2,3)*zb(j2,6)*XM2 + 2*za(3,5)*za(j1,3)*za(j1,4)*
     &    zb(4,6)*zb(j2,3)**2 )

      m1(r,l,l,l) =
     &  + s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j2,5)
     &    **2*zb(j1,5)*zb(j1,6)*XM2**2 + 2*za(3,5)*za(j2,4)*za(j2,5)*
     &    zb(4,6)*zb(j1,3)*zb(j1,5)*XM2 )
     &  + s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * ( 2*za(j2,3)*
     &    za(j2,5)*zb(j1,3)*zb(j1,6)*XM2 + 2*za(3,5)*za(j2,3)*za(j2,4)*
     &    zb(4,6)*zb(j1,3)**2 )

      m1(l,r,l,r) =
     &  + XM*s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(
     &    3,5)*za(j1,5)*za(j1,6)*zb(j2,3)*zb(j2,5)*XM2 - 2*za(4,6)*za(
     &    j1,5)**2*zb(j2,4)*zb(j2,5)*XM2 )
     &  + XM*s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(3,
     &    5)*za(j1,3)*za(j1,6)*zb(j2,3)**2 - 2*za(4,6)*za(j1,3)*za(j1,5
     &    )*zb(j2,3)*zb(j2,4) )

      m1(r,l,l,r) =
     &  + XM*s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(
     &    3,5)*za(j2,5)*za(j2,6)*zb(j1,3)*zb(j1,5)*XM2 - 2*za(4,6)*za(
     &    j2,5)**2*zb(j1,4)*zb(j1,5)*XM2 )
     &  + XM*s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(3,
     &    5)*za(j2,3)*za(j2,6)*zb(j1,3)**2 - 2*za(4,6)*za(j2,3)*za(j2,5
     &    )*zb(j1,3)*zb(j1,4) )

      m1(l,r,r,l) =
     &  + XM*s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(
     &    j1,3)*za(j1,5)*zb(3,5)*zb(j2,5)*zb(j2,6)*XM2 - 2*za(j1,4)*za(
     &    j1,5)*zb(4,6)*zb(j2,5)**2*XM2 )
     &  + XM*s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j1
     &    ,3)**2*zb(3,5)*zb(j2,3)*zb(j2,6) - 2*za(j1,3)*za(j1,4)*zb(4,6
     &    )*zb(j2,3)*zb(j2,5) )

      m1(r,l,r,l) =
     &  + XM*s12**(-1)*s13**(-1)*s35**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(
     &    j2,3)*za(j2,5)*zb(3,5)*zb(j1,5)*zb(j1,6)*XM2 - 2*za(j2,4)*za(
     &    j2,5)*zb(4,6)*zb(j1,5)**2*XM2 )
     &  + XM*s12**(-1)*s13**(-1)*rt2O3T3*rt2O4T4 * (  - 2*za(j2
     &    ,3)**2*zb(3,5)*zb(j1,3)*zb(j1,6) - 2*za(j2,3)*za(j2,4)*zb(4,6
     &    )*zb(j1,3)*zb(j1,5) )


      if (loop.eq.1) then
         do l3=1,2
         do l4=1,2
             m2(1,1,l3,l4)=m1(1,1,l3,l4)
             m2(2,2,l3,l4)=m1(2,2,l3,l4)
             m2(1,2,l3,l4)=m1(2,1,l3,l4)
             m2(2,1,l3,l4)=m1(1,2,l3,l4)
         enddo
         enddo
      endif

      enddo

      return

      end
 
