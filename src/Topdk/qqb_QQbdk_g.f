      subroutine qqb_QQbdk_g(p,msq)
      implicit none
C***********************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                     *
*     calculate the element squared                                    *
*     for the process                                                  *
*----My notation                                                       *
*     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
*                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
*                                                                      * 
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer b,j,k,h1,h2,nu,jj,kk
      double precision t(4),r(4),
     . msq(-nf:nf,-nf:nf),p(mxpart,4),ampsq(2,2),ps(mxpart,4)
      double precision wt(-1:1,-1:1),fac,al4,al7,propw,
     . s34,s35,s68,s78,rDp7,tDp4,gggampsq,qqgampsq

      do nu=1,4
      do j=1,mxpart
      ps(j,nu)=zip
      enddo
      do j=1,2
      ps(j,nu)=p(j,nu)
      enddo
      enddo
      
      s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)-p(3,2)*p(4,2)-p(3,1)*p(4,1)) 
      s78=2d0*(p(7,4)*p(8,4)-p(7,3)*p(8,3)-p(7,2)*p(8,2)-p(7,1)*p(8,1)) 
      s35=2d0*(p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1))
      s68=2d0*(p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1)) 
      propw=((s34-wmass**2)**2+(wmass*wwidth)**2)
     .     *((s78-wmass**2)**2+(wmass*wwidth)**2)

      do nu=1,4
      t(nu)=p(3,nu)+p(4,nu)+p(5,nu)
      r(nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo      
c      we will have no further need for p3 and p5 
c      we will have no further need for p6 and p8 
      
      tDp4=t(4)*p(4,4)-t(3)*p(4,3)-t(2)*p(4,2)-t(1)*p(4,1) 
      rDp7=r(4)*p(7,4)-r(3)*p(7,3)-r(2)*p(7,2)-r(1)*p(7,1)             
      al4=0.5d0*mt**2/tDp4
      al7=0.5d0*mt**2/rDp7

      do nu=1,4
C  Fold 9 into 3
      ps(3,nu)=p(9,nu)
C  4 is the positron momentum which is rescaled
      ps(4,nu)=al4*p(4,nu)
C  5 is the massless residue of the top momentum
      ps(5,nu)=t(nu)-ps(4,nu)

      
C  7 is the electron momentum which is rescaled
      ps(7,nu)=al7*p(7,nu)
C  6 is the massless residue of the top-bar momentum
      ps(6,nu)=r(nu)-ps(7,nu)

      enddo

      call spinoru(7,ps,za,zb)

C---set all matrix elements to zero.

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      fac=gwsq**4*gsq**3/(al7*al4*propw*(mt*twidth)**4)*s35*s68

      wt(+1,-1)=fac*aveqq*qqgampsq(1,2,3)
      wt(-1,+1)=fac*aveqq*qqgampsq(2,1,3)

      wt(+1, 0)=fac*aveqg*qqgampsq(1,3,2)
      wt(-1, 0)=fac*aveqg*qqgampsq(3,1,2)

      wt( 0,+1)=fac*aveqg*qqgampsq(2,3,1)
      wt( 0,-1)=fac*aveqg*qqgampsq(3,2,1)

      wt( 0, 0)=fac*avegg*gggampsq()

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if (j .gt. 0) then
          if (k.eq.-j) then
          msq(j,k)=wt(1,-1)
          elseif (k .eq .0) then
          msq(j,k)=wt(1,0)
          endif
      elseif (j .lt. 0) then
          if (k.eq.-j) then
          msq(j,k)=wt(-1,1)
          elseif (k .eq .0) then
          msq(j,k)=wt(-1,0)
          endif
      elseif (j .eq. 0) then
          if (k.gt.0) then
          msq(j,k)=wt(0,1)
         elseif (k.eq.0) then
          msq(j,k)=wt(0,0)
          elseif (k.lt.0) then
          msq(j,k)=wt(0,-1)
          endif
      endif
      enddo
      enddo
      return
      end

