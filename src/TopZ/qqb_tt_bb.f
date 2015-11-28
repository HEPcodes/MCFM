      subroutine qqb_tt_bb(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
C                       +t~(b~(p6)+e^-(p7)+nu(p8))+b(p9)+bbar(p10)     *
C                                                                      * 
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprodx.f'
      integer nu,np,j,k,h1,h2
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),ampsq(2,2)
      double precision wtgg,wtqqb,wtqbq
      double complex T1(2,2),T2(2,2),T3(2,2),T4(2,2),T5(2,2),
     . T6(2,2),T7(2,2),x1(2,2),x2(2,2),x3(2,2),x4(2,2),x12(2,2)
      call spinoru(10,p,za,zb)
      call diag1(1,2,3,4,5,6,7,8,9,10,one,T1)      
      call diag2(1,2,3,4,5,6,7,8,9,10,one,T2)      
      call diag3(1,2,3,4,5,6,7,8,9,10,one,T3)      
      call diag4(1,2,3,4,5,6,7,8,9,10,one,T4)      
      call diag5(1,2,3,4,5,6,7,8,9,10,one,T5)      
      call diag6(1,2,3,4,5,6,7,8,9,10,one,T6)      
      call diag7(1,2,3,4,5,6,7,8,9,10,one,T7)      

      wtqqb=0d0
      do h1=1,2
      do h2=1,2
      x1(h1,h2)=0.25d0*(t1(h1,h2)+t2(h1,h2))
      x2(h1,h2)=0.25d0*(t3(h1,h2)+t4(h1,h2)+t5(h1,h2)+t6(h1,h2))
      x12(h1,h2)=x1(h1,h2)+x2(h1,h2)
      x3(h1,h2)=0.25d0*(t1(h1,h2)+t3(h1,h2)+t5(h1,h2)-t7(h1,h2))
      x4(h1,h2)=0.25d0*(t2(h1,h2)+t4(h1,h2)+t6(h1,h2)+t7(h1,h2))
      ampsq(h1,h2)=
     . +2d0*V/xn*dble(x1(h1,h2)*Dconjg(x2(h1,h2)))
     . +2d0*V*xn*dble(x3(h1,h2)*Dconjg(x4(h1,h2)))
     . +(xn**4+6d0*xn**2-3d0)/xn*abs(x12(h1,h2))**2
      wtqqb=wtqqb+ampsq(h1,h2)
      enddo
      enddo

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo


C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=wtqqb
      elseif (j .eq. 0) then
          msq(j,j)=wtgg
      elseif (j .gt. 0) then
          msq(j,-j)=wtqbq
      endif
      enddo
      return
      end

