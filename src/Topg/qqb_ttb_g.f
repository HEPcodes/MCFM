      subroutine qqb_ttb_g(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=t(nu(p3)+e^+(p4)+b(p5))                        *
C                       +t~(b~(p6)+e^-(p7)+nu(p8))+g(p9)               *
C                                                                      * 
C     Only seven diagrams including leading to 2 on-shell top quarks   *
C***********************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'prods.f'
      integer b,j,k,h1,h2
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),ampsq(2,2)
      double precision wtgg,wtqqb,wtqbq,wtqg,fac
      double complex T1(2,2),T2(2,2),T3(2,2),T4(2,2),T5(2,2),
     . T6(2,2),T7(2,2),x1(2,2),x2(2,2),x3(2,2),x4(2,2),x12(2,2)

      
      call spinoru(9,p,za,zb)

      b=1
      call diagg1(1,2,9,b,T1)      
      call diagg2(1,2,9,b,T2)      
      call diagg3(1,2,9,b,T3)      
      call diagg4(1,2,9,b,T4)      
      call diagg5(1,2,9,b,T5)      
      call diagg6(1,2,9,b,T6)      
      call diagg7(1,2,9,b,T7)      
      wtqqb=0d0
     
      do h1=1,2
      do h2=1,2
      x1(h1,h2)=-0.5d0*(t1(h1,h2)+t2(h1,h2))/xn
      x2(h1,h2)=-0.5d0*(t3(h1,h2)+t4(h1,h2)+t5(h1,h2)+t6(h1,h2))/xn
      x12(h1,h2)=x1(h1,h2)+x2(h1,h2)
      x3(h1,h2)=0.5d0*(t1(h1,h2)+t3(h1,h2)+t5(h1,h2)-t7(h1,h2))
      x4(h1,h2)=0.5d0*(t2(h1,h2)+t4(h1,h2)+t6(h1,h2)+t7(h1,h2))
      ampsq(h1,h2)=0.5d0*V/xn*(abs(x1(h1,h2))**2+abs(x2(h1,h2))**2
     .  +abs(x3(h1,h2))**2+abs(x4(h1,h2))**2-2d0*abs(x12(h1,h2))**2)
      wtqqb=wtqqb+ampsq(h1,h2)
      enddo
      enddo

      wtqbq=0d0


      b=1
      call diagg1(2,1,9,b,T1)      
      call diagg2(2,1,9,b,T2)      
      call diagg3(2,1,9,b,T3)      
      call diagg4(2,1,9,b,T4)      
      call diagg5(2,1,9,b,T5)      
      call diagg6(2,1,9,b,T6)      
      call diagg7(2,1,9,b,T7)      


      wtqqb=0d0


      do h1=1,2
      do h2=1,2
      x1(h1,h2)=-0.5d0*(t1(h1,h2)+t2(h1,h2))/xn
      x2(h1,h2)=-0.5d0*(t3(h1,h2)+t4(h1,h2)+t5(h1,h2)+t6(h1,h2))/xn
      x12(h1,h2)=x1(h1,h2)+x2(h1,h2)
      x3(h1,h2)=0.5d0*(t1(h1,h2)+t3(h1,h2)+t5(h1,h2)-t7(h1,h2))
      x4(h1,h2)=0.5d0*(t2(h1,h2)+t4(h1,h2)+t6(h1,h2)+t7(h1,h2))
      ampsq(h1,h2)=0.5d0*V/xn*(abs(x1(h1,h2))**2+abs(x2(h1,h2))**2
     .  +abs(x3(h1,h2))**2+abs(x4(h1,h2))**2-2d0*abs(x12(h1,h2))**2)
      wtqbq=wtqbq+ampsq(h1,h2)
      enddo
      enddo

      wtqg=0d0

      call diagg1(1,9,2,b,T1)      
      call diagg2(1,9,2,b,T2)      
      call diagg3(1,9,2,b,T3)      
      call diagg4(1,9,2,b,T4)      
      call diagg5(1,9,2,b,T5)      
      call diagg6(1,9,2,b,T6)      
      call diagg7(1,9,2,b,T7)      


      wtqg=0d0

      do h1=1,2
      do h2=1,2
      x1(h1,h2)=-0.5d0*(t1(h1,h2)+t2(h1,h2))/xn
      x2(h1,h2)=-0.5d0*(t3(h1,h2)+t4(h1,h2)+t5(h1,h2)+t6(h1,h2))/xn
      x12(h1,h2)=x1(h1,h2)+x2(h1,h2)
      x3(h1,h2)=0.5d0*(t1(h1,h2)+t3(h1,h2)+t5(h1,h2)-t7(h1,h2))
      x4(h1,h2)=0.5d0*(t2(h1,h2)+t4(h1,h2)+t6(h1,h2)+t7(h1,h2))
      ampsq(h1,h2)=0.5d0*V/xn*(abs(x1(h1,h2))**2+abs(x2(h1,h2))**2
     .  +abs(x3(h1,h2))**2+abs(x4(h1,h2))**2-2d0*abs(x12(h1,h2))**2)
      wtqg=wtqg+ampsq(h1,h2)
      enddo
      enddo



C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      wtgg=0d0

      fac=(gwsq/2d0)**4*gsq**3
C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=aveqq*fac*wtqqb
      elseif (j .eq. 0) then
          msq(j,j)=avegg*fac*wtgg
      elseif (j .gt. 0) then
          msq(j,-j)=aveqq*fac*wtqbq
      endif
      enddo
      return
      end

