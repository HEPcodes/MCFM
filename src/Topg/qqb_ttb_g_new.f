      subroutine qqb_ttb_g_new(p,msq)
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
*     Only seven diagrams including leading to 2 on-shell top quarks   *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'prods.f'
      integer b,j,k,h1,h2,nu
      double precision t(4),r(4),
     . msq(-nf:nf,-nf:nf),p(mxpart,4),ampsq(2,2),ps(mxpart,4)
      double precision wtgg,wtqqb,wtqbq,wtqg,fac
c      double complex T1(2,2),T2(2,2),T3(2,2),T4(2,2),T5(2,2),
c     . T6(2,2),T7(2,2),x1(2,2),x2(2,2),x3(2,2),x4(2,2),x12(2,2)
      double complex a,
     . ttbgggppp,ttbgggmpp,ttbgggpmp,ttbgggppm,
     . ttbgggmmm,ttbgggpmm,ttbgggmpm,ttbgggmmp
      double precision p3Dp5,p6Dp8,rDp7,tDp4
c      integer q,a,t1,t2,r1,r2,ep,em,g
c      parameter(q=1,a=2,t1=3,t2=5,r1=6,r2=8,ep=4,em=7,g=9)

      do nu=1,4
      do j=1,mxpart
      ps(j,nu)=p(j,nu)
      enddo
      enddo
      
      p3Dp5=p(3,4)*p(5,4)-p(3,3)*p(5,3)-p(3,2)*p(5,2)-p(3,1)*p(5,1) 
      p6Dp8=p(6,4)*p(8,4)-p(6,3)*p(8,3)-p(6,2)*p(8,2)-p(6,1)*p(8,1) 

c      we will have no further need for p3 and p5 
c      we will have no further need for p6 and p8 
      
      do nu=1,4
      t(nu)=p(3,nu)+p(4,nu)+p(5,nu)
      r(nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo      
      tDp4=t(4)*p(4,4)-t(3)*p(4,3)-t(2)*p(4,2)-t(1)*p(4,1) 
      rDp7=r(4)*p(7,4)-r(3)*p(7,3)-r(2)*p(7,2)-r(1)*p(7,1)             
      do nu=1,4
c---t2
      ps(5,nu)=0.5d0*mt**2/tDp4*p(4,nu)
c---t1
      ps(3,nu)=t(nu)-ps(5,nu)

c---t2
      ps(8,nu)=0.5d0*mt**2/rDp7*p(7,nu)
c---t1
      ps(6,nu)=r(nu)-ps(8,nu)

      enddo
c      call writeout(ps)
c      pause
      call spinoru(9,ps,za,zb)

      write(6,*) 'ppp'
      do j=3,8
      a=ttbgggppp(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo
      write(6,*) 'mmm'
      do j=3,8
      a=ttbgggmmm(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo
      write(6,*) 'mpp'
      do j=3,8
      a=ttbgggmpp(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo
      write(6,*) 'pmm'
      do j=3,8
      a=ttbgggpmm(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo

      write(6,*) 'mpm'
      do j=3,8
      a=ttbgggmpm(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo

      write(6,*) 'mmp'
      do j=3,8
      a=ttbgggmmp(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo
       
      write(6,*) 'pmp'
      do j=3,8
      a=ttbgggpmp(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo

      write(6,*) 'ppm'
      do j=3,8
      a=ttbgggppm(1,2,3,4,5,6,7,8,9,j)
      write(6,*) j,a,cdabs(a)
      enddo
      pause
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

