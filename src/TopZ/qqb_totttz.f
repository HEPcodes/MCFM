      subroutine qqb_totttz(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared 
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=t(p3)+t(p4)+z(p5)
C  
************************************************************************
      include 'constants.f'
      
      integer j,k,nu,np
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision ddb,dbd,uub,ubu,
     . SDDB_TTBZ,SUUB_TTBZ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      
      do nu=1,4
      np=nu
      if (nu .eq. 4) np=0
      p1(np)=-p(1,nu)
      p2(np)=-p(2,nu)
      p3(np)=+p(3,nu)
      p4(np)=+p(4,nu)
      p5(np)=+p(5,nu)
      enddo

      ddb=SDDB_TTBZ(P1, P2, P3, P4, P5)
      uub=SuuB_TTBZ(P1, P2, P3, P4, P5)
      dbd=SDDB_TTBZ(P2, P1, P3, P4, P5)
      ubu=SuuB_TTBZ(P2, P1, P3, P4, P5)
C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j .lt. 0) then
          msq(j,-j)=0d0
      elseif ((j .eq. +1) .or. (j .eq. +3)  .or. (j .eq. +5)) then
          msq(j,j)=ddb
      elseif ((j .eq. -1) .or. (j .eq. -3)  .or. (j .eq. -5)) then
          msq(j,j)=dbd
      elseif ((j .eq. +2) .or. (j .eq. +4)) then
          msq(j,j)=uub
      elseif ((j .eq. -2) .or. (j .eq. -4)) then
          msq(j,j)=ubu
      endif
      enddo
      return
      end
 
