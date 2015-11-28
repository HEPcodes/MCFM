      subroutine qqb_ttbgg(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=t(p3)+tbar(p4)+g(p5)+g(p6)
C  
************************************************************************
      include 'constants.f'
      integer nu,np,j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision pu(0:3),pub(0:3),pt(0:3),ptb(0:3),pb(0:3),
     . pbb(0:3),wtqqb,wtqbq,wtgg,SUUB_TTBgg,SGG_TTBgg
      logical first
      data first/.true./

      if (first) then
      first=.false.
      call initialize
      endif
      
      do nu=1,4
      np=nu
      if (nu.eq.4) np=0
      pu(np)=-p(1,nu)
      pub(np)=-p(2,nu)
      pt(np)= p(3,nu)
      ptb(np)=p(4,nu)
      pb(np)= p(5,nu)
      pbb(np)=p(6,nu)
      enddo


      wtqqb=SUUB_TTBgg(pu, pub, pt, ptb, pb, pbb)
      wtqbq=SUUB_TTBgg(pub, pu, pt, ptb, pb, pbb)
      wtgg=SGG_TTBgg(pub, pu, pt, ptb, pb, pbb)

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

