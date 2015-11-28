      subroutine qqb_wbbm_g(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 2003.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+bb(p6)
c   with mass for the b and the bbar
c   positively charged W only
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j,k,nu
      double precision p(mxpart,4),q(mxpart,4),dot,
     . msq(-nf:nf,-nf:nf),redmsqm,fac,al(5:6)
      double precision qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb
      double precision mbbwbbj
      common/wbbj/mbbwbbj

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      fac=gsq**3*gw**4/4d0

c      write(6,*) 'mb',mb

c--- q-qb and qb-q
      al(5)=mb**2/(2d0*dot(p,5,7))
      al(6)=mb**2/(2d0*dot(p,6,7))
      do nu=1,4
      do j=1,7
         if ((j.eq.5).or.(j.eq.6)) then
             q(j,nu)=p(j,nu)-al(j)*p(7,nu)
         else
             q(j,nu)=p(j,nu)
         endif
      enddo
      enddo
      call spinoru(7,q,za,zb)

      mbbwbbj=dsqrt(s(5,6))

      qqbWbbg =+redmsqm(1,2,7,5,6,3,4,mb)*fac*aveqq
      qbqWbbg =+redmsqm(2,1,7,5,6,3,4,mb)*fac*aveqq

c--- q-g and qb-g
      al(5)=mb**2/(2d0*dot(p,5,2))
      al(6)=mb**2/(2d0*dot(p,6,2))
      do nu=1,4
      do j=1,7
         if ((j.eq.5).or.(j.eq.6)) then
             q(j,nu)=p(j,nu)-al(j)*p(2,nu)
         else
             q(j,nu)=p(j,nu)
         endif
      enddo
      enddo
      call spinoru(7,q,za,zb)

      qgWbbq  =+redmsqm(1,7,2,5,6,3,4,mb)*fac*aveqg
      qbgWbbqb=+redmsqm(7,1,2,5,6,3,4,mb)*fac*aveqg

c--- g-q and g-qb
      al(5)=mb**2/(2d0*dot(p,5,1))
      al(6)=mb**2/(2d0*dot(p,6,1))
      do nu=1,4
      do j=1,7
         if ((j.eq.5).or.(j.eq.6)) then
             q(j,nu)=p(j,nu)-al(j)*p(1,nu)
         else
             q(j,nu)=p(j,nu)
         endif
      enddo
      enddo
      call spinoru(7,q,za,zb)

      gqWbbq  =+redmsqm(2,7,1,5,6,3,4,mb)*fac*aveqg
      gqbWbbqb=+redmsqm(7,2,1,5,6,3,4,mb)*fac*aveqg

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)

      msq(j,k)=0d0

      if     ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsq(j,k)*qqbWbbg

      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsq(j,k)*qbqWbbg

      elseif ((j .gt. 0) .and. (k .eq. 0)) then
      msq(j,k)=Vsum(j)*qgWbbq 
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
      msq(j,k)=Vsum(j)*qbgWbbqb

      elseif ((j .eq. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsum(k)*gqWbbq

      elseif ((j .eq. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsum(k)*gqbWbbqb
      endif

      enddo
      enddo

      return
      end


      double precision function redmsqm(j1,j2,j3,j4,j5,j6,j7,bmass)
      implicit none
c matrix element squared summed over colors and spins
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jh,sb,sc
      double complex qcda(2,2,2),qcdb(2,2,2),qedi(2,2,2),qedf(2,2,2)
      double precision prop,bmass
c---calculate the W propagator
      prop=((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2)
            
      call Wbbgmassa(j1,j2,j3,j4,j5,j6,j7,bmass,qcda)
      call Wbbgmassb(j1,j2,j3,j4,j5,j6,j7,bmass,qcdb)
      call Wbbgmassi(j1,j2,j3,j4,j5,j6,j7,bmass,qedi)
      call Wbbgmassf(j1,j2,j3,j4,j5,j6,j7,bmass,qedf)

      redmsqm=zip
      do jh=1,2
      do sb=1,2
      do sc=1,2

      redmsqm=redmsqm+
     & V*xn/eight*(cdabs(qcda(jh,sb,sc))**2+cdabs(qcdb(jh,sb,sc))**2)
     &+V/(eight*xn)*(cdabs(qedi(jh,sb,sc))**2+cdabs(qedf(jh,sb,sc))**2
     &-two*(cdabs(qedi(jh,sb,sc)+qedf(jh,sb,sc)))**2)
      enddo
      enddo
      enddo

      redmsqm=redmsqm/prop
      return
      end
