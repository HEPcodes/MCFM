      subroutine qqb_wbb_v(P,msqv)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     calculate the virtual matrix element squared and subtraction terms
*     for the process                                                  *
*     q(-p1) +Q(-p6)+ l(-p4) -->   q(p2)+Q(p5) +l(p3)                  *
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'sprodx.f'
      include 'epinv.f'
      include 'scale.f'
      include 'hardscale.f'
      logical msbar
      common/msbar/msbar
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . p(mxpart,4),q(mxpart,4),srke(mxpart,mxpart),faclo,fac,
     . qqb,qbq,
     . xl12,xl15,xl16,xl25,xl26,xl56,subqbq,subqqb
      double complex atrLLL,atrLRL,a61LLL,a61LRL
      double complex tLLL,tLRL,fLLL,fLRL
      integer nu,j,k
      double precision ii_qg,ii_gq,if_qg,fi_qg,ff_qg
      common/twopij/srke

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

c---calculate the lowest order matrix element and fill the common block
c---twopij with srke_{ij}
      call qqb_wbb(p,msq)
      if (
     .      (srke(5,6) .lt. four*hscalesq) 
     . .or. (srke(1,5)*srke(2,5)/srke(1,2) .lt. hscalesq) 
     . .or. (srke(1,6)*srke(2,6)/srke(1,2) .lt. hscalesq) ) return 

c---add result of integrating subtraction terms
      xl12=log(srke(1,2)/musq)
      xl15=log(-srke(1,5)/musq)
      xl16=log(-srke(1,6)/musq)
      xl25=log(-srke(2,5)/musq)
      xl26=log(-srke(2,6)/musq)
      xl56=log(srke(5,6)/musq)

c      subqqb=
c     & xn*(two*epinv**2+epinv*(three-xl15-xl26)
c     & +2.5d0-two*pisqo6-0.75d0*(xl15+xl26)
c     & +half*(xl15**2+xl26**2))
c     & +(-two*epinv**2
c     & +epinv*(two*(xl15-xl16+xl26-xl25)+xl12+xl56-three)
c     & +4*(pisqo6-one)+1.5d0*(xl15-xl16+xl26-xl25+xl56)
c     & -half*(xl12**2+xl56**2)
c     & +xl16**2-xl15**2+xl25**2-xl26**2)/xn      
c      if (msbar) subqqb=subqqb+2*cf
c      subqqb=subqqb*ason2pi

c      subqbq=
c     & xn*(two*epinv**2+epinv*(three-xl25-xl16)
c     & +2.5d0-two*pisqo6-0.75d0*(xl25+xl16)
c     & +half*(xl25**2+xl16**2))
c     & +(-two*epinv**2
c     & +epinv*(two*(xl25-xl26+xl16-xl15)+xl12+xl56-three)
c     & +4*(pisqo6-one)+1.5d0*(xl25-xl26+xl16-xl15+xl56)
c     & -half*(xl12**2+xl56**2)
c     & +xl26**2-xl25**2+xl15**2-xl16**2)/xn      
c      if (msbar) subqbq=subqbq+2*cf
c      subqbq=subqbq*ason2pi

c      write(*,*) 'subqqb (old) ',subqqb
c      write(*,*) 'subqbq (old) ',subqbq
      
      subqqb=0.5d0*((xn-two/xn)*(if_qg(one,xl15,1)+fi_qg(one,xl15,1))
     .                 +two/xn *(if_qg(one,xl16,1)+fi_qg(one,xl16,1))
     .                 -one/xn *(ii_qg(one,xl12,1)+ff_qg(one,xl56,1)))
     .      +0.5d0*((xn-two/xn)*(if_qg(one,xl26,1)+fi_qg(one,xl26,1))
     .                 +two/xn *(if_qg(one,xl25,1)+fi_qg(one,xl25,1))
     .                 -one/xn *(ii_qg(one,xl12,1)+ff_qg(one,xl56,1)))
     
      subqbq=0.5d0*((xn-two/xn)*(if_qg(one,xl25,1)+fi_qg(one,xl25,1))
     .                 +two/xn *(if_qg(one,xl26,1)+fi_qg(one,xl26,1))
     .                 -one/xn *(ii_qg(one,xl12,1)+ff_qg(one,xl56,1)))
     .      +0.5d0*((xn-two/xn)*(if_qg(one,xl16,1)+fi_qg(one,xl16,1))
     .                 +two/xn *(if_qg(one,xl15,1)+fi_qg(one,xl15,1))
     .                 -one/xn *(ii_qg(one,xl12,1)+ff_qg(one,xl56,1)))

      if (msbar) then
        subqbq=subqbq+2*cf
        subqqb=subqqb+2*cf
      endif

      subqqb=subqqb*ason2pi
      subqbq=subqbq*ason2pi

c      write(*,*) 'subqqb (new) ',subqqb
c      write(*,*) 'subqbq (new) ',subqbq
c      pause
      
c---  Now transform momenta into a notation 
c---  suitable for calling the BDKW function with notation which is 
c---  q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
      do nu=1,4
      q(1,nu)=p(2,nu)
      q(2,nu)=p(6,nu)
      q(3,nu)=p(5,nu)
      q(4,nu)=p(1,nu)
      q(5,nu)=p(4,nu)
      q(6,nu)=p(3,nu)
      enddo      
      call spinoru(6,q,za,zb)
      faclo=V*aveqq*gw**4*gsq**2
      fac=faclo*xn*0.5d0*ason2pi

c----do whatever needs to be done q-qb case
      tLLL=atrLLL(1,2,3,4,5,6,za,zb)
      fLLL=a61LLL(1,2,3,4,5,6,za,zb)
      tLRL=atrLRL(1,2,3,4,5,6,za,zb)
      fLRL=a61LRL(1,2,3,4,5,6,za,zb)

      qqb=fac*dble(tLLL*dconjg(fLLL)+fLLL*dconjg(tLLL)
     .            +tLRL*dconjg(fLRL)+fLRL*dconjg(tLRL))
      
c----now look at qb-q case, swap the momenta
      tLLL=atrLLL(4,2,3,1,5,6,za,zb)
      fLLL=a61LLL(4,2,3,1,5,6,za,zb)
      tLRL=atrLRL(4,2,3,1,5,6,za,zb)
      fLRL=a61LRL(4,2,3,1,5,6,za,zb)

      qbq=fac*dble(tLLL*dconjg(fLLL)+fLLL*dconjg(tLLL)
     .            +tLRL*dconjg(fLRL)+fLRL*dconjg(tLRL))

      do j=-nf,nf
      do k=-nf,nf
      if (Vsq(j,k) .eq. 0d0) goto 20
      if     ((j .gt. 0) .and. (k .lt. 0)) then
               msqv(j,k)=subqqb*msq(j,k)+Vsq(j,k)*qqb
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
               msqv(j,k)=subqbq*msq(j,k)+Vsq(j,k)*qbq
      else
      msqv(j,k)=0d0
      endif
   20 enddo
      enddo

      return
      end
     
