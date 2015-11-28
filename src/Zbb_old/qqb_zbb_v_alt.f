      subroutine qqb_zbb_v(P,msqv)
      implicit none
************************************************************************
*     Author: J.M. Campbell                                            *
*     March, 1999.                                                     *
*     calculate the virtual matrix element squared and subtraction terms
*     for the process                                                  *
*     q(-p1)+qb(-p2) --> b(p4)+bb(p5)+e^-(p6)+e^+(p7)                  *
*     q(-p1) +Q(-p5)+ l(-p7) -->   q(p2)+Q(p4) +l(p6)                  *
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'epinv.f'
      include 'scale.f'
      include 'hardscale.f'
      logical msbar
      common/msbar/msbar
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . p(mxpart,4),q_bdkw(mxpart,4),faclo,dot,
     . fac,xl12,xl15,xl16,xl25,xl26,xl56,subqbq,subqqb,v2(2),vQ(nf,2)
      double complex tamp,lamp,atreez,a61z
      integer nu,j,k,polq,polb,polz

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

c---calculate the lowest order matrix element and fill the common block
c---twopij with s_{ij} (in rke notation)
      call qqb_zbb(p,msq)

      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return 

c---add result of integrating subtraction terms
      xl12=log(s(1,2)/musq)
      xl15=log(-s(1,5)/musq)
      xl16=log(-s(1,6)/musq)
      xl25=log(-s(2,5)/musq)
      xl26=log(-s(2,6)/musq)
      xl56=log(s(5,6)/musq)

      subqqb=
     & xn*(two*epinv**2+epinv*(three-xl15-xl26)
     & +2.5d0-two*pisqo6-0.75d0*(xl15+xl26)
     & +half*(xl15**2+xl26**2))
     & +(-two*epinv**2
     & +epinv*(two*(xl15-xl16+xl26-xl25)+xl12+xl56-three)
     & +4*(pisqo6-one)+1.5d0*(xl15-xl16+xl26-xl25+xl56)
     & -half*(xl12**2+xl56**2)
     & +xl16**2-xl15**2+xl25**2-xl26**2)/xn      
      if (msbar) subqqb=subqqb+2*cf
      subqqb=subqqb*ason2pi

      subqbq=
     & xn*(two*epinv**2+epinv*(three-xl25-xl16)
     & +2.5d0-two*pisqo6-0.75d0*(xl25+xl16)
     & +half*(xl25**2+xl16**2))
     & +(-two*epinv**2
     & +epinv*(two*(xl25-xl26+xl16-xl15)+xl12+xl56-three)
     & +4*(pisqo6-one)+1.5d0*(xl25-xl26+xl16-xl15+xl56)
     & -half*(xl12**2+xl56**2)
     & +xl26**2-xl25**2+xl15**2-xl16**2)/xn      
      if (msbar) subqbq=subqbq+2*cf
      subqbq=subqbq*ason2pi



c---  Now transform momenta into a notation 
c---  suitable for calling the BDKW function with notation which is 
c---  q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
      do nu=1,4
      q_bdkw(1,nu)=p(2,nu)
      q_bdkw(2,nu)=p(6,nu)
      q_bdkw(3,nu)=p(5,nu)
      q_bdkw(4,nu)=p(1,nu)
      q_bdkw(5,nu)=p(4,nu)
      q_bdkw(6,nu)=p(3,nu)
      enddo      
      call spinoru(6,q_bdkw,za,zb)

      faclo=4d0*V*aveqq*esq**2*gsq**2
      fac=faclo*xn*0.5d0*ason2pi

      v2(1)=l1
      v2(2)=r1

      do j=1,nf
        vQ(j,1)=L(j)
        vQ(j,2)=R(j)
      enddo

      do j=-nf,nf
      k=-j

      do polq=1,2
      do polb=1,2
      do polz=1,2
        if     ((j .eq. 0) .and. (k .eq. 0)) then
          tamp=0d0
          lamp=0d0
        elseif ((j .gt. 0) .and. (k .lt. 0)) then
          tamp=atreez(polq,polb,polz,1,2,3,4,5,6,za,zb)
     .         *vQ(j,polq)*v2(polz)
          lamp=a61z(polq,polb,polz,1,2,3,4,5,6,za,zb)
     .         *vQ(j,polq)*v2(polz)
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
          tamp=atreez(3-polq,polb,polz,1,2,3,4,5,6,za,zb)
     .         *vQ(k,polq)*v2(polz)
          lamp=a61z(3-polq,polb,polz,1,2,3,4,5,6,za,zb)
     .         *vQ(k,polq)*v2(polz)
        endif
c        write(*,*) tamp,lamp,tamp*dconjg(lamp)+lamp*dconjg(tamp)
        msqv(j,k)=msqv(j,k)+fac*2d0*dble(tamp*dconjg(lamp))
      enddo
      enddo
      enddo

      if    ((j .gt. 0) .and. (k .lt. 0)) then
        msqv(j,k)=msqv(j,k)+subqqb*msq(j,k)          
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
        msqv(j,k)=msqv(j,k)+subqbq*msq(j,k)    
      endif
c      pause
      
      enddo

      return
      end
     
