      subroutine qqb_zbb_v(P,msqv)
      implicit none
************************************************************************
*     Author: J.M. Campbell                                            *
*     March, 1999.                                                     *
*     Updated February, 2000                                           *
*     calculate the virtual matrix element squared and subtraction     *
*     terms for the process                                            *
*     q(-p1)+qb(-p2) --> e^-(p3)+e^+(p4)+b(p5)+bb(p6)                  *
*     q(-p1) +Q(-p5)+ l(-p7) -->   q(p2)+Q(p4) +l(p6)                  *
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'prods.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'epinv.f'
      include 'scale.f'
      include 'hardscale.f'
      logical msbar
      common/msbar/msbar
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . p(mxpart,4),q_bdkw(mxpart,4),faclo,dot,
     . fac,xl12,xl15,xl16,xl25,xl26,xl56,subqbq,subqqb,v2(2),vQ(nf,2),
     . mmsq(2,2),mmsq_vec(2,2),mmsq_ax(2,2),pswap(mxpart,4),subgg
      double complex tamp,lamp,atreez,a61z,prop
      integer nu,j,k,polq,polb,polz
      double precision ii_qg,ii_gq,if_qg,fi_qg,ff_qg,ii_gg,if_gg
      double precision msq_cs(0:2,-nf:nf,-nf:nf),mmsq_cs(0:2,2,2)
      common/msq_cols/msq_cs,mmsq_cs

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

      prop=s(3,4)/(s(3,4)-zmass**2+im*zmass*zwidth)

c---add result of integrating subtraction terms
      xl12=log(s(1,2)/musq)
      xl15=log(-s(1,5)/musq)
      xl16=log(-s(1,6)/musq)
      xl25=log(-s(2,5)/musq)
      xl26=log(-s(2,6)/musq)
      xl56=log(s(5,6)/musq)

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
        subqqb=subqqb+2*cf
        subqbq=subqbq+2*cf
      endif
      
      subqqb=subqqb*ason2pi
      subqbq=subqbq*ason2pi

c--- calculate the gg terms
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(6,nu)
      pswap(5,nu)=p(3,nu)
      pswap(6,nu)=p(4,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq,mmsq_vec,mmsq_ax)
      

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
      do polz=1,2
      do polb=1,2
        if     ((j .eq. 0) .and. (k .eq. 0)) then
          tamp=0d0
          lamp=0d0
        elseif ((j .gt. 0) .and. (k .lt. 0)) then
          tamp=atreez(polq,polb,polz,1,2,3,4,5,6,za,zb)
     .         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .        -atreez(3-polb,3-polq,polz,2,1,4,3,5,6,za,zb)
     .         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
          lamp=a61z(polq,polb,polz,1,2,3,4,5,6,za,zb)
     .         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     .        -a61z(3-polb,3-polq,polz,2,1,4,3,5,6,za,zb)
     .         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
          tamp=atreez(polq,polb,polz,4,2,3,1,5,6,za,zb)
     .         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .        -atreez(3-polb,3-polq,polz,2,4,1,3,5,6,za,zb)
     .         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
          lamp=a61z(polq,polb,polz,4,2,3,1,5,6,za,zb)
     .         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     .        -a61z(3-polb,3-polq,polz,2,4,1,3,5,6,za,zb)
     .         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
        endif
        msqv(j,k)=msqv(j,k)+fac*2d0*dble(tamp*dconjg(lamp))
      enddo
      if ((j .eq. 0) .and. (k .eq. 0)) then
        msqv(j,k)=msqv(j,k)+mmsq(polq,polz)*(
     .             cdabs(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)**2)
     .                     +mmsq_vec(polq,polz)*dreal(
     .             (Q(1)*q1+vQ(1,polq)*v2(polz)*prop)
     .            *(Q(1)*q1+0.5d0*(vQ(1,1)+vQ(1,2))*v2(polz)*prop))
     .                     +mmsq_ax(polq,polz)*dreal(
     .             (Q(1)*q1+vQ(1,polq)*v2(polz)*prop)
     .            *(v2(polz)*prop)/sin2w)
      endif
      enddo
      enddo
      
      if    ((j .gt. 0) .and. (k .lt. 0)) then
        msqv(j,k)=msqv(j,k)+subqqb*msq(j,k) 
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        msqv(j,k)=msqv(j,k)+subqbq*msq(j,k) 
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
        subgg=xn*(if_gg(one,xl15,1)+fi_qg(one,xl15,1)
     .           +if_gg(one,xl26,1)+fi_qg(one,xl26,1)
     .           +two*ii_gg(one,xl12,1))/2d0*msq_cs(1,j,k)
     .       +xn*(if_gg(one,xl25,1)+fi_qg(one,xl25,1)
     .           +if_gg(one,xl16,1)+fi_qg(one,xl16,1)
     .           +two*ii_gg(one,xl12,1))/2d0*msq_cs(2,j,k)
     .     -one/xn*(two*ff_qg(one,xl56,1))/2d0*msq_cs(0,j,k)
     .     -one/xn*(two*ff_qg(one,xl56,1))/2d0*msq_cs(1,j,k)
     .     -one/xn*(two*ff_qg(one,xl56,1))/2d0*msq_cs(2,j,k)
     .     +xn*(if_gg(one,xl15,1)+if_gg(one,xl26,1))/2d0*msq_cs(0,j,k)
     .     +xn*(fi_qg(one,xl15,1)+fi_qg(one,xl26,1))/2d0*msq_cs(0,j,k)
     .     +xn*(if_gg(one,xl16,1)+if_gg(one,xl25,1))/2d0*msq_cs(0,j,k)
     .     +xn*(fi_qg(one,xl16,1)+fi_qg(one,xl25,1))/2d0*msq_cs(0,j,k)
     .     -xn*(two*ff_qg(one,xl56,1))/2d0*msq_cs(0,j,k)
c--- add in UV counter-term here
        subgg=subgg-
     .    xn*(epinv-log(fourpi))*(11d0-2d0*dble(nf)/xn)/3d0*msq(j,k)
        if (msbar) then 
          subgg=subgg-cf*msq(j,k)
        endif
        subgg=subgg*ason2pi
c      write(*,*) 'msqv(j,k),subgg,msqv(j,k)+subgg,1d0+subgg/msqv(j,k)'
c      write(*,*) msqv(j,k),subgg,msqv(j,k)+subgg,1d0+subgg/msqv(j,k)
c      pause
        msqv(j,k)=msqv(j,k)+subgg
      endif
      
      enddo

      return
      end
     
