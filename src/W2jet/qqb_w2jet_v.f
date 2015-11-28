      subroutine qqb_w2jet_v(p,msqv)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     January 2001.                                                    *
*     Additions Aug. 2001, for the 4Q piece                            *
*      (singularities cancelled)                                       *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> W + j(p5) + j(p6)                         *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
*                                                                      *
*     where the partons are either q(p5) and qbar(p6) [Qflag = .true.] *
*                               or g(p5) and g(p6)    [Gflag = .true.] *
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the terms for the QQGG piece                     *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'prods.f'
      include 'epinv.f'
      include 'scale.f'
      include 'hardscale.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      logical msbar
      common/msbar/msbar
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . mmsq_qqb,mmsq_qbq,mmsq_gq,mmsq_gqb,mmsq_qg,mmsq_qbg,mmsq_gg,
     . p(mxpart,4),q(mxpart,4),pswap(mxpart,4),
     . fac,qqb,qbq,
     . xl12,xl15,xl16,xl25,xl26,xl56
      double complex atrLLL,atrLRL,a61LLL,a61LRL
      double complex tLLL,tLRL,fLLL,fLRL,prop
      integer nu,j,k,n1,n2,cs
      double precision ii_qg,ii_gq,ii_gg,if_qg,if_gg,fi_qg,fi_gg,
     .                 ff_qg,ff_gg,fi_gq
      double precision subqqb(0:2),subqbq(0:2),subgg(0:2),subgq(0:2),
     .                 subqg(0:2),subqbg(0:2),subgqb(0:2),subuv(0:2),
     .                 subqq(0:2),subqbqb(0:2)
      double precision facqq,facgg,Vfac,temp
      double precision mqq(0:2,fn:nf,fn:nf)
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     .                 qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji,
     .                 qbq_ijkk,qbq_iikl,qbq_ijkj,qbq_ijik,
     .                 qbq_ijii,qbq_ijjj,qbq_iiij,qbq_iiji,
     .                 qq_ijkk,qq_iikl,qq_ijkj,qq_ijik,
     .                 qq_ijii,qq_ijjj,qq_iiij,qq_iiji,
     .                 qbqb_ijkk,qbqb_iikl,qbqb_ijkj,qbqb_ijik,
     .                 qbqb_ijii,qbqb_ijjj,qbqb_iiij,qbqb_iiji
      logical first
      common/mqq/mqq
      data first/.true./
      save first

      if (first) then
        first=.false.
        if (Gflag) then
          write(*,*) 'Using QQGG (VIRTUAL) matrix elements'
          write(*,*) '[LC is     N   ]'
          write(*,*) '[SLC is   1/N  ]'
          write(*,*) '[SSLC is 1/N**3]'
        endif
        if (Qflag) then
          write(*,*) 'Using QQBQQB (VIRTUAL) matrix elements'
          write(*,*) '[LC is   1 ]'
          write(*,*) '[SLC is 1/N]'
        endif
        if     (colourchoice .eq. 1) then
          write(*,*) 'Leading colour only in VIRTUAL'
        elseif (colourchoice .eq. 2) then
          write(*,*) 'Sub-leading colour only in VIRTUAL'
        elseif (colourchoice .eq. 3) then
          write(*,*) 'Sub-sub-leading colour only in VIRTUAL'
        elseif (colourchoice .eq. 0) then
          write(*,*) 'Total of all colour structures in VIRTUAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
      endif

c--- calculate the lowest order matrix element and fill the
c--- common block twopij with s_{ij}
      call qqb_w2jetx(p,msq,mqq,msqx,msqx_cs)
      
      if  ( (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return 

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

      prop=s(3,4)/(s(3,4)-wmass**2+im*wmass*wwidth)
      
c--- endpoint contributions from subtraction terms
c--- Note that, compared to the _z terms, we must add additional
c--- contributions representing the final-final singularities
     
      xl12=log(+s(1,2)/musq)
      xl15=log(-s(1,5)/musq)
      xl16=log(-s(1,6)/musq)
      xl25=log(-s(2,5)/musq)
      xl26=log(-s(2,6)/musq)
      xl56=log(+s(5,6)/musq)

      do cs=0,2
      subgg(cs)=0d0
      subqqb(cs)=0d0
      subqbq(cs)=0d0
      subqq(cs)=0d0
      subqbqb(cs)=0d0
      subqg(cs)=0d0
      subgq(cs)=0d0
      subqbg(cs)=0d0
      subgqb(cs)=0d0
      subuv(cs)=0d0
      enddo

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************            
      if (Gflag) then
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      subgg(1)=xn*(if_gg(one,xl15,1)+fi_qg(one,xl15,1)
     .            +if_gg(one,xl26,1)+fi_qg(one,xl26,1)
     .            +two*ii_gg(one,xl12,1))/2d0
      subgg(2)=xn*(if_gg(one,xl16,1)+fi_qg(one,xl16,1)
     .            +if_gg(one,xl25,1)+fi_qg(one,xl25,1)
     .            +two*ii_gg(one,xl12,1))/2d0
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      subgg(0)=subgg(0)+
     .         xn*(if_gg(one,xl15,1)+fi_qg(one,xl15,1)
     .            +if_gg(one,xl16,1)+fi_qg(one,xl16,1)
     .            +if_gg(one,xl25,1)+fi_qg(one,xl25,1)
     .            +if_gg(one,xl26,1)+fi_qg(one,xl26,1))/2d0
      subgg(1)=subgg(1)-ff_qg(one,xl56,1)/xn
      subgg(2)=subgg(2)-ff_qg(one,xl56,1)/xn
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      subgg(0)=subgg(0)-ff_qg(one,xl56,1)*(xn+1d0/xn)
      endif

c--- subqqb and subqbq should really include the _gq terms arising
c--- from the q_gq subtraction contributions, but these are zero  
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      subqqb(1)=xn*(if_qg(one,xl15,1)+fi_gg(one,xl15,1)/2d0
     .             +if_qg(one,xl26,1)+fi_gg(one,xl26,1)/2d0
     .             +ff_gg(one,xl56,1))/2d0
      subqqb(2)=xn*(if_qg(one,xl16,1)+fi_gg(one,xl16,1)/2d0
     .             +if_qg(one,xl25,1)+fi_gg(one,xl25,1)/2d0
     .             +ff_gg(one,xl56,1))/2d0
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      subqqb(0)=subqqb(0)+
     .          xn*(if_qg(one,xl15,1)+fi_gg(one,xl15,1)/2d0
     .             +if_qg(one,xl16,1)+fi_gg(one,xl16,1)/2d0
     .             +if_qg(one,xl25,1)+fi_gg(one,xl25,1)/2d0
     .             +if_qg(one,xl26,1)+fi_gg(one,xl26,1)/2d0)/2d0
      subqqb(1)=subqqb(1)-ii_qg(one,xl12,1)/xn
      subqqb(2)=subqqb(2)-ii_qg(one,xl12,1)/xn
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      subqqb(0)=subqqb(0)-ii_qg(one,xl12,1)*(xn+1d0/xn)
      endif
     
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      subqbq(1)=xn*(if_qg(one,xl16,1)+fi_gg(one,xl16,1)/2d0
     .             +if_qg(one,xl25,1)+fi_gg(one,xl25,1)/2d0
     .             +ff_gg(one,xl56,1))/2d0
      subqbq(2)=xn*(if_qg(one,xl15,1)+fi_gg(one,xl15,1)/2d0
     .             +if_qg(one,xl26,1)+fi_gg(one,xl26,1)/2d0
     .             +ff_gg(one,xl56,1))/2d0
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      subqbq(0)=subqbq(0)+
     .          xn*(if_qg(one,xl15,1)+fi_gg(one,xl15,1)/2d0
     .             +if_qg(one,xl16,1)+fi_gg(one,xl16,1)/2d0
     .             +if_qg(one,xl25,1)+fi_gg(one,xl25,1)/2d0
     .             +if_qg(one,xl26,1)+fi_gg(one,xl26,1)/2d0)/2d0
      subqbq(1)=subqbq(1)-ii_qg(one,xl12,1)/xn
      subqbq(2)=subqbq(2)-ii_qg(one,xl12,1)/xn
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      subqbq(0)=subqbq(0)-ii_qg(one,xl12,1)*(xn+1d0/xn)
      endif
     
c--- subqg etc. should also include _gq terms from the gq_g piece
c--- but once again these are zero anyway
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      subqg(1)=xn*(2d0*ii_qg(one,xl12,1)+2d0*if_gg(one,xl26,1)
     .                +fi_gg(one,xl26,1)
     .            +2d0*ii_gg(one,xl12,1)
     .            +2d0*ff_qg(one,xl56,1)+ff_gg(one,xl56,1))/4d0
      subqg(2)=xn*(2d0*if_qg(one,xl16,1)+fi_gg(one,xl16,1)
     .            +2d0*if_gg(one,xl26,1)+fi_gg(one,xl26,1)
     .            +2d0*if_gg(one,xl25,1)
     .            +2d0*fi_qg(one,xl25,1))/4d0
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      subqg(0)=subqg(0)+
     .         xn*(2d0*if_qg(one,xl16,1)+fi_gg(one,xl16,1)
     .            +2d0*if_gg(one,xl25,1)+2d0*fi_qg(one,xl25,1)
     .            +2d0*ii_qg(one,xl12,1)+2d0*ii_gg(one,xl12,1)
     .            +2d0*ff_qg(one,xl56,1)+ff_gg(one,xl56,1))/4d0
      subqg(1)=subqg(1)-(if_qg(one,xl15,1)+fi_qg(one,xl15,1))/xn/2d0
      subqg(2)=subqg(2)-(if_qg(one,xl15,1)+fi_qg(one,xl15,1))/xn/2d0     
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      subqg(0)=subqg(0)-(xn+1d0/xn)*
     .        (if_qg(one,xl15,1)+fi_qg(one,xl15,1))/2d0
      endif
      
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      subgq(1)=xn*(2d0*ii_gg(one,xl12,1)
     .            +2d0*if_gg(one,xl16,1)+fi_gg(one,xl16,1)
     .            +2d0*ii_qg(one,xl12,1)
     .            +2d0*ff_qg(one,xl56,1)+ff_gg(one,xl56,1))/4d0
      subgq(2)=xn*(2d0*if_gg(one,xl16,1)+fi_gg(one,xl16,1)
     .            +2d0*if_gg(one,xl15,1)
     .            +2d0*fi_qg(one,xl15,1)+2d0*if_qg(one,xl26,1)
     .                +fi_gg(one,xl26,1))/4d0
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      subgq(0)=subgq(0)+
     .         xn*(2d0*if_qg(one,xl26,1)+fi_gg(one,xl26,1)
     .            +2d0*if_gg(one,xl15,1)+2d0*fi_qg(one,xl15,1)
     .            +2d0*ii_qg(one,xl12,1)+2d0*ii_gg(one,xl12,1)
     .            +2d0*ff_qg(one,xl56,1)+ff_gg(one,xl56,1))/4d0
      subgq(1)=subgq(1)-(if_qg(one,xl25,1)+fi_qg(one,xl25,1))/xn/2d0
      subgq(2)=subgq(2)-(if_qg(one,xl25,1)+fi_qg(one,xl25,1))/xn/2d0     
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      subgq(0)=subgq(0)-(xn+1d0/xn)*
     .        (if_qg(one,xl25,1)+fi_qg(one,xl25,1))/2d0
      endif
     
      subgqb(0)=subgq(0)
      subgqb(1)=subgq(2)
      subgqb(2)=subgq(1)
      
      subqbg(0)=subqg(0)
      subqbg(1)=subqg(2)
      subqbg(2)=subqg(1)
     
c--- UV counter-term
      
      if     (colourchoice .eq. 1) then
        subuv(1)=xn*(epinv-log(fourpi))*(11d0-2d0*dble(nf)/xn)/3d0
        subuv(2)=subuv(1)
      elseif (colourchoice .eq. 2) then
        subuv(0)=xn*(epinv-log(fourpi))*(11d0-2d0*dble(nf)/xn)/3d0
      elseif (colourchoice .eq. 3) then
c--- all zero already
      elseif (colourchoice .eq. 0) then
        subuv(1)=xn*(epinv-log(fourpi))*(11d0-2d0*dble(nf)/xn)/3d0
        subuv(2)=subuv(1)
        subuv(0)=subuv(1)
      endif
      
c--- End of endpoint contributions

c--- Now calculate the relevant lowest-order matrix elements
c--- for each possible initial state from the QQGG contribution
      
c---  calculate the qqb terms
CALL    0--> q(p2)+g(p5)+g(p6)+qb(p1)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qqb)

c---  calculate the qbq terms
CALL    0--> q(p1)+g(p5)+g(p6)+qb(p2)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(1,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qbq)

c---  calculate the gq terms
CALL    0--> q(p5)+g(p1)+g(p6)+qb(p2)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gq)

c---  calculate the qg terms
CALL    0--> q(p5)+g(p2)+g(p6)+qb(p1)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qg)

c---  calculate the gqb terms
CALL    0--> q(p2)+g(p1)+g(p6)+qb(p5)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gqb)

c---  calculate the qbg terms
CALL    0--> q(p1)+g(p2)+g(p6)+qb(p5)+l(p3)+lbar(p4) 
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(1,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qbg)

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
      call xwqqgg_v(mmsq_gg)
      endif
      
************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************            
      if (Qflag) then
      subqq(0)=(
     . -(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*one/xn
     . -(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*one/xn
     . +ii_qg(one,xl12,1)*(xn+one/xn)
     . +ff_qg(one,xl56,1)*(xn+one/xn)
     . -(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*one/xn
     . -(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*one/xn
     . +ii_qg(one,xl12,1)*(xn+one/xn)
     . +ff_qg(one,xl56,1)*(xn+one/xn))/2d0
      subqq(1)=(
     . -(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*one/xn
     . +(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*(xn-two/xn)
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn
     . +(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*(xn-two/xn)
     . -(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*one/xn
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn)/2d0
      subqq(2)=(
     . +(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*(xn-two/xn)
     . -(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*one/xn
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn
     . -(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*one/xn
     . +(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*(xn-two/xn)
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn)/2d0
     
      subqbqb(0)=(
     . -(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*one/xn
     . -(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*one/xn
     . +ii_qg(one,xl12,1)*(xn+one/xn)
     . +ff_qg(one,xl56,1)*(xn+one/xn)
     . -(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*one/xn
     . -(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*one/xn
     . +ii_qg(one,xl12,1)*(xn+one/xn)
     . +ff_qg(one,xl56,1)*(xn+one/xn))/2d0
      subqbqb(1)=(
     . -(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*one/xn
     . +(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*(xn-two/xn)
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn
     . +(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*(xn-two/xn)
     . -(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*one/xn
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn)/2d0
      subqbqb(2)=(
     . +(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*(xn-two/xn)
     . -(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*one/xn
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn
     . -(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*one/xn
     . +(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*(xn-two/xn)
     . +ii_qg(one,xl12,1)*two/xn
     . +ff_qg(one,xl56,1)*two/xn)/2d0
 
      subqqb(0)=(
     . -(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*one/xn
     . +(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*(xn+one/xn)
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn
     . +(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*(xn+one/xn)
     . -(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*one/xn
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn)/2d0
      subqqb(1)=(
     . +(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*(xn-two/xn)
     . +(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*two/xn
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn
     . +(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*two/xn
     . +(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*(xn-two/xn)
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn)/2d0
      subqqb(2)=(
     . -(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*one/xn
     . +(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*two/xn
     . +ii_qg(one,xl12,1)*(xn-two/xn)
     . +ff_qg(one,xl56,1)*(xn-two/xn)
     . +(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*two/xn
     . -(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*one/xn
     . +ii_qg(one,xl12,1)*(xn-two/xn)
     . +ff_qg(one,xl56,1)*(xn-two/xn))/2d0
      
      subqbq(0)=(
     . +(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*(xn+one/xn)
     . -(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*one/xn
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn
     . -(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*one/xn
     . +(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*(xn+one/xn)
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn)/2d0
      subqbq(1)=(
     . +(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*two/xn
     . +(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*(xn-two/xn)
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn
     . +(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*(xn-two/xn)
     . +(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*two/xn
     . -ii_qg(one,xl12,1)*one/xn
     . -ff_qg(one,xl56,1)*one/xn)/2d0
      subqbq(2)=(
     . +(if_qg(one,xl15,1)+fi_qg(one,xl15,1))*two/xn
     . -(if_qg(one,xl16,1)+fi_qg(one,xl16,1))*one/xn
     . +ii_qg(one,xl12,1)*(xn-two/xn)
     . +ff_qg(one,xl56,1)*(xn-two/xn)
     . -(if_qg(one,xl25,1)+fi_qg(one,xl25,1))*one/xn
     . +(if_qg(one,xl26,1)+fi_qg(one,xl26,1))*two/xn
     . +ii_qg(one,xl12,1)*(xn-two/xn)
     . +ff_qg(one,xl56,1)*(xn-two/xn))/2d0

c--- UV counter-term is already included in a6routine.f
      subuv(1)=0d0
      subuv(2)=subuv(1)
      subuv(0)=subuv(1)
      
c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=mqq(cs,j,k)
        enddo
        enddo
      enddo
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
      fac=V*xn*gw**4*gsq**2*ason2pi

c--- set-up qqb matrix elements
      call qqbw2j_loop(1,2,3,4,5,6,qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     .                             qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji)
c--- qbq
      call qqbw2j_loop(4,2,3,1,5,6,qbq_ijkk,qbq_iikl,qbq_ijkj,qbq_ijik,
     .                             qbq_ijii,qbq_ijjj,qbq_iiij,qbq_iiji)
c--- qq (note that roles of iiij and ijii are reversed)
      call qqbw2j_loop(2,1,3,4,5,6,qq_ijkk,qq_iikl,qq_ijkj,qq_ijik,
     .                             qq_iiij,qq_ijjj,qq_ijii,qq_iiji)
c--- qbqb (note that roles of ijjj and iiji are reversed)
      call qqbw2j_loop(1,2,4,3,5,6,qbqb_ijkk,qbqb_iikl,qbqb_ijkj,
     .         qbqb_ijik,qbqb_ijii,qbqb_iiji,qbqb_iiij,qbqb_ijjj)

      endif

c--- Add VIRTUAL terms

      do j=-nf,nf
      do k=-nf,nf

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************            
      if (Gflag) then
      if     ((j .gt. 0) .and. (k .lt. 0)) then
        msqv(j,k)=msqv(j,k)+Vsq(j,k)*mmsq_qqb*cdabs(prop)**2
     .           *half*(aveqq/avegg)*(gwsq**2/4d0/esq**2)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        msqv(j,k)=msqv(j,k)+Vsq(j,k)*mmsq_qbq*cdabs(prop)**2
     .           *half*(aveqq/avegg)*(gwsq**2/4d0/esq**2)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_qg*cdabs(prop)**2*
     .           (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_qbg*cdabs(prop)**2*
     .           (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_gq*cdabs(prop)**2*
     .           (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_gqb*cdabs(prop)**2*
     .           (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))
     .           *(aveqg/avegg)*(gwsq**2/4d0/esq**2)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
        Vfac=0d0
        do n1=1,nf
          do n2=-nf,-1
            Vfac=Vfac+Vsq(n1,n2)
          enddo
        enddo
        msqv(j,k)=msqv(j,k)
     .           +mmsq_gg*Vfac*cdabs(prop)**2*(gwsq**2/4d0/esq**2)
      endif
      endif
      
      if (Qflag) then
      if     ((j .gt. 0) .and. (k .lt. 0)) then
        if (j .ne. -k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsq(j,k)*(qqb_ijii+qqb_ijjj)
     .     +Vsq(j,k)*dfloat(nf-2)*qqb_ijkk
     .     +(Vsum(j)-Vsq(j,k))*qqb_ijkj
     .     +(Vsum(k)-Vsq(j,k))*qqb_ijik)
        else
          Vfac=0d0
          do n1=1,nf
          do n2=-nf,-1
          if ((n1 .ne. j) .and. (n2 .ne. k)) then
            Vfac=Vfac+Vsq(n1,n2)
          endif
          enddo
          enddo
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsum(k)*qqb_iiij
     .     +Vsum(j)*qqb_iiji
     .     +Vfac*qqb_iikl)
        endif
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        if (j .ne. -k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsq(j,k)*(qbq_ijii+qbq_ijjj)
     .     +Vsq(j,k)*dfloat(nf-2)*qbq_ijkk
     .     +(Vsum(k)-Vsq(j,k))*qbq_ijkj
     .     +(Vsum(j)-Vsq(j,k))*qbq_ijik)
        else
          Vfac=0d0
          do n1=-nf,-1
          do n2=1,nf
          if ((n1 .ne. j) .and. (n2 .ne. k)) then
            Vfac=Vfac+Vsq(n1,n2)
          endif
          enddo
          enddo
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsum(j)*qbq_iiij
     .     +Vsum(k)*qbq_iiji
     .     +Vfac*qbq_iikl)
        endif
      elseif ((j .gt. 0) .and. (k .gt. 0)) then
        if (j .ne. k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsq(j,-k)*half*qq_ijjj
     .     +Vsq(k,-j)*half*qq_ijii
     .     +(Vsum(j)-Vsq(j,-k))*qq_ijkj
     .     +(Vsum(k)-Vsq(k,-j))*qq_ijik)
        else
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsum(j)*qq_iiji)
        endif
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
        if (j .ne. k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsq(j,-k)*half*qbqb_ijjj
     .     +Vsq(k,-j)*half*qbqb_ijii
     .     +(Vsum(j)-Vsq(j,-k))*qbqb_ijkj
     .     +(Vsum(k)-Vsq(k,-j))*qbqb_ijik)
        else
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     .      Vsum(j)*qbqb_iiji)
        endif
      endif
      endif
      
      enddo
      enddo
  
c--- Add ENDPOINT terms for QQGG contributions

      do j=-nf,nf
      do k=-nf,nf

      temp=msqv(j,k)

      if     ((j .gt. 0) .and. (k .lt. 0)) then
        do cs=0,2
        if (subqqb(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subqqb(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        do cs=0,2
        if (subqbq(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subqbq(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .gt. 0) .and. (k .gt. 0)) then
        do cs=0,2
        if (subqq(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subqq(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
        do cs=0,2
        if (subqbqb(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subqbqb(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
        do cs=0,2
        if (subqg(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subqg(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
        do cs=0,2
        if (subqbg(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subqbg(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        do cs=0,2
        if (subgq(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subgq(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        do cs=0,2
        if (subgqb(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subgqb(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
        do cs=0,2
        if (subgg(cs) .ne. 0d0) then
          msqv(j,k)=msqv(j,k)+
     .      ason2pi*(subgg(cs)-subuv(cs))*msq_cs(cs,j,k)
        endif
        enddo
      endif

c      if ((temp .ne. 0d0)) then
c      write(*,91) 'j, k, msqv, sub',j,k,temp,msqv(j,k)-temp
c     .      ,abs(temp/(msqv(j,k)-temp))-1d0,msqv(j,k)
c      endif

      enddo
      enddo

c      write(*,*) 'Using epinv =',epinv
c      pause
 
   91 format(a16,2i3,2e16.7,f12.8,e18.10)
      
      return
      end
     
     
 
