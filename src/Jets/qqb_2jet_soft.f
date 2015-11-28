      subroutine qqb_2jet_soft(p,msq)
************************************************************************
*                                                                      *
*  Author: J.M. Campbell, October 2002                                 *
*                                                                      *
*  Soft approximation to matrix elements for 3 jet production          *
*     f(-p1) + f(-p2) --> f(p3) + f(p4) + f(p5)                        *
*                                                                      *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'jetlabel.f'
      integer i,j,k,l,m,n,nqcdjets,nqcdstart,a(6),b(6),c(6),qq
      integer perm(3:5,3:5),nid
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),pjet(mxpart,4),
     . msq0(-nf:nf,-nf:nf),s(mxpart,mxpart),Rcut,dum(36)
      double precision eik1a_b(6),eikba_1(6),eikab_c(6),
     .                 eikcb_a(6),eikbc_2(6),eik2c_b(6),
     .                 eik1b_2(6),eik2b_1(6),eikf,eikc,
     .                 msq_ac(6,-nf:nf,-nf:nf),p_ac(mxpart,4),
     .                 msq_acx(6,0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     .                 msq_ac_1d(6,-nf:nf,-nf:nf),
     .                 msq_ac_2d(6,-nf:nf,-nf:nf),
     .                 msq_ac_3d(6,-nf:nf,-nf:nf),
     .                 msq_ac_1b(6,-nf:nf,-nf:nf),
     .                 msq_ac_1c(6,-nf:nf,-nf:nf),
     .                 msq_ac_2c(6,-nf:nf,-nf:nf)
      double precision mqq(0:2,fn:nf,fn:nf)
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision dot,ss,tt,uu,smallb_1,smallc_1,smallc_2,smalld_1,
     . xmsq
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/eikeik/eik1a_b,eikba_1,eik2c_b,eikbc_2,eikab_c,eikcb_a,
     . eik1b_2,eik2b_1,msq_ac,msq_ac_1b,msq_ac_1d,msq_ac_2d,msq_ac_3d,
     . msq_acx
      parameter(qq=1)
      data a/3,3,5,4,4,5/
      data b/4,5,3,5,3,4/
      data c/5,4,4,3,5,3/
      data perm/0,4,6,2,0,3,1,5,0/
      
c      call dip2j_qrqr(1,1,1,1,1,dum,1,1,1,1)
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do j=1,2
      do k=1,4
      p_ac(j,k)=p(j,k)
      enddo
      enddo

      call dotem(5,p,s)

c---- gg --> ggg singularities
      do n=1,6
      eik1a_b(n)=+s(1,b(n))/(-s(1,a(n))+s(a(n),b(n)))/s(1,a(n))*two*gsq
      eikba_1(n)=-s(b(n),1)/(s(b(n),a(n))-s(a(n),1))
     .             /s(b(n),a(n))*two*gsq
      eikab_c(n)=+s(a(n),c(n))/(s(a(n),b(n))+s(b(n),c(n)))
     .             /s(a(n),b(n))*two*gsq
      eikcb_a(n)=+s(c(n),a(n))/(s(c(n),b(n))+s(b(n),a(n)))
     .             /s(c(n),b(n))*two*gsq
      eikbc_2(n)=-s(b(n),2)/(s(b(n),c(n))-s(c(n),2))
     .             /s(b(n),c(n))*two*gsq
      eik2c_b(n)=+s(2,b(n))/(-s(2,c(n))+s(c(n),b(n)))/s(2,c(n))*two*gsq
      eik1b_2(n)=+s(1,2)/(s(1,b(n))+s(2,b(n)))/s(1,b(n))*two*gsq
      eik2b_1(n)=+s(1,2)/(s(1,b(n))+s(2,b(n)))/s(2,b(n))*two*gsq

      do k=1,4
        p_ac(3,k)=p(a(n),k)
        p_ac(4,k)=p(c(n),k)
        p_ac(5,k)=0d0
      enddo

      call qqb_2jetx(p_ac,msq0,msqx)
      ss=2d0*dot(p_ac,1,2)
      tt=2d0*dot(p_ac,1,3)
      uu=2d0*dot(p_ac,2,3)
      call genclust2(p_ac,rcut,pjet,0)
      
      do j=-nf,nf
      do k=-nf,nf
        if (jets .ne. nqcdjets) then
          msq_ac(n,j,k)=0d0
          do l=-nf,nf
          do m=-nf,nf
          do i=0,2
            msqx(i,j,k,l,m)=0d0
          enddo
          enddo
          enddo          
        else
          msq_ac(n,j,k)=msq0(j,k)
          do l=-nf,nf
          do m=-nf,nf
          do i=0,2
            msq_acx(n,i,j,k,l,m)=msqx(i,j,k,l,m)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo
      
      enddo

c--- sum up contributions      
      do j=-nf,nf
      do k=-nf,nf
      
c--- gluon-gluon
      if     ((j .eq. 0) .and. (k.eq.0)) then
        call eik2jet_gggg(1,2,3,4,5,xmsq)
        msq(j,k)=xmsq
        call eik2jet_qqgg(4,3,2,1,5,xmsq,0,0,1,-1)
        msq(j,k)=msq(j,k)+dfloat(nf)*xmsq

c--- quark-quark or antiquark-antiquark
      elseif (((j .gt. 0) .and. (k. gt. 0)) .or.
     .        ((j .lt. 0) .and. (k. lt. 0))) then
c---     (identical)
        if (j .eq. k) then
          call eik2jet_qqqq(1,2,3,4,5,xmsq,1,1,1,1)
          msq(j,k)=xmsq
c---     (non-identical)
        else
          call eik2jet_qrqr(1,2,3,4,5,xmsq,1,2,1,2)
          msq(j,k)=xmsq
        endif

c--- quark-antiquark
      elseif ((j .gt. 0) .and. (k. lt. 0)) then
c---     (identical)
        if (j .eq. -k) then
          call eik2jet_qqqq(1,4,3,2,5,xmsq,2,-2,2,-2)
          msq(j,k)=xmsq
          call eik2jet_qrqr(1,4,2,3,5,xmsq,2,-2,1,-1)
          msq(j,k)=msq(j,k)+dfloat(nf-1)*xmsq
          call eik2jet_qqgg(2,1,4,3,5,xmsq,1,-1,0,0)
          msq(j,k)=msq(j,k)+xmsq/3d0
c---     (non-identical)
        else
          call eik2jet_qrqr(1,4,3,2,5,xmsq,1,2,1,2)
          msq(j,k)=xmsq
        endif

c--- antiquark-quark
      elseif ((j .lt. 0) .and. (k. gt. 0)) then
c---     (identical)
        if (j .eq. -k) then
          call eik2jet_qqqq(2,4,3,1,5,xmsq,-2,2,2,-2)
          msq(j,k)=xmsq
          call eik2jet_qrqr(2,4,1,3,5,xmsq,-2,2,1,-1)
          msq(j,k)=msq(j,k)+dfloat(nf-1)*xmsq
          call eik2jet_qqgg(1,2,4,3,5,xmsq,-1,1,0,0)
          msq(j,k)=msq(j,k)+xmsq/3d0
c---     (non-identical)
        else
          call eik2jet_qrqr(2,4,3,1,5,xmsq,1,2,1,2)
          msq(j,k)=xmsq
        endif

c--- quark-gluon
      elseif ((j .gt. 0) .and. (k. eq. 0)) then
        call eik2jet_qqgg(3,1,2,4,5,xmsq,2,0,2,0)
        msq(j,k)=msq(j,k)+xmsq/2d0

      elseif ((j .lt. 0) .and. (k. eq. 0)) then
        call eik2jet_qqgg(1,3,2,4,5,xmsq,-2,0,-2,0)
        msq(j,k)=msq(j,k)+xmsq/2d0

      elseif ((j .eq. 0) .and. (k. gt. 0)) then
        call eik2jet_qqgg(2,3,1,4,5,xmsq,0,2,0,2)
        msq(j,k)=msq(j,k)+xmsq/2d0

      elseif ((j .eq. 0) .and. (k. lt. 0)) then
        call eik2jet_qqgg(3,2,1,4,5,xmsq,0,-2,0,-2)
        msq(j,k)=msq(j,k)+xmsq/2d0

      endif  

      enddo
      enddo          
          
      return
      end
      
