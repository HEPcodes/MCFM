      subroutine qqb_z2jet_gs(p,msq)
************************************************************************
*     Author: R.K.Ellis                                                *
*     April, 2001.                                                     *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z + f(p5) + f(p6) + g(p7)
c                           |
c                           --> l(p3) + lbar(p4)
c     where the fermions are either q(p5) and qbar(p6) [Qflag = .true.]
c                                or g(p5) and g(p6)    [Gflag = .true.]
c--- all momenta are incoming

      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'msq_cs.f'
      include 'flags.f'
      integer j,k,mj,mk,n,nd,ii,np6,np12,np18,np21
      integer a(6),b(6),c(6),pntr(5:7,5:7)
      double precision p(mxpart,4),msq(maxd,fn:nf,fn:nf),dot,temp
      double precision msumsq(fn:nf,fn:nf),
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq16_2(fn:nf,fn:nf),msq26_1(fn:nf,fn:nf),
     & msq17_2(fn:nf,fn:nf),msq27_1(fn:nf,fn:nf),
     & msq17_5(fn:nf,fn:nf),msq17_5v(fn:nf,fn:nf),
     & msq17_6(fn:nf,fn:nf),msq17_6v(fn:nf,fn:nf),
     & msq15_6(fn:nf,fn:nf),msq15_6v(fn:nf,fn:nf),
     & msq27_5(fn:nf,fn:nf),msq27_5v(fn:nf,fn:nf),
     & msq27_6(fn:nf,fn:nf),msq27_6v(fn:nf,fn:nf),
     & msq57_6(fn:nf,fn:nf),msq67_5(fn:nf,fn:nf),
     & msq57_1(fn:nf,fn:nf),msq67_1(fn:nf,fn:nf),
     & msq57_2(fn:nf,fn:nf),msq67_2(fn:nf,fn:nf),
     & msq65_1(fn:nf,fn:nf),msq65_1v(fn:nf,fn:nf),
     & msq17_2v(fn:nf,fn:nf),msq27_1v(fn:nf,fn:nf),
     & msq57_6v(fn:nf,fn:nf),msq67_5v(fn:nf,fn:nf),
     & msq57_1v(fn:nf,fn:nf),msq67_1v(fn:nf,fn:nf),
     & msq57_2v(fn:nf,fn:nf),msq67_2v(fn:nf,fn:nf),
     & dummy(fn:nf,fn:nf),dummyv(fn:nf,fn:nf),dummy1(fn:nf,fn:nf),
     & sub15_2(4),sub25_1(4),
     & sub16_2(4),sub26_1(4),
     & sub17_2(4),sub27_1(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),
     & sub15_6(4),sub65_1(4),
     & sub57_6(4),sub67_5(4),dsubv,dsub(4),
     & sub17_2v,sub27_1v,sub17_5v,sub17_6v,sub27_5v,sub27_6v,
     & sub57_6v,sub67_5v,sub57_1v,sub57_2v,sub67_1v,sub67_2v
      double precision
     & mgg17_2(0:2),mgg27_1(0:2),mgg57_6(0:2),mgg67_5(0:2),
     & mgg56_7(0:2),mgg76_5(0:2),mgg17_5(0:2),mgg57_1(0:2),
     & mgg17_6(0:2),mgg67_1(0:2),mgg15_6(0:2),mgg65_1(0:2),
     & mgg27_5(0:2),mgg57_2(0:2),mgg27_6(0:2),mgg67_2(0:2)
      double precision
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m16_2(0:2,fn:nf,fn:nf),m26_1(0:2,fn:nf,fn:nf),
     & m17_2(0:2,fn:nf,fn:nf),m27_1(0:2,fn:nf,fn:nf),
     & m57_6(0:2,fn:nf,fn:nf),m67_5(0:2,fn:nf,fn:nf),
     & m17_5(0:2,fn:nf,fn:nf),
     & m17_6(0:2,fn:nf,fn:nf),
     & m27_5(0:2,fn:nf,fn:nf),
     & m27_6(0:2,fn:nf,fn:nf)
      double precision
     & isq15_2(fn:nf,fn:nf),isq25_1(fn:nf,fn:nf),
     & isq16_2(fn:nf,fn:nf),isq26_1(fn:nf,fn:nf),
     & isq17_2(fn:nf,fn:nf),isq27_1(fn:nf,fn:nf),
     & isq57_6(fn:nf,fn:nf),isq67_5(fn:nf,fn:nf),
     & isq17_5(fn:nf,fn:nf),
     & isq17_6(fn:nf,fn:nf),
     & isq27_5(fn:nf,fn:nf),
     & isq27_6(fn:nf,fn:nf)
      double precision
     & mgg17_2v(0:2),mgg27_1v(0:2),mgg57_6v(0:2),mgg67_5v(0:2),
     & mgg56_7v(0:2),mgg76_5v(0:2),mgg17_5v(0:2),mgg57_1v(0:2),
     & mgg17_6v(0:2),mgg67_1v(0:2),mgg15_6v(0:2),mgg65_1v(0:2),
     & mgg27_5v(0:2),mgg57_2v(0:2),mgg27_6v(0:2),mgg67_2v(0:2)
      double precision 
     . msq1a_b(6,0:2,fn:nf,fn:nf),msqba_1(6,0:2,fn:nf,fn:nf),
     . msqab_c(6,0:2,fn:nf,fn:nf),msqcb_a(6,0:2,fn:nf,fn:nf),
     . msqbc_2(6,0:2,fn:nf,fn:nf),msq2c_b(6,0:2,fn:nf,fn:nf),
     . msq56_1(6,0:2,fn:nf,fn:nf),msq56_2(6,0:2,fn:nf,fn:nf),
     . msq56_7(6,0:2,fn:nf,fn:nf),
     . sub1a_b(6,4),subba_1(6,4),
     . subab_c(6,4),subcb_a(6,4),
     . subbc_2(6,4),sub2c_b(6,4),
     . sub56_1(6,4),sub56_2(6,4),
     . sub56_7(6,4),
     . msq1a_bv(6,0:2,fn:nf,fn:nf),msqba_1v(6,0:2,fn:nf,fn:nf),
     . msqab_cv(6,0:2,fn:nf,fn:nf),msqcb_av(6,0:2,fn:nf,fn:nf),
     . msqbc_2v(6,0:2,fn:nf,fn:nf),msq2c_bv(6,0:2,fn:nf,fn:nf),
     . sub1a_bv(6),subba_1v(6),
     . subab_cv(6),subcb_av(6),
     . subbc_2v(6),sub2c_bv(6),
     . msq1b_2(6,0:2,fn:nf,fn:nf),msq2b_1(6,0:2,fn:nf,fn:nf),
     . sub1b_2(6,4),sub2b_1(6,4),
     . msq1b_2v(6,0:2,fn:nf,fn:nf),msq2b_1v(6,0:2,fn:nf,fn:nf),
     . sub1b_2v(6),sub2b_1v(6),
     . msq56_1v(6,0:2,fn:nf,fn:nf),msq56_2v(6,0:2,fn:nf,fn:nf),
     . msq56_7v(6,0:2,fn:nf,fn:nf),
     . sub56_1v(6),sub56_2v(6),sub56_7v(6)
      external qqb_z2jet,qqb_z2jet_gvec,donothing_gvec
      common/msumsq/msumsq
      data a/5,5,7,6,6,7/
      data b/6,7,5,7,5,6/
      data c/7,6,6,5,7,5/
      data pntr/0,2,2,1,0,2,1,1,0/

      ndmax=24
      if (Gflag) then
c--- new final-final
      do n=1,6
      call dips(n,p,a(n),b(n),c(n),dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msqab_c,msqab_cv,dsub,dsubv,subab_c,subab_cv,n)
      enddo

c--- new initial-final/final-initial
      do n=1,6
      np6=n+6
      np12=n+12
      call dips(np6,p,1,a(n),b(n),dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq1a_b,msq1a_bv,dsub,dsubv,sub1a_b,sub1a_bv,n)
      call dips(np6,p,b(n),a(n),1,dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msqba_1,msqba_1v,dsub,dsubv,subba_1,subba_1v,n)
      call dips(np12,p,2,c(n),b(n),dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq2c_b,msq2c_bv,dsub,dsubv,sub2c_b,sub2c_bv,n)
      call dips(np12,p,b(n),c(n),2,dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msqbc_2,msqbc_2v,dsub,dsubv,subbc_2,subbc_2v,n)
      enddo

c--- new initial-initial
      do n=1,3
      np18=n+18
      np21=n+21
      call dips(np18,p,1,b(n),2,dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq1b_2,msq1b_2v,dsub,dsubv,sub1b_2,sub1b_2v,n)
      call dips(np21,p,2,b(n),1,dsub,dsubv,dummy,dummyv,
     . qqb_z2jet,qqb_z2jet_gvec)
      call storedip(msq2b_1,msq2b_1v,dsub,dsubv,sub2b_1,sub2b_1v,n)
      enddo

c-- fill the matrix elements
      do j=-nf,nf
      do k=-nf,nf
      
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo

c--- QUARK-ANTIQUARK contributions
      if    ((j.gt.0).and.(k.lt.0)) then
          do n=1,6   
          msq(n,j,k)   =subab_c(n,gg)/2d0
     .                  *(msqab_c(n,pntr(a(n),c(n)),j,k)
     .                   +msqab_c(n,pntr(c(n),a(n)),j,k))*xn/3d0
     .                 +subab_cv(n)/2d0
     .                  *(msqab_cv(n,pntr(a(n),c(n)),j,k)
     .                   +msqab_cv(n,pntr(c(n),a(n)),j,k))*xn/3d0
         msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2d0)
     .                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/3d0
     .                 +subba_1v(n)/2d0
     .                  *msqba_1v(n,pntr(b(n),c(n)),j,k)*xn/3d0
          msq(12+n,j,k)=(subbc_2(n,gg)/2d0+sub2c_b(n,qq))
     .                  *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn/3d0
     .                 +subbc_2v(n)/2d0
     .                  *msqbc_2v(n,pntr(a(n),b(n)),j,k)*xn/3d0
          enddo
c--- ANTIQUARK-QUARK contributions
      elseif((j.lt.0).and.(k.gt.0)) then
          do n=1,6   
          msq(n,j,k)   =subab_c(n,gg)/2d0
     .                  *(msqab_c(n,pntr(a(n),c(n)),j,k)
     .                   +msqab_c(n,pntr(c(n),a(n)),j,k))*xn/3d0
     .                 +subab_cv(n)/2d0
     .                  *(msqab_cv(n,pntr(a(n),c(n)),j,k)
     .                   +msqab_cv(n,pntr(c(n),a(n)),j,k))*xn/3d0
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2d0)
     .                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/3d0
     .                 +subba_1v(n)/2d0
     .                  *msqba_1v(n,pntr(c(n),b(n)),j,k)*xn/3d0
          msq(12+n,j,k)=(subbc_2(n,gg)/2d0+sub2c_b(n,qq))
     .                  *msq2c_b(n,pntr(b(n),a(n)),j,k)*xn/3d0
     .                 +subbc_2v(n)/2d0
     .                  *msqbc_2v(n,pntr(b(n),a(n)),j,k)*xn/3d0
          enddo
c--- GLUON-GLUON contributions
      elseif ((j.eq.0) .and. (k.eq.0)) then
c--- choose n=2 which is (5,7,6)
          do n=2,2
          msq(18+n,j,k)=sub1b_2(n,gg)
     .                  *(msq1b_2(n,1,j,k)+msq1b_2(n,2,j,k))*xn
     .                 +sub1b_2v(n)
     .                  *(msq1b_2v(n,1,j,k)+msq1b_2v(n,2,j,k))*xn
          msq(21+n,j,k)=sub2b_1(n,gg)
     .                  *(msq2b_1(n,1,j,k)+msq2b_1(n,2,j,k))*xn
     .                 +sub2b_1v(n)
     .                  *(msq2b_1v(n,1,j,k)+msq2b_1v(n,2,j,k))*xn
          enddo
c--- choose n=3,6 which is (7,5,6) and (7,6,5)
          do n=3,6,3
          msq(6+n,j,k)=(sub1a_b(n,gg)+subba_1(n,qq))
     .                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn
     .                 +sub1a_bv(n)
     .                  *msq1a_bv(n,pntr(b(n),c(n)),j,k)*xn
          enddo
c--- choose n=1,5 which is (5,6,7) and (6,5,7)
          do n=1,5,4
          msq(12+n,j,k)=(sub2c_b(n,gg)+subbc_2(n,qq))
     .                   *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn
     .                  +sub2c_bv(n)
     .                   *msq2c_bv(n,pntr(a(n),b(n)),j,k)*xn
         enddo
c--- additional (qg) collinear contributions
c--- choose n=3 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     .      (msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     .      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     .      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     .      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))*xn*(avegg/aveqg)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     .      (msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     .      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     .      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     .      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))*xn*(avegg/aveqg)
          enddo
c--- choose n=1 which is (5,6,7)
          do n=1,1
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     .      (msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     .      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     .      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     .      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))*xn*(avegg/aveqg)
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     .      (msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     .      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     .      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     .      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))*xn*(avegg/aveqg)
          enddo
c--- QUARK-GLUON contributions
      elseif ((j.gt.0) .and. (k.eq.0)) then
          do n=1,2  
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     .                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(12+n,j,k)=(subbc_2(n,gg)/2d0+sub2c_b(n,gg))
     .                  *(msq2c_b(n,pntr(a(n),b(n)),j,k)
     .                   +msq2c_b(n,pntr(b(n),a(n)),j,k))*xn/2d0
     .                 +subbc_2v(n)/2d0
     .                  *(msqbc_2v(n,pntr(a(n),b(n)),j,k)
     .                   +msqbc_2v(n,pntr(b(n),a(n)),j,k))*xn/2d0
     .                 +sub2c_bv(n)
     .                  *(msq2c_bv(n,pntr(a(n),b(n)),j,k)
     .                   +msq2c_bv(n,pntr(b(n),a(n)),j,k))*xn/2d0
          msq(18+n,j,k)=sub1b_2(n,qq)
     .                  *msq1b_2(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(21+n,j,k)=sub2b_1(n,gg)
     .                  *msq2b_1(n,pntr(a(n),c(n)),j,k)*xn/2d0
     .                 +sub2b_1v(n)
     .                  *msq2b_1v(n,pntr(a(n),c(n)),j,k)*xn/2d0
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2d0
     .                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2d0
     .                 +subab_cv(n)/2d0
     .                  *msqab_cv(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2d0)
     .                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/2d0
     .                 +subba_1v(n)/2d0
     .                  *msqba_1v(n,pntr(b(n),c(n)),j,k)*xn/2d0
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=(subbc_2(n,qq)+sub2c_b(n,gg))
     .                  *msq2c_b(n,pntr(a(n),b(n)),j,k)*xn/2d0
     .                 +sub2c_bv(n)
     .                  *msq2c_bv(n,pntr(a(n),b(n)),j,k)*xn/2d0
          enddo      
c--- additional (qg) collinear contributions
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     .      (msq2b_1(n,1,j,-1)+msq2b_1(n,1,j,-2)+msq2b_1(n,1,j,-3)
     .      +msq2b_1(n,1,j,-4)+msq2b_1(n,1,j,-5)
     .      +msq2b_1(n,2,j,-1)+msq2b_1(n,2,j,-2)+msq2b_1(n,2,j,-3)
     .      +msq2b_1(n,2,j,-4)+msq2b_1(n,2,j,-5))*xn*(aveqg/aveqq)
          enddo
c--- GLUON-QUARK contributions
      elseif ((j.eq.0) .and. (k.gt.0)) then
          do n=1,2  
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     .                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(12+n,j,k) =(sub2c_b(n,qq)+subbc_2(n,gg)/2d0)
     .                    *msq2c_b(n,pntr(c(n),a(n)),j,k)*xn/2d0
     .                   +subbc_2v(n)/2d0
     .                    *msqbc_2v(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(18+n,j,k)=sub1b_2(n,gg)
     .                  *msq1b_2(n,pntr(a(n),c(n)),j,k)*xn/2d0
     .                 +sub1b_2v(n)
     .                  *msq1b_2v(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(21+n,j,k)=sub2b_1(n,qq)
     .                  *msq2b_1(n,pntr(a(n),c(n)),j,k)*xn/2d0
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2d0
     .                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2d0
     .                 +subab_cv(n)/2d0
     .                  *msqab_cv(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(6+n,j,k)=(subba_1(n,gg)/2d0+sub1a_b(n,gg))
     .                 *(msq1a_b(n,pntr(c(n),b(n)),j,k)
     .                  +msq1a_b(n,pntr(b(n),c(n)),j,k))*xn/2d0
     .                +subba_1v(n)/2d0
     .                 *(msqba_1v(n,pntr(c(n),b(n)),j,k)
     .                  +msqba_1v(n,pntr(b(n),c(n)),j,k))*xn/2d0
     .                +sub1a_bv(n)
     .                 *(msq1a_bv(n,pntr(c(n),b(n)),j,k)
     .                  +msq1a_bv(n,pntr(b(n),c(n)),j,k))*xn/2d0
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,5,7) and (7,5,6)
          msq(6+n,j,k)=(subba_1(n,qq)+sub1a_b(n,gg))
     .                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/2d0
     .                 +sub1a_bv(n)
     .                  *msq1a_bv(n,pntr(c(n),b(n)),j,k)*xn/2d0
          enddo              
c--- additional (qg) collinear contributions
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     .      (msq1b_2(n,1,-1,k)+msq1b_2(n,1,-2,k)+msq1b_2(n,1,-3,k)
     .      +msq1b_2(n,1,-4,k)+msq1b_2(n,1,-5,k)
     .      +msq1b_2(n,2,-1,k)+msq1b_2(n,2,-2,k)+msq1b_2(n,2,-3,k)
     .      +msq1b_2(n,2,-4,k)+msq1b_2(n,2,-5,k))*xn*(aveqg/aveqq)
          enddo
c--- GLUON-ANTIQUARK contributions
      elseif ((j.eq.0) .and. (k.lt.0)) then
          do n=1,2  
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     .                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(12+n,j,k) =(sub2c_b(n,qq)+subbc_2(n,gg)/2d0)
     .                    *msq2c_b(n,pntr(a(n),c(n)),j,k)*xn/2d0
     .                   +subbc_2v(n)/2d0
     .                    *msqbc_2v(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(18+n,j,k)=sub1b_2(n,gg)
     .                  *msq1b_2(n,pntr(c(n),a(n)),j,k)*xn/2d0
     .                 +sub1b_2v(n)
     .                  *msq1b_2v(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(21+n,j,k)=sub2b_1(n,qq)
     .                  *msq2b_1(n,pntr(c(n),a(n)),j,k)*xn/2d0
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2d0
     .                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2d0
     .                 +subab_cv(n)/2d0
     .                  *msqab_cv(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(6+n,j,k)=(subba_1(n,gg)/2d0+sub1a_b(n,gg))
     .                 *(msq1a_b(n,pntr(b(n),c(n)),j,k)
     .                  +msq1a_b(n,pntr(c(n),b(n)),j,k))*xn/2d0
     .                +subba_1v(n)/2d0
     .                 *(msqba_1v(n,pntr(b(n),c(n)),j,k)
     .                  +msqba_1v(n,pntr(c(n),b(n)),j,k))*xn/2d0
     .                +sub1a_bv(n)
     .                 *(msq1a_bv(n,pntr(b(n),c(n)),j,k)
     .                  +msq1a_bv(n,pntr(c(n),b(n)),j,k))*xn/2d0
          enddo
          do n=3,5,2
c-- (a,b,c) = (6,5,7) and (7,5,6)
          msq(6+n,j,k)=(subba_1(n,qq)+sub1a_b(n,gg))
     .                  *msq1a_b(n,pntr(b(n),c(n)),j,k)*xn/2d0
     .                 +sub1a_bv(n)
     .                  *msq1a_bv(n,pntr(b(n),c(n)),j,k)*xn/2d0
          enddo              
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(18+n,j,k)=msq(18+n,j,k)+sub1b_2(n,qg)*
     .      (msq1b_2(n,1,+1,k)+msq1b_2(n,1,+2,k)+msq1b_2(n,1,+3,k)
     .      +msq1b_2(n,1,+4,k)+msq1b_2(n,1,+5,k)
     .      +msq1b_2(n,2,+1,k)+msq1b_2(n,2,+2,k)+msq1b_2(n,2,+3,k)
     .      +msq1b_2(n,2,+4,k)+msq1b_2(n,2,+5,k))*xn*(aveqg/aveqq)
          enddo
c--- ANTIQUARK-GLUON contributions
      elseif ((j.lt.0) .and. (k.eq.0)) then
          do n=1,2  
c-- (a,b,c) = (5,6,7) and (5,7,6)
          msq(n,j,k)   =subab_c(n,qq)
     .                  *msqab_c(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(12+n,j,k)=(subbc_2(n,gg)/2d0+sub2c_b(n,gg))
     .                  *(msq2c_b(n,pntr(b(n),a(n)),j,k)
     .                   +msq2c_b(n,pntr(a(n),b(n)),j,k))*xn/2d0
     .                 +subbc_2v(n)/2d0
     .                  *(msqbc_2v(n,pntr(b(n),a(n)),j,k)
     .                   +msqbc_2v(n,pntr(a(n),b(n)),j,k))*xn/2d0
     .                 +sub2c_bv(n)
     .                  *(msq2c_bv(n,pntr(b(n),a(n)),j,k)
     .                   +msq2c_bv(n,pntr(a(n),b(n)),j,k))*xn/2d0
          msq(18+n,j,k)=sub1b_2(n,qq)
     .                  *msq1b_2(n,pntr(c(n),a(n)),j,k)*xn/2d0
          msq(21+n,j,k)=sub2b_1(n,gg)
     .                  *msq2b_1(n,pntr(c(n),a(n)),j,k)*xn/2d0
     .                 +sub2b_1v(n)
     .                  *msq2b_1v(n,pntr(c(n),a(n)),j,k)*xn/2d0
          enddo
          do n=4,6,2
c-- (a,b,c) = (6,7,5) and (7,6,5)
          msq(n,j,k)   =subab_c(n,gg)/2d0
     .                  *msqab_c(n,pntr(a(n),c(n)),j,k)*xn/2d0
     .                 +subab_cv(n)/2d0
     .                  *msqab_cv(n,pntr(a(n),c(n)),j,k)*xn/2d0
          msq(6+n,j,k) =(sub1a_b(n,qq)+subba_1(n,gg)/2d0)
     .                  *msq1a_b(n,pntr(c(n),b(n)),j,k)*xn/2d0
     .                 +subba_1v(n)/2d0
     .                  *msqba_1v(n,pntr(c(n),b(n)),j,k)*xn/2d0
          enddo
          do n=3,5,2
c-- (a,b,c) = (7,5,6) and (6,5,7)
          msq(12+n,j,k)=(subbc_2(n,qq)+sub2c_b(n,gg))
     .                  *msq2c_b(n,pntr(b(n),a(n)),j,k)*xn/2d0
     .                 +sub2c_bv(n)
     .                  *msq2c_bv(n,pntr(b(n),a(n)),j,k)*xn/2d0
          enddo      
c--- additional (qg) collinear contributions
c--- choose n=1 which is (7,5,6)
          do n=3,3
          msq(21+n,j,k)=msq(21+n,j,k)+sub2b_1(n,qg)*
     .      (msq2b_1(n,1,j,+1)+msq2b_1(n,1,j,+2)+msq2b_1(n,1,j,+3)
     .      +msq2b_1(n,1,j,+4)+msq2b_1(n,1,j,+5)
     .      +msq2b_1(n,2,j,+1)+msq2b_1(n,2,j,+2)+msq2b_1(n,2,j,+3)
     .      +msq2b_1(n,2,j,+4)+msq2b_1(n,2,j,+5))*xn*(aveqg/aveqq)
          enddo
      endif

      enddo
      enddo

c--debug
c
c      do j=-nf,nf
c      do k=-nf,nf
c      msumsq(j,k)=0d0
c      do nd=1,ndmax
c      msumsq(j,k)=msumsq(j,k)+msq(nd,j,k)
c      enddo
c      enddo
c      enddo
         if (Qflag) then
c---rename terms required for q-dipoles
            do ii=1,4
              sub57_6(ii)=subab_c(2,ii)
              sub67_5(ii)=subab_c(4,ii)

              sub17_5(ii)=sub1a_b(3,ii)
              sub17_6(ii)=sub1a_b(6,ii)
              sub57_1(ii)=subba_1(3,ii)
              sub67_1(ii)=subba_1(6,ii)
              sub27_6(ii)=sub2c_b(1,ii)
              sub27_5(ii)=sub2c_b(5,ii)
              sub67_2(ii)=subbc_2(1,ii)
              sub57_2(ii)=subbc_2(5,ii)

              sub16_2(ii)=sub1b_2(1,ii)
              sub17_2(ii)=sub1b_2(2,ii)
              sub15_2(ii)=sub1b_2(3,ii)

              sub26_1(ii)=sub2b_1(1,ii)
              sub27_1(ii)=sub2b_1(2,ii)
              sub25_1(ii)=sub2b_1(3,ii)

            enddo
         endif
      endif

      if (qflag) then
      if (.not. gflag) then
      do nd=1,ndmax
      do j=-nf,nf
      do k=-nf,nf
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo

C --if we have not called dips before we need to do it now
c--- final-final
      call dips(2,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m57_6,isq57_6)
      call dips(4,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m67_5,isq67_5)
c--- now the basic initial final and final initial
c--- second call for dipole 9,12 etc  only supplies new values for
c--- sub..
      call dips(9,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m17_5,isq17_5)
      call dips(9,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     . qqb_z2jet,donothing_gvec)

      call dips(12,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m17_6,isq17_6)
      call dips(12,p,6,7,1,sub67_1,dsubv,dummy,dummyv,
     . qqb_z2jet,donothing_gvec)

      call dips(13,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m27_6,isq27_6)
      call dips(13,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     . qqb_z2jet,donothing_gvec)

      call dips(17,p,2,7,5,sub27_5,dsubv,msq27_5,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m27_5,isq27_5)
      call dips(17,p,5,7,2,sub57_2,dsubv,dummy,dummyv,
     . qqb_z2jet,donothing_gvec)
c--- calculate all the initial-initial dipoles
      call dips(19,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m16_2,isq16_2)
      call dips(21,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m15_2,isq15_2)
      call dips(20,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m17_2,isq17_2)
      call dips(22,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m26_1,isq26_1)
      call dips(23,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m27_1,isq27_1)
      call dips(24,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     . qqb_z2jet,donothing_gvec)
      call storeqq(m25_1,isq25_1)
      endif

c      do ii=1,4
c      write(6,*)'sub57_6(ii)-subab_c(2,ii)',ii,sub57_6(ii)-subab_c(2,ii)
c      write(6,*)'sub67_5(ii)-subab_c(4,ii)',ii,sub67_5(ii)-subab_c(4,ii)
c
c      write(6,*)'sub17_5(ii)-sub1a_b(3,ii)',ii,sub17_5(ii)-sub1a_b(3,ii)
c      write(6,*)'sub17_6(ii)-sub1a_b(6,ii)',ii,sub17_6(ii)-sub1a_b(6,ii)
c      write(6,*)'sub57_1(ii)-subba_1(3,ii)',ii,sub57_1(ii)-subba_1(3,ii)
c      write(6,*)'sub67_1(ii)-subba_1(6,ii)',ii,sub67_1(ii)-subba_1(6,ii)
c
c      write(6,*)'sub27_6(ii)-sub2c_b(1,ii)',ii,sub27_6(ii)-sub2c_b(1,ii)
c      write(6,*)'sub27_5(ii)-sub2c_b(5,ii)',ii,sub27_5(ii)-sub2c_b(5,ii)
c      write(6,*)'sub67_2(ii)-subbc_2(1,ii)',ii,sub67_2(ii)-subbc_2(1,ii)
c      write(6,*)'sub57_2(ii)-subbc_2(5,ii)',ii,sub57_2(ii)-subbc_2(5,ii)
c
c      write(6,*)'sub16_2(ii)-sub1b_2(1,ii)',ii,sub16_2(ii)-sub1b_2(1,ii)
c      write(6,*)'sub17_2(ii)-sub1b_2(2,ii)',ii,sub17_2(ii)-sub1b_2(2,ii)
c      write(6,*)'sub15_2(ii)-sub1b_2(3,ii)',ii,sub15_2(ii)-sub1b_2(3,ii)

c      write(6,*)'sub26_1(ii)-sub2b_1(1,ii)',ii,sub26_1(ii)-sub2b_1(1,ii)
c      write(6,*)'sub27_1(ii)-sub2b_1(2,ii)',ii,sub27_1(ii)-sub2b_1(2,ii)
c      write(6,*)'sub25_1(ii)-sub2b_1(3,ii)',ii,sub25_1(ii)-sub2b_1(3,ii)
c      enddo



      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .gt. 0) .and. (k.gt.0)) then
c---QQ
      msq(2,j,k)=msq(2,j,k)
     . +sub57_6(qq)
     . *((xn+1d0/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     . +sub67_5(qq)
     . *((xn+1d0/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(9,j,k)=msq(9,j,k)
     . +sub17_5(qq)
     . *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     . +sub57_1(qq)
     . *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     . +sub17_6(qq)
     . *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     . +sub67_1(qq)
     . *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     . +sub27_6(qq)
     . *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     . +sub67_2(qq)
     . *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     . +sub27_5(qq)
     . *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     . +sub57_2(qq)
     . *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     . +sub17_2(qq)
     . *((xn+1d0/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     . +sub27_1(qq)
     . *((xn+1d0/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)

      elseif ((j .lt. 0) .and. (k.lt.0)) then
c---QbarQbar
      msq(2,j,k)=msq(2,j,k)
     . +sub57_6(qq)
     . *((xn+1d0/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(4,j,k)=msq(4,j,k)
     . +sub67_5(qq)
     . *((xn+1d0/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(9,j,k)=msq(9,j,k)
     . +sub17_5(qq)
     . *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     . +sub57_1(qq)
     . *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     . +sub17_6(qq)
     . *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     . +sub67_1(qq)
     . *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     . +sub27_6(qq)
     . *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     . +sub67_2(qq)
     . *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     . +sub27_5(qq)
     . *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     . +sub57_2(qq)
     . *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     . +sub17_2(qq)
     . *((xn+1d0/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(23,j,k)=msq(23,j,k)
     . +sub27_1(qq)
     . *((xn+1d0/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      elseif ((j .gt. 0) .and. (k.lt.0)) then
c---QQbar
      msq(2,j,k)=msq(2,j,k)
     . +sub57_6(qq)
     . *((xn-two/xn)*m57_6(2,j,k)-m57_6(0,j,k)/xn-m57_6(1,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     . +sub67_5(qq)
     . *((xn-two/xn)*m67_5(2,j,k)-m67_5(0,j,k)/xn-m67_5(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     . +sub17_5(qq)
     . *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
     . +sub57_1(qq)
     . *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
      msq(12,j,k)=msq(12,j,k)
     . +sub17_6(qq)
     . *((xn+1d0/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
     . +sub67_1(qq)
     . *((xn+1d0/xn)*m17_6(0,j,k)+two*(m17_6(1,j,k)+m17_6(2,j,k))/xn)
      msq(13,j,k)=msq(13,j,k)
     . +sub27_6(qq)
     . *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
     . +sub67_2(qq)
     . *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
      msq(17,j,k)=msq(17,j,k)
     . +sub27_5(qq)
     . *((xn+1d0/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
     . +sub57_2(qq)
     . *((xn+1d0/xn)*m27_5(0,j,k)+two*(m27_5(1,j,k)+m27_5(2,j,k))/xn)
      msq(20,j,k)=msq(20,j,k)
     . +sub17_2(qq)
     . *((xn-two/xn)*m17_2(2,j,k)-m17_2(0,j,k)/xn-m17_2(1,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     . +sub27_1(qq)
     . *((xn-two/xn)*m27_1(2,j,k)-m27_1(0,j,k)/xn-m27_1(1,j,k)/xn)
      elseif ((j .lt. 0) .and. (k.gt.0)) then
c---QbarQ
      msq(2,j,k)=msq(2,j,k)
     . +sub57_6(qq)
     . *((xn-two/xn)*m57_6(2,j,k)-m57_6(0,j,k)/xn-m57_6(1,j,k)/xn)
      msq(4,j,k)=msq(4,j,k)
     . +sub67_5(qq)
     . *((xn-two/xn)*m67_5(2,j,k)-m67_5(0,j,k)/xn-m67_5(1,j,k)/xn)
      msq(9,j,k)=msq(9,j,k)
     . +sub17_5(qq)
     . *((xn+1d0/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
     . +sub57_1(qq)
     . *((xn+1d0/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
      msq(12,j,k)=msq(12,j,k)
     . +sub17_6(qq)
     . *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     . +sub67_1(qq)
     . *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(13,j,k)=msq(13,j,k)
     . +sub27_6(qq)
     . *((xn+1d0/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
     . +sub67_2(qq)
     . *((xn+1d0/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
      msq(17,j,k)=msq(17,j,k)
     . +sub27_5(qq)
     . *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     . +sub57_2(qq)
     . *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(20,j,k)=msq(20,j,k)
     . +sub17_2(qq)
     . *((xn-two/xn)*m17_2(2,j,k)-m17_2(0,j,k)/xn-m17_2(1,j,k)/xn)
      msq(23,j,k)=msq(23,j,k)
     . +sub27_1(qq)
     . *((xn-two/xn)*m27_1(2,j,k)-m27_1(0,j,k)/xn-m27_1(1,j,k)/xn)

      elseif (j.eq.0) then
c---------g-q and g-qbar
        msq(20,j,k)=msq(20,j,k)+sub17_2(qg)*(msq17_2(-5,k)+msq17_2(+5,k)
     .  +msq17_2(-1,k)+msq17_2(-2,k)+msq17_2(-3,k)+msq17_2(-4,k)
     .  +msq17_2(+1,k)+msq17_2(+2,k)+msq17_2(+3,k)+msq17_2(+4,k))
      elseif (k.eq.0) then
c---------q-g and qbar-g
c        msq(23,j,k)=msq(23,j,k)+sub27_1(qg)*(msq27_1(j,-5)
C+msq27_1(j,+5)
c     .  +msq27_1(j,-1)+msq27_1(j,-2)+msq27_1(j,-3)+msq27_1(j,-4)
C     .  +msq27_1(j,+1)+msq27_1(j,+2)+msq27_1(j,+3)+msq27_1(j,+4)
c     . )
        msq(23,j,k)=msq(23,j,k)+sub27_1(qg)*(msq27_1(j,+1)+msq27_1(j,+2)
     . +msq27_1(j,+3)+msq27_1(j,+4)+msq27_1(j,+5))

        msq(22,j,k)=msq(22,j,k)
     . +sub26_1(qg)*(
     . +msq26_1(j,-1)+msq26_1(j,-2)
     . +msq26_1(j,-3)+msq26_1(j,-4)+msq26_1(j,-5)
     . -isq26_1(j,-1)-isq26_1(j,-2)
     . -isq26_1(j,-3)-isq26_1(j,-4)-isq26_1(j,-5))

       write(6,*) 'half*sub26_1(qg)*isq26_1(j,-1)',
     . half*sub26_1(qg)*isq26_1(j,-1)
       write(6,*) 'othyer 26 dutt',
     . +sub26_1(qg)*(
     . +msq26_1(j,-1)+msq26_1(j,-2)
     . +msq26_1(j,-3)+msq26_1(j,-4)+msq26_1(j,-5)
     . -isq26_1(j,-1)-isq26_1(j,-2)
     . -isq26_1(j,-3)-isq26_1(j,-4)-isq26_1(j,-5))


c       write(6,*) 'half*sub25_1(qg)*isq25_1(j,-1)',
c     . half*sub25_1(qg)*isq25_1(j,-1)
c       write(6,*) 'other 25 stuff',
c     . +sub25_1(qg)*(msq25_1(j,-1)-isq25_1(j,-1))

        msq(24,j,k)=msq(24,j,k)
     . +sub25_1(qg)*(msq25_1(j,-1)-isq25_1(j,-1))


c        msq(22,j,k)=half*isq25_1(j,-1)


c          write(6,*) 'sub27_1(qg)',sub26_1(qg)
c          write(6,*) 'j,msq27_1(j,2)',j,msq26_1(j,-2)
c          write(6,*) 'j,k,msq(22,j,k)',j,k,msq(22,j,k)
      endif


      enddo
      enddo

      endif 


      return
      end
      
      
