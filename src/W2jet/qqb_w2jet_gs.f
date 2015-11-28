      subroutine qqb_w2jet_gs(p,msq)
************************************************************************
*     Author: J.M. Campbell                                            *
*     August, 1999.                                                    *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W + parton(p5) + parton(p6) + parton(p7)
c                           |
c                            --> l(p3)+a(p4)
c   positively charged W only

      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,mj,mk
c --- remember: nd will count the dipoles
      integer nd,nu,fn
      parameter (fn=-nf)
c--- slightly obtuse notation, to simplify declaration lines      
      double precision p(mxpart,4),ptrans(mxpart,4),
     & msq(maxd,-nf:nf,-nf:nf),x,omx,y,u,z,gamqq,gamqg,
     & facq,facg,color
      double precision Q(mxpart,4),s(mxpart,mxpart)
      double precision
     & sub15_2(4),sub15_2v,msq15_2(fn:nf,fn:nf),msq15_2v(fn:nf,fn:nf),
     & sub25_1(4),sub25_1v,msq25_1(fn:nf,fn:nf),msq25_1v(fn:nf,fn:nf),
     & sub16_2(4),sub16_2v,msq16_2(fn:nf,fn:nf),msq16_2v(fn:nf,fn:nf),
     & sub26_1(4),sub26_1v,msq26_1(fn:nf,fn:nf),msq26_1v(fn:nf,fn:nf),
     & sub17_2(4),sub17_2v,msq17_2(fn:nf,fn:nf),msq17_2v(fn:nf,fn:nf),
     & sub27_1(4),sub27_1v,msq27_1(fn:nf,fn:nf),msq27_1v(fn:nf,fn:nf),
     & sub15_6(4),sub15_6v,msq15_6(fn:nf,fn:nf),msq15_6v(fn:nf,fn:nf),
     & sub15_7(4),sub15_7v,msq15_7(fn:nf,fn:nf),msq15_7v(fn:nf,fn:nf),
     & sub16_7(4),sub16_7v,msq16_7(fn:nf,fn:nf),msq16_7v(fn:nf,fn:nf),
     & sub25_6(4),sub25_6v,msq25_6(fn:nf,fn:nf),msq25_6v(fn:nf,fn:nf),
     & sub25_7(4),sub25_7v,msq25_7(fn:nf,fn:nf),msq25_7v(fn:nf,fn:nf),
     & sub26_7(4),sub26_7v,msq26_7(fn:nf,fn:nf),msq26_7v(fn:nf,fn:nf),
     & sub16_5(4),sub16_5v,sub17_5(4),sub17_5v,sub17_6(4),sub17_6v,
     & sub26_5(4),sub26_5v,sub27_5(4),sub27_5v,sub27_6(4),sub27_6v,
     & sub56_1(4),sub56_1v,sub57_1(4),sub57_1v,sub67_1(4),sub67_1v,
     & msq56_1v(fn:nf,fn:nf),msq56_2v(fn:nf,fn:nf),
     & msq57_1v(fn:nf,fn:nf),msq57_2v(fn:nf,fn:nf),
     & msq67_1v(fn:nf,fn:nf),msq67_2v(fn:nf,fn:nf),
     & sub56_2(4),sub56_2v,sub57_2(4),sub57_2v,sub67_2(4),sub67_2v,
     & sub56_7(4),sub56_7v,msq56_7(fn:nf,fn:nf),msq56_7v(fn:nf,fn:nf),
     & sub57_6(4),sub57_6v,msq57_6(fn:nf,fn:nf),msq57_6v(fn:nf,fn:nf),
     & sub67_5(4),sub67_5v,msq67_5(fn:nf,fn:nf),msq67_5v(fn:nf,fn:nf),
     & dsub,dsubv,dummy(fn:nf,fn:nf),dummyv(fn:nf,fn:nf)
      external qqb_w2jet,qqb_w2jet_gvec

      ndmax=15

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(4,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(5,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(6,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)
c--- now the initial final and final initial ones
      call dips(7,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(7,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(7,p,1,6,5,sub16_5,sub16_5v,dummy,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)

      call dips(8,p,1,5,7,sub15_7,sub15_7v,msq15_7,msq15_7v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(8,p,5,7,1,sub57_1,sub57_1v,dummy,msq57_1v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(8,p,1,7,5,sub17_5,sub17_5v,dummy,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)


      call dips(9,p,1,6,7,sub16_7,sub16_7v,msq16_7,msq16_7v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(9,p,6,7,1,sub67_1,sub67_1v,dummy,msq67_1v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(9,p,1,7,6,sub17_6,sub17_6v,dummy,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)

      call dips(10,p,2,5,6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(10,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(10,p,2,6,5,sub26_5,sub26_5v,dummy,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)

      call dips(11,p,2,5,7,sub25_7,sub25_7v,msq25_7,msq25_7v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(11,p,2,7,5,sub27_5,sub27_5v,dummy,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(11,p,5,7,2,sub57_2,sub57_2v,dummy,msq57_2v,
     . qqb_w2jet,qqb_w2jet_gvec)

      call dips(12,p,2,7,6,sub27_6,sub27_6v,dummy,dummyv,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(12,p,6,7,2,sub67_2,sub67_2v,dummy,msq67_2v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(12,p,2,6,7,sub26_7,sub26_7v,msq26_7,msq26_7v,
     . qqb_w2jet,qqb_w2jet_gvec)



c--- final-final dipoles
      call dips(13,p,5,6,7,sub56_7,sub56_7v,msq56_7,msq56_7v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(14,p,5,7,6,sub57_6,sub57_6v,msq57_6,msq57_6v,
     . qqb_w2jet,qqb_w2jet_gvec)
      call dips(15,p,6,7,5,sub67_5,sub67_5v,msq67_5,msq67_5v,
     . qqb_w2jet,qqb_w2jet_gvec)
           
      do j=-nf,nf
      do k=-nf,nf
      
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo

c--- do only q-qb and qb-q cases      
      if (  ((j .gt. 0).and.(k .lt. 0))
     . .or. ((j .lt. 0).and.(k .gt. 0))) then
      msq(1,j,k)=msq15_2(j,k)*sub15_2(qq)/xn**4
      msq(2,j,k)=msq25_1(j,k)*sub25_1(qq)/xn**4
      msq(3,j,k)=msq16_2(j,k)*sub16_2(qq)/xn**4
      msq(4,j,k)=msq26_1(j,k)*sub26_1(qq)/xn**4
      msq(5,j,k)=msq17_2(j,k)*sub17_2(qq)/xn**4
      msq(6,j,k)=msq27_1(j,k)*sub27_1(qq)/xn**4
      msq(7,j,k)=(one-ninth)*
     & (msq15_6(j,k)*(sub15_6(qq)+sub56_1(gg))+msq56_1v(j,k)*sub56_1v
     & +msq15_6(j,k)*(sub16_5(qq)+sub56_1(gg))+msq56_1v(j,k)*sub56_1v)
      msq(8,j,k)=(one-ninth)*
     & (msq15_7(j,k)*(sub15_7(qq)+sub57_1(gg))+msq57_1v(j,k)*sub57_1v
     & +msq15_7(j,k)*(sub17_5(qq)+sub57_1(gg))+msq57_1v(j,k)*sub57_1v)
      msq(9,j,k)=(one-ninth)*
     & (msq16_7(j,k)*(sub16_7(qq)+sub67_1(gg))+msq67_1v(j,k)*sub67_1v
     & +msq16_7(j,k)*(sub17_6(qq)+sub67_1(gg))+msq67_1v(j,k)*sub67_1v)
      msq(10,j,k)=(one-ninth)*
     & (msq25_6(j,k)*(sub25_6(qq)+sub56_2(gg))+msq56_2v(j,k)*sub56_2v
     & +msq25_6(j,k)*(sub26_5(qq)+sub56_2(gg))+msq56_2v(j,k)*sub56_2v)
      msq(11,j,k)=(one-ninth)*
     & (msq25_7(j,k)*(sub25_7(qq)+sub57_2(gg))+msq57_2v(j,k)*sub57_2v
     & +msq25_7(j,k)*(sub27_5(qq)+sub57_2(gg))+msq57_2v(j,k)*sub57_2v)
      msq(12,j,k)=(one-ninth)*
     & (msq26_7(j,k)*(sub26_7(qq)+sub67_2(gg))+msq67_2v(j,k)*sub67_2v
     & +msq26_7(j,k)*(sub27_6(qq)+sub67_2(gg))+msq67_2v(j,k)*sub67_2v)
      msq(13,j,k)=msq56_7(j,k)*sub56_7(gg)+msq56_7v(j,k)*sub56_7v
      msq(14,j,k)=msq57_6(j,k)*sub57_6(gg)+msq57_6v(j,k)*sub57_6v
      msq(15,j,k)=msq67_5(j,k)*sub67_5(gg)+msq67_5v(j,k)*sub67_5v
      elseif ((k .eq. 0).and.((j.gt.0) .or. (j.lt.0))) then
c--- q-g and qb-g cases
      elseif ((j .eq. 0).and.((k.gt.0) .or. (k.lt.0))) then
c--- g-q and g-qb cases
      elseif ((j .eq. 0).and.(k .eq. 0)) then
c--- g-g case
      endif
      
      enddo
      enddo

      return
      end
      
