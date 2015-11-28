
      subroutine qqb_zbb_gs(P,msq)
************************************************************************
*     Author: J. M. Campbell                                           *
*     March, 2000.                                                     *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1) + qbar(-p2) -->  g(p5) + g(p6) + Z +g(p7)
c                                             |
c                                             --> e^-(p3)+e^+(p4)
c                            
c--all momenta incoming so some of the signs may look odd. 
c
c  This routine implements the 2-quark, 3-gluon subtractions and
c  is adapted from the gg -> Zbbg routine in Zbb
c  All dipoles are obtained from those by crossings
c  PLUS WE NEED ADDITIONAL DIPOLES TO INCLUDE NON-bbar FINAL STATES
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd
c --- remember: nd will count the dipoles
      
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),dot
      double precision 
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq17_5(-nf:nf,-nf:nf),msq17_5v(-nf:nf,-nf:nf),
     & msq17_6(-nf:nf,-nf:nf),msq17_6v(-nf:nf,-nf:nf),
     & msq27_5(-nf:nf,-nf:nf),msq27_5v(-nf:nf,-nf:nf),
     & msq27_6(-nf:nf,-nf:nf),msq27_6v(-nf:nf,-nf:nf),
     & msq57_6(-nf:nf,-nf:nf),msq67_5(-nf:nf,-nf:nf),
     & msq17_2v(-nf:nf,-nf:nf),msq27_1v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),
     & sub57_6(4),sub67_5(4),dsubv,
     & sub17_2v,sub27_1v,sub17_5v,sub17_6v,sub27_5v,sub27_6v
      double precision mgg17_2(0:2),mgg27_1(0:2),
     &                 mgg57_6(0:2),mgg67_5(0:2),
     &                 mgg17_5(0:2),mgg57_1(0:2),
     &                 mgg17_6(0:2),mgg67_1(0:2),
     &                 mgg27_5(0:2),mgg57_2(0:2),
     &                 mgg27_6(0:2),mgg67_2(0:2)
      double precision mgg17_2v(0:2),mgg27_1v(0:2),
     &                 mgg57_6v(0:2),mgg67_5v(0:2),
     &                 mgg17_5v(0:2),mgg57_1v(0:2),
     &                 mgg17_6v(0:2),mgg67_1v(0:2),
     &                 mgg27_5v(0:2),mgg57_2v(0:2),
     &                 mgg27_6v(0:2),mgg67_2v(0:2)
      double precision msq_cs(0:2,-nf:nf,-nf:nf),mmsq_cs(0:2,2,2)
      common/msq_cols/msq_cs,mmsq_cs

      external qqb_zbb,donothing_gvec
      external qqb_zbb_gvec,qqb_z2jet

      ndmax=8

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     . qqb_zbb,qqb_zbb_gvec)
      call storegg(mgg17_2,mgg17_2v)
      call dips(2,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     . qqb_zbb,qqb_zbb_gvec)
      call storegg(mgg27_1,mgg27_1v)

c--- final-final
      call dips(3,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     . qqb_zbb,donothing_gvec)
      call storegg(mgg57_6,mgg57_6v)
      call dips(4,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     . qqb_zbb,donothing_gvec)
      call storegg(mgg67_5,mgg67_5v)

c--- now the basic initial final ones
      call dips(5,p,1,7,5,sub17_5,sub17_5v,msq17_5,msq17_5v,
     . qqb_zbb,qqb_zbb_gvec)
      call storegg(mgg17_5,mgg17_5v)
      call dips(5,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     . qqb_zbb,donothing_gvec)
      call dips(6,p,1,7,6,sub17_6,sub17_6v,msq17_6,msq17_6v,
     . qqb_zbb,qqb_zbb_gvec)
      call storegg(mgg17_6,mgg17_6v)
      call dips(6,p,6,7,1,sub67_1,dsubv,dummy,dummyv,
     . qqb_zbb,donothing_gvec)
      call dips(7,p,2,7,5,sub27_5,sub27_5v,msq27_5,msq27_5v,
     . qqb_zbb,qqb_zbb_gvec)
      call storegg(mgg27_5,mgg27_5v)
      call dips(7,p,5,7,2,sub57_2,dsubv,dummy,dummyv,
     . qqb_zbb,donothing_gvec)
      call dips(8,p,2,7,6,sub27_6,sub27_6v,msq27_6,msq27_6v,
     . qqb_zbb,qqb_zbb_gvec)
      call storegg(mgg27_6,mgg27_6v)
      call dips(8,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     . qqb_zbb,donothing_gvec)
      
      do j=-nf,nf
      do k=-nf,nf
      
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo

      if ((j.gt.0).and.(k.lt.0)) then
c----------q-qb
          msq(1,j,k)=-sub17_2(qq)*mgg17_2(0)*(xn+one/xn)/2d0
     .               -sub17_2(qq)*(mgg17_2(1)+mgg17_2(2))/xn/2d0
          msq(2,j,k)=-sub27_1(qq)*mgg27_1(0)*(xn+one/xn)/2d0
     .               -sub27_1(qq)*(mgg27_1(1)+mgg27_1(2))/xn/2d0     
          msq(3,j,k)=(sub57_6(gg)*(mgg57_6(1)+mgg57_6(2))
     .               +sub57_6v*(mgg57_6v(1)+mgg57_6v(2)))*xn/2d0
          msq(4,j,k)=(sub67_5(gg)*(mgg67_5(1)+mgg67_5(2))
     .               +sub67_5v*(mgg67_5v(1)+mgg67_5v(2)))*xn/2d0
          msq(5,j,k)=((sub57_1(gg)+sub17_5(qq))*(mgg17_5(1)+mgg17_5(0))
     .                +sub57_1v*(mgg17_5v(1)+mgg17_5v(0)))*xn/2d0
          msq(6,j,k)=((sub67_1(gg)+sub17_6(qq))*(mgg17_6(2)+mgg17_6(0))
     .                +sub67_1v*(mgg17_6v(2)+mgg17_6v(0)))*xn/2d0
          msq(7,j,k)=((sub57_2(gg)+sub27_5(qq))*(mgg27_5(2)+mgg27_5(0))
     .                +sub57_2v*(mgg27_5v(2)+mgg27_5v(0)))*xn/2d0
          msq(8,j,k)=((sub67_2(gg)+sub27_6(qq))*(mgg27_6(1)+mgg27_6(0))
     .                +sub67_2v*(mgg27_6v(1)+mgg27_6v(0)))*xn/2d0

      elseif((j.lt.0).and.(k.gt.0)) then
c----------qb-q
          msq(1,j,k)=-sub17_2(qq)*mgg17_2(0)*(xn+one/xn)/2d0
     .               -sub17_2(qq)*(mgg17_2(1)+mgg17_2(2))/xn/2d0
          msq(2,j,k)=-sub27_1(qq)*mgg27_1(0)*(xn+one/xn)/2d0
     .               -sub27_1(qq)*(mgg27_1(1)+mgg27_1(2))/xn/2d0     
          msq(3,j,k)=(sub57_6(gg)*(mgg57_6(1)+mgg57_6(2))
     .               +sub57_6v*(mgg57_6v(1)+mgg57_6v(2)))*xn/2d0
          msq(4,j,k)=(sub67_5(gg)*(mgg67_5(1)+mgg67_5(2))
     .               +sub67_5v*(mgg67_5v(1)+mgg67_5v(2)))*xn/2d0
          msq(5,j,k)=((sub57_1(gg)+sub17_5(qq))*(mgg17_5(2)+mgg17_5(0))
     .                +sub57_1v*(mgg17_5v(2)+mgg17_5v(0)))*xn/2d0
          msq(6,j,k)=((sub67_1(gg)+sub17_6(qq))*(mgg17_6(1)+mgg17_6(0))
     .                +sub67_1v*(mgg17_6v(1)+mgg17_6v(0)))*xn/2d0
          msq(7,j,k)=((sub57_2(gg)+sub27_5(qq))*(mgg27_5(1)+mgg27_5(0))
     .                +sub57_2v*(mgg27_5v(1)+mgg27_5v(0)))*xn/2d0
          msq(8,j,k)=((sub67_2(gg)+sub27_6(qq))*(mgg27_6(2)+mgg27_6(0))
     .                +sub67_2v*(mgg27_6v(2)+mgg27_6v(0)))*xn/2d0

      elseif ((j.eq.0) .and. (k.eq.0)) then
c---------g-g
          msq(1,j,k)=(sub17_2(gg)*(mgg17_2(1)+mgg17_2(2))
     .               +sub17_2v*(mgg17_2v(1)+mgg17_2v(2)))*xn/2d0
          msq(2,j,k)=(sub27_1(gg)*(mgg27_1(1)+mgg27_1(2))
     .               +sub27_1v*(mgg27_1v(1)+mgg27_1v(2)))*xn/2d0
          msq(3,j,k)=-sub57_6(qq)*mgg57_6(0)*(xn+one/xn)/2d0
     .               -sub57_6(qq)*(mgg57_6(1)+mgg57_6(2))/xn/2d0
          msq(4,j,k)=-sub67_5(qq)*mgg67_5(0)*(xn+one/xn)/2d0
     .               -sub67_5(qq)*(mgg67_5(1)+mgg67_5(2))/xn/2d0
          msq(5,j,k)=((sub17_5(gg)+sub57_1(qq))*(mgg17_5(1)+mgg17_5(0))
     .                +sub17_5v*(mgg17_5v(1)+mgg17_5v(0)))*xn/2d0
          msq(6,j,k)=((sub17_6(gg)+sub67_1(qq))*(mgg17_6(2)+mgg17_6(0))
     .                +sub17_6v*(mgg17_6v(2)+mgg17_6v(0)))*xn/2d0
          msq(7,j,k)=((sub27_5(gg)+sub57_2(qq))*(mgg27_5(2)+mgg27_5(0))
     .                +sub27_5v*(mgg27_5v(2)+mgg27_5v(0)))*xn/2d0
          msq(8,j,k)=((sub27_6(gg)+sub67_2(qq))*(mgg27_6(1)+mgg27_6(0))
     .                +sub27_6v*(mgg27_6v(1)+mgg27_6v(0)))*xn/2d0

      elseif (j.eq.0) then
          if (k.gt.0) then
c---------g-q
          msq(1,j,k)=sub17_2(qg)
     &    *(msq17_2(-1,k)+msq17_2(-2,k)+msq17_2(-3,k)+msq17_2(-4,k)
     .     +msq17_2(-5,k))

          elseif (k.lt.0) then
c---------g-qbar
          msq(1,j,k)=sub17_2(qg)
     &    *(msq17_2(+1,k)+msq17_2(+2,k)+msq17_2(+3,k)+msq17_2(+4,k)
     .    +msq17_2(+5,k))
          endif 

      elseif (k.eq.0) then
          if (j.gt.0) then
c---------q-g
          msq(1,j,k)=(sub17_2(gg)*(mgg17_2(1)+mgg17_2(2))
     .               +sub17_2v*(mgg17_2v(1)+mgg17_2v(2)))*xn/2d0
          msq(2,j,k)=(sub27_1(gg)*(mgg27_1(1)+mgg27_1(2))
     .               +sub27_1v*(mgg27_1v(1)+mgg27_1v(2)))*xn/2d0
          msq(3,j,k)=-sub57_6(qq)*mgg57_6(0)*(xn+one/xn)/2d0
     .               -sub57_6(qq)*(mgg57_6(1)+mgg57_6(2))/xn/2d0
          msq(4,j,k)=-sub67_5(qq)*mgg67_5(0)*(xn+one/xn)/2d0
     .               -sub67_5(qq)*(mgg67_5(1)+mgg67_5(2))/xn/2d0
          msq(5,j,k)=((sub17_5(gg)+sub57_1(qq))*(mgg17_5(1)+mgg17_5(0))
     .                +sub17_5v*(mgg17_5v(1)+mgg17_5v(0)))*xn/2d0
          msq(6,j,k)=((sub17_6(gg)+sub67_1(qq))*(mgg17_6(2)+mgg17_6(0))
     .                +sub17_6v*(mgg17_6v(2)+mgg17_6v(0)))*xn/2d0
          msq(7,j,k)=((sub27_5(gg)+sub57_2(qq))*(mgg27_5(2)+mgg27_5(0))
     .                +sub27_5v*(mgg27_5v(2)+mgg27_5v(0)))*xn/2d0
          msq(8,j,k)=((sub27_6(gg)+sub67_2(qq))*(mgg27_6(1)+mgg27_6(0))
     .                +sub27_6v*(mgg27_6v(1)+mgg27_6v(0)))*xn/2d0

          msq(2,j,k)=sub27_1(qg)
     &    *(msq27_1(j,-1)+msq27_1(j,-2)+msq27_1(j,-3)+msq27_1(j,-4)
     .    +msq27_1(j,-5))

          elseif (j.lt.0) then
c---------qbar-g
          msq(2,j,k)=sub27_1(qg)
     &    *(msq27_1(j,+1)+msq27_1(j,+2)+msq27_1(j,+3)+msq27_1(j,+4)
     .     +msq27_1(j,+5))
          endif 
      endif

      enddo
      enddo

      return
      end

      
      subroutine storegg(mgg,mggv)
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration
      implicit none
      include 'constants.f'
      integer i
      double precision mgg(0:2),mggv(0:2)
      double precision msq_cs(0:2,-nf:nf,-nf:nf),mmsq_cs(0:2,2,2),
     .                msqv_cs(0:2,-nf:nf,-nf:nf),mmsqv_cs(0:2,2,2)
      common/msq_cols/msq_cs,mmsq_cs
      common/msqv_cols/msqv_cs,mmsqv_cs
      
      do i=0,2
        mgg(i)=msq_cs(i,0,0)
c--- need to add part for mggv here
        mggv(i)=msqv_cs(i,0,0)
      enddo
      
      return
      end
      

