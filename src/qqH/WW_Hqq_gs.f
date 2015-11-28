      subroutine WW_Hqq_gs(p,msq)
      implicit none 
************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2002.                                                      *
*                                                                      *
*     Weak Boson Fusion by W-W exchange only                           *
*     This routine calculates the dipole subtraction terms             *
*     for the process:                                                 *
*     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,m,n,msgn,nsgn,nd,pn(-nf:nf)

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq17_5(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     . msq17_6(-nf:nf,-nf:nf),msq27_5(-nf:nf,-nf:nf),
     . sub17_5(4),sub27_6(4),sub57_1(4),sub67_2(4),
     . sub17_6(4),sub27_5(4),sub67_1(4),sub57_2(4),
     . msq15_7(-nf:nf,-nf:nf),msq26_7(-nf:nf,-nf:nf),
     . sub15_7(4),sub26_7(4),
     . msq16_7(-nf:nf,-nf:nf),
     . sub16_7(4),
     . dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),dsubv,stat
      external WW_Hqq,donothing_gvec

      data pn/-1,-2,-1,-2,-1,0,1,2,1,2,1/

      ndmax=5

c---- calculate the dipoles: initial-final and final-initial
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     . WW_Hqq,donothing_gvec)
      call dips(1,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     . WW_Hqq,donothing_gvec)

      call dips(2,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     . WW_Hqq,donothing_gvec)
      call dips(2,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     . WW_Hqq,donothing_gvec)

      call dips(3,p,1,5,7,sub15_7,dsubv,msq15_7,dummyv,
     . WW_Hqq,donothing_gvec)
      call dips(4,p,2,6,7,sub26_7,dsubv,msq26_7,dummyv,
     . WW_Hqq,donothing_gvec)
      call dips(5,p,1,6,7,sub16_7,dsubv,msq16_7,dummyv,
     . WW_Hqq,donothing_gvec)
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if  ((j .eq. 0) .and. (k .eq. 0)) then
         goto 20
      elseif ((j .ne. 0) .and. (k .ne. 0)) then
         if     ((j .gt. 0) .and. (k .lt. 0)) then
c--- q-qb case
         msq(1,j,k)=2d0*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2d0*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
         elseif ((j .lt. 0) .and. (k .gt. 0)) then
c--- qb-q case
         msq(1,j,k)=2d0*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2d0*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
         else
c--- q-q and qb-qb cases
         msq(1,j,k)=2d0*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2d0*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
         endif
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub27_6(qg)*(
     .              +msq27_6(j,+1)+msq27_6(j,+2)+msq27_6(j,+3)
     .              +msq27_6(j,+4)+msq27_6(j,+5))
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     .              +msq26_7(j,-1)+msq26_7(j,-2)+msq26_7(j,-3)
     .              +msq26_7(j,-4)+msq26_7(j,-5))
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub27_6(qg)*(
     .              +msq27_6(j,-5)+msq27_6(j,-4)+msq27_6(j,-3)
     .              +msq27_6(j,-2)+msq27_6(j,-1))
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     .              +msq26_7(j,+1)+msq26_7(j,+2)+msq26_7(j,+3)
     .              +msq26_7(j,+4)+msq26_7(j,+5))
       elseif ((j .eq. 0) .and. (k .gt. 0)) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     .              +msq16_7(+5,k)+msq16_7(+4,k)+msq16_7(+3,k)
     .              +msq16_7(+2,k)+msq16_7(+1,k))
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     .              +msq15_7(-1,k)+msq15_7(-2,k)+msq15_7(-3,k)
     .              +msq15_7(-4,k)+msq15_7(-5,k))
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     .              +msq16_7(-5,k)+msq16_7(-4,k)+msq16_7(-3,k)
     .              +msq16_7(-2,k)+msq16_7(-1,k))
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     .              +msq15_7(+1,k)+msq15_7(+2,k)+msq15_7(+3,k)
     .              +msq15_7(+4,k)+msq15_7(+5,k))
      endif
 20   continue
      enddo
      enddo

      return      
      end


