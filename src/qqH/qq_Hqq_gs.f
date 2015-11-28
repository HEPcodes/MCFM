      subroutine qq_Hqq_gs(p,msq)
      implicit none 
************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2002.                                                      *
*                                                                      *
*     This routine calculates the dipole subtraction terms             *
*     for the process:                                                 *
*     q(-p1) + q(-p2) --> H(p34) + q(p5) + q(p6)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,m,n,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq17_5(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     . sub17_5(4),sub27_6(4),sub57_1(4),sub67_2(4),
     . msq15_7(-nf:nf,-nf:nf),msq26_7(-nf:nf,-nf:nf),
     . sub15_7(4),sub26_7(4),
     . msq16_7(-nf:nf,-nf:nf),msq25_7(-nf:nf,-nf:nf),
     . sub16_7(4),sub25_7(4),
     . msq15_7x(-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . msq26_7x(-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . msq16_7x(-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . msq25_7x(-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . msqx(-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),dsubv
      external qq_Hqq,donothing_gvec
      common/msq_all/msqx

      ndmax=6

c---- calculate the dipoles: initial-final and final-initial
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     . qq_Hqq,donothing_gvec)
      call dips(1,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     . qq_Hqq,donothing_gvec)
      call dips(2,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     . qq_Hqq,donothing_gvec)
      call dips(2,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     . qq_Hqq,donothing_gvec)

      call dips(3,p,1,5,7,sub15_7,dsubv,msq15_7,dummyv,
     . qq_Hqq,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
        msq15_7x(j,k,m,n)=msqx(j,k,m,n)
      enddo
      enddo
      enddo
      enddo

      call dips(4,p,2,6,7,sub26_7,dsubv,msq26_7,dummyv,
     . qq_Hqq,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
        msq26_7x(j,k,m,n)=msqx(j,k,m,n)
      enddo
      enddo
      enddo
      enddo

      call dips(5,p,1,6,7,sub16_7,dsubv,msq16_7,dummyv,
     . qq_Hqq,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
        msq16_7x(j,k,m,n)=msqx(j,k,m,n)
      enddo
      enddo
      enddo
      enddo

      call dips(6,p,2,5,7,sub25_7,dsubv,msq25_7,dummyv,
     . qq_Hqq,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
        msq25_7x(j,k,m,n)=msqx(j,k,m,n)
      enddo
      enddo
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if  ((j .eq. 0) .and. (k .eq. 0)) then
         goto 20
      elseif ((j .ne. 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*cf*(sub17_5(qq)+sub57_1(qq))*msq17_5(j,k)
         msq(2,j,k)=2d0*cf*(sub27_6(qq)+sub67_2(qq))*msq27_6(j,k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub27_6(qg)*(
     .              +msq27_6(j,+1)+msq27_6(j,+2)+msq27_6(j,+3)
     .              +msq27_6(j,+4)+msq27_6(j,+5))
         if     (j .eq. 1) then
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(1,-1,1,-1)/2d0+msq26_7x(1,-1,2,-2)+msq26_7x(1,-1,3,-3)
     . +msq26_7x(1,-1,4,-4)+msq26_7x(1,-1,5,-5)+msq26_7x(1,-2,3,-4))
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(1,-1,1,-1)/2d0+msq25_7x(1,-2,1,-2)+msq25_7x(1,-3,1,-3)
     . +msq25_7x(1,-4,1,-4)+msq25_7x(1,-5,1,-5)+msq25_7x(1,-3,2,-4))
         elseif (j .eq. 2) then
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(2,-1,2,-1)+msq26_7x(2,-2,2,-2)/2d0+msq26_7x(2,-2,3,-3)
     . +msq26_7x(2,-1,4,-4)+msq26_7x(2,-2,4,-4)+msq26_7x(2,-2,5,-5))
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(2,-2,1,-1)+msq25_7x(2,-2,2,-2)/2d0+msq25_7x(2,-3,2,-3)
     . +msq25_7x(2,-4,1,-3)+msq25_7x(2,-4,2,-4)+msq25_7x(2,-5,2,-5))
         elseif (j .eq. 3) then
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(3,-1,3,-1)+msq26_7x(3,-2,3,-2)+msq26_7x(3,-3,3,-3)/2d0
     . +msq26_7x(3,-1,4,-2)+msq26_7x(3,-3,4,-4)+msq26_7x(3,-3,5,-5))
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(3,-3,1,-1)+msq25_7x(3,-3,2,-2)+msq25_7x(3,-3,3,-3)/2d0
     . +msq25_7x(3,-4,1,-2)+msq25_7x(3,-4,3,-4)+msq25_7x(3,-5,3,-5))
         elseif (j .eq. 4) then
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(4,-2,3,-1)+msq26_7x(4,-1,4,-1)+msq26_7x(4,-2,4,-2)
     . +msq26_7x(4,-3,4,-3)+msq26_7x(4,-4,4,-4)/2d0+msq26_7x(4,-4,5,-5))
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(4,-3,2,-1)+msq25_7x(4,-4,1,-1)+msq25_7x(4,-4,2,-2)
     . +msq25_7x(4,-4,3,-3)+msq25_7x(4,-4,4,-4)/2d0+msq25_7x(4,-5,4,-5))
         elseif (j .eq. 5) then
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(5,-1,5,-1)+msq26_7x(5,-2,5,-2)+msq26_7x(5,-3,5,-3)
     . +msq26_7x(5,-4,5,-4)+msq26_7x(5,-5,5,-5)/2d0)
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(5,-5,1,-1)+msq25_7x(5,-5,2,-2)/2d0+msq25_7x(5,-5,3,-3)
     . +msq25_7x(5,-5,4,-4)+msq25_7x(5,-5,5,-5)/2d0)
         endif
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub27_6(qg)*(
     .              +msq27_6(j,-5)+msq27_6(j,-4)+msq27_6(j,-3)
     .              +msq27_6(j,-2)+msq27_6(j,-1))
         if     (j .eq. -1) then
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(-1,1,-1,1)/2d0+msq25_7x(-1,1,-2,2)+msq25_7x(-1,1,-3,3)
     . +msq25_7x(-1,1,-4,4)+msq25_7x(-1,1,-5,5)+msq25_7x(-1,2,-3,4))
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(-1,1,-1,1)/2d0+msq26_7x(-1,2,-1,2)+msq26_7x(-1,3,-1,3)
     . +msq26_7x(-1,4,-1,4)+msq26_7x(-1,5,-1,5)+msq26_7x(-1,3,-2,4))
         elseif (j .eq. -2) then
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(-2,1,-2,1)+msq25_7x(-2,2,-2,2)/2d0+msq25_7x(-2,2,-3,3)
     . +msq25_7x(-2,1,-4,4)+msq25_7x(-2,2,-4,4)+msq25_7x(-2,2,-5,5))
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(-2,2,-1,1)+msq26_7x(-2,2,-2,2)/2d0+msq26_7x(-2,3,-2,3)
     . +msq26_7x(-2,4,-1,3)+msq26_7x(-2,4,-2,4)+msq26_7x(-2,5,-2,5))
         elseif (j .eq. -3) then
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(-3,1,-3,1)+msq25_7x(-3,2,-3,2)+msq25_7x(-3,3,-3,3)/2d0
     . +msq25_7x(-3,1,-4,2)+msq25_7x(-3,3,-4,4)+msq25_7x(-3,3,-5,5))
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(-3,3,-1,1)+msq26_7x(-3,3,-2,2)+msq26_7x(-3,3,-3,3)/2d0
     . +msq26_7x(-3,4,-1,2)+msq26_7x(-3,4,-3,4)+msq26_7x(-3,5,-3,5))
         elseif (j .eq. -4) then
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(-4,2,-3,1)+msq25_7x(-4,1,-4,1)+msq25_7x(-4,2,-4,2)
     . +msq25_7x(-4,3,-4,3)+msq25_7x(-4,4,-4,4)/2d0+msq25_7x(-4,4,-5,5))
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(-4,3,-2,1)+msq26_7x(-4,4,-1,1)+msq26_7x(-4,4,-2,2)
     . +msq26_7x(-4,4,-3,3)+msq26_7x(-4,4,-4,4)/2d0+msq26_7x(-4,5,-4,5))
         elseif (j .eq. -5) then
         msq(6,j,k)=2d0*tr*sub25_7(qg)*(
     . +msq25_7x(-5,1,-5,1)+msq25_7x(-5,2,-5,2)+msq25_7x(-5,3,-5,3)
     . +msq25_7x(-5,4,-5,4)+msq25_7x(-5,5,-5,5)/2d0)
         msq(4,j,k)=2d0*tr*sub26_7(qg)*(
     . +msq26_7x(-5,5,-1,1)+msq26_7x(-5,5,-2,2)/2d0+msq26_7x(-5,5,-3,3)
     . +msq26_7x(-5,5,-4,4)+msq26_7x(-5,5,-5,5)/2d0)
         endif
       elseif ((j .eq. 0) .and. (k .gt. 0)) then
         msq(1,j,k)=2d0*tr*sub17_5(qg)*(
     .              +msq17_5(+5,k)+msq17_5(+4,k)+msq17_5(+3,k)
     .              +msq17_5(+2,k)+msq17_5(+1,k))
         if     (k .eq. 1) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(-1,1,1,-1)/2d0+msq16_7x(-1,1,2,-2)+msq16_7x(-1,1,3,-3)
     . +msq16_7x(-1,1,4,-4)+msq16_7x(-1,1,5,-5)+msq16_7x(-2,1,3,-4))
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(-1,1,1,-1)/2d0+msq15_7x(-2,1,1,-2)+msq15_7x(-3,1,1,-3)
     . +msq15_7x(-4,1,1,-4)+msq15_7x(-5,1,1,-5)+msq15_7x(-3,1,2,-4))
         elseif (k .eq. 2) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(-1,2,2,-1)+msq16_7x(-2,2,2,-2)/2d0+msq16_7x(-2,2,3,-3)
     . +msq16_7x(-1,2,4,-4)+msq16_7x(-2,2,4,-4)+msq16_7x(-2,2,5,-5))
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(-2,2,1,-1)+msq15_7x(-2,2,2,-2)/2d0+msq15_7x(-3,2,2,-3)
     . +msq15_7x(-4,2,1,-3)+msq15_7x(-4,2,2,-4)+msq15_7x(-5,2,2,-5))
         elseif (k .eq. 3) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(-1,3,3,-1)+msq16_7x(-2,3,3,-2)+msq16_7x(-3,3,3,-3)/2d0
     . +msq16_7x(-1,3,4,-2)+msq16_7x(-3,3,4,-4)+msq16_7x(-3,3,5,-5))
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(-3,3,1,-1)+msq15_7x(-3,3,2,-2)+msq15_7x(-3,3,3,-3)/2d0
     . +msq15_7x(-4,3,1,-2)+msq15_7x(-4,3,3,-4)+msq15_7x(-5,3,3,-5))
         elseif (k .eq. 4) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(-2,4,3,-1)+msq16_7x(-1,4,4,-1)+msq16_7x(-2,4,4,-2)
     . +msq16_7x(-3,4,4,-3)+msq16_7x(-4,4,4,-4)/2d0+msq16_7x(-4,4,5,-5))
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(-3,4,2,-1)+msq15_7x(-4,4,1,-1)+msq15_7x(-4,4,2,-2)
     . +msq15_7x(-4,4,3,-3)+msq15_7x(-4,4,4,-4)/2d0+msq15_7x(-5,4,4,-5))
         elseif (k .eq. 5) then
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(-1,5,5,-1)+msq16_7x(-2,5,5,-2)+msq16_7x(-3,5,5,-3)
     . +msq16_7x(-4,5,5,-4)+msq16_7x(-5,5,5,-5)/2d0)
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(-5,5,1,-1)+msq15_7x(-5,5,2,-2)/2d0+msq15_7x(-5,5,3,-3)
     . +msq15_7x(-5,5,4,-4)+msq15_7x(-5,5,5,-5)/2d0)
         endif
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
         msq(1,j,k)=2d0*tr*sub17_5(qg)*(
     .              +msq17_5(-5,k)+msq17_5(-4,k)+msq17_5(-3,k)
     .              +msq17_5(-2,k)+msq17_5(-1,k))
         if     (k .eq. -1) then
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(1,-1,-1,1)/2d0+msq15_7x(1,-1,-2,2)+msq15_7x(1,-1,-3,3)
     . +msq15_7x(1,-1,-4,4)+msq15_7x(1,-1,-5,5)+msq15_7x(2,-1,-3,4))
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(1,-1,-1,1)/2d0+msq16_7x(2,-1,-1,2)+msq16_7x(3,-1,-1,3)
     . +msq16_7x(4,-1,-1,4)+msq16_7x(5,-1,-1,5)+msq16_7x(3,-1,-2,4))
         elseif (k .eq. -2) then
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(1,-2,-2,1)+msq15_7x(2,-2,-2,2)/2d0+msq15_7x(2,-2,-3,3)
     . +msq15_7x(1,-2,-4,4)+msq15_7x(2,-2,-4,4)+msq15_7x(2,-2,-5,5))
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(2,-2,-1,1)+msq16_7x(2,-2,-2,2)/2d0+msq16_7x(3,-2,-2,3)
     . +msq16_7x(4,-2,-1,3)+msq16_7x(4,-2,-2,4)+msq16_7x(5,-2,-2,5))
         elseif (k .eq. -3) then
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(1,-3,-3,1)+msq15_7x(2,-3,-3,2)+msq15_7x(3,-3,-3,3)/2d0
     . +msq15_7x(1,-3,-4,2)+msq15_7x(3,-3,-4,4)+msq15_7x(3,-3,-5,5))
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(3,-3,-1,1)+msq16_7x(3,-3,-2,2)+msq16_7x(3,-3,-3,3)/2d0
     . +msq16_7x(4,-3,-1,2)+msq16_7x(4,-3,-3,4)+msq16_7x(5,-3,-3,5))
         elseif (k .eq. -4) then
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(2,-4,-3,1)+msq15_7x(1,-4,-4,1)+msq15_7x(2,-4,-4,2)
     . +msq15_7x(3,-4,-4,3)+msq15_7x(4,-4,-4,4)/2d0+msq15_7x(4,-4,-5,5))
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(3,-4,-2,1)+msq16_7x(4,-4,-1,1)+msq16_7x(4,-4,-2,2)
     . +msq16_7x(4,-4,-3,3)+msq16_7x(4,-4,-4,4)/2d0+msq16_7x(5,-4,-4,5))
         elseif (k .eq. -5) then
         msq(3,j,k)=2d0*tr*sub15_7(qg)*(
     . +msq15_7x(1,-5,-5,1)+msq15_7x(2,-5,-5,2)+msq15_7x(3,-5,-5,3)
     . +msq15_7x(4,-5,-5,4)+msq15_7x(5,-5,-5,5)/2d0)
         msq(5,j,k)=2d0*tr*sub16_7(qg)*(
     . +msq16_7x(5,-5,-1,1)+msq16_7x(5,-5,-2,2)/2d0+msq16_7x(5,-5,-3,3)
     . +msq16_7x(5,-5,-4,4)+msq16_7x(5,-5,-5,5)/2d0)
         endif
      endif
 20   continue
      enddo
      enddo

      return      
      end


