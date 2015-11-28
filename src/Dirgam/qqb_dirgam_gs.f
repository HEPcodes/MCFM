      subroutine qqb_dirgam_gs(p,msq)
************************************************************************
*     Author: J.M. Campbell                                            *
*     October, 2002.                                                   *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
*     f(-p1) + f(-p2) -->  gamma(p3) + parton(p4) + parton(p5)         *
************************************************************************

      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd
c --- remember: nd will count the dipoles
      
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision 
     & msq14_2(-nf:nf,-nf:nf),msq24_1(-nf:nf,-nf:nf),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq14_5(-nf:nf,-nf:nf),msq25_4(-nf:nf,-nf:nf),
     & msq45_1v(-nf:nf,-nf:nf),msq45_2v(-nf:nf,-nf:nf),
     & msq25_4v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     & msq14_5v(-nf:nf,-nf:nf),msq15_2v(-nf:nf,-nf:nf),
     & msq24_1v(-nf:nf,-nf:nf),
     & msq14_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),dummyv(-nf:nf,-nf:nf),
     & sub14_2(4),sub24_1(4),sub15_2(4),sub25_1(4),
     & sub14_5(4),sub15_4(4),sub25_4(4),sub24_5(4),
     & sub45_1(4),sub45_2(4),sub45_1v,sub45_2v,
     & sub25_4v,sub24_1v,sub25_1v,sub15_4v,sub15_2v,sub14_2v,sub14_5v,
     & sub24_5v
      external qqb_dirgam,qqb_dirgam_gvec,donothing_gvec
      ndmax=6

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,4,2,sub14_2,sub14_2v,msq14_2,msq14_2v,
     . qqb_dirgam,qqb_dirgam_gvec)
      call dips(2,p,2,4,1,sub24_1,sub24_1v,msq24_1,msq24_1v,
     . qqb_dirgam,qqb_dirgam_gvec)
      call dips(3,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     . qqb_dirgam,qqb_dirgam_gvec)
      call dips(4,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     . qqb_dirgam,qqb_dirgam_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,4,5,sub14_5,sub14_5v,msq14_5,msq14_5v,
     . qqb_dirgam,qqb_dirgam_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,4,5,1,sub45_1,sub45_1v,dummy,msq45_1v,
     . qqb_dirgam,qqb_dirgam_gvec)
      call dips(5,p,1,5,4,sub15_4,sub15_4v,dummy,dummyv,
     . qqb_dirgam,donothing_gvec)

      call dips(6,p,2,5,4,sub25_4,sub25_4v,msq25_4,msq25_4v,
     . qqb_dirgam,qqb_dirgam_gvec)
      call dips(6,p,4,5,2,sub45_2,sub45_2v,dummy,msq45_2v,
     . qqb_dirgam,qqb_dirgam_gvec)
      call dips(6,p,2,4,5,sub24_5,sub24_5v,dummy,dummyv,
     . qqb_dirgam,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf      
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo
      
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

c--- do only q-qb and qb-q cases      
      if (  ((j .gt. 0).and.(k .lt. 0))
     . .or. ((j .lt. 0).and.(k .gt. 0))) then
C-----half=statistical factor
      msq(1,j,k)=-half*msq14_2(j,k)*sub14_2(qq)/xn
      msq(2,j,k)=-half*msq24_1(j,k)*sub24_1(qq)/xn
      msq(3,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
      msq(4,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
      msq(5,j,k)=half*xn*(
     .  msq14_5(j,k)*(sub14_5(qq)+0.5d0*sub45_1(gg))
     . +0.5d0*msq45_1v(j,k)*sub45_1v
     . +msq14_5(j,k)*(sub15_4(qq)+0.5d0*sub45_1(gg))
     . +0.5d0*msq45_1v(j,k)*sub45_1v)
      msq(6,j,k)=half*xn*(
     .  msq25_4(j,k)*(sub25_4(qq)+0.5d0*sub45_2(gg))
     . +0.5d0*msq45_2v(j,k)*sub45_2v
     . +msq25_4(j,k)*(sub24_5(qq)+0.5d0*sub45_2(gg))
     . +0.5d0*msq45_2v(j,k)*sub45_2v)
      elseif ((k .eq. 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2d0*tr*msq24_1(j,-j)*sub24_1(qg)
      msq(3,j,k)=xn*msq15_2(j,k)*sub15_2(qq)
      msq(4,j,k)=xn*(msq25_1(j,k)*sub25_1(gg)+msq25_1v(j,k)*sub25_1v)
      msq(5,j,k)=-(msq14_5(j,k)*sub15_4(qq)+msq14_5(j,k)*sub45_1(qq))/xn
      msq(6,j,k)=xn*(msq25_4(j,k)*sub25_4(gg)+msq25_4v(j,k)*sub25_4v
     .              +msq25_4(j,k)*sub45_2(qq))
 
      elseif ((j .eq. 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2d0*tr*msq14_2(-k,k)*sub14_2(qg)
      msq(3,j,k)=xn*(msq15_2(j,k)*sub15_2(gg)+msq15_2v(j,k)*sub15_2v)
      msq(4,j,k)=xn*msq25_1(j,k)*sub25_1(qq)
      msq(5,j,k)=xn*(msq14_5(j,k)*sub15_4(gg)+msq14_5v(j,k)*sub15_4v
     .              +msq14_5(j,k)*sub45_1(qq))
      msq(6,j,k)=-(msq25_4(j,k)*sub25_4(qq)+msq25_4(j,k)*sub45_2(qq))/xn

      elseif ((j .eq. 0).and.(k .eq. 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> gamma(p3)+q(p4)+qb(p5)
c---Hence 14 split multiplies qb(14)+g(p2) --> gamma(p3)+qb(p5)
c---Hence 24 split multiplies g(p1)+qb(p24) --> gamma(p3)+qb(p5)
      msq(1,j,k)=(msq14_2(-1,k)+msq14_2(-2,k)+msq14_2(-3,k)
     .           +msq14_2(-4,k)+msq14_2(-5,k))*sub14_2(qg)*2d0*tr
      msq(2,j,k)=(msq24_1(k,-1)+msq24_1(k,-2)+msq24_1(k,-3)
     .           +msq24_1(k,-4)+msq24_1(k,-5))*sub24_1(qg)*2d0*tr
      msq(3,j,k)=(msq15_2(+5,k)+msq15_2(+4,k)+msq15_2(+3,k)
     .           +msq15_2(+2,k)+msq15_2(+1,k))*sub15_2(qg)*2d0*tr
      msq(4,j,k)=(msq25_1(k,+5)+msq25_1(k,+4)+msq25_1(k,+3)
     .           +msq25_1(k,+2)+msq25_1(k,+1))*sub25_1(qg)*2d0*tr

      endif

 19   continue
      enddo
      enddo


      do j=-nf,nf
      do k=-nf,nf      

      if (((j .gt. 0).and.(k .gt. 0)) .or. 
     .    ((j .lt. 0).and.(k .lt. 0))) then
c--q-q or qb-qb
      if (j.eq.k) then
      msq(1,j,k)=msq(1,j,k)+0.5d0*(xn-1d0/xn)
     .  *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
      msq(2,j,k)=msq(2,j,k)+0.5d0*(xn-1d0/xn)
     .  *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
      msq(3,j,k)=msq(3,j,k)+0.5d0*(xn-1d0/xn)
     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+0.5d0*(xn-1d0/xn)
     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      else
      msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     .  *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      endif
      elseif ((j .gt. 0).and.(k .lt. 0)) then
c q-qbar
      if (j.eq.-k) then
      msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     .  *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(6,j,k)=msq(6,j,k)+2d0*tr*dfloat(nf)
     . *(msq25_4(j,k)*sub45_2(gq)-msq45_2v(j,k)*sub45_2v)
      else 
      msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     .  *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      endif
c--qbar-q
      elseif ((j .lt. 0).and.(k .gt. 0)) then
      if (j.eq.-k) then
      msq(2,j,k)=msq(2,j,k)+(xn-1d0/xn)
     .  *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1d0/xn)
     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(6,j,k)=msq(6,j,k)+2d0*tr*dfloat(nf)
     . *(msq25_4(j,k)*sub45_2(gq)-msq45_2v(j,k)*sub45_2v)
      else 
      msq(2,j,k)=msq(2,j,k)+(xn-1d0/xn)
     .  *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1d0/xn)
     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)

      endif
      endif


      enddo
      enddo

      return
      end
      
