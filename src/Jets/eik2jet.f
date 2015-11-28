      subroutine eik2jet_gggg(i1,i2,i3,i4,i5,xmsq)
************************************************************************
*   General eikole functions                                           *
*   g(i1) + g(i2) --> g(i3) + g(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
      integer i1,i2,i3,i4,i5,di(42),dk(42),dj(42)
      integer i,j,k,pntr(5,5,5),n,perm(3:5,3:5),pab,pba
c--- di(n), dk(n), dj(n) tell the eikole (ik_j) for a pointer n  
      double precision eik1a_b(6),eikba_1(6),eikab_c(6),
     .                 eikcb_a(6),eikbc_2(6),eik2c_b(6),
     .                 eik1b_2(6),eik2b_1(6),eik(42),
     .                 msq_ac(6,-nf:nf,-nf:nf),
     .                 msq_acx(6,0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     .                 msq_ac_1b(6,-nf:nf,-nf:nf),
     .                 msq_ac_1d(6,-nf:nf,-nf:nf),
     .                 msq_ac_2d(6,-nf:nf,-nf:nf),
     .                 msq_ac_3d(6,-nf:nf,-nf:nf)
      double precision eik13,eik14,eik23,eik24,eik12,eik34,xmsq
      common/eikeik/eik1a_b,eikba_1,eik2c_b,eikbc_2,eikab_c,eikcb_a,
     . eik1b_2,eik2b_1,msq_ac,msq_ac_1b,msq_ac_1d,msq_ac_2d,msq_ac_3d,
     . msq_acx
      data di/1,1,1,1,1,1,4,5,3,5,3,4,4,5,3,5,3,4,2,2,2,2,2,2,
     . 3,3,5,4,4,5,5,4,4,3,5,3,1,1,1,2,2,2/
      data dj/3,3,5,4,4,5,3,3,5,4,4,5,5,4,4,3,5,3,5,4,4,3,5,3,
     . 4,5,3,5,3,4,4,5,3,5,3,4,4,5,3,4,5,3/
      data dk/4,5,3,5,3,4,1,1,1,1,1,1,2,2,2,2,2,2,4,5,3,5,3,4,
     . 5,4,4,3,5,3,3,3,5,4,4,5,2,2,2,1,1,1/
      data pntr/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,42, 0, 7, 8, 0,40,11,
     . 0,10, 0,41, 9,12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,39, 0, 0,18,
     . 16,37, 0,15, 0,14,38, 0,17,13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 5,21, 0, 0,31, 3,23, 0,32, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 1,24, 0, 0,35, 0, 0, 0, 0, 0, 6,19,34, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 0, 0, 0, 2,22, 0,33, 0, 4,20,36, 0, 0, 0, 0, 0,
     . 0, 0/
      data perm/0,4,6,2,0,3,1,5,0/
      
c--- fill the "general" eikonal array
      do n=1,6
      eik(   n)=eik1a_b(n)
      eik( 6+n)=eikba_1(n)
      eik(12+n)=eikbc_2(n)
      eik(18+n)=eik2c_b(n)
      eik(24+n)=eikab_c(n)
      eik(30+n)=eikcb_a(n)
      enddo      
      do n=1,3
      eik(36+n)=eik1b_2(n)
      eik(39+n)=eik2b_1(n)
      enddo      

      eik14=eik(pntr(i1,i5,i4))+eik(pntr(i4,i5,i1))
      eik23=eik(pntr(i2,i5,i3))+eik(pntr(i3,i5,i2))
      eik13=eik(pntr(i1,i5,i3))+eik(pntr(i3,i5,i1))
      eik24=eik(pntr(i2,i5,i4))+eik(pntr(i4,i5,i2))
      eik12=eik(pntr(i1,i5,i2))+eik(pntr(i2,i5,i1))
      eik34=eik(pntr(i3,i5,i4))+eik(pntr(i4,i5,i3))
      
      pab=perm(3,4)
      pba=perm(4,3)
      
      xmsq=0d0
      xmsq=xmsq+eik14*(msq_acx(pab,2,0,0,0,0)
     . +half*msq_acx(pab,0,0,0,0,0)+half*msq_acx(pab,1,0,0,0,0))
      xmsq=xmsq+eik23*(msq_acx(pab,2,0,0,0,0)
     . +half*msq_acx(pab,0,0,0,0,0)+half*msq_acx(pab,1,0,0,0,0))
      xmsq=xmsq+eik13*(msq_acx(pab,1,0,0,0,0)
     . +half*msq_acx(pab,0,0,0,0,0)+half*msq_acx(pab,2,0,0,0,0))
      xmsq=xmsq+eik24*(msq_acx(pab,1,0,0,0,0)
     . +half*msq_acx(pab,0,0,0,0,0)+half*msq_acx(pab,2,0,0,0,0))
      xmsq=xmsq+eik12*(msq_acx(pab,0,0,0,0,0)
     . +half*msq_acx(pab,1,0,0,0,0)+half*msq_acx(pab,2,0,0,0,0))
      xmsq=xmsq+eik34*(msq_acx(pab,0,0,0,0,0)
     . +half*msq_acx(pab,1,0,0,0,0)+half*msq_acx(pab,2,0,0,0,0))

      return
      end
      

      subroutine eik2jet_qrqr(i1,i2,i3,i4,i5,xmsq,jp,kp,lp,mp)
************************************************************************
*   General eikole functions                                           *
*   q(i1) + r(i2) --> q(i3) + r(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
      integer i1,i2,i3,i4,i5,di(42),dk(42),dj(42)
      integer jp,kp,lp,mp
      integer i,j,k,pntr(5,5,5),n,perm(3:5,3:5),pab,pba
c--- di(n), dk(n), dj(n) tell the eikole (ik_j) for a pointer n  
      double precision eik1a_b(6),eikba_1(6),eikab_c(6),
     .                 eikcb_a(6),eikbc_2(6),eik2c_b(6),
     .                 eik1b_2(6),eik2b_1(6),eik(42),
     .                 msq_ac(6,-nf:nf,-nf:nf),
     .                 msq_acx(6,0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     .                 msq_ac_1b(6,-nf:nf,-nf:nf),
     .                 msq_ac_1d(6,-nf:nf,-nf:nf),
     .                 msq_ac_2d(6,-nf:nf,-nf:nf),
     .                 msq_ac_3d(6,-nf:nf,-nf:nf)
      double precision eik13,eik14,eik23,eik24,eik12,eik34,xmsq
      common/eikeik/eik1a_b,eikba_1,eik2c_b,eikbc_2,eikab_c,eikcb_a,
     . eik1b_2,eik2b_1,msq_ac,msq_ac_1b,msq_ac_1d,msq_ac_2d,msq_ac_3d,
     . msq_acx
      data di/1,1,1,1,1,1,4,5,3,5,3,4,4,5,3,5,3,4,2,2,2,2,2,2,
     . 3,3,5,4,4,5,5,4,4,3,5,3,1,1,1,2,2,2/
      data dj/3,3,5,4,4,5,3,3,5,4,4,5,5,4,4,3,5,3,5,4,4,3,5,3,
     . 4,5,3,5,3,4,4,5,3,5,3,4,4,5,3,4,5,3/
      data dk/4,5,3,5,3,4,1,1,1,1,1,1,2,2,2,2,2,2,4,5,3,5,3,4,
     . 5,4,4,3,5,3,3,3,5,4,4,5,2,2,2,1,1,1/
      data pntr/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,42, 0, 7, 8, 0,40,11,
     . 0,10, 0,41, 9,12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,39, 0, 0,18,
     . 16,37, 0,15, 0,14,38, 0,17,13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 5,21, 0, 0,31, 3,23, 0,32, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 1,24, 0, 0,35, 0, 0, 0, 0, 0, 6,19,34, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 0, 0, 0, 2,22, 0,33, 0, 4,20,36, 0, 0, 0, 0, 0,
     . 0, 0/
      data perm/0,4,6,2,0,3,1,5,0/
      
c--- fill the "general" eikonal array
      do n=1,6
      eik(   n)=eik1a_b(n)
      eik( 6+n)=eikba_1(n)
      eik(12+n)=eikbc_2(n)
      eik(18+n)=eik2c_b(n)
      eik(24+n)=eikab_c(n)
      eik(30+n)=eikcb_a(n)
      enddo      
      do n=1,3
      eik(36+n)=eik1b_2(n)
      eik(39+n)=eik2b_1(n)
      enddo      

      eik14=eik(pntr(i1,i5,i4))+eik(pntr(i4,i5,i1))
      eik23=eik(pntr(i2,i5,i3))+eik(pntr(i3,i5,i2))
      eik13=eik(pntr(i1,i5,i3))+eik(pntr(i3,i5,i1))
      eik24=eik(pntr(i2,i5,i4))+eik(pntr(i4,i5,i2))
      eik12=eik(pntr(i1,i5,i2))+eik(pntr(i2,i5,i1))
      eik34=eik(pntr(i3,i5,i4))+eik(pntr(i4,i5,i3))
      
c--- this is hacked at the moment
      if (i1 .lt. i4) then      
        pab=perm(3,4)
        pba=perm(4,3)
      else
        pab=perm(4,3)
        pba=perm(3,4)
      endif        
      
      xmsq=0d0
      xmsq=xmsq+eik14*(2d0*cf-1d0/xn)*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+eik23*(2d0*cf-1d0/xn)*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+eik13*(      -1d0/xn)*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+eik24*(      -1d0/xn)*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+eik12*(      +2d0/xn)*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+eik34*(      +2d0/xn)*msq_acx(pab,0,jp,kp,lp,mp)
      
      return
      end

      
      subroutine eik2jet_qqqq(i1,i2,i3,i4,i5,xmsq,jp,kp,lp,mp)
************************************************************************
*   General eikole functions                                           *
*   q(i1) + q(i2) --> q(i3) + q(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
      integer i1,i2,i3,i4,i5,di(42),dk(42),dj(42)
      integer jp,kp,lp,mp
      integer i,j,k,pntr(5,5,5),n,perm(3:5,3:5),pab,pba
c--- di(n), dk(n), dj(n) tell the eikole (ik_j) for a pointer n  
      double precision eik1a_b(6),eikba_1(6),eikab_c(6),
     .                 eikcb_a(6),eikbc_2(6),eik2c_b(6),
     .                 eik1b_2(6),eik2b_1(6),eik(42),
     .                 msq_ac(6,-nf:nf,-nf:nf),
     .                 msq_acx(6,0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     .                 msq_ac_1b(6,-nf:nf,-nf:nf),
     .                 msq_ac_1d(6,-nf:nf,-nf:nf),
     .                 msq_ac_2d(6,-nf:nf,-nf:nf),
     .                 msq_ac_3d(6,-nf:nf,-nf:nf)
      double precision eik13,eik14,eik23,eik24,eik12,eik34,xmsq
      common/eikeik/eik1a_b,eikba_1,eik2c_b,eikbc_2,eikab_c,eikcb_a,
     . eik1b_2,eik2b_1,msq_ac,msq_ac_1b,msq_ac_1d,msq_ac_2d,msq_ac_3d,
     . msq_acx
      data di/1,1,1,1,1,1,4,5,3,5,3,4,4,5,3,5,3,4,2,2,2,2,2,2,
     . 3,3,5,4,4,5,5,4,4,3,5,3,1,1,1,2,2,2/
      data dj/3,3,5,4,4,5,3,3,5,4,4,5,5,4,4,3,5,3,5,4,4,3,5,3,
     . 4,5,3,5,3,4,4,5,3,5,3,4,4,5,3,4,5,3/
      data dk/4,5,3,5,3,4,1,1,1,1,1,1,2,2,2,2,2,2,4,5,3,5,3,4,
     . 5,4,4,3,5,3,3,3,5,4,4,5,2,2,2,1,1,1/
      data pntr/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,42, 0, 7, 8, 0,40,11,
     . 0,10, 0,41, 9,12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,39, 0, 0,18,
     . 16,37, 0,15, 0,14,38, 0,17,13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 5,21, 0, 0,31, 3,23, 0,32, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 1,24, 0, 0,35, 0, 0, 0, 0, 0, 6,19,34, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 0, 0, 0, 2,22, 0,33, 0, 4,20,36, 0, 0, 0, 0, 0,
     . 0, 0/
      data perm/0,4,6,2,0,3,1,5,0/
      
c--- fill the "general" eikonal array
      do n=1,6
      eik(   n)=eik1a_b(n)
      eik( 6+n)=eikba_1(n)
      eik(12+n)=eikbc_2(n)
      eik(18+n)=eik2c_b(n)
      eik(24+n)=eikab_c(n)
      eik(30+n)=eikcb_a(n)
      enddo      
      do n=1,3
      eik(36+n)=eik1b_2(n)
      eik(39+n)=eik2b_1(n)
      enddo      

      eik14=eik(pntr(i1,i5,i4))+eik(pntr(i4,i5,i1))
      eik23=eik(pntr(i2,i5,i3))+eik(pntr(i3,i5,i2))
      eik13=eik(pntr(i1,i5,i3))+eik(pntr(i3,i5,i1))
      eik24=eik(pntr(i2,i5,i4))+eik(pntr(i4,i5,i2))
      eik12=eik(pntr(i1,i5,i2))+eik(pntr(i2,i5,i1))
      eik34=eik(pntr(i3,i5,i4))+eik(pntr(i4,i5,i3))
      
c--- this is hacked at the moment
      if (i1 .lt. i4) then      
        pab=perm(3,4)
        pba=perm(4,3)
      else
        pab=perm(4,3)
        pba=perm(3,4)
      endif        
      
      xmsq=0d0
      xmsq=xmsq+2d0*Cf*eik12*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+2d0*Cf*eik34*msq_acx(pab,0,jp,kp,lp,mp)
      xmsq=xmsq+2d0*Cf*eik13*msq_acx(pab,1,jp,kp,lp,mp)
      xmsq=xmsq+2d0*Cf*eik24*msq_acx(pab,1,jp,kp,lp,mp)
      xmsq=xmsq+2d0*Cf*eik14*msq_acx(pab,2,jp,kp,lp,mp)
      xmsq=xmsq+2d0*Cf*eik23*msq_acx(pab,2,jp,kp,lp,mp)
     
      return
      end
      

      subroutine eik2jet_qqgg(i1,i2,i3,i4,i5,xmsq,jp,kp,lp,mp)
************************************************************************
*   General eikole functions                                           *
*   q(i1) + q(i2) --> g(i3) + g(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
      integer i1,i2,i3,i4,i5,di(42),dk(42),dj(42)
      integer jp,kp,lp,mp
      integer i,j,k,pntr(5,5,5),n,perm(3:5,3:5),pab,pba
c--- di(n), dk(n), dj(n) tell the eikole (ik_j) for a pointer n  
      double precision eik1a_b(6),eikba_1(6),eikab_c(6),
     .                 eikcb_a(6),eikbc_2(6),eik2c_b(6),
     .                 eik1b_2(6),eik2b_1(6),eik(42),
     .                 msq_ac(6,-nf:nf,-nf:nf),
     .                 msq_acx(6,0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     .                 msq_ac_1b(6,-nf:nf,-nf:nf),
     .                 msq_ac_1d(6,-nf:nf,-nf:nf),
     .                 msq_ac_2d(6,-nf:nf,-nf:nf),
     .                 msq_ac_3d(6,-nf:nf,-nf:nf)
      double precision eik13,eik14,eik23,eik24,eik12,eik34,xmsq
      common/eikeik/eik1a_b,eikba_1,eik2c_b,eikbc_2,eikab_c,eikcb_a,
     . eik1b_2,eik2b_1,msq_ac,msq_ac_1b,msq_ac_1d,msq_ac_2d,msq_ac_3d,
     . msq_acx
      data di/1,1,1,1,1,1,4,5,3,5,3,4,4,5,3,5,3,4,2,2,2,2,2,2,
     . 3,3,5,4,4,5,5,4,4,3,5,3,1,1,1,2,2,2/
      data dj/3,3,5,4,4,5,3,3,5,4,4,5,5,4,4,3,5,3,5,4,4,3,5,3,
     . 4,5,3,5,3,4,4,5,3,5,3,4,4,5,3,4,5,3/
      data dk/4,5,3,5,3,4,1,1,1,1,1,1,2,2,2,2,2,2,4,5,3,5,3,4,
     . 5,4,4,3,5,3,3,3,5,4,4,5,2,2,2,1,1,1/
      data pntr/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,42, 0, 7, 8, 0,40,11,
     . 0,10, 0,41, 9,12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,39, 0, 0,18,
     . 16,37, 0,15, 0,14,38, 0,17,13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 5,21, 0, 0,31, 3,23, 0,32, 0, 0, 0, 0, 0, 0, 0,
     . 0, 0, 0, 0, 1,24, 0, 0,35, 0, 0, 0, 0, 0, 6,19,34, 0, 0, 0, 0,
     . 0, 0, 0, 0, 0, 0, 0, 0, 2,22, 0,33, 0, 4,20,36, 0, 0, 0, 0, 0,
     . 0, 0/
      data perm/0,4,6,2,0,3,1,5,0/
      
c--- fill the "general" eikonal array
      do n=1,6
      eik(   n)=eik1a_b(n)
      eik( 6+n)=eikba_1(n)
      eik(12+n)=eikbc_2(n)
      eik(18+n)=eik2c_b(n)
      eik(24+n)=eikab_c(n)
      eik(30+n)=eikcb_a(n)
      enddo      
      do n=1,3
      eik(36+n)=eik1b_2(n)
      eik(39+n)=eik2b_1(n)
      enddo      

      eik14=eik(pntr(i1,i5,i4))+eik(pntr(i4,i5,i1))
      eik23=eik(pntr(i2,i5,i3))+eik(pntr(i3,i5,i2))
      eik13=eik(pntr(i1,i5,i3))+eik(pntr(i3,i5,i1))
      eik24=eik(pntr(i2,i5,i4))+eik(pntr(i4,i5,i2))
      eik12=eik(pntr(i1,i5,i2))+eik(pntr(i2,i5,i1))
      eik34=eik(pntr(i3,i5,i4))+eik(pntr(i4,i5,i3))
     

c--- this is hacked at the moment - checks that 1 comes before 2
      if ( (i1 .eq. 2) .or. (i4 .eq. 1)
     . .or. ((i2 .eq. 2) .and. (i1 .ne. 1)) 
     . .or. ((i3 .eq. 1) .and. (i4 .ne. 2)) ) then      
        pab=perm(4,3)
      else
        pab=perm(3,4)
      endif        

      xmsq=0d0
      xmsq=xmsq+xn*eik14*(
     . msq_acx(pab,0,jp,kp,lp,mp)+msq_acx(pab,2,jp,kp,lp,mp))
      xmsq=xmsq+xn*eik23*(
     . msq_acx(pab,0,jp,kp,lp,mp)+msq_acx(pab,2,jp,kp,lp,mp))
      xmsq=xmsq+xn*eik13*(
     . msq_acx(pab,0,jp,kp,lp,mp)+msq_acx(pab,1,jp,kp,lp,mp))
      xmsq=xmsq+xn*eik24*(
     . msq_acx(pab,0,jp,kp,lp,mp)+msq_acx(pab,1,jp,kp,lp,mp))
      xmsq=xmsq+xn*eik34*(
     . msq_acx(pab,1,jp,kp,lp,mp)+msq_acx(pab,2,jp,kp,lp,mp))
      xmsq=xmsq+eik12*(
     .     -xn*msq_acx(pab,0,jp,kp,lp,mp)
     . -1d0/xn*msq_acx(pab,0,jp,kp,lp,mp)
     . -1d0/xn*msq_acx(pab,1,jp,kp,lp,mp)
     . -1d0/xn*msq_acx(pab,2,jp,kp,lp,mp))
      
      return
      end
      

c--- code used to generate pntr      
c      do i=1,5
c      do j=1,5
c      do k=1,5     
c      pntr(i,j,k)=0
c      do n=1,42
c      if ((di(n) .eq. i).and.(dj(n) .eq. j).and. (dk(n) .eq. k)) then
c        pntr(i,j,k)=n
c      endif
c      enddo
c      
c      enddo
c      enddo
c      enddo
c
c      write(6,99) pntr
c
c      pause
c   
c   99 format(125(i2,',')) 
            
      
