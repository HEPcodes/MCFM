************************************************************************
*   Subtractions for g --> q qbar (final) or q --> g + qbar (initial)  *
************************************************************************
       subroutine coll_gq(i1,i3,i2,xmsq,jx,kx,lx,mx,fac)
       implicit none
       include 'constants.f'
       include 'qqgg.f'
      
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36),coll_born,coll_corr
      integer i,pntr(5,5,5),xntr(5,5,5)
      integer i1,i2,i3,jx,kx,lx,mx
      double precision xmsq(36),fac   
    
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
    
c      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
c     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
c      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
c     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
c      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
c     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/
      
c--- pntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k
      data pntr/0,0,0,0,0,0,0,0,0,0,0,19,0,20,21,0,25,26,0,27,0,31,32,
     . 33,0,0,0,0,0,0,0,0,0,0,0,1,0,0,22,23,7,0,28,0,29,13,0,34,35,0,0,
     . 0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,
     . 0,0,0,0,0,2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,
     . 3,5,0,6,0,9,11,12,0,0,0,0,0,0,0/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

      coll_born=
     . +subij_k(pntr(i1,i3,i2),gq)*msqij_k(xntr(i1,i3,i2),0,jx,kx,lx,mx)
     . +subij_k(pntr(i1,i3,i2),gq)*msqij_k(xntr(i1,i3,i2),1,jx,kx,lx,mx)
     . +subij_k(pntr(i1,i3,i2),gq)*msqij_k(xntr(i1,i3,i2),2,jx,kx,lx,mx)

      coll_corr=
     . +subij_kv(pntr(i1,i3,i2))*msqij_kv(pntr(i1,i3,i2),0,jx,kx,lx,mx)
     . +subij_kv(pntr(i1,i3,i2))*msqij_kv(pntr(i1,i3,i2),1,jx,kx,lx,mx)
     . +subij_kv(pntr(i1,i3,i2))*msqij_kv(pntr(i1,i3,i2),2,jx,kx,lx,mx)

      xmsq(xntr(i1,i3,i2))=xmsq(xntr(i1,i3,i2))+fac*coll_born
 
c--- both the final-initial and the final-final (gq) have a minus sign      
c--- in the correlation piece, relative to (gg)   
      if (i1.gt.2) then
      xmsq(xntr(i1,i3,i2))=xmsq(xntr(i1,i3,i2))-fac*coll_corr
      else
      xmsq(xntr(i1,i3,i2))=xmsq(xntr(i1,i3,i2))+fac*coll_corr      
      endif
      
      return
      end
      
      
************************************************************************
*   Subtractions for g --> q + qbar (initial)                          *
************************************************************************
       subroutine coll_qg(i1,i3,i2,xmsq,jx,kx,lx,mx,fac)
       implicit none
       include 'constants.f'
       include 'qqgg.f'
      
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36)
      integer i,pntr(5,5,5),xntr(5,5,5)
      integer i1,i2,i3,jx,kx,lx,mx
      double precision xmsq(36),fac   
    
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
    
c      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
c     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
c      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
c     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
c      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
c     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/
      
c--- pntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k
      data pntr/0,0,0,0,0,0,0,0,0,0,0,19,0,20,21,0,25,26,0,27,0,31,32,
     . 33,0,0,0,0,0,0,0,0,0,0,0,1,0,0,22,23,7,0,28,0,29,13,0,34,35,0,0,
     . 0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,
     . 0,0,0,0,0,2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,
     . 3,5,0,6,0,9,11,12,0,0,0,0,0,0,0/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

      xmsq(xntr(i1,i3,i2))=xmsq(xntr(i1,i3,i2))+fac*(
     . +subij_k(pntr(i1,i3,i2),qg)*msqij_k(xntr(i1,i3,i2),0,jx,kx,lx,mx)
     . +subij_k(pntr(i1,i3,i2),qg)*msqij_k(xntr(i1,i3,i2),1,jx,kx,lx,mx)
     . +subij_k(pntr(i1,i3,i2),qg)*msqij_k(xntr(i1,i3,i2),2,jx,kx,lx,mx)
     . )
      
      return
      end
      
      
      subroutine dip2jet_qrqr(i1,i2,i3,i4,i5,xmsq,jx,kx,lx,mx)
************************************************************************
*   General subtraction functions                                      *
*   q(i1) + r(i2) --> q(i3) + r(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
       include 'qqgg.f'
      
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36)
      integer i,pntr(5,5,5),xntr(5,5,5)
      integer i1,i2,i3,i4,i5,jx,kx,lx,mx,jg,kg,lg,mg
      double precision xmsq(36),fac
      
      double precision 
     . sub15_4(4),sub25_3(4),sub15_3(4),sub25_4(4),
     . sub45_1(4),sub35_2(4),sub35_1(4),sub45_2(4),
     . sub15_2(4),sub25_1(4),sub35_4(4),sub45_3(4)
      double precision 
     . msq15_4(0:2),msq25_3(0:2),msq15_3(0:2),msq25_4(0:2),
     . msq45_1(0:2),msq35_2(0:2),msq35_1(0:2),msq45_2(0:2),
     . msq15_2(0:2),msq25_1(0:2),msq35_4(0:2),msq45_3(0:2)
      
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
    
c      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
c     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
c      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
c     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
c      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
c     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/
      
c--- pntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k
      data pntr/0,0,0,0,0,0,0,0,0,0,0,19,0,20,21,0,25,26,0,27,0,31,32,
     . 33,0,0,0,0,0,0,0,0,0,0,0,1,0,0,22,23,7,0,28,0,29,13,0,34,35,0,0,
     . 0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,
     . 0,0,0,0,0,2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,
     . 3,5,0,6,0,9,11,12,0,0,0,0,0,0,0/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

c--- initialize the dipole contributions to zero
      do i=1,36
        xmsq(i)=0d0
      enddo

c--- fill up the subtraction terms and reduced matrix elements
      do i=qq,qq
      sub15_4(i)=subij_k(pntr(i1,i5,i4),i)
      sub45_1(i)=subij_k(pntr(i4,i5,i1),i)
      sub25_3(i)=subij_k(pntr(i2,i5,i3),i)
      sub35_2(i)=subij_k(pntr(i3,i5,i2),i)
      sub15_3(i)=subij_k(pntr(i1,i5,i3),i)
      sub35_1(i)=subij_k(pntr(i3,i5,i1),i)
      sub25_4(i)=subij_k(pntr(i2,i5,i4),i)
      sub45_2(i)=subij_k(pntr(i4,i5,i2),i)
      sub15_2(i)=subij_k(pntr(i1,i5,i2),i)
      sub25_1(i)=subij_k(pntr(i2,i5,i1),i)
      sub35_4(i)=subij_k(pntr(i3,i5,i4),i)
      sub45_3(i)=subij_k(pntr(i4,i5,i3),i)
      enddo
      
      do i=0,2
      msq15_4(i)=msqij_k(xntr(i1,i5,i4),i,jx,kx,lx,mx)
      msq45_1(i)=msqij_k(xntr(i4,i5,i1),i,jx,kx,lx,mx)
      msq25_3(i)=msqij_k(xntr(i2,i5,i3),i,jx,kx,lx,mx)
      msq35_2(i)=msqij_k(xntr(i3,i5,i2),i,jx,kx,lx,mx)
      msq15_3(i)=msqij_k(xntr(i1,i5,i3),i,jx,kx,lx,mx)
      msq35_1(i)=msqij_k(xntr(i3,i5,i1),i,jx,kx,lx,mx)
      msq25_4(i)=msqij_k(xntr(i2,i5,i4),i,jx,kx,lx,mx)
      msq45_2(i)=msqij_k(xntr(i4,i5,i2),i,jx,kx,lx,mx)
      msq15_2(i)=msqij_k(xntr(i1,i5,i2),i,jx,kx,lx,mx)
      msq25_1(i)=msqij_k(xntr(i2,i5,i1),i,jx,kx,lx,mx)
      msq35_4(i)=msqij_k(xntr(i3,i5,i4),i,jx,kx,lx,mx)
      msq45_3(i)=msqij_k(xntr(i4,i5,i3),i,jx,kx,lx,mx)
      enddo
      
c--- fill up the dipole contributions - soft terms
      xmsq(xntr(i1,i5,i4))=xmsq(xntr(i1,i5,i4))+
     . sub15_4(qq)*(2d0*cf-1d0/xn)*msq15_4(0)
      xmsq(xntr(i4,i5,i1))=xmsq(xntr(i4,i5,i1))+
     . sub45_1(qq)*(2d0*cf-1d0/xn)*msq45_1(0)
      xmsq(xntr(i2,i5,i3))=xmsq(xntr(i2,i5,i3))+
     . sub25_3(qq)*(2d0*cf-1d0/xn)*msq25_3(0)
      xmsq(xntr(i3,i5,i2))=xmsq(xntr(i3,i5,i2))+
     . sub35_2(qq)*(2d0*cf-1d0/xn)*msq35_2(0)
      xmsq(xntr(i1,i5,i3))=xmsq(xntr(i1,i5,i3))+
     . sub15_3(qq)*(      -1d0/xn)*msq15_3(0)
      xmsq(xntr(i3,i5,i1))=xmsq(xntr(i3,i5,i1))+
     . sub35_1(qq)*(      -1d0/xn)*msq35_1(0)
      xmsq(xntr(i2,i5,i4))=xmsq(xntr(i2,i5,i4))+
     . sub25_4(qq)*(      -1d0/xn)*msq25_4(0)
      xmsq(xntr(i4,i5,i2))=xmsq(xntr(i4,i5,i2))+
     . sub45_2(qq)*(      -1d0/xn)*msq45_2(0)
      xmsq(xntr(i1,i5,i2))=xmsq(xntr(i1,i5,i2))+
     . sub15_2(qq)*(      +2d0/xn)*msq15_2(0)
      xmsq(xntr(i2,i5,i1))=xmsq(xntr(i2,i5,i1))+
     . sub25_1(qq)*(      +2d0/xn)*msq25_1(0)
      xmsq(xntr(i3,i5,i4))=xmsq(xntr(i3,i5,i4))+
     . sub35_4(qq)*(      +2d0/xn)*msq35_4(0)
      xmsq(xntr(i4,i5,i3))=xmsq(xntr(i4,i5,i3))+
     . sub45_3(qq)*(      +2d0/xn)*msq45_3(0)
      
      return
      end



      subroutine dip2jet_qqqq(i1,i2,i3,i4,i5,xmsq,jx,kx,lx,mx)
************************************************************************
*   General subtraction functions                                      *
*   q(i1) + q(i2) --> q(i3) + q(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
       include 'qqgg.f'
      
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36)
      integer i,pntr(5,5,5),xntr(5,5,5)
      integer i1,i2,i3,i4,i5,jx,kx,lx,mx
      double precision xmsq(36),tmp
      
      double precision 
     . sub15_4(4),sub25_3(4),sub15_3(4),sub25_4(4),
     . sub45_1(4),sub35_2(4),sub35_1(4),sub45_2(4),
     . sub15_2(4),sub25_1(4),sub35_4(4),sub45_3(4)
      double precision 
     . msq15_4(0:2),msq25_3(0:2),msq15_3(0:2),msq25_4(0:2),
     . msq45_1(0:2),msq35_2(0:2),msq35_1(0:2),msq45_2(0:2),
     . msq15_2(0:2),msq25_1(0:2),msq35_4(0:2),msq45_3(0:2)
      
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
    
c      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
c     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
c      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
c     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
c      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
c     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/
      
c--- pntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k
      data pntr/0,0,0,0,0,0,0,0,0,0,0,19,0,20,21,0,25,26,0,27,0,31,32,
     . 33,0,0,0,0,0,0,0,0,0,0,0,1,0,0,22,23,7,0,28,0,29,13,0,34,35,0,0,
     . 0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,
     . 0,0,0,0,0,2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,
     . 3,5,0,6,0,9,11,12,0,0,0,0,0,0,0/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

c--- initialize the dipole contributions to zero
      do i=1,36
        xmsq(i)=0d0
      enddo

c--- fill up the subtraction terms and reduced matrix elements
      do i=qq,qq
      sub15_4(i)=subij_k(pntr(i1,i5,i4),i)
      sub45_1(i)=subij_k(pntr(i4,i5,i1),i)
      sub25_3(i)=subij_k(pntr(i2,i5,i3),i)
      sub35_2(i)=subij_k(pntr(i3,i5,i2),i)
      sub15_3(i)=subij_k(pntr(i1,i5,i3),i)
      sub35_1(i)=subij_k(pntr(i3,i5,i1),i)
      sub25_4(i)=subij_k(pntr(i2,i5,i4),i)
      sub45_2(i)=subij_k(pntr(i4,i5,i2),i)
      sub15_2(i)=subij_k(pntr(i1,i5,i2),i)
      sub25_1(i)=subij_k(pntr(i2,i5,i1),i)
      sub35_4(i)=subij_k(pntr(i3,i5,i4),i)
      sub45_3(i)=subij_k(pntr(i4,i5,i3),i)
      enddo
      
      do i=0,2
      msq15_4(i)=msqij_k(xntr(i1,i5,i4),i,jx,kx,lx,mx)
      msq45_1(i)=msqij_k(xntr(i4,i5,i1),i,jx,kx,lx,mx)
      msq25_3(i)=msqij_k(xntr(i2,i5,i3),i,jx,kx,lx,mx)
      msq35_2(i)=msqij_k(xntr(i3,i5,i2),i,jx,kx,lx,mx)
      msq15_3(i)=msqij_k(xntr(i1,i5,i3),i,jx,kx,lx,mx)
      msq35_1(i)=msqij_k(xntr(i3,i5,i1),i,jx,kx,lx,mx)
      msq25_4(i)=msqij_k(xntr(i2,i5,i4),i,jx,kx,lx,mx)
      msq45_2(i)=msqij_k(xntr(i4,i5,i2),i,jx,kx,lx,mx)
      msq15_2(i)=msqij_k(xntr(i1,i5,i2),i,jx,kx,lx,mx)
      msq25_1(i)=msqij_k(xntr(i2,i5,i1),i,jx,kx,lx,mx)
      msq35_4(i)=msqij_k(xntr(i3,i5,i4),i,jx,kx,lx,mx)
      msq45_3(i)=msqij_k(xntr(i4,i5,i3),i,jx,kx,lx,mx)
      enddo      
      
c--- fill up the dipole contributions
      xmsq(xntr(i1,i5,i4))=xmsq(xntr(i1,i5,i4))+
     . sub15_4(qq)*2d0*cf*msq15_4(2)
      xmsq(xntr(i4,i5,i1))=xmsq(xntr(i4,i5,i1))+
     . sub45_1(qq)*2d0*cf*msq45_1(2)
      xmsq(xntr(i2,i5,i3))=xmsq(xntr(i2,i5,i3))+
     . sub25_3(qq)*2d0*cf*msq25_3(2)
      xmsq(xntr(i3,i5,i2))=xmsq(xntr(i3,i5,i2))+
     . sub35_2(qq)*2d0*cf*msq35_2(2)
      xmsq(xntr(i1,i5,i3))=xmsq(xntr(i1,i5,i3))+
     . sub15_3(qq)*2d0*cf*msq15_3(1)
      xmsq(xntr(i3,i5,i1))=xmsq(xntr(i3,i5,i1))+
     . sub35_1(qq)*2d0*cf*msq35_1(1)
      xmsq(xntr(i2,i5,i4))=xmsq(xntr(i2,i5,i4))+
     . sub25_4(qq)*2d0*cf*msq25_4(1)
      xmsq(xntr(i4,i5,i2))=xmsq(xntr(i4,i5,i2))+
     . sub45_2(qq)*2d0*cf*msq45_2(1)
      xmsq(xntr(i1,i5,i2))=xmsq(xntr(i1,i5,i2))+
     . sub15_2(qq)*2d0*cf*msq15_2(0)
      xmsq(xntr(i2,i5,i1))=xmsq(xntr(i2,i5,i1))+
     . sub25_1(qq)*2d0*cf*msq25_1(0)
      xmsq(xntr(i3,i5,i4))=xmsq(xntr(i3,i5,i4))+
     . sub35_4(qq)*2d0*cf*msq35_4(0)
      xmsq(xntr(i4,i5,i3))=xmsq(xntr(i4,i5,i3))+
     . sub45_3(qq)*2d0*cf*msq45_3(0)
      
      return
      end



      subroutine dip2jet_qqgg(i1,i2,i3,i4,i5,xmsq,jx,kx,lx,mx,permfac)
************************************************************************
*   General subtraction functions                                      *
*   q(i1) + q(i2) --> g(i3) + g(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
       include 'qqgg.f'
      
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36)
      integer i,j,icrit,pntr(5,5,5),xntr(5,5,5)
      integer i1,i2,i3,i4,i5,jx,kx,lx,mx
      double precision xmsq(36),tmp,permfac,xfac(5,5)
      
      double precision 
     . sub15_4(4),sub25_3(4),sub15_3(4),sub25_4(4),
     . sub45_1(4),sub35_2(4),sub35_1(4),sub45_2(4),
     . sub15_2(4),sub25_1(4),sub35_4(4),sub45_3(4),
     . sub45_1v,sub35_2v,sub35_1v,sub45_2v,
     . sub35_4v,sub45_3v
      double precision 
     . msq15_4(0:2),msq25_3(0:2),msq15_3(0:2),msq25_4(0:2),
     . msq45_1(0:2),msq35_2(0:2),msq35_1(0:2),msq45_2(0:2),
     . msq15_2(0:2),msq25_1(0:2),msq35_4(0:2),msq45_3(0:2),
     . msq45_1v(0:2),msq35_2v(0:2),msq35_1v(0:2),msq45_2v(0:2),
     . msq35_4v(0:2),msq45_3v(0:2)
      
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
    
c      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
c     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
c      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
c     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
c      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
c     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/
      
c--- pntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k
      data pntr/0,0,0,0,0,0,0,0,0,0,0,19,0,20,21,0,25,26,0,27,0,31,32,
     . 33,0,0,0,0,0,0,0,0,0,0,0,1,0,0,22,23,7,0,28,0,29,13,0,34,35,0,0,
     . 0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,
     . 0,0,0,0,0,2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,
     . 3,5,0,6,0,9,11,12,0,0,0,0,0,0,0/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

c--- initialize the dipole contributions to zero
      do i=1,36
        xmsq(i)=0d0
      enddo

      if     (i3+i4+i5 .eq. 7) then
c--- gg --> qqg case, so no factors of half
        icrit=6
      elseif (i3+i4+i5 .eq. 12) then
c--- qq --> ggg case, so all of 34,35 and 45 need factor of half
        icrit=3
      else
c--- qg --> qgg case, so only 45 needs factor of half
        icrit=4
      endif

      do i=1,5
      do j=1,5
        if ((i .ge. icrit) .and. (j .ge.icrit)) then
          xfac(i,j)=permfac
        else
          xfac(i,j)=1d0
        endif  
      enddo
      enddo
      
c--- fill up the subtraction terms and reduced matrix elements
c--- NOTE: final-final singularities have an extra factor of "permfac"
c---  since they may be counted twice when summing over permutations
      do i=qq,gg,gg-qq
      sub15_4(i)=xfac(i1,i5)*subij_k(pntr(i1,i5,i4),i)
      sub45_1(i)=xfac(i4,i5)*subij_k(pntr(i4,i5,i1),i)
      sub25_3(i)=xfac(i2,i5)*subij_k(pntr(i2,i5,i3),i)
      sub35_2(i)=xfac(i3,i5)*subij_k(pntr(i3,i5,i2),i)
      sub15_3(i)=xfac(i1,i5)*subij_k(pntr(i1,i5,i3),i)
      sub35_1(i)=xfac(i3,i5)*subij_k(pntr(i3,i5,i1),i)
      sub25_4(i)=xfac(i2,i5)*subij_k(pntr(i2,i5,i4),i)
      sub45_2(i)=xfac(i4,i5)*subij_k(pntr(i4,i5,i2),i)
      sub15_2(i)=xfac(i1,i5)*subij_k(pntr(i1,i5,i2),i)
      sub25_1(i)=xfac(i2,i5)*subij_k(pntr(i2,i5,i1),i)
      sub35_4(i)=xfac(i3,i5)*subij_k(pntr(i3,i5,i4),i)
      sub45_3(i)=xfac(i4,i5)*subij_k(pntr(i4,i5,i3),i)
      enddo
      sub45_1v=xfac(i4,i5)*subij_kv(pntr(i4,i5,i1))
      sub35_2v=xfac(i3,i5)*subij_kv(pntr(i3,i5,i2))
      sub35_1v=xfac(i3,i5)*subij_kv(pntr(i3,i5,i1))
      sub45_2v=xfac(i4,i5)*subij_kv(pntr(i4,i5,i2))
      sub35_4v=xfac(i3,i5)*subij_kv(pntr(i3,i5,i4))
      sub45_3v=xfac(i4,i5)*subij_kv(pntr(i4,i5,i3))
            
      do i=0,2
      msq15_4(i)=msqij_k(xntr(i1,i5,i4),i,jx,kx,lx,mx)
      msq45_1(i)=msqij_k(xntr(i4,i5,i1),i,jx,kx,lx,mx)
      msq25_3(i)=msqij_k(xntr(i2,i5,i3),i,jx,kx,lx,mx)
      msq35_2(i)=msqij_k(xntr(i3,i5,i2),i,jx,kx,lx,mx)
      msq15_3(i)=msqij_k(xntr(i1,i5,i3),i,jx,kx,lx,mx)
      msq35_1(i)=msqij_k(xntr(i3,i5,i1),i,jx,kx,lx,mx)
      msq25_4(i)=msqij_k(xntr(i2,i5,i4),i,jx,kx,lx,mx)
      msq45_2(i)=msqij_k(xntr(i4,i5,i2),i,jx,kx,lx,mx)
      msq15_2(i)=msqij_k(xntr(i1,i5,i2),i,jx,kx,lx,mx)
      msq25_1(i)=msqij_k(xntr(i2,i5,i1),i,jx,kx,lx,mx)
      msq35_4(i)=msqij_k(xntr(i3,i5,i4),i,jx,kx,lx,mx)
      msq45_3(i)=msqij_k(xntr(i4,i5,i3),i,jx,kx,lx,mx)

      msq45_1v(i)=msqij_kv(pntr(i4,i5,i1),i,jx,kx,lx,mx)
      msq35_2v(i)=msqij_kv(pntr(i3,i5,i2),i,jx,kx,lx,mx)
      msq35_1v(i)=msqij_kv(pntr(i3,i5,i1),i,jx,kx,lx,mx)
      msq45_2v(i)=msqij_kv(pntr(i4,i5,i2),i,jx,kx,lx,mx)
      msq35_4v(i)=msqij_kv(pntr(i3,i5,i4),i,jx,kx,lx,mx)
      msq45_3v(i)=msqij_kv(pntr(i4,i5,i3),i,jx,kx,lx,mx)
      enddo      
      
c--- fill up the dipole contributions
      xmsq(xntr(i1,i5,i4))=xmsq(xntr(i1,i5,i4))
     . +xn*sub15_4(qq)*(msq15_4(0)+msq15_4(2))
      xmsq(xntr(i4,i5,i1))=xmsq(xntr(i4,i5,i1))
     . +xn*sub45_1(gg)*(msq45_1(0)+msq45_1(2))
     . +xn*sub45_1v*(msq45_1v(0)+msq45_1v(2))
      xmsq(xntr(i2,i5,i3))=xmsq(xntr(i2,i5,i3))
     . +xn*sub25_3(qq)*(msq25_3(0)+msq25_3(2))
      xmsq(xntr(i3,i5,i2))=xmsq(xntr(i3,i5,i2))
     . +xn*sub35_2(gg)*(msq35_2(0)+msq35_2(2))
     . +xn*sub35_2v*(msq35_2v(0)+msq35_2v(2))
      xmsq(xntr(i1,i5,i3))=xmsq(xntr(i1,i5,i3))
     . +xn*sub15_3(qq)*(msq15_3(0)+msq15_3(1))
      xmsq(xntr(i3,i5,i1))=xmsq(xntr(i3,i5,i1))
     . +xn*sub35_1(gg)*(msq35_1(0)+msq35_1(1))
     . +xn*sub35_1v*(msq35_1v(0)+msq35_1v(1))
      xmsq(xntr(i2,i5,i4))=xmsq(xntr(i2,i5,i4))
     . +xn*sub25_4(qq)*(msq25_4(0)+msq25_4(1))
      xmsq(xntr(i4,i5,i2))=xmsq(xntr(i4,i5,i2))
     . +xn*sub45_2(gg)*(msq45_2(0)+msq45_2(1))
     . +xn*sub45_2v*(msq45_2v(0)+msq45_2v(1))
      xmsq(xntr(i1,i5,i2))=xmsq(xntr(i1,i5,i2))
     . +sub15_2(qq)*(-(xn+1d0/xn)*msq15_2(0)
     .              -1d0/xn*(msq15_2(1)+msq15_2(2)))
      xmsq(xntr(i2,i5,i1))=xmsq(xntr(i2,i5,i1))
     . +sub25_1(qq)*(-(xn+1d0/xn)*msq25_1(0)
     .              -1d0/xn*(msq25_1(1)+msq25_1(2)))
      xmsq(xntr(i3,i5,i4))=xmsq(xntr(i3,i5,i4))
     . +xn*sub35_4(gg)*(msq35_4(1)+msq35_4(2))
     . +xn*sub35_4v*(msq35_4v(1)+msq35_4v(2))
      xmsq(xntr(i4,i5,i3))=xmsq(xntr(i4,i5,i3))
     . +xn*sub45_3(gg)*(msq45_3(1)+msq45_3(2))
     . +xn*sub45_3v*(msq45_3v(1)+msq45_3v(2))
      
      return
      end



      subroutine dip2jet_gggg(i1,i2,i3,i4,i5,xmsq,jx,kx,lx,mx)
************************************************************************
*   General subtraction functions                                      *
*   g(i1) + g(i2) --> g(i3) + g(i4) + g(i5)                            *
************************************************************************
       implicit none
       include 'constants.f'
       include 'qqgg.f'
      
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36)
      integer i,pntr(5,5,5),xntr(5,5,5)
      integer i1,i2,i3,i4,i5,jx,kx,lx,mx
      double precision xmsq(36),tmp
      
      double precision 
     . sub15_4(4),sub25_3(4),sub15_3(4),sub25_4(4),
     . sub45_1(4),sub35_2(4),sub35_1(4),sub45_2(4),
     . sub15_2(4),sub25_1(4),sub35_4(4),sub45_3(4),
     . sub15_4v,sub25_3v,sub15_3v,sub25_4v,
     . sub45_1v,sub35_2v,sub35_1v,sub45_2v,
     . sub15_2v,sub25_1v,sub35_4v,sub45_3v
      double precision 
     . msq15_4(0:2),msq25_3(0:2),msq15_3(0:2),msq25_4(0:2),
     . msq45_1(0:2),msq35_2(0:2),msq35_1(0:2),msq45_2(0:2),
     . msq15_2(0:2),msq25_1(0:2),msq35_4(0:2),msq45_3(0:2),
     . msq15_4v(0:2),msq25_3v(0:2),msq15_3v(0:2),msq25_4v(0:2),
     . msq45_1v(0:2),msq35_2v(0:2),msq35_1v(0:2),msq45_2v(0:2),
     . msq15_2v(0:2),msq25_1v(0:2),msq35_4v(0:2),msq45_3v(0:2)
      
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
    
c      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
c     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
c      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
c     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
c      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
c     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/
      
c--- pntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k
      data pntr/0,0,0,0,0,0,0,0,0,0,0,19,0,20,21,0,25,26,0,27,0,31,32,
     . 33,0,0,0,0,0,0,0,0,0,0,0,1,0,0,22,23,7,0,28,0,29,13,0,34,35,0,0,
     . 0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,
     . 0,0,0,0,0,2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,
     . 3,5,0,6,0,9,11,12,0,0,0,0,0,0,0/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

c--- initialize the dipole contributions to zero
      do i=1,36
        xmsq(i)=0d0
      enddo

c--- fill up the subtraction terms and reduced matrix elements
c--- NOTE: final-final singularities have an extra factor of 1/2
c---  since they are counted twice when summing over permutations
      do i=gg,gg
      sub15_4(i)=subij_k(pntr(i1,i5,i4),i)
      sub45_1(i)=half*subij_k(pntr(i4,i5,i1),i)
      sub25_3(i)=subij_k(pntr(i2,i5,i3),i)
      sub35_2(i)=half*subij_k(pntr(i3,i5,i2),i)
      sub15_3(i)=subij_k(pntr(i1,i5,i3),i)
      sub35_1(i)=half*subij_k(pntr(i3,i5,i1),i)
      sub25_4(i)=subij_k(pntr(i2,i5,i4),i)
      sub45_2(i)=half*subij_k(pntr(i4,i5,i2),i)
      sub15_2(i)=subij_k(pntr(i1,i5,i2),i)
      sub25_1(i)=subij_k(pntr(i2,i5,i1),i)
      sub35_4(i)=half*subij_k(pntr(i3,i5,i4),i)
      sub45_3(i)=half*subij_k(pntr(i4,i5,i3),i)
      enddo
      sub15_4v=subij_kv(pntr(i1,i5,i4))
      sub45_1v=half*subij_kv(pntr(i4,i5,i1))
      sub25_3v=subij_kv(pntr(i2,i5,i3))
      sub35_2v=half*subij_kv(pntr(i3,i5,i2))
      sub15_3v=subij_kv(pntr(i1,i5,i3))
      sub35_1v=half*subij_kv(pntr(i3,i5,i1))
      sub25_4v=subij_kv(pntr(i2,i5,i4))
      sub45_2v=half*subij_kv(pntr(i4,i5,i2))
      sub15_2v=subij_kv(pntr(i1,i5,i2))
      sub25_1v=subij_kv(pntr(i2,i5,i1))
      sub35_4v=half*subij_kv(pntr(i3,i5,i4))
      sub45_3v=half*subij_kv(pntr(i4,i5,i3))
      
      do i=0,2
      msq15_4(i)=msqij_k(xntr(i1,i5,i4),i,jx,kx,lx,mx)
      msq45_1(i)=msqij_k(xntr(i4,i5,i1),i,jx,kx,lx,mx)
      msq25_3(i)=msqij_k(xntr(i2,i5,i3),i,jx,kx,lx,mx)
      msq35_2(i)=msqij_k(xntr(i3,i5,i2),i,jx,kx,lx,mx)
      msq15_3(i)=msqij_k(xntr(i1,i5,i3),i,jx,kx,lx,mx)
      msq35_1(i)=msqij_k(xntr(i3,i5,i1),i,jx,kx,lx,mx)
      msq25_4(i)=msqij_k(xntr(i2,i5,i4),i,jx,kx,lx,mx)
      msq45_2(i)=msqij_k(xntr(i4,i5,i2),i,jx,kx,lx,mx)
      msq15_2(i)=msqij_k(xntr(i1,i5,i2),i,jx,kx,lx,mx)
      msq25_1(i)=msqij_k(xntr(i2,i5,i1),i,jx,kx,lx,mx)
      msq35_4(i)=msqij_k(xntr(i3,i5,i4),i,jx,kx,lx,mx)
      msq45_3(i)=msqij_k(xntr(i4,i5,i3),i,jx,kx,lx,mx)

      msq15_4v(i)=msqij_kv(pntr(i1,i5,i4),i,jx,kx,lx,mx)
      msq45_1v(i)=msqij_kv(pntr(i4,i5,i1),i,jx,kx,lx,mx)
      msq25_3v(i)=msqij_kv(pntr(i2,i5,i3),i,jx,kx,lx,mx)
      msq35_2v(i)=msqij_kv(pntr(i3,i5,i2),i,jx,kx,lx,mx)
      msq15_3v(i)=msqij_kv(pntr(i1,i5,i3),i,jx,kx,lx,mx)
      msq35_1v(i)=msqij_kv(pntr(i3,i5,i1),i,jx,kx,lx,mx)
      msq25_4v(i)=msqij_kv(pntr(i2,i5,i4),i,jx,kx,lx,mx)
      msq45_2v(i)=msqij_kv(pntr(i4,i5,i2),i,jx,kx,lx,mx)
      msq15_2v(i)=msqij_kv(pntr(i1,i5,i2),i,jx,kx,lx,mx)
      msq25_1v(i)=msqij_kv(pntr(i2,i5,i1),i,jx,kx,lx,mx)
      msq35_4v(i)=msqij_kv(pntr(i3,i5,i4),i,jx,kx,lx,mx)
      msq45_3v(i)=msqij_kv(pntr(i4,i5,i3),i,jx,kx,lx,mx)
      enddo      
      
c--- fill up the dipole contributions
      xmsq(xntr(i1,i5,i4))=xmsq(xntr(i1,i5,i4))
     . +sub15_4(gg)*(msq15_4(2)+half*msq15_4(0)+half*msq15_4(1))
     . +sub15_4v*(msq15_4v(2)+half*msq15_4v(0)+half*msq15_4v(1))
      xmsq(xntr(i4,i5,i1))=xmsq(xntr(i4,i5,i1))
     . +sub45_1(gg)*(msq45_1(2)+half*msq45_1(0)+half*msq45_1(1))
     . +sub45_1v*(msq45_1v(2)+half*msq45_1v(0)+half*msq45_1v(1))
      xmsq(xntr(i2,i5,i3))=xmsq(xntr(i2,i5,i3))
     . +sub25_3(gg)*(msq25_3(2)+half*msq25_3(0)+half*msq25_3(1))
     . +sub25_3v*(msq25_3v(2)+half*msq25_3v(0)+half*msq25_3v(1))
      xmsq(xntr(i3,i5,i2))=xmsq(xntr(i3,i5,i2))
     . +sub35_2(gg)*(msq35_2(2)+half*msq35_2(0)+half*msq35_2(1))
     . +sub35_2v*(msq35_2v(2)+half*msq35_2v(0)+half*msq35_2v(1))
      xmsq(xntr(i1,i5,i3))=xmsq(xntr(i1,i5,i3))
     . +sub15_3(gg)*(msq15_3(1)+half*msq15_3(0)+half*msq15_3(2))
     . +sub15_3v*(msq15_3v(1)+half*msq15_3v(0)+half*msq15_3v(2))
      xmsq(xntr(i3,i5,i1))=xmsq(xntr(i3,i5,i1))
     . +sub35_1(gg)*(msq35_1(1)+half*msq35_1(0)+half*msq35_1(2))
     . +sub35_1v*(msq35_1v(1)+half*msq35_1v(0)+half*msq35_1v(2))
      xmsq(xntr(i2,i5,i4))=xmsq(xntr(i2,i5,i4))
     . +sub25_4(gg)*(msq25_4(1)+half*msq25_4(0)+half*msq25_4(2))
     . +sub25_4v*(msq25_4v(1)+half*msq25_4v(0)+half*msq25_4v(2))
      xmsq(xntr(i4,i5,i2))=xmsq(xntr(i4,i5,i2))
     . +sub45_2(gg)*(msq45_2(1)+half*msq45_2(0)+half*msq45_2(2))
     . +sub45_2v*(msq45_2v(1)+half*msq45_2v(0)+half*msq45_2v(2))
      xmsq(xntr(i1,i5,i2))=xmsq(xntr(i1,i5,i2))
     . +sub15_2(gg)*(msq15_2(0)+half*msq15_2(1)+half*msq15_2(2))
     . +sub15_2v*(msq15_2v(0)+half*msq15_2v(1)+half*msq15_2v(2))
      xmsq(xntr(i2,i5,i1))=xmsq(xntr(i2,i5,i1))
     . +sub25_1(gg)*(msq25_1(0)+half*msq25_1(1)+half*msq25_1(2))
     . +sub25_1v*(msq25_1v(0)+half*msq25_1v(1)+half*msq25_1v(2))
      xmsq(xntr(i3,i5,i4))=xmsq(xntr(i3,i5,i4))
     . +sub35_4(gg)*(msq35_4(0)+half*msq35_4(1)+half*msq35_4(2))
     . +sub35_4v*(msq35_4v(0)+half*msq35_4v(1)+half*msq35_4v(2))
      xmsq(xntr(i4,i5,i3))=xmsq(xntr(i4,i5,i3))
     . +sub45_3(gg)*(msq45_3(0)+half*msq45_3(1)+half*msq45_3(2))
     . +sub45_3v*(msq45_3v(0)+half*msq45_3v(1)+half*msq45_3v(2))
            
      return
      end






c      do i=1,5
c      do j=1,5
c      do k=1,5
c      pntr(i,j,k)=0
c      enddo
c      enddo
c      enddo
c      
c      do i=1,36
c        pntr(ip(i),jp(i),kp(i))=i
c      enddo
c      
c      write(6,*) pntr
c      stop
      
c      do i=1,5
c      do j=1,5
c      do k=1,5
c        if ((i .ge. 3) .and. (k .le. 2)) then
c          xntr(i,j,k)=pntr(k,j,i)
c        else
c          xntr(i,j,k)=pntr(i,j,k)
c        endif
c      enddo
c      enddo
c      enddo
c      
c      write(6,*) xntr
c      stop

