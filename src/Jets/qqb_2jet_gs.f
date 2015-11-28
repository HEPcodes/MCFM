      subroutine qqb_2jet_gs(p,msq)
************************************************************************
*                                                                      *
*  Author: J.M. Campbell, October 2002                                 *
*                                                                      *
*  Subtraction terms for the matrix elements for 3 jet production      *
*     f(-p1) + f(-p2) --> f(p3) + f(p4) + f(p5)                        *
*                                                                      *
*                                                                      *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
     
      integer j,k,n
      integer nd
c--- slightly obtuse notation, to simplify declaration lines      
      double precision p(mxpart,4),msq(maxd,fn:nf,fn:nf),third
      
      double precision sub_st(4),subv_st,
     . msqx_st(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . msqvx_st(0:2,-1:1,-1:1,-1:1,-1:1)
    
      double precision msqij_k(36,0:2,fn:nf,fn:nf,fn:nf,fn:nf),
     .                msqij_kv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     .                 subij_k(36,4),subij_kv(36),xmsq(36)
      integer ip(36),jp(36),kp(36),xntr(5,5,5),dipn,m(-5:5)
      
      common/dipij_k/msqij_k,msqij_kv,subij_k,subij_kv
        
      external qqb_2jetx,qqb_2jet_gvecx

      data m/-1,-2,-1,-2,-1,0,1,2,1,2,1/ 

      data ip/1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3,
     .        2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4/
      data jp/3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,
     .        3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5/
      data kp/2,4,5,4,5,5,2,3,5,3,5,5,2,3,4,3,4,4,
     .        1,1,1,2,2,4,1,1,1,2,2,3,1,1,1,2,2,3/

c--- xntr(i,j,k) -> n such that ip(n)=i,jp(n)=j,kp(n)=k for i<3 or k>2
c--- xntr(i,j,k) -> n such that ip(n)=k,jp(n)=j,kp(n)=i for i>2 and k<3
      data xntr/0,0,0,0,0,0,0,0,0,0,0,19,0,2,3,0,25,8,0,9,0,31,14,15,0,
     . 0,0,0,0,0,0,0,0,0,0,1,0,0,4,5,7,0,10,0,11,13,0,16,17,0,0,0,0,0,0,
     . 0,0,0,0,0,0,0,0,0,0,8,10,0,0,30,14,16,0,36,0,0,0,0,0,0,0,0,0,0,0,
     . 2,4,0,0,24,0,0,0,0,0,15,17,18,0,0,0,0,0,0,0,0,0,0,0,0,3,5,0,6,0,
     . 9,11,12,0,0,0,0,0,0,0/

      ndmax=36

c-- initialize the matrix elements to zero
      do j=-nf,nf
      do k=-nf,nf      
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo
            
c---arguments of dipsx:
c---    1  dipole number
c---    2  momentum
c---    3  emitter
c---    4  emitted
c---    5  spectator
c---    6  interference-free subtraction, equiv. to AP kernel for qq,qg
c---    7  correlation piece of subtraction, relevant only for gg,gq
c---    8  lowest order matrix elements at rescaled momentum, msq(j,k)
c---    9  lowest order matrix elements at rescaled momentum
c---        with emitter contracted with appropriate vector, msqv(j,k)
c---   10  appropriate lowest order calculating routine for msqx
c---   11  appropriate lowest order calculating routine for msqvx
c---   12  dummy variable, dimension (0:2,-nf:nf,-nf:nf)
c---   13  lowest order matrix elements with colour and 4 indices,
c---          msqx(i,j,k,l,m)  -  i runs over 0:2, j,k,l,m are flavours
c---   14  dummy variable, dimension (0:2,-nf:nf,-nf:nf)
c---   15  lowest order matrix elements squared, separated by colours, 
c---        contracted with appropriate vector
c---   16  lowest order matrix elements with 4 indices and
c----        contracted with appropriate vector, msqvx(j,k,l,m)

c--- calculate the generic dipoles
      do n=1,36
c--- dipn=n if i<3 or k>2, but if i>2 and k<3 then points to
c--- n that corresponds to kj_i
        dipn=xntr(ip(n),jp(n),kp(n))
        call dipsxx(dipn,p,ip(n),jp(n),kp(n),sub_st,subv_st,
     .                qqb_2jetx,qqb_2jet_gvecx,msqx_st,msqvx_st)
        call storedipx(msqx_st,msqvx_st,msqij_k,msqij_kv,
     .                sub_st,subv_st,subij_k,subij_kv,n)
      enddo

      third=1d0/3d0

c--- sum up contributions      
      do j=-2,2
      do k=-2,2
      
c--- gluon-gluon
      if     ((j .eq. 0) .and. (k.eq.0)) then
        call dip2jet_gggg(1,2,3,4,5,xmsq,0,0,0,0)
        call msqupdate(j,k,msq,xmsq,1d0)
        call dip2jet_gggg(1,2,5,3,4,xmsq,0,0,0,0)
        call msqupdate(j,k,msq,xmsq,1d0)
        call dip2jet_gggg(1,2,4,5,3,xmsq,0,0,0,0)
        call msqupdate(j,k,msq,xmsq,1d0)
        call dip2jet_qqgg(3,4,1,2,5,xmsq,0,0,1,-1,1d0)
        call coll_gq(3,4,1,xmsq,0,0,0,0,2d0)
        call coll_qg(1,3,2,xmsq,-1,0,-1,0,2d0*cf*avegg/aveqg)
        call coll_qg(1,4,2,xmsq,1,0,1,0,2d0*cf*avegg/aveqg)
        call coll_qg(2,3,1,xmsq,0,-1,-1,0,2d0*cf*avegg/aveqg)
        call coll_qg(2,4,1,xmsq,0,1,1,0,2d0*cf*avegg/aveqg)
        call msqupdate(j,k,msq,xmsq,dfloat(nf))

c--- quark-quark or antiquark-antiquark
      elseif (((j .gt. 0) .and. (k. gt. 0)) .or.
     .        ((j .lt. 0) .and. (k. lt. 0))) then
c---     (identical)
        if (j .eq. k) then
         call dip2jet_qqqq(1,2,3,4,5,xmsq,1,1,1,1)
         call coll_gq(1,3,2,xmsq,0,1,1,0,half*aveqq/aveqg)
         call coll_gq(2,4,1,xmsq,1,0,1,0,half*aveqq/aveqg)
         call coll_gq(1,4,2,xmsq,0,1,1,0,half*aveqq/aveqg)
         call coll_gq(2,3,1,xmsq,1,0,1,0,half*aveqq/aveqg)
         call msqupdate(j,k,msq,xmsq,1d0)
c---     (non-identical)
        else
          call dip2jet_qrqr(1,2,3,4,5,xmsq,1,2,1,2)
          call coll_gq(1,3,2,xmsq,0,1,1,0,aveqq/aveqg)
          call coll_gq(2,4,1,xmsq,1,0,1,0,aveqq/aveqg)
          call msqupdate(j,k,msq,xmsq,1d0)          
        endif

c--- quark-antiquark
      elseif ((j .gt. 0) .and. (k. lt. 0)) then
c---     (identical)
        if (j .eq. -k) then
         call dip2jet_qqqq(1,4,3,2,5,xmsq,-2,2,-2,2)
         call coll_gq(3,4,1,xmsq,1,-1,0,0,2d0)
         call coll_gq(1,3,2,xmsq,0,1,1,0,aveqq/aveqg)
         call coll_gq(2,4,1,xmsq,1,0,1,0,aveqq/aveqg)
         call msqupdate(j,k,msq,xmsq,1d0)
         call dip2jet_qrqr(1,4,2,3,5,xmsq,1,-1,2,-2)
         call coll_gq(3,4,1,xmsq,1,-1,0,0,2d0)
         call msqupdate(j,k,msq,xmsq,dfloat(nf-1))
         call dip2jet_qqgg(2,1,3,4,5,xmsq,1,-1,0,0,half)
         call msqupdate(j,k,msq,xmsq,third)
         call dip2jet_qqgg(2,1,3,5,4,xmsq,1,-1,0,0,half)
         call msqupdate(j,k,msq,xmsq,third)
         call dip2jet_qqgg(2,1,4,5,3,xmsq,1,-1,0,0,half)
         call msqupdate(j,k,msq,xmsq,third)
c---     (non-identical)
        else
         call dip2jet_qrqr(1,4,3,2,5,xmsq,1,2,1,2)
         call coll_gq(1,3,2,xmsq,0,1,1,0,aveqq/aveqg)
         call coll_gq(2,4,1,xmsq,1,0,1,0,aveqq/aveqg)
         call msqupdate(j,k,msq,xmsq,1d0)
        endif

c--- antiquark-quark
      elseif ((j .lt. 0) .and. (k. gt. 0)) then
c---     (identical)
        if (j .eq. -k) then
         call dip2jet_qqqq(2,3,4,1,5,xmsq,2,-2,2,-2)
         call coll_gq(3,4,1,xmsq,-1,1,0,0,2d0)
         call coll_gq(1,3,2,xmsq,0,1,1,0,aveqq/aveqg)
         call coll_gq(2,4,1,xmsq,1,0,1,0,aveqq/aveqg)
         call msqupdate(j,k,msq,xmsq,1d0)
         call dip2jet_qrqr(2,3,1,4,5,xmsq,-1,1,-2,2)
         call coll_gq(3,4,1,xmsq,-1,1,0,0,2d0)
         call msqupdate(j,k,msq,xmsq,dfloat(nf-1))
         call dip2jet_qqgg(1,2,3,4,5,xmsq,-1,1,0,0,half)
         call msqupdate(j,k,msq,xmsq,third)
         call dip2jet_qqgg(1,2,3,5,4,xmsq,-1,1,0,0,half)
         call msqupdate(j,k,msq,xmsq,third)
         call dip2jet_qqgg(1,2,4,5,3,xmsq,-1,1,0,0,half)
         call msqupdate(j,k,msq,xmsq,third)
c---     (non-identical)
        else
         call dip2jet_qrqr(2,3,4,1,5,xmsq,1,2,1,2)
         call coll_gq(1,3,2,xmsq,0,1,1,0,aveqq/aveqg)
         call coll_gq(2,4,1,xmsq,1,0,1,0,aveqq/aveqg)
         call msqupdate(j,k,msq,xmsq,1d0)
        endif

c--- quark-gluon
      elseif ((j .gt. 0) .and. (k. eq. 0)) then
        call dip2jet_qqgg(1,3,2,4,5,xmsq,1,0,1,0,half)
        call msqupdate(j,k,msq,xmsq,half)
        call dip2jet_qqgg(1,3,2,5,4,xmsq,1,0,1,0,half)
        call coll_gq(1,3,2,xmsq,0,0,0,0,2d0*aveqg/avegg)
        call coll_qg(2,3,1,xmsq,1,-1,0,0,
     .               2d0*cf*2d0*aveqg/aveqq) ! end of qqgg
        call coll_gq(1,3,2,xmsq,0,0,1,-1,2d0*dfloat(nf-1)*aveqg/avegg)
        call coll_gq(4,5,1,xmsq,1,0,1,0,2d0*dfloat(nf-1))
        call coll_qg(2,5,1,xmsq,2,1,2,1,
     .               2d0*2d0*cf*aveqg/aveqq)
        call coll_qg(2,3,1,xmsq,1,-1,2,-2,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq)
        call coll_qg(2,4,1,xmsq,2,-1,2,-1,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq) ! end of qrqr
        call coll_gq(1,3,2,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(1,4,2,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(3,5,1,xmsq,1,0,0,1,1d0)
        call coll_gq(4,5,1,xmsq,1,0,1,0,1d0)
        call coll_qg(2,5,1,xmsq,1,1,1,1,
     .               2d0*2d0*cf*aveqg/aveqq)
        call coll_qg(2,3,1,xmsq,1,-1,1,-1,
     .               2d0*cf*aveqg/aveqq)
        call coll_qg(2,4,1,xmsq,1,-1,1,-1,
     .               2d0*cf*aveqg/aveqq) ! end of qqqq
        call msqupdate(j,k,msq,xmsq,half)

c--- antiquark-gluon
      elseif ((j .lt. 0) .and. (k. eq. 0)) then
        call dip2jet_qqgg(3,1,2,4,5,xmsq,-1,0,-1,0,half)
        call msqupdate(j,k,msq,xmsq,half)
        call dip2jet_qqgg(3,1,2,5,4,xmsq,-1,0,-1,0,half)
        call coll_gq(1,3,2,xmsq,0,0,0,0,2d0*aveqg/avegg)
        call coll_qg(2,3,1,xmsq,-1,1,0,0,
     .               2d0*cf*2d0*aveqg/aveqq) ! end of qqgg
        call coll_gq(1,5,2,xmsq,0,0,1,-1,2d0*dfloat(nf-1)*aveqg/avegg)
        call coll_gq(3,4,1,xmsq,-1,0,0,-1,2d0*dfloat(nf-1))
        call coll_qg(2,5,1,xmsq,1,-1,2,-2,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq)
        call coll_qg(2,4,1,xmsq,-2,1,1,-2,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq) 
        call coll_qg(2,3,1,xmsq,-2,-1,-1,-2,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq) ! end of qrqr
        call coll_gq(1,5,2,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(1,4,2,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(3,5,1,xmsq,-1,0,0,-1,1d0)
        call coll_gq(3,4,1,xmsq,-1,0,0,-1,1d0)
        call coll_qg(2,5,1,xmsq,-1,1,1,-1,
     .               2d0*cf*aveqg/aveqq)
        call coll_qg(2,3,1,xmsq,-1,-1,-1,-1,
     .               2d0*2d0*cf*aveqg/aveqq)
        call coll_qg(2,4,1,xmsq,-1,1,1,-1,
     .               2d0*cf*aveqg/aveqq) ! end of qqqq
        call msqupdate(j,k,msq,xmsq,half)

c--- gluon-quark
      elseif ((j .eq. 0) .and. (k. gt. 0)) then
        call dip2jet_qqgg(2,3,1,4,5,xmsq,0,1,1,0,half)
        call msqupdate(j,k,msq,xmsq,half)
        call dip2jet_qqgg(2,3,1,5,4,xmsq,0,1,1,0,half)
        call coll_gq(2,3,1,xmsq,0,0,0,0,2d0*aveqg/avegg)
        call coll_qg(1,3,2,xmsq,1,-1,0,0,
     .               2d0*cf*2d0*aveqg/aveqq) ! end of qqgg
        call coll_gq(2,3,1,xmsq,0,0,1,-1,2d0*dfloat(nf-1)*aveqg/avegg)
        call coll_gq(4,5,1,xmsq,0,1,1,0,2d0*dfloat(nf-1))
        call coll_qg(1,5,2,xmsq,1,2,2,1,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq)
        call coll_qg(1,3,2,xmsq,-2,2,1,-1,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq)
        call coll_qg(1,4,2,xmsq,-1,2,2,-1,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq) ! end of qrqr
        call coll_gq(2,3,1,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(2,4,1,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(3,5,1,xmsq,0,1,0,1,1d0)
        call coll_gq(4,5,1,xmsq,0,1,1,0,1d0)
        call coll_qg(1,5,2,xmsq,1,1,1,1,
     .               2d0*2d0*cf*aveqg/aveqq)
        call coll_qg(1,3,2,xmsq,-1,1,1,-1,
     .               2d0*cf*aveqg/aveqq)
        call coll_qg(1,4,2,xmsq,-1,1,1,-1,
     .               2d0*cf*aveqg/aveqq) ! end of qqqq
        call msqupdate(j,k,msq,xmsq,half)

c--- gluon-antiquark
      elseif ((j .eq. 0) .and. (k. lt. 0)) then
        call dip2jet_qqgg(3,2,1,4,5,xmsq,0,-1,-1,0,half)
        call msqupdate(j,k,msq,xmsq,half)
        call dip2jet_qqgg(3,2,1,5,4,xmsq,0,-1,-1,0,half)
        call coll_gq(2,3,1,xmsq,0,0,0,0,2d0*aveqg/avegg)
        call coll_qg(1,3,2,xmsq,-1,1,0,0,
     .               2d0*cf*2d0*aveqg/aveqq) ! end of qqgg
        call coll_gq(2,5,1,xmsq,0,0,1,-1,2d0*dfloat(nf-1)*aveqg/avegg)
        call coll_gq(3,4,1,xmsq,0,-1,0,-1,2d0*dfloat(nf-1))
        call coll_qg(1,5,2,xmsq,2,-2,1,-1,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq)
        call coll_qg(1,4,2,xmsq,1,-2,1,-2,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq) 
        call coll_qg(1,3,2,xmsq,-1,-2,-1,-2,
     .               2d0*cf*2d0*dfloat(nf-1)*aveqg/aveqq) ! end of qrqr
        call coll_gq(2,5,1,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(2,4,1,xmsq,0,0,1,-1,aveqg/avegg)
        call coll_gq(3,5,1,xmsq,0,-1,0,-1,1d0)
        call coll_gq(3,4,1,xmsq,0,-1,0,-1,1d0)
        call coll_qg(1,5,2,xmsq,1,-1,1,-1,
     .               2d0*cf*aveqg/aveqq)
        call coll_qg(1,3,2,xmsq,-1,-1,-1,-1,
     .               2d0*2d0*cf*aveqg/aveqq)
        call coll_qg(1,4,2,xmsq,1,-1,1,-1,
     .               2d0*cf*aveqg/aveqq) ! end of qqqq
        call msqupdate(j,k,msq,xmsq,half)

      endif

      enddo
      enddo                    

c--- the other flavour combinations are easy now
      do j=-nf,nf
      do k=-nf,nf
        if ((abs(j) .gt. 2) .or. (abs(k) .gt. 2)) then
          do nd=1,ndmax
            msq(nd,j,k)=msq(nd,m(j),m(k))
          enddo
        endif
      enddo
      enddo

      return
      end
      
      
c--- incorporates the dipole contribution from a single sub-process
c--- into the total for flavours (j,k), weighted by a given factor, fac
      subroutine msqupdate(j,k,msq,xmsq,fac)
      implicit none
      include 'constants.f'
      include 'ptilde.f'
      double precision msq(maxd,-nf:nf,-nf:nf),xmsq(36),fac
      integer n,j,k 
      
      do n=1,36
        msq(n,j,k)=msq(n,j,k)+xmsq(n)*fac
      enddo
      
      return
      end      
      
      
c      subroutine qqb_2jet_gvecx(p,n,in,msqv,mvg,mvxg)
c      implicit none
c      include 'constants.f'
c      integer j,k,in
c      double precision msqv(-nf:nf,-nf:nf),p(mxpart,4),n(4)
c      double precision mvg(0:2,-nf:nf,-nf:nf)                                        
c      double precision mvxg(-nf:nf,-nf:nf,-nf:nf,-nf:nf)                        
c      do j=-nf,nf
c         do k=-nf,nf
c           msqv(j,k)=0d0
c         enddo
c      enddo
c      return
c      end                                                                       
      
c      subroutine donothing_gvecx(p,n,in,msq,mvg,mvxg)
c      implicit none
c      include 'constants.f'
c      integer j,k,in
c      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),n(4)
c      double precision mvg(0:2,-nf:nf,-nf:nf)                                        
c      double precision mvxg(-nf:nf,-nf:nf,-nf:nf,-nf:nf)                        
c      do j=-nf,nf
c         do k=-nf,nf
c            msq(j,k)=0d0
c         enddo
c      enddo
c      return
c      end                                                                       

