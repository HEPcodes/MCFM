      subroutine qqb_zbb(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> b(p5)+bb(p6)+e^-(p3)+e^+(p4)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      include 'zcouple.f'
c      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'hardscale.f'
      integer j,k,nu,ics
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),msqzbb,mmsq(2,2),
     . prop,qqb,qbq,pswap(mxpart,4)
      double precision msq_cs(0:2,-nf:nf,-nf:nf),mmsq_cs(0:2,2,2)
      common/msq_cols/msq_cs,mmsq_cs

c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(6,nu)
      pswap(5,nu)=p(3,nu)
      pswap(6,nu)=p(4,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_alt(mmsq)
      
C---Fill spinor products
      call spinoru(6,p,za,zb)
      prop=s(3,4)/sqrt((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)

c ensure that we have a hard process
      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return

      qqb=msqzbb(1,2,5,6)
      qbq=msqzbb(2,1,5,6)

      do j=-nf,nf
      k=-j
          if ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=
     .      +(Q(1)*q1+L(1)*l1*prop)**2*mmsq(1,1)
     .      +(Q(1)*q1+R(1)*l1*prop)**2*mmsq(2,1)
     .      +(Q(1)*q1+L(1)*r1*prop)**2*mmsq(1,2)
     .      +(Q(1)*q1+R(1)*r1*prop)**2*mmsq(2,2)
            do ics=0,2
            msq_cs(ics,j,k)=
     .      +(Q(1)*q1+L(1)*l1*prop)**2*mmsq_cs(ics,1,1)
     .      +(Q(1)*q1+R(1)*l1*prop)**2*mmsq_cs(ics,2,1)
     .      +(Q(1)*q1+L(1)*r1*prop)**2*mmsq_cs(ics,1,2)
     .      +(Q(1)*q1+R(1)*r1*prop)**2*mmsq_cs(ics,2,2)
            enddo
          elseif ((j .gt. 0) .and. (k .lt. 0)) then
       msq(j,k)=((Q(j)*q1+L(j)*l1*prop)**2
     .          +(Q(j)*q1+R(j)*r1*prop)**2)*qqb
     .         +((Q(j)*q1+L(j)*r1*prop)**2
     .          +(Q(j)*q1+R(j)*l1*prop)**2)*qbq
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
       msq(j,k)=((Q(k)*q1+L(k)*l1*prop)**2
     .          +(Q(k)*q1+R(k)*r1*prop)**2)*qbq
     .         +((Q(k)*q1+L(k)*r1*prop)**2
     .          +(Q(k)*q1+R(k)*l1*prop)**2)*qqb
          endif
      enddo
      return
      end







