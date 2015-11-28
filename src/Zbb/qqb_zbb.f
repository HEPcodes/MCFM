      subroutine qqb_zbb(p,msq)
c---  Matrix elements squared 
c     q(-p1)+qb(-p2) --> bbar(p5)+b(p6)+e^-(p3)+e^+(p4)
c---  averaged(summed) over initial(final) colours and spins
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      include 'mmsq_cs.f'
      integer j,k,nu,ics,j1,j2,j3,swap(2)
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),mmsq(2,2),
     . pswap(mxpart,4),faclo
      double complex tamp,prop
c      double complex qqb5,qbq5,qqb6,qbq6,qqb7,qbq7,qqb8,qbq8
      double complex qqb_a(2,2,2),qqb_b(2,2,2)
      double complex qbq_a(2,2,2),qbq_b(2,2,2)
      data swap/2,1/
      save swap
      
c--initialize to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c ---Call the two gluon process which is defined in xzqqgg 
C ---in the notation
C     0 ---> q(p1)+g(p2)+g(p3)+qbar(p4)+a(p5)+  l(p6)
C ---compared with ours which is:-
c     0 ---> b(p6)+g(p1)+g(p2)+bb(p5)+e^+(p4)+e^-(p3)

      do nu=1,4
      pswap(1,nu)=p(6,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg(mmsq)
      
C---Fill spinor products
      call spinoru(6,p,za,zb)
      prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)

c ensure that we have a hard process
      if (
     .      (s(5,6) .lt. four*mbsq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. mbsq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. mbsq) ) return

c--- qqb
      call ampqqb_qqb(1,2,6,5,qqb_a,qqb_b)
C Instead of calling ampqqb_qqb(2,1,5,6,qbq_a,qbq_b)
c--- qbq from symmetries
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbq_a(j1,j2,j3)=-qqb_a(swap(j1),j2,j3)
      qbq_b(j1,j2,j3)=+qqb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

      faclo=4d0*V*gsq**2*esq**2*aveqq 

      do j=-(nf-1),(nf-1)
      k=-j
          if ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=
     .      +abs(Q(1)*q1+L(1)*l1*prop)**2*mmsq(1,1)
     .      +abs(Q(1)*q1+R(1)*l1*prop)**2*mmsq(2,1)
     .      +abs(Q(1)*q1+L(1)*r1*prop)**2*mmsq(1,2)
     .      +abs(Q(1)*q1+R(1)*r1*prop)**2*mmsq(2,2)
            do ics=0,2
            msq_cs(ics,j,k)=
     .      +abs(Q(1)*q1+L(1)*l1*prop)**2*mmsq_cs(ics,1,1)
     .      +abs(Q(1)*q1+R(1)*l1*prop)**2*mmsq_cs(ics,2,1)
     .      +abs(Q(1)*q1+L(1)*r1*prop)**2*mmsq_cs(ics,1,2)
     .      +abs(Q(1)*q1+R(1)*r1*prop)**2*mmsq_cs(ics,2,2)
            enddo
          elseif ((j .gt. 0) .and. (k .lt. 0)) then
            tamp=(Q(j)*q1+L(j)*l1*prop)*qqb_a(1,1,1)
     .          +(Q(1)*q1+L(1)*l1*prop)*qqb_b(1,1,1)
            msq(j,k)=faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(1)*q1+R(1)*l1*prop)*qqb_b(1,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*r1*prop)*qqb_a(1,1,2)
     .          +(Q(1)*q1+L(1)*r1*prop)*qqb_b(1,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(1)*q1+R(1)*r1*prop)*qqb_b(1,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(1)*q1+L(1)*l1*prop)*qqb_b(2,1,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*l1*prop)*qqb_a(2,2,1)
     .          +(Q(1)*q1+R(1)*l1*prop)*qqb_b(2,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(1)*q1+L(1)*r1*prop)*qqb_b(2,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*r1*prop)*qqb_a(2,2,2)
     .          +(Q(1)*q1+R(1)*r1*prop)*qqb_b(2,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            tamp=(Q(k)*q1+L(k)*l1*prop)*qbq_a(1,1,1)
     .          +(Q(1)*q1+L(1)*l1*prop)*qbq_b(1,1,1)
            msq(j,k)=faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(1)*q1+R(1)*l1*prop)*qbq_b(1,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*r1*prop)*qbq_a(1,1,2)
     .          +(Q(1)*q1+L(1)*r1*prop)*qbq_b(1,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(1)*q1+R(1)*r1*prop)*qbq_b(1,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(1)*q1+L(1)*l1*prop)*qbq_b(2,1,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*l1*prop)*qbq_a(2,2,1)
     .          +(Q(1)*q1+R(1)*l1*prop)*qbq_b(2,2,1)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(1)*q1+L(1)*r1*prop)*qbq_b(2,1,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*r1*prop)*qbq_a(2,2,2)
     .          +(Q(1)*q1+R(1)*r1*prop)*qbq_b(2,2,2)
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
          endif
      enddo
      return
      end





