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
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer j,k,nu,ics
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),mmsq(2,2),
     . pswap(mxpart,4)
      double precision msq_cs(0:2,-nf:nf,-nf:nf),mmsq_cs(0:2,2,2),faclo
      double complex qqb1,qbq1,qqb2,qbq2,qqb3,qbq3,qqb4,qbq4,tamp,prop
      double complex qqb5,qbq5,qqb6,qbq6,qqb7,qbq7,qqb8,qbq8
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
      call xzqqgg(mmsq)
      
C---Fill spinor products
      call spinoru(6,p,za,zb)
      prop=s(3,4)/(s(3,4)-zmass**2+im*zmass*zwidth)

c ensure that we have a hard process
      if (
     .      (s(5,6) .lt. four*hscalesq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. hscalesq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. hscalesq) ) return

c--- qqb
      call msqzbb(1,2,5,6,qqb1,qqb2,qqb3,qqb4,qqb5,qqb6,qqb7,qqb8)
c--- qbq from symmetries
      qbq1 = -qqb5
      qbq2 = +qqb6
      qbq3 = -qqb7
      qbq4 = +qqb8
      qbq5 = -qqb1
      qbq6 = +qqb2
      qbq7 = -qqb3
      qbq8 = +qqb4

      faclo=4d0*V*gsq**2*esq**2*aveqq 

      do j=-nf,nf
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
            tamp=(Q(j)*q1+L(j)*l1*prop)*qqb1
     .          +(Q(1)*q1+L(1)*l1*prop)*qqb2
            msq(j,k)=faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*l1*prop)*qqb3
     .          +(Q(1)*q1+R(1)*l1*prop)*qqb4
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*r1*Dconjg(prop))*qqb7
     .          +(Q(1)*q1+L(1)*r1*Dconjg(prop))*qqb8
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+L(j)*r1*Dconjg(prop))*qqb5
     .          +(Q(1)*q1+R(1)*r1*Dconjg(prop))*qqb6
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*l1*prop)*qqb5
     .          +(Q(1)*q1+L(1)*l1*prop)*qqb6
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*l1*prop)*qqb7
     .          +(Q(1)*q1+R(1)*l1*prop)*qqb8
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*r1*Dconjg(prop))*qqb3
     .          +(Q(1)*q1+L(1)*r1*Dconjg(prop))*qqb4
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(j)*q1+R(j)*r1*Dconjg(prop))*qqb1
     .          +(Q(1)*q1+R(1)*r1*Dconjg(prop))*qqb2
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
          elseif ((j .lt. 0) .and. (k. gt. 0)) then
            tamp=(Q(k)*q1+L(k)*l1*prop)*qbq1
     .          +(Q(1)*q1+L(1)*l1*prop)*qbq2
            msq(j,k)=faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*l1*prop)*qbq3
     .          +(Q(1)*q1+R(1)*l1*prop)*qbq4
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*r1*Dconjg(prop))*qbq7
     .          +(Q(1)*q1+L(1)*r1*Dconjg(prop))*qbq8
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+L(k)*r1*Dconjg(prop))*qbq5
     .          +(Q(1)*q1+R(1)*r1*Dconjg(prop))*qbq6
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*l1*prop)*qbq5
     .          +(Q(1)*q1+L(1)*l1*prop)*qbq6
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*l1*prop)*qbq7
     .          +(Q(1)*q1+R(1)*l1*prop)*qbq8
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*r1*Dconjg(prop))*qbq3
     .          +(Q(1)*q1+L(1)*r1*Dconjg(prop))*qbq4
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
            tamp=(Q(k)*q1+R(k)*r1*Dconjg(prop))*qbq1
     .          +(Q(1)*q1+R(1)*r1*Dconjg(prop))*qbq2
            msq(j,k)=msq(j,k)+faclo*abs(tamp)**2
          endif
      enddo
      return
      end





      subroutine msqzbb(i1,i2,i5,i6,amp1,amp2,amp3,amp4,
     .                                    amp5,amp6,amp7,amp8)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer i1,i2,i5,i6
      double complex aqqb_zbb_new,amp1,amp2,amp3,amp4
      double complex amp5,amp6,amp7,amp8

c--- also include diagrams where the Z is attached to b-bbar line  
c--- we swap 1 <--> 5, 2 <-->6 to obtain the correct diagram
c--- parity (charge conjugation) fixes helicities for the quarks
c--- and lepton helicity is recovered by 3 <--> 4
c--- overall minus sign for parity
c---
c--- notation: (qqb hel, bbb hel) - lepton helicity (p4) is right-handed
c--- Z to qqb, L L
      amp1=+aqqb_zbb_new(i1,i6,i5,i2,3,4)
c--- Z to bbb, L L
      amp2=-Dconjg(aqqb_zbb_new(i5,i2,i1,i6,4,3))
c--- Z to qqb, L R
      amp3=-aqqb_zbb_new(i1,i5,i6,i2,3,4)
c--- Z to bbb, L R
      amp4=+(aqqb_zbb_new(i5,i1,i2,i6,3,4))
c--- Z to qqb, R L
      amp5=+Dconjg(aqqb_zbb_new(i1,i5,i6,i2,4,3))    
c--- Z to bbb, R L
      amp6=-Dconjg(aqqb_zbb_new(i5,i1,i2,i6,4,3))    
c--- Z to qqb, R R
      amp7=-Dconjg(aqqb_zbb_new(i1,i6,i5,i2,4,3))
c--- Z to bbb, R R
      amp8=+aqqb_zbb_new(i5,i2,i1,i6,3,4)
      
c--- note that qbq amplitudes can be obtained from these by the
c--- following symmetries, since 1 <--> 2 is equivalent to
c--- parity (c.c.) followed by 3 <--> 4 and 5 <--> 6
c--- qbq1 = -qqb5
c--- qbq2 = +qqb6
c--- qbq3 = -qqb7
c--- qbq4 = +qqb8
c--- qbq5 = -qqb1
c--- qbq6 = +qqb2
c--- qbq7 = -qqb3
c--- qbq8 = +qqb4
      
      return 
      end


