      subroutine qqb_zbb_soft(P,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  bb(p5)+b(p6)+ Z +g(p7)
c                                         |
c                                          --> mu(p3)+nubar(p4)
c                            
c   positively charged W only
c--all momenta incoming
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      integer j,k

      double precision P(mxpart,4),msq(-nf:nf,-nf:nf),
     . msq0(-nf:nf,-nf:nf),s(mxpart,mxpart),facqqb,facqbq
      double precision eik17_2,eik17_6,eik17_5,eik27_6,eik27_5
      double precision eik27_1,eik67_1,eik57_1,eik67_2,eik57_2
      double precision eik67_5,eik57_6
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      call dotem(7,p,s)

      eik17_2=+s(1,2)/(s(1,7)+s(2,7))/s(1,7)*two*gsq
      eik27_1=+s(1,2)/(s(1,7)+s(2,7))/s(2,7)*two*gsq

      eik17_5=+s(1,5)/(-s(1,7)+s(5,7))/s(1,7)*two*gsq
      eik57_1=-s(1,5)/(-s(1,7)+s(5,7))/s(5,7)*two*gsq

      eik27_6=+s(2,6)/(-s(2,7)+s(6,7))/s(2,7)*two*gsq
      eik67_2=-s(2,6)/(-s(2,7)+s(6,7))/s(6,7)*two*gsq

      eik27_5=+s(2,5)/(-s(2,7)+s(5,7))/s(2,7)*two*gsq
      eik57_2=-s(2,5)/(-s(2,7)+s(5,7))/s(5,7)*two*gsq

      eik17_6=+s(1,6)/(-s(1,7)+s(6,7))/s(1,7)*two*gsq
      eik67_1=-s(1,6)/(-s(1,7)+s(6,7))/s(6,7)*two*gsq

      eik67_5=+s(6,5)/(s(6,7)+s(5,7))/s(6,7)*two*gsq
      eik57_6=+s(6,5)/(s(6,7)+s(5,7))/s(5,7)*two*gsq

      facqqb=(eik27_5+eik57_2+eik17_6+eik67_1)*(xn-two/xn)
     .      +(eik17_5+eik57_1+eik27_6+eik67_2)*(two/xn)
     .      +(eik27_1+eik17_2+eik67_5+eik57_6)*(-one/xn)
      facqbq=(eik17_5+eik57_1+eik27_6+eik67_2)*(xn-two/xn)
     .      +(eik27_5+eik57_2+eik17_6+eik67_1)*(two/xn)
     .      +(eik17_2+eik27_1+eik67_5+eik57_6)*(-one/xn)

      call qqb_zbb(p,msq0)

      
      do j=-nf,nf
      do k=-nf,nf
            if ((j .gt. 0) .and. (k.lt.0)) then
             msq(j,k)=facqqb*msq0(j,k)
            elseif ((j .lt. 0) .and. (k.gt.0)) then
               msq(j,k)=facqbq*msq0(j,k)
            elseif ((j .eq. 0) .and. (k .eq. 0)) then
             if (p(7,4) .lt. 2d0) then
          msq(j,k)=+(eik17_6+eik67_1)*msq_cs(1,j,k)*xn
     .             +(eik27_5+eik57_2)*msq_cs(1,j,k)*xn
     .             +(eik17_2+eik27_1)*msq_cs(1,j,k)*xn
     .             +(eik27_6+eik67_2)*msq_cs(2,j,k)*xn
     .             +(eik17_5+eik57_1)*msq_cs(2,j,k)*xn
     .             +(eik17_2+eik27_1)*msq_cs(2,j,k)*xn
c--- done with leading colour
     .             -(eik67_5+eik57_6)*msq_cs(1,j,k)/xn
     .             -(eik67_5+eik57_6)*msq_cs(2,j,k)/xn
     .             +(eik67_1+eik17_6)*msq_cs(0,j,k)*xn
     .             +(eik57_2+eik27_5)*msq_cs(0,j,k)*xn
     .             +(eik57_1+eik17_5)*msq_cs(0,j,k)*xn
     .             +(eik67_2+eik27_6)*msq_cs(0,j,k)*xn
c--- done with one gluon QED-like (1/N**2 colour-suppressed)
     .             -(eik67_5+eik57_6)*msq_cs(0,j,k)*(xn+one/xn)
c--- done with all gluons QED-like
             endif
            endif
      enddo
      enddo
      return
      end

 
