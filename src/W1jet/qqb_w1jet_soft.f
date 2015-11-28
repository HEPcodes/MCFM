      subroutine qqb_w1jet_soft(P,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> W +g(p5) +g(p6)
c                          |
c                          --> nu(p3)+e^+(p4)
c                            
c   positively charged W only with p6 soft
c--all momenta incoming
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      integer j,k
      logical hhhh
      common/hhhh/hhhh
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf),
     . msq0(-nf:nf,-nf:nf),s(mxpart,mxpart),
     . facqqb,facqbq,facgqb,facqbg,facgq,facqg
      double precision eik12,eik15,eik25
      double precision eik16_2,eik26_1,eik16_5,eik56_1
      double precision eik26_5,eik56_2
      
      call dotem(6,p,s)

      eik12=s(1,2)/(s(1,6)*s(2,6))
      eik15=s(1,5)/(s(1,6)*s(5,6))
      eik25=s(2,5)/(s(2,6)*s(5,6))


      eik16_2=two*gsq*s(1,2)/(s(1,6)+s(2,6))/s(1,6)
      eik26_1=two*gsq*s(1,2)/(s(1,6)+s(2,6))/s(2,6)
      eik16_5=two*gsq*s(1,5)/(s(5,6)-s(1,6))/s(1,6)
      eik56_1=-two*gsq*s(1,5)/(s(5,6)-s(1,6))/s(5,6)
      eik26_5=two*gsq*s(2,5)/(s(5,6)-s(2,6))/s(2,6)
      eik56_2=-two*gsq*s(2,5)/(s(5,6)-s(2,6))/s(5,6)

c      write(6,*) 'eik16_2',eik16_2
c      write(6,*) 'eik26_1',eik26_1
c      write(6,*) 'eik16_5',eik16_5
c      write(6,*) 'eik56_1',eik56_1
c      write(6,*) 'eik26_5',eik26_5
c      write(6,*) 'eik56_2',eik56_2

c      write(6,*) 'eik12',eik16_2+eik26_1
c      write(6,*) 'eik15',eik16_5+eik56_1
c      write(6,*) 'eik25',eik26_5+eik56_2

      facqqb=two*gsq*(xn*(eik15+eik25)-eik12/xn)
      facqg=two*gsq*(xn*(eik12+eik25)-eik15/xn)
      facgq=two*gsq*(xn*(eik12+eik15)-eik25/xn)
      facqbq=facqqb
      facqbg=facqg
      facgqb=facgq

c--- SUB-LEADING ONLY
c      facqqb=two*gsq*(-eik12/xn)
c      facqbq=facqbq
c      facqg=two*gsq*(-eik15/xn)
c      facgq=two*gsq*(-eik25/xn)
c      facgqb=two*gsq*(-eik25/xn)
c      facqbg=two*gsq*(-eik15/xn)

      call qqb_w_g(p,msq0)
 

c      write(6,*) 
c      write(6,*) 'In soft'
c      write(6,*) 'msq0(-1,0)',msq0(-1,0)
      
c--debug
      do j=-nf,nf
      do k=-nf,nf
            if ((j .gt. 0) .and. (k.lt.0)) then
               msq(j,k)=facqqb*msq0(j,k)
c               write(6,*) 'eik16_2*msq0(j,k)/xn',-eik16_2*msq0(j,k)/xn
c               write(6,*) 'eik26_1*msq0(j,k)/xn',-eik26_1*msq0(j,k)/xn
c               write(6,*) 'eik16_5*msq0(j,k)*xn',eik16_5*msq0(j,k)*xn
c               write(6,*) 'eik56_1*msq0(j,k)*xn',eik56_1*msq0(j,k)*xn

c               write(6,*) 'eik16_5*msq0(j,k)*xn+eik56_1*msq0(j,k)*xn',
c     . eik16_5*msq0(j,k)*xn+eik56_1*msq0(j,k)*xn

c               write(6,*) 'eik56_2*msq0(j,k)*xn',eik56_2*msq0(j,k)*xn
c               write(6,*) 'eik26_5*msq0(j,k)*xn',eik26_5*msq0(j,k)*xn

c               write(6,*) 'eik56_2*msq0(j,k)*xn+eik26_5*msq0(j,k)*xn',
c     .         eik56_2*msq0(j,k)*xn+eik26_5*msq0(j,k)*xn
c      pause
            elseif ((j .lt. 0) .and. (k.gt.0)) then
               msq(j,k)=facqbq*msq0(j,k)
            elseif ((j .gt. 0) .and. (k.eq.0)) then
               msq(j,k)=facqg*msq0(j,k)
            elseif ((j .lt. 0) .and. (k.eq.0)) then
               msq(j,k)=facqbg*msq0(j,k)
            elseif ((j .eq. 0) .and. (k.gt.0)) then
               msq(j,k)=facgq*msq0(j,k)
            elseif ((j .eq. 0) .and. (k.lt.0)) then
               msq(j,k)=facgqb*msq0(j,k)
            endif
      enddo
      enddo
      return
      end

 
