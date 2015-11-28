      subroutine gg_hgg_soft(p,msq)
      implicit none
c---- matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H(-->  b(p3)+b~(p4))+g(p5)
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j,k,iglue
      double precision msq(-nf:nf,-nf:nf),msqskip(fn:nf,fn:nf),
     . p(mxpart,4),pskip(mxpart,4)
      double precision sh,ss,tt,uu,hdecay,s(mxpart,mxpart)
      double precision qa,qg,gq,gg,s34,Asq
      double precision 
     . eik16_2,eik16_5,eik26_5,
     . eik26_1,eik56_1,eik56_2,
     . eik15_2,eik15_6,eik25_6,
     . eik25_1,eik65_1,eik65_2

      parameter(iglue=5)

      do k=1,4
      do j=1,4
      pskip(j,k)=p(j,k)
      enddo
      pskip(5,k)=p(6,k)
      enddo

      call gg_hg(pskip,msqskip)
      call gg_hg(p,msq)

      call dotem(iglue+1,p,s)

C   Deal withb Higgs decay to b-bbar
      s34=s(3,4)+2d0*mb**2
      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s34-4d0*mb**2) 
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)


      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)
      sh=ss+tt+uu
      Asq=(as/(3d0*pi))**2/vevsq

      eik16_2=+two*gsq*s(1,2)/(s(1,6)+s(2,6))/s(1,6)
      eik16_5=+two*gsq*s(1,5)/(-s(1,6)+s(5,6))/s(1,6)
      eik26_5=+two*gsq*s(2,5)/(-s(2,6)+s(5,6))/s(2,6)

      eik26_1=+two*gsq*s(1,2)/(s(1,6)+s(2,6))/s(2,6)
      eik56_1=-two*gsq*s(1,5)/(-s(1,6)+s(5,6))/s(5,6)
      eik56_2=-two*gsq*s(2,5)/(-s(2,6)+s(5,6))/s(5,6)

      eik15_2=+two*gsq*s(1,2)/(s(1,5)+s(2,5))/s(1,5)
      eik15_6=+two*gsq*s(1,6)/(-s(1,5)+s(5,6))/s(1,5)
      eik25_6=+two*gsq*s(2,6)/(-s(2,5)+s(5,6))/s(2,5)

      eik25_1=+two*gsq*s(1,2)/(s(1,5)+s(2,5))/s(2,5)
      eik65_1=-two*gsq*s(1,6)/(-s(1,5)+s(5,6))/s(6,5)
      eik65_2=-two*gsq*s(2,6)/(-s(2,5)+s(5,6))/s(6,5)

      gg=
     .  msq(0,0)*half*xn
     . *(eik16_2+eik26_1+eik16_5+eik56_1+eik26_5+eik56_2)
     . +msqskip(0,0)*half*xn
     . *(eik15_2+eik25_1+eik15_6+eik65_1+eik25_6+eik65_2)
      qa=
     . +msq(1,-1)*half*(
     . +xn*(eik16_5+eik56_1+eik26_5+eik56_2)-(eik16_2+eik26_1)/xn)
     . +msqskip(1,-1)*half*(
     . +xn*(eik15_6+eik65_1+eik25_6+eik65_2)-(eik15_2+eik25_1)/xn)
      gq=msq(0,1)
     . *(xn*(eik56_1+eik16_5+eik26_1+eik16_2)-(eik56_2+eik26_5)/xn)
      qg=msq(1,0)
     . *(xn*(eik56_2+eik26_5+eik26_1+eik16_2)-(eik56_1+eik16_5)/xn)


      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initialize
      msq(j,k)=0d0
      if ((k .eq. -j) .and. (j .ne. 0)) then
      msq(j,k)=qa
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
      msq(j,k)=gq
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
      msq(j,k)=qg
      elseif ((k .eq. 0) .and. (j .eq. 0)) then
      msq(j,k)=gg
      endif
      enddo
      enddo

      return
      end
