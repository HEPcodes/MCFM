       subroutine chkC5_46(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat3zz,Shat4zz,S0000,N0)
      implicit none
      include 'constants.f'  
      include 'Cnames.f'  
      include 'Cv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      include 'weenumber.f'  
      integer ep,N0,k,l,n,m,np
      parameter(np=2)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),
     . f(np)
      double complex Shat3zz(np,-2:0),S0000(-2:0),
     . Shat4zz(np,z1max,-2:0),bit,pole,diff

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat3zz(m,ep)-f(n)*f(m)*Cv(cc00+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Cv(cc0000+N0,ep-1)

      diff=
     . 10*Gtwiddle(k,l)*Cv(cc0000+N0,ep)
     . +pole
     . +DetGr*Cv(czzii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S0000(ep)
     . -Gtwiddle(1,l)*Shat4zz(1,k,ep)
     . -Gtwiddle(2,l)*Shat4zz(2,k,ep)
     . +Gtwiddle(k,l)*(Shat4zz(1,1,ep)+Shat4zz(2,2,ep))
     . +bit
      if (abs(diff) .gt. weenumber) write(6,*) 'chkC5_46',k,l,diff
      enddo


      return
      end
  



