      subroutine chkG_0000(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shati00,Si00m,S0000,N0)
      implicit none
C-----DD 5.46
      include 'constants.f'  
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,k,l,n,m,np
      parameter(np=3)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),
     . f(np)
      double complex Shati00(np,-2:0),S0000(-2:0),Si00m(np,np,-2:0),
     . bit,pole,diff

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shati00(m,ep)-f(n)*f(m)*Dv(dd00+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Dv(dd0000+N0,ep-1)

      diff=
     . 8*Gtwiddle(k,l)*Dv(dd0000+N0,ep)
     . +pole
     . +DetGr*Dv(dzzii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S0000(ep)
     . -Gtwiddle(1,l)*Si00m(1,k,ep)
     . -Gtwiddle(2,l)*Si00m(2,k,ep)
     . -Gtwiddle(3,l)*Si00m(3,k,ep)
     . +Gtwiddle(k,l)*(Si00m(1,1,ep)+Si00m(2,2,ep)+Si00m(3,3,ep))
     . +bit
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.))
     . write(6,*) 'chkG_0000',k,l,diff
      enddo


      return
      end
  



