      subroutine chkG_00(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat1,Shat2,S00,N0)
      implicit none
C-----DD Eq. 5.42
      include 'constants.f' 
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,n,m,np
      parameter(np=3)
      double precision DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      double complex S00(-2:0),Shat1(np,-2:0),Shat2(np,np,-2:0),bit,diff

      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat1(m,ep)-f(n)*f(m)*Dv(dd0+N0,ep))
      enddo
      enddo
      do n=1,np
      bit=bit
      enddo
      diff=
     . 4*Gtwiddle(k,l)*Dv(dd00+N0,ep)
     . +DetGr*Dv(dii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S00(ep)
     . +Gtwiddle(k,l)*Shat2(1,1,ep)-Gtwiddle(1,l)*Shat2(1,k,ep)
     . +Gtwiddle(k,l)*Shat2(2,2,ep)-Gtwiddle(2,l)*Shat2(2,k,ep)
     . +Gtwiddle(k,l)*Shat2(3,3,ep)-Gtwiddle(3,l)*Shat2(3,k,ep)
     . +bit

      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.)) 
     . write(6,*) 'chkG_00',k,l,diff
      enddo


      return
      end
  



