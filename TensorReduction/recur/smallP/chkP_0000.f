      subroutine chkP_0000(Gr,S0000,N0)
      implicit none
      include 'constants.f'
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex S0000(-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(dzzii(z2(n,m))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dd0000+N0,ep-1)

c--- note: we have simplified the recursion relation by using DD (5.9)
c--- so that S0000 = 2C00(0)+2*m0sq*D00      
      diff=Dv(dd0000+N0,ep)*12d0-
     . (pole
     . +S0000(ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_0000',diff
     
      enddo
      
      return
      end
