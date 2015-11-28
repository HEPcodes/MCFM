      subroutine chkP_0(k,f,Gr,Shat1,N0)
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat1(np,-2:0),diff
       
      do ep=-2,0
      diff=Dv(dd0+N0,ep)*f(k)-
     . (Shat1(k,ep)
     . -Gr(k,1)*Dv(di(1)+N0,ep) 
     . -Gr(k,2)*Dv(di(2)+N0,ep) 
     . -Gr(k,3)*Dv(di(3)+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_0',k,diff
     
      enddo
      
      return
      end
