      subroutine chkCP_0(k,f,Gr,Shat1,N0)
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat1(np,-2:0),diff
       
      do ep=-2,0
      diff=Cv(cc0+N0,ep)*f(k)-
     . (Shat1(k,ep)
     . -Gr(k,1)*Cv(ci(1)+N0,ep) 
     . -Gr(k,2)*Cv(ci(2)+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_0',k,diff
     
      enddo
      
      return
      end
