      subroutine chkCP_i(k,i1,f,Gr,Shat2,N0)
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat2(np,np,-2:0),diff
       
      do ep=-2,0
      diff=Cv(ci(i1)+N0,ep)*f(k)-
     . (Shat2(k,i1,ep)
     . -2d0*delta(k,i1)*Cv(cc00+N0,ep)
     . -Gr(k,1)*Cv(cii(z2(1,i1))+N0,ep) 
     . -Gr(k,2)*Cv(cii(z2(2,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_i',k,i1,diff
     
      enddo
      
      return
      end
