      subroutine chkCP_ii(k,i1,i2,f,Gr,Shat3,N0)
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat3(np,z2max,-2:0),diff
       
      do ep=-2,0
      diff=Cv(cii(z2(i1,i2))+N0,ep)*f(k)-
     . (Shat3(k,z2(i1,i2),ep)
     . -2d0*delta(k,i1)*Cv(czzi(i2)+N0,ep)
     . -2d0*delta(k,i2)*Cv(czzi(i1)+N0,ep)
     . -Gr(k,1)*Cv(ciii(z3(1,i1,i2))+N0,ep) 
     . -Gr(k,2)*Cv(ciii(z3(2,i1,i2))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_ii',k,i1,i2,diff
     
      enddo
      
      return
      end
