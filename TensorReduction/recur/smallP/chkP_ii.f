      subroutine chkP_ii(k,i1,i2,f,Gr,Shat3,N0)
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat3(np,z2max,-2:0),diff
       
      do ep=-2,0
      diff=Dv(dii(z2(i1,i2))+N0,ep)*f(k)-
     . (Shat3(k,z2(i1,i2),ep)
     . -2d0*delta(k,i1)*Dv(dzzi(i2)+N0,ep)
     . -2d0*delta(k,i2)*Dv(dzzi(i1)+N0,ep)
     . -Gr(k,1)*Dv(diii(z3(1,i1,i2))+N0,ep) 
     . -Gr(k,2)*Dv(diii(z3(2,i1,i2))+N0,ep)
     . -Gr(k,3)*Dv(diii(z3(3,i1,i2))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_ii',k,i1,i2,diff
     
      enddo
      
      return
      end
