      subroutine chkP_iii(k,i1,i2,i3,f,Gr,Shat4,N0)
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat4(np,z3max,-2:0),diff
       
      do ep=-2,0
      diff=Dv(diii(z3(i1,i2,i3))+N0,ep)*f(k)-
     . (Shat4(k,z3(i1,i2,i3),ep)
     . -2d0*delta(k,i1)*Dv(dzzii(z2(i2,i3))+N0,ep)
     . -2d0*delta(k,i2)*Dv(dzzii(z2(i1,i3))+N0,ep)
     . -2d0*delta(k,i3)*Dv(dzzii(z2(i1,i2))+N0,ep)
     . -Gr(k,1)*Dv(diiii(z4(1,i1,i2,i3))+N0,ep) 
     . -Gr(k,2)*Dv(diiii(z4(2,i1,i2,i3))+N0,ep)
     . -Gr(k,3)*Dv(diiii(z4(3,i1,i2,i3))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_iii',k,i1,i2,i3,diff
     
      enddo
      
      return
      end
