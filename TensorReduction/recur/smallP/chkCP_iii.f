      subroutine chkCP_iii(k,i1,i2,i3,f,Gr,Shat4,N0)
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat4(np,z3max,-2:0),diff
       
      do ep=-2,0
      diff=Cv(ciii(z3(i1,i2,i3))+N0,ep)*f(k)-
     . (Shat4(k,z3(i1,i2,i3),ep)
     . -2d0*delta(k,i1)*Cv(czzii(z2(i2,i3))+N0,ep)
     . -2d0*delta(k,i2)*Cv(czzii(z2(i1,i3))+N0,ep)
     . -2d0*delta(k,i3)*Cv(czzii(z2(i1,i2))+N0,ep)
     . -Gr(k,1)*Cv(ciiii(z4(1,i1,i2,i3))+N0,ep) 
     . -Gr(k,2)*Cv(ciiii(z4(2,i1,i2,i3))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_iii',k,i1,i2,i3,diff
     
      enddo
      
      return
      end
