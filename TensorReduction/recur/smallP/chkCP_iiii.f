      subroutine chkCP_iiii(k,i1,i2,i3,i4,f,Gr,Shat5,N0)
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,i4,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat5(np,z4max,-2:0),diff
       
      do ep=-2,0
      diff=Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)*f(k)-
     . (Shat5(k,z4(i1,i2,i3,i4),ep)
     . -2d0*delta(k,i1)*Cv(czziii(z3(i2,i3,i4))+N0,ep)
     . -2d0*delta(k,i2)*Cv(czziii(z3(i1,i3,i4))+N0,ep)
     . -2d0*delta(k,i3)*Cv(czziii(z3(i1,i2,i4))+N0,ep)
     . -2d0*delta(k,i4)*Cv(czziii(z3(i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Cv(ciiiii(z5(1,i1,i2,i3,i4))+N0,ep) 
     . -Gr(k,2)*Cv(ciiiii(z5(2,i1,i2,i3,i4))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_iiii',k,i1,i2,i3,i4,diff
     
      enddo
      
      return
      end
