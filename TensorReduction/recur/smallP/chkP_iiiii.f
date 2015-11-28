      subroutine chkP_iiiii(k,i1,i2,i3,i4,i5,f,Gr,Shat6,N0)
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,i4,i5,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat6(np,z5max,-2:0),diff
       
      do ep=-2,0
      diff=Dv(diiiii(z5(i1,i2,i3,i4,i5))+N0,ep)*f(k)-
     . (Shat6(k,z5(i1,i2,i3,i4,i5),ep)
     . -2d0*delta(k,i1)*Dv(dzziiii(z4(i2,i3,i4,i5))+N0,ep)
     . -2d0*delta(k,i2)*Dv(dzziiii(z4(i1,i3,i4,i5))+N0,ep)
     . -2d0*delta(k,i3)*Dv(dzziiii(z4(i1,i2,i4,i5))+N0,ep)
     . -2d0*delta(k,i4)*Dv(dzziiii(z4(i1,i2,i3,i5))+N0,ep)
     . -2d0*delta(k,i5)*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep)
     . -Gr(k,1)*Dv(diiiiii(z6(1,i1,i2,i3,i4,i5))+N0,ep) 
     . -Gr(k,2)*Dv(diiiiii(z6(2,i1,i2,i3,i4,i5))+N0,ep)
     . -Gr(k,3)*Dv(diiiiii(z6(3,i1,i2,i3,i4,i5))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_iiiii',k,i1,i2,i3,i4,i5,diff
     
      enddo
      
      return
      end
