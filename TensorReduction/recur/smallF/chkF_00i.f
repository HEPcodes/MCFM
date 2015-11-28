      subroutine chkF_00i(i1,f,Gr,Shat3,N0)
C---  Expression for rearrangement of Eq. 5.68
C---  Checks D00i
C---  Small terms of order f(i)*Dij,Gr(i,j)*Dijk
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat3(np,z2max,-2:0),diff
       
      do ep=-2,0
      diff=Dv(dzzi(i1)+N0,ep)*4d0-
     . (Shat3(i1,z2(i1,i1),ep)
     . -f(i1)*Dv(dii(z2(i1,i1))+N0,ep)
     . -Gr(i1,1)*Dv(diii(z3(1,i1,i1))+N0,ep) 
     . -Gr(i1,2)*Dv(diii(z3(2,i1,i1))+N0,ep)
     . -Gr(i1,3)*Dv(diii(z3(3,i1,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkF_00i',i1,diff
     
      enddo
      
      return
      end
