      subroutine chkCF_00i(i1,f,Gr,Shat3,N0)
C---  Expression for rearrangement of Eq. 5.68
C---  Checks C00i
C---  Small terms of order f(i)*Cij,Gr(i,j)*Cijk
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat3(np,z2max,-2:0),diff
       
      do ep=-2,0
      diff=Cv(czzi(i1)+N0,ep)*4d0-
     . (Shat3(i1,z2(i1,i1),ep)
     . -f(i1)*Cv(cii(z2(i1,i1))+N0,ep)
     . -Gr(i1,1)*Cv(ciii(z3(1,i1,i1))+N0,ep) 
     . -Gr(i1,2)*Cv(ciii(z3(2,i1,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCF_00i',i1,diff
     
      enddo
      
      return
      end
