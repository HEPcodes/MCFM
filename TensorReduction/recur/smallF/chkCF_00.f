      subroutine chkCF_00(i1,f,Gr,Shat2,N0)
C---  Expression for rearrangement of Eq. 5.66
C---  Checks C00
C---  Small terms of order f(i)*Ci,Gr(i,j)*Cij
      implicit none
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,np
      parameter(np=2)
      double precision f(np),Gr(np,np)
      double complex Shat2(np,np,-2:0),diff
       
      do ep=-2,0
      diff=Cv(cc00+N0,ep)*2d0-
     . (Shat2(i1,i1,ep)
     . -f(i1)*Cv(ci(i1)+N0,ep)
     . -Gr(i1,1)*Cv(cii(z2(1,i1))+N0,ep) 
     . -Gr(i1,2)*Cv(cii(z2(2,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCF_00',i1,diff
     
      enddo
      
      return
      end
