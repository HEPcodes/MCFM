      subroutine chkF_00(i1,f,Gr,Shat2,N0)
C---  Expression for rearrangement of Eq. 5.66
C---  Checks D00
C---  Small terms of order f(i)*Di,Gr(i,j)*Dij
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,np
      parameter(np=3)
      double precision f(np),Gr(np,np)
      double complex Shat2(np,np,-2:0),diff
       
      do ep=-2,0
      diff=Dv(dd00+N0,ep)*2d0-
     . (Shat2(i1,i1,ep)
     . -f(i1)*Dv(di(i1)+N0,ep)
     . -Gr(i1,1)*Dv(dii(z2(1,i1))+N0,ep) 
     . -Gr(i1,2)*Dv(dii(z2(2,i1))+N0,ep)
     . -Gr(i1,3)*Dv(dii(z2(3,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkF_00',i1,diff
     
      enddo
      
      return
      end
