      subroutine chkF_00iii(i1,i2,i3,f,Gr,Shat5,N0)
C---  Expression for rearrangement of extension of Eq. 5.70
C---  Checks D00iii
C---  Small terms of order f(i)*Dijkl,Gr(i,j)*Dijklm
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,np
      parameter(np=3)
      double precision f(np),Gr(np,np),den
      double complex Shat5(np,z4max,-2:0),diff
       
      do ep=-2,0
      if     ((i1 .eq. i2) .and. (i1 .eq. i3)) then
        den=8d0
	k=i1
      elseif (i1 .eq. i2) then
        den=6d0
	k=i1
      elseif (i1 .eq. i3) then
        den=6d0
	k=i1
      elseif (i2 .eq. i3) then
        den=6d0
	k=i2
      else
        den=4d0
	k=i1
      endif      
      diff=Dv(dzziii(z3(i1,i2,i3))+N0,ep)*den-
     . (Shat5(k,z4(k,i1,i2,i3),ep)
     . -f(k)*Dv(diiii(z4(k,i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Dv(diiiii(z5(1,k,i1,i2,i3))+N0,ep) 
     . -Gr(k,2)*Dv(diiiii(z5(2,k,i1,i2,i3))+N0,ep)
     . -Gr(k,3)*Dv(diiiii(z5(3,k,i1,i2,i3))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkF_00iii',i1,i2,i3,diff
     
      enddo
      
      return
      end
