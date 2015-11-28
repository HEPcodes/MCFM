      subroutine chkCY_0000l(k,l,Xtwiddle,Gtwiddle,Shat5zz,N0)
      implicit none
C---  Expression for C0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adcing adcitional "00" pairs
C---  (similar to Eq. 5.56a but with "00" added) 
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5zz(np,z2max,-2:0),diff

      do ep=-2,0

      diff=
     .Cv(czzzzi(l)+N0,ep)*(4d0*Gtwiddle(k,l))-( -(
     . -(Gtwiddle(k,1)*Shat5zz(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat5zz(2,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Cv(czzii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Cv(czziii(z3(k,l,l))+N0,ep))))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0000l',k,l,diff
     
      enddo

      return
      end
  

