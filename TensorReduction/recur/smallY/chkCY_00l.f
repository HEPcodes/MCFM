      subroutine chkCY_00l(k,l,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for Eq. 5.56a
C---  Checks C00l
C---  Small terms of order Xtwiddle(0,k)*Cii,Xtwiddle(0,0)*Ciii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),diff

      do ep=-2,0
      diff=
     .Cv(czzi(l)+N0,ep)*(4*Gtwiddle(k,l))-
     . (Gtwiddle(k,1)*Shat3(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat3(2,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Cv(cii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Cv(ciii(z3(k,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_00l',k,l,diff
     
      enddo

      return
      end
  



