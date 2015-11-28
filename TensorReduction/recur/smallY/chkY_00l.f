      subroutine chkY_00l(k,l,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for DD Eq. 5.56a
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),diff

      do ep=-2,0

      diff=
     . 4*Gtwiddle(k,l)*Dv(dzzi(l)+N0,ep)
     . -(Gtwiddle(k,1)*Shat3(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat3(2,z2(l,l),ep)
     .  +Gtwiddle(k,3)*Shat3(3,z2(l,l),ep)
     .  +Xtwiddle(k,0)*Dv(dii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(diii(z3(k,l,l))+N0,ep))
 
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chk5_56a',k,l,diff

      enddo

      return
      end
  



