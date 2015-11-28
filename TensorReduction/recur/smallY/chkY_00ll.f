      subroutine chkY_00ll(k,l,Xtwiddle,Gtwiddle,Shat4,N0)
      implicit none
C---  Expression for Eq. 5.58a
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0),diff

      do ep=-2,0

      diff=
     . 6*Gtwiddle(k,l)*Dv(dzzii(z2(l,l))+N0,ep)
     . -(Gtwiddle(k,1)*Shat4(1,z3(l,l,l),ep)
     .  +Gtwiddle(k,2)*Shat4(2,z3(l,l,l),ep)
     .  +Gtwiddle(k,3)*Shat4(3,z3(l,l,l),ep)
     .  +Xtwiddle(k,0)*Dv(diii(z3(l,l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(diiii(z4(k,l,l,l))+N0,ep))
 
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chk5_58a',k,l,diff
      enddo


      return
      end
  



