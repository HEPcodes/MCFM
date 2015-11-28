      subroutine chkY_00(k,l,Xtwiddle,Gtwiddle,Shat2,N0)
      implicit none
C---  Expression for DD Eq. 5.54
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat2(np,z1max,-2:0),diff

      do ep=-2,0

      diff=
     . 2*Gtwiddle(k,l)*Dv(dd00+N0,ep)
     . -(Gtwiddle(k,1)*Shat2(1,l,ep)
     .  +Gtwiddle(k,2)*Shat2(2,l,ep)
     .  +Gtwiddle(k,3)*Shat2(3,l,ep)
     .  +Xtwiddle(k,0)*Dv(di(l)+N0,ep)
     .  -Xtwiddle(0,0)*Dv(dii(z2(k,l))+N0,ep))
 
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chk5_54',k,l,diff

      enddo

      return
      end
  



