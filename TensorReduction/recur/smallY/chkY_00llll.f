      subroutine chkY_00llll(k,l,Xtwiddle,Gtwiddle,Shat6,N0)
      implicit none
C---  Expression for extension of Eq. 5.60a
C---  Checks D00llll
C---  Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---  Denominator Gtwiddle(k,l)
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0),diff

      do ep=-2,0
      diff=
     . 10d0*Gtwiddle(k,l)*Dv(dzziiii(z4(l,l,l,l))+N0,ep)-
     . (Gtwiddle(k,1)*Shat6(1,z5(l,l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat6(2,z5(l,l,l,l,l),ep)
     . +Gtwiddle(k,3)*Shat6(3,z5(l,l,l,l,l),ep)
     . +Xtwiddle(k,0)*Dv(diiiii(z5(l,l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiiii(z6(k,l,l,l,l,l))+N0,ep))
     
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_00llll',k,l,diff,
     . Xtwiddle(k,0)*Dv(diiiii(z5(l,l,l,l,l))+N0,ep),
     . Xtwiddle(0,0),Dv(diiiiii(z6(k,l,l,l,l,l))+N0,ep),
     . 10d0*Gtwiddle(k,l)*Dv(dzziiii(z4(l,l,l,l))+N0,ep)
     
      enddo

      return
      end
  



