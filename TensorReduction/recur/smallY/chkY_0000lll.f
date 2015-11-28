      subroutine chkY_0000lll(k,l,Xtwiddle,Gtwiddle,Shat7zz,N0)
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
      double complex Shat7zz(np,z4max,-2:0),diff

      do ep=-2,0
      diff=8d0*Gtwiddle(k,l)*Dv(dzzzziii(z3(l,l,l))+N0,ep)-
     . (Gtwiddle(k,1)*Shat7zz(1,z4(l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,l,l,l),ep)
     . +Gtwiddle(k,3)*Shat7zz(3,z4(l,l,l,l),ep)
     . +Xtwiddle(k,0)*Dv(dzziiii(z4(l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzziiiii(z5(k,l,l,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_0000lll',k,l,diff

      enddo

      return
      end
  



