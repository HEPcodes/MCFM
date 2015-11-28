      subroutine chkY_000000l(k,l,Xtwiddle,Gtwiddle,Shat7zzzz,N0)
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
      double complex Shat7zzzz(np,z2max,-2:0),diff

      do ep=-2,0
      diff=4d0*Gtwiddle(k,l)*Dv(dzzzzzzi(l)+N0,ep)-
     . (Gtwiddle(k,1)*Shat7zzzz(1,z2(l,l),ep)
     . +Gtwiddle(k,2)*Shat7zzzz(2,z2(l,l),ep)
     . +Gtwiddle(k,3)*Shat7zzzz(3,z2(l,l),ep)
     . +Xtwiddle(k,0)*Dv(dzzzzii(z2(l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(dzzzziii(z3(k,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_000000l',k,l,diff
     
      enddo

      return
      end
  



