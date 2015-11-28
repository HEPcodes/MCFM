      subroutine chkCY_0000ll(k,l,Xtwiddle,Gtwiddle,Shat6zz,N0)
      implicit none
C---  Expression for extension of Eq. 5.60a
C---  Checks C00llll
C---  Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6zz(np,z3max,-2:0),diff

      do ep=-2,0
      diff=
     .Cv(czzzzii(z2(l,l))+N0,ep)*(6d0*Gtwiddle(k,l))-
     . (Gtwiddle(k,1)*Shat6zz(1,z3(l,l,l),ep)
     . +Gtwiddle(k,2)*Shat6zz(2,z3(l,l,l),ep)
     . +Xtwiddle(k,0)*Cv(czziii(z3(l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiii(z4(k,l,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0000ll',k,l,diff
     
      enddo

      return
      end
  



