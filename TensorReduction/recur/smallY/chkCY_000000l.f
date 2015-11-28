      subroutine chkCY_000000l(k,l,Xtwiddle,Gtwiddle,Shat7zzzz,N0)
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
      double complex Shat7zzzz(np,z2max,-2:0),diff

      do ep=-2,0
      diff=
     .Cv(czzzzzzi(l)+N0,ep)*(4d0*Gtwiddle(k,l))-
     . (Gtwiddle(k,1)*Shat7zzzz(1,z2(l,l),ep)
     . +Gtwiddle(k,2)*Shat7zzzz(2,z2(l,l),ep)
     . +Xtwiddle(k,0)*Cv(czzzzii(z2(l,l))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czzzziii(z3(k,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_000000l',k,l,diff
     
      enddo

      return
      end
  



