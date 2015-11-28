      subroutine chkCY_00lll(k,l,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for Eq. 5.60a
C---  Checks C00lll
C---  Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0),diff

      do ep=-2,0

      diff=
     .Cv(czziii(z3(l,l,l))+N0,ep)*(8*Gtwiddle(k,l))-
     . (Gtwiddle(k,1)*Shat5(1,z4(l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat5(2,z4(l,l,l,l),ep)
     . +Xtwiddle(0,k)*Cv(ciiii(z4(l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiii(z5(k,l,l,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_00lll',k,l,diff
     
      enddo

      return
      end
  



