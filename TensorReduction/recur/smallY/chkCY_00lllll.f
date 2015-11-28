      subroutine chkCY_00lllll(k,l,Xtwiddle,Gtwiddle,Shat7,N0)
      implicit none
C---  Expression for extension of Eq. 5.60a
C---  Checks C00llll
C---  Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7(np,z6max,-2:0),diff

      do ep=-2,0
      diff=
     .Cv(czziiiii(z5(l,l,l,l,l))+N0,ep)*(12d0*Gtwiddle(k,l))-
     . (Gtwiddle(k,1)*Shat7(1,z6(l,l,l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat7(2,z6(l,l,l,l,l,l),ep)
     . +Xtwiddle(k,0)*Cv(ciiiiii(z6(l,l,l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiiiii(z7(k,l,l,l,l,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_00lllll',k,l,diff
     
      enddo

      return
      end
  



