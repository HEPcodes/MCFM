      subroutine chkY_00lllll(k,l,Xtwiddle,Gtwiddle,Shat7,N0)
      implicit none
C---  Expression for extension of Eq. 5.60a
C---  Checks D00llll
C---  Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---  Denominator Gtwiddle(k,l)
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7(np,z6max,-2:0),diff

      do ep=-2,0
      diff=12d0*Gtwiddle(k,l)*Dv(dzziiiii(z5(l,l,l,l,l))+N0,ep)-
     . (Gtwiddle(k,1)*Shat7(1,z6(l,l,l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat7(2,z6(l,l,l,l,l,l),ep)
     . +Gtwiddle(k,3)*Shat7(3,z6(l,l,l,l,l,l),ep)
     . +Xtwiddle(k,0)*Dv(diiiiii(z6(l,l,l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiiiii(z7(k,l,l,l,l,l,l))+N0,ep))

      if ((abs(diff) .gt. weenumber))  then
      write(6,*) 'chkY_00lllll',k,l,diff
      endif
      
      enddo

      return
      end
  



