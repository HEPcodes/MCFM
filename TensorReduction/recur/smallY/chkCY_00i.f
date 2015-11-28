      subroutine chkCY_00i(k,l,i1,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for Eq. 5.56b
C---  Checks C00i1
C---  Requires C00l
C---  Small terms of order Xtwiddle(0,k)*Dli1,Xtwiddle(0,0)*Dkli1
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),diff

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      diff=
     .Cv(czzi(i1)+N0,ep)*(2*Gtwiddle(k,l))-
     . (-2*Gtwiddle(k,i1)*Cv(czzi(l)+N0,ep)
     . +Gtwiddle(k,1)*Shat3(1,z2(l,i1),ep)
     . +Gtwiddle(k,2)*Shat3(2,z2(l,i1),ep)
     . +Xtwiddle(0,k)*Cv(cii(z2(l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciii(z3(k,l,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_00i',k,l,i1,diff
     
      enddo

      return
      end
  



