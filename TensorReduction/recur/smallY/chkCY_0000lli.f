      subroutine chkCY_0000lli(k,l,i1,Xtwiddle,Gtwiddle,Shat7zz,N0)
      implicit none
C---  Expression for Eq. 5.60b
C---  Checks C0000lli, requires C0000lll
C---  Small terms of order Xtwiddle(0,k)*C00iiii,Xtwiddle(0,0)*C00iiiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zz(np,z4max,-2:0),diff

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      diff=
     .Cv(czzzziii(z3(l,l,i1))+N0,ep)*(6d0*Gtwiddle(k,l))-
     . (-2d0*Gtwiddle(k,i1)*Cv(czzzziii(z3(l,l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat7zz(1,z4(l,l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,l,l,i1),ep)
     . +Xtwiddle(k,0)*Cv(czziiii(z4(l,l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiiii(z5(k,l,l,l,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0000lll',k,l,diff
      
      enddo

      return
      end
  



