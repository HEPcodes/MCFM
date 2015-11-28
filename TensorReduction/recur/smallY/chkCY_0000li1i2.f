      subroutine chkCY_0000li1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat7zz,N0)
      implicit none
C---  Expression for Eq. 5.60c
C---  Checks C0000li1i2, requires C0000lli1
C---  Small terms of order Xtwiddle(0,k)*C00iiii,Xtwiddle(0,0)*C00iiiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat7zz(np,z4max,-2:0),diff

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0
      diff=
     .Cv(czzzziii(z3(l,i1,i2))+N0,ep)*(4d0*Gtwiddle(k,l))-
     .(-2d0*Gtwiddle(k,i1)*Cv(czzzziii(z3(l,l,i2))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Cv(czzzziii(z3(l,l,i1))+N0,ep)
     . +Gtwiddle(k,1)*Shat7zz(1,z4(l,l,i1,i2),ep)
     . +Gtwiddle(k,2)*Shat7zz(2,z4(l,l,i1,i2),ep)
     . +Xtwiddle(0,k)*Cv(czziiii(z4(l,l,i1,i2))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiiii(z5(k,l,l,i1,i2))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0000li1i2',k,l,i1,i2,diff
      
      enddo

      return
      end
  



