      subroutine chkCY_00i1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat4,N0)
      implicit none
C---  Expression for Eq. 5.58c
C---  Checks C00i1i2, requires C00li1,C00li2
C---  Small terms of order Xtwiddle(0,k)*Ciii,Xtwiddle(0,0)*Ciiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0),diff

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0

      diff=
     .Cv(czzii(z2(i1,i2))+N0,ep)*(2*Gtwiddle(k,l))-
     .(-2*Gtwiddle(k,i1)*Cv(czzii(z2(l,i2))+N0,ep)
     . -2*Gtwiddle(k,i2)*Cv(czzii(z2(l,i1))+N0,ep)
     . +Gtwiddle(k,1)*Shat4(1,z3(l,i1,i2),ep)
     . +Gtwiddle(k,2)*Shat4(2,z3(l,i1,i2),ep)
     . +Xtwiddle(0,k)*Cv(ciii(z3(l,i1,i2))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiii(z4(k,l,i1,i2))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_00i1i2',k,l,i1,i2,diff
     
      enddo

      return
      end
  



