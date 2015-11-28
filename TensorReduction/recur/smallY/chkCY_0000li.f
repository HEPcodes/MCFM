      subroutine chkCY_0000li(k,l,i1,Xtwiddle,Gtwiddle,Shat6zz,N0)
      implicit none
C---  Expression for Eq. 5.58b
C---  Checks C0000li, requires C0000ll
C---  Small terms of order Xtwiddle(0,k)*C00iii,Xtwiddle(0,0)*C00iiii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6zz(np,z3max,-2:0),diff

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      diff=
     .Cv(czzzzii(z2(l,i1))+N0,ep)*(4d0*Gtwiddle(k,l))-
     . (-2d0*Gtwiddle(k,i1)*Cv(czzzzii(z2(l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat6zz(1,z3(l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat6zz(2,z3(l,l,i1),ep)
     . +Xtwiddle(0,k)*Cv(czziii(z3(l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiii(z4(k,l,l,i1))+N0,ep)
     . )

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0000li',k,l,i1,diff
      
      enddo

      return
      end
  



