      subroutine chkY_00li1i2i3(k,l,i1,i2,i3,Xtwiddle,Gtwiddle,Shat6,N0)
      implicit none
C---  Expression for extension of Eq. 5.60c
C---  Checks D00li1i2i3, requires D00lli1i2
C---  Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---  Denominator Gtwiddle(k,l)
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,i2,i3,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0),diff

      if (  (i1 .eq. l) .or. (i2 .eq. l) .or. (i3 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0) .or. (i3 .eq. 0)) then
      return
      endif

      do ep=-2,0
      diff=4d0*Gtwiddle(k,l)*Dv(dzziiii(z4(l,i1,i2,i3))+N0,ep)-
     .(-2d0*Gtwiddle(k,i1)*Dv(dzziiii(z4(l,l,i2,i3))+N0,ep)
     . -2d0*Gtwiddle(k,i2)*Dv(dzziiii(z4(l,l,i1,i3))+N0,ep)
     . -2d0*Gtwiddle(k,i3)*Dv(dzziiii(z4(l,l,i1,i2))+N0,ep)
     . +Gtwiddle(k,1)*Shat6(1,z5(l,l,i1,i2,i3),ep)
     . +Gtwiddle(k,2)*Shat6(2,z5(l,l,i1,i2,i3),ep)
     . +Gtwiddle(k,3)*Shat6(3,z5(l,l,i1,i2,i3),ep)
     . +Xtwiddle(k,0)*Dv(diiiii(z5(l,l,i1,i2,i3))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiiii(z6(k,l,l,i1,i2,i3))+N0,ep))
 
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_00li1i2i3',k,l,i1,i2,i3,diff
     
      enddo


      return
      end
  



