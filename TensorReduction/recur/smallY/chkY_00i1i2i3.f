      subroutine chkY_00i1i2i3(k,l,i1,i2,i3,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for DD Eq. 5.60d
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,i2,i3,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0),diff

      if (  (i1 .eq. l) .or. (i2 .eq. l).or. (i3 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0).or. (i3 .eq. 0)) then
      return
      endif

      do ep=-2,0

      diff=
     .  2d0*Gtwiddle(k,l)*Dv(dzziii(z3(i1,i2,i3))+N0,ep)
     . -(-2d0*Gtwiddle(k,i1)*Dv(dzziii(z3(l,i2,i3))+N0,ep)
     .   -2d0*Gtwiddle(k,i2)*Dv(dzziii(z3(l,i1,i3))+N0,ep)
     .   -2d0*Gtwiddle(k,i3)*Dv(dzziii(z3(l,i1,i2))+N0,ep)
     . +Gtwiddle(k,1)*Shat5(1,z4(l,i1,i2,i3),ep)
     . +Gtwiddle(k,2)*Shat5(2,z4(l,i1,i2,i3),ep)
     . +Gtwiddle(k,3)*Shat5(3,z4(l,i1,i2,i3),ep)
     . +Xtwiddle(k,0)*Dv(diiii(z4(l,i1,i2,i3))+N0,ep)
     . -Xtwiddle(0,0)*Dv(diiiii(z5(k,l,i1,i2,i3))+N0,ep))
 
      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_00i1i2i3',k,l,i1,i2,i3,diff

      enddo

      return
      end
  



