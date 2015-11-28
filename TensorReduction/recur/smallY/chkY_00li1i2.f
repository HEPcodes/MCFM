      subroutine chkY_00li1i2(k,l,i1,i2,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for DD Eq. 5.60c
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,i1,i2,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0),diff,res

      if (  (i1 .eq. l) .or. (i2 .eq. l)
     . .or. (i1 .eq. 0) .or. (i2 .eq. 0)) then
      return
      endif

      do ep=-2,0
      res =          
     . (-2d0*Gtwiddle(k,i1)*Dv(dzziii(z3(l,l,i2))+N0,ep)
     .  -2d0*Gtwiddle(k,i2)*Dv(dzziii(z3(l,l,i1))+N0,ep)
     .  +Gtwiddle(k,1)*Shat5(1,z4(l,l,i1,i2),ep)
     .  +Gtwiddle(k,2)*Shat5(2,z4(l,l,i1,i2),ep)
     .  +Gtwiddle(k,3)*Shat5(3,z4(l,l,i1,i2),ep)
     .  +Xtwiddle(k,0)*Dv(diiii(z4(l,l,i1,i2))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(diiiii(z5(k,l,l,i1,i2))+N0,ep))
     .        /(4d0*Gtwiddle(k,l))
      diff=
     . 4d0*Gtwiddle(k,l)*(Dv(dzziii(z3(l,i1,i2))+N0,ep)-res) 
 
      if (((abs(diff)     .gt. weenumber ) 
     .   .or. (abs(diff/res).gt.weenumber .and.abs(res).gt.weenumber)) 
     .  ) write(6,10) 'chkY_00li1i2',ep,k,l,
     .   i1,i2,diff,res,diff/res

      enddo
 10   format(A10,5I2,6d20.10)

      return
      end
  



