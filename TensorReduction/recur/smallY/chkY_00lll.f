      subroutine chkY_00lll(k,l,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for DD Eq. 5.60a
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0),diff,res

      do ep=-2,0

      res =
     .  (Gtwiddle(k,1)*Shat5(1,z4(l,l,l,l),ep)
     .  +Gtwiddle(k,2)*Shat5(2,z4(l,l,l,l),ep)
     .  +Gtwiddle(k,3)*Shat5(3,z4(l,l,l,l),ep)
     .  +Xtwiddle(k,0)*Dv(diiii(z4(l,l,l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(diiiii(z5(k,l,l,l,l))+N0,ep))
     . /(8d0*Gtwiddle(k,l))
      diff=
     . 8d0*Gtwiddle(k,l)*(Dv(dzziii(z3(l,l,l))+N0,ep)-res)
 
      if (((abs(diff)     .gt. weenumber ) 
     .    .or. (abs(diff/res).gt.weenumber .and.abs(res).gt.weenumber)) 
     .   )  write(6,10) 'chkY_00lll',ep,k,l,
     .    diff,res,diff/res

      enddo
      
 10   format(A10,3I2,7d20.10)

      return
      end
  



