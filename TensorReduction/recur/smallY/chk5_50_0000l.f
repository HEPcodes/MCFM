      subroutine chk5_50_0000l(k,l,Xtwiddle,Gtwiddle,Shat5zz,N0)
      implicit none
C---  Expression for D0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adding additional "00" pairs
C---  (similar to Eq. 5.56a but with "00" added) 
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat5zz(np,z2max,-2:0),diff,res

      do ep=-2,0

      res = -(
     . -(Gtwiddle(k,1)*Shat5zz(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat5zz(2,z2(l,l),ep)
     .  +Gtwiddle(k,3)*Shat5zz(3,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Dv(dzzii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Dv(dzziii(z3(k,l,l))+N0,ep)))/(4d0*Gtwiddle(k,l))

      diff=
     . 4*Gtwiddle(k,l)*(Dv(dzzzzi(l)+N0,ep)-res)
 
            if (((abs(diff)     .gt. weenumber ) 
     .     .or. (abs(diff/res).gt.weenumber .and.abs(res).gt.weenumber)) 
     .     .and. (Gsing .eqv. .false.))  write(6,10) 'chk5_50_0000i',
     .     ep,k,l,diff,res,diff/res

      enddo
 10   format(A10,3I2,6d20.10)

      return
      end
  

