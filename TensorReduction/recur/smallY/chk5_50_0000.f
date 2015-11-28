      subroutine chk5_50_0000(k,l,Xtwiddle,Gtwiddle,Shat4zz,N0)
      implicit none
C---  Expression for D0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adding additional "00" pairs
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(3,3)
      double complex Shat4zz(np,z1max,-2:0),diff,res

      do ep=-2,0
      res =  
     . -(-(Gtwiddle(k,1)*Shat4zz(1,l,ep)
     .    +Gtwiddle(k,2)*Shat4zz(2,l,ep)
     .    +Gtwiddle(k,3)*Shat4zz(3,l,ep)
     .    +Xtwiddle(k,0)*Dv(dzzi(l)+N0,ep)
     .    -Xtwiddle(0,0)*Dv(dzzii(z2(k,l))+N0,ep)
     .        ))/(2d0*Gtwiddle(k,l))

      diff=
     . 2d0*Gtwiddle(k,l)*(Dv(dd0000+N0,ep)-res)
 
      if (((abs(diff)     .gt. weenumber ) 
     .     .or. (abs(diff/res).gt.weenumber .and.abs(res).gt.weenumber)) 
     .     .and. (Gsing .eqv. .false.))  write(6,10) 'chk5_50_0000',
     .     ep,k,l,diff,res,diff/res

      enddo
 10   format(A10,3I2,6d20.10)

      return
      end
  
