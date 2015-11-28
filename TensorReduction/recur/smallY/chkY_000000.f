      subroutine chkY_000000(k,l,Xtwiddle,Gtwiddle,Shat6zzzz,N0)
      implicit none
C---  Expression for D0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adding additional "00" pairs
C---  Checks D0000
C---  Small terms of order Xtwiddle(0,k)*D00i,Xtwiddle(0,0)*D00ii
C---  Denominator Gtwiddle(k,l)
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(3,3)
      double complex Shat6zzzz(np,z1max,-2:0),diff

      do ep=-2,0
      diff=2d0*Gtwiddle(k,l)*Dv(dd000000+N0,ep)-(  
     . -(-(Gtwiddle(k,1)*Shat6zzzz(1,l,ep)
     .    +Gtwiddle(k,2)*Shat6zzzz(2,l,ep)
     .    +Gtwiddle(k,3)*Shat6zzzz(3,l,ep)
     .    +Xtwiddle(k,0)*Dv(dzzzzi(l)+N0,ep)
     .    -Xtwiddle(0,0)*Dv(dzzzzii(z2(k,l))+N0,ep)
     .        ))
     .  )

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_000000',k,l,diff
     
      enddo

      return
      end
  
