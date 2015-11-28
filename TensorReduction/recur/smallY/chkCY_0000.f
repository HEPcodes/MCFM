      subroutine chkCY_0000(k,l,Xtwiddle,Gtwiddle,Shat4zz,N0)
      implicit none
C---  Expression for C0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adcing adcitional "00" pairs
C---  Checks C0000
C---  Small terms of order Xtwiddle(0,k)*C00i,Xtwiddle(0,0)*C00ii
C---  Denominator Gtwiddle(k,l)
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      double complex Shat4zz(np,z1max,-2:0),diff

      do ep=-2,0
      diff=
     .   Cv(cc0000+N0,ep)*(2d0*Gtwiddle(k,l))-(
     .    -(-(Gtwiddle(k,1)*Shat4zz(1,l,ep)
     .    +Gtwiddle(k,2)*Shat4zz(2,l,ep)
     .    +Xtwiddle(k,0)*Cv(czzi(l)+N0,ep)
     .    -Xtwiddle(0,0)*Cv(czzii(z2(k,l))+N0,ep)
     .        )))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0000',k,l,diff
     
      enddo

      return
      end
  
