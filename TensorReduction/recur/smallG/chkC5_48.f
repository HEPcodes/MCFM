      subroutine chkC5_48(j,i1,i2,i3,DetGr,Xtwiddle0,Gtwiddle,
     . Shat4,N0)
      implicit none
      include 'constants.f'  
      include 'Cnames.f'  
      include 'Cv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,i3,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0),Shat4s(np,z3max,-2:0),diff
       
      do ep=-2,0
      do n=1,np
      Shat4s(n,z3(i1,i2,i3),ep)=Shat4(n,z3(i1,i2,i3),ep)
     .   -2d0*(delta(n,i1)*Cv(N0+czzii(z2(i2,i3)),ep)
     .        +delta(n,i2)*Cv(N0+czzii(z2(i1,i3)),ep)
     .        +delta(n,i3)*Cv(N0+czzii(z2(i1,i2)),ep))
      enddo 
      diff=
     . +Xtwiddle0(j)*Cv(ciii(z3(i1,i2,i3))+N0,ep)
     . +Gtwiddle(j,1)*Shat4s(1,z3(i1,i2,i3),ep)
     . +Gtwiddle(j,2)*Shat4s(2,z3(i1,i2,i3),ep)
     . -DetGr*Cv(ciiii(z4(j,i1,i2,i3))+N0,ep)
      if (abs(diff) .gt. weenumber) 
     . write(6,*) 'chkC5_48',j,i1,i2,i3,diff
      enddo

      return
      end
