      subroutine chkG_iii(j,i1,i2,i3,DetGr,Xtwiddle0,Gtwiddle,
     . Shat4,N0)
C     DD 5.48
      implicit none
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,i3,n,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat4(np,z3max,-2:0),Shat4s(np,z3max,-2:0),diff
       
      do ep=-2,0
      do n=1,np
      Shat4s(n,z3(i1,i2,i3),ep)=Shat4(n,z3(i1,i2,i3),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzzii(z2(i2,i3)),ep)
     .        +delta(n,i2)*Dv(N0+dzzii(z2(i1,i3)),ep)
     .        +delta(n,i3)*Dv(N0+dzzii(z2(i1,i2)),ep))
      enddo 
      diff=
     . +Xtwiddle0(j)*Dv(diii(z3(i1,i2,i3))+N0,ep)
     . +Gtwiddle(j,1)*Shat4s(1,z3(i1,i2,i3),ep)
     . +Gtwiddle(j,2)*Shat4s(2,z3(i1,i2,i3),ep)
     . +Gtwiddle(j,3)*Shat4s(3,z3(i1,i2,i3),ep)
     . -DetGr*Dv(diiii(z4(j,i1,i2,i3))+N0,ep)
      if ((abs(diff) .gt. weenumber)  .and. (Gsing .eqv. .false.))
     . write(6,*) 'chkG_iii',j,i1,i2,i3,diff
      enddo

      return
      end
