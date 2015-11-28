      subroutine chkG_ii(j,i1,i2,DetGr,Xtwiddle0,Gtwiddle,Shat3,N0)
      implicit none
C     DD 5.45
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,n,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat3(np,z2max,-2:0),Shat3s(np,z2max,-2:0),diff
       
      do ep=-2,0
      do n=1,np
      Shat3s(n,z2(i1,i2),ep)=Shat3(n,z2(i1,i2),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzzi(i2),ep)
     .        +delta(n,i2)*Dv(N0+dzzi(i1),ep))
      enddo 
         
      diff=
     . +Xtwiddle0(j)*Dv(dii(z2(i1,i2))+N0,ep)
     . +Gtwiddle(j,1)*Shat3s(1,z2(i1,i2),ep)
     . +Gtwiddle(j,2)*Shat3s(2,z2(i1,i2),ep)
     . +Gtwiddle(j,3)*Shat3s(3,z2(i1,i2),ep)
     . -DetGr*Dv(diii(z3(j,i1,i2))+N0,ep)
      if ((abs(diff) .gt. weenumber)  .and. (Gsing .eqv. .false.))
     . write(6,*) 'chkG_ii',j,i1,i2,diff
      enddo

      return
      end
