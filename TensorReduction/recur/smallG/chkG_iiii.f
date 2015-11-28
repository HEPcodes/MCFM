      subroutine chkG_iiii(j,i1,i2,i3,i4,DetGr,Xtwiddle0,Gtwiddle,
     . Shat5,N0)
      implicit none
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,i3,i4,n,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat5(np,z4max,-2:0),Shat5s(np,z4max,-2:0),diff
       

      do ep=-2,0
      do n=1,np
      Shat5s(n,z4(i1,i2,i3,i4),ep)=Shat5(n,z4(i1,i2,i3,i4),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzziii(z3(i2,i3,i4)),ep)
     .        +delta(n,i2)*Dv(N0+dzziii(z3(i1,i3,i4)),ep)
     .        +delta(n,i3)*Dv(N0+dzziii(z3(i1,i2,i4)),ep)
     .        +delta(n,i4)*Dv(N0+dzziii(z3(i1,i2,i3)),ep))
      enddo 
          
      diff=
     . +Xtwiddle0(j)*Dv(diiii(z4(i1,i2,i3,i4))+N0,ep)
     . +Gtwiddle(j,1)*Shat5s(1,z4(i1,i2,i3,i4),ep)
     . +Gtwiddle(j,2)*Shat5s(2,z4(i1,i2,i3,i4),ep)
     . +Gtwiddle(j,3)*Shat5s(3,z4(i1,i2,i3,i4),ep)
     . -DetGr*Dv(diiiii(z5(j,i1,i2,i3,i4))+N0,ep)
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.)) 
     . write(6,*) 'chkG_iiii',j,i1,i2,i3,i4,diff
      enddo

      return
      end
