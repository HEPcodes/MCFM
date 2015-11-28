      subroutine chkG_iiiii(j,i1,i2,i3,i4,i5,DetGr,Xtwiddle0,
     . Gtwiddle,Shat6,N0)
      implicit none
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,i3,i4,i5,n,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0),Shat6s(np,z5max,-2:0),diff,
     . res
       

      do ep=-2,0
      do n=1,np
      Shat6s(n,z5(i1,i2,i3,i4,i5),ep)=Shat6(n,z5(i1,i2,i3,i4,i5),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzziiii(z4(i2,i3,i4,i5)),ep)
     .        +delta(n,i2)*Dv(N0+dzziiii(z4(i1,i3,i4,i5)),ep)
     .        +delta(n,i3)*Dv(N0+dzziiii(z4(i1,i2,i4,i5)),ep)
     .        +delta(n,i4)*Dv(N0+dzziiii(z4(i1,i2,i3,i5)),ep)
     .        +delta(n,i5)*Dv(N0+dzziiii(z4(i1,i2,i3,i4)),ep))

      enddo 
      res=
     . -(Gtwiddle(j,1)*Shat6s(1,z5(i1,i2,i3,i4,i5),ep)
     .  +Gtwiddle(j,2)*Shat6s(2,z5(i1,i2,i3,i4,i5),ep)
     .  +Gtwiddle(j,3)*Shat6s(3,z5(i1,i2,i3,i4,i5),ep)
     .  -DetGr*Dv(diiiiii(z6(j,i1,i2,i3,i4,i5))+N0,ep))
     . /Xtwiddle0(j)
      diff=
     . +Xtwiddle0(j)*Dv(diiiii(z5(i1,i2,i3,i4,i5))+N0,ep)
     . +Gtwiddle(j,1)*Shat6s(1,z5(i1,i2,i3,i4,i5),ep)
     . +Gtwiddle(j,2)*Shat6s(2,z5(i1,i2,i3,i4,i5),ep)
     . +Gtwiddle(j,3)*Shat6s(3,z5(i1,i2,i3,i4,i5),ep)
     . -DetGr*Dv(diiiiii(z6(j,i1,i2,i3,i4,i5))+N0,ep)
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.))
c     . .and. (abs(res) .gt. 1d-14))
     . write(6,*) 'chkG_iiiii,ep',ep,j,i1,i2,i3,i4,i5,diff,res,
     . diff/res
  
      enddo

      return
      end
