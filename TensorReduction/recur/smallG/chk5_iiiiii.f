      subroutine chk5_iiiiii(j,i1,i2,i3,i4,i5,i6,DetGr,Xtwiddle0,
     . Gtwiddle,Shat7,N0)
      implicit none
      include 'Dnames.f'  
      include 'Dv.f'  
      include 'Darraydef.f'  
      include 'Darrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,i3,i4,i5,n,np
      parameter(np=3)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat7(np,z6max,-2:0),Shat6s(np,z6max,-2:0),diff,
     . res
       

      do ep=-2,0
      do n=1,np
      Shat7s(n,z5(i1,i2,i3,i4,i5,i6),ep)=
     . Shat7(n,z5(i1,i2,i3,i4,i5,i6),ep)
     .   -2d0*(delta(n,i1)*Dv(N0+dzziiiii(z5(i2,i3,i4,i5,i6)),ep)
     .        +delta(n,i2)*Dv(N0+dzziiiii(z5(i1,i3,i4,i5,i6)),ep)
     .        +delta(n,i3)*Dv(N0+dzziiiii(z5(i1,i2,i4,i5,i6)),ep)
     .        +delta(n,i4)*Dv(N0+dzziiiii(z5(i1,i2,i3,i5,i6)),ep)
     .        +delta(n,i5)*Dv(N0+dzziiiii(z5(i1,i2,i3,i4,i6)),ep)
     .        +delta(n,i6)*Dv(N0+dzziiiii(z5(i1,i2,i3,i4,i5)),ep))

      enddo 
      res=
     . -(Gtwiddle(j,1)*Shat6s(1,z6(i1,i2,i3,i4,i5,i6),ep)
     .  +Gtwiddle(j,2)*Shat6s(2,z6(i1,i2,i3,i4,i5,i6),ep)
     .  +Gtwiddle(j,3)*Shat6s(3,z6(i1,i2,i3,i4,i5,i6),ep)
c     .  -0*DetGr*Dv(diiiiiii(z7(j,i1,i2,i3,i4,i5))+N0,ep)
      )/Xtwiddle0(j)
      diff=
     . +Xtwiddle0(j)*Dv(diiiiii(z6(i1,i2,i3,i4,i5,i6))+N0,ep)
     . -Gtwiddle(j,1)*Shat7s(1,z5(i1,i2,i3,i4,i5,i6),ep)
     . +Gtwiddle(j,2)*Shat7s(2,z5(i1,i2,i3,i4,i5,i6),ep)
     . +Gtwiddle(j,3)*Shat7s(3,z5(i1,i2,i3,i4,i5,i6),ep)
c     . -DetGr*Dv(diiiiii(z7(j,i1,i2,i3,i4,i5,i6))+N0,ep)
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.))
c     . .and. (abs(res) .gt. 1d-14))
     . write(6,*) 'chk5_iiiii,ep',ep,j,i1,i2,i3,i4,i5,i6,diff,res,
     . diff/res
  
      enddo

      return
      end
