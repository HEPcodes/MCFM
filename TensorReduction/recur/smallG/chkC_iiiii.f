      subroutine chkC_iiiii(j,i1,i2,i3,i4,i5,DetGr,Xtwiddle0,
     . Gtwiddle,Shat6,N0)
      implicit none
      include 'Cnames.f'  
      include 'Cv.f'  
      include 'Carraydef.f'  
      include 'Carrays.f'  
      include 'weenumber.f'  
      integer ep,N0,j,i1,i2,i3,i4,i5,n,np
      parameter(np=2)
      double precision DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      double complex Shat6(np,z5max,-2:0),Shat6s(np,z5max,-2:0),diff
       

      do ep=-2,0
      do n=1,np
      Shat6s(n,z5(i1,i2,i3,i4,i5),ep)=Shat6(n,z5(i1,i2,i3,i4,i5),ep)
     .   -2d0*(delta(n,i1)*Cv(N0+czziiii(z4(i2,i3,i4,i5)),ep)
     .        +delta(n,i2)*Cv(N0+czziiii(z4(i1,i3,i4,i5)),ep)
     .        +delta(n,i3)*Cv(N0+czziiii(z4(i1,i2,i4,i5)),ep)
     .        +delta(n,i4)*Cv(N0+czziiii(z4(i1,i2,i3,i5)),ep)
     .        +delta(n,i5)*Cv(N0+czziiii(z4(i1,i2,i3,i4)),ep))
      enddo 

      diff=
     . Xtwiddle0(j)*Cv(ciiiii(z5(i1,i2,i3,i4,i5))+N0,ep)
     . +Gtwiddle(j,1)*Shat6s(1,z5(i1,i2,i3,i4,i5),ep)
     . +Gtwiddle(j,2)*Shat6s(2,z5(i1,i2,i3,i4,i5),ep)
     . -DetGr*Cv(ciiiiii(z6(j,i1,i2,i3,i4,i5))+N0,ep)
 
      if ((abs(diff) .gt. weenumber) .and. (Gsing .eqv. .false.))
     . write(6,*) 'chkC_iiiii,ep',ep,j,i1,i2,i3,i4,i5,diff
  
      enddo

      return
      end
