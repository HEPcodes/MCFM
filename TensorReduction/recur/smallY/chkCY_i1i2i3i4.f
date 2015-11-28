      subroutine chkCY_i1i2i3i4(i,j,i1,i2,i3,i4,f,Xtwiddle,Gtt,Gtwiddle,
     . Shat5,Bzero4,N0)
      implicit none
C---  Expression for Eq. 5.61
      include 'constants.f' 
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,i,j,i1,i2,i3,i4,n,m,np
      parameter(np=2)
      double precision 
     . Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat5(np,z4max,-2:0),Bzero4(z4max,-2:0),
     . bit,pole,diff


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat5(m,z4(i1,i2,i3,i4),ep)
     . -2d0*delta(m,i1)*Cv(czziii(z3(i2,i3,i4))+N0,ep)
     . -2d0*delta(m,i2)*Cv(czziii(z3(i1,i3,i4))+N0,ep)
     . -2d0*delta(m,i3)*Cv(czziii(z3(i1,i2,i4))+N0,ep)
     . -2d0*delta(m,i4)*Cv(czziii(z3(i1,i2,i3))+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep-1)

      diff=
     .Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)*Xtwiddle(i,j)-
     .  (Gtwiddle(i,j)
     . *(12d0*Cv(czziiii(z4(i1,i2,i3,i4))+N0,ep)
     .  +pole-Bzero4(z4(i1,i2,i3,i4),ep))
     .  +bit+Xtwiddle(0,j)*Cv(ciiiii(z5(i,i1,i2,i3,i4))+N0,ep)
     .   )

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_i1i2i3i4',i,j,i1,i2,i3,i4,diff
     
      enddo
      return
      end
  



