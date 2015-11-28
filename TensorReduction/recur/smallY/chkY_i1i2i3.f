      subroutine chkY_i1i2i3(i,j,i1,i2,i3,f,Xtwiddle,Gtt,Gtwiddle,Shat4,
     . Czero3,N0)
      implicit none
C---  Expression for DD Eq. 5.61
      include 'constants.f' 
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,i,j,i1,i2,i3,n,m,np
      parameter(np=3)
      double precision 
     . Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat4(np,z3max,-2:0),Czero3(z3max,-2:0),
     . bit,diff,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat4(m,z3(i1,i2,i3),ep)
     . -2d0*delta(m,i1)*Dv(dzzii(z2(i2,i3))+N0,ep)
     . -2d0*delta(m,i2)*Dv(dzzii(z2(i1,i3))+N0,ep)
     . -2d0*delta(m,i3)*Dv(dzzii(z2(i1,i2))+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Dv(dzziii(z3(i1,i2,i3))+N0,ep-1)
      diff=
     .  Xtwiddle(i,j)*Dv(diii(z3(i1,i2,i3))+N0,ep)
     .  -(Gtwiddle(i,j)
     . *(8d0*Dv(dzziii(z3(i1,i2,i3))+N0,ep)
     . +pole-Czero3(z3(i1,i2,i3),ep))
     . +bit+Xtwiddle(0,j)*Dv(diiii(z4(i,i1,i2,i3))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkY_i1i2i3',i,j,i1,i2,i3,diff

      enddo
      
      return
      end
  



