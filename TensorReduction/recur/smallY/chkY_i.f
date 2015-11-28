      subroutine chkY_i(i,j,i1,f,Xtwiddle,Gtt,Gtwiddle,Shat2,Czero1,N0)
      implicit none
C---  Expression for Eq. 5.57
      include 'constants.f' 
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,i,j,i1,n,m,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),
     . Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat2(np,z1max,-2:0),Czero1(z1max,-2:0),
     . bit,diff,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat2(m,i1,ep)-2d0*delta(m,i1)*Dv(dd00+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) 
     . pole=-4d0*Dv(dzzi(i1)+N0,ep-1)
      diff=
     .  Xtwiddle(i,j)*Dv(di(i1)+N0,ep)
     .  -(Gtwiddle(i,j)*(4d0*Dv(dzzi(i1)+N0,ep)+pole-Czero1(i1,ep))
     . +bit+Xtwiddle(0,j)*Dv(dii(z2(i,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chk5_57',i,j,i1,diff

      enddo

      return
      end
  



