      subroutine chkY_0(i,j,f,Xtwiddle,Gtwiddle,Gtt,Shat1,Czero0,N0)
      implicit none
C---  Expression for Eq. 5.55
      include 'constants.f' 
      include 'Dnames.f' 
      include 'Dv.f' 
      include 'Darraydef.f' 
      include 'Darrays.f' 
      include 'weenumber.f' 
      integer ep,N0,i,j,n,m,np
      parameter(np=3)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),
     . Gtt(np,np,np,np)
      double complex Shat1(np,-2:0),Czero0(-2:0),
     . bit,diff,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)*Shat1(m,ep)
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4d0*Dv(dd00+N0,ep-1)
      diff=
     .  Xtwiddle(i,j)*Dv(dd0+N0,ep)
     .  -(Gtwiddle(i,j)*(2d0*Dv(dd00+N0,ep)+pole-Czero0(ep))
     . +bit+Xtwiddle(0,j)*Dv(di(i)+N0,ep))

      if ((abs(diff) .gt. weenumber))
     . write(6,*) 'chkY_0',i,j,diff,Xtwiddle(i,j)*Dv(dd0+N0,ep)
     
      enddo
      
      return
      end
  



