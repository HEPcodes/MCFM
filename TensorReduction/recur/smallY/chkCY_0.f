      subroutine chkCY_0(i,j,f,Xtwiddle,Gtwiddle,Gtt,Shat1,Bzero0,N0)
      implicit none
C---  Expression for Eq. 5.55
C---  Checks C0, requires C00
C---  Small terms of order Xtwiddle(0,j)*Ci
C---  Denominator Xtwiddle(i,j)
      include 'constants.f' 
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,i,j,n,m,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),
     . Gtt(np,np,np,np)
      double complex Shat1(np,-2:0),Bzero0(-2:0),
     . bit,pole,diff

      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)*Shat1(m,ep)
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4*Cv(cc00+N0,ep-1)

      diff=
     .Cv(cc0+N0,ep)*Xtwiddle(i,j)-
     .  (Gtwiddle(i,j)*(4d0*Cv(cc00+N0,ep)+pole-Bzero0(ep))
     . +bit+Xtwiddle(0,j)*Cv(ci(i)+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_0',i,j,diff
     
      enddo

      return
      end
  



