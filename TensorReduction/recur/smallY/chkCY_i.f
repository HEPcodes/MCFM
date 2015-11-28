      subroutine chkCY_i(i,j,i1,f,Xtwiddle,Gtt,Gtwiddle,Shat2,Bzero1,N0)
      implicit none
C---  Expression for Eq. 5.57
C---  Checks Ci, requires C00i
C---  Small terms of order Xtwiddle(0,j)*Cii
C---  Denominator Xtwiddle(i,j)

      include 'constants.f' 
      include 'Cnames.f' 
      include 'Cv.f' 
      include 'Carraydef.f' 
      include 'Carrays.f' 
      include 'weenumber.f' 
      integer ep,N0,i,j,i1,n,m,np
      parameter(np=2)
      double precision Xtwiddle(0:np,0:np),
     . Gtwiddle(np,np),f(np),Gtt(np,np,np,np)
      double complex Shat2(np,z1max,-2:0),Bzero1(z1max,-2:0),
     . bit,pole,diff


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat2(m,i1,ep)-2d0*delta(m,i1)*Cv(cc00+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4*Cv(czzi(i1)+N0,ep-1)

      diff=
     .Cv(ci(i1)+N0,ep)*Xtwiddle(i,j)- 
     .  (Gtwiddle(i,j)*(6d0*Cv(czzi(i1)+N0,ep)+pole-Bzero1(i1,ep))
     . +bit+Xtwiddle(0,j)*Cv(cii(z2(i,i1))+N0,ep))

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCY_i',i,j,i1,diff
     
      enddo
      return
      end
  



