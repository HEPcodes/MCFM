      subroutine chkF_0(m0sq,Gr,Czero0,N0)
C---  Expression for rearrangement of Eq. 5.65
C---  Checks D0, requires D00
C---  Small terms of order Gr(i,j)*Dij
C---  Denominator m0sq
      implicit none
      include 'constants.f'
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex Czero0(-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(dii(z2(n,m))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dd00+N0,ep-1)
      
      diff=Dv(dd0+N0,ep)*2d0*m0sq-
     . (8d0*Dv(dd00+N0,ep)
     . -pole
     . -2d0*Czero0(ep)
     . +bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkF_0',diff
     
      enddo
      
      return
      end
