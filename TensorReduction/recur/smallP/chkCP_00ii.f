      subroutine chkCP_00ii(i1,i2,m0sq,Gr,Bzero2,N0)
      implicit none
      include 'constants.f'
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,i2,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero2(z2max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiii(z4(n,m,i1,i2))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czzii(z2(i1,i2))+N0,ep-1)
      
      diff=Cv(czzii(z2(i1,i2))+N0,ep)*16d0-
     . (pole
     . +2d0*Bzero2(z2(i1,i2),ep)
     . +2d0*m0sq*Cv(cii(z2(i1,i2))+N0,ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_00ii',i1,i2,diff
     
      enddo
      
      return
      end
