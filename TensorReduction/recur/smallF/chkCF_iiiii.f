      subroutine chkCF_iiiii(i1,i2,i3,i4,i5,m0sq,Gr,Bzero5,N0)
C---  Expression for rearrangement of extension of Eq. 5.69
C---  Checks Ciiii, requires C00iiiii
C---  Small terms of order Gr(i,j)*Cijklmno
C---  Denominator m0sq
      implicit none
      include 'constants.f'
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,i2,i3,i4,i5,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero5(z5max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiiiiii(z7(n,m,i1,i2,i3,i4,i5))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czziiiii(z5(i1,i2,i3,i4,i5))+N0,ep-1)
      
      diff=Cv(ciiiii(z5(i1,i2,i3,i4,i5))+N0,ep)*2d0*m0sq-
     . (28d0*Cv(czziiiii(z5(i1,i2,i3,i4,i5))+N0,ep)
     . -pole
     . -2d0*Bzero5(z5(i1,i2,i3,i4,i5),ep)
     . +bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCF_iiiii',i1,i2,i3,i4,i5,diff
     
      enddo
      
      return
      end
