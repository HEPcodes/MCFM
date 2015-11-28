      subroutine chkCF_i(i1,m0sq,Gr,Bzero1,N0)
C---  Expression for rearrangement of Eq. 5.67
C---  Checks Ci, requires C00i
C---  Small terms of order Gr(i,j)*Cijk
C---  Denominator m0sq
      implicit none
      include 'constants.f'
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero1(z1max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czzi(i1)+N0,ep-1)
      
      diff=Cv(ci(i1)+N0,ep)*2d0*m0sq-      
     . (12d0*Cv(czzi(i1)+N0,ep)
     . -pole
     . -2d0*Bzero1(i1,ep)
     . +bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCF_i',i1,diff
     
      enddo
      
      return
      end
