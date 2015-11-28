      subroutine chkF_i(i1,m0sq,Gr,Czero1,N0)
C---  Expression for rearrangement of Eq. 5.67
C---  Checks Di, requires D00i
C---  Small terms of order Gr(i,j)*Dijk
C---  Denominator m0sq
      implicit none
      include 'constants.f'
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex Czero1(z1max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzzi(i1)+N0,ep-1)
      
      diff=Dv(di(i1)+N0,ep)*2d0*m0sq-      
     . (12d0*Dv(dzzi(i1)+N0,ep)
     . -pole
     . -2d0*Czero1(i1,ep)
     . +bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkF_i',i1,diff
     
      enddo
      
      return
      end
