      subroutine chkP_00ii(i1,i2,m0sq,Gr,Czero2,N0)
      implicit none
      include 'constants.f'
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,i2,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex Czero2(z2max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diiii(z4(n,m,i1,i2))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzzii(z2(i1,i2))+N0,ep-1)
      
      diff=Dv(dzzii(z2(i1,i2))+N0,ep)*16d0-
     . (pole
     . +2d0*Czero2(z2(i1,i2),ep)
     . +2d0*m0sq*Dv(dii(z2(i1,i2))+N0,ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_00ii',i1,i2,diff
     
      enddo
      
      return
      end
