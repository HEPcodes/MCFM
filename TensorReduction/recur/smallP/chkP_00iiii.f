      subroutine chkP_00iiii(i1,i2,i3,i4,m0sq,Gr,Czero4,N0)
      implicit none
      include 'constants.f'
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,i2,i3,i4,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex Czero4(z4max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diiiiii(z6(n,m,i1,i2,i3,i4))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep-1)
      
      diff=Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep)*24d0-
     . (pole
     . +2d0*Czero4(z4(i1,i2,i3,i4),ep)
     . +2d0*m0sq*Dv(diiii(z4(i1,i2,i3,i4))+N0,ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_00iiii',i1,i2,i3,i4,diff
     
      enddo
      
      return
      end
