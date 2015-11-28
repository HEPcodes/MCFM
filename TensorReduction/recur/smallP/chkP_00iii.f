      subroutine chkP_00iii(i1,i2,i3,m0sq,Gr,Czero3,N0)
      implicit none
      include 'constants.f'
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,i2,i3,m,n,np
      parameter(np=3)
      double precision m0sq,Gr(np,np)
      double complex Czero3(z3max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diiiii(z5(n,m,i1,i2,i3))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzziii(z3(i1,i2,i3))+N0,ep-1)
      
      diff=Dv(dzziii(z3(i1,i2,i3))+N0,ep)*20d0-
     . (pole
     . +2d0*Czero3(z3(i1,i2,i3),ep)
     . +2d0*m0sq*Dv(diii(z3(i1,i2,i3))+N0,ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkP_00iii',i1,i2,i3,diff
     
      enddo
      
      return
      end
