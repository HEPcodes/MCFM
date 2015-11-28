      subroutine chkCP_00iii(i1,i2,i3,m0sq,Gr,Bzero3,N0)
      implicit none
      include 'constants.f'
      include 'Cnames.f'
      include 'Cv.f'
      include 'Carraydef.f'
      include 'Carrays.f'
      include 'weenumber.f' 
      integer ep,N0,i1,i2,i3,m,n,np
      parameter(np=2)
      double precision m0sq,Gr(np,np)
      double complex Bzero3(z3max,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiiii(z5(n,m,i1,i2,i3))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czziii(z3(i1,i2,i3))+N0,ep-1)
      
      diff=Cv(czziii(z3(i1,i2,i3))+N0,ep)*20d0-
     . (pole
     . +2d0*Bzero3(z3(i1,i2,i3),ep)
     . +2d0*m0sq*Cv(ciii(z3(i1,i2,i3))+N0,ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_00iii',i1,i2,i3,diff
     
      enddo
      
      return
      end
