      subroutine chkCP_0000i(i1,Gr,S0000i,N0)
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
      double complex S0000i(np,-2:0),bit,pole,diff
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(czziii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czzzzi(i1)+N0,ep-1)

c--- note: we have simplified the recursion relation by using DD (5.9)
c--- so that S0000i = 2B00i(0)+2*m0sq*C00i      
      diff=Cv(czzzzi(i1)+N0,ep)*16d0-
     . (pole
     . +S0000i(i1,ep)
     . -bit)

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkCP_0000i',i1,diff
     
      enddo
      
      return
      end
