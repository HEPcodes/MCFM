      subroutine wzamps(j1,j2,j3,j4,j5,j6,j7,za,zb,f)
c  -first label of fs,ft,fu, is gluon polarization, second is zdecay line

      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'zerowidth.f'
      integer j1,j2,j3,j4,j5,j6,j7,mplus,minus,jtype,j,k
      double complex A7treea,A7treeb,B7treea,B7treeb
      double complex f(9,2,2)
      data minus,mplus/1,2/

c----initialize to zero
      do jtype=4,9
      do j=1,2
      do k=1,2
      f(jtype,j,k)=dcmplx(0d0)
      enddo
      enddo
      enddo

      f(1,mplus,minus)=+A7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)    !fs
      f(2,mplus,minus)=+A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)    !ft
      f(3,mplus,minus)=+A7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)    !fu

      f(1,mplus,mplus)=+A7treeb(j1,j2,j3,j4,j6,j5,j7,za,zb)
      f(2,mplus,mplus)=+A7treea(j1,j2,j3,j4,j6,j5,j7,za,zb)
      f(3,mplus,mplus)=+A7treea(j1,j2,j5,j6,j4,j3,j7,za,zb)

      f(1,minus,minus)=-A7treeb(j2,j1,j5,j6,j3,j4,j7,zb,za)
      f(2,minus,minus)=-A7treea(j2,j1,j5,j6,j3,j4,j7,zb,za)
      f(3,minus,minus)=-A7treea(j2,j1,j4,j3,j6,j5,j7,zb,za)

      f(1,minus,mplus)=-A7treeb(j2,j1,j6,j5,j3,j4,j7,zb,za)
      f(2,minus,mplus)=-A7treea(j2,j1,j6,j5,j3,j4,j7,zb,za)
      f(3,minus,mplus)=-A7treea(j2,j1,j4,j3,j5,j6,j7,zb,za)

      if (zerowidth) return

c+++ the extra amplitude 7 accounts for W+ production instead of W-
c+++ amplitudes 8 and 9 are needed when considering ub+d instead of
c+++ d+ub and are obtained from 6 and 7 by simply swapping 1 <+> 2

      f(4,mplus,minus)=+B7treeb(j2,j1,j6,j5,j4,j3,j7,za,zb)    
      f(5,mplus,minus)=+B7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)    
      f(6,mplus,minus)=+B7treeb(j2,j1,j3,j4,j5,j6,j7,za,zb)    
      f(7,mplus,minus)=+B7treea(j1,j2,j6,j5,j4,j3,j7,za,zb)    
      f(8,mplus,minus)=+B7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)    
      f(9,mplus,minus)=+B7treea(j2,j1,j6,j5,j4,j3,j7,za,zb)    

      f(4,mplus,mplus)=+B7treeb(j2,j1,j5,j6,j4,j3,j7,za,zb)
      f(5,mplus,mplus)=+B7treea(j1,j2,j3,j4,j6,j5,j7,za,zb)
         
      f(4,minus,minus)=-B7treea(j2,j1,j4,j3,j6,j5,j7,zb,za)
      f(5,minus,minus)=-B7treeb(j1,j2,j5,j6,j3,j4,j7,zb,za)
      f(6,minus,minus)=-B7treea(j2,j1,j5,j6,j3,j4,j7,zb,za)
      f(7,minus,minus)=-B7treeb(j1,j2,j4,j3,j6,j5,j7,zb,za)
      f(8,minus,minus)=-B7treea(j1,j2,j5,j6,j3,j4,j7,zb,za)
      f(9,minus,minus)=-B7treeb(j2,j1,j4,j3,j6,j5,j7,zb,za)

      f(4,minus,mplus)=-B7treea(j2,j1,j3,j4,j6,j5,j7,zb,za)
      f(5,minus,mplus)=-B7treeb(j1,j2,j5,j6,j4,j3,j7,zb,za)
     
      return
      end
