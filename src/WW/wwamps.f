      subroutine wwamps(j1,j2,j3,j4,j5,j6,j7,za,zb,f)
c  -first label of fs,ft is gluon polarization, second is qqb line

      implicit none
      include 'constants.f'
      include 'dprodx.f'
      include 'sprodx.f'
      include 'zerowidth.f'
      include 'masses.f'
      integer j,k,jtype,j1,j2,j3,j4,j5,j6,j7,mplus,minus
      double complex A7treea,A7treeb,B7treea,B7treeb
      double complex f(4,2,2)
      double complex prop34,prop56,propboth
      data minus,mplus/1,2/
      
c----initialize to zero
      do jtype=3,4
      do j=1,2
      do k=1,2
      f(jtype,j,k)=czip
      enddo
      enddo
      enddo
            
      if     (zerowidth  .eqv. .true.) then
      prop34=dcmplx(s(3,4))/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
      prop56=dcmplx(s(5,6))/dcmplx(s(5,6)-wmass**2,wmass*wwidth)
      elseif (zerowidth .neqv. .true.) then
      prop34=dcmplx(s(3,4)/(s(3,4)-wmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-wmass**2)) 
      endif
      propboth=prop34*prop56

      f(1,mplus,mplus)= czip
      f(2,mplus,mplus)=-A7treeb(j2,j1,j3,j4,j5,j6,j7,za,zb)*propboth

      f(1,mplus,minus)=+A7treea(j1,j2,j3,j4,j5,j6,j7,za,zb)*propboth
      f(2,mplus,minus)=+A7treeb(j1,j2,j3,j4,j5,j6,j7,za,zb)*propboth

      f(1,minus,mplus)= czip
      f(2,minus,mplus)=+A7treeb(j1,j2,j5,j6,j3,j4,j7,zb,za)*propboth

      f(1,minus,minus)=-A7treea(j2,j1,j5,j6,j3,j4,j7,zb,za)*propboth
      f(2,minus,minus)=-A7treeb(j2,j1,j5,j6,j3,j4,j7,zb,za)*propboth

      if (zerowidth) return   ! Done all amplitudes needed for zerowidth


      f(3,mplus,mplus)=-B7treeb(j1,j2,3,4,5,6,j7,za,zb)*prop34
     .                 -B7treea(j2,j1,3,4,5,6,j7,za,zb)*prop56
      f(4,mplus,mplus)=-B7treea(j2,j1,6,5,4,3,j7,za,zb)*prop34
     .                 -B7treeb(j1,j2,6,5,4,3,j7,za,zb)*prop56

      f(3,mplus,minus)=+B7treeb(j2,j1,3,4,5,6,j7,za,zb)*prop34
     .                 +B7treea(j1,j2,3,4,5,6,j7,za,zb)*prop56
      f(4,mplus,minus)=+B7treea(j1,j2,6,5,4,3,j7,za,zb)*prop34
     .                 +B7treeb(j2,j1,6,5,4,3,j7,za,zb)*prop56

      f(3,minus,mplus)=+B7treea(j1,j2,5,6,3,4,j7,zb,za)*prop34
     .                 +B7treeb(j2,j1,5,6,3,4,j7,zb,za)*prop56
      f(4,minus,mplus)=+B7treeb(j2,j1,4,3,6,5,j7,zb,za)*prop34
     .                 +B7treea(j1,j2,4,3,6,5,j7,zb,za)*prop56

      f(3,minus,minus)=-B7treea(j2,j1,5,6,3,4,j7,zb,za)*prop34
     .                 -B7treeb(j1,j2,5,6,3,4,j7,zb,za)*prop56
      f(4,minus,minus)=-B7treeb(j1,j2,4,3,6,5,j7,zb,za)*prop34
     .                 -B7treea(j2,j1,4,3,6,5,j7,zb,za)*prop56

      return
      end

