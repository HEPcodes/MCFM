      subroutine chkF_00iiiii(i1,i2,i3,i4,i5,f,Gr,Shat7,N0)
      implicit none
      include 'Dnames.f'
      include 'Dv.f'
      include 'Darraydef.f'
      include 'Darrays.f'
      include 'weenumber.f' 
      integer ep,N0,k,i1,i2,i3,i4,i5,np
      parameter(np=3)
      double precision f(np),Gr(np,np),den
      double complex Shat7(np,z6max,-2:0),diff
       
      do ep=-2,0
      if     ((i1 .eq. i2) .and. (i1 .eq. i3) .and. (i1 .eq. i4)
     .  .and. (i1 .eq. i5)) then
        den=12d0
	k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i3) .and. (i1 .eq. i4)) then
        den=10d0
	k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i3) .and. (i1 .eq. i5)) then
        den=10d0
	k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i4) .and. (i1 .eq. i5)) then
        den=10d0
	k=i1
      elseif ((i1 .eq. i3) .and. (i1 .eq. i4) .and. (i1 .eq. i5)) then
        den=10d0
	k=i1
      elseif ((i2 .eq. i3) .and. (i2 .eq. i4) .and. (i2 .eq. i5)) then
        den=10d0
	k=i2
      elseif ((i1 .eq. i2) .and. (i1 .eq. i3)) then
        den=8d0
	k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i4)) then
        den=8d0
	k=i1
      elseif ((i1 .eq. i2) .and. (i1 .eq. i5)) then
        den=8d0
	k=i1
      elseif ((i2 .eq. i3) .and. (i2 .eq. i4)) then
        den=8d0
	k=i2
      elseif ((i2 .eq. i3) .and. (i2 .eq. i5)) then
        den=8d0
	k=i2
      elseif ((i3 .eq. i4) .and. (i3 .eq. i5)) then
        den=8d0
	k=i3
      elseif ((i1 .eq. i2) .or. (i1 .eq. i3) .or. (i1 .eq. i4)
     .   .or. (i1 .eq. i5)) then
        den=6d0
	k=i1
      elseif ((i2 .eq. i3) .or. (i2 .eq. i4) .or. (i2 .eq. i5)) then
        den=6d0
	k=i2
      elseif ((i3 .eq. i4) .or. (i3 .eq. i5)) then
        den=6d0
	k=i3
      elseif (i4 .eq. i5) then
        den=6d0
	k=i4
      else
        den=4d0
	k=i1
      endif      
      diff=Dv(dzziiiii(z5(i1,i2,i3,i4,i5))+N0,ep)*den-
     . (Shat7(k,z6(k,i1,i2,i3,i4,i5),ep)
     . -f(k)*Dv(diiiiii(z6(k,i1,i2,i3,i4,i5))+N0,ep)
     . -Gr(k,1)*Dv(diiiiiii(z7(1,k,i1,i2,i3,i4,i5))+N0,ep) 
     . -Gr(k,2)*Dv(diiiiiii(z7(2,k,i1,i2,i3,i4,i5))+N0,ep)
     . -Gr(k,3)*Dv(diiiiiii(z7(3,k,i1,i2,i3,i4,i5))+N0,ep)) 

      if ((abs(diff) .gt. weenumber)) 
     . write(6,*) 'chkF_00iiiii',i1,i2,i3,i4,i5,diff,
     . Dv(dzziiiii(z5(i1,i2,i3,i4,i5))+N0,ep)*den,
     . (Shat7(k,z6(k,i1,i2,i3,i4,i5),ep)
     . -f(k)*Dv(diiiiii(z6(k,i1,i2,i3,i4,i5))+N0,ep)
     . -Gr(k,1)*Dv(diiiiiii(z7(1,k,i1,i2,i3,i4,i5))+N0,ep) 
     . -Gr(k,2)*Dv(diiiiiii(z7(2,k,i1,i2,i3,i4,i5))+N0,ep)) 
     
      enddo
      
      return
      end
