      double complex function amp_qqgg(i1,h1,i2,h2,i3,h3,i4,lh)
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C--Appendix A
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'prods.f'
      integer i1,i2,i3,i4,i5,i6,i7,j,k,h1,h2,h3,h4,lh
      double complex t2,t3,tx
      double complex xa(mxpart,mxpart),xb(mxpart,mxpart)
      double precision t123,t234
C---A(hq,h2,h3,h4,lh)      
C---h=1 LH
C---h=2 RH

      if ((h1.eq.2) .and. (lh.eq.1)) then
      i6=7
      i7=6
      lh=2  
      elseif ((h1.eq.1) .and. (lh.eq.2)) then
      i6=7
      i7=6
      lh=1  
      endif 

      if (h1.eq.2) then
      do j=1,mxpart
      do k=1,mxpart
      xa(j,k)=za(j,k)
      xb(j,k)=zb(j,k)
      enddo
      enddo
      elseif (h1.eq.1) then
      h1=2
      h2=h2-(-1)**h2
      h3=h3-(-1)**h3
      lh=lh-(-1)**lh
      do j=1,mxpart
      do k=1,mxpart
      xa(j,k)=zb(k,j)
      xb(j,k)=za(k,j)
      enddo
      enddo
       
      endif


      t123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t234=s(i2,i3)+s(i3,i4)+s(i4,i2)
C--A37
c      A(2,2,2,2) =
      if ((h1.eq.2).and.(h1.eq.2).and.(h1.eq.2).and.(h1.eq.2)) then
      amp_qqgg=
     .-za(i4,i5)**2*zb(i5,i6)/(za(i1,i2)*za(i2,i3)*za(i3,i4))
 
C--A38
c      A(2,1,1,2) =
      elseif ((h1.eq.2).and.(h1.eq.1).and.(h1.eq.1).and.(h1.eq.2)) then
      amp_qqgg=
     . -zb(i1,i6)**2*za(i5,i6)/(zb(i1,i2)*zb(i2,i3)*zb(i3,i4))


C--A39
c      A(2,2,1,2)=
      elseif ((h1.eq.2).and.(h1.eq.2).and.(h1.eq.1).and.(h1.eq.2)) then
      amp_qqgg=
     . -za(i3,i1)*zb(i1,i2)*za(i4,i5)*t2(i3,i1,i2,i6)
     . /(za(i1,i2)*s(i2,i3)*t123) 
     . +za(i3,i4)*zb(i4,i2)*zb(i1,i6)*t2(i5,i3,i4,i2)
     . /(zb(i3,i4)*s(i2,i3)*t234)
     . +t2(i5,i3,i4,i2)*t2(i3,i1,i2,i6)
     . /(za(i1,i2)*zb(i3,i4)*s(i2,i3)) 
 
C---A40
c      A(2,1,2,2) =
      elseif ((h1.eq.2).and.(h1.eq.1).and.(h1.eq.2).and.(h1.eq.2)) then
      amp_qqgg=
     . +zb(i1,i3)**2*zb(i4,i5)*t2(i2,i1,i3,i6)
     . /(zb(i1,i2)*s(i2,i3)*t123)
     . -za(i2,i4)**2*zb(i1,i6)*t2(i5,i2,i4,i3)
     . /(za(i3,i4)*s(i2,i3)*t234)
     . -zb(i1,i3)*za(i2,i4)*zb(i1,i6)*zb(i4,i5)
     . /(zb(i1,i2)*za(i3,i4)*s(i2,i3))
      endif
      return 
      end

