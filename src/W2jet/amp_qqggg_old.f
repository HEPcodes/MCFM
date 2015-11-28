      double complex function amp_qqggg(i1,h1,i2,h2,i3,h3,i4,h4,i5,lh)
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C--Appendix A
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'dprodx.f'
      include 'sprodx.f'
      integer i1,i2,i3,i4,i5,i6,i7,j,k,h1,h2,h3,h4,lh
      double complex t2,t3,tx,A(2,2,2,2,2)
      double complex xa(mxpart,mxpart),xb(mxpart,mxpart),t2a
      double precision t123,t234,t345,t567,t167

      t2a(i1,i2,i3,i4)=xa(i1,i2)*xb(i2,i4)+xa(i1,i3)*xb(i3,i4)
      
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
      h4=h4-(-1)**h4
      lh=lh-(-1)**lh
      do j=1,mxpart
      do k=1,mxpart
      xa(j,k)=zb(k,j)
      xb(j,k)=za(k,j)
      enddo
      enddo
       
      endif
      t123=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t167=s(i1,i6)+s(i6,i7)+s(i7,i1)
      t567=s(i5,i6)+s(i6,i7)+s(i7,i5)
      t345=s(i3,i4)+s(i4,i5)+s(i5,i3)
      t2a34=s(i2,i3)+s(i3,i4)+s(i4,i2)

C--A43 (2,2,2,2,2)
      if (
     . (h1.eq.2).and.(h2.eq.2).and.(h3.eq.2).and.(h4.eq.2).and.(lh.eq.2)
     . ) then 
      amp_qqggg =
     .-xa(i6,i5)**2*xb(i6,i7)/(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xa(i4,i5))

C--A49
C        A(2,1,1,1,2) =
      elseif 
     .((h1.eq.2).and.(h2.eq.1).and.(h3.eq.1).and.(h4.eq.1).and.(lh.eq.2)
     . ) then
      amp_qqggg =
     .-xb(i5,i7)**2*xa(i6,i7)/(xb(i2,i1)*xb(i3,i2)*xb(i4,i3)*xb(i5,i4))
 
C--A44
C      A(2,2,2,1,2) =
      elseif 
     .((h1.eq.2).and.(h2.eq.2).and.(h3.eq.2).and.(h4.eq.1).and.(lh.eq.2)
     . ) then
      amp_qqggg =
     . +xa(i6,i5)*t2a(i4,i5,i6,i7)/(xa(i2,i3)*xa(i3,i4)*xb(i3,i4)*t567)
     .*(xb(i2,i3)*t2a(i4,i2,i3,i1)/t2a34+t2a(i4,i1,i2,i3)/xa(i1,i2)) 
     . +t2a(i4,i5,i6,i7)*t2a(i6,i5,i4,i3)
     . /(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xb(i3,i4)*xb(i4,i5))
     . -t2a(i4,i1,i2,i7)*t2a(i6,i5,i4,i3)*xa(i4,i5)*xb(i5,i3)
     . /(xa(i1,i2)*xa(i2,i4)*xb(i4,i5)*s(i3,i4)*t345) 
     . -xb(i1,i7)*t2a(i6,i1,i7,i2)*xa(i4,i5)**2*xb(i5,i3)**2
     . /(xb(i4,i5)*xa(i4,i2)*s(i3,i4)*t345*t167) 
     . +(xb(i1,i7)*xa(i4,i5)*xb(i5,i3))
     . /(xa(i2,i3)*xb(i3,i4)*xb(i4,i5)*t167)
     . *(t2a(i6,i1,i7,i2)/xa(i3,i4)+t2a(i6,i1,i7,i3)/xa(i2,i4))
     . -xb(i1,i7)*xa(i4,i5)*xb(i2,i3)/(xa(i2,i3)*xb(i3,i4)*t2a34*t167)
     . *(t2a(i6,i1,i7,i2)*xa(i2,i4)/xa(i3,i4)+t2a(i6,i1,i7,i3)) 


C--A45
C      A(2,2,1,2,2)=
      elseif 
     .((h1.eq.2).and.(h2.eq.2).and.(h3.eq.1).and.(h4.eq.2).and.(lh.eq.2)
     . ) then
      amp_qqggg =
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i4)*t2a(i3,i5,i6,i7)*xa(i6,i5)
     . /(xa(i1,i2)*xa(i3,i4)*s(i2,i3)*t123*t567) 
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i7)*xa(i6,i5)*xa(i3,i5)
     . /(xa(i1,i2)*xa(i3,i4)*xa(i4,i5)*s(i2,i3)*t123)
     . -t2a(i3,i1,i2,i7)*t2a(i6,i1,i7,i2)*xa(i3,i5)**2
     . /(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xa(i4,i5)*xb(i2,i3)*t345) 
     . +t2a(i3,i1,i2,i7)*t2a(i6,i5,i3,i4)*xa(i3,i5)*xb(i4,i2)
     . /(xb(i2,i3)*xa(i1,i2)*xa(i2,i3)*s(i3,i4)*t345)
     . -(xb(i4,i2)**2*xb(i1,i7)*t3(i6,i1,i7,i2,i4,i3)*xa(i3,i5))
     . /(s(i2,i3)*s(i3,i4)*t2a34*t167) 
     . +(xb(i4,i2)*t2a(i3,i5,i6,i7)*xa(i6,i5))/(s(i2,i3)*s(i3,i4)*t567)
     . *((xb(i4,i2)*t2a(i3,i2,i4,i1))/t2a34-t2a(i3,i1,i2,i4)/xa(i1,i2))
     . +xb(i1,i7)*t2a(i6,i1,i7,i2)*xa(i3,i5)**2/(s(i2,i3)*t345*t167)
     . *(t2a(i3,i5,i4,i2)/(xa(i3,i4)*xa(i4,i5))
     . +xb(i4,i2)*xb(i5,i4)/s(i3,i4))
 
C---A46
C      A(2,1,2,2,2) =
      elseif 
     .((h1.eq.2).and.(h2.eq.1).and.(h3.eq.2).and.(h4.eq.2).and.(lh.eq.2)
     . ) then
      amp_qqggg =
     . +xb(i4,i3)**2*t2a(i2,i3,i4,i1)*t2a(i2,i5,i6,i7)*xa(i6,i5)
     . /(s(i2,i3)*s(i3,i4)*t2a34*t567)
     . +xb(i1,i3)*t2a(i2,i3,i4,i1)*t2a(i2,i5,i6,i7)*xa(i6,i5)
     . /(xb(i1,i2)*xa(i3,i4)*xa(i4,i2)*s(i2,i3)*t567)
     . -xb(i1,i3)**2*t2a(i2,i1,i3,i4)*t2a(i2,i5,i6,i7)*xa(i6,i5)
     . /(xb(i1,i2)*xa(i2,i4)*s(i2,i3)*t123*t567)
     . -xb(i1,i3)**2*t2a(i2,i1,i3,i7)*xa(i6,i5)*xa(i2,i5)
     . /(xb(i1,i2)*xa(i2,i4)*xa(i4,i5)*s(i2,i3)*t123)
     . -xb(i1,i7)*xa(i2,i5)*xb(i4,i3)**2*t3(i6,i1,i7,i3,i4,i2)
     . /(s(i2,i3)*s(i3,i4)*t2a34*t167) 
     . +xb(i1,i3)*xb(i1,i7)*xa(i2,i5)
     . /(xb(i1,i2)*s(i2,i3)*xa(i2,i4)*t345)
     . *(t2a(i6,i5,i4,i3)*(xa(i2,i5)/xa(i4,i5)-xa(i3,i2)/xa(i3,i4))
     . - t2a(i6,i5,i3,i4)*(xa(i4,i2))/xa(i3,i4)) 
     . +xb(i1,i7)*t2a(i6,i1,i7,i3)*xa(i2,i5)
     . /(xa(i2,i4)*s(i2,i3)*t345*t167)
     . *(t2a(i2,i5,i4,i3)
     . *((xa(i2,i5))/(xa(i4,i5))-(xa(i3,i2))/(xa(i3,i4)))
     . -t2a(i2,i5,i3,i4)*xa(i4,i2)/xa(i3,i4))
 
C--A47
C      A(2,2,1,1,2) =
      elseif 
     .((h1.eq.2).and.(h2.eq.2).and.(h3.eq.1).and.(h4.eq.1).and.(lh.eq.2)
     . ) then
      amp_qqggg =
     . +xb(i1,i2)*tx(i2,i3,i4,i5,i6,i7)*xa(i6,i5)/(s(i2,i3)*t567)
     . *(xa(i4,i3)**2/(s(i3,i4)*t2a34)
     . -xa(i3,i1)/(xb(i3,i4)*xb(i4,i2)*xa(i1,i2)))
     . -xa(i3,i1)**2*xb(i1,i2)**2*t2a(i4,i5,i6,i7)*xa(i6,i5)
     . /(xa(i1,i2)*xb(i2,i4)*s(i2,i3)*t123*t567)
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i7)*t2a(i6,i5,i4,i2)
     . /(xa(i1,i2)*xb(i2,i4)*xb(i4,i5)*s(i2,i3)*t123) 
     . +t2a(i3,i1,i2,i7)*t2a(i6,i1,i7,i2)
     . /(xa(i1,i2)*xb(i3,i4)*xb(i4,i5)*s(i2,i3))
     . +xb(i1,i7)*t2a(i6,i1,i7,i2)*t2a(i3,i5,i4,i2)
     . /(xb(i3,i4)*xb(i4,i5)*s(i2,i3)*t167)
     . -xb(i1,i7)*t2a(i6,i1,i7,i2)*t2a(i5,i3,i4,i2)*xa(i4,i3)**2
     . /(s(i2,i3)*s(i3,i4)*t2a34*t167)

C--A48
C      A(2,1,1,2,2) =
      elseif 
     .((h1.eq.2).and.(h2.eq.1).and.(h3.eq.1).and.(h4.eq.2).and.(lh.eq.2)
     . ) then
      amp_qqggg =
     . +xa(i2,i3)**2*xb(i1,i4)*tx(i4,i2,i3,i5,i6,i7)*xa(i6,i5)
     . /(s(i2,i3)*s(i3,i4)*t2a34*t567)
     . -xb(i1,i4)*t2a(i3,i5,i6,i7)*xa(i6,i5)
     . /(xb(i4,i2)*s(i3,i4)*t123*t567)
     . *((xb(i4,i2)*t2a(i2,i1,i3,i4)+xb(i4,i3)*t2a(i3,i1,i2,i4))
     . /xb(i2,i3) 
     . -xb(i1,i4)*t2a(i3,i1,i2,i4)/xb(i1,i2)) 
     . -xb(i1,i4)*xa(i6,i5)*xa(i3,i5)
     . /(xb(i4,i2)*xa(i4,i5)*s(i3,i4)*t123)
     . *((xb(i4,i2)*t2a(i2,i1,i3,i7)+xb(i4,i3)*t2a(i3,i1,i2,i7))
     . /xb(i2,i3)
     .  -xb(i1,i4)*t2a(i3,i1,i2,i7)/xb(i1,i2)) 
     . +xb(i1,i4)*xb(i1,i7)*t2a(i6,i5,i3,i4)*xa(i3,i5)**2
     . /(xb(i1,i2)*xb(i2,i4)*xa(i4,i5)*s(i3,i4)*t345)
     . -xb(i1,i7)*t2a(i6,i1,i7,i4)*t2a(i2,i5,i3,i4)*xa(i3,i5)**2
     . /(xb(i4,i2)*xa(i4,i5)*s(i3,i4)*t345*t167)
     . -xb(i1,i7)*t2a(i6,i1,i7,i4)*t2a(i5,i2,i3,i4)*xa(i3,i5)
     . /(xa(i4,i5)*xb(i4,i2)*xb(i2,i3)*s(i3,i4)*t167)
     . -xb(i1,i7)*t2a(i6,i1,i7,i4)*t2a(i5,i2,i3,i4)*xa(i2,i3)**2
     . /(s(i2,i3)*s(i3,i4)*t2a34*t167)
 

      elseif 
     .((h1.eq.2).and.(h2.eq.1).and.(h3.eq.2).and.(h4.eq.1).and.(lh.eq.2)
     . ) then

C--A49
C----Constructed from A45 with P operation (above) and 6 and 7 exchanged
C      A(2,1,2,1,2)=
      amp_qqggg =
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i4)*t2a(i3,i5,i7,i6)*xa(i7,i5)
     . /(xa(i1,i2)*xa(i3,i4)*s(i2,i3)*t123*t567) 
     . +xa(i3,i1)*xb(i1,i2)*t2a(i3,i1,i2,i6)*xa(i7,i5)*xa(i3,i5)
     . /(xa(i1,i2)*xa(i3,i4)*xa(i4,i5)*s(i2,i3)*t123)
     . -t2a(i3,i1,i2,i6)*t2a(i7,i1,i6,i2)*xa(i3,i5)**2
     . /(xa(i1,i2)*xa(i2,i3)*xa(i3,i4)*xa(i4,i5)*xb(i2,i3)*t345) 
     . +t2a(i3,i1,i2,i6)*t2a(i7,i5,i3,i4)*xa(i3,i5)*xb(i4,i2)
     . /(xb(i2,i3)*xa(i1,i2)*xa(i2,i3)*s(i3,i4)*t345)
     . -(xb(i4,i2)**2*xb(i1,i6)*t3(i7,i1,i6,i2,i4,i3)*xa(i3,i5))
     . /(s(i2,i3)*s(i3,i4)*t2a34*t167) 
     . +(xb(i4,i2)*t2a(i3,i5,i7,i6)*xa(i7,i5))/(s(i2,i3)*s(i3,i4)*t567)
     . *((xb(i4,i2)*t2a(i3,i2,i4,i1))/t2a34-t2a(i3,i1,i2,i4)/xa(i1,i2))
     . +xb(i1,i6)*t2a(i7,i1,i6,i2)*xa(i3,i5)**2/(s(i2,i3)*t345*t167)
     . *(t2a(i3,i5,i4,i2)/(xa(i3,i4)*xa(i4,i5))
     . +xb(i4,i2)*xb(i5,i4)/s(i3,i4))
 
      endif

      return 
      end

