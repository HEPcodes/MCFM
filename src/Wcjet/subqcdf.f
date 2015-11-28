      subroutine subqcdf(i1,i2,i3,i4,i5,i6,i7,p156,p256,za,zb,amp)
c*******************************************************************
c     the matrix elements of the
C     helicity amplitudes for the QCD process
c     s(-p1)+cbar(-p2) --> l(p3)+abar(p4)+g(p5)+g(p6)
c     multiplied by ((a+l)^2-M**2)/g^4/gwsq^2/2
c*******************************************************************
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'masses.f'
      integer i1,i2,i3,i4,i5,i6,i7,j1,j2,j3
      double precision p156,p256,s34,s56
C     p156=2*p1.p5+2*p1.p6+2*p5.p6
C     p256=2*p2.p5+2*p2.p6+2*p5.p6
      double complex amp(2,2,2)

      s34=dble(za(i3,i4)*zb(i4,i3))
      s56=dble(za(i5,i6)*zb(i6,i5))

      do j1=1,2
      do j2=1,2
      do j3=1,2
      amp(j1,j2,j3)=0d0
      enddo
      enddo
      enddo

  
      amp(1,1,2)= four*mc/p256/za(i3,i4)/zb(i3,i4) 
     . /zb(i5,i6)/zb(i1,i6)/zb(i2,i5)**2*zb(i1,i4)*(
     . -za(i2,i5)*za(i3,i5)*zb(i1,i2)*zb(i2,i5)
     . -za(i2,i5)*za(i3,i6)*zb(i1,i2)*zb(i2,i6)
     . +za(i3,i4)*za(i5,i6)*zb(i1,i4)*zb(i2,i6))

      amp(1,1,1)=four/p256*(zb(i1,i4)**2/zb(i3,i4)
     . *(za(i5,i6)*zb(i6,i5)+za(i5,i2)*zb(i2,i5)+za(i6,i2)*zb(i2,i6))
     . +mc**2/s34*zb(i1,i4)*zb(i1,i5)/zb(i2,i5)
     . *(za(i3,i5)*zb(i5,i2)+za(i3,i6)*zb(i6,i2)))
     . /zb(i5,i6)/zb(i1,i6)/zb(i2,i5)

      amp(2,2,1)=four*za(i2,i3)*(za(i2,i3)/za(i3,i4)
     . -mc**2/s34/p256/zb(i2,i5)*(za(i1,i6)*zb(i1,i4)*zb(i5,i6))
     . -mc**2/s34*zb(i4,i5)/zb(i2,i5)
     . )/za(i5,i6)/za(i1,i6)/za(i2,i5)

      amp(2,2,2)=four*mc*((za(i3,i5)*za(i2,i3)/za(i3,i4)
     . -mc**2/(s34*zb(i2,i5))*za(i3,i5)*zb(i4,i5))
     . /za(i1,i6)/za(i2,i5)/za(i5,i6)/zb(i2,i5)
     . +za(i2,i3)*zb(i1,i4)*(
     . +za(i2,i5)*(za(i3,i4)*zb(i4,i5)+za(i3,i1)*zb(i1,i5))
     . +za(i2,i6)*(za(i3,i1)*zb(i1,i6)+za(i3,i4)*zb(i4,i6))
     . -za(i2,i3)
     . *(za(i2,i5)*zb(i2,i5)+za(i2,i6)*zb(i2,i6)+za(i5,i6)*zb(i5,i6)))
     . /p256/za(i5,i6)/za(i2,i6)/za(i2,i3)
     . /zb(i2,i5)/zb(i3,i4)/za(i3,i4))

      amp(1,2,1)=four*(mc**2/s34/p256*(za(i3,i5)*zb(i1,i4)*zb(i2,i6))
     . /za(i5,i6)/zb(i2,i5)/zb(i2,i5)
     . +za(i2,i5)*zb(i2,i6)*zb(i1,i4)
     . *(za(i3,i1)*zb(i1,i6)+za(i3,i4)*zb(i4,i6))/zb(i2,i5)/p256/s34/s56
     . -za(i1,i5)*za(i2,i3)*zb(i1,i6)        
     . *(za(i2,i5)*zb(i2,i4)+zb(i3,i4)*za(i3,i5))/za(i1,i6)/p156/s34/s56
     . -(zb(i1,i4)/za(i5,i6)
     . *(-za(i1,i5)*za(i3,i5)+za(i1,i5)*za(i2,i3)*zb(i2,i6)/zb(i5,i6))
     . -(zb(i5,i6)*za(i3,i5)-za(i2,i3)*zb(i2,i6))*zb(i4,i6)/zb(i5,i6)
     . )/za(i1,i6)/zb(i2,i5)/s34);



      amp(1,2,2)=mc*four*za(i3,i5)*zb(i2,i6)/zb(i2,i5)/s34/s56
     . *(+za(i2,i5)*zb(i1,i4)*zb(i2,i6)/zb(i2,i5)/p256
     . -za(i1,i5)*zb(i1,i6)/za(i1,i6)/zb(i2,i6)
     . *(za(i2,i5)*zb(i2,i4)+za(i3,i5)*zb(i3,i4))/p156
     . +(za(i1,i5)*zb(i1,i4)+(za(i5,i6))*zb(i4,i6))
     . /za(i1,i6)/zb(i2,i5))


      amp(2,1,2)=mc*four*(-za(i2,i6)*zb(i1,i4)/s34/p256*(
     .  +(za(i1,i3)*zb(i1,i5)+za(i4,i3)*zb(i4,i5))
     .  /za(i2,i5)/zb(i2,i5)/zb(i5,i6)
     .  +za(i3,i6)/za(i5,i6)/zb(i5,i6))


     .  +zb(i1,i4)*zb(i1,i5)/s34/zb(i1,i6)/zb(i2,i5)/zb(i5,i6)
     .  *(za(i2,i3)*za(i5,i6)-za(i3,i6)*za(i2,i5))/za(i5,i6)/za(i2,i5)

     .  +za(i3,i5)*zb(i1,i5)**2
     .  /p156/s34/s56/zb(i1,i6)/zb(i2,i5)
     .  *(-za(i2,i6)*zb(i2,i4)-zb(i3,i4)*za(i3,i6))

     .  -za(i3,i5)*zb(i1,i5)**2
     .  /p156/s34/s56/zb(i1,i6)/zb(i2,i5)
     .  *mc**2*zb(i4,i5)*za(i2,i5)*za(i5,i6)/zb(i2,i5))

      amp(2,1,2)=mc*four*(za(i2,i6)*zb(i1,i4)/s34/s56/p256*(
     . +(za(i1,i3)*zb(i1,i5)+za(i4,i3)*zb(i4,i5))
     . /za(i2,i5)/zb(i2,i5)*za(i5,i6)+za(i3,i6))

     . -za(i3,i5)*zb(i1,i5)**2/p156/s34
     . /za(i5,i6)/zb(i1,i6)/zb(i2,i5)/zb(i5,i6)
     . *(-za(i2,i6)*zb(i2,i4)-zb(i3,i4)*za(i3,i6))

     . +zb(i1,i4)*zb(i1,i5)/s34/zb(i1,i6)/zb(i2,i5)/zb(i5,i6)
     . *(za(i2,i3)*za(i5,i6)-za(i3,i6)*za(i2,i5))/za(i5,i6)/za(i2,i5)

     . +mc**2/p156/s34*za(i3,i5)*zb(i1,i5)**2*zb(i4,i5)
     . /za(i2,i5)/zb(i1,i6)/zb(i2,i5)/zb(i2,i5)/zb(i5,i6))


      amp(2,1,1)=four*(mc**2*za(i2,i3)*zb(i1,i5)**2*zb(i4,i5)/p156/s34
     . /za(i2,i5)/zb(i2,i5)/zb(i1,i6)/zb(i5,i6)
     . -za(i2,i3)*zb(i1,i5)**2/p156/s34/s56/zb(i1,i6)
     . *(za(i2,i6)*zb(i2,i4)+za(i3,i6)*zb(i3,i4))
     . -za(i2,i6)**2*zb(i1,i4)/za(i2,i5)/p256/s34/s56
     . *(za(i1,i3)*zb(i1,i5)-zb(i4,i5)*za(i3,i4))
     . +za(i2,i3)*za(i2,i6)*zb(i1,i4)*zb(i1,i5)
     . /s34/s56/za(i2,i5)/zb(i1,i6));

      return
      end
