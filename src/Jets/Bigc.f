      double precision function Bigc(ia,ib,i1,i2,i3)
      implicit none
C     Expressions taken from Berends,Kleiss,DeC,Gastmans,Wu
C     Physics Letters 103B 124, (1981)
C     Eqn (16)
      include 'constants.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      integer ia,ib,i1,i2,i3,j
      double precision a(3),b(3)
      double precision xnum

      a(1)=0.5d0*s(ia,i1)
      a(2)=0.5d0*s(ia,i2)
      a(3)=0.5d0*s(ia,i3)
      b(1)=0.5d0*s(ib,i1)
      b(2)=0.5d0*s(ib,i2)
      b(3)=0.5d0*s(ib,i3)
      Bigc=0d0
      do j=1,3
      Bigc=Bigc+a(j)*b(j)*(a(j)**2+b(j)**2)
      enddo
      Bigc=Bigc/(a(1)*a(2)*a(3)*b(1)*b(2)*b(3))
      xnum=0.5d0*s(ia,ib)
     . +xn**2*(0.5d0*s(ia,ib)
     . -2d0*(a(1)*b(2)+a(2)*b(1))/s(i1,i2)
     . -2d0*(a(2)*b(3)+a(3)*b(2))/s(i2,i3)
     . -2d0*(a(3)*b(1)+a(1)*b(3))/s(i3,i1))
     . +2d0*xn**4/s(ia,ib)*(
     . +4d0*a(3)*b(3)*(a(1)*b(2)+a(2)*b(1))/(s(i2,i3)*s(i3,i1))
     . +4d0*a(1)*b(1)*(a(2)*b(3)+a(3)*b(2))/(s(i3,i1)*s(i1,i2))
     . +4d0*a(2)*b(2)*(a(3)*b(1)+a(1)*b(3))/(s(i1,i2)*s(i2,i3)))
      Bigc=gsq**3*V/xn**2*xnum*Bigc
      return
      end

