      double complex function ntm(j1,j2,j3,j4,j5,j6,j7,za,zb)
      implicit none
      include 'constants.f'
      include 'prods.f'
      double complex function 
     . A234(2,2,2,2,2),A342(2,2,2,2,2),A423(2,2,2,2,2),
     . A432(2,2,2,2,2),A324(2,2,2,2,2),A243(2,2,2,2,2)
      integer j1,j2,j3,j4,j5,j6,j7
      call trocs(j1,j2,j3,j4,j5,j6,j7,za,zb,A234)
      call trocs(j1,j3,j4,j2,j5,j6,j7,za,zb,A342)
      call trocs(j1,j4,j2,j3,j5,j6,j7,za,zb,A423)
      call trocs(j1,j4,j3,j2,j5,j6,j7,za,zb,A432)
      call trocs(j1,j3,j2,j4,j5,j6,j7,za,zb,A324)
      call trocs(j1,j2,j4,j3,j5,j6,j7,za,zb,A243)
      return 
      end


      double precsion function ntM5sq(j1,j2,j3,j4,j5,j6,j7)
      implicit none
      include 'constants.f'
      integer j1,j2,j3,j4,j5,j6,j7
      double precision ntM0,ntM1,ntM2,x
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C  Eqn B30
      x=xn/cf
      ntm5sq=CF**3*XN*(ntM0-0.5d0*(ntM1+2d0*ntM0)
     . +0.25d0*x**2*(ntM0+ntM1+ntM2)
      return 
      end

      double precision function ntM0()
      implicit none
      include 'constants.f'
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C  Eqn B31
      double complex sum,ntm
      integer i2(6),i3(6),i4(6)
      data i2/2,3,4,4,3,2/
      data i2/2,3,4,3,2,4/
      data i2/2,3,4,2,4,3/
      sum=czip
      do k=1,6
      sum=sum+ntm(i2(k),i3(k),i4(k))
      enddo
      ntM0=dabs(sum)
      return 
      end
 

      double precision function ntM2()
      implicit none
      include 'constants.f'
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C  Eqn B32
      double complex ntm
      double precision sum
      integer i2(6),i3(6),i4(6)
      data i2/2,3,4,4,3,2/
      data i3/2,3,4,3,2,4/
      data i4/2,3,4,2,4,3/
      sum=zip
      do k=1,6
      sum=sum+dabs(ntm(1,i2(k),i3(k),i4(k)))**2
      enddo
      ntM2=sum
      return 
      end
 


      double precision function ntM1()
      implicit none
      include 'constants.f'
C--Results taken from Nagy and Trocsanyi-hep-ph/9806317
C  Eqn B33
      double complex ntm,sum
      double precision ntM2
      integer i2(3),i3(3),i4(3)
      data i2/2,3,4/
      data i3/3,4,2/
      data i4/4,2,3/
      sum=czip
      do k=1,3
      sum=sum+dconjg(ntm(i2(k),i3(k),i4(k)))
     . *(ntm(i2(k),i4(k),i3(k))
     .  +ntm(i3(k),i2(k),i4(k))
     .  -ntm(i4(k),i3(k),i2(k)))
      enddo
      ntM1=-2d0*(ntM2+dble(sum))
      return 
      end
 
