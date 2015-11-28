      double complex function t1(j1,j2,j3)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3
      t1=zb(j1,j2)*za(j3,j2)
      return
      end

      double complex function t2(j1,j2,j3,j4)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4
      t2=zb(j1,j2)*za(j4,j2)+zb(j1,j3)*za(j4,j3)
      return
      end

      double complex function t3(j1,j2,j3,j4,j5)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4,j5
      t3=zb(j1,j2)*za(j5,j2)+zb(j1,j3)*za(j5,j3)+zb(j1,j4)*za(j5,j4)
      return
      end

      double complex function t4(j1,j2,j3,j4,j5,j6)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4,j5,j6
      t4=zb(j1,j2)*za(j6,j2)+zb(j1,j3)*za(j6,j3)
     &  +zb(j1,j4)*za(j6,j4)+zb(j1,j5)*za(j6,j5)
      return
      end
