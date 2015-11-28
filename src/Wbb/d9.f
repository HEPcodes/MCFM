      double complex function d9mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s2345,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d9mp=fourrt2*za(j5,j2)*zb(j7,j1)*zb(j4,jb)
     &*t2(j4,j1,j7,j6)/(zb(j3,jb)*zb(j4,j3)*s345*s2345)
      return
      end

      double complex function d9pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2345,s345
      double complex t2
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d9pm=fourrt2*za(jb,j4)*za(j2,j4)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,jb)*za(j3,j4)*s345*s2345)
      return
      end

      double complex function d9mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s2345,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d9mm=-fourrt2*t2(jb,j3,j4,j2)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,jb)*zb(j3,j4)*s345*s2345)
      return
      end

      double complex function d9pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2345,s345
      double complex t2
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d9pp=-fourrt2*za(j5,j2)*zb(j7,j1)
     & *(za(j3,jb)*t2(j3,j1,j7,j6)+za(j4,jb)*t2(j4,j1,j7,j6))
     & /(za(j3,jb)*za(j4,j3)*s345*s2345)
      return
      end

