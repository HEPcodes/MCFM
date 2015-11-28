      double complex function d10pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2345,s345
      double complex t2
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d10pm=fourrt2*za(j2,j4)*zb(j7,j1)
     & *(za(j6,j1)*t2(j1,j3,j5,jb)+za(j6,j7)*t2(j7,j3,j5,jb))
     & /(za(j3,jb)*za(j3,j5)*s345*s2345)
      return
      end

      double complex function d10mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2345,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d10mm=-fourrt2*za(j2,j4)*zb(jb,j5)*zb(j7,j1)
     & *(za(j6,j1)*zb(j1,j5)+za(j6,j7)*zb(j7,j5))
     & /(zb(j3,jb)*zb(j3,j5)*s345*s2345)
      return
      end

      double complex function d10pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2345,s345
      double complex t2
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d10pp=fourrt2*za(j5,j2)*za(j5,jb)*zb(j7,j1)
     & *t2(j4,j1,j7,j6)/(za(j3,jb)*za(j5,j3)*s345*s2345)
      return
      end

      double complex function d10mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s2345,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      d10mp=-fourrt2*zb(j7,j1)*t2(j4,j1,j7,j6)
     & *(za(j5,j2)*zb(j5,jb)+za(j3,j2)*zb(j3,jb))
     & /(zb(j3,jb)*zb(j5,j3)*s345*s2345)
      return
      end


