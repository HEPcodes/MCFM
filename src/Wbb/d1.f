      double complex function d1pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
C     DENOMS (dg),(bc),(neu)
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2345
      double complex t2
      s2345=(s(j2,j3)+s(j2,j4)+s(j2,j5)+s(j3,j4)+s(j3,j5)+s(j4,j5))
      d1pm=-fourrt2*za(jb,j2)*za(j4,j2)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,jb)*za(j3,j2)*s(j4,j5)*s2345)
      return
      end

      double complex function d1mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s2345
      s2345=(s(j2,j3)+s(j2,j4)+s(j2,j5)+s(j3,j4)+s(j3,j5)+s(j4,j5))
      d1mm=fourrt2*t2(jb,j2,j3,j4)*t2(j5,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,jb)*zb(j3,j2)*s(j4,j5)*s2345)
      return
      end

      double complex function d1pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d1pm
      d1pp=d1pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d1mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d1mm
      d1mp=d1mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end









