      double complex function d2pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
C     DENOMS (dg) (bc) (nedg)
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s2367
      double complex t2
      s2367=(s(j2,j3)+s(j2,j6)+s(j2,j7)+s(j3,j6)+s(j3,j7)+s(j6,j7))
      d2pm=-fourrt2*za(jb,j2)*za(j6,j2)*t2(j7,j1,j5,j4)
     & *zb(j5,j1)/(za(j3,jb)*za(j3,j2)*s(j4,j5)*s2367)
      return
      end

      double complex function d2mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s2367
      s2367=(s(j2,j3)+s(j2,j6)+s(j2,j7)+s(j3,j6)+s(j3,j7)+s(j6,j7))
      d2mm=fourrt2*t2(jb,j2,j3,j6)*t2(j7,j1,j5,j4)
     & *zb(j5,j1)/(zb(j3,jb)*zb(j3,j2)*s(j4,j5)*s2367)
      return
      end

      double complex function d2pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d2pm
      d2pp=d2pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d2mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d2mm
      d2mp=d2mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end









