      double complex function d6pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
C     DENOMS (ug),(bc),(neug)
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s245
      double complex t2
      s245=s(j2,j4)+s(j2,j5)+s(j4,j5)
      d6pm=-fourrt2*za(j4,j2)*t2(j5,j2,j4,j6)*t2(j7,j1,j3,jb)
     & /(za(j3,jb)*za(j3,j1)*s(j4,j5)*s245)
      return
      end

      double complex function d6mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s245
      s245=s(j2,j4)+s(j2,j5)+s(j4,j5)
      d6mm=fourrt2*za(j4,j2)*t2(j5,j2,j4,j6)*zb(j7,j1)
     & *zb(jb,j1)/(zb(j3,jb)*zb(j3,j1)*s(j4,j5)*s245)
      return
      end



      double complex function d6pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d6pm
      d6pp=d6pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d6mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d6mm
      d6mp=d6mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end



