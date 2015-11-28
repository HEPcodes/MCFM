      double complex function d4pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
C     DENOMS (bc),(neug),(neu)
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s245,s2345
      double complex t2
      s245=s(j2,j4)+s(j2,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j3,j4)+s(j3,j5)+s245
      d4pm=fourrt2*za(j4,j2)*t2(j5,j2,j4,jb)*t2(j3,j1,j7,j6)
     & *zb(j7,j1)/(za(j3,jb)*s(j4,j5)*s245*s2345)
      return
      end

      double complex function d4mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s245,s2345
      s245=s(j2,j4)+s(j2,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j3,j4)+s(j3,j5)+s245
      d4mm=-fourrt2*za(j4,j2)*t2(j5,j2,j4,j3)*t2(jb,j1,j7,j6)
     & *zb(j7,j1)/(zb(j3,jb)*s(j4,j5)*s245*s2345)
      return
      end
      
      double complex function d4pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d4pm
      d4pp=d4pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d4mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d4mm
      d4mp=d4mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end


