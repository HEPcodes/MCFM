      double complex function d11pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
C     DENOMS 345,ned,bg
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s267,s345
      double complex t2
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d11pp=fourrt2*za(j6,j2)*t2(j7,j2,j6,j5)
     & *(za(j3,jb)*zb(j3,j1)+za(j4,jb)*zb(j4,j1))
     & /(za(j3,jb)*za(j4,j3)*s345*s267)
      return
      end

      double complex function d11mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s267,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d11mp=-fourrt2*za(j6,j2)*zb(j4,jb)*zb(j4,j1)
     &*t2(j7,j2,j6,j5)/(zb(j3,jb)*zb(j4,j3)*s345*s267)
      return
      end

      double complex function d11pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s267,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d11pm=fourrt2*za(j6,j2)*za(jb,j4)*zb(j5,j1)
     & *t2(j7,j2,j6,j4)/(za(j3,jb)*za(j3,j4)*s345*s267)
      return
      end

      double complex function d11mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s267,s345
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d11mm=-fourrt2*za(j6,j2)*zb(j5,j1)
     & *(zb(jb,j3)*t2(j7,j2,j6,j3)+zb(jb,j4)*t2(j7,j2,j6,j4))
     & /(zb(j3,jb)*zb(j3,j4)*s345*s267)
      return
      end

