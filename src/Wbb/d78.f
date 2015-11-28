      double complex function d78pm(j1,j2,j3,j4,j5,j6,j7,jb)
c---sum of d7+d8
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2,d8,d7
      double precision s45,s267,s2345,s345
      s45=s(j4,j5)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)

      d7=-t2(j3,j1,j7,j6)*za(j3,j2)*zb(j3,j5)*za(jb,j4)
     &   +t2(j3,j1,j7,j6)*za(jb,j2)*zb(j3,j5)*za(j3,j4)
     &   -t2(j5,j1,j7,j6)*za(j4,j2)
     & *(za(jb,j4)*zb(j3,j4)+za(jb,j5)*zb(j3,j5))

      d8=+t2(j7,j2,j6,j3)*zb(j3,j1)*zb(j3,j5)*za(jb,j4)
     &   -t2(j7,j2,j6,jb)*zb(j3,j1)*zb(j3,j5)*za(j3,j4)
     &   +t2(j7,j2,j6,j4)*zb(j5,j1)
     & *(za(jb,j4)*zb(j3,j4)+za(jb,j5)*zb(j3,j5))
      d7=-fourrt2*zb(j7,j1)/(za(j3,jb)*s45*s345*s2345)*d7
      d8=-fourrt2*za(j6,j2)/(za(j3,jb)*s45*s345*s267)*d8
      d78pm=d7+d8
      return
      end


      double complex function d78mm(j1,j2,j3,j4,j5,j6,j7,jb)
c---sum of d7+d8
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2,d8,d7
      double precision s45,s267,s2345,s345
c      double precision x7sq,x8sq,x78sq
c      common/tempbit/x7sq,x8sq,x78sq
      s45=s(j4,j5)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)
      s2345=s(j2,j3)+s(j2,j4)+s(j2,j5)+s345
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d7=-t2(j3,j1,j7,j6)*za(j3,j2)*zb(jb,j5)*za(j3,j4)
     &   +t2(jb,j1,j7,j6)*za(j3,j2)*zb(j3,j5)*za(j3,j4)
     &   -t2(j5,j1,j7,j6)*za(j4,j2)
     & *(zb(jb,j4)*za(j3,j4)+zb(jb,j5)*za(j3,j5))

      d8=+t2(j7,j2,j6,j3)*zb(j3,j1)*zb(jb,j5)*za(j3,j4)
     &   -t2(j7,j2,j6,j3)*zb(jb,j1)*zb(j3,j5)*za(j3,j4)
     &   +t2(j7,j2,j6,j4)*zb(j5,j1)
     & *(zb(jb,j4)*za(j3,j4)+zb(jb,j5)*za(j3,j5))
      d8=fourrt2*za(j6,j2)/(zb(j3,jb)*s45*s345*s267)*d8
      d7=fourrt2*zb(j7,j1)/(zb(j3,jb)*s45*s345*s2345)*d7
      d78mm=d7+d8
      return
      end


      double complex function d78pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d78pm
      d78pp=d78pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d78mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d78mm
      d78mp=d78mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

