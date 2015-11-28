      double complex function d5pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s267
      double complex t2
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d5pm=fourrt2*za(j6,j2)*t2(j7,j2,j6,j4)*t2(j5,j1,j3,jb)
     & *zb(j3,j1)/(za(j3,jb)*s(j1,j3)*s(j4,j5)*s267)
      return
      end

      double complex function d5mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s267
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      d5mm=fourrt2*za(j6,j2)*t2(j7,j2,j6,j4)*zb(j5,j1)
     & *zb(jb,j1)/(zb(j3,jb)*zb(j3,j1)*s(j4,j5)*s267)
      return
      end


      double complex function d5pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d5pm
      d5pp=d5pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d5mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d5mm
      d5mp=d5mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end


