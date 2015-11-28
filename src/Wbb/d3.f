      double complex function d3pm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double precision s267,s2367
      double complex t2
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      s2367=s(j2,j3)+s(j3,j6)+s(j3,j7)+s267
      d3pm=fourrt2*za(j6,j2)*t2(j7,j2,j6,jb)*t2(j3,j1,j5,j4)
     & *zb(j5,j1)/(za(j3,jb)*s(j4,j5)*s267*s2367)
      return
      end

      double complex function d3mm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      include 'prods.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex t2
      double precision s267,s2367
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      s2367=s(j2,j3)+s(j3,j6)+s(j3,j7)+s267
      d3mm=-fourrt2*za(j6,j2)*t2(j7,j2,j6,j3)*t2(jb,j1,j5,j4)
     & *zb(j5,j1)/(zb(j3,jb)*s(j4,j5)*s267*s2367)
      return
      end
      
      double complex function d3pp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d3pm
      d3pp=d3pm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end

      double complex function d3mp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      include 'constants.f'
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d3mm
      d3mp=d3mm(j1,j2,j3,j5,j4,j6,j7,jb)
      return
      end


