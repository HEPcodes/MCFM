      double complex function qcdapp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d1pp,d2pp,d3pp,d10pp,d12pp,d78pp
      qcdapp=
     & +d1pp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d2pp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d3pp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d10pp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d12pp(j1,j2,j3,j4,j5,j6,j7,jb)
     & -d78pp(j1,j2,j3,j4,j5,j6,j7,jb)
      return
      end

      double complex function qcdamp(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      double complex d1mp,d2mp,d3mp,d10mp,d12mp,d78mp
      integer j1,j2,j3,j4,j5,j6,j7,jb
      qcdamp=
     & +d1mp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d2mp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d3mp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d10mp(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d12mp(j1,j2,j3,j4,j5,j6,j7,jb)
     & -d78mp(j1,j2,j3,j4,j5,j6,j7,jb)
      return
      end

      double complex function qcdapm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      integer j1,j2,j3,j4,j5,j6,j7,jb
      double complex d1pm,d2pm,d3pm,d10pm,d12pm,d78pm
      qcdapm=
     & +d1pm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d2pm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d3pm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d10pm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d12pm(j1,j2,j3,j4,j5,j6,j7,jb)
     & -d78pm(j1,j2,j3,j4,j5,j6,j7,jb)
      return
      end

      double complex function qcdamm(j1,j2,j3,j4,j5,j6,j7,jb)
      implicit none
      double complex d1mm,d2mm,d3mm,d10mm,d12mm,d78mm
      integer j1,j2,j3,j4,j5,j6,j7,jb
      qcdamm=
     & +d1mm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d2mm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d3mm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d10mm(j1,j2,j3,j4,j5,j6,j7,jb)
     & +d12mm(j1,j2,j3,j4,j5,j6,j7,jb)
     & -d78mm(j1,j2,j3,j4,j5,j6,j7,jb)
      return
      end


