      subroutine gggghn(p1,p2,p3,p4,p,n,c1234,c1243,c1423)
      implicit none
C     This is the reduced matrix element squared
C     for the process
c     g(p1)+g(p2) --> H((p5+p6)+g(p3)+g(p4)
c     calculated by the program gggghnn.frm
c     with p1 contracted with the vector n
c     Returned are: c1234, c1243, c1324 which are the 3 colour
c     orderings necessary for the subtraction. The sum of these
c     is the full matrix element
      include 'constants.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4
      double precision p(mxpart,4),n(4)
      double precision c1234,c1243,c1423
      call dotem(6,p,s)
      call gggghn1(p1,p2,p3,p4,p,n,c1234)
      call gggghn1(p1,p2,p4,p3,p,n,c1243)
      call gggghn1(p1,p4,p2,p3,p,n,c1423)
      c1234=0.5d0*V*xnsq*c1234
      c1243=0.5d0*V*xnsq*c1243
      c1423=0.5d0*V*xnsq*c1423

c--- This is the total (no longer returned)
c      gggghn=0.5d0*V*xnsq*(c1234+c1243+c1324)

      return
      end
