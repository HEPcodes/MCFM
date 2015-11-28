      subroutine gggghn1amp(p1,p2,p3,p4,p,n,square)
      implicit none
C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p5+p6)+g(p3)+g(p4)
c   with momentum 4 contracted with the vector n(mu)
c     calculated by the program gggghn1amp.frm
      include 'constants.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer p1,p2,p3,p4
      double precision ABSq,BAsq,sumsq,square
      double precision p(mxpart,4),n(4),nDn,nDp1,nDp2,nDp3,nDp4
      double precision s123,s124,s234,s134
      double complex vcm(mxpart,mxpart),vcp(mxpart,mxpart)
      double complex za(mxpart,mxpart),zb(mxpart,mxpart)
      double complex n4g1234(2,2,2),n4g1342(2,2,2),n4g1423(2,2,2)
      double complex ABpp,ABpm,ABmp,ABmm,BApp,BApm,BAmp,BAmm
      call ndveccur(6,n,p,vcm)
      call spinoru(6,p,za,zb)
      nDp2=n(4)*p(p2,4)-n(3)*p(p2,3)-n(2)*p(p2,2)-n(1)*p(p2,1)
      nDp3=n(4)*p(p3,4)-n(3)*p(p3,3)-n(2)*p(p3,2)-n(1)*p(p3,1)
      nDp4=n(4)*p(p4,4)-n(3)*p(p4,3)-n(2)*p(p4,2)-n(1)*p(p4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      s123=s(p1,p2)+s(p2,p3)+s(p3,p1)
      s124=s(p1,p2)+s(p2,p4)+s(p4,p1)
      s134=s(p1,p3)+s(p3,p4)+s(p4,p1)
      s234=s(p2,p3)+s(p3,p4)+s(p4,p2)
c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp4).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'gggghn1amp:Error for :',p1,p2,p3,p4
         write(*,*) 'cutoff',1d-3*abs(p(1,4))
         write(6,*) 'nDp4',nDp4
         call flush(6)
         stop
      endif
      n4g1234(2,2,2)= Hgggg1.frm Line 120 ==> gg1234ppp is not an expression
      n4g1234(1,2,2)= Hgggg1.frm Line 121 ==> gg1234mpp is not an expression
      n4g1234(2,1,2)= Hgggg1.frm Line 122 ==> gg1234pmp is not an expression
      n4g1234(2,2,1)= Hgggg1.frm Line 123 ==> gg1234ppm is not an expression
      n4g1234(2,1,1)= Hgggg1.frm Line 124 ==> gg1234pmm is not an expression
      n4g1234(1,2,1)= Hgggg1.frm Line 125 ==> gg1234mpm is not an expression
      n4g1234(1,1,2)= Hgggg1.frm Line 126 ==> gg1234mmp is not an expression
      n4g1234(1,1,1)= Hgggg1.frm Line 127 ==> gg1234mmm is not an expression
C
      n4g1234(2,2,2)= Hgggg1.frm Line 129 ==> gg1234ppp is not an expression
      n4g1234(1,2,2)= Hgggg1.frm Line 130 ==> gg1234mpp is not an expression
      n4g1234(2,1,2)= Hgggg1.frm Line 131 ==> gg1234pmp is not an expression
      n4g1234(2,2,1)= Hgggg1.frm Line 132 ==> gg1234ppm is not an expression
      n4g1234(2,1,1)= Hgggg1.frm Line 133 ==> gg1234pmm is not an expression
      n4g1234(1,2,1)= Hgggg1.frm Line 134 ==> gg1234mpm is not an expression
      n4g1234(1,1,2)= Hgggg1.frm Line 135 ==> gg1234mmp is not an expression
      n4g1234(1,1,1)= Hgggg1.frm Line 136 ==> gg1234mmm is not an expression
C
      n4g1234(2,2,2)= Hgggg1.frm Line 138 ==> gg1234ppp is not an expression
      n4g1234(1,2,2)= Hgggg1.frm Line 139 ==> gg1234mpp is not an expression
      n4g1234(2,1,2)= Hgggg1.frm Line 140 ==> gg1234pmp is not an expression
      n4g1234(2,2,1)= Hgggg1.frm Line 141 ==> gg1234ppm is not an expression
      n4g1234(2,1,1)= Hgggg1.frm Line 142 ==> gg1234pmm is not an expression
      n4g1234(1,2,1)= Hgggg1.frm Line 143 ==> gg1234mpm is not an expression
      n4g1234(1,1,2)= Hgggg1.frm Line 144 ==> gg1234mmp is not an expression
      n4g1234(1,1,1)= Hgggg1.frm Line 145 ==> gg1234mmm is not an expression
      write(6,*) n4g1234(1,1,1)+n4g1342(1,1,1)+n4g1423(1,1,1)
      write(6,*) n4g1234(2,1,1)+n4g1342(2,1,1)+n4g1423(2,1,1)
      write(6,*) n4g1234(1,2,1)+n4g1342(1,2,1)+n4g1423(1,2,1)
      write(6,*) n4g1234(1,1,2)+n4g1342(1,1,2)+n4g1423(1,1,2)
      write(6,*) n4g1234(1,2,2)+n4g1342(1,2,2)+n4g1423(1,2,2)
      write(6,*) n4g1234(2,1,2)+n4g1342(2,1,2)+n4g1423(2,1,2)
      write(6,*) n4g1234(2,2,1)+n4g1342(2,2,1)+n4g1423(2,2,1)
      return
      end
