      subroutine phase5_bgd(r,p1,p2,p3,p4,p5,p6,p7,wt)
c----phase space for bgd
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision p123(4),p12(4),p45(4),p67(4),s4567min
      double precision wt,wt123,wt4567,wt45,wt67,wt0

      integer j,n2,n3
      double precision mass2,width2,mass3,width3 
      common/breit/n2,n3,mass2,width2,mass3,width3 

      parameter(wt0=1d0/twopi**3)
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      s4567min=mb**2
      n2=0
      n3=0
      call phi1_2m(zip,r(13),r(12),r(11),s4567min,p12,p3,p123,wt123,*99)
      n2=0
      n3=1
      mass3=wmass
      width3=wwidth
      call phi1_2(r(1),r(2),r(3),r(4),p123,p45,p67,wt4567,*99)
      call phi3m0(r(5),r(6),p45,p4,p5,wt45,*99)
      call phi3m0(r(7),r(8),p67,p6,p7,wt67,*99)
      wt=wt0*wt123*wt4567*wt45*wt67
      return
 99   continue
      wt=0d0
      return
      end

