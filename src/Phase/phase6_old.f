      subroutine phase6(r,p1,p2,p6,p7,p4,p5,p8,p9,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
c******* generate phase space for 2-->4 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p4+p5+p6+p7+p8+p9
c---- with all 2 pi's (ie 1/(2*pi)^14)
      logical oldzerowidth
      integer n2,n3
      double precision r(mxdim)
      double precision p1(4),p2(4),p4(4),p5(4),p6(4),p7(4),p8(4),p9(4)
      double precision p12(4),p467(4),p589(4),p89(4),p67(4),smin
      double precision wt,wt0,wt12,wt589,wt467,wt67,wt89
      double precision mass2,width2,mass3,width3
      character*6 case
      common/process/case
      common/breit/n2,n3,mass2,width2,mass3,width3 

      integer j
      parameter(wt0=1d0/twopi**4)
      wt=0d0
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=mb**2
c---- calculate momenta of top and bbbar
      
      n2=1
      n3=1
      if (
     .        (case .eq. 'tt_bbh')
     .   .or. (case .eq. 'tt_bbl')
     .   .or. (case .eq. 'Httbar')
     .   .or. (case .eq. 'vlchk6')) then
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        call phi1_2(r(1),r(2),r(3),r(4),p12,p467,p589,wt12,*99)
      elseif (case .eq. 'tautau') then
        mass2=mtau
        width2=tauwidth
        mass3=mtau
        width3=tauwidth
        oldzerowidth=zerowidth
        zerowidth=.true.
        call phi1_2(r(1),r(2),r(3),r(4),p12,p467,p589,wt12,*99)
        zerowidth=oldzerowidth
      else
        write(*,*) 'Case not supported in phase6.f'
        stop
      endif

c      write(6,*) '467',(p467(4)**2-p467(3)**2-p467(2)**2-p467(1)**2)
c      write(6,*) '589',(p589(4)**2-p589(3)**2-p589(2)**2-p589(1)**2)
            
      mass3=wmass
      width3=wwidth
      if ((case .eq. 'tt_bbh')  .or. (case .eq. 'tt_bbl')
     .     .or. (case .eq. 'vlchk6')) then
        n3=1
        call phi1_2m(mb,r(5),r(6),r(7),smin,p467,p4,p67,wt467,*99)
        call phi1_2m(mb,r(8),r(11),r(12),smin,p589,p5,p89,wt589,*99)
      elseif (case .eq. 'tautau') then
        n3=0
        call phi1_2m(smin,r(5),r(6),r(7),smin,p467,p4,p67,wt467,*99)
        call phi1_2m(smin,r(8),r(11),r(12),smin,p589,p5,p89,wt589,*99)
      endif

c      write(6,*) '67',(p67(4)**2-p67(3)**2-p67(2)**2-p67(1)**2)
c      write(6,*) '89',(p89(4)**2-p89(3)**2-p89(2)**2-p89(1)**2)

      call phi3m0(r(13),r(14),p67,p6,p7,wt67,*99)
      call phi3m0(r(15),r(16),p89,p8,p9,wt89,*99)

      wt=wt0*wt12*wt467*wt589*wt67*wt89
      
      return
 99   wt=0d0
      return 1
      end

