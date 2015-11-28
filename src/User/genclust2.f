      subroutine genclust2(q,R,qfinal,isub)
c--- this is a wrapper routine for the jet clustering algorithm
c--- either re-route to:
c---  genclust_kt.f     for kt clustering
c---  genclust_cone.f   for cone algorithm
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'jetcuts.f'
      include 'bbproc.f'
      include 'process.f'
      double precision q(mxpart,4),qfinal(mxpart,4),R,Rbbmin
      integer nqcdjets,nqcdstart,isub
      logical first
      character*4 part
      common/part/part
      common/nqcdjets/nqcdjets,nqcdstart
      common/Rbbmin/Rbbmin
      data first/.true./
      save first
      
      if ((first) .and. ((nqcdjets .gt. 0).or.(part .eq. 'real'))) then
        first=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
      write(6,*)
      write(6,*) '*********** Basic jet-defining parameters **********'
      if     (algorithm .eq. 'ktal') then
      write(6,*) '*          (Run II kT clustering algorithm)        *'
      elseif (algorithm .eq. 'cone') then
      write(6,*) '*              (Run II cone algorithm)             *'
      elseif (algorithm .eq. 'hqrk') then
      write(6,*) '*        (Simple cone algorithm for W/Z+Q+j)       *'
      else
      write(6,*)
      write(6,*) 'Invalid selection of algorithm in input file.'
      write(6,*) 'Please select either ktal, cone or hqrk'
      stop
      endif
      write(6,*) '*                                                  *'
      write(6,79) ' *     pt(jet)         > ',ptjetmin
      write(6,79) ' *   |pseudo-rap(jet)| > ',etajetmin   
      write(6,79) ' *   |pseudo-rap(jet)| < ',etajetmax   
      if (bbproc) then
        ptbjetmin=max(ptjetmin,ptbjetmin)
        etabjetmax=min(etajetmax,etabjetmax)
      write(6,79) ' *   pt(b-jet)         > ',ptbjetmin
      write(6,79) ' * |pseudo-rap(b-jet)| < ',etabjetmax   
      endif
      if (algorithm .eq. 'hqrk') then
      write(6,79) ' *   b-bbar separation : ',Rbbmin
      write(6,79) ' *        cone size, R : ',R      
      else
      write(6,79) ' * pseudo-cone size, R : ',R
      endif
      write(6,*) '*                                                  *'
      if ((case .eq. 'W_twdk') .or. (case .eq. 'Wtdkay')) then
      write(6,79) ' *   pt(b-jet @ NLO)   < ',ptbjetmin
      write(6,*) '*                                                  *'
      endif
      if (inclusive) then
      write(6,*) '*        Jet cross-section is INCLUSIVE            *'
      else
      write(6,*) '*        Jet cross-section is EXCLUSIVE            *'
      endif
      write(6,*) '****************************************************'
      call flush(6)
      endif
   79 format(a25,f8.4,'                   *')

      if     (algorithm .eq. 'ktal') then
        call genclust_kt(q,R,qfinal,isub)
      elseif (algorithm .eq. 'cone') then
        call genclust_cone(q,R,qfinal,isub)
      elseif (algorithm .eq. 'hqrk') then
        call genclust_hqrk(q,R,qfinal,isub)
      else
        write(6,*) 'Invalid choice for clustering algorithm: ',algorithm
        stop
      endif

      return
      end
      
