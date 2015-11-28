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
      double precision q(mxpart,4),qfinal(mxpart,4),R
      integer nqcdjets,nqcdstart,isub
      logical first
      character*4 part
      common/part/part
      common/nqcdjets/nqcdjets,nqcdstart
      data first/.true./
      save first
      
      if ((first) .and. ((nqcdjets .gt. 0).or.(part .eq. 'real'))) then
        first=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
      write(6,*)
      write(6,*) '*********** Basic jet-defining parameters **********'
      if (algorithm .eq. 'ktal') then
      write(6,*) '*          (Run II kT clustering algorithm)        *'
      else
      write(6,*) '*              (Run II cone algorithm)             *'
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
      write(6,79) ' * pseudo-cone size, R : ',R
      write(6,*) '*                                                  *'
      if (inclusive) then
      write(6,*) '*        Jet cross-section is INCLUSIVE            *'
      else
      write(6,*) '*        Jet cross-section is EXCLUSIVE            *'
      endif
      write(6,*) '****************************************************'
      call flush(6)
      endif
   79 format(a25,f6.3,'                     *')

      if (algorithm .eq. 'ktal') then
        call genclust_kt(q,R,qfinal,isub)
      else
        call genclust_cone(q,R,qfinal,isub)
      endif

      return
      end
      
