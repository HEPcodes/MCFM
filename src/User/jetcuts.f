
*********************************************************
*   Switch to routines suitable for various jet cuts    *
*********************************************************
      logical function madejetcuts(p,pjet,njets)
      implicit none
      include 'constants.f'
      logical cutsGGK0,cutsGGK1,cutsPM_W,cutsPM_Z
      double precision p(mxpart,4),pjet(mxpart,4)
      integer njets,nqcdjets,nqcdstart
      character*6 case
      common/nqcdjets/nqcdjets,nqcdstart
      common/process/case
                 
c---- we branch on the number of jets required for the process (nqcdjets), 
c---- not the number of jets actually clustered so far (njets)
      if     (nqcdjets .eq. 0) then
        madejetcuts=cutsGGK0(p,pjet,njets)
      elseif (nqcdjets .eq. 1) then
        madejetcuts=cutsGGK1(p,pjet,njets)
      elseif (nqcdjets .eq. 2) then
        if     (case .eq. 'W_2jet') then
          madejetcuts=cutsPM_W(p,pjet,njets)
        elseif (case .eq. 'Z_2jet') then
          madejetcuts=cutsPM_Z(p,pjet,njets)
        endif
      endif
      
      return
      end
      
*********************************************************
c--- Implements cuts in Z + 2 jets a la Parke and Mangano
*********************************************************
      logical function cutsPM_Z(p,pjet,njets)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4)
      double precision ptmin,ymax,ptjet,yrap,pt,deltarj
      integer njets,i
      logical passed3,passed4
      parameter (ptmin=15d0,ymax=2d0)
      character*6 case
      common/process/case

c--- jet pt and rapidity cuts
      
      if (njets .ne. 2) then
        cutsPM_Z=.false.
        return
      else
        cutsPM_Z=.true.
      endif

      do i=5,4+njets
        if (ptjet(i,p,pjet) .lt. ptmin .or.
     .      abs(yrap(i,pjet)) .gt. ymax) then
          cutsPM_Z=.false.
        endif
      enddo

c--- electron and positron jet isolation

      do i=5,4+njets
        if (deltarj(3,i,p,pjet) .lt. 0.4d0) then 
          cutsPM_Z=.false.
        endif
        if (deltarj(4,i,p,pjet) .lt. 0.4d0) then 
          cutsPM_Z=.false.
        endif
      enddo
           
      passed3=.true.     
      passed4=.true.     
            
      if ((pt(3,p). lt. 20d0) .or.
     .    (abs(yrap(3,p)) .gt. 1d0)) then
          passed3=.false.
      endif
      if ((pt(4,p). lt. 20d0) .or.
     .    (abs(yrap(4,p)) .gt. 1d0)) then
          passed4=.false.
      endif
       
      if ((passed3 .eqv. .false.).and.(passed4 .eqv. .false.)) 
     .  cutsPM_Z=.false.
       
      return
      end
      
*********************************************************
c--- Implements cuts in W + 2 jets a la Parke and Mangano
*********************************************************
c--- Returns FALSE if cuts are FAILED
*********************************************************

      logical function cutsPM_W(p,pjet,njets)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),r
      double precision ptmin,ymax,ptjet,yrap,pt
      integer njets,i
      parameter (ptmin=15d0,ymax=2.4d0)

c--- Cuts are implemented for W+ --> nu(p3) + e-(p4) ONLY

c--- jet pt and rapidity cuts
      
      if (njets .ne. 2) then
        cutsPM_W=.false.
        return
      else
        cutsPM_W=.true.
      endif

      do i=5,4+njets
        if ((ptjet(i,p,pjet) .lt. ptmin) .or.
     .      (abs(yrap(i,pjet)) .gt. ymax)) then
          cutsPM_W=.false.
        endif
      enddo

c--- electron and neutrino cuts

      if ((pt(4,p). lt. 20d0) .or.
     .    (abs(yrap(4,p)) .gt. 1d0)) then
          cutsPM_W=.false.
      endif
       
      do i=5,4+njets
        if (r(pjet,4,i) .lt. 0.4d0) then 
          cutsPM_W=.false.
        endif
      enddo

      if (r(pjet,5,6) .lt. 0.7d0) then 
          cutsPM_W=.false.
        endif
            
      if (pt(3,p). lt. 20d0) cutsPM_W=.false.
       
      return
      end
      
**************************************************************
c--- Implements cuts in W + 1 jets a la Giele, Glover, Kosower
**************************************************************
      logical function cutsGGK1(p,pjet,njets)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      character*2 plabel(mxpart)
      double precision p(mxpart,4),pjet(mxpart,4),etvec(4)
      double precision ptjet,ptmin,ymax,yrap,pt,pttwo,etmiss,misset
      integer njets,j
c--- these parameters define the jet observability
      logical first
      data first/.true./
      save first
      common/plabel/plabel
      parameter (ptmin=15d0,ymax=2d0)

c--- special cut for comparing with Arnold and Reno, PRD 1989
c--- only applied if CLUSTERING is FALSE
      if (clustering .eqv. .false.) then
        if (first) then
          first=.false.
      write(6,*)
      write(6,*) '*** Cut for W/Z+1 jet process with no clustering ***'
      write(6,*) '*                                                  *'
      write(6,*) '*   cf. Arnold and Reno, PRD 1989                  *'
      write(6,*) '*                                                  *'
      write(6,*) '*      pt(W)   >   20 GeV                          *'
      write(6,*) '****************************************************'
        endif
        if (pttwo(3,4,pjet) .lt. 20d0)  then
          cutsGGK1=.false.
        else
          cutsGGK1=.true.
        endif
        return
      endif
       
      if (first) then
        first=.false.
      write(6,*)
      write(6,*) '******** Special cuts for W/Z+1 jet process ********'
      write(6,*) '*                                                  *'
      write(6,*) '*   cf. Giele, Glover, Kosower - hep-ph/9302225    *'
      write(6,*) '*                                                  *'
      write(6,*) '*    pt(lepton)      >   20 GeV                    *'
      write(6,*) '*   |rap(lepton)|    <    1                        *'
      write(6,*) '*    pt(missing)     >   20 GeV                    *'
      write(6,*) '****************************************************'
      endif

c--- jet pt and rapidity cuts
      
      if(((ptjet(6,p,pjet). lt. ptmin).or.(abs(yrap(6,pjet)) .gt. ymax))
     .    .and. (njets .eq. 2)) njets=njets-1   
      
      if(((ptjet(5,p,pjet). lt. ptmin).or.(abs(yrap(5,pjet)) .gt. ymax))
     .    .and. (njets .gt. 0)) njets=njets-1   
      
      if (njets .ne. 1) then
        cutsGGK1=.false.
        return
      else
        cutsGGK1=.true.
      endif

c--- electron and missing-Et cuts

      do j=3,mxpart
        if ((plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')) then
          if ((pt(j,p). lt. 20d0) .or.
     .    (abs(yrap(j,p)) .gt. 1d0)) then
            cutsGGK1=.false.
           endif
        endif
      enddo
      
      misset=etmiss(p,etvec)
      if ((misset. lt. 20d0) .and. (misset .ne. 0d0)) cutsGGK1=.false.
      
      return
      end
      
**************************************************************
c--- Implements cuts in W + 0 jets a la Giele, Glover, Kosower
**************************************************************
      logical function cutsGGK0(p,pjet,njets)
      implicit none
      include 'constants.f'
      character*2 plabel(mxpart)
      double precision p(mxpart,4),pjet(mxpart,4),etvec(4)
      double precision ptjet,ptmin,ymax,yrap,pt,etmiss,misset
      integer njets,j
c--- these parameters define the jet observability
      logical first
      data first/.true./
      save first
      common/plabel/plabel
      parameter (ptmin=50d0,ymax=2d0)

      if (first) then
        first=.false.
      write(6,*)
      write(6,*) '******** Special cuts for W/Z+0 jet process ********'
      write(6,*) '*                                                  *'
      write(6,*) '*   cf. Giele, Glover, Kosower - hep-ph/9302225    *'
      write(6,*) '*                                                  *'
      write(6,*) '*    pt(lepton)      >   20 GeV                    *'
      write(6,*) '*   |rap(lepton)|    <    1                        *'
      write(6,*) '*    pt(missing)     >   20 GeV                    *'
      write(6,*) '****************************************************'
      endif

c--- jet pt and rapidity cuts
      if (njets .ne. 0) then
        if ((ptjet(5,p,pjet). lt. ptmin) .or. 
     .      (abs(yrap(5,pjet)) .gt. ymax)) then
          cutsGGK0=.true.    
          njets=0   
        else
          cutsGGK0=.false.
          return
        endif
      else
        cutsGGK0=.true.
      endif

c--- electron and missing-Et cuts

      do j=3,mxpart
        if ((plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')) then
          if ((pt(j,p). lt. 20d0) .or.
     .    (abs(yrap(j,p)) .gt. 1d0)) then
            cutsGGK0=.false.
           endif
        endif
      enddo
       
      misset=etmiss(p,etvec)
      if ((misset. lt. 20d0) .and. (misset .ne. 0d0)) cutsGGK0=.false.
      
      return
      end
      
