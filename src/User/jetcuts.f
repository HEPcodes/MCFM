
*********************************************************
*   Switch to routines suitable for various jet cuts    *
*********************************************************
      logical function madejetcuts(p,pjet,njets)
      implicit none
      include 'constants.f'
      logical cutsW_3j,cutsW_2j,cutsW_1j
      logical cutsGGK0,cutsPM_Z
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
        madejetcuts=cutsW_1j(p,pjet,njets)
      elseif (nqcdjets .eq. 2) then
        if     (case .eq. 'W_2jet') then
          madejetcuts=cutsW_2j(p,pjet,njets)
        elseif (case .eq. 'Z_2jet') then
          madejetcuts=cutsW_2j(p,pjet,njets)
        endif
      elseif (nqcdjets .eq. 3) then
          madejetcuts=cutsW_3j(p,pjet,njets)
      endif
      
      return
      end
      
**************************************************************
c--- Example of implementing cuts in W + 3 jets
**************************************************************
      logical function cutsW_3j(p,pjet,njets)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      character*2 plabel(mxpart)
      double precision p(mxpart,4),pjet(mxpart,4),etvec(4)
      double precision ptjet,ptmin,ymax,yrap,pt,pttwo,etmiss,misset
      double precision r,rcut
      integer njets,j
c--- these parameters define the jet observability
      logical first
      common/plabel/plabel
      common/rcut/rcut
      parameter (ptmin=15d0,ymax=2d0)
      data first/.true./
      save first
       
      if (first) then
        first=.false.
      write(6,*)
      write(6,*) '********* Special cuts for W+3 jet process *********'
      write(6,*) '*                                                  *'
      write(6,*) '*    pt(lepton)      >   20 GeV                    *'
      write(6,*) '*   |rap(lepton)|    <    1                        *'
      write(6,*) '*    pt(missing)     >   20 GeV                    *'
      write(6,*) '****************************************************'
      endif

c--- jet pt and rapidity cuts
      
      if(((ptjet(7,p,pjet). lt. ptmin).or.(abs(yrap(7,pjet)) .gt. ymax))
     .    .and. (njets .gt. 2)) njets=njets-1   
      
      if(((ptjet(6,p,pjet). lt. ptmin).or.(abs(yrap(6,pjet)) .gt. ymax))
     .    .and. (njets .gt. 1)) njets=njets-1   
      
      if(((ptjet(5,p,pjet). lt. ptmin).or.(abs(yrap(5,pjet)) .gt. ymax))
     .    .and. (njets .gt. 0)) njets=njets-1   
      
      if (njets .ne. 3) then
        cutsW_3j=.false.
        return
      else
        cutsW_3j=.true.
      endif

c--- special cut for comparing with VECBOS
c--- only applied if CLUSTERING is FALSE and only applicable at LO
      if (clustering .eqv. .false.) then
        if ((r(pjet,5,6) .lt. rcut) .or. (r(pjet,5,7) .lt. rcut)
     . .or. (r(pjet,6,7) .lt. rcut)) then
          cutsW_3j=.false.
          return
        endif
      endif
      
c--- electron and missing-Et cuts

      do j=3,mxpart
        if ((plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')) then
          if ((pt(j,p). lt. 20d0) .or.
     .    (abs(yrap(j,p)) .gt. 1d0)) then
            cutsW_3j=.false.
           endif
        endif
      enddo
      
      misset=etmiss(p,etvec)
      if ((misset. lt. 20d0) .and. (misset .ne. 0d0)) cutsW_3j=.false.
      
      return
      end
      
*********************************************************
c--- Example of implementing cuts in W + 2 jets
*********************************************************

      logical function cutsW_2j(p,pjet,njets)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      character*2 plabel(mxpart)
      double precision p(mxpart,4),pjet(mxpart,4),etvec(4)
      double precision ptjet,ptmin,ymax,yrap,pt,pttwo,etmiss,misset
      double precision r,rcut
      integer njets,j
c--- these parameters define the jet observability
      logical first
      data first/.true./
      save first
      common/plabel/plabel
      common/rcut/rcut
      parameter (ptmin=15d0,ymax=2d0)
       
      if (first) then
        first=.false.
      write(6,*)
      write(6,*) '******** Special cuts for W/Z+2 jet process ********'
      write(6,*) '*                                                  *'
      write(6,*) '*    pt(lepton)      >   20 GeV                    *'
      write(6,*) '*   |rap(lepton)|    <    1                        *'
      write(6,*) '*    pt(missing)     >   20 GeV                    *'
      write(6,*) '****************************************************'
      endif

c--- jet pt and rapidity cuts
      
      if(((ptjet(7,p,pjet). lt. ptmin).or.(abs(yrap(7,pjet)) .gt. ymax))
     .    .and. (njets .gt. 2)) njets=njets-1   
      
      if(((ptjet(6,p,pjet). lt. ptmin).or.(abs(yrap(6,pjet)) .gt. ymax))
     .    .and. (njets .gt. 1)) njets=njets-1   
      
      if(((ptjet(5,p,pjet). lt. ptmin).or.(abs(yrap(5,pjet)) .gt. ymax))
     .    .and. (njets .gt. 0)) njets=njets-1   
      
      if (njets .ne. 2) then
        cutsW_2j=.false.
        return
      else
        cutsW_2j=.true.
      endif

c--- special cut for comparing with VECBOS
c--- only applied if CLUSTERING is FALSE and only applicable at LO
      if (clustering .eqv. .false.) then
        if (r(pjet,5,6) .lt. rcut) then
          cutsW_2j=.false.
          return
        endif
      endif
      
c--- electron and missing-Et cuts

      do j=3,mxpart
        if ((plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')) then
          if ((pt(j,p). lt. 20d0) .or.
     .    (abs(yrap(j,p)) .gt. 1d0)) then
            cutsW_2j=.false.
           endif
        endif
      enddo
      
      misset=etmiss(p,etvec)
      if ((misset. lt. 20d0) .and. (misset .ne. 0d0)) cutsW_2j=.false.
      
      return
      end
      
**************************************************************
c--- Example of implementing cuts in W/Z + 1 jet
**************************************************************
      logical function cutsW_1j(p,pjet,njets)
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
c      if (clustering .eqv. .false.) then
c        if (first) then
c          first=.false.
c      write(6,*)
c      write(6,*) '*** Cut for W/Z+1 jet process with no clustering ***'
c      write(6,*) '*                                                  *'
c      write(6,*) '*   cf. Arnold and Reno, PRD 1989                  *'
c      write(6,*) '*                                                  *'
c      write(6,*) '*      pt(W)   >   20 GeV                          *'
c      write(6,*) '****************************************************'
c        endif
c        if (pttwo(3,4,pjet) .lt. 20d0)  then
c          cutsW_1j=.false.
c        else
c          cutsW_1j=.true.
c        endif
c        return
c      endif
       
c--- Otherwise, implements cuts in W+1 jet a la Giele, Glover, Kosower
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
     .    .and. (njets .gt. 1)) njets=njets-1   
      
      if(((ptjet(5,p,pjet). lt. ptmin).or.(abs(yrap(5,pjet)) .gt. ymax))
     .    .and. (njets .gt. 0)) njets=njets-1   

      cutsW_1j=.true.

c--- if inclusive, njets=1 or 2 is ok      
      if (inclusive) then      
        if (njets .lt. 1) then
          cutsW_1j=.false.
          return
        endif
c--- if exclusive, we need only njets=1
      else
        if (njets .ne. 1) then
          cutsW_1j=.false.
          return
        endif
      endif

c--- electron and missing-Et cuts

      do j=3,mxpart
        if ((plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')) then
          if ((pt(j,p). lt. 20d0) .or.
     .    (abs(yrap(j,p)) .gt. 1d0)) then
            cutsW_1j=.false.
           endif
        endif
      enddo
      
      misset=etmiss(p,etvec)
      if ((misset. lt. 20d0) .and. (misset .ne. 0d0)) cutsW_1j=.false.
      
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
      
