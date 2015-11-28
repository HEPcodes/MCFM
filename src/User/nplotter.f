      subroutine bookplot(n,tag,titlex,var,wt,xmin,xmax,dx,llplot) 
      implicit none
      include 'nplot.f'
      integer n
      character*(*) titlex
      character llplot*3,tag*4
      double precision var,wt,xmin,xmax,dx
      logical creatent,dswhisto
      common/outputflags/creatent,dswhisto

      if (tag.eq.'book') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call mbook(n,titlex,dx,xmin,xmax)
        else
c--- DSW histograms - call hbook booking routine
        call dswhbook(n,titlex,dx,xmin,xmax)
        endif
      elseif (tag .eq. 'plot') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call mfill(n,var,wt)
        else
c--- DSW histograms - call hbook filling routine
          call dswhfill(n,var,wt)
        endif
        linlog(n)=llplot
        titlearray(n)=titlex
      endif
      return
      end

      subroutine ebookplot(n,tag,var,wt) 
      implicit none
      include 'PDFerrors.f'
      integer n
      double precision var,wt
      character tag*4
      logical creatent,dswhisto
      common/outputflags/creatent,dswhisto

      if (PDFerrors .eqv. .false.) return

      if (tag.eq.'book') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call ebook(n)
        endif
      elseif (tag .eq. 'plot') then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call efill(n,var,wt)
        endif
      endif

      return
      end

      subroutine nplotter(p,wt,switch)
      implicit none
      include 'bbproc.f'
      include 'clustering.f'
      include 'constants.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'process.f'
      integer n,switch,i5,i6,i7,nu
      character tag*4
      double precision DETABB
     & ,DPHIBB
     & ,ETAB1
     & ,ETAB2
     & ,ETANOB
     & ,ETARAP
     & ,ETARAPTWO
     & ,ETBIN
     & ,ETDOUBLEBIN
     & ,GETET
     & ,HT_STOP
     & ,M56CLUST
     & ,MBB
     & ,MERECON
     & ,MLBNU
     & ,PT
     & ,PTB1
     & ,PTB2
     & ,PTNOB
     & ,PTQ1
     & ,ETAQ1
     & ,PTOTHER
     & ,ETAOTHER
     & ,PTTWO
     & ,QETA
     & ,R
     & ,R57
     & ,R67
     & ,RBB
     & ,RECONCORR
     & ,SWAP
     & ,WT
     & ,YRAPTWO


      double precision m34,etmiss,misset,
     . p(mxpart,4),HT
      double precision eta3,eta4,eta5,eta6,eta7,eta8,eta34,eta56,y34
      double precision r34,r35,r45,r36,r46,r56
      double precision pt3,pt4,pt5,pt6,pt7,pt8,pt34,pt56
      double precision bclustmass,etvec(4)
      integer nproc,eventpart,ib1,ib2,nqcdjets,nqcdstart
      logical first,jetmerge
      logical creatent,dswhisto
      character*30 runstring
      common/runstring/runstring
      common/outputflags/creatent,dswhisto
      common/nproc/nproc
      common/nqcdjets/nqcdjets,nqcdstart
      common/jetmerge/jetmerge
      common/stopvars/ht_stop,qeta,mlbnu,merecon,reconcorr
      data first/.true./
      save first
       if (first) then
        first=.false.
        tag='book'
c--- ensure we initialize all possible histograms
        eventpart=npart+3
        eta3=0d0
        pt3=0d0
        eta4=0d0
        pt4=0d0
        eta5=0d0
        pt5=0d0
        eta6=0d0
        pt6=0d0
        eta7=0d0
        pt7=0d0
        eta8=0d0
        pt8=0d0
        eta34=0d0
        y34=0d0
        pt34=0d0
        r34=0d0
        r35=0d0
        r45=0d0
        eta56=0d0
        pt56=0d0        
        m56clust=0d0
        r36=0d0
        r46=0d0
        r56=0d0
        r57=0d0
        r67=0d0
        misset=0d0
        etbin=0d0
        mbb=0d0
        etab1=0d0
        etab2=0d0
        etanob=0d0
        ptb1=0d0
        ptb2=0d0
        ptQ1=0d0
        etaQ1=0d0
        ptother=0d0
        etaother=0d0
        ptnob=0d0
        rbb=0d0
        detabb=0d0
        dphibb=0d0
        HT=0d0
        jetmerge=.true.
        jets=nqcdjets
        goto 99
      else
        tag='plot'
      endif

c--- 'eventpart' will contain the number of actual particles that have
c--- a defined momentum. For most processes, this is calculated as follows:
c----  for lowest order and virtual terms switch=0 and eventpart=npart+2
c---   for real events switch=0 and eventpart=npart+2
c---   for real counter-events switch=1 and eventpart=npart+1
c--- There are some processes for which this is not correct and these
c---  are handled with reference to nproc  
      eventpart=npart-switch+2
      if (jets .gt. 0) eventpart=4+jets
      if (    (case .eq. 'WWqqbr') .or. (case .eq. 'WWnpol')
     .    .or.(case .eq. 'WZbbar') .or. (case .eq. 'ZZlept') ) then
        eventpart=6+jets
      elseif ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh')) then
        eventpart=6+jets
      endif
      if (nproc .eq. 73) eventpart=4+jets
      
c--- re-order jets according to pt, for a W/Z+jet event        
      if ((case .eq. 'W_1jet') .or. (case .eq. 'Z_1jet') .or.
     .    (case .eq. 'W_2jet') .or. (case .eq. 'Z_2jet') .or.
     .    (case .eq. 'W_3jet') .or. (case .eq. 'Z_3jet')) then
        if (algorithm .eq. 'cone') then
          if (jets .gt. 0) pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))
          if (jets .gt. 1) pt6=getet(p(6,4),p(6,1),p(6,2),p(6,3))
          if (jets .gt. 2) pt7=getet(p(7,4),p(7,1),p(7,2),p(7,3))
        else
          if (jets .gt. 0) pt5=pt(5,p)
          if (jets .gt. 1) pt6=pt(6,p)
          if (jets .gt. 2) pt7=pt(7,p)        
        endif
        i5=5
        i6=6
        i7=7
c--- sort for 2 jets 
        if (jets .eq. 2) then          
          if (pt6 .gt. pt5) then
            i5=6
            i6=5
          endif
        endif
c--- sort for 3 jets 
        if (jets .eq. 3) then
          if ((pt5 .gt. pt6) .and. (pt5 .gt. pt7)) then
             i5=5
            if (pt6 .gt. pt7) then
              i6=6
              i7=7
            else
              i6=7
              i7=6
            endif
          endif
          if ((pt6 .gt. pt5) .and. (pt6 .gt. pt7)) then
             i5=6
            if (pt5 .gt. pt7) then
              i6=5
              i7=7
            else
              i6=7
              i7=5
            endif
          endif
          if ((pt7 .gt. pt5) .and. (pt7 .gt. pt6)) then
             i5=7
            if (pt5 .gt. pt6) then
              i6=5
              i7=6
            else
              i6=6
              i7=5
            endif
          endif
        endif
c--- perform exchange
        do nu=1,4
          p(5,nu)=p(i5,nu)
          p(6,nu)=p(i6,nu)
          p(7,nu)=p(i7,nu)
        enddo
      endif      

      m34=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)

      if (bbproc .and. clustering) then
c--- returns zero cluster mass if two b's are in one jet
        m56clust=dsqrt(bclustmass(p))
      else
        m56clust=dsqrt((p(5,4)+p(6,4))**2
     .                -(p(5,1)+p(6,1))**2
     .                -(p(5,2)+p(6,2))**2
     .                -(p(5,3)+p(6,3))**2)
      endif  

      eta3=etarap(3,p)
      pt3=pt(3,p)
      eta4=etarap(4,p)
      pt4=pt(4,p)        
      eta34=etaraptwo(3,4,p)
      y34=yraptwo(3,4,p)
      pt34=pttwo(3,4,p)
      r34=R(p,3,4)
      HT=pt3+pt4
              
      if (eventpart .gt. 4) then        
      eta5=etarap(5,p)
      pt5=pt(5,p)
      r35=R(p,3,5)
      r45=R(p,4,5)
      HT=HT+pt5
 

      endif
      
      if (eventpart .gt. 5) then        
      eta6=etarap(6,p)
      pt6=pt(6,p)
      eta56=etaraptwo(5,6,p)
      eta56=etaraptwo(5,6,p)
      pt56=pttwo(5,6,p)
      r36=R(p,3,6)
      r46=R(p,4,6)
      r56=R(p,5,6)
      HT=HT+pt6
      endif

      if (eventpart .gt. 6) then        
      eta7=etarap(7,p)
      pt7=pt(7,p)
      r57=R(p,5,7)
      r67=R(p,6,7)
      HT=HT+pt7
      endif

      if (eventpart .gt. 7) then        
      eta8=etarap(8,p)
      pt8=pt(8,p)
      HT=HT+pt8
      endif
       
      misset=etmiss(p,etvec)

c--- set-up variables to catch b's
        if (bbproc) then
          if     (jets .eq. 1) then
            write(6,*) 'Error: bbproc set, but only 1 jet in nplotter.f'
            stop
          elseif (jets .eq. 2) then
            mbb=m56clust
            ptb1=pt5
            ptb2=pt6
            etab1=eta5
            etab2=eta6
            rbb=r56
            if (ptb2 .gt. ptb1) then
              swap=ptb1
              ptb1=ptb2
              ptb2=swap
              swap=etab1
              etab1=etab2
              etab2=swap
            endif
          elseif (jets .eq. 3) then
            call getbs(p,ib1,ib2)
            if     (ib1 .eq. 5) then
              ptb1=pt5
              etab1=eta5
            elseif (ib1 .eq. 6) then
              ptb1=pt6
              etab1=eta6
            elseif (ib1 .eq. 7) then
              ptb1=pt7
              etab1=eta7
            endif
            if     (ib2 .eq. 5) then
              ptb2=pt5
              etab2=eta5
            elseif (ib2 .eq. 6) then
              ptb2=pt6
              etab2=eta6
            elseif (ib2 .eq. 7) then
              ptb2=pt7
              etab2=eta7
            endif
            if (ptb2 .gt. ptb1) then
              swap=ptb1
              ptb1=ptb2
              ptb2=swap
              swap=etab1
              etab1=etab2
              etab2=swap
            endif
            if     (ib1+ib2 .eq. 11) then
              mbb=dsqrt((p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     .                 -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2)
              ptnob=pt7
              etanob=eta7
            elseif (ib1+ib2 .eq. 12) then
              mbb=dsqrt((p(5,4)+p(7,4))**2-(p(5,1)+p(7,1))**2
     .                 -(p(5,2)+p(7,2))**2-(p(5,3)+p(7,3))**2)
              ptnob=pt6
              etanob=eta6
            elseif (ib1+ib2 .eq. 13) then
              mbb=dsqrt((p(6,4)+p(7,4))**2-(p(6,1)+p(7,1))**2
     .                 -(p(6,2)+p(7,2))**2-(p(6,3)+p(7,3))**2)
              ptnob=pt5
              etanob=eta5
            endif
            rbb=r(p,ib1,ib2)
          else
            if     (nproc .eq. 152) then
              call getbs(p,ib1,ib2)
              ptb1=pt(ib1,p)
              ptb2=pt(ib2,p)
              etab1=etarap(ib1,p)
              etab2=etarap(ib2,p)
              rbb=r(p,ib1,ib2)
            elseif (nproc .eq. 103) then
c--- set b-quarks to 5 and 6 arbitrarily: if you're interested in
c--- this process, these should be set up properly 
              ib1=5
              ib2=6
              ptb1=pt(ib1,p)
              ptb2=pt(ib2,p)
              etab1=etarap(ib1,p)
              etab2=etarap(ib2,p)
              rbb=r(p,ib1,ib2)
            else
            write(6,*) 'Unforeseen # of jets and b-quarks in nplotter.f'
            stop
            endif
          endif
          detabb=dabs(etab1-etab2)
          dphibb=dsqrt(dabs(rbb**2-detabb**2))
        endif

      if     (case .eq. 'gQ__ZQ') then
c--- for the Z+b process, catch the highest pt heavy quark
        call getptQ1(pt5,pt6,eta5,eta6,ptQ1,etaQ1,1)
      elseif (case .eq. 'H_1jet') then
c--- for the H+b process, catch the most central heavy quark
        call getptQ1(pt5,pt6,eta5,eta6,ptQ1,etaQ1,2)
        if (jets .gt. 1) then
          ptother=pt5+pt6-ptQ1
          etaother=eta5+eta6-etaQ1
        else
          ptother=-1d0
          etaother=99d0
        endif
      else
        ptQ1=-1d0
        etaQ1=99d0
        ptother=-1d0
        etaother=99d0
      endif

      if (nproc .eq. 61) then
        etbin=etdoublebin(pt4,pt5)
      endif

   99 continue

c--- Book and fill ntuple if that option is set
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt)  
        return    
      endif

c--- Otherwise, fill the histograms 
      n=1                  
      
      if (runstring(1:4) .eq. 'stop') then
c--- Special histograms for single top search
        if (tag .eq. 'plot') then
c        call stopcuts(p,eventpart,ht,qeta,mlbnu,merecon,reconcorr)
c these variables should already be in the common-block
        else
          ht_stop=-1d0
          qeta=99d0
          mlbnu=-1d0
          merecon=-1d0
          reconcorr=-99d0
        endif
      call bookplot(n,tag,'    HT',HT_stop,wt,100d0,500d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'  qeta',qeta,wt,-2.8d0,3.2d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,' mlbnu',mlbnu,wt,0d0,300d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'reconcorr?',reconcorr,wt,-1d0,1d0,1d0,'lin')
      n=n+1
      call bookplot(n,tag,'MEn recon',merecon,wt,0d0,300d0,10d0,'lin')
      n=n+1

      endif
      
      call bookplot(n,tag,'      HT',HT,wt,0d0,500d0,10d0,'lin')
      n=n+1
c --- Histograms to monitor exclusive/inclusive cross-sections:
      if ((jets .eq. nqcdjets)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'= #LO j ',0.5d0,wt,0d0,1d0,1d0,'lin')
      endif
      n=n+1
      if ((jets .ge. nqcdjets)  .or. (tag .eq. 'book')) then
      call bookplot(n,tag,'>= #LO j',0.5d0,wt,0d0,1d0,1d0,'lin')
      endif
      n=n+1
      call bookplot(n,tag,'    eta3',eta3,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt3',pt3,wt,0d0,150d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'    eta4',eta4,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt4',pt4,wt,0d0,150d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'   eta34',eta34,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'   eta34',eta34,wt,-1d0,1d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'     y34',y34,wt,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'    pt34',pt34,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'    pt34',pt34,wt,0d0,200d0,10d0,'log')
      call ebookplot(n,tag,pt34,wt)
      n=n+1
      call bookplot(n,tag,'     m34',m34,wt,10d0,200d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r34',r34,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'  misset',misset,wt,0d0,100d0,2d0,'lin')
      n=n+1
      if (eventpart .gt. 4) then
      call bookplot(n,tag,'    eta5',eta5,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt5',pt5,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'     r35',r35,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r45',r45,wt,0d0,4d0,0.1d0,'lin')
      n=n+1

c--- The selected heavy quark jet for Z+b and H+b
      call bookplot(n,tag,' ptQ_central',ptQ1,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'etaQ_central',etaQ1,wt,-5d0,5d0,0.4d0,'lin')
      n=n+1
c--- The other jet for H+b (which is a heavy quark for process 143)
      call bookplot(n,tag,' pt other',ptother,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'eta other',etaother,wt,-5d0,5d0,0.4d0,'lin')
      n=n+1

      endif
      
      if (eventpart .gt. 5) then
      call bookplot(n,tag,'    eta6',eta6,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt6',pt6,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'   eta56',eta56,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'    pt56',pt56,wt,10d0,150d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'     r36',r36,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r46',r46,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r56',r56,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,0d0,200d0,5d0,'log')
      n=n+1
      endif

      if (bbproc) then
c--- Leading b Et, 5 GeV bins from 15 to 200
      call bookplot(n,tag,'bjet1 Et',ptb1,wt,15d0,200d0,5d0,'log')
      n=n+1
c--- Second b Et, 5 GeV bins from 15 to 200
      call bookplot(n,tag,'bjet2 Et',ptb2,wt,15d0,200d0,5d0,'log')
      n=n+1
c--- Dijet invariant mass, 2 b-jets, 10 GeV bins from 0 to 250
      call bookplot(n,tag,' bb mass',mbb,wt,0d0,350d0,10d0,'log')
      n=n+1
c--- Delta_R(b,b), 2 b-jets, bins of 0.2 from 0.35 to 4.95
      call bookplot(n,tag,'deltaRbb',rbb,wt,0.35d0,4.95d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'   etab1',etab1,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'    ptb1',ptb1,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'   etab2',etab2,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'    ptb2',ptb2,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'     mbb',mbb,wt,0d0,200d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'     mbb',mbb,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'     rbb',rbb,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'  detabb',detabb,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'  dphibb',dphibb,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      endif

      if (nproc .eq. 61) then
      call bookplot(n,tag,'   etbin',etbin,wt,0.5d0,25.5d0,1d0,'lin')
      n=n+1
      endif
      
      if (eventpart .gt. 6) then

      if (bbproc) then
c--- Non b-jet Et, 5 GeV bins from 15 to 200
      call bookplot(n,tag,'non-b Et',ptnob,wt,15d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'  etanob',etanob,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'   ptnob',ptnob,wt,0d0,200d0,5d0,'log')
      n=n+1
      endif

      call bookplot(n,tag,'    eta7',eta7,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt7',pt7,wt,0d0,100d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r57',r57,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r67',r67,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      endif      

      if (eventpart .gt. 7) then
      call bookplot(n,tag,'    eta8',eta8,wt,-5d0,5d0,0.5d0,'lin')
      n=n+1
      call bookplot(n,tag,'     pt8',pt8,wt,0d0,100d0,5d0,'lin')
      n=n+1
      endif      

      n=n-1

      if (n .gt. 100) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > 100, which is the built-in maximum'
        stop
      endif
      
      return 
      end
      
      subroutine cross(p,i,j,r)
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),r(3)
      
      r(1)=p(i,2)*p(j,3)-p(j,2)*p(i,3)
      r(2)=p(i,3)*p(j,1)-p(j,3)*p(i,1)
      r(3)=p(i,1)*p(j,2)-p(j,1)*p(i,2)
      
      return
      end
            
      double precision function fphi(n1,n2,p)
      implicit none
      include 'constants.f'
      integer n1,n2
      double precision p(mxpart,4)
    
      fphi=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)
      fphi=fphi/dsqrt(p(n1,1)**2+p(n1,2)**2)
      fphi=fphi/dsqrt(p(n2,1)**2+p(n2,2)**2)
      if (fphi .gt. 1d0) then
        fphi=0d0
      elseif (fphi .lt. -1d0) then
        fphi=pi
      else
        fphi=dacos(fphi)
      endif

      return
      end
                     
      double precision function coslpairet(n1,n2,nm1,nm2,p)
      implicit none
      include 'constants.f'
      integer n1,n2,nm1,nm2,i
      double precision p(mxpart,4),misset(4),pp(4)
      
      if (nm2 .eq. 0) then
        do i=1,4
          misset(i)=p(nm1,i)
        enddo
      else
        do i=1,4
          misset(i)=p(nm1,i)+p(nm2,i)
        enddo
      endif
      
      do i=1,4
        pp(i)=p(n1,i)+p(n2,i)
      enddo
      
      coslpairet=pp(1)*misset(1)+pp(2)*misset(2)
      coslpairet=coslpairet/dsqrt(pp(1)**2+pp(2)**2)
      coslpairet=coslpairet/dsqrt(misset(1)**2+misset(2)**2)
      
      return
      end
            
      double precision function deltar(i,j,p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),phi1,phi2,etarap,dphi
      integer i,j
      
      phi1=atan2(p(i,1),p(i,2))
      phi2=atan2(p(j,1),p(j,2))
      dphi=phi1-phi2
      if (dphi .gt. pi) dphi=twopi-dphi
      if (dphi .lt. -pi) dphi=twopi+dphi
      deltar=(etarap(i,p)-etarap(j,p))**2+dphi**2
      deltar=dsqrt(deltar)
      
      return
      end
      
