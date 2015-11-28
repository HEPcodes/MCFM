      logical function gencuts(pjet,njets)
************************************************************************
*   Author: J.M. Campbell, 5th December 2001                           *
*                                                                      *
*   This routine imposes a generic set of cuts that can be applied     *
*   to all processes in process.DAT, using the parton momenta in pjet  *
*   which have already passed through the jet clustering algorithm     *
*                                                                      *
*   Only a basic set of variables is tested:                           *
*     pt(lepton) > leptpt, eta(lepton) < leptrap, missing Et > misspt  *
*                                                                      *
*   If leptpt2 or leptrap2 is not zero, then any leptons beyond the    *
*   leading-pt one must instead satisfy:                               *
*     pt(lepton) > leptpt2, eta(lepton) < leptrap2                     *
*                                                                      *
*   For processes where one would like to apply jet-like cuts, but     *
*   no clustering has been performed, additional cuts apply:           *
*     pt(jet) > jetpt, eta(jet) < jetrap, R(jet1,jet2) > Rcut          *
*                                                                      *
*   Finally, if further process-specific cuts are necessary,           *
*   an appropriate second routine may be called                        *
*                                                                      *
*   Return TRUE if this point FAILS the cuts                           *
*                                                                      *
************************************************************************
      implicit none
      include 'bbproc.f'
      include 'constants.f'
      include 'jetcuts.f'
      include 'process.f'
      logical first,passedlept
      character*2 plabel(mxpart)
      integer njets,j,k,countb,bindex(mxpart),jindex,kindex,ib1,ib2
      integer countlept,leptindex(mxpart),countgamm,gammindex(mxpart),
     . countjet,jetindex(mxpart),pntr,lbjscheme,maxparts,notag
      double precision pjet(mxpart,4),etvec(4)
      double precision leptpt,leptrap,misspt,jetpt,jetrap,gammpt,gammrap
      double precision pt,etarap,etmiss,evtmisset,R,Rcut,gammcone,
     . gammcut,etaj,etak,etalept,mll
      double precision Rjlmin,Rllmin,delyjjmin,leptpt2,leptrap2
      double precision delta(mxpart),discr,ptjet(mxpart),etabuffer
c    . ,MJJ
      logical jetsopphem,passed,hwwjetcuts
c      integer nu
c      double precision sumjetpt(2)
      double precision ht,qeta,mlbnu,merecon,reconcorr
      double precision dphi_ll,m_ll,mtrans,scut1,scut2
      double precision pt34,pttwo,phill,phillcut,etajet2cut,mllcut
      character*30 runstring
      common/runstring/runstring
      common/stopvars/ht,qeta,mlbnu,merecon,reconcorr
      common/hwwvars/dphi_ll,m_ll,mtrans,scut1,scut2
      common/leptcuts/leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut,
     . lbjscheme,jetsopphem
      common/plabel/plabel
      common/rcut/Rcut
      common/notag/notag
************************************************************************
*     Set-up the jet-like cut parameters here                          *
      parameter (jetpt=15d0,jetrap=2d0)
************************************************************************
      parameter(phillcut=1.2d0,etajet2cut=2.5d0,mllcut=75d0)
      data first/.true./
      
      gencuts=.false.

c--- perform extra H(->WW)+jet search cuts if there is a minimum
c---  jet rapidity, as is usually done
      if (abs(etajetmin) .gt. 1d-6) then
        hwwjetcuts=.true.
      else
        hwwjetcuts=.false.
      endif 

c--- extra transverse mass cut in W+jets for CDF
      if (runstring(1:7) .eq. 'cdfjoey') then
        mtrans=
     .   (pjet(3,1)*pjet(4,1)+pjet(3,2)*pjet(4,2))
     .   /dsqrt((pjet(3,1)**2+pjet(3,2)**2)
     .         *(pjet(4,1)**2+pjet(4,2)**2))
c---    transverse mass calculation
        mtrans=2d0*dsqrt(pjet(3,1)**2+pjet(3,2)**2)
     .   *dsqrt(pjet(4,1)**2+pjet(4,2)**2)*(1d0-mtrans)
        mtrans=dsqrt(max(mtrans,0d0))
        if (mtrans .lt. 20d0) then
	  gencuts=.true.
	  return
	endif       
      endif
      
      if (runstring(1:4) .eq. 'stop') then
c--- do single-top search cuts instead
        maxparts=4+njets
        if ((case .eq. 'tt_bbl') .or. (case .eq. 'tt_bbh'))
     .    maxparts=6+njets
        call stopcuts(pjet,maxparts,ht,qeta,mlbnu,merecon,reconcorr)  
        if (ht .lt. 0d0) gencuts=.true.
        return
      endif
       
      if (runstring(1:3) .eq. 'hww') then
c--- do H->WW search cuts instead
        maxparts=6+njets
        call hwwcuts(pjet,maxparts,dphi_ll,m_ll,mtrans,scut1,scut2)  
        if (mtrans .lt. 0d0) gencuts=.true.
        return
      endif
       
      if (runstring(1:3) .eq. 'wbf') then
c--- do WBF search cuts instead
        maxparts=4+njets
	if (runstring(4:8) .eq. 'jeppe') then
          call wbfcuts_jeppe(pjet,maxparts,passed)  
        else
	  call wbfcuts(pjet,maxparts,passed)  
	endif
        if (passed .eqv. .false.) gencuts=.true.
        return
      endif
       
c--- Look for particles that should be treated as jets,
c--- so far only b decays from Z-bosons and
c--- hadronic decay of the W in diboson processes
c--- If these are present, we will do additional cuts
      countb=0
      do j=3,mxpart
         if ((plabel(j) .eq. 'qb') .or. (plabel(j) .eq. 'ab')
     .  .or. (plabel(j) .eq. 'qq') .or. (plabel(j) .eq. 'qa')) then
           countb=countb+1
           bindex(countb)=j
         endif
      enddo

c--- write-out the cuts we are using
      if (first) then
      first=.false.
      
      write(6,*)
      write(6,*)  '****************** Generic cuts ********************'
      write(6,*)  '*                                                  *'
      write(6,99) '*        pt(lepton)      >   ',leptpt,
     .                ' GeV            *'
      write(6,99) '*      |eta(lepton)|     <   ',leptrap,
     .                '                *'
      write(6,99) '*       pt(missing)      >   ',misspt,
     .                ' GeV            *'
      if ((leptpt2 .ne. 0d0) .or. (leptrap2 .ne. 0d0)) then
      write(6,99) '*     pt(2nd+ lepton)    >   ',leptpt2,
     .                ' GeV            *'
      write(6,99) '*   |eta(2nd+ lepton)|   <   ',leptrap2,
     .                '                *'
      endif
      write(6,99) '*      R(jet,lepton)     >   ',Rjlmin,
     .                '                *'
      write(6,99) '*     R(lepton,lepton)   >   ',Rllmin,
     .                '                *'
      write(6,99) '* |eta(jet1)-eta(jet2)|  >   ',delyjjmin,
     .                '                *'
      if (jetsopphem) then
      write(6,*) '*           eta(jet1) . eta(jet2)  <  0            *' 
      endif
      if     (lbjscheme .eq. 1) then
      write(6,*) '*        eta(jet1)  <  eta(lept)  <  eta(jet2)     *' 
      elseif (lbjscheme .ge. 2) then
      write(6,*) '*  eta(jet1)+Rcut  <  eta(lept)  <  eta(jet2)-Rcut *' 
      lbjscheme=2
      endif
      write(6,99) '*   pt(photon)           >   ',gammpt,
     .                '                *'
      write(6,99) '*   eta(photon)          <   ',gammrap,
     .                '                *'
      write(6,99) '*   pt(hadronic)         <   ',gammcut,
     .                '    pt(photon)  *'
      write(6,99) '*   (in cone around photon of',gammcone,
     .                ')               *'
      if (countb .gt. 0) then
      write(6,*)  '*                                                  *'
      write(6,99) '*      pt(jet)       >   ',jetpt,
     .                ' GeV                *'
      write(6,99) '*    |eta(jet)|      <   ',jetrap,
     .                '                    *'
      write(6,99) '*   R(jet1,jet2)     >   ',Rcut,
     .                '                    *'
      endif
      if (hwwjetcuts) then
      write(6,99) '*   phi(lepton,lepton)   <   ',phillcut,
     .                '                *'
      write(6,99) '*  |eta(2nd jet)|        >   ',etajet2cut,
     .                '                *'
      write(6,99) '*  m(lepton,lepton)      <   ',mllcut,
     .                '                *'
      endif
      write(6,*)  '****************************************************'
      endif

      countlept=0
      countgamm=0
c--- lepton pt and rapidity cuts
      do j=3,mxpart
        if (     (plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')
     .      .or. (plabel(j) .eq. 'ml') .or. (plabel(j) .eq. 'ma')
     .      .or. (plabel(j) .eq. 'tl') .or. (plabel(j) .eq. 'ta')) then
          countlept=countlept+1
          leptindex(countlept)=j
        endif
        if (plabel(j) .eq. 'ga') then
          countgamm=countgamm+1
          gammindex(countgamm)=j
        endif
      enddo
      
C     Basic pt and rapidity cuts for lepton
      if     (countlept .eq. 1) then
          if (     (pt(leptindex(1),pjet) .lt. leptpt) .or.
     .      (abs(etarap(leptindex(1),pjet)) .gt. leptrap)) then
            gencuts=.true.
            return
          endif
      elseif (countlept .gt. 1) then
c--- loop over all the lepton possibilities for lepton 1 (j)
          j=0
  77      continue
          j=j+1
          passedlept=.true.
          if (     (pt(leptindex(j),pjet) .lt. leptpt) .or.
     .      (abs(etarap(leptindex(j),pjet)) .gt. leptrap)) then
            passedlept=.false.
            goto 78
          endif
          do k=1,countlept
              if (k .ne. j) then
                  if (     (pt(leptindex(k),pjet) .lt. leptpt2) .or.
     .              (abs(etarap(leptindex(k),pjet)) .gt. leptrap2)) then
                    passedlept=.false.
                   endif
              endif 
          enddo         
  78      continue
c--- return to beginning if we failed and there are more leptons to try
          if ((passedlept .eqv. .false.).and.(j .lt. countlept)) goto 77
          gencuts=.not.(passedlept)
      endif

C     Basic pt and rapidity cuts for photon
      if (countgamm .gt. 0) then
            do j=1,countgamm
              if (     (pt(gammindex(j),pjet) .lt. gammpt) .or.
     .          (abs(etarap(gammindex(j),pjet)) .gt. gammrap)) then
                gencuts=.true.
                return
              endif
             enddo
      endif

c--- missing energy cut
      evtmisset=etmiss(pjet,etvec)
      if ((evtmisset .lt. misspt) .and. (evtmisset .ne. 0d0)) then
        gencuts=.true.
        return
      endif
      
c--- lepton-lepton separation (if there are 2 or more leptons)
      if ((countlept .gt. 1)) then
        do j=1,countlept
        do k=j+1,countlept
          if (R(pjet,leptindex(j),leptindex(k)) .lt. Rllmin) then
            gencuts=.true.
            return
          endif
        enddo
        enddo
c--- extra cut on phi(lept,lept) for H(->WW)+jet search
	if (hwwjetcuts) then
          phill=
     .       (pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .       +pjet(leptindex(1),2)*pjet(leptindex(2),2))
     .       /dsqrt((pjet(leptindex(1),1)**2+pjet(leptindex(1),2)**2)
     .             *(pjet(leptindex(2),1)**2+pjet(leptindex(2),2)**2))
          if (phill .lt. -0.999999999D0) phill=-1d0
          phill=dacos(phill) 
	  if (phill .gt. phillcut) then
	    gencuts=.true.
	    return
	  endif
	endif
c--- extra cut on m(lept,lept) for H(->WW)+jet search
	if (hwwjetcuts) then
          mll=dsqrt(2d0*(
     .       +pjet(leptindex(1),4)*pjet(leptindex(2),4)
     .       -pjet(leptindex(1),1)*pjet(leptindex(2),1)
     .       -pjet(leptindex(1),2)*pjet(leptindex(2),2)
     .       -pjet(leptindex(1),3)*pjet(leptindex(2),3)))
	  if (mll .gt. mllcut) then
	    gencuts=.true.
	    return
	  endif
	endif
      endif
 
c--- lepton-photon separation 
      if ((countlept .ge. 1) .and. (countgamm .ge. 1)) then
        do j=1,countgamm
        do k=1,countlept
          if (R(pjet,gammindex(j),leptindex(k)) .lt. Rllmin) then
            gencuts=.true.
            return
          endif
        enddo
        enddo
      endif
 
c--- if there are no cuts on the jets - or no jets - we are done      
      if ((Rjlmin .le. 0d0) .and. (delyjjmin .le. 0d0)
     . .and. (gammcone .le. 0d0) .and. (gammcut .le. 0d0)) return
      if ((njets .eq. 0) .and. (countb .eq. 0)) return
      
c--- identify the jets
      countjet=0      
      do j=3,mxpart
        if (     (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .      .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo

c--- countjet will pick up the extra 'pp' needed for the real piece,
c--- therefore we should subtract 1 from this number     
      if (countjet .gt. njets) countjet=countjet-1

      if ((njets .ne. countjet) .and. (notag .eq. 0)) then
        write(6,*) 'Something is wrong in gencuts.f -'
        write(6,*) 'countjet = ',countjet,' BUT njets = ',njets
        stop
      endif
      
c--- extra cut on eta(2nd jet) for H(->WW)+jet search
      if ((hwwjetcuts) .and. (countjet .ge. 2)) then
        if (pt(jetindex(1),pjet) .gt. pt(jetindex(2),pjet)) then
          if (abs(etarap(jetindex(2),pjet)) .lt. etajet2cut) then
            gencuts=.true.
            return
          endif
        else
          if (abs(etarap(jetindex(1),pjet)) .lt. etajet2cut) then
            gencuts=.true.
            return
          endif
        endif
      endif

c--- jet-lepton separation (if there are 1 or more jets and leptons)
      if ((njets .gt. 0) .and. (countlept .gt. 0)) then
        do j=1,countlept
        do k=1,njets
          if (R(pjet,leptindex(j),jetindex(k)) .lt. Rjlmin) then
            gencuts=.true.
            return
          endif
        enddo
        enddo
      endif
      
c--- jet-photon separation (if there are 1 or more jets and photons)
      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
        do j=1,countgamm
        do k=1,njets
          if (R(pjet,gammindex(j),jetindex(k)) .lt. Rjlmin) then
            gencuts=.true.
            return
          endif
        enddo
        enddo
      endif
      
c--- photon/hadron isolation     
      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
        do j=1,countgamm
c--- Frixione cut, hep-ph/9801442
          if (njets .gt. 2) then
            write(6,*) 'Photon-hadron isolation not coded for njets > 2'
            stop
          endif
          do k=1,njets
            delta(k)=R(pjet,gammindex(j),jetindex(k))
            ptjet(k)=pt(jetindex(k),pjet)
          enddo
          pntr=1
          if (delta(2) .lt. delta(1)) pntr=2
          if (njets .eq. 1) delta(2)=gammcone   
          discr=(1d0-dcos(delta(pntr)))/(1d0-dcos(gammcone))
          if (ptjet(pntr) .gt. discr*pt(gammindex(j),pjet)) then
            gencuts=.true.
            return
          endif 
          if (njets .ge. 2) then
            discr=(1d0-dcos(delta(3-pntr)))/(1d0-dcos(gammcone))
            if (ptjet(1)+ptjet(2) .gt. discr*pt(gammindex(j),pjet))then
              gencuts=.true.
              return
            endif
          endif
c--- optional jet-veto
c          if (   (pt(4+countgamm+1,pjet) .gt. 50d0)    
c     .     .and. (abs(etarap(4+countgamm+1,pjet)) .lt. 2.5d0)) then
c            gencuts=.true.
c          endif                     
c--- de-Florian,Signer cut
c          do nu=1,2
c            sumjetpt(nu)=0d0
c          enddo
c          do k=1,njets
c            if (R(pjet,gammindex(j),4+countgamm+k) .lt. gammcone) then
c              do nu=1,2
c              sumjetpt(nu)=sumjetpt(nu)+pjet(4+countgamm+k,nu)
c              enddo
c            endif
c          enddo
c          if ( dsqrt(sumjetpt(1)**2+sumjetpt(2)**2) 
c     .    .gt. gammcut*pt(gammindex(j),pjet) ) then
c            gencuts=.true.
c          endif
        enddo
      endif  
      
c--- WBF-style cuts (if there are 2 or more jets)
      if ((njets .gt. 1)) then
c--- jet-jet rapidity separation 
c--- j and k point to the two highest pt ('tagging') jets
        j=1
        k=2
        if (njets .eq. 3) then
          if     ( pt(jetindex(1),pjet) .lt.
     .      min(pt(jetindex(2),pjet),pt(jetindex(3),pjet)) ) then
            j=2
            k=3
          elseif ( pt(jetindex(2),pjet) .lt.
     .      min(pt(jetindex(1),pjet),pt(jetindex(3),pjet)) ) then
            j=1
            k=3
          endif
        endif
        if (abs(etarap(jetindex(j),pjet)-etarap(jetindex(k),pjet))
     .         .lt. delyjjmin) then
          gencuts=.true.
          return
        endif

c--- Requirement that jets be in opposite hemispheres
        if (jetsopphem) then
          if(etarap(jetindex(j),pjet)*etarap(jetindex(k),pjet) .gt. 0d0)
     .       then
            gencuts=.true.
            return
          endif
        endif

        if (lbjscheme .gt. 0) then
c--- Cut to require lepton to be between jets
          etabuffer=dble(lbjscheme-1)*Rcut
          etaj=etarap(jetindex(j),pjet)
          etak=etarap(jetindex(k),pjet)
          do pntr=1,countlept
            etalept=etarap(leptindex(pntr),pjet)
            if ( (etalept .lt. min(etaj,etak)+etabuffer) .or.
     .           (etalept .gt. max(etaj,etak)-etabuffer) ) then 
              gencuts=.true.
              return
            endif

c--- these lines impose an alternative WBF selection
c             if ( ((etaj .gt. 1.6d0) .and. (etaj .lt. 4.4d0))
c     .        .or.((etak .gt. 1.6d0) .and. (etak .lt. 4.4d0)) ) then
c             else
c               gencuts=.true.
c             endif             
          enddo
          
          
        endif

      endif
      
c-- cuts on b-quarks
      if (bbproc) then
        call getbs(pjet,ib1,ib2)
        if ( (abs(etarap(ib1,pjet)) .gt. etabjetmax)
     .  .or. (pt(ib1,pjet) .lt. ptbjetmin) ) gencuts=.true. 
        if ( (abs(etarap(ib2,pjet)) .gt. etabjetmax)
     .  .or. (pt(ib2,pjet) .lt. ptbjetmin) ) gencuts=.true. 
      endif

c--- completed basic cuts
C--- if there are jet-like particles (see above), do more cuts
      if (countb .gt. 0) then
        do j=1,countb
          jindex=bindex(j)
          if (          (pt(jindex,pjet) .lt. jetpt) .or.
     .           (abs(etarap(jindex,pjet)) .gt. jetrap)) then
            gencuts=.true.
            return
          endif
          do k=j+1,countb
            kindex=bindex(k)
            if ((r(pjet,jindex,kindex) .lt. rcut)) then
              gencuts=.true.
              return
            endif
          enddo
        enddo
      endif
 
      return

   99 format(1x,a29,f6.2,a17)
      
  999 write(6,*) 'Error reading gencuts.DAT'
      call flush(6)
      stop      

      end
 
 
 
 
