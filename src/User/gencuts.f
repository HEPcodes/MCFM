      logical function gencuts(p,pjet,njets)
************************************************************************
*   Author: J.M. Campbell, 5th December 2001                           *
*                                                                      *
*   This routine imposes a generic set of cuts that can be applied     *
*   to all processes in process.DAT                                    *
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
      include 'constants.f'
      logical first,passedlept
      character*2 plabel(mxpart)
      integer njets,j,k,nu,countb,bindex(mxpart),jindex,kindex
      integer countlept,leptindex(mxpart),countgamm,gammindex(mxpart),
     . countjet,jetindex(mxpart),pntr
      double precision p(mxpart,4),pjet(mxpart,4),etvec(4),sumjetpt(2)
      double precision leptpt,leptrap,misspt,jetpt,jetrap,gammpt,gammrap
      double precision pt,etarap,etmiss,evtmisset,R,Rcut,gammcone,
     . gammcut
      double precision Rjlmin,Rllmin,delyjjmin,leptpt2,leptrap2
      double precision rapj,delta(mxpart),discr,ptjet(mxpart)
      logical newinput
      common/newinput/newinput
      common/leptcuts/leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut
      common/plabel/plabel
      common/rcut/Rcut
************************************************************************
*     Set-up the jet-like cut parameters here                          *
      parameter (jetpt=15d0,jetrap=2d0)
************************************************************************
      data first/.true./

      gencuts=.false.
      
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
      
      if (newinput .eqv. .false.) then
c--- read-in cuts from file gencuts.DAT
        open(unit=21,file='gencuts.DAT',status='old',err=999)
        call checkversion(21,'gencuts.DAT')
        read(21,*) leptpt
        read(21,*) leptrap
        read(21,*) misspt
        read(21,*) leptpt2
        read(21,*) leptrap2
        read(21,*) Rjlmin
        read(21,*) Rllmin
        read(21,*) delyjjmin
        read(21,*) gammpt
        read(21,*) gammrap
        read(21,*) gammcone
        read(21,*) gammcut
        close(21)
      endif
      
      write(6,*)
      write(6,*)  '****************** Generic cuts ********************'
      write(6,*)  '*                                                  *'
      write(6,99) '*        pt(lepton)      >   ',leptpt,
     .                ' GeV            *'
      write(6,99) '*      |rap(lepton)|     <   ',leptrap,
     .                '                *'
      write(6,99) '*       pt(missing)      >   ',misspt,
     .                ' GeV            *'
      if ((leptpt2 .ne. 0d0) .or. (leptrap2 .ne. 0d0)) then
      write(6,99) '*     pt(2nd+ lepton)    >   ',leptpt2,
     .                ' GeV            *'
      write(6,99) '*   |rap(2nd+ lepton)|   <   ',leptrap2,
     .                '                *'
      endif
      write(6,99) '*      R(jet,lepton)     >   ',Rjlmin,
     .                '                *'
      write(6,99) '*     R(lepton,lepton)   >   ',Rllmin,
     .                '                *'
      write(6,99) '*     |y(jet1)-y(jet2)|  >   ',delyjjmin,
     .                '                *'
      write(6,99) '*   pt(photon)           >   ',gammpt,
     .                '                *'
      write(6,99) '*   rap(photon)          <   ',gammrap,
     .                '                *'
      write(6,99) '*   pt(hadronic)         <   ',gammcut,
     .                '    pt(photon)  *'
      write(6,99) '*   (in cone around photon of',gammcone,
     .                ')               *'
      if (countb .gt. 0) then
      write(6,*)  '*                                                  *'
      write(6,99) '*      pt(jet)       >   ',jetpt,
     .                ' GeV                *'
      write(6,99) '*    |rap(jet)|      <   ',jetrap,
     .                '                    *'
      write(6,99) '*   R(jet1,jet2)     >   ',Rcut,
     .                '                    *'
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
        if (     (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'pj')
     .      .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          countjet=countjet+1
          jetindex(countjet)=j
        endif
      enddo

c--- countjet will pick up the extra 'pp' needed for the real piece,
c--- therefore we should subtract 1 from this number     
      if (countjet .gt. njets) countjet=countjet-1

      if (njets .ne. countjet) then
        write(6,*) 'Something is wrong in gencuts.f -'
        write(6,*) 'countjet = ',countjet,' BUT njets = ',njets
        stop
      endif
      
c-- check that momenta 5...4+njets are indeed jets
c---(shifted by 1 for the number of photons)
c      if (njets .gt. 0) then
c        do j=5+countgamm,4+countgamm+njets
c          if ((plabel(j) .ne. 'pp') .and. (plabel(j) .ne. 'pj')
c     .   .and.(plabel(j) .ne. 'bq') .and. (plabel(j) .ne. 'ba')) then
c            write(6,*) 'Cannot identify jets in gencuts.f'
c            stop
c          endif
c        enddo
c      endif     

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
      
c--- jet-jet rapidity separation (if there are 2 or more jets)
      if ((njets .gt. 1)) then
c--- j and k point to the two highest pt jets
        j=1
        k=2
        if (njets. eq. 3) then
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
 
 
 
 
