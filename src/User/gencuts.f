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
      logical first
      character*2 plabel(mxpart)
      integer njets,j,k,nu,countb,bindex(mxpart),jindex,kindex
      integer countlept,leptindex(mxpart),countgamm,gammindex(mxpart),
     . maxindex
      double precision p(mxpart,4),pjet(mxpart,4),etvec(4),sumjetpt(2)
      double precision leptpt,leptrap,misspt,jetpt,jetrap,gammpt,gammrap
      double precision pt,yrap,etmiss,evtmisset,R,Rcut,gammcone,gammcut
      double precision Rjlmin,Rllmin,delyjjmin,leptpt2,leptrap2
      double precision ptj,maxpt
      common/plabel/plabel
      common/rcut/Rcut
************************************************************************
*     Set-up the jet-like cut parameters here                          *
      parameter (jetpt=15d0,jetrap=2d0)
************************************************************************
      data first/.true./
      save first,leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut

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
      maxindex=0
      maxpt=-1d0
c--- lepton pt and rapidity cuts
      do j=3,mxpart
        if (     (plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')
     .      .or. (plabel(j) .eq. 'ml') .or. (plabel(j) .eq. 'ma')
     .      .or. (plabel(j) .eq. 'tl') .or. (plabel(j) .eq. 'ta')) then
          countlept=countlept+1
          leptindex(countlept)=j
          ptj=pt(j,pjet)
          if (ptj .gt. maxpt) then
            maxindex=j
            maxpt=ptj
          endif
        endif
        if (plabel(j) .eq. 'ga') then
          countgamm=countgamm+1
          gammindex(countgamm)=j
        endif
      enddo
      
C     Basic pt and rapidity cuts for lepton
      if (countlept .gt. 0) then
          if (     (pt(maxindex,pjet) .lt. leptpt) .or.
     .      (abs(yrap(maxindex,pjet)) .gt. leptrap)) then
            gencuts=.true.
          endif
          if (countlept .gt. 1) then
            do j=1,countlept
              if (leptindex(j) .eq. maxindex) goto 77
              if (     (pt(leptindex(j),pjet) .lt. leptpt2) .or.
     .          (abs(yrap(leptindex(j),pjet)) .gt. leptrap2)) then
                gencuts=.true.
              endif
  77          continue
            enddo
          endif
      endif

C     Basic pt and rapidity cuts for photon
      if (countgamm .gt. 0) then
            do j=1,countgamm
              if (     (pt(gammindex(j),pjet) .lt. gammpt) .or.
     .          (abs(yrap(gammindex(j),pjet)) .gt. gammrap)) then
                gencuts=.true.
              endif
             enddo
      endif

c--- missing energy cut
      evtmisset=etmiss(pjet,etvec)
      if ((evtmisset. lt. misspt) .and. (evtmisset .ne. 0d0)) then
        gencuts=.true.
      endif
      
c--- lepton-lepton separation (if there are 2 or more leptons)
      if ((countlept .gt. 1)) then
        do j=1,countlept
        do k=j+1,countlept
          if (R(pjet,leptindex(j),leptindex(k)) .lt. Rllmin) then
            gencuts=.true.
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
          endif
        enddo
        enddo
      endif
 
c--- if there are no cuts on the jets, we are done      
      if ((Rjlmin .le. 0d0) .and. (delyjjmin .le. 0d0)
     . .and. (gammcone .le. 0d0) .and. (gammcut .le. 0d0)) return
      
c-- check that momenta 5...4+njets are indeed jets
c---(shifted by 1 for the number of photons)
      if (njets .gt. 0) then
        do j=5+countgamm,4+countgamm+njets
          if ((plabel(j) .ne. 'pp') .and. (plabel(j) .ne. 'pj')
     .   .and.(plabel(j) .ne. 'bq') .and. (plabel(j) .ne. 'ba')) then
            write(6,*) 'Cannot identify jets in gencuts.DAT'
            stop
          endif
        enddo
      endif     

c--- jet-lepton separation (if there are 1 or more jets and leptons)
      if ((njets .gt. 0) .and. (countlept .gt. 0)) then
        do j=1,countlept
        do k=1,njets
          if (R(pjet,leptindex(j),4+countgamm+k) .lt. Rjlmin) then
            gencuts=.true.
          endif
        enddo
        enddo
      endif
      
c--- jet-photon separation (if there are 1 or more jets and photons)
      if ((njets .gt. 0) .and. (countgamm .gt. 0)) then
        do j=1,countgamm
        do k=1,njets
          if (R(pjet,gammindex(j),4+countgamm+k) .lt. Rjlmin) then
            gencuts=.true.
          endif
        enddo
        enddo
      endif
      
c--- photon/hadron isolation     
      if ((njets. gt. 0) .and. (countgamm .gt. 0)) then
        do j=1,countgamm
          do nu=1,2
            sumjetpt(nu)=0d0
          enddo
          do k=1,njets
            if (R(pjet,gammindex(j),4+countgamm+k) .lt. gammcone) then
              do nu=1,2
              sumjetpt(nu)=sumjetpt(nu)+pjet(4+countgamm+k,nu)
              enddo
            endif
          enddo
          if ( dsqrt(sumjetpt(1)**2+sumjetpt(2)**2) 
     .    .gt. gammcut*pt(gammindex(j),pjet) ) then
            gencuts=.true.
          endif
        enddo
      endif  
      
c--- jet-jet rapidity separation (if there are 2 or more jets)
      if ((njets .gt. 1)) then
        do j=1,njets
        do k=j+1,njets
          if (abs(yrap(4+countgamm+j,pjet)-yrap(4+countgamm+k,pjet))
     .           .lt. delyjjmin) then
            gencuts=.true.
          endif
        enddo
        enddo
      endif
      
c--- completed basic cutsc
C--- if there are jet-like particles (see above), do more cuts
      if (countb .gt. 0) then
        do j=1,countb
          jindex=bindex(j)
          if (          (pt(jindex,pjet) .lt. jetpt) .or.
     .           (abs(yrap(jindex,pjet)) .gt. jetrap)) then
            gencuts=.true.
          endif
          do k=j+1,countb
            kindex=bindex(k)
            if ((r(pjet,jindex,kindex) .lt. rcut)) then
              gencuts=.true.
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
 
 
 
 
