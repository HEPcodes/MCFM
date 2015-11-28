      subroutine genclust_kt(q,Rmin,njet,qfinal,jetlabel)
c--- clusters momenta using plabel to determine which 
c--- particles should be clustered. Forms njet jets according to
c--- the standard kT algorithm with cone size Rmin.
c--- Furthermore, the clustered jets are only observed if
c--- pT(jet) > ptjetmin and eta(jet) < etajetmax
c--- 
c--- qfinal is the final vector q1,.... q(4+njets)
c--- where non-jet four vectors are set equal to the incoming q 
      implicit none
      include 'constants.f'
      include 'bbproc.f'
      include 'limits.f'
      include 'npart.f'
      include 'jetcuts.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision Rmin,dijmin,dkmin,aetarap
      double precision bclustmass,ptjet,m56,m57,m67
      integer njet,i,nu,iter,nmin1,nmin2,maxjet,nk,
     . ajet,jetindex(mxpart),nproc,countb
      character jetlabel(mxpart)*2,plabel(mxpart)*2
      common/plabel/plabel
      common/nproc/nproc

      njet=0
      maxjet=0

      do i=1,mxpart
        do nu=1,4
        qfinal(i,nu)=0d0
        enddo
      enddo

c--- pick out jets: note that we search to npart+2 and to prevent
c--- extra gluons identified in 'chooser' contaminating lowest order
c--- we check that the energy of the particle is not zero     
      do i=3,npart+2
      if (((plabel(i) .eq. 'pp') .or. (plabel(i) .eq. 'pj')
     ..or.(plabel(i) .eq. 'bq') .or. (plabel(i) .eq. 'ba')
     ..or.(plabel(i) .eq. 'qj'))
     ..and. (q(i,4) .gt. 1d-10)) then
        maxjet=maxjet+1
        jetindex(maxjet)=i
        jetlabel(maxjet)=plabel(i)
        do nu=1,4
          qjet(maxjet,nu)=q(i,nu)
        enddo
      endif
      enddo

c--- for no partons, just switch q into qfinal
      if (maxjet .eq. 0) then
        do i=1,mxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        njet=0
        return
      endif

c--- skip clustering if we only have one parton  
      if (maxjet .eq. 1) goto 2

      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

c      write(*,*) 'iter ',iter
c      write(*,*) 'njet ',njet
c      write(*,*) 'maxjet ',maxjet

c--- step1: find (i,j) pair with lowest measure of all non-jets so far
      call findmind(q,qjet,iter,maxjet,dijmin,nmin1,nmin2)
 
c--- step2: find jet K with lowest Et
      call findminet(q,qjet,iter,maxjet,dkmin,nk)
      dkmin=dkmin*Rmin

c      write(*,*) 'Comparing pair (',nmin1,',',nmin2,') value of'
c      write(*,*) 'dijmin = ',dijmin,' with ',nk,' value of dk = ',dkmin
      
c--- step3: compare the two ...      
      if (dijmin .lt. dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2
        call combine(q,qjet,nmin1,nmin2,jetlabel)
c--- combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet(qjet,jetlabel,jetindex,nmin2,maxjet)        
        maxjet=maxjet-1
        iter=iter-1
c        do i=1,maxjet
c          do j=1,4
c            write(*,*) 'qjet(',i,',',nu,') = ',qjet(i,nu)
c          enddo
c        enddo
      else
c---  ... we've finished a jet
        njet=njet+1
c        write(*,*) 'Now swapping ',njet,' and ',nk
        call swapjet(qjet,jetlabel,jetindex,njet,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter .lt. maxjet-1) goto 1
      
 2    continue      
      njet=njet+1

c--- restore incoming partons
      do i=1,2
        do nu=1,4
          qfinal(i,nu)=q(i,nu)
        enddo
      enddo
c--- set all other momenta to zero and restore leptons
      do i=3,npart+2
        do nu=1,4
          qfinal(i,nu)=0d0
          if ((plabel(i) .ne. 'pp') .and. (plabel(i) .ne. 'pj')
     .   .and.(plabel(i) .ne. 'bq') .and. (plabel(i) .ne. 'ba')
     .   .and.(plabel(i) .ne. 'qj')) then
            qfinal(i,nu)=q(i,nu)
          endif
        enddo
      enddo
      
      
c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'AFTER CLUSTERING: Obtained ',njet,' jets'

c--- restore jets
      ajet=0
      do i=1,njet
c        write(*,*) 'Jet ',i,'(',jetlabel(i),')',jetindex(i)
c        write(*,*) 'pt: ',ptjet(i,q,qjet),' vs min. ',ptjetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs min. ',etajetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs max. ',etajetmax
        if ((ptjet(i,q,qjet) .gt. ptjetmin) .and.
     .      (aetarap(i,qjet)   .gt. etajetmin) .and.
     .      (aetarap(i,qjet)   .lt. etajetmax)) then     
        ajet=ajet+1
        do nu=1,4
          qfinal(jetindex(ajet),nu)=qjet(i,nu)
        enddo
        jetlabel(ajet)=jetlabel(i)
        endif
      enddo
      
c--- if no jets are removed by eta and pt cuts, then njet=ajet
      if (ajet .lt. njet) then
        do i=ajet+1,njet
          do nu=1,4
            qfinal(jetindex(i),nu)=0d0
          enddo
        enddo
        njet=ajet
      endif
      
c      write(*,*) '... and ',njet,' jets after pt and eta cuts'
c      do i=1,njet
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause

c--- check that 5 is a b (for bH process)
      if ((nproc/10 .eq. 14)) then
        countb=0
        if ((njet.ge. 1) .and. ((jetlabel(1) .eq. 'bq')
     .    .or. (jetlabel(1) .eq. 'ba'))) countb=1
        if ((njet.eq. 2) .and. ((jetlabel(2) .eq. 'bq')
     .    .or. (jetlabel(2) .eq. 'ba'))) countb=countb+1
        if ((njet .eq. 1) .and. (countb .eq. 0)) njet=-1
        if ((nproc .eq. 142) .and. (njet .eq. 2)
     .      .and. (countb .ne. 1)) njet=-1
        if ((nproc .eq. 145) .and. (njet .eq. 2)
     .      .and. (countb .ne. 2)) njet=-1
      endif
      
c--- check that 5 and 6 are b and b-bar (if appropriate)
      if ((bbproc) .and. (bclustmass(njet,qfinal,jetlabel) .eq. 0d0))
     .  njet=-1    

c--- check that 5 and 6 are not b and b-bar if there are 2 jets 
      if (((nproc .eq. 24) .or. (nproc .eq. 29)) .and. (njet .eq. 2) 
     . .and. (bclustmass(njet,qfinal,jetlabel) .ne. 0d0))
     .  njet=-1    

c--- perform m56 mass cut if there are 2 or more jets
      if (njet .ge. 2) then
        m56=(qfinal(5,4)+qfinal(6,4))**2
     .     -(qfinal(5,1)+qfinal(6,1))**2
     .     -(qfinal(5,2)+qfinal(6,2))**2
     .     -(qfinal(5,3)+qfinal(6,3))**2
        if (njet .ge. 3) then
        m57=(qfinal(5,4)+qfinal(7,4))**2
     .     -(qfinal(5,1)+qfinal(7,1))**2
     .     -(qfinal(5,2)+qfinal(7,2))**2
     .     -(qfinal(5,3)+qfinal(7,3))**2
        m67=(qfinal(6,4)+qfinal(7,4))**2
     .     -(qfinal(6,1)+qfinal(7,1))**2
     .     -(qfinal(6,2)+qfinal(7,2))**2
     .     -(qfinal(6,3)+qfinal(7,3))**2
        m56=max(m56,max(m57,m67))
        endif
        if ((m56 .lt. bbsqmin) .or. (m56 .gt. bbsqmax)) then
          njet=-1
        endif
      endif

c      if (nproc .eq. 140) then
c      write(6,*) njet
c      write(6,*) 'jetlabel(1)',jetlabel(1)
c      write(6,*) 'jetlabel(2)',jetlabel(2)
c      write(6,*) 'jetlabel(3)',jetlabel(3)
c      write(6,*) 'jetlabel(4)',jetlabel(4)
c      write(6,*) 'jetlabel(5)',jetlabel(5)
c      write(6,*) 'jetlabel(6)',jetlabel(6)
c      write(6,*) 'jetlabel(7)',jetlabel(7)
c      write(6,*) 'threebee(njet,jetlabel)',threebee(njet,jetlabel)

c      if (threebee(njet,jetlabel) .eqv. .false.) njet=-1    

c      endif

       
            
c      do i=1,njet
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause
c      write(*,*) 'Started with ',njet+nremoved,' and now have ',njet
c      write(*,*) 'Found ',njet,' jets'

      return
      end
