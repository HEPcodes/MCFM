      subroutine genclust_kt(q,Rmin,qfinal,isub)
c--- clusters momenta using plabel to determine which 
c--- particles should be clustered. Forms 'jets' jets according to
c--- the standard kT algorithm with cone size Rmin.
c--- Furthermore, the clustered jets are only observed if
c--- pT(jet) > ptjetmin and eta(jet) < etajetmax
c--- 
c--- qfinal is the final vector q1,.... q(4+jets)
c--- where non-jet four vectors are set equal to the incoming q 
      implicit none
      include 'constants.f'
      include 'bbproc.f'
      include 'limits.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision Rmin,dijmin,dkmin,aetarap
      double precision ptjet,m56,m57,m67
      integer i,nu,iter,nmin1,nmin2,maxjet,nk,
     . ajet,jetindex(mxpart),nproc,countb,nbq,nba,isub
      character*2 plabel(mxpart)
      common/plabel/plabel
      common/nproc/nproc

      jets=0
      maxjet=0

      do i=1,mxpart
        do nu=1,4
        qfinal(i,nu)=0d0
        enddo
      enddo

c--- pick out jets: note that we search to npart+2-isub, to get the
c--- number of particles right. Note that isub=0 for all calls except
c--- the dipole contributions, where isub=1.   
      do i=3,npart+2-isub
      if ( (plabel(i) .eq. 'pp') .or. (plabel(i) .eq. 'pj')
     . .or.(plabel(i) .eq. 'bq') .or. (plabel(i) .eq. 'ba')
     . .or.(plabel(i) .eq. 'qj') ) then
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
        jets=0
        return
      endif

c--- skip clustering if we only have one parton  
      if (maxjet .eq. 1) goto 2

      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

c      write(*,*) 'iter ',iter
c      write(*,*) 'jets ',jets
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
        call combine(qjet,nmin1,nmin2)
c--- combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet(qjet,jetindex,nmin2,maxjet)        
        maxjet=maxjet-1
        iter=iter-1
c        do i=1,maxjet
c          do j=1,4
c            write(*,*) 'qjet(',i,',',nu,') = ',qjet(i,nu)
c          enddo
c        enddo
      else
c---  ... we've finished a jet
        jets=jets+1
c        write(*,*) 'Now swapping ',jets,' and ',nk
        call swapjet(qjet,jetindex,jets,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter .lt. maxjet-1) goto 1
      
 2    continue      
      jets=jets+1

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
     
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'

c--- restore jets
      ajet=0
      do i=1,jets
c        write(*,*) 'Jet ',i,'(',jetlabel(i),')',jetindex(i)
c        write(*,*) 'pt: ',ptjet(i,q,qjet),' vs min. ',ptjetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs min. ',etajetmin
c        write(*,*) 'aeta: ',aetarap(i,qjet),' vs max. ',etajetmax
        if ((ptjet(i,q,qjet) .gt. ptjetmin) .and.
     .      (aetarap(i,qjet) .gt. etajetmin) .and.
     .      (aetarap(i,qjet) .lt. etajetmax)) then     
        ajet=ajet+1
        do nu=1,4
          qfinal(jetindex(ajet),nu)=qjet(i,nu)
        enddo
        jetlabel(ajet)=jetlabel(i)
        endif
      enddo
      
c--- if no jets are removed by eta and pt cuts, then jets=ajet
      if (ajet .lt. jets) then
        do i=ajet+1,jets
          do nu=1,4
            qfinal(jetindex(i),nu)=0d0
          enddo
        enddo
        jets=ajet
      endif
      
c      write(*,*) '... and ',jets,' jets after pt and eta cuts'
c      do i=1,jets
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause

c--- check that 5 is a b (for bH process)
      if ((nproc/10 .eq. 14)) then
        countb=0
        if ((jets.ge. 1) .and. ((jetlabel(1) .eq. 'bq')
     .    .or. (jetlabel(1) .eq. 'ba'))) countb=1
        if ((jets.eq. 2) .and. ((jetlabel(2) .eq. 'bq')
     .    .or. (jetlabel(2) .eq. 'ba'))) countb=countb+1
        if ((jets .eq. 1) .and. (countb .eq. 0)) jets=-1
        if ((nproc .eq. 142) .and. (jets .eq. 2)
     .      .and. (countb .ne. 1)) jets=-1
        if ((nproc .eq. 145) .and. (jets .eq. 2)
     .      .and. (countb .ne. 2)) jets=-1
      endif
      
c--- check that 5 and 6 are b and b-bar (if appropriate)
      if (bbproc) then
        call getbs(qfinal,nbq,nba)
        if ((nbq .eq. 0) .or. (nba .eq. 0)) jets=-1    
      endif

c--- check that 5 and 6 are not b and b-bar if there are 2 jets 
c      if (((nproc .eq. 24) .or. (nproc .eq. 29)) .and. (jets .eq. 2) 
c     . .and. (bclustmass(qfinal) .ne. 0d0))
c     .  jets=-1    

c--- perform m56 mass cut if there are 2 or more jets
      if (jets .ge. 2) then
        m56=(qfinal(5,4)+qfinal(6,4))**2
     .     -(qfinal(5,1)+qfinal(6,1))**2
     .     -(qfinal(5,2)+qfinal(6,2))**2
     .     -(qfinal(5,3)+qfinal(6,3))**2
        if (jets .ge. 3) then
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
          jets=-1
        endif
      endif

c      if (nproc .eq. 140) then
c      write(6,*) jets
c      write(6,*) 'jetlabel(1)',jetlabel(1)
c      write(6,*) 'jetlabel(2)',jetlabel(2)
c      write(6,*) 'jetlabel(3)',jetlabel(3)
c      write(6,*) 'jetlabel(4)',jetlabel(4)
c      write(6,*) 'jetlabel(5)',jetlabel(5)
c      write(6,*) 'jetlabel(6)',jetlabel(6)
c      write(6,*) 'jetlabel(7)',jetlabel(7)
c      write(6,*) 'threebee()',threebee()

c      if (threebee() .eqv. .false.) jets=-1    

c      endif

       
            
c      do i=1,jets
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause
c      write(*,*) 'Started with ',jets+nremoved,' and now have ',jets
c      write(*,*) 'Found ',jets,' jets '

      return
      end
