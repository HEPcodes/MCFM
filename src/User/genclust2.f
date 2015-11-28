      subroutine genclust2(q,R,njet,qfinal,jetlabel)
c--- clusters momenta using plabel to determine which 
c--- particles should be clustered. Forms njet jets according to
c--- the standard kT algorithm with cone size R.
c--- Furthermore, the clustered jets are only observed if
c--- pT(jet) > ptjetmin and y(jet) < yjetmax
c--- 
c--- qfinal is the final vector q1,.... q(4+njets)
c--- where non-jet four vectors are set equal to the incoming q 
      implicit none
      include 'constants.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision R,dijmin,dkmin,pt,ayrap,ptjetmin,yjetmax
      double precision bclustmass,ptjet
      integer nparti,npartf,njet,i,nu,iter,nmin1,nmin2,maxjet,nk,
     . npart,nqcdjets,nqcdstart,ajet,nproc,jetindex(mxpart)
      character jetlabel(mxpart)*2,plabel(mxpart)*2
      logical bbproc
      common/plabel/plabel
      common/npart/npart
      common/nqcdjets/nqcdjets,nqcdstart
      common/nproc/nproc
      parameter (ptjetmin=15d0,yjetmax=2.5d0)
      
      njet=0
      maxjet=0

      do i=1,mxpart
        do nu=1,4
        qfinal(i,nu)=0d0
        enddo
      enddo

c--- pick out jets: note that we search to npart+2 to prevent
c--- extra gluons identified in 'chooser' contaminating lowest order      
      do i=3,npart+2
      if ((plabel(i) .eq. 'pp') .or. (plabel(i) .eq. 'pj')
     ..or.(plabel(i) .eq. 'bq') .or. (plabel(i) .eq. 'ba')) then
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
      dkmin=dkmin*R
      
c      write(*,*) 'Comparing pair (',nmin1,',',nmin2,') value of'
c      write(*,*) 'dijmin = ',dijmin,' with ',nk,' value of dk = ',dkmin
      
c--- step3: compare the two ...      
      if (dijmin .lt. dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2
        call combine(q,qjet,nmin1,nmin2,jetlabel)
c--- combined object goes into nmin1, now shuffle nmin2 off the end 
        call swap(qjet,jetlabel,nmin2,maxjet)        
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
        call swap(qjet,jetlabel,njet,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter. lt. maxjet-1) goto 1
      
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
     .   .and.(plabel(i) .ne. 'bq') .and. (plabel(i) .ne. 'ba')) then
            qfinal(i,nu)=q(i,nu)
          endif
        enddo
      enddo
      
c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'Obtained ',njet,' jets'

c--- restore jets
      ajet=0
      do i=1,njet
        if ((ptjet(i,q,qjet) .gt. ptjetmin) .and.
     .      (ayrap(i,qjet) .lt. yjetmax)) then     
        ajet=ajet+1
        do nu=1,4
          qfinal(jetindex(ajet),nu)=qjet(i,nu)
        enddo
        jetlabel(ajet)=jetlabel(i)
        endif
      enddo

c--- if any jets are removed by eta and pt cuts, then njet=ajet
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

c--- set bbproc to TRUE if the process involves two b-jets
      if (
     .      (nproc .eq.  21)
     . .or. (nproc .eq.  26)
     . .or. (nproc .eq.  51)
     . .or. (nproc .eq.  52)
     . .or. (nproc .eq.  53)
     . .or. (nproc .eq.  73)
     . .or. (nproc .eq.  78)
     . .or. (nproc .eq.  84)
     . .or. (nproc .eq.  89)
     . .or. (nproc .eq.  91)
     . .or. (nproc .eq.  96)
     . .or. (nproc .eq.  101)
     . .or. (nproc .eq.  102)
     . .or. (nproc .eq.  123)
     . .or. (nproc .eq.  151)
     . .or. (nproc .eq.  152)
     . .or. (nproc .eq.  161)
     . .or. (nproc .eq.  171)
     . ) then
        bbproc=.true.
      else
        bbproc=.false.
      endif

c--- check that 5 and 6 are b and b-bar
      if ((bbproc) .and. (bclustmass(njet,qfinal,jetlabel) .eq. 0d0))
     .  njet=-1    
      
c      do i=1,njet
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause
c      write(*,*) 'Started with ',njet+nremoved,' and now have ',njet
c      write(*,*) 'Found ',njet,' jets'

      return
      end
            
