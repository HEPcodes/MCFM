      subroutine genclust(q,nparti,npartf,R,njet,qfinal,jetlabel)
c--- clusters momenta from momentum number nparti till 
c--- momentum number npartf in q into njet jets according to
c--- standard kT algorithm with cone size R
c--- furthermore, the clustered jets are only observed if
c--- pT(jet) > ptjetmin and y(jet) < etajetmax
c--- 
c--- Quantities returned are qfinal, njet, parts
c--- parton content of jets is encoded in parts(njet) 
c--- 10**n representing original parton label n
c--- qfinal is the final vector q1,q2,q3,q4,.... q(4+njets)
c--- first four vectors are set equal to the incoming q 
      implicit none
      include 'constants.f'
      include 'bbproc.f'
      include 'npart.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision R,dijmin,dkmin,pt,aetarap,ptjetmin,etajetmax
      double precision bclustmass
      integer nparti,npartf,njet,i,nu,iter,nmin1,nmin2,maxjet,nk,
     . nqcdjets,nqcdstart,nremoved,nproc
      character*2 plabel(mxpart)
      character*2 jetlabel(mxpart)
      common/plabel/plabel
      common/nqcdjets/nqcdjets,nqcdstart
      common/nproc/nproc
      parameter (ptjetmin=15d0,etajetmax=2.5d0)
      
      maxjet=npartf-nparti+1
      njet=0
      
C---calculate the four-jet momenta from i to maxjet and label each in parts
      do i=1,maxjet
        jetlabel(i)=plabel(nparti+i-1)
        do nu=1,4
          qjet(i,nu)=q(nparti+i-1,nu)
c          write(*,*) 'qjet(',i,',',nu,') = ',qjet(i,nu)
        enddo
      enddo

c--- for no partons, just switch q into qfinal
      if (npartf .lt. nparti) then
        do i=1,mxpart
          do nu=1,4
            qfinal(i,nu)=q(i,nu)
          enddo
        enddo
        njet=0
        return
      endif

c--- skip clustering if we only have one parton  
      if (npartf .eq. nparti) goto 2

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
        call combine(qjet,nmin1,nmin2)
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

      if (iter .lt. maxjet-1) goto 1
      
 2    continue      
      njet=njet+1

C---store total result in
       do nu=1,4
          do i=1,mxpart
          qfinal(i,nu)=0d0
          enddo
c----set first (nparti-1) momenta qfinal equal to original q's
          do i=1,nparti-1
          qfinal(i,nu)=q(i,nu)
          enddo
          do i=1,njet
          qfinal(nparti+i-1,nu)=qjet(i,nu)
          enddo
c----make sure we also set lepton momenta for processes such as 151 and 181.
c----(npart+2) is the total number of partons in the event, which can be
c----greater than (4+nqcdjets) only for real events (when it is a gluon, in
c----which case plabel='pp' and it is ignored) and events with trailing leptons
c----such as 151 and 181. This should be generic enough to handle future
c----extensions of process.DAT.
          if (npart+2 .gt. nparti+nqcdjets) then
            do i=nparti+nqcdjets,npart+2
              if (plabel(i) .ne. 'pp') then
                 qfinal(i,nu)=q(i,nu)
              endif
            enddo
          endif
      enddo

c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'Tried to cluster ',nparti,' to ',npartf
c      write(*,*) 'Obtained ',njet,' jets'

      nremoved=0
      do i=1,njet
        if ((pt(nparti+i-1,qfinal) .lt. ptjetmin) .or.
     .      (aetarap(nparti+i-1,qfinal) .gt. etajetmax)) then
          nremoved=nremoved+1
c          write(*,*) ' pt ',pt(nparti+i-1,qfinal),ptjetmin
c          write(*,*) 'eta ',aetarap(nparti+i-1,qfinal),etajetmax
        else
          jetlabel(i-nremoved)=jetlabel(i)
          do nu=1,4
            qfinal(nparti+i-1-nremoved,nu)=qfinal(nparti+i-1,nu)
          enddo
        endif
      enddo
      njet=njet-nremoved

c      write(*,*) '... and ',njet,' jets after pt and eta cuts'
c      do i=1,njet
c        write(*,*) i,jetlabel(i)
c      enddo
c      pause

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
            
