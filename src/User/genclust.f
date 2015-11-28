      subroutine genclust(q,nparti,npartf,R,njet,qfinal,jetlabel)
c--- clusters momenta from momentum number nparti till 
c--- momentum number npartf in q into njet jets according to
c--- standard kT algorithm with cone size R
c--- furthermore, the clustered jets are only observed if
c--- pT(jet) > ptjetmin and y(jet) < yjetmax
c--- 
c--- Quantities returned are qfinal, njet, parts
c--- parton content of jets is encoded in parts(njet) 
c--- 10**n representing original parton label n
c--- qfinal is the final vector q1,q2,q3,q4,.... q(4+njets)
c--- first four vectors are set equal to the incoming q 
      implicit none
      include 'constants.f'
      include 'bbproc.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision R,dijmin,dkmin,pt,ayrap,ptjetmin,yjetmax
      double precision bclustmass
      integer nparti,npartf,njet,i,nu,iter,nmin1,nmin2,maxjet,nk,
     . npart,nqcdjets,nqcdstart,nremoved,nproc
      character jetlabel(mxpart)*2,plabel(mxpart)*2
      common/plabel/plabel
      common/npart/npart
      common/nqcdjets/nqcdjets,nqcdstart
      common/nproc/nproc
      parameter (ptjetmin=15d0,yjetmax=2.5d0)
      
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
          if (npart+2 . gt. nparti+nqcdjets) then
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
     .      (ayrap(nparti+i-1,qfinal) .gt. yjetmax)) then
          nremoved=nremoved+1
c          write(*,*) ' pt ',pt(nparti+i-1,qfinal),ptjetmin
c          write(*,*) 'eta ',ayrap(nparti+i-1,qfinal),yjetmax
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
            
      subroutine findmind(p,pjet,pjetmin,pjetmax,dijmin,nmin1,nmin2)
c--- this finds the minimum dij for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),dijmin,dij,d
      integer pjetmin,pjetmax,nmin1,nmin2,i,j
      
      dijmin=1d9

      do i=pjetmin,pjetmax
        do j=i+1,pjetmax
          if (i .ne. j) then
            d=dij(p,pjet,i,j)
            if (d .lt. dijmin) then
              dijmin=d
              nmin1=i
              nmin2=j
            endif
          endif
        enddo
      enddo
      
      if (dijmin. eq. 1d9) then
        write(*,*) 'Error in dij minimum-finding routine'
        stop
      endif
      
      return
      end
      
      subroutine findminet(p,pjet,pjetmin,pjetmax,dkmin,nk)
c--- this finds the minimum dkmin for pjet indices pjetmin through pjetmax
c--- returns dijmin and indices of minimum in (nmin1,nmin2)
C--- calculate the beam proto-jet separation see NPB406(1993)187, Eqn. 7
C--- in  practice this is just the minimum ptsq of protojets       
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),dkmin,dk,ptjet
      integer pjetmin,pjetmax,nk,i
      
      dkmin=1d9
      
      do i=pjetmin,pjetmax
        dk=ptjet(i,p,pjet)
        if (dk .lt. dkmin) then
          dkmin=dk
          nk=i
        endif
      enddo
      
      if (dkmin. eq. 1d9) then
        write(*,*) 'Error in dk minimum-finding routine'
        stop
      endif
      
      return
      end
      
      double precision function dij(p,pjet,i,j)
C---calculate the proto-jet separation see NPB406(1993)187, Eqn. 7
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),pjet(mxpart,4),pti,ptj,phii,phij,
     . yi,yj,ptjet,yrap,r
      
      pti=ptjet(i,p,pjet)
      ptj=ptjet(j,p,pjet)

c--- old method - bad because (phii-phij) can be > pi       
c      yi=yrap(i,pjet)
c      yj=yrap(j,pjet)

c      phii=atan2(pjet(i,1),pjet(i,2))
c      phij=atan2(pjet(j,1),pjet(j,2))
      
c      dij=dsqrt((yi-yj)**2+(phii-phij)**2)

c--- new method - r() calculates true value of 0 < (phi-phij) < pi
      dij=r(pjet,i,j)
            
      dij=dij*min(pti,ptj)
      
      return
      end
      
      subroutine combine(p,pjet,i,j,jetlabel)
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),pjet(mxpart,4)
      character jetlabel(mxpart)*2
      
c--Run II prescription
      pjet(i,1)=pjet(i,1)+pjet(j,1)
      pjet(i,2)=pjet(i,2)+pjet(j,2)
      pjet(i,3)=pjet(i,3)+pjet(j,3)
      pjet(i,4)=pjet(i,4)+pjet(j,4)


      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'qj'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'qj'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'qj'))) then
        jetlabel(i)='ba'
        return
      endif

      return
      end
      
      subroutine combine_snowmass(p,pjet,i,j,jetlabel)
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),pjet(mxpart,4),ptjetij,yjet,phijet,
     . ejet,ptjet,yrap,pti,ptj,yi,yj,phii,phij
      character jetlabel(mxpart)*2
      
C----Snowmass style prescripton 
      pti=ptjet(i,p,pjet)
      ptj=ptjet(j,p,pjet)
      
      yi=yrap(i,pjet)
      yj=yrap(j,pjet)
      
      phii=atan2(pjet(i,1),pjet(i,2))
      phij=atan2(pjet(j,1),pjet(j,2))
c

      ptjetij=pti+ptj
      yjet=(pti*yi+ptj*yj)/ptjetij
      phijet=(pti*phii+ptj*phij)/ptjetij
      ejet=exp(yjet)
      
      pjet(i,1)=ptjetij*dsin(phijet)
      pjet(i,2)=ptjetij*dcos(phijet)
      pjet(i,3)=ptjetij*(ejet-1d0/ejet)/2d0
      pjet(i,4)=ptjetij*(ejet+1d0/ejet)/2d0


      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='bq'
        return
      endif
      if (((jetlabel(i) .eq. 'ba') .and. (jetlabel(j) .eq. 'pp'))
     ..or.((jetlabel(j) .eq. 'ba') .and. (jetlabel(i) .eq. 'pp'))) then
        jetlabel(i)='ba'
        return
      endif
      if (((jetlabel(i) .eq. 'bq') .and. (jetlabel(j) .eq. 'ba'))
     ..or.((jetlabel(j) .eq. 'bq') .and. (jetlabel(i) .eq. 'ba'))) then
        jetlabel(i)='pp'
        return
      endif

      return
      end
      
      subroutine shuffle(pjet,nmin,nmax)
c--- shuffles jets nmin..nmax-1 in pjet down by 1 index
      implicit none
      include 'constants.f'
      integer i,j,nmin,nmax
      double precision pjet(mxpart,4)
      
      if (nmin .eq. nmax) return
      
      do i=nmin,nmax-1
        do j=1,4
          pjet(i,j)=pjet(i+1,j)
        enddo
      enddo
      
      return
      end
      
      subroutine swap(pjet,jetlabel,i,j)
c--- swaps jets i..j in pjet
      implicit none
      include 'constants.f'
      integer i,j,k
      double precision pjet(mxpart,4),tmp
      character jetlabel(mxpart)*2,chartmp*2
      
      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
      
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp
      
      return
      end
      
      double precision function ptjet(j,p,pjet)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4),pjet(mxpart,4),dotjet
      ptjet=sqrt(2d0*dotjet(p,1,pjet,j)*dotjet(p,2,pjet,j)
     . /dotjet(p,1,p,2))
      return
      end

      double precision function dotjet(p,i,pjet,j)
C---Dot the ith vector p with the jth vector pjet
      implicit none
      include 'constants.f'
      integer i,j
      double precision p(mxpart,4),pjet(mxpart,4)
       
      dotjet=p(i,4)*pjet(j,4)-p(i,1)*pjet(j,1)
     .      -p(i,2)*pjet(j,2)-p(i,3)*pjet(j,3)

      return
      end
      
      double precision function bclustmass(jets,pjet,jetlabel)
      implicit none
      include 'constants.f'
      integer i,jets,nbq,nba
      double precision pjet(mxpart,4)
      character jetlabel(mxpart)*2
      
c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      bclustmass=0d0
      nbq=0
      nba=0
      
      do i=1,jets
        if (jetlabel(i) .eq. 'bq') nbq=i+4
        if (jetlabel(i) .eq. 'ba') nba=i+4
      enddo

      if ((nbq .eq. 0) .or. (nba .eq. 0)) return
      
      bclustmass=(pjet(nbq,4)+pjet(nba,4))**2
      do i=1,3
        bclustmass=bclustmass-(pjet(nbq,i)+pjet(nba,i))**2
      enddo

      return
      end
       
