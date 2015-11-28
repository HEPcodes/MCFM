      subroutine genclust_cone(q,R,qfinal,isub)
c---  clusters momenta using plabel to determine which 
c---  particles should be clustered. Forms 'jets' jets according to
c---  the Run II cone algorithm with cone size R.
c---  Furthermore, the clustered jets are only observed if
c---  pT(jet) > ptjetmin and y(jet) < etajetmax
c--- 
c---  qfinal is the final vector q1,.... q(4+jets)
c---  where non-jet four vectors are set equal to the incoming q 
      implicit none
      include 'constants.f'
      include 'bbproc.f'
      include 'limits.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      double precision q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4)
      double precision R,aetarap
      double precision bclustmass,ptjet,m56,m57,m67
      integer i,j,k,nu,iter,maxjet,
     . ajet,jetindex(mxpart),nproc,countb,nbq,nba,isub
      character*2 plabel(mxpart)
      double precision protoq(mxpart,4),deltar,deltarj,et,etmax,net,
     . qshared(4),sharedet
      integer maxproto,protoc(mxpart,0:mxpart),eti,shared,
     . sharedc(mxpart),ni
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
      if (maxjet .eq. 1) then
        jets=1
        do nu=1,4
          qfinal(1,nu)=qjet(1,nu)
        enddo
        goto 2
      endif
      
c--- set up the proto-jets
      maxproto=0
      do i=1,maxjet
        maxproto=maxproto+1
        protoc(maxproto,0)=1
        protoc(maxproto,1)=i
        do nu=1,4
          protoq(maxproto,nu)=qjet(i,nu)
        enddo
      enddo
      do i=1,maxjet
        do j=i+1,maxjet
          maxproto=maxproto+1
          protoc(maxproto,0)=2
          protoc(maxproto,1)=i
          protoc(maxproto,2)=j
          do nu=1,4
            protoq(maxproto,nu)=qjet(i,nu)+qjet(j,nu)
          enddo
          if (  (deltar(maxproto,i,protoq) .gt. R)
     .     .or. (deltar(maxproto,j,protoq) .gt. R)) then
            maxproto=maxproto-1
          endif
        enddo
      enddo
      if (maxjet .gt. 2) then
      do i=1,maxjet
        do j=i+1,maxjet
          do k=j+1,maxjet
            maxproto=maxproto+1
            protoc(maxproto,0)=3
            protoc(maxproto,1)=i
            protoc(maxproto,2)=j
            protoc(maxproto,3)=k
            do nu=1,4
              protoq(maxproto,nu)=qjet(i,nu)+qjet(j,nu)+qjet(k,nu)
            enddo
            if (  (deltar(maxproto,i,protoq) .gt. R)
     .       .or. (deltar(maxproto,j,protoq) .gt. R)
     .       .or. (deltar(maxproto,k,protoq) .gt. R)) then
              maxproto=maxproto-1
            endif
          enddo
        enddo
      enddo
      endif
      if (maxjet .gt. 3) then
       write(6,*) 'Too many jets for this version of the cone algorithm'
       stop
      endif
      
c      write(6,*) 'Found ',maxproto,' proto-jets'
c      stop
      
      jets=0
      
      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

      if (maxproto .eq. 0) goto 2

c--- find the highest Et proto-jet
      eti=0
      etmax=-1d0
      do i=1,maxproto
        et=dsqrt(protoq(i,1)**2+protoq(i,2)**2)
        if (et .gt. etmax) then
          eti=i
          etmax=et
        endif
      enddo
      
c      write(6,*) 'Max Et proto-jet is ',eti
      
c--- check to see if any partons are shared by this proto-jet
      shared=0
      do i=1,maxproto
        sharedc(i)=0
        if (i .ne. eti) then
          do j=1,protoc(i,0)
            do k=1,protoc(eti,0)
              if (protoc(i,j) .eq. protoc(eti,k)) then
                shared=shared+1
                sharedc(i)=1
              endif
            enddo
          enddo
        endif
      enddo
      
      if (shared .eq. 0) then
c-- proto-jet does not share any partons - move it to qfinal and repeat
        jets=jets+1
        do nu=1,4
          qfinal(jets,nu)=protoq(eti,nu)
        enddo
c--- shuffle down the proto-jets
        do i=eti+1,maxproto
          do nu=1,4
            protoq(i-1,nu)=protoq(i,nu)
          enddo
          do j=0,mxpart
            protoc(i-1,j)=protoc(i,j)
          enddo
        enddo
        maxproto=maxproto-1
c        write(6,*) 'Found jet number ',jets
        goto 1
      endif

c--- a parton is shared: perform split/merge procedure
c      write(6,*) 'Need to do split/merge'

c--- calculate which proto-jet that shares has the highest Et      
      ni=0
      net=-1d0
      do i=1,maxproto
        et=dsqrt(protoq(i,1)**2+protoq(i,2)**2)
        if ((sharedc(i) .eq. 1) .and. (et .gt. net)) then
          ni=i
          net=et
        endif
      enddo
     
c--- calculate the shared Et
      do j=1,protoc(eti,0)
        do k=1,protoc(ni,0)
          if (protoc(eti,j) .eq. protoc(ni,k)) then
            do nu=1,4
              qshared(nu)=qshared(nu)+qjet(protoc(eti,j),nu)
            enddo
          endif
        enddo
      enddo
      sharedet=dsqrt(qshared(1)**2+qshared(2)**2)
      
c      write(6,*) 'Proto-jet is',eti
c      write(6,*) 'Highest et neighbour is',ni
c      write(6,*) 'Shared Et is',sharedet
c      write(6,*) 'Neighbour Et is',net

      if (sharedet/net .gt. 0.5d0) then
c---  we should merge the proto-jets
        do i=1,protoc(ni,0)
          shared=0
          do j=1,protoc(eti,0)
            if (protoc(ni,i) .eq. protoc(eti,j)) shared=1
          enddo
c--- add cells that are not shared
          if (shared .eq. 0) then
            protoc(eti,0)=protoc(eti,0)+1
            protoc(eti,protoc(eti,0))=protoc(ni,i)
            do nu=1,4
              protoq(eti,nu)=protoq(eti,nu)+qjet(protoc(ni,i),nu)
            enddo
          endif
        enddo
c--- shuffle down the proto-jets
        do i=ni+1,maxproto
          do nu=1,4
            protoq(i-1,nu)=protoq(i,nu)
          enddo
          do j=0,mxpart
            protoc(i-1,j)=protoc(i,j)
          enddo
        enddo
        maxproto=maxproto-1
c        write(6,*) 'Merged proto-jets',eti,' and ',ni
      else
c---  we should split the proto-jets
        do i=1,protoc(ni,0)
          shared=0
          do j=1,protoc(eti,0)
            if (protoc(ni,i) .eq. protoc(eti,j)) shared=j
          enddo
c--- add cells that are not shared
          if (shared .gt. 0) then
            if (deltarj(protoc(ni,i),ni,qjet,protoq)
     .     .lt. deltarj(protoc(ni,i),eti,qjet,protoq)) then
c--- shared cell is closer to neighbour, ni
              do j=shared+1,protoc(eti,0)
              protoc(eti,j-1)=protoc(eti,j)
              enddo
              protoc(eti,0)=protoc(eti,0)-1     
            else            
c--- shared cell is closer to original proto-jet, eti     
              do j=i+1,protoc(ni,0)
              protoc(ni,j-1)=protoc(ni,j)
              enddo
              protoc(ni,0)=protoc(ni,0)-1     
            endif
          endif
        enddo
c        write(6,*) 'Split proto-jets',eti,' and ',ni
      endif
      
c      pause
      goto 1                   
      
 2    continue 
 
c---- transfer qfinal --> qjet
      do i=1,jets
        do nu=1,4
          qjet(i,nu)=qfinal(i,nu)
        enddo
      enddo
            
c      write(6,*) 'Finished finding jets: got ',jets
c      pause
      
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
c        write(*,*) 'ay: ',aetarap(i,qjet),' vs max. ',etajetmax
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
c      write(*,*) 'Found ',jets,' jets'

      return
      end
