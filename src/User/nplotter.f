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

      subroutine nplotter(vector,s,p,wt,switch)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'clustering.f'
      include 'bbproc.f'
      include 'jetlabel.f'
      integer idum_gasdev,n,switch,i5,i6,i7
      character tag*4
      double precision m56,m56_5,m56_10,m56_11,m56_12,m56_13,m56_15,
     . sigma,m34,m345,m346,m3456,m678,m47,etmiss,misset,m35,m45,
     . s(mxpart,mxpart),p(mxpart,4),eta,root,wt1,mtw
      double precision eta3,eta4,eta5,eta6,eta7,eta8,eta34,eta56
      double precision r34,r35,r45,r36,r46,r56,pt345
      double precision pt3,pt4,pt5,pt6,pt7,pt8,pt34,pt56,pt34a,pt34b
      double precision ptbbsq,ptbbpair,m56smw,gasdev,
     . pt,etarap,chi,cosphi,phi,var,vector(mxdim),wt,r,
     . kt12,kt14,kt15,kt24,kt25,kt56,costh,cosnew1,cosnew2,cosnew3
      double precision transm,transcm,dot,pttwo,etaraptwo,sdot30m,mttbar
      double precision phill,thetall,fphi,ftheta,mtsqlet,mt1,mt2,mll
      double precision c4,cosnchi,nchi,m56psm20,m56psm40,smearp,swap
      double precision clustermass,m56clust,costhdd,cosllet,coslpairet
      double precision pjet(mxpart,4),bclustmass,rn,etvec(4)
      double precision dsigdy,dsigdytmp,transm345,transm435
      double precision es17,es27,es56,es57,es67,etbin,etdoublebin
      double precision deta53,deta54,r57,r67
      double precision mbb,etab1,etab2,etanob,ptb1,ptb2,ptnob,rbb
      double precision detatags,dphitags
      integer nproc,eventpart,ib1,ib2
      logical first
      logical creatent,dswhisto
      common/outputflags/creatent,dswhisto
      common/nproc/nproc
      data idum_gasdev/34265765/
      data first/.true./
      save first,eta
      save es17,es27,es56,es57,es67
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
        pt34=0d0
        r34=0d0
        r35=0d0
        r45=0d0
        transm345=0d0
        transm435=0d0
        deta53=0d0
        deta54=0d0
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
        ptnob=0d0
        rbb=0d0
        goto 99
      else
        tag='plot'
      endif

c--- eventpart will contain the number of actual particles that have
c--- a defined momentum
c--- for lowest order and virtual terms switch=0 and eventpart=npart+2
c--- for real events switch=0 and eventpart=npart+2
c--- for real counter-events switch=1 and eventpart=npart+1
 
      m345=dsqrt(2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5)))
      m678=dsqrt(2d0*(dot(p,6,4)+dot(p,6,8)+dot(p,7,8)))
      m346=dsqrt(2d0*(dot(p,3,4)+dot(p,3,6)+dot(p,4,6)))
      m3456=dsqrt(2d0*(dot(p,3,4)+dot(p,3,5)+dot(p,3,6)
     .                +dot(p,4,5)+dot(p,4,6)+dot(p,5,6)))
      mttbar=2d0*(
     . +dot(p,3,4)+dot(p,3,5)+dot(p,3,6)+dot(p,3,7)+dot(p,3,8)
     .            +dot(p,4,5)+dot(p,4,6)+dot(p,4,7)+dot(p,4,8)
     .                       +dot(p,5,6)+dot(p,5,7)+dot(p,5,8)
     .                                  +dot(p,6,7)+dot(p,6,8)
     .                                             +dot(p,7,8))

      mttbar=dsqrt(mttbar)
c      write(6,*) 'm345',m345
c      write(6,*) 'm678',m678
c      write(6,*) 'mttbar',mttbar
c      pause
        eventpart=npart-switch+2
        if (jets .gt. 0) eventpart=4+jets
        if ((nproc .ge. 60) .and. (nproc .le. 89)) eventpart=6+jets
        if ((nproc .ge. 230) .and. (nproc .le. 239)) eventpart=2+jets
        
        if (switch .eq. 0) then
          es17=2d0*dot(p,1,7)  
          es27=2d0*dot(p,2,7)  
          es56=2d0*dot(p,5,6)  
          es57=2d0*dot(p,5,7)  
          es67=2d0*dot(p,6,7)  
        endif

        m34=dsqrt((p(3,4)+p(4,4))**2
     .           -(p(3,1)+p(4,1))**2
     .           -(p(3,2)+p(4,2))**2
     .           -(p(3,3)+p(4,3))**2)

        if (bbproc .and. clustering) then
c--- returns zero cluster mass if two b's are in one jet
          m56clust=dsqrt(bclustmass(p))
        else
          m56clust=dsqrt((p(5,4)+p(6,4))**2
     .                  -(p(5,1)+p(6,1))**2
     .                  -(p(5,2)+p(6,2))**2
     .                  -(p(5,3)+p(6,3))**2)
        endif  

c        r56=r(p,5,6)
c        r35=r(p,3,5)
c        r36=r(p,3,6)

c--generate a gaussian kick only for event
        if (switch .eq. 0) eta=gasdev(idum_gasdev)
        sigma=5d0
        m56_5=m56clust+sigma*eta
        sigma=10d0
        m56_10=m56clust+sigma*eta
        sigma=11d0
        m56_11=m56clust+sigma*eta
        sigma=12d0
        m56_12=m56clust+sigma*eta
        sigma=13d0
        m56_13=m56clust+sigma*eta
        sigma=15d0
        m56_15=m56clust+sigma*eta

        eta3=etarap(3,p)
        pt3=pt(3,p)
        eta4=etarap(4,p)
        pt4=pt(4,p)        
        eta34=etaraptwo(3,4,p)
        pt34=pttwo(3,4,p)
        r34=R(p,3,4)
        if ((eta34 .lt. -0.5d0) .or. (eta34 .gt. 0.5d0)) then
          pt34a=-1d0
        else
          pt34a=pt34
        endif
        if ((eta34 .lt. -0.1d0) .or. (eta34 .gt. 0.1d0)) then
          pt34b=-1d0
        else
          pt34b=pt34*5d0
        endif
                
        if (eventpart .gt. 4) then        
        eta5=etarap(5,p)
        pt5=pt(5,p)
        r35=R(p,3,5)
        r45=R(p,4,5)
        m35=dsqrt((p(3,4)+p(5,4))**2
     .           -(p(3,1)+p(5,1))**2
     .           -(p(3,2)+p(5,2))**2
     .           -(p(3,3)+p(5,3))**2)
        m45=dsqrt((p(4,4)+p(5,4))**2
     .           -(p(4,1)+p(5,1))**2
     .           -(p(4,2)+p(5,2))**2
     .           -(p(4,3)+p(5,3))**2)
        if (eventpart .gt. 5) then
          pt345=-pt(6,p)
        else
          pt345=0d0
        endif
 
        transm345=dsqrt((dsqrt(m45**2+pttwo(4,5,p)**2)+pt3)**2-pt345**2) 
        transm435=dsqrt((dsqrt(m35**2+pttwo(3,5,p)**2)+pt4)**2-pt345**2) 
        
        deta53=eta5-eta3
        deta54=eta5-eta4

c--- debug: transverse mass cut, this should be removed        
c        if (transm345 .lt. 90d0) wt=0d0
c--- debug: transverse mass cut, this should be removed        
        
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
        transm=2d0*dsqrt(pttwo(4,5,p)**2+2d0*dot(p,4,5))
        transcm=dsqrt(pttwo(4,5,p)**2+2d0*dot(p,4,5))+pttwo(3,6,p)
        endif

        if (eventpart .gt. 6) then        
        eta7=etarap(7,p)
        pt7=pt(7,p)
        r57=R(p,5,7)
        r67=R(p,6,7)
        endif

        if (eventpart .gt. 7) then        
        eta8=etarap(8,p)
        pt8=pt(8,p)
        endif
         
        misset=etmiss(p,etvec)

c--- if we're doing W/Z + jets then order the jets according to pt's:
c--- JET 5 = highest Et
c--- JET 6 = next-highest Et
c--- JET 7 = lowest Et
c--- one-jet processes are 11,16 and 41,42,43
c--- two-jet processes are 22 and 44
c--- three-jet process is 23

c--- case where we have 2 jets to order
      if ( (((nproc .eq. 22)
     .    .or.(nproc .eq. 27)
     .    .or.(nproc .eq. 44)
     .    .or.(nproc .eq. 217)
     .    .or.(nproc .eq. 219)
     .    .or.(nproc .eq. 272)
     .    .or.(nproc .eq. 273)
     .    .or.(nproc .eq. 145))
     .    .and. (jets .eq. 2))
     . .or.(((nproc .eq. 11).or.(nproc .eq. 46).or.(nproc .eq. 41)
     .                      .or.(nproc .eq. 42).or.(nproc .eq. 43))
     .    .and. (jets .eq. 2)) ) then
        if (pt6 .gt. pt5) then
          swap=pt5
          pt5=pt6
          pt6=swap
          swap=eta5
          eta5=eta6
          eta6=swap
        endif
        detatags = abs(eta6-eta5)
        dphitags = fphi(5,6,p)
      endif
      
c--- case where we have 3 jets to order
      if ( (((nproc .eq. 22)
     .    .or.(nproc .eq. 27)
     .    .or.(nproc .eq. 44)
     .    .or.(nproc .eq. 217)
     .    .or.(nproc .eq. 219)
     .    .or.(nproc .eq. 272)
     .    .or.(nproc .eq. 273))
     .    .and. (jets .eq. 3))
     . .or.  ((nproc .eq. 23).or.(nproc .eq. 28)
     . .or.  (nproc .eq. 218)) ) then
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
        eta5=etarap(i5,p)
        pt5=pt(i5,p)
        eta6=etarap(i6,p)
        pt6=pt(i6,p)
        eta56=etaraptwo(i5,i6,p)
        pt56=pttwo(i5,i6,p)
        r56=r(p,i5,i6)
        eta7=etarap(i7,p)
        pt7=pt(i7,p)
        m56clust=dsqrt((p(i5,4)+p(i6,4))**2
     .                -(p(i5,1)+p(i6,1))**2
     .                -(p(i5,2)+p(i6,2))**2
     .                -(p(i5,3)+p(i6,3))**2)
        if ( abs(eta6-eta5).gt.abs(eta7-eta6) .and.
     .       abs(eta6-eta5).gt.abs(eta7-eta5) ) then
          detatags = abs(eta6-eta5)
          dphitags = fphi(5,6,p)
        elseif ( abs(eta7-eta5).gt.abs(eta7-eta6) .and.
     .           abs(eta7-eta5).gt.abs(eta6-eta5) ) then
          detatags = abs(eta7-eta5)
          dphitags = fphi(5,7,p)
        else
          detatags = abs(eta7-eta6)
          dphitags = fphi(6,7,p)
        endif
      endif

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
            write(6,*) 'Unforeseen # of jets and b-quarks in nplotter.f'
            stop
          endif
        endif

      if ((nproc .eq. 60) .or. (nproc .eq. 61)) then
        etbin=etdoublebin(pt4,pt5)
      endif

   99 continue

c--- only fill the histograms if we're not creating n-tuples
      if (creatent .eqv. .false.) then
      n=1                  

c --- Histograms to monitor the weight distributions :
      call bookplot(n,tag,'      wt',wt,1.0d0,-1d-2,1d-2,2d-4,'lin')
      n=n+1
c      call bookplot(n,tag,' log10wt',dlog10(dabs(wt)),
c     . 1.0d0,-2d0,0d0,0.1d0,'lin')
     
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
      call bookplot(n,tag,'   eta34',eta34,wt,-5d0,5d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'    pt34',pt34,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'    pt34',pt34,wt,0d0,200d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'d/pt34^2',pt34,wt/2d0*pt34**3,
     .     20d0,200d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'pt34,eta=0',pt34a,wt,10d0,150d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'pt34,eta=0',pt34a,wt,20d0,480d0,40d0,'log')
      n=n+1
      call bookplot(n,tag,'     m34',m34,wt,10d0,140d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'    mT34',mtw,wt,20d0,120d0,2d0,'lin')
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
      call bookplot(n,tag,'trnsm345',transm345,wt,0d0,200d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'trnsm345',transm345,wt,90d0,990d0,1d2,'lin')
      n=n+1
      call bookplot(n,tag,'trnsm435',transm435,wt,0d0,200d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'trnsm435',transm435,wt,90d0,990d0,1d2,'lin')
      n=n+1
      call bookplot(n,tag,'  deta53',deta53,wt,-5d0,5d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'  deta54',deta54,wt,-5d0,5d0,0.4d0,'lin')
      n=n+1

c--- Kramer comparison
      call bookplot(n,tag,' pth(mk)',pt34,wt,9.2d0,101.2d0,9.2d0,'log')
      n=n+1
      call bookplot(n,tag,' ptb(mk)',pt5,wt,9.2d0,101.2d0,9.2d0,'log')
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
      call bookplot(n,tag,'|y5-y6|',
     .                    abs(eta5-eta6),wt,0d0,10d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'dphi(56)',fphi(5,6,p),wt,0d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'detatags',detatags,wt,0d0,10d0,0.4d0,'lin')
      n=n+1
      call bookplot(n,tag,'dphitags',dphitags,wt,0d0,4d0,0.2d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r36',r36,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r46',r46,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     r56',r56,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,0d0,200d0,5d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,0d0,1000d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,0d0,1000d0,50d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,0d0,200d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,0d0,400d0,10d0,'log')
      n=n+1
      call bookplot(n,tag,'     m56',m56clust,wt,71d0,111d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'   m56_5',m56_5,wt,20d0,160d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'  m56_10',m56_10,wt,20d0,160d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'  m56_15',m56_15,wt,20d0,160d0,5d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_10_a',m56_10,wt,86d0,114d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_11_b',m56_11,wt,94d0,126d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_12_c',m56_12,wt,103d0,137d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'m56_13_d',m56_13,wt,112d0,148d0,2d0,'lin')
      n=n+1
      call bookplot(n,tag,'    m345',m345,wt,50d0,250d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'    m346',m346,wt,50d0,250d0,10d0,'lin')
      n=n+1
      call bookplot(n,tag,'   m3456',m3456,wt,50d0,250d0,10d0,'lin')
      n=n+1
c      call bookplot(n,tag,'  mttbar',mttbar,wt,300d0,1d3,20d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'  transm',transm,wt,20d0,200d0,20d0,'lin')
c      n=n+1
c      call bookplot(n,tag,' transcm',transcm,wt,20d0,200d0,20d0,'lin')
c      n=n+1
      endif

      if (bbproc) then
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
      call bookplot(n,tag,'     rbb',rbb,wt,0d0,4d0,0.1d0,'lin')
      n=n+1
      endif

      if ((nproc .eq. 60) .or. (nproc .eq. 61)) then
      call bookplot(n,tag,'   etbin',etbin,wt,0.5d0,25.5d0,1d0,'lin')
      n=n+1
      endif

      if (eventpart .gt. 6) then

      if (bbproc) then
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

      else
c--- Book and fill ntuple if we're not doing histograms       
         call bookfill(tag,p,wt)       
      endif

      n=n-1

      return 
      end

c--- this is trying to be decay angle, but need to unboost i back to (jk)
c--- rest frame before doing this
c      double precision function sdot30m(i,j,k,p)
c      implicit none
c      include 'constants.f'
c      integer i,j,k,n
c      double precision p(mxpart,4),n1,n2
      
c      sdot30m=0d0
c      n1=0d0
c      n2=0d0
c      do n=1,3
c      sdot30m=p(i,n)*(p(j,n)+p(k,n))
c      n1=n1+p(i,n)**2
c      n2=n2+(p(j,n)+p(k,n))**2
c      enddo
      
c      sdot30m=sdot30m/sqrt(n1*n2)
c      
c      return
c      end
     
c      double precision function cosnchi(i,j,k,l,p)
c      implicit none
c      include 'constants.f'
c      integer i,j,k,l,n
c      double precision p(mxpart,4),cr1(3),cr2(3),n1,n2
      
c      call cross(p,i,j,cr1)
c      call cross(p,k,l,cr2)
      
c      n1=0d0
c      n2=0d0
c      cosnchi=0d0
c      do n=1,3
c        n1=n1+cr1(n)**2
c        n2=n2+cr2(n)**2
c        cosnchi=cosnchi+cr1(n)*cr2(n)
c      enddo
c
c      cosnchi=cosnchi/dsqrt(n1*n2)
      
c      return
c      end
      
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
            
c      double precision function smearp(i,j,p,sd)
c      implicit none
c      include 'constants.f'
c      include 'masses.f'     
c      integer i,j,k,idum_gasdev
c      double precision p(mxpart,4),r1(4),r2(4),gasdev,sm1,sm2,sd
c      data idum_gasdev/56735345/

c      sm1=1d0+gasdev(idum_gasdev)/sd
c      sm2=1d0+gasdev(idum_gasdev)/sd

c      do k=1,4
c        r1(k)=p(i,k)*sm1
c        r2(k)=p(j,k)*sm2
c      enddo
      
c      smearp=sqrt(2d0*(r1(4)*r2(4)-r1(1)*r2(1)-r1(2)*r2(2)-r1(3)*r2(3))
c     . +mb**2)
      
c      return
c      end
      
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
          
c      double precision function ftheta(n1,n2,p)
c      implicit none
c      include 'constants.f'
c      integer n1,n2
c      double precision p(mxpart,4)
c    
c      ftheta=p(n1,1)*p(n2,1)+p(n1,2)*p(n2,2)+p(n1,3)*p(n2,3)
c      ftheta=ftheta/dsqrt(p(n1,1)**2+p(n1,2)**2+p(n1,3)**2)
c      ftheta=ftheta/dsqrt(p(n2,1)**2+p(n2,2)**2+p(n2,3)**2)
c      if (ftheta .gt. 1d0) then
c        ftheta=0d0
c      elseif (ftheta .lt. -1d0) then
c        ftheta=pi
c      else
c        ftheta=dacos(ftheta)
c      endif
c   
c      return
c      end
      
c      double precision function mtsqlet(n,nm1,nm2,p)
c      implicit none
c      include 'constants.f'
c      integer n,nm1,nm2,i
c      double precision p(mxpart,4),misset(4),etmiss,coslet,pt,pttwo
      
c      if (nm2 .eq. 0) then
c        etmiss=pt(nm1,p)
c        do i=1,4
c          misset(i)=p(nm1,i)
c        enddo
c      else
c        etmiss=pttwo(nm1,nm2,p)
c        do i=1,4
c          misset(i)=p(nm1,i)+p(nm2,i)
c        enddo
c      endif
      
c      coslet=p(n,1)*misset(1)+p(n,2)*misset(2)
c      coslet=coslet/dsqrt(p(n,1)**2+p(n,2)**2)
c      coslet=coslet/dsqrt(misset(1)**2+misset(2)**2)
c      mtsqlet=2d0*pt(n,p)*etmiss*(1d0-coslet)
c      
c      return
c      end
            
      double precision function coslpairet(n1,n2,nm1,nm2,p)
      implicit none
      include 'constants.f'
      integer n1,n2,nm1,nm2,i
      double precision p(mxpart,4),misset(4),pp(4),coslet,pt,pttwo
      
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
      
c      double precision function deltarpjet(i,j,p,pjet)
c      implicit none
c      include 'constants.f'
c      double precision p(mxpart,4),pjet(mxpart,4),phi1,phi2,etarap,dphi
c      integer i,j
      
c      phi1=atan2(p(i,1),p(i,2))
c      phi2=atan2(pjet(j,1),pjet(j,2))
c      dphi=phi1-phi2
c      if (dphi .gt. pi) dphi=twopi-dphi
c      if (dphi .lt. -pi) dphi=twopi+dphi
c      deltarpjet=(etarap(i,p)-etarap(j,p))**2+dphi**2
c      deltarpjet=dsqrt(deltarpjet)
      
c      return
c      end
      
c      double precision function ptqfour(q,j,k,l,m)
c      implicit none
c      include 'constants.f'
c      integer j,k,l,m
c      double precision q(mxpart,4)
c      ptqfour=dsqrt((q(j,1)+q(k,1)+q(l,1)+q(m,1))**2
c     .             +(q(j,2)+q(k,2)+q(l,2)+q(m,2))**2)
c      return
c      end

          
