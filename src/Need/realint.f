      double precision function realint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'realwt.f'
      include 'clustering.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'flags.f'
      include 'efficiency.f'
      integer ih1,ih2,j,k,nd,nmax,nmin,jets,nproc,nvec,jj,kk     
      double precision vector(mxdim),W,val,xint
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,msqs(-nf:nf,-nf:nf)
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4),rcut
      double precision debugsmall
      character*32 debugmsg
      character*6 case
      integer nqcdjets,nqcdstart,i,kx
      character pdlabel*7,jetlabel(mxpart)*2
      common/nqcdjets/nqcdjets,nqcdstart
      common/parts/jets,jetlabel
      common/process/case
      common/xreal/xreal,xreal2
      logical bin,makecuts,first,madejetcuts,failed,cuts,bbproc,test1
      common/density/ih1,ih2
      common/energy/sqrts
      common/pdlabel/pdlabel
      common/bin/bin
      common/makecuts/makecuts
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/nproc/nproc
      common/rcut/rcut
      data first/.true./
      
      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      

      W=sqrts**2
      
      if (first) then
         write(6,*)
         write(6,*) 'nmin=',nmin,',nmax=',nmax
         write(6,*)
         first=.false.
      endif
      if     (
     .       (case .eq. 'W_only')
     .  .or. (case .eq. 'Z_only')
     .  .or. (case .eq. 'Hbbbar')
     . ) then
          npart=3
          if (new_pspace) then
          call gen3a(vector,p,pswt,*999)      
          else
          call gen3(vector,p,pswt,*999)      
          endif
      elseif ((case .eq. 'W_1jet') .or. (case .eq. 'Z_1jet')) then
          npart=4
          if (new_pspace) then
          call gen4a(vector,p,pswt,*999)      
          else
          call gen4(vector,p,pswt,*999)
          endif
      else
          npart=5
          if (new_pspace) then
          call gen5a(vector,p,pswt,*999)
          else
          call gen5(vector,p,pswt,*999)      
          endif
      endif
      nvec=npart+2
      
c--- DEBUG to test limit where p5 is soft
c      do j=1,4
c        p(8,j)=p(5,j)
c        p(5,j)=p(7,j)
c        p(7,j)=p(8,j)
c        p(8,j)=0d0
c      enddo     
c--- DEBUG      
      
      call dotem(nvec,p,s)

c      write(6,*) 's(1,5),s(1,6),s(2,5),s(2,6),s(5,6),s(3,4)',
c     . s(1,5),s(1,6),s(2,5),s(2,6),s(5,6),s(3,4)
c---impose cuts on final state
      call masscuts(s,*999)
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      
c---- DEBUG 
c---- generate collinear points that satisy the jet cuts    
c      call genclust2(p,rcut,jets,pjet,jetlabel)
c      if ((-s(1,7) .lt. 1d-1) .and. (jets .eq. nqcdjets)
c     . .and. (-p(1,4) .gt. 1d1) .and. (p(7,4) .gt. 1d1)) then
c      do j=1,7
c      write(6,77) j,p(j,1),p(j,2)
c      write(6,78) p(j,3),p(j,4)
c      enddo
c      pause
c      endif      
c      goto 998
c   77 format('      data p',i1,'/ ',f20.14,'d0,',f20.14,'d0,')
c   78 format('     .  ',f20.14,'d0,',f20.14,'d0/')
c---- DEBUG      
      
      
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 .gt. 1) .or. (xx2 .gt. 1)) then
         realint=0d0
         return
      endif

      call fdist(pdlabel,ih1,xx1,scale,fx1)
      call fdist(pdlabel,ih2,xx2,scale,fx2)

      flux=fbGeV2/(two*xx1*xx2*W)

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
     . .or. (nproc .eq.  151)
     . .or. (nproc .eq.  152)
     . .or. (nproc .eq.  161)
     . .or. (nproc .eq.  171)
     . ) then
        bbproc=.true.
      else
        bbproc=.false.
      endif

      if     (case .eq. 'Wbbbar') then
        call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (case .eq. 'Zbbbar') then
      call qqb_Zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (case .eq. 'WHbbar') then
        call qqb_wh_gs(p,msqc)     
        call qqb_wh_g(p,msq)      
      elseif (case .eq. 'ZHbbar') then
        call qqb_zh_gs(p,msqc)     
        call qqb_zh_g(p,msq)      
      elseif (case .eq. 'W_2jet') then
        call qqb_w2jet_g(p,msq)  
        call qqb_w2jet_gs(p,msqc)
c---- include this section if you want to check subtraction cancellation
        if (1 .eq. 2) then
        do j=1,4
        do k=1,9
           call coll2(p,k,j)
           write(*,*) 'Point ',j,':'
           call qqb_w2jet_g(p,msq)  
           call qqb_w2jet_gs(p,msqc) 
           do jj=-nf,nf
           do kk=-nf,nf
              msqs(jj,kk)=0d0
           enddo
           enddo
           do nd=1,ndmax
              do jj=1,7
              do kk=1,4
                 q(jj,kk)=ptilde(nd,jj,kk)
              enddo
              enddo
              call dotem(6,q,s)
              call genclust2(q,rcut,jets,pjet,jetlabel)
              if (jets .eq. nqcdjets) then
                 do jj=-nf,nf
                 do kk=-nf,nf
                    msqs(jj,kk)=msqs(jj,kk)+msqc(nd,jj,kk)
                 enddo
                 enddo
              endif
           enddo

c--- find smallest value of msq
           debugsmall=1d0
           do jj=-nf,nf
           do kk=-nf,nf
             if ((msq(jj,kk) .lt. debugsmall)
     .       .and. (msq(jj,kk) .gt. 0d0)) debugsmall=msq(jj,kk)        
           enddo
           enddo                  
        
           do jj=-nf,nf
           do kk=-nf,nf
           if (msq(jj,kk) .eq. 0d0) then
              if (msqs(jj,kk) .eq. 0d0) then
                 debugmsg='   OK   zero msq and subtraction'
              else
                 debugmsg=' FAILED subtraction with msq=0'
              endif
              write(*,65) jj,kk,'    n/a   ',msq(jj,kk),msqs(jj,kk),
     .        debugmsg
              goto 64
              endif

           if ((msq(jj,kk)/debugsmall) .lt. 1d4) then
              debugmsg='   OK   not singular'
              write(*,63) jj,kk,msq(jj,kk)/msqs(jj,kk),
     .                      msq(jj,kk),msqs(jj,kk),debugmsg
           else
               if (msqs(jj,kk) .eq. 0d0) then
                  debugmsg=' FAILED singular, no subtraction'
                  write(*,65) jj,kk,'    n/a   ',msq(jj,kk),msqs(jj,kk),
     . debugmsg
                  goto 64
               endif
               if (abs(abs(msq(jj,kk)/msqs(jj,kk))-1d0) .gt. 0.03d0)then
                  debugmsg=' FAILED singularity uncancelled'
               else
                  debugmsg='   OK   cancelled'
               endif
               write(*,63) jj,kk,msq(jj,kk)/msqs(jj,kk),
     .                    msq(jj,kk),msqs(jj,kk),debugmsg
           endif
        
   64      continue
           enddo
           enddo

        pause
        
        enddo
        enddo 
        endif

      elseif (case .eq. 'Z_2jet') then
        call compare(p)
        pause
        if (1 .eq. 2) then
        do j=1,5
        do k=1,9
        call coll2(p,k,j)
        write(*,*) 'Point ',j,':'
c        call dotem(7,p,s)
c        call compare(p)
c        call qqb_w2jet_soft(p,msqs)  
c        write(*,*) 'Expecting --->',nqcdjets,' jets'
c        write(*,*) 'Clustering -->',jets,' jets'
        call qqb_z2jet_g(p,msq)  
        call qqb_z2jet_gs(p,msqc) 
        do jj=-nf,nf
        do kk=-nf,nf
          msqs(jj,kk)=0d0
        enddo
        enddo
        do nd=1,ndmax
          do jj=1,7
          do kk=1,4
            q(jj,kk)=ptilde(nd,jj,kk)
          enddo
          enddo
          call dotem(6,q,s)
          call genclust2(q,rcut,jets,pjet,jetlabel)
          if (jets .eq. nqcdjets) then
            do jj=-nf,nf
            do kk=-nf,nf
              msqs(jj,kk)=msqs(jj,kk)+msqc(nd,jj,kk)
            enddo
            enddo
          else
          endif
        enddo
        do jj=-nf,nf
        do kk=-nf,nf
        if (msq(jj,kk) .ne. 0d0) then
c          write(*,*) 'real  ',jj,kk,msq(jj,kk)
c          write(*,*) 'dips  ',jj,kk,msqs(jj,kk)
          write(*,*) 'ratio ',jj,kk,msq(jj,kk)/msqs(jj,kk),msq(jj,kk)
        endif
        enddo
        enddo
        pause
        enddo
        enddo
        endif 
      elseif (case .eq. 'HWW_4l') then
        call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)     
      elseif (case .eq. 'Hbbbar') then
        call qqb_hbbbar_g(p,msq)      
        call qqb_hbbbar_gs(p,msqc)     
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)     
      elseif (case .eq. 'W_only') then
        call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (case .eq. 'Z_only') then
        call qqb_z_g(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (case .eq. 'W_1jet') then
        call qqb_w2jet(p,msq)      
        call qqb_w1jet_gs(p,msqc)  
        if (1 .eq. 2) then
        do i=1,1
        do kx=2,6
        call coll6(p,kx,i)
        call dotem(6,p,s)
        write(6,*)
        if (abs(s(1,6)) .lt. 1d0) write(6,*) 's(1,6)',s(1,6)
        if (abs(s(2,6)) .lt. 1d0) write(6,*) 's(2,6)',s(2,6)
        if (abs(s(1,5)) .lt. 1d0) write(6,*) 's(1,5)',s(1,5)
        if (abs(s(2,5)) .lt. 1d0) write(6,*) 's(2,5)',s(2,5)
        if (abs(s(5,6)) .lt. 1d0) write(6,*) 's(5,6)',s(5,6)
        call qqb_w2jet(p,msq)      
        call qqb_w1jet_gs(p,msqc)  
         do jj=-nf,nf
         do kk=-nf,nf
           msqs(jj,kk)=0d0
         enddo
         enddo
         do nd=1,ndmax
           do jj=1,6
           do kk=1,4
             q(jj,kk)=ptilde(nd,jj,kk)
           enddo
           enddo
           call dotem(6,q,s)
         if (i.le.3) 
     .   test1=((abs(s(1,5)) .gt. 10d0) .and. (abs(s(2,5)) .gt. 10d0))
         if (i.gt.3) 
     .   test1=((abs(s(1,6)) .gt. 10d0) .and. (abs(s(2,6)) .gt. 10d0))
         if (test1) then
             do jj=-nf,nf
             do kk=-nf,nf
               msqs(jj,kk)=msqs(jj,kk)+msqc(nd,jj,kk)
c            if ((jj.eq.2).and.(kk.eq.1)) write(*,*) nd,msqc(nd,jj,kk)
             enddo
             enddo
         else
         endif
         enddo

         do jj=-nf,nf
         do kk=-nf,nf
         if (abs(msq(jj,kk)/msqs(jj,kk))-1 .gt. 0.01d0) then
           write(*,*) jj,kk,msq(jj,kk)/msqs(jj,kk),msq(jj,kk)
c           write(*,*) 'real  ',jj,kk,msq(jj,kk)
c           write(*,*) 'dips ',jj,kk,msqs(jj,kk)
c           write(*,*) 'ratio ',jj,kk,msq(jj,kk)/msqs(jj,kk)
         endif
         enddo
         enddo

         enddo
         pause
         enddo
         
         endif
      elseif (case .eq. 'Z_1jet') then
        call qqb_z1jet_gs(p,msqc)  
        call qqb_z2jet(p,msq)      
      elseif ((case .eq. 'tt_bbl') .or .(case .eq. 'tt_bbl')) then
        write(6,*) 'No real correction to t_bbar yet'
        stop
      elseif (case .eq. 't_bbar') then
        call qqb_tbb(p,msq)
      elseif (case .eq. 'WZbbar') then
         call qqb_wz_g(p,msq)      
         call qqb_wz_gs(p,msqc)      
      elseif (case .eq. 'WWqqbr') then
         call qqb_ww_g(p,msq)      
         call qqb_ww_gs(p,msqc)      
      elseif (case .eq. 'ZZlept') then
         call qqb_zz_g(p,msq)      
         call qqb_zz_gs(p,msqc)      
      endif
      
      do j=0,ndmax
      xmsq(j)=0d0
      enddo
      
      do j=-nf,nf
      do k=-nf,nf

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (ggexcl) then
      if ((j.eq.0) .and. (k.eq.0)) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (realonly) then 
        xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        do nd=1,ndmax
        xmsq(nd)=0d0
        enddo
      elseif (virtonly) then
         xmsq(0)=0d0
         do nd=1,ndmax
           xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
         enddo
      else
         xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
         do nd=1,ndmax
           xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
         enddo
      endif
 20   continue

      enddo
      enddo

      realint=0d0
      xint=0d0

c --- TESTING THIS
      if (xmsq(0) .eq. 0d0) goto 999
c --- TESTING THIS

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat

        if (nd .eq. 0) then
c---call clustering for event
c--- cluster partons (nqcdstart) to (nqcdstart+nqcdjets)
c--- if nqcdjets=0, no clustering is performed and pjet=p     
        if (clustering .eqv. .false.) then
          do j=1,7
          do k=1,4
          pjet(j,k)=p(j,k)
          enddo
          enddo
          jets=3
        else
          do j=1,7
          do k=1,4
          q(j,k)=p(j,k)
          enddo
          enddo
          call genclust2(q,rcut,jets,pjet,jetlabel)
c          write(*,*) 'should have ',nqcdjets
c          write(*,*) 'actually have ',jets
c          do j=1,jets
c            write(*,*) j,' ',jetlabel(j),' ',
c     .         ((pt(4+j,pjet).lt.15d0) .or. (ayrap(4+j,pjet).gt.2.5d0))
c            if (((jetlabel(j).eq.'pj') .or. (jetlabel(j).eq.'pp')) .and.
c     .         ((pt(4+j,pjet).lt.15d0) .or. (ayrap(4+j,pjet).gt.2.5d0)))
c     .        jets=jets-1
c          enddo  
c          write(*,*) 'now have ',jets
c          pause
          if (jets .ne. nqcdjets) then
c--- DEBUG
          njetzero=njetzero+1
          ncutzero=ncutzero-1
          if (nqcdjets .eq. 0) then
            do j=0,ndmax
            xmsq(j)=0d0
            enddo       
            goto 998
          else
            failed=.true.
            goto 996
          endif
c--- two lines above are added
c            njetzero=njetzero+1
c            ntotzero=ntotzero+1
c            do j=0,ndmax
c            xmsq(j)=0d0
c            enddo       
c            goto 998
c--- DEBUG
          endif
        endif
        else
c---call clustering for counter-event
        if (clustering .eqv. .false.) then
          do j=1,7
          do k=1,4
           pjet(j,k)=ptilde(nd,j,k)
          enddo
          enddo
          jets=2
        else
          do j=1,7
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo
c--- cluster partons (nqcdstart) to (nqcdstart+nqcdjets-1)
c--- if nqcdjets=0, no clustering is performed and pjet=p      
          call genclust2(q,rcut,jets,pjet,jetlabel)
        endif
        endif

        call dotem(nvec,pjet,s)

c---check to see if (counter-)event passes cuts
        failed=.false.
        if (clustering .and. (jets .ne. nqcdjets)) then
          failed=.true.
          goto 996
        endif
        if (makecuts) then
          if (bbproc) then 
            failed=cuts(pjet)
          else
            if (madejetcuts(q,pjet,jets) .eqv. .false.) failed=.true.
          endif
        endif
 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997         
        endif
c---if it does, add to total
        xint=xint+xmsq(nd)
c---if we're binning, add to histo too
        if (bin) then
          val=xmsq(nd)*wgt/dfloat(itmx)
          if (nd .eq. 0) then
            call nplotter(vector,s,pjet,val,0)
          else
            call nplotter(vector,s,pjet,val,1)
          endif
        endif
c---otherwise, skip contribution
 997    continue
      enddo

      call dotem(nvec,p,s)

 998  continue

      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+xint**2*wgt/dfloat(itmx)

      if (abs(realint) .gt. 1d99) then
        write(*,*) realint
        call dotem(nvec,p,s)
        write(*,*) 's12',s(1,2)
        write(*,*) 's34',s(3,4)
        write(*,*) 's15',s(1,5)
        write(*,*) 's25',s(2,5)
        write(*,*) 's16',s(1,6)
        write(*,*) 's26',s(2,6)
        write(*,*) 's17',s(1,7)
        write(*,*) 's27',s(2,7)
        write(*,*) 's56',s(5,6)
        write(*,*) 's57',s(5,7)
        write(*,*) 's67',s(6,7)
        write(*,*) 'pt5',dsqrt(s(1,5)*s(2,5)/s(1,2))
        write(*,*) 'pt6',dsqrt(s(1,6)*s(2,6)/s(1,2))
        write(*,*) 'pt7',dsqrt(s(1,7)*s(2,7)/s(1,2))
        write(*,*) 'y5',0.5d0*dabs(log((p(5,4)+p(5,3))/(p(5,4)-p(5,3))))
        write(*,*) 'y6',0.5d0*dabs(log((p(6,4)+p(6,3))/(p(6,4)-p(6,3))))
        write(*,*) 'y7',0.5d0*dabs(log((p(7,4)+p(7,3))/(p(7,4)-p(7,3))))
        do nd=1,ndmax
c          write(*,*) 'nd',nd,xmsq(nd)
        enddo
        write(*,*) 'dipoles',realint-xmsq(0)
        write(*,*) 'real   ',xmsq(0)
        do j=1,7
        do k=1,4
c          write(*,*) '     p(',j,',',k,')=',p(j,k),'d0'
        enddo
        enddo
c        write(*,*) '     pswt=',pswt,'d0'
        pause
      endif
      
      return

 999  realint=0d0
      ntotzero=ntotzero+1
 
   63 format(1x,2i3,f10.6,2e14.6,1x,a32)
   65 format(1x,2i3,a10,2e14.6,1x,a32)
      return
      end
















