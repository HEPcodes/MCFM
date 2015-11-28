      double precision function lowint(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'dprodx.f'
      include 'scale.f'
      include 'clustering.f'
      include 'noglue.f'
      include 'process.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer ih1,ih2,j,k,jets,nproc,nvec
      double precision r(mxdim),W,sqrts,xmsq,val
      double precision fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),pnew(mxpart,4)
      double precision pswt,wbb,zz,alphas,pttwo
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5
      double precision xx(2),flux,vol,taumin,BrnRat,ptjet,rcut,amz
      integer nqcdjets,nqcdstart
      common/parts/jets,jetlabel
      common/nqcdjets/nqcdjets,nqcdstart
      logical bin,makecuts,madejetcuts,cuts,bbproc
      character pdlabel*7,jetlabel(mxpart)*2
      common/density/ih1,ih2
      common/energy/sqrts
      common/pdlabel/pdlabel
      common/bin/bin
      common/x1x2/xx
      common/taumin/taumin
      common/makecuts/makecuts
      common/BrnRat/BrnRat
      common/nproc/nproc
      common/rcut/rcut
      common/couple/amz
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                     
     . ,SUUB_VEVEBBB,SDDB_VEVEBBB,SUUB_EMEPBBB,SDDB_EMEPBBB
     . ,SUUB_VEVE,SDDB_VEVE,SUUB_EMEP,SDDB_EMEP
      lowint=0d0

      W=sqrts**2


      if ((case .eq. 'tt_bbh')  .or. (case .eq. 'tt_bbl')
     .      .or. (case .eq. 'tautau') .or. (case .eq. 'vlchk6')) then
          npart=6
          call gen6(r,p,pswt,*999)
      elseif (
     .        (case .eq. 'qq_tth') 
     .   .or. (case .eq. 'qq_ttz') 
     .   .or. (case .eq. 'vlchk8')
     .        ) then
          npart=8
          call gen8(r,p,pswt,*999)
      elseif (
     .        (case .eq. 'qq_ttg') 
     .   .or. (case .eq. 'hlljet') 
     .        ) then
          npart=7
          call gen7(r,p,pswt,*999)
      elseif (
     .        (case .eq. 'qg_tbb') 
     .   .or. (case .eq. 'vlchk5') 
     .        ) then
          npart=5
          call gen5(r,p,pswt,*999)
c          call gen5a(r,p,pswt,*999)
      elseif ((case .eq. 'vlchk3') 
     .   .or. (case .eq. 'httjet')
     .   .or. (case .eq. 'Z_1jet')
     .   .or. (case .eq. 'W_1jet')) then
          npart=3
          if (new_pspace) then
          call gen3a(r,p,pswt,*999)
          else
          call gen3b(r,p,pswt,*999)
          endif
      elseif (case .eq. 'tottth') then
          npart=3
          m3=mt
          m4=mt
          m5=hmass
          call gen3m(r,p,m3,m4,m5,pswt,*999)
      elseif (case .eq. 'totttz') then
          npart=3
          m3=mt
          m4=mt
          m5=zmass
          call gen3m(r,p,m3,m4,m5,pswt,*999)
      elseif (
     .        (case .eq. 'W_only')
     .   .or. (case .eq. 'Z_only')
     .   .or. (case .eq. 'vlchk2')
     .        ) then
c          wsqmin=0d0
c          wsqmax=sqrts**2
          npart=2
          if (new_pspace) then
          call gen2a(r,p,pswt,*999)
          else
          call gen2(r,p,pswt,*999)
          endif
      else
          npart=4
          call gen4(r,p,pswt,*999)      
      endif
      nvec=npart+2
      call dotem(nvec,p,s)

      if (case .ne. 'vlchk6' .and. case .ne. 'tautau') then
        call masscuts(s,*999)
      endif

c--debug 
c      shat=s(1,2)
c      scale=sqrt(shat)
c      if (scale .gt. 1000d0) scale=1000d0
c--debug

      if   ((case .eq. 'httjet')
     . .or. (case .eq. 'hlljet')) then
c--- scale for processes 200 and 201 only
        scale=dsqrt(hmass**2+pttwo(3,4,p)**2)
        if (scale .gt. 1000d0) scale=1000d0
        as=alphas(scale,amz,2)
        ason2pi=as/twopi
        gsq=fourpi*as
      endif
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

      call fdist(pdlabel,ih1,xx(1),scale,fx1)
      call fdist(pdlabel,ih2,xx(2),scale,fx2)

 
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)

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
      call qqb_Wbb(p,msq)
      elseif (case .eq. 'Zbbbar') then
c      call qqb_Zbb(p,msq)
c      write(*,*) 'OLD u ',msq(2,-2)
c      write(*,*) 'OLD ub',msq(-2,2)
c      write(*,*) 'OLD d ',msq(1,-1)
c      write(*,*) 'OLD db',msq(-1,1)
      call qqb_Zbb(p,msq)
c      write(*,*) 'NEW u ',msq(2,-2)
c      write(*,*) 'NEW ub',msq(-2,2)
c      write(*,*) 'NEW d ',msq(1,-1)
c      write(*,*) 'NEW db',msq(-1,1)
c--- call MADGRAPH routines     
c      do k=1,4
c        if (k.lt.4) then
c          j=k
c        else
c          j=0
c        endif 
c        p1(j)=-p(1,k)
c        p2(j)=-p(2,k)
c        p3(j)=p(3,k)
c        p4(j)=p(4,k)
c        p5(j)=p(5,k)
c        p6(j)=p(6,k)
c      enddo            
c      call initialize
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      enddo
c      enddo
c      do j=-nf,nf
c      if     ((j.gt.0) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,-j)=gsq**2*esq**2/(fourpi/128d0)**2*
c     .             suub_emepbbb(p1,p2,p3,p4,p5,p6)
c      elseif ((j.lt.0) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,-j)=gsq**2*esq**2/(fourpi/128d0)**2*
c     .             suub_emepbbb(p2,p1,p3,p4,p5,p6)
c      elseif ((j.gt.0) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,-j)=gsq**2*esq**2/(fourpi/128d0)**2*
c     .             sddb_emepbbb(p1,p2,p3,p4,p5,p6)
c      elseif ((j.lt.0) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,-j)=gsq**2*esq**2/(fourpi/128d0)**2*
c     .             sddb_emepbbb(p2,p1,p3,p4,p5,p6)
c      endif
c      enddo
c      write(*,*) 'MAD u ',msq(2,-2)
c      write(*,*) 'MAD ub',msq(-2,2)
c      write(*,*) 'MAD d ',msq(1,-1)
c      write(*,*) 'MAD db',msq(-1,1)
c      pause
c      call qqb_Zbb_alt(p,msq)
      elseif (case .eq. 'Z_only') then
      call qqb_Z(p,msq)
c      write(*,*) 'MCFM u ',msq(2,-2)
c      write(*,*) 'MCFM ub',msq(-2,2)
c      write(*,*) 'MCFM d ',msq(1,-1)
c      write(*,*) 'MCFM db',msq(-1,1)
c--- call MADGRAPH routines     
c      do k=1,4
c        if (k.lt.4) then
c          j=k
c        else
c          j=0
c        endif 
c        p1(j)=-p(1,k)
c        p2(j)=-p(2,k)
c        p3(j)=p(3,k)
c        p4(j)=p(4,k)
c      enddo            
c      call initialize
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      enddo
c      enddo
c      do j=-nf,nf
c      if     ((j.gt.0) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,-j)=3d0*esq**2/(fourpi/128d0)**2*
c     .             suub_veve(p1,p2,p3,p4)
c      elseif ((j.lt.0) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,-j)=3d0*esq**2/(fourpi/128d0)**2*
c     .             suub_veve(p2,p1,p3,p4)
c      elseif ((j.gt.0) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,-j)=3d0*esq**2/(fourpi/128d0)**2*
c     .             sddb_veve(p1,p2,p3,p4)
c      elseif ((j.lt.0) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,-j)=3d0*esq**2/(fourpi/128d0)**2*
c     .             sddb_veve(p2,p1,p3,p4)
c      endif
c      enddo
c      write(*,*) 'MAD u ',msq(2,-2)
c      write(*,*) 'MAD ub',msq(-2,2)
c      write(*,*) 'MAD d ',msq(1,-1)
c      write(*,*) 'MAD db',msq(-1,1)
c      pause
      elseif (case .eq. 'W_1jet') then
      call qqb_w_g(p,msq)
      elseif (case .eq. 'Z_1jet') then
      call qqb_z_g(p,msq)
      elseif (case .eq. 'W_2jet') then
      call qqb_w2jet(p,msq)
      elseif (case .eq. 'Z_2jet') then
      call qqb_z2jet(p,msq)
      elseif (case .eq. 'W_only') then
      call qqb_W(p,msq)
      elseif (case .eq. 't_bbar') then
      call qqb_tbb(p,msq)
      elseif (case .eq. 'WmZbbr') then
      call qqb_wz(p,msq)
      elseif (case .eq. 'WZbbar') then
      call qqb_wz(p,msq)
      elseif (case .eq. 'ZZlept') then
      call qqb_ZZ(p,msq)
      elseif (case .eq. 'WWqqbr') then
      call qqb_WW(p,msq)
      elseif (case .eq. 'WHbbar') then
      call qqb_WH(p,msq)
      elseif (case .eq. 'ZHbbar') then
      call qqb_ZH(p,msq)
      elseif (case .eq. 'HWW_4l') then
      call qqb_HWW(p,msq)
      elseif (case .eq. 'HZZ_4l') then
      call qqb_HZZ(p,msq)
      elseif (case .eq. 'qg_tbb') then
      call qg_tbb(p,msq)
      elseif ((case .eq. 'tt_bbh') .or. (case .eq. 'tt_bbl')) then
      call qqb_ttb(p,msq)
      elseif (case .eq. 'qq_ttg') then
      call qqb_ttb_g(p,msq)
      elseif (case .eq. 'tautau') then
      call qqb_tautau(p,msq)
      elseif (case .eq. 'tottth') then
      call qqb_tottth(p,msq)
      elseif (case .eq. 'qq_tth') then
      call qqb_tth(p,msq)
      elseif (case .eq. 'qq_ttz') then
      call qqb_ttz(p,msq)
      elseif (case .eq. 'httjet') then
      call qqb_higgs(p,msq)
      elseif (case .eq. 'hlljet') then
      call Htautau(p,msq)
      elseif (case .eq. 'qqttbb') then
c      call qqb_ttbbbb(p,msq)
      write(*,*) 'MADGRAPH files need to be included'
      stop
      elseif (case .eq. 'qqttgg') then
c      call qqb_ttbgg(p,msq)
      write(*,*) 'MADGRAPH files need to be included'
      stop
      elseif (case .eq. 'vlchk2') then
      call qqb_vol(p,msq)
      flux=one/vol(W,2)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=1d0
      fx2(-1)=1d0
      elseif (case .eq. 'vlchk3') then
      call qqb_vol(p,msq)
      flux=one/vol(W,3)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=2d0
      fx2(-1)=2d0
      elseif (case .eq. 'vlchk4') then
      taumin=0.0001d0
      bbsqmax=W
      bbsqmin=0d0
      call qqb_vol(p,msq)
      flux=one/vol(W,4)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=4d0*xx(1)
      fx2(-1)=2d0/xx(2)
      elseif (case .eq. 'vlchk5') then
      taumin=0.0001d0
      bbsqmax=W
      bbsqmin=0d0
      call qqb_vol(p,msq)
      flux=one/vol(W,5)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=4d0
      fx2(-1)=4d0
      elseif (case .eq. 'vlchk6') then
      call qqb_vol(p,msq)
      flux=one/vol(W,6)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=6d0*xx(1)
      fx2(-1)=4d0/xx(2)

      elseif (case .eq. 'vlchk8') then
      call qqb_vol(p,msq)
      flux=one/vol(W,8)
      bbsqmax=W
      bbsqmin=0d0
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=6d0/xx(1)
      fx2(-1)=6d0/xx(2)
      endif
      
      xmsq=0d0
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

      xmsq=xmsq+fx1(j)*fx2(k)*msq(j,k)
    
 20   continue
      enddo
      enddo

c--- cluster partons (nqcdstart) to (nqcdstart+nqcdjets-1)
c--- if nqcdjets=0, no clustering is performed and pjet=p  
      if (clustering .eqv. .false.) then
        do j=1,mxpart
        do k=1,4
          pjet(j,k)=p(j,k)
        enddo
        enddo
      else
        call genclust2(p,rcut,jets,pjet,jetlabel)
        if((nproc .eq. 152) .or. (nproc .eq. 161))then
          if (jets .ne. 2) goto 999
        else
          if ((jets .ne. nqcdjets) .and. (nqcdjets .gt. 0)
     .        .and. (clustering)) then
            goto 999
          endif
        endif       
      endif
      
      call dotem(nvec,pjet,s)

c--- make some cuts ....
      if (makecuts) then
        if (bbproc) then
          if (cuts(pjet)) goto 999
        else
          if (madejetcuts(p,pjet,jets) .eqv. .false.) goto 999
        endif
      endif
      
      lowint=flux*pswt*xmsq/BrnRat

      if (bin) then
      val=lowint*wgt/dfloat(itmx)
      call nplotter(r,s,pjet,val,0)
      endif

 999  continue
      
      return
      end


