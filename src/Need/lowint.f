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
      include 'efficiency.f'
      include 'facscale.f'
      include 'maxwt.f'
      integer ih1,ih2,j,k,jets,nproc,nvec,sgnj,sgnk
      double precision r(mxdim),W,sqrts,xmsq,val,temp
      double precision fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4)
      double precision pswt,scalestart,pdfscale
      double precision wgt,msq(-nf:nf,-nf:nf),m3,m4,m5,xmsqjk
      double precision msqa(-nf:nf,-nf:nf),n(4)
      double precision xx(2),flux,vol,vol_mass,vol3_mass,
     . taumin,BrnRat,rcut
      double precision xmsq_bypart(-1:1,-1:1),lord_bypart(-1:1,-1:1)
      integer nqcdjets,nqcdstart
      common/parts/jets,jetlabel
      common/nqcdjets/nqcdjets,nqcdstart
      logical bin,makecuts,gencuts,first
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
      common/bypart/lord_bypart
      data p/48*0d0/
      data first/.true./
      save first,scalestart
      if (first) then
         first=.false.
         scalestart=scale
      endif

      ntotshot=ntotshot+1
      lowint=0d0

      W=sqrts**2


      if ((case .eq. 'tt_bbh')  .or. (case .eq. 'tt_bbl')
     .      .or. (case .eq. 'tautau') .or. (case .eq. 'vlchk6')) then
          npart=6
          call gen6(r,p,pswt,*999)
      elseif ((case .eq. 'ttbdkl')  .or. (case .eq. 'ttbdkh')) then
          npart=6
          call gen6_rap(r,p,pswt,*999)                                          
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
     .   .or. (case .eq. 'W_3jet') 
     .   .or. (case .eq. 'Z_3jet') 
     .   .or. (case .eq. 'vlchk5') 
     .        ) then
          npart=5
          call gen_njets(r,3,p,pswt,*999)      
c          call gen5(r,p,pswt,*999)
c          call gen5a(r,p,pswt,*999)
      elseif ((case .eq. 'vlchk3') 
     .   .or. (case .eq. 'httjet')
     .   .or. (case .eq. 'attjet')
     .   .or. (case .eq. 'H_1jet')
     .   .or. (case .eq. 'Z_1jet')
     .   .or. (case .eq. 'Wgamma')
     .   .or. (case .eq. 'W_1jet')) then
          if (case .eq. 'vlchk3') then
            wsqmin=0d0
            wsqmax=sqrts**2
          endif
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
      elseif (case .eq. 'vlchm3') then
          taumin=(2d0*mt/sqrts)**2
          npart=3
          m3=mt
          m4=mt
          m5=0d0
          call gen3m_rap(r,p,m3,m4,m5,pswt,*999)
      elseif (case .eq. 'threeb') then
          npart=3
          m3=zip
          m4=zip
          m5=zip
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
     .   .or. (case .eq. 'Hbbbar')
     .   .or. (case .eq. 'vlchk2')
     .        ) then
          if (case .eq. 'vlchk2') then
            wsqmin=0d0
            wsqmax=sqrts**2
          endif
          npart=2
          if (new_pspace) then
          call gen2a(r,p,pswt,*999)
          else
          call gen2(r,p,pswt,*999)
          endif
      elseif ((case .eq. 'tt_tot')
     .   .or. (case .eq. 'bb_tot')
     .   .or. (case .eq. 'cc_tot')) then
          npart=2
          call gen2m(r,p,pswt,*999)
      elseif (
     .        (case .eq. 'W_2jet')
     .   .or. (case .eq. 'Z_2jet')
     .        ) then
          npart=4
          call gen_njets(r,2,p,pswt,*999)      
c          call gen4(r,p,pswt,*999)      
      else
          if ((case .eq. 'vlchk4') .or. (case .eq. 'vlchkm')) then
            wsqmin=0d0
            wsqmax=sqrts**2
          endif
          npart=4
          call gen4(r,p,pswt,*999)      
      endif
      nvec=npart+2
      call dotem(nvec,p,s)

c---impose cuts on final state
      if (case .ne. 'vlchk6' .and. case .ne. 'tautau') then
        call masscuts(s,*999)
      endif
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)                                                 

      if (scalestart .lt. 0d0) call scaleset(scalestart,p)
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx(1),xx(2)      

      if (case .eq. 'H_1jet') then
        pdfscale=facscale
      else
        pdfscale=scale
      endif   
      
      if (case(1:4) .ne. 'vlch') then      
      call fdist(pdlabel,ih1,xx(1),pdfscale,fx1)
      call fdist(pdlabel,ih2,xx(2),pdfscale,fx2)
      endif

      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)

      if     (case .eq. 'Wbbbar') then
      call qqb_Wbb(p,msq)
      elseif (case .eq. 'Wbbmas') then
      call qqb_Wbbm(p,msq)
c      call compare_wbb_mad(p)
      elseif (case .eq. 'Zbbbar') then
      call qqb_Zbb(p,msq)
      elseif (case .eq. 'Zbbmas') then
      call qqb_Zbbm(p,msq)
c      call compare_zbb_mad(p)
      elseif (case .eq. 'Z_only') then
      call qqb_Z(p,msq)
      elseif (case .eq. 'H_1jet') then
      call qqb_Hg(p,msq)
      elseif (case .eq. 'W_1jet') then
      call qqb_w_g(p,msq)
      elseif (case .eq. 'Wgamma') then
      call qqb_wgam(p,msq)
      elseif (case .eq. 'Z_1jet') then
      call qqb_z1jet(p,msq)
      elseif (case .eq. 'W_2jet') then
      call qqb_w2jet(p,msq)
      elseif (case .eq. 'W_3jet') then
      call qqb_w2jet_g(p,msq)
      elseif (case .eq. 'Z_2jet') then
      call qqb_z2jet(p,msq)
      elseif (case .eq. 'Z_3jet') then
      call qqb_z2jet_g(p,msq)
      elseif (case .eq. 'W_only') then
      call qqb_W(p,msq)
      elseif (case .eq. 't_bbar') then
      call qqb_tbb(p,msq)
      elseif (case .eq. 'WmZbbr') then
      call qqb_wz(p,msq)
      elseif (case .eq. 'WZbbar') then
c      call comparewz(p)
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
      elseif (case .eq. 'Hbbbar') then
      call qqb_Hbbbar(p,msq)
      elseif (case .eq. 'qg_tbb') then
      call qg_tbb(p,msq)
c      elseif ((case .eq. 'tt_bbh') .or. (case .eq. 'tt_bbl')) then
c      call qqb_ttb(p,msq)
c      elseif ((case .eq. 'ttbdkl') .or. (case .eq. 'ttbdkh')) then
c      call qqb_QQbdk(p,msq)                                                     
c      elseif (case .eq. 'tt_tot') then
c      call qqb_QQb(p,msq)
c      elseif (case .eq. 'bb_tot') then
c      call qqb_QQb(p,msq)
c      elseif (case .eq. 'cc_tot') then
c      call qqb_QQb(p,msq)
      elseif (case .eq. 'qq_ttg') then
      call qqb_ttb_g(p,msq)
      elseif (case .eq. 'tautau') then
      call qqb_tautau(p,msq)
      elseif (case .eq. 'tottth') then
      call qqb_tottth(p,msq)
      elseif (case .eq. 'threeb') then
      call gq_qqqb(p,msq)
      elseif (case .eq. 'qq_tth') then
      call qqb_tth(p,msq)
      elseif (case .eq. 'qq_ttz') then
      call qqb_ttz(p,msq)
      elseif (case .eq. 'httjet') then
      call qqb_higgs(p,msq)
      elseif (case .eq. 'attjet') then
      call qqb_higgs_odd(p,msq)
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

      elseif (case .eq. 'vlchkm') then
      taumin=0.0001d0
      bbsqmax=W
      bbsqmin=0d0
      call qqb_vol(p,msq)
      flux=one/vol_mass(mb,W)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=4d0*xx(1)
      fx2(-1)=2d0/xx(2)
      
      elseif (case .eq. 'vlchm3') then
      taumin=(2d0*mt/sqrts)**2
      call qqb_vol(p,msq)
      flux=one/vol3_mass(mt,W)
      do j=-nf,nf
      fx1(j)=0d0
      fx2(j)=0d0
      enddo
      fx1(2)=2d0
      fx2(-1)=2d0
      endif
      
      xmsq=0d0
      do j=-1,1
      do k=-1,1
      xmsq_bypart(j,k)=0d0
      enddo
      enddo
      
      do j=-nf,nf
      do k=-nf,nf    
c      do j=-5,5,5
c      do k=-5,5,5
c      if (abs(j+k) .ne. 5) goto 20   

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (ggexcl) then
      if ((j.eq.0) .and. (k.eq.0)) goto 20
      endif
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      xmsqjk=fx1(j)*fx2(k)*msq(j,k)
      xmsq=xmsq+xmsqjk

      if     (j .gt. 0) then
        sgnj=+1
      elseif (j .lt. 0) then
        sgnj=-1
      else
        sgnj=0
      endif
      if     (k .gt. 0) then
        sgnk=+1
      elseif (k .lt. 0) then
        sgnk=-1
      else
        sgnk=0
      endif

      xmsq_bypart(sgnj,sgnk)=xmsq_bypart(sgnj,sgnk)+xmsqjk
      
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
        jets=nqcdjets
      else
        call genclust2(p,rcut,jets,pjet,jetlabel)
        if((nproc .eq. 152) .or. (nproc .eq. 161))then
          if (jets .ne. 2) then
            njetzero=njetzero+1
            goto 999
          endif
        else
          if ((jets .ne. nqcdjets) .and. (nqcdjets .gt. 0)
     .        .and. clustering) then
            njetzero=njetzero+1
            goto 999
          endif
        endif       
      endif
      
      call dotem(nvec,pjet,s)

c--- Apply the event cuts
      if (makecuts) then
        if (gencuts(p,pjet,jets)) then
          ncutzero=ncutzero+1
          goto 999
         endif
      endif
      
      lowint=flux*pswt*xmsq/BrnRat
      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=lord_bypart(j,k)+
     .       flux*pswt*xmsq_bypart(j,k)/BrnRat
      enddo
      enddo

      val=lowint*wgt
c--- update the maximum weight so far, if necessary
      if (val .gt. wtmax) then
        wtmax=val
      endif

      if (bin) then
        val=val/dfloat(itmx)
        call nplotter(r,s,pjet,val,0)
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


