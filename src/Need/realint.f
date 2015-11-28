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
      include 'qcdcouple.f'
      include 'efficiency.f'
      include 'bbproc.f'
      include 'facscale.f'
      include 'maxwt.f'
      include 'process.f'
      include 'jetlabel.f'
      include 'pdlabel.f'
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec     
      double precision vector(mxdim),W,val,xint,reweight,n(4)
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,scalestart,pdfscale
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqa(-nf:nf,-nf:nf),msqs(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4),rcut,pttwo,pt
      double precision m3,m4,m5
      double precision s19,s29,s3459,s6789
      integer nqcdjets,nqcdstart
      integer n2,n3
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/nqcdjets/nqcdjets,nqcdstart
      common/xreal/xreal,xreal2
      logical bin,makecuts,first,failed,gencuts
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs,
     . qqb_w2jet,qqb_w1jet_gs,qqb_z2jet,qqb_z1jet_gs,qqb_Hg_g,qqb_Hg_gs,
     . qqb_hww_g,qqb_hww_gs,qqb_zbb_g,qqb_zbb_gs,qqb_wbb_g,qqb_wbb_gs,
     . qqb_w_g,qqb_w_gs,qqb_z1jet,qqb_z_gs,qqb_ww_g,qqb_ww_gs,
     . qqb_wz_g,qqb_wz_gs,qqb_zz_g,qqb_zz_gs,qqb_wgam_g,qqb_wgam_gs,
     . qqb_zgam_g,qqb_zgam_gs,qqb_dirgam_g,qqb_dirgam_gs,
     . qq_Hqq_g,qq_Hqq_gs,WW_Hqq_g,WW_Hqq_gs,ZZ_Hqq_g,ZZ_Hqq_gs,
     . gg_Hg,gg_H_gs,gg_Hgg,gg_Hg_gs
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/makecuts/makecuts
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/rcut/rcut
      data p/48*0d0/
      data first/.true./
      save first,scalestart
      if (first) then
         first=.false.
         scalestart=scale
      endif
      
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
      if (
     .       (case .eq. 'W_only')
     .  .or. (case .eq. 'Z_only')
     .  .or. (case .eq. 'ggfus0')
     .  .or. (case .eq. 'Hbbbar')
     . ) then
          npart=3
          if (new_pspace) then
          call gen3a(vector,p,pswt,*999)      
          else
          call gen3(vector,p,pswt,*999)      
          endif
      elseif (
     .       (case .eq. 'twojet')
     .  .or. (case .eq. 'dirgam')
     . ) then
          npart=3
          call gen3jet(vector,p,pswt,*999)      
      elseif ((case .eq. 'tt_tot') .or. (case .eq. 'cc_tot')
     .   .or. (case .eq. 'bb_tot')) then
          npart=3
          m3=mass2
          m4=mass2
          m5=0d0
c          call gen3m_rap(vector,p,m3,m4,pswt,*999)
          call gen3m(vector,p,m3,m4,m5,pswt,*999)
      elseif ((case .eq. 'W_1jet')
     .   .or. (case .eq. 'Z_1jet')
     .   .or. (case .eq. 'H_1jet')
     .   .or. (case .eq. 'ggfus1')) then
          npart=4
          if (new_pspace) then
          call gen4a(vector,p,pswt,*999)      
          else
          call gen_njets(vector,2,p,pswt,*999)
c          call gen4(vector,p,pswt,*999)
          endif
      elseif ((case .eq. 'W_2jet')
     .   .or. (case .eq. 'Z_2jet')
     .   .or. (case .eq. 'WW_Hqq')
     .   .or. (case .eq. 'ZZ_Hqq')
     .   .or. (case .eq. 'VV_Hqq')
     .   .or. (case .eq. 'Wbbbar')
     .   .or. (case .eq. 'Zbbbar')
     .       ) then
          npart=5
          call gen_njets(vector,3,p,pswt,*999)
      elseif ((case .eq. 'Wgamma') .or. (case .eq. 'Zgamma')) then
          npart=4
          call gen4(vector,p,pswt,*999)
      elseif ((case .eq. 'tt_bbh')  .or. (case .eq. 'tt_bbl')) then
          npart=7
          call gen7(vector,p,pswt,*999)
      elseif ((case .eq. 'ttbdkl')  .or. (case .eq. 'ttbdkh')) then
          npart=7
          call gen7_rap(vector,p,pswt,*999)
      else
          npart=5
          if (new_pspace) then
          call gen5a(vector,p,pswt,*999)
          else
          call gen5(vector,p,pswt,*999)      
          endif
      endif
      nvec=npart+2

      call dotem(nvec,p,s)

c---impose cuts on final state
      call masscuts(s,*999)
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      
      if ((scalestart .lt. 0d0) .or. (case .eq. 'H_1jet'))
     . call scaleset(scalestart,p)
      
c---- generate collinear points that satisy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2
      
      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         realint=0d0
         return
      endif

      if (case .eq. 'H_1jet') then
        pdfscale=facscale
      else
        pdfscale=scale
      endif   
            
      call fdist(ih1,xx1,pdfscale,fx1)
      call fdist(ih2,xx2,pdfscale,fx2)
      
      flux=fbGeV2/(two*xx1*xx2*W)

      if     (case .eq. 'W_only') then
c        call singcheck(qqb_w_g,qqb_w_gs,p)         ! Checked 11/30/01
        call qqb_w_g(p,msq)      
        call qqb_w_gs(p,msqc)     
      elseif (case .eq. 'W_1jet') then
c        call singcheck(qqb_w2jet,qqb_w1jet_gs,p)   ! Checked 11/16/01
        call qqb_w2jet(p,msq)      
        call qqb_w1jet_gs(p,msqc)  
      elseif (case .eq. 'Wgamma') then
c        call singcheck(qqb_wgam_g,qqb_wgam_gs,p)   ! Checked 08/27/02
        call qqb_wgam_g(p,msq)      
        call qqb_wgam_gs(p,msqc)  
      elseif (case .eq. 'Zgamma') then
c        call singcheck(qqb_zgam_g,qqb_zgam_gs,p)
        call qqb_zgam_g(p,msq)      
        call qqb_zgam_gs(p,msqc)  
      elseif (case .eq. 'Wbbbar') then
c        call singcheck(qqb_wbb_g,qqb_wbb_gs,p)     ! Checked 11/30/01
        call qqb_wbb_g(p,msq)      
        call qqb_wbb_gs(p,msqc)      
      elseif (case .eq. 'W_2jet') then
c        call singcheck(qqb_w2jet_g,qqb_w2jet_gs,p) ! Checked 11/16/01
        call qqb_w2jet_g(p,msq)  
        call qqb_w2jet_gs(p,msqc)
      elseif (case .eq. 'Z_only') then
c        call singcheck(qqb_z1jet,qqb_z_gs,p)         ! Checked 11/30/01
        call qqb_z1jet(p,msq)      
        call qqb_z_gs(p,msqc)     
      elseif (case .eq. 'Z_1jet') then
c        call singcheck(qqb_z2jet,qqb_z1jet_gs,p)   ! Checked 11/16/01
        call qqb_z2jet(p,msq)      
        call qqb_z1jet_gs(p,msqc)  
      elseif (case .eq. 'Z_2jet') then
c        call singcheck(qqb_z2jet_g,qqb_z2jet_gs,p) ! Checked 11/16/01
        call qqb_z2jet_g(p,msq)  
        call qqb_z2jet_gs(p,msqc) 
      elseif (case .eq. 'Zbbbar') then
c        call singcheck(qqb_zbb_g,qqb_zbb_gs,p)     ! Checked 11/30/01
        call qqb_zbb_g(p,msq)
        call qqb_zbb_gs(p,msqc) 
      elseif (case .eq. 'WWqqbr') then
c        call singcheck(qqb_ww_g,qqb_ww_gs,p)       ! Checked 11/30/01
        call qqb_ww_g(p,msq)      
        call qqb_ww_gs(p,msqc)      
      elseif (case .eq. 'WZbbar') then
c        call singcheck(qqb_wz_g,qqb_wz_gs,p)       ! Checked 12/05/01
        call qqb_wz_g(p,msq)      
        call qqb_wz_gs(p,msqc)      
      elseif (case .eq. 'ZZlept') then
c        call singcheck(qqb_zz_g,qqb_zz_gs,p)       ! Checked 12/05/01
        call qqb_zz_g(p,msq)      
        call qqb_zz_gs(p,msqc)      
      elseif (case .eq. 'WHbbar') then
        call qqb_wh_g(p,msq)      
        call qqb_wh_gs(p,msqc)     
      elseif (case .eq. 'ZHbbar') then
        call qqb_zh_g(p,msq)      
        call qqb_zh_gs(p,msqc)     
      elseif (case .eq. 'dirgam') then
c        call singcheck(qqb_dirgam_g,qqb_dirgam_gs,p)
c        call compare_dirgam_g_mad(p)     
c        pause
c        call qqb_dirgam_g(p,msq)      
        call qqb_dirgam_g(p,msq)      
        call qqb_dirgam_gs(p,msqc)      
      elseif (case .eq. 'qq_Hqq') then
c        call singcheck(qq_Hqq_g,qq_Hqq_gs,p)
c        call compareqqh_g(p)        ! Done now
        call qq_Hqq_g(p,msq)
        call qq_Hqq_gs(p,msqc)
      elseif (case .eq. 'WW_Hqq') then
c        call singcheck(WW_Hqq_g,WW_Hqq_gs,p)   ! Checked 3/03
        call WW_Hqq_g(p,msq)
        call WW_Hqq_gs(p,msqc)
      elseif (case .eq. 'ZZ_Hqq') then
c        call singcheck(ZZ_Hqq_g,ZZ_Hqq_gs,p)   ! Checked 3/03
        call ZZ_Hqq_g(p,msq)
        call ZZ_Hqq_gs(p,msqc)
      elseif (case .eq. 'VV_Hqq') then
        call VV_Hqq_g(p,msq)
        call VV_Hqq_gs(p,msqc)
       elseif (case .eq. 'ggfus0') then
c         call singcheck(gg_hg,gg_h_gs,p)       ! Checked 28/02/03
         call gg_hg(p,msq)
         call gg_h_gs(p,msqc)
       elseif (case .eq. 'ggfus1') then
c         call singcheck(gg_hgg,gg_hg_gs,p)
         call gg_hgg(p,msq)
         call gg_hg_gs(p,msqc)
      elseif (case .eq. 'HWW_4l') then
c        call singcheck(qqb_hww_g,qqb_hww_gs,p)
        call qqb_hww_g(p,msq)      
        call qqb_hww_gs(p,msqc)     
      elseif (case .eq. 'Hbbbar') then
        call qqb_hbbbar_g(p,msq)      
        call qqb_hbbbar_gs(p,msqc)     
      elseif (case .eq. 'HZZ_4l') then
        call qqb_hzz_g(p,msq)      
        call qqb_hzz_gs(p,msqc)  
      elseif (case .eq. 'H_1jet') then
c        call singcheck(qqb_Hg_g,qqb_Hg_gs,p)       ! Checked 19/02/02
        call qqb_Hg_g(p,msq)  
        call qqb_Hg_gs(p,msqc) 
      elseif ((case .eq. 'tt_bbh') .or .(case .eq. 'tt_bbl')) then
        write(6,*) 'No real correction to t_bbar yet'
        write(6,*) 'case=',case
        stop
      elseif ((case .eq. 'ttbdkl') .or .(case .eq. 'ttbdkh')) then
        write(6,*) 'No real correction to t_bbar yet'
        write(6,*) 'case=',case
        stop
      elseif ((case .eq. 'tt_tot') .or. (case .eq. 'cc_tot')
     .   .or. (case .eq. 'bb_tot')) then
       write(6,*) 'Real corrections not yet included!'
       stop
      elseif (case .eq. 't_bbar') then
        call qqb_tbb(p,msq)
      endif
      
      do j=0,ndmax
      xmsq(j)=0d0
      enddo
      
      do j=-nf,nf
      do k=-nf,nf
     
      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
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

      if (xmsq(0) .eq. 0d0) goto 999

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat

c--- if this dipole has no contribution, go to end of loop
        if (xmsq(nd) .eq. 0d0) goto 997

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
          jets=nqcdjets+1
        else
          do j=1,7
          do k=1,4
          q(j,k)=p(j,k)
          enddo
          enddo
          call genclust2(q,rcut,pjet,0)
          if (((jets .ne. nqcdjets) .and. (inclusive .eqv. .false.))
     .    .or.((jets .lt. nqcdjets) .and. (inclusive .eqv. .true.)))then
            njetzero=njetzero+1
            ncutzero=ncutzero-1
            failed=.true.
            goto 996
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
          jets=nqcdjets
        else
          do j=1,7
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo
c--- cluster partons (nqcdstart) to (nqcdstart+nqcdjets-1)
c--- if nqcdjets=0, no clustering is performed and pjet=p      
          call genclust2(q,rcut,pjet,1)
        endif
        endif

        call dotem(nvec,pjet,s)

c--- Check that there is the right number of jets
        failed=.false.
        if ((clustering .and. (jets .ne. nqcdjets)
     .         .and. (inclusive .eqv. .false.)) .or.
     .      (clustering .and. (jets .lt. nqcdjets)
     .         .and. (inclusive .eqv. .true.))) then
          failed=.true.
          goto 996
        endif

c--- Apply the (counter-)event cuts
        if (makecuts) then
          failed=gencuts(p,pjet,jets)
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

        val=xmsq(nd)*wgt
c--- update the maximum weight so far, if necessary
        if (dabs(val) .gt. wtmax) then
          wtmax=dabs(val)
        endif

c---if we're binning, add to histo too
        if (bin) then
          val=val/dfloat(itmx)
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
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      
c--- re-weight to improve m56 distribution
c      do j=1,7
c      do k=1,4
c      q(j,k)=p(j,k)
c      enddo
c      enddo
c      call genclust2(q,rcut,pjet,0)
c      realint=realint*reweight(pjet)
c--- end re-weight


      return

 999  realint=0d0
      ntotzero=ntotzero+1
 
      return
      end
















