      double precision function virtint(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'dprodx.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'PR_cs_new.f'
      include 'msq_cs.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'clustering.f'
      include 'efficiency.f'
      include 'bbproc.f'
      include 'lc.f'
      include 'facscale.f'
      include 'maxwt.f'
      double precision mqq(0:2,fn:nf,fn:nf)
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)
      double precision AP(-1:1,-1:1,3)

      integer ih1,ih2,j,k,cs,jets,nvec,is,ia,ib,ic
      double precision p(mxpart,4),pjet(mxpart,4),r(mxdim),W,sqrts,xmsq,
     . val,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf)
      double precision pswt,xjac,scalestart,pdfscale,
     . wgt,msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq
      double precision xx(2),z,x1onz,x2onz,flux,vol,taumin,omz,
     . BrnRat,xmsq_old
      double precision bbsqmin,bbsqmax,wsqmin,wsqmax,rcut     
      integer nqcdjets,nqcdstart,nshot,rvcolourchoice
      character*6 case
      character pdlabel*7,jetlabel(mxpart)*2
      common/nqcdjets/nqcdjets,nqcdstart
      common/parts/jets,jetlabel
      common/process/case
      logical bin,makecuts,gencuts,first
      common/limits/bbsqmin,bbsqmax,wsqmin,wsqmax      
      common/density/ih1,ih2
      common/energy/sqrts
      common/pdlabel/pdlabel
      common/bin/bin
      common/x1x2/xx
      common/taumin/taumin
      common/BrnRat/BrnRat
      common/makecuts/makecuts
      common/rcut/rcut
      common/rvcolourchoice/rvcolourchoice
      data p/48*0d0/
      data nshot/1/
      data first/.true./
      save first,scalestart
      if (first) then
         first=.false.
         scalestart=scale
      endif

      ntotshot=ntotshot+1
      virtint=0d0

      W=sqrts**2

      if (
     .        (case .eq. 'tt_bbl')
     .   .or. (case .eq. 'tt_bbh')
     .   .or. (case .eq. 'tautau')
     .   .or. (case .eq. 'vlchk6')) then
          write(6,*) 'Not implemented for virtual yet'
          write(6,*) 'case=',case
          stop
      elseif ((case .eq. 'ttbdkl')
     .   .or. (case .eq. 'ttbdkh')) then
          npart=6
          call gen6_rap(r,p,pswt,*999)
      elseif ((case .eq. 'W_only')
     .   .or. (case .eq. 'Z_only')
     .   .or. (case .eq. 'Hbbbar')) then
          npart=2
          call gen2(r,p,pswt,*999)
      elseif ((case .eq. 'tt_tot')
     .   .or. (case .eq. 'bb_tot')
     .   .or. (case .eq. 'cc_tot')) then
          npart=2
          call gen2m(r,p,pswt,*999)
      elseif  (case .eq. 'Wgamma') then
        npart=3
        call gen3b(r,p,pswt,*999)
      elseif ((case .eq. 'H_1jet')
     .   .or. (case .eq. 'Z_1jet')
     .   .or. (case .eq. 'W_1jet')) then
        npart=3
        call gen_njets(r,1,p,pswt,*999)
      elseif (
     .        (case .eq. 'W_2jet')
     .   .or. (case .eq. 'Z_2jet')
     .        ) then
        npart=4
        call gen_njets(r,2,p,pswt,*999)
c        call gen4(r,p,pswt,*999)
      elseif (
     .        (case .eq. 'tt_bbh')
     .   .or. (case .eq. 'tt_bbl')) then
          npart=6
          call gen6(r,p,pswt,*999)
      else
        npart=4
        call gen4(r,p,pswt,*999)      
      endif
      nvec=npart+2
      call dotem(nvec,p,s)
c---impose mass cuts on final state
      call masscuts(s,*999)
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)

      if (scalestart .lt. 0d0) call scaleset(scalestart,p)      

      do j=-nf,nf
      fx1z(j)=0d0
      fx2z(j)=0d0
      enddo
            
      if (case .eq. 'H_1jet') then
        pdfscale=facscale
      else
        pdfscale=scale
      endif   
            
      call fdist(pdlabel,ih1,xx(1),pdfscale,fx1)
      call fdist(pdlabel,ih2,xx(2),pdfscale,fx2)

      z=r(ndim)**2
      if (nshot .eq. 1) z=0.95d0
      xjac=two*sqrt(z)
      if (z .gt. xx(1)) then
         x1onz=xx(1)/z
         call fdist(pdlabel,ih1,x1onz,pdfscale,fx1z)
      endif
      if (z .gt. xx(2)) then
         x2onz=xx(2)/z
         call fdist(pdlabel,ih2,x2onz,pdfscale,fx2z)
      endif         
      
      omz=1d0-z

      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)

c--- to test poles, we need colourchoice=0, but save real value
      if (nshot .eq. 1) then
        rvcolourchoice=colourchoice
        colourchoice=0
      endif
      
   12 continue
c--- point to restart from when checking epsilon poles

      AP(q,q,1)=+ason2pi*Cf*1.5d0*epinv
      AP(q,q,2)=+ason2pi*Cf*(-1d0-z)*epinv
      AP(q,q,3)=+ason2pi*Cf*2d0/omz*epinv
      AP(a,a,1)=+ason2pi*Cf*1.5d0*epinv
      AP(a,a,2)=+ason2pi*Cf*(-1d0-z)*epinv
      AP(a,a,3)=+ason2pi*Cf*2d0/omz*epinv

      AP(q,g,1)=0d0
      AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epinv
      AP(q,g,3)=0d0
      AP(a,g,1)=0d0
      AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epinv
      AP(a,g,3)=0d0

      AP(g,q,1)=0d0
      AP(g,q,2)=ason2pi*Cf*(1d0+omz**2)/z*epinv
      AP(g,q,3)=0d0
      AP(g,a,1)=0d0
      AP(g,a,2)=ason2pi*Cf*(1d0+omz**2)/z*epinv
      AP(g,a,3)=0d0

      AP(g,g,1)=+ason2pi*b0*epinv
      AP(g,g,2)=+ason2pi*xn*2d0*(1d0/z+z*omz-2d0)*epinv
      AP(g,g,3)=+ason2pi*xn*2d0/omz*epinv

      do ia=-1,+1
      do ib=-1,+1
      do ic=-1,+1
      do is=1,3
        Q1(ia,ib,ic,is)=0d0
        Q2(ia,ib,ic,is)=0d0
      do cs=0,2
        R1(ia,ib,ic,cs,is)=0d0
        R2(ia,ib,ic,cs,is)=0d0
      enddo
      enddo
      enddo
      enddo
      enddo
     
      if     (case .eq. 'Wbbbar') then
         call qqb_Wbb(p,msq)
         call qqb_Wbb_v(p,msqv)
         call qqb_Wbb_z(p,z)
      elseif (case .eq. 'Zbbbar') then
         call qqb_Zbb(p,msq)
         call qqb_Zbb_v(p,msqv)
         call qqb_Zbb_z(p,z)
      elseif (case .eq. 't_bbar') then
         call qqb_tbb(p,msq)
      elseif (case .eq. 'WZbbar') then
         call qqb_wz(p,msq)
         call qqb_wz_v(p,msqv)
         call qqb_wz_z(p,z)
      elseif (case .eq. 'WmZbbr') then
         call qqb_wz(p,msq)
         call qqb_wz_v(p,msqv)
         call qqb_wz_z(p,z)
      elseif (case .eq. 'WWqqbr') then
         call qqb_WW(p,msq)
         call qqb_WW_v(p,msqv)
         call qqb_WW_z(p,z)
      elseif (case .eq. 'ZZlept') then
         call qqb_ZZ(p,msq)
         call qqb_ZZ_v(p,msqv)
         call qqb_ZZ_z(p,z)
      elseif (case .eq. 'WHbbar') then
         call qqb_WH(p,msq)
         call qqb_WH_v(p,msqv)
         call qqb_WH_z(p,z)
      elseif (case .eq. 'H_1jet') then
         call qqb_Hg(p,msq)
         call qqb_Hg_v(p,msqv)
         call qqb_Hg_z(p,z)
      elseif (case .eq. 'HWW_4l') then
         call qqb_HWW(p,msq)
         call qqb_HWW_v(p,msqv)
         call qqb_HWW_z(p,z)
      elseif (case .eq. 'Hbbbar') then
         call qqb_Hbbbar(p,msq)
         call qqb_Hbbbar_v(p,msqv)
         call qqb_Hbbbar_z(p,z)
      elseif (case .eq. 'HZZ_4l') then
         call qqb_HZZ(p,msq)
         call qqb_HZZ_v(p,msqv)
         call qqb_HZZ_z(p,z)
      elseif (case .eq. 'W_only') then
         call qqb_w(p,msq)
         call qqb_w_v(p,msqv)
         call qqb_w_z(p,z)
      elseif (case .eq. 'W_1jet') then
         call qqb_w_g(p,msq)
         call qqb_w1jet_v(p,msqv)
         call qqb_w1jet_z(p,z)
      elseif (case .eq. 'Wgamma') then
         call qqb_wgam(p,msq)
         call qqb_wgam_v(p,msqv)
         call qqb_wgam_z(p,z)
      elseif (case .eq. 'W_2jet') then         
         call qqb_w2jetx(p,msq,mqq,msqx,msqx_cs)
         call qqb_w2jet_v(p,msqv)
         call qqb_w2jet_z(p,z)
      elseif (case .eq. 'Z_only') then
         call qqb_z(p,msq)
         call qqb_z_v(p,msqv)
         call qqb_z_z(p,z)
      elseif (case .eq. 'Z_1jet') then
         call qqb_z1jet(p,msq)
         call qqb_z1jet_v(p,msqv)
         call qqb_z1jet_z(p,z)
      elseif (case .eq. 'Z_2jet') then
         call qqb_z2jetx(p,msq,mqq,msqx,msqx_cs)
         call qqb_z2jet_v(p,msqv)
         call qqb_z2jet_z(p,z)
      elseif (case .eq. 'ZHbbar') then
         call qqb_ZH(p,msq)
         call qqb_ZH_v(p,msqv)
         call qqb_ZH_z(p,z)
c      elseif ((case .eq. 'tt_tot')
c     .   .or. (case .eq. 'bb_tot')
c     .   .or. (case .eq. 'cc_tot')) then
c         call qqb_QQb(p,msq)
c         call qqb_QQb_v(p,msqv)
c      elseif ((case .eq. 'ttbdkl') .or. (case .eq. 'ttbdkh')) then
c         call qqb_QQbdk(p,msq)
c         write(6,*) 'msq(1,-1)',msq(1,-1)
c         write(6,*) 'msq(0,0)',msq(0,0)
c         write(6,*) 'msq(-1,1)',msq(-1,1)
c         call qqb_QQbdk_v(p,msqv)
c         write(6,*) 'msqv(1,-1)',msqv(1,-1)
c         write(6,*) 'msqv(0,0)',msqv(0,0)
c         write(6,*) 'msqv(-1,1)',msqv(-1,1)
c         pause
      elseif (case .eq. 'vlchk4') then
         taumin=0.0001d0
         bbsqmax=W
         bbsqmin=0d0
         call qqb_vol(p,msq)
         flux=9d0/vol(W,4)

         do j=-nf,nf
         fx1(j)=0d0
         fx2(j)=0d0
         enddo
         fx1(2)=1d0
         fx2(-1)=1d0

      endif

C---initialize to zero
      xmsq=0d0
C----------------------**************************************

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

c--- The variables R1 and R2 provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (R1(a,b,c,cs,is)) and leg 2 (R2(a,b,c,cs,is))
c--- In each case the parton labelling is using the normal QM notation of 
c--- putting everything backward
c---       emitted line after emission =    a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      if ((case .eq. 'W_2jet') .or. (case .eq. 'Z_2jet')) then
c--- SUM BY COLOUR STRUCTURES: W/Z + 2 jet only

      xmsq=xmsq+fx1(j)*fx2(k)*(
     . msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))

      if ((j .gt. 0) .and. (k.gt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(q,q,q,cs,1)-R1(q,q,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(q,q,q,cs,1)-R2(q,q,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R1(q,q,q,cs,2)+R1(q,q,q,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,q,2)+R1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R2(q,q,q,cs,2)+R2(q,q,q,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,q,2)+R2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      elseif ((j .lt. 0) .and. (k.lt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                 +R1(a,a,a,cs,1)-R1(a,a,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,a,cs,1)-R2(a,a,a,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,a,cs,2)+R1(a,a,a,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,a,2)+R1(g,a,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,a,cs,2)+R2(a,a,a,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,a,2)+R2(g,a,a,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      elseif     ((j .gt. 0) .and. (k.lt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(q,q,a,cs,1)-R1(q,q,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,q,cs,1)-R2(a,a,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R1(q,q,a,cs,2)+R1(q,q,a,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,q,2)+R1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,q,cs,2)+R2(a,a,q,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,a,2)+R2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      elseif ((j .lt. 0) .and. (k.gt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                 +R1(a,a,q,cs,1)-R1(a,a,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(q,q,a,cs,1)-R2(q,q,a,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,q,cs,2)+R1(a,a,q,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,a,2)+R1(g,a,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R2(q,q,a,cs,2)+R2(q,q,a,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,q,2)+R2(g,q,a,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      
      elseif ((j .eq. g) .and. (k .eq. g)) then
      do cs=0,2
      msq_qg=msq_cs(cs,+5,g)+msq_cs(cs,+4,g)+msq_cs(cs,+3,g)
     &      +msq_cs(cs,+2,g)+msq_cs(cs,+1,g)
     &      +msq_cs(cs,-5,g)+msq_cs(cs,-4,g)+msq_cs(cs,-3,g)
     &      +msq_cs(cs,-2,g)+msq_cs(cs,-1,g)
      msq_gq=msq_cs(cs,g,+5)+msq_cs(cs,g,+4)+msq_cs(cs,g,+3)
     &      +msq_cs(cs,g,+2)+msq_cs(cs,g,+1)
     &      +msq_cs(cs,g,-5)+msq_cs(cs,g,-4)+msq_cs(cs,g,-3)
     &      +msq_cs(cs,g,-2)+msq_cs(cs,g,-1)
      xmsq=xmsq
     & +msq_cs(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3)
     &                 +R2(g,g,g,cs,1)-R2(g,g,g,cs,3))*fx1(g)*fx2(g)
     & +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,g,cs,3)+R1(g,g,g,cs,2))
     & + msq_qg*(AP(q,g,2)+R1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     & +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R2(g,g,g,cs,3)+R2(g,g,g,cs,2))
     & + msq_gq*(AP(q,g,2)+R2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
      enddo
      
      elseif ((j .eq. g) .and. (k .gt. 0)) then
      do cs=0,2
      msq_aq=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      msq_qq=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      xmsq=xmsq     
     & +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,q,cs,1)-R1(g,g,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(q,q,g,cs,1)-R2(q,q,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,q,cs,2)+R1(g,g,q,cs,3))
     & + msq_aq*(AP(a,g,2)+R1(a,g,q,cs,2))
     & + msq_qq*(AP(q,g,2)+R1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(q,q,2)+AP(q,q,3)
     &                +R2(q,q,g,cs,2)+R2(q,q,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,q,2)+R2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
      enddo

      elseif ((j .eq. g) .and. (k .lt. 0)) then
      do cs=0,2
      msq_qa=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      msq_aa=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      xmsq=xmsq     
     & +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,a,cs,1)-R1(g,g,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,g,cs,1)-R2(a,a,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,a,cs,2)+R1(g,g,a,cs,3))
     & + msq_qa*(AP(q,g,2)+R1(q,g,a,cs,2))
     & + msq_aa*(AP(a,g,2)+R1(a,g,a,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,g,cs,2)+R2(a,a,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,a,2)+R2(g,a,g,cs,2)))*fx1(g)*fx2z(k)/z
      enddo

      elseif ((j .gt. 0) .and. (k .eq. g)) then
      do cs=0,2
      msq_qa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
      msq_qq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
       xmsq=xmsq     
     &+ msq_cs(cs,j,g)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(q,q,g,cs,1)-R1(q,q,g,cs,3)
     &                 +R2(g,g,q,cs,1)-R2(g,g,q,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &+(msq_cs(cs,j,g)*(AP(q,q,2)+AP(q,q,3)
     &                 +R1(q,q,g,cs,2)+R1(q,q,g,cs,3))
     &+ msq_cs(cs,g,g)*(AP(g,q,2)+R1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                 +R2(g,g,q,cs,2)+R2(g,g,q,cs,3))
     &+ msq_qa*(AP(a,g,2)+R2(a,g,q,cs,2))
     &+ msq_qq*(AP(q,g,2)+R2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z
      enddo
      
      elseif ((j .lt. 0) .and. (k .eq. g)) then
      do cs=0,2
      msq_aq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
      msq_aa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
       xmsq=xmsq     
     & + msq_cs(cs,j,g)*(AP(a,a,1)-AP(a,a,3)
     &                  +R1(a,a,g,cs,1)-R1(a,a,g,cs,3)
     &                  +AP(g,g,1)-AP(g,g,3)
     &                  +R2(g,g,a,cs,1)-R2(g,g,a,cs,3))*fx1(j)*fx2(g)
     & +(msq_cs(cs,j,g)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,g,cs,2)+R1(a,a,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,a,2)+R1(g,a,g,cs,2)))*fx1z(j)/z*fx2(g)
     & +(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R2(g,g,a,cs,2)+R2(g,g,a,cs,3))
     & + msq_aq*(AP(q,g,2)+R2(q,g,a,cs,2))
     & + msq_aa*(AP(a,g,2)+R2(a,g,a,cs,2)))*fx1(j)*fx2z(g)/z
      enddo
      endif

      else

c--- SUM BY TOTAL MATRIX ELEMENTS: everything else
C--QQ
      if     ((j .gt. 0) .and. (k.gt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,q,1)-Q1(q,q,q,3)
     &                +AP(q,q,1)-AP(q,q,3)+Q2(q,q,q,1)-Q2(q,q,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,q,2)+Q1(q,q,q,3))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,q,2)+Q2(q,q,q,3))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,q,2)))*fx1(j)*fx2z(k)/z
C--QbarQbar
      elseif ((j .lt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,a,1)-Q1(a,a,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,a,1)-Q2(a,a,a,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,a,2)+Q1(a,a,a,3))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,a,2)+Q2(a,a,a,3))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,a,2)))*fx1(j)*fx2z(k)/z

C--QQbar
      elseif ((j .gt. 0) .and. (k.lt.0)) then
      xmsq=xmsq+(msqv(j,k)
     & + msq(j,k)*(one+AP(q,q,1)-AP(q,q,3)+Q1(q,q,a,1)-Q1(q,q,a,3)
     &                +AP(a,a,1)-AP(a,a,3)+Q2(a,a,q,1)-Q2(a,a,q,3)))
     &                *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,a,3)+Q1(q,q,a,2))
     & + msq(g,k)*(AP(g,q,2)+Q1(g,q,a,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,q,3)+Q2(a,a,q,2))
     & + msq(j,g)*(AP(g,a,2)+Q2(g,a,q,2)))*fx1(j)*fx2z(k)/z

      elseif ((j .lt. 0) .and. (k.gt.0)) then
C--QbarQ
      xmsq=xmsq+(msqv(j,k)
     & +msq(j,k)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,q,1)-Q1(a,a,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,a,1)-Q2(q,q,a,3)))
     &               *fx1(j)*fx2(k)
     & +(msq(j,k)*(AP(a,a,3)+AP(a,a,2)+Q1(a,a,q,3)+Q1(a,a,q,2))
     & + msq(g,k)*(AP(g,a,2)+Q1(g,a,q,2)))*fx1z(j)/z*fx2(k)
     & +(msq(j,k)*(AP(q,q,3)+AP(q,q,2)+Q2(q,q,a,3)+Q2(q,q,a,2))
     & + msq(j,g)*(AP(g,q,2)+Q2(g,q,a,2)))*fx1(j)*fx2z(k)/z

      elseif ((j .eq. g) .and. (k.eq.g)) then
C--gg
      if (case .eq. 'Zbbbar') then
         xmsq=xmsq+msqv(g,g)*fx1(g)*fx2(g)
         do cs=0,2
      xmsq=xmsq+(
     &+msq_cs(cs,g,g)*(one+two*(AP(g,g,1)-AP(g,g,3))
     &                    +R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &                    +R2(g,g,g,cs,1)-R2(g,g,g,cs,3)))*fx1(g)*fx2(g)
     &+msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                +R1(g,g,g,cs,2)+R1(g,g,g,cs,3))*fx1z(g)/z*fx2(g)
     &+msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                +R2(g,g,g,cs,2)+R2(g,g,g,cs,3))*fx1(g)*fx2z(g)/z
        enddo
      else
      msq_qg=msq(+5,g)+msq(+4,g)+msq(+3,g)+msq(+2,g)+msq(+1,g)
     &      +msq(-5,g)+msq(-4,g)+msq(-3,g)+msq(-2,g)+msq(-1,g)
      msq_gq=msq(g,+5)+msq(g,+4)+msq(g,+3)+msq(g,+2)+msq(g,+1)
     &      +msq(g,-5)+msq(g,-4)+msq(g,-3)+msq(g,-2)+msq(g,-1)

      xmsq=xmsq+(msqv(g,g)
     & +msq(g,g)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,g,1)-Q1(g,g,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,g,1)-Q2(g,g,g,3)))
     &               *fx1(g)*fx2(g)
     & +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,g,2)+Q1(g,g,g,3))
     & + msq_qg*(AP(q,g,2)+Q1(q,g,g,2)))*fx1z(g)/z*fx2(g)
     & +(msq(g,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,g,2)+Q2(g,g,g,3))
     & + msq_gq*(AP(q,g,2)+Q2(q,g,g,2)))*fx1(g)*fx2z(g)/z
      endif

      elseif (j .eq. g) then
C--gQ
       if    (k .gt. 0) then
       msq_aq=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,q,1)-Q1(g,g,q,3)
     &               +AP(q,q,1)-AP(q,q,3)+Q2(q,q,g,1)-Q2(q,q,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,q,2)+Q1(g,g,q,3))
     & +   msq_aq*(AP(a,g,2)+Q1(a,g,q,2)))*fx1z(g)/z*fx2(k)
     & +(msq(g,k)*(AP(q,q,2)+AP(q,q,3)+Q2(q,q,g,2)+Q2(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q2(g,q,g,2)))*fx1(g)*fx2z(k)/z
C--gQbar

       elseif (k.lt.0) then
       msq_qa=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       xmsq=xmsq+(msqv(g,k)
     & +msq(g,k)*(one+AP(g,g,1)-AP(g,g,3)+Q1(g,g,a,1)-Q1(g,g,a,3)
     &               +AP(a,a,1)-AP(a,a,3)+Q2(a,a,g,1)-Q2(a,a,g,3)))
     &               *fx1(g)*fx2(k)
     & +(msq(g,k)*(AP(g,g,2)+AP(g,g,3)+Q1(g,g,a,2)+Q1(g,g,a,3))
     & +   msq_qa*(AP(q,g,2)+Q1(q,g,a,2)))*fx1z(g)/z*fx2(k)
     & +(msq(g,k)*(AP(a,a,2)+AP(a,a,3)+Q2(a,a,g,2)+Q2(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q2(g,a,g,2)))*fx1(g)*fx2z(k)/z
       endif
C--Qg
      elseif (k .eq. g) then
       if     (j.gt.0) then
       msq_qa=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(one
     &               +AP(q,q,1)-AP(q,q,3)+Q1(q,q,g,1)-Q1(q,q,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,q,1)-Q2(g,g,q,3)))
     &               *fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(q,q,2)+AP(q,q,3)+Q1(q,q,g,2)+Q1(q,q,g,3))
     & + msq(g,g)*(AP(g,q,2)+Q1(g,q,g,2)))*fx1z(j)/z*fx2(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,q,2)+Q2(g,g,q,3))
     & +   msq_qa*(AP(a,g,2)+Q2(a,g,q,2)))*fx1(j)*fx2z(g)/z
C--Qbarg
       elseif (j.lt.0) then
       msq_aq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       xmsq=xmsq+(msqv(j,g)
     & +msq(j,g)*(one+AP(a,a,1)-AP(a,a,3)+Q1(a,a,g,1)-Q1(a,a,g,3)
     &               +AP(g,g,1)-AP(g,g,3)+Q2(g,g,a,1)-Q2(g,g,a,3)))
     &                *fx1(j)*fx2(g)
     & +(msq(j,g)*(AP(a,a,2)+AP(a,a,3)+Q1(a,a,g,2)+Q1(a,a,g,3))
     & + msq(g,g)*(AP(g,a,2)+Q1(g,a,g,2)))*fx1z(j)/z*fx2(g)
     & +(msq(j,g)*(AP(g,g,2)+AP(g,g,3)+Q2(g,g,a,3)+Q2(g,g,a,2))
     & + msq_aq*(AP(q,g,2)+Q2(q,g,a,2)))*fx1(j)*fx2z(g)/z
       endif
      endif
      
      endif
      
 20   continue

      enddo
      enddo

c--- code to check that epsilon poles cancel      
      if (nshot .eq. 1) then
        if (xmsq .eq. 0d0) goto 999
        xmsq_old=xmsq
        nshot=nshot+1
        epinv=0d0
        epinv2=0d0
        goto 12
      elseif (nshot .eq. 2) then
        nshot=nshot+1
        if (abs(xmsq_old/xmsq-1d0) .gt. 1d-7) then
          write(6,*) 'epsilon fails to cancel'
          write(6,*) 'xmsq (epinv=large) = ',xmsq_old
          write(6,*) 'xmsq (epinv=zero ) = ',xmsq
          stop
        else
          write(6,*) 'Poles cancelled!'
          colourchoice=rvcolourchoice
        endif
      endif
      
c--- cluster partons 5 to (4+nqcdjets)
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
        if ((jets .ne. nqcdjets) .and. (nqcdjets .gt. 0)) then
          njetzero=njetzero+1
          goto 999
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
         
      virtint=flux*xjac*pswt*xmsq/BrnRat

      val=virtint*wgt 
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


