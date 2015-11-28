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
      integer ih1,ih2,j,k,nd,nmax,nmin,jets,nproc,nvec     
      double precision vector(mxdim),W,val,xint
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,wbb,bit,sum,pt,ayrap,msqs(-nf:nf,-nf:nf)
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4),ptjet,bclustmass,rcut
      character*6 case
      integer nqcdjets,nqcdstart
      character pdlabel*7,jetlabel(mxpart)*2
      common/nqcdjets/nqcdjets,nqcdstart
      common/parts/jets,jetlabel
      common/process/case
      common/xreal/xreal,xreal2
      logical bin,makecuts,first,madejetcuts,failed,cuts,bbproc
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
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                     
     . ,SUUB_VEVEBBBG,SDDB_VEVEBBBG,SUUB_EMEPBBBG,SDDB_EMEPBBBG
     . ,SUG_VEVEBBBU,SUBG_VEVEBBBUB,SGU_VEVEBBBU,SGUB_VEVEBBBUB
     . ,SUUB_VEVEG,SDDB_VEVEG,SUUB_EMEPG
     . ,SUG_VEVEU,SUBG_VEVEUB,SGU_VEVEU,SGUB_VEVEUB
      data first/.true./
      
      pswt=0d0
      realint=0d0
      

      W=sqrts**2
      
      if (first) then
      write(6,*)
      write(6,*) 'nmin=',nmin,',nmax=',nmax
      write(6,*)
      first=.false.
      endif
 30   continue
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
      call dotem(nvec,p,s)

c      write(6,*) 's(1,5),s(1,6),s(2,5),s(2,6),s(5,6),s(3,4)',
c     . s(1,5),s(1,6),s(2,5),s(2,6),s(5,6),s(3,4)
c---impose cuts on final state
      call masscuts(s,*999)
c----reject event if any  s(i,j)  is too small
      call smalls(s,npart,*999)
      
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
c      call qqb_Zbb_g(p,msq)
c      write(*,*) 'OLD u ',msq(2,-2)
c      write(*,*) 'OLD ub',msq(-2,2)
c      write(*,*) 'OLD d ',msq(1,-1)
c      write(*,*) 'OLD db',msq(-1,1)
c      write(*,*) 'OLD ug ',msq(2,0)
c      write(*,*) 'OLD ubg',msq(-2,0)
c      write(*,*) 'OLD gu',msq(0,2)
c      write(*,*) 'OLD gub',msq(0,-2)
      call qqb_Zbb_g(p,msq)
c      write(*,*) 'NEW u ',msq(2,-2)
c      write(*,*) 'NEW ub',msq(-2,2)
c      write(*,*) 'NEW d ',msq(1,-1)
c      write(*,*) 'NEW db',msq(-1,1)      
c      write(*,*) 'NEW ug ',msq(2,0)
c      write(*,*) 'NEW ubg',msq(-2,0)
c      write(*,*) 'NEW gu',msq(0,2)
c      write(*,*) 'NEW gub',msq(0,-2)
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
c        p7(j)=p(7,k)
c      enddo            
c      call initialize
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      enddo
c      enddo
c      do j=-nf,nf
c      do k=-nf,nf
c      if     ((j.gt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             suub_vevebbbg(p1,p2,p3,p4,p5,p6,p7)
c      elseif ((j.lt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             suub_vevebbbg(p2,p1,p3,p4,p5,p6,p7)
c      elseif ((j.gt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             sddb_vevebbbg(p1,p2,p3,p4,p5,p6,p7)
c      elseif ((j.lt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             sddb_vevebbbg(p2,p1,p3,p4,p5,p6,p7)
c      elseif ((j.gt.0) .and. (k.eq.0) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             sug_vevebbbu(p1,p2,p3,p4,p5,p6,p7)
c      elseif ((j.lt.0) .and. (k.eq.0) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             subg_vevebbbub(p1,p2,p3,p4,p5,p6,p7)
c      elseif ((k.gt.0) .and. (j.eq.0) .and. (k/2 .eq. (k+1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             sgu_vevebbbu(p1,p2,p3,p4,p5,p6,p7)
c      elseif ((k.lt.0) .and. (j.eq.0) .and. (k/2 .eq. (k-1)/2)) then
c        msq(j,k)=3d0*gsq**3*esq**2/(fourpi/128d0)**2*
c     .             sgub_vevebbbub(p1,p2,p3,p4,p5,p6,p7)
c      endif
c      enddo
c      enddo
c      write(*,*) 'MAD u ',msq(2,-2)
c      write(*,*) 'MAD ub',msq(-2,2)
c      write(*,*) 'MAD d ',msq(1,-1)
c      write(*,*) 'MAD db',msq(-1,1)
c      write(*,*) 'MAD ug ',msq(2,0)
c      write(*,*) 'MAD ubg',msq(-2,0)
c      write(*,*) 'MAD gu',msq(0,2)
c      write(*,*) 'MAD gub',msq(0,-2)
c      pause
        call qqb_zbb_gs(p,msqc) 
      elseif (case .eq. 'WHbbar') then
        call qqb_wh_gs(p,msqc)     
        call qqb_wh_g(p,msq)      
      elseif (case .eq. 'ZHbbar') then
        call qqb_zh_gs(p,msqc)     
        call qqb_zh_g(p,msq)      
      elseif (case .eq. 'W_2jet') then
        call dotem(7,p,s)
        call qqb_w2jet_g(p,msq)      
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
c      write(*,*) 'NEW u ',msq(2,-2)
c      write(*,*) 'NEW ub',msq(-2,2)
c      write(*,*) 'NEW d ',msq(1,-1)
c      write(*,*) 'NEW db',msq(-1,1)      
c      write(*,*) 'NEW ug ',msq(2,0)
c      write(*,*) 'NEW ubg',msq(-2,0)
c      write(*,*) 'NEW gu',msq(0,2)
c      write(*,*) 'NEW gub',msq(0,-2)
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
c      enddo            
c      call initialize
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      enddo
c      enddo
c      do j=-nf,nf
c      do k=-nf,nf
c      if     ((j.gt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             suub_veveg(p1,p2,p3,p4,p5)
c        msq(j,k)=gsq*esq**2/(fourpi/128d0)**2*
c     .             suub_emepg(p1,p2,p3,p4,p5)
c      elseif ((j.lt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             suub_veveg(p2,p1,p3,p4,p5)
c      elseif ((j.gt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             sddb_veveg(p1,p2,p3,p4,p5)
c      elseif ((j.lt.0) .and. (k.eq.-j) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             sddb_veveg(p2,p1,p3,p4,p5)
c      elseif ((j.gt.0) .and. (k.eq.0) .and. (j/2 .eq. (j+1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             sug_veveu(p1,p2,p3,p4,p5)
c      elseif ((j.lt.0) .and. (k.eq.0) .and. (j/2 .eq. (j-1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             subg_veveub(p1,p2,p3,p4,p5)
c      elseif ((k.gt.0) .and. (j.eq.0) .and. (k/2 .eq. (k+1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             sgu_veveu(p1,p2,p3,p4,p5)
c      elseif ((k.lt.0) .and. (j.eq.0) .and. (k/2 .eq. (k-1)/2)) then
c        msq(j,k)=3d0*gsq*esq**2/(fourpi/128d0)**2*
c     .             sgub_veveub(p1,p2,p3,p4,p5)
c      endif
c      enddo
c      enddo
c      write(*,*) 'MAD u ',msq(2,-2)
c      write(*,*) 'MAD ub',msq(-2,2)
c      write(*,*) 'MAD d ',msq(1,-1)
c      write(*,*) 'MAD db',msq(-1,1)
c      write(*,*) 'MAD ug ',msq(2,0)
c      write(*,*) 'MAD ubg',msq(-2,0)
c      write(*,*) 'MAD gu',msq(0,2)
c      write(*,*) 'MAD gub',msq(0,-2)
c      pause
        call qqb_z_gs(p,msqc)     
      elseif (case .eq. 'W_1jet') then
        call qqb_w2jet(p,msq)      
        call qqb_w1jet_gs(p,msqc)  
c        call qqb_w1jet_soft(p,msqs)  
      elseif (case .eq. 'Z_1jet') then
        call qqb_z1jet_gs(p,msqc)  
        call qqb_z2jet(p,msq)      
c        call qqb_z1jet_soft(p,msqs)  
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
      
c      if (p(7,4) .lt. 2d0) then
c      if (-s(1,7) .lt. 20d0) then
c      write(*,*) 's17',s(1,7),'   s27',s(2,7)
c      write(*,*) 's57',s(5,7),'   s67',s(6,7)
c      write(*,*) 's56',s(5,6)
c      write(*,*) 'p(7,4)',p(7,4)
c      realint=0d0
c      do nd=1,ndmax
c      if (msqc(nd,0,0).ne.0d0) write(*,*) 'nd',nd,msqc(nd,0,0)
c      realint=realint+msqc(nd,0,0)
c      enddo
c      write(*,*) 'c1+2',msqc(1,0,0)+msqc(2,0,0)
c      write(*,*) 'c3+4',msqc(3,0,0)+msqc(4,0,0)
c      write(*,*) 'c5+8',msqc(5,0,0)+msqc(8,0,0)
c      write(*,*) 'c6+7',msqc(6,0,0)+msqc(7,0,0)
c      write(*,*) 'c5-8',msqc(5,0,0)+msqc(6,0,0)
c     .                 +msqc(7,0,0)+msqc(8,0,0)
c      write(*,*) 'msqc',realint
c      write(*,*) 'soft',msqs(0,0)
c      write(*,*) 'msq ',msq(0,0)
c      write(*,*) 'DIFF',1d0-msq(0,0)/realint
c      pause
c      endif
      
c--------debug
       if (debug) then
       if ((case .eq. 'W_1jet') .or. (case .eq. 'HWW_4l')) then
        sum=+msqc(1,2,-1)+msqc(2,2,-1)+msqc(3,2,-1) 
     .      +msqc(4,2,-1)+msqc(5,2,-1)+msqc(6,2,-1) 
        write(6,*) 'msq(2,-1),sum,ratio',msq(2,-1),sum,msq(2,-1)/sum 
        write(6,*) 

        sum=+msqc(1,-1,2)+msqc(2,-1,2)+msqc(3,-1,2) 
     .      +msqc(4,-1,2)+msqc(5,-1,2)+msqc(6,-1,2) 
        write(6,*) 'msq(-1,2),sum,ratio',msq(-1,2),sum,msq(-1,2)/sum 
        write(6,*) 

        elseif (case .eq. 'Z_1jet') then

        sum=+msqc(1,1,-1)+msqc(2,1,-1)+msqc(3,1,-1) 
     .      +msqc(4,1,-1)+msqc(5,1,-1)+msqc(6,1,-1) 
        write(6,*) 'msq(1,-1),sum,ratio',msq(1,-1),sum,msq(1,-1)/sum 
        write(6,*) 

        sum=+msqc(1,-1,1)+msqc(2,-1,1)+msqc(3,-1,1) 
     .      +msqc(4,-1,1)+msqc(5,-1,1)+msqc(6,-1,1) 
        write(6,*) 'msq(-1,1),sum,ratio',msq(-1,1),sum,msq(-1,1)/sum 
        write(6,*) 

        endif

        sum=+msqc(1,-1,0)+msqc(2,-1,0)+msqc(3,-1,0) 
     .      +msqc(4,-1,0)+msqc(5,-1,0)+msqc(6,-1,0) 
        write(6,*) 'msq(-1,0),sum,ratio',msq(-1,0),sum,msq(-1,0)/sum 
        write(6,*) 

        sum=+msqc(1,0,-1)+msqc(2,0,-1)+msqc(3,0,-1) 
     .      +msqc(4,0,-1)+msqc(5,0,-1)+msqc(6,0,-1) 
        write(6,*) 'msq(0,-1),sum,ratio',msq(0,-1),sum,msq(0,-1)/sum 
        write(6,*) 

        sum=+msqc(1,2,0)+msqc(2,2,0)+msqc(3,2,0) 
     .      +msqc(4,2,0)+msqc(5,2,0)+msqc(6,2,0) 
        write(6,*) 'msq(2,0),sum,ratio',msq(2,0),sum,msq(2,0)/sum 
        write(6,*) 

        sum=+msqc(1,0,2)+msqc(2,0,2)+msqc(3,0,2) 
     .      +msqc(4,0,2)+msqc(5,0,2)+msqc(6,0,2) 
        write(6,*) 'msq(0,2),sum,ratio',msq(0,2),sum,msq(0,2)/sum 
        write(6,*) 

        sum=+msqc(1,0,0)+msqc(2,0,0)+msqc(3,0,0) 
     .      +msqc(4,0,0)+msqc(5,0,0)+msqc(6,0,0) 
        write(6,*) 'msq(0,0),sum,ratio',msq(0,0),sum,msq(0,0)/sum 
        write(6,*) 

        pause
        goto 30
        endif 
c--------debug

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
          call genclust2(p,rcut,jets,pjet,jetlabel)
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
            do j=0,ndmax
            xmsq(j)=0d0
            enddo       
            goto 998
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
          call genclust(q,nqcdstart,nqcdstart+nqcdjets-1,
     .                   rcut,jets,pjet,jetlabel)
        endif
        endif

        call dotem(nvec,pjet,s)

c---check to see if (counter-)event passes cuts
        failed=.false.
        if ((clustering) .and. (jets .ne. nqcdjets)) then
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
c          call dotem(nvec,p,s)
c          if (( s(6,7) .lt. 50d0).and.(xmsq(0).ne.0d0)) then
c            write(*,*) 'dropping ',nd,xmsq(nd),' with ',jets,' jets'
c          endif
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

c      if (abs(xmsq(0)). gt. 2000d0) then
c      if ((xmsq(0) .ne. 0d0) .and.
c     .   ((-s(1,7) .lt. 1d0) .or. (-s(2,7) .lt. 1d0)
c     ..or.( s(5,7) .lt. 1d0) .or. ( s(6,7) .lt. 1d0))) then
c      if (p(7,4) .lt. 1d0) then

c      if (( -s(1,7) .lt. 100d0).and.( -s(2,7) .lt. 100d0)
c     . .and. (xmsq(0).ne.0d0)) then
c      write(*,*) 's17',s(1,7),'   s27',s(2,7)
c      write(*,*) 's57',s(5,7),'   s67',s(6,7)
c      write(*,*) 's56',s(5,6)
c      do j=1,7
c        write(*,*) 'part ',j,p(j,4),p(j,1),p(j,2),p(j,3)
c      enddo
c      write(*,*)
c      call genclust(p,nqcdstart,nqcdstart+nqcdjets,
c     .               rcut,jets,pjet,jetlabel)
c      do j=1,6
c        write(*,*) ' jet ',j,pjet(j,4),pjet(j,1),pjet(j,2),pjet(j,3)
c      enddo
c      write(*,*)
c      realint=0d0
c      do nd=1,ndmax
c      write(*,*) 'nd',nd,xmsq(nd)/flux/pswt*BrnRat
c      realint=realint+xmsq(nd)
c      enddo
c      write(*,*) 'msqc',realint/flux/pswt*BrnRat
c      write(*,*) 'msq ',xmsq(0)/flux/pswt*BrnRat
c      write(*,*) 'DIFF',1d0+xmsq(0)/realint
c      pause
c      endif

c      if ((-s(1,7) .lt. 50d0) .or. (-s(2,7) .lt. 50d0)) then
c      call compare(p)
c      if (xmsq(0) .eq. 0d0) then
c      write(*,*) s(5,6),s(1,5)*s(2,5)/s(1,2),s(1,6)*s(2,6)/s(1,2)
c      write(*,*) s(1,5),s(2,5)
c      write(*,*) 's17',s(1,7),'   s27',s(2,7)
c      if ((s(5,7) .lt. 50d0) .or. (s(6,7) .lt. 50d0)) then
c      write(*,*) 'msq ',xmsq(0)
c      realint=0d0
c      do nd=1,ndmax
c      write(*,*) 'nd',nd,xmsq(nd)
c      realint=realint+xmsq(nd)
c      enddo
c      write(*,*) 'msqc',realint 
c      write(*,*) 'DIFF',1d0+xmsq(0)/realint
c      if (xmsq(0) .eq. 0d0) then
c      if ((s(5,7) .lt. 50d0) .or. (s(6,7) .lt. 50d0)) then
c      write(*,*) 'nd     s56         pt5       pt6      xmsq'
c      write(*,*) '0',s(5,6),s(1,5)*s(2,5)/s(1,2),s(1,6)*s(2,6)/s(1,2),
c     . xmsq(0)
c      do nd=1,ndmax
c        do j=1,7
c        do k=1,4
c        q(j,k)=ptilde(nd,j,k)
c        enddo
c        enddo
c      call dotem(nvec,q,s)
c      write(*,*) nd,s(5,6),s(1,5)*s(2,5)/s(1,2),s(1,6)*s(2,6)/s(1,2),
c     . xmsq(nd)
c      enddo
c      write(*,*) 'total ',realint,'    with diff ',1d0+xmsq(0)/realint
c      do j=1,6
c      write(*,*) j,p(j,4),p(j,1),p(j,2),p(j,3)
c      enddo
c      do nd=1,2
c      do j=1,6
c      write(*,*) j,ptilde(nd,j,4),ptilde(nd,j,1),
c     .ptilde(nd,j,2),ptilde(nd,j,3)
c      enddo
c      enddo
c      call dotem(p,s)
c      write(*,*) 's12',s(1,2)
c      write(*,*) 's34',s(3,4)
c      write(*,*) 's15',s(1,5)
c      write(*,*) 's16',s(1,6)
c      write(*,*) 's25',s(2,5)
c      write(*,*) 's26',s(2,6)
c      write(*,*) 's56',s(5,6)
c      pause
c      endif            

 998  continue

      if (realwt) then
      realint=xmsq(0)
      else
      realint=xint
      endif
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+xint**2*wgt/dfloat(itmx)

      if (debug) write(6,*) 'flux',flux
      if (debug) write(6,*) 'pswt',pswt
      if (debug) write(6,*) 'xmsq(0)',xmsq(0)
      if (debug) write(6,*) 'xmsq(1)',xmsq(1)
      if (debug) write(6,*) 'xmsq(2)',xmsq(2)
      if (debug) write(6,*) 'xmsq(3)',xmsq(3)
      if (debug) write(6,*) 'xmsq(4)',xmsq(4)
      if (debug) write(6,*) 'xmsq(5)',xmsq(5)
      if (debug) write(6,*) 'xmsq(6)',xmsq(6)
      if (debug) pause

      return

 999  realint=0d0
      return
      end
















