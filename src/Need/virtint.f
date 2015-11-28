      double precision function virtint(r,wgt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'dprodx.f'
      include 'npart.f'
      include 'PR.f'
      include 'PR_cs.f'
      include 'msq_cs.f'
      include 'scale.f'
      include 'clustering.f'
      include 'flags.f'
      include 'efficiency.f'
      integer ih1,ih2,j,k,cs,jets,nproc,nvec
      double precision p(mxpart,4),pjet(mxpart,4),r(mxdim),W,sqrts,xmsq,
     . val,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf)
      double precision pswt,wbb,zz,xjac,
     . wgt,msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),msq_qq,msq_qg,msq_gq
      double precision xx(2),z,x1onz,x2onz,flux,vol,taumin,
     . xlog,BrnRat
      double precision bbsqmin,bbsqmax,wsqmin,wsqmax,ptjet,rcut     
      integer nqcdjets,nqcdstart
      character*6 case
      character pdlabel*7,jetlabel(mxpart)*2
      common/nqcdjets/nqcdjets,nqcdstart
      common/parts/jets,jetlabel
      common/process/case
      logical bin,makecuts,madejetcuts,cuts,bbproc
      common/limits/bbsqmin,bbsqmax,wsqmin,wsqmax      
      common/density/ih1,ih2
      common/energy/sqrts
      common/pdlabel/pdlabel
      common/bin/bin
      common/x1x2/xx
      common/taumin/taumin
      common/BrnRat/BrnRat
      common/xlog/xlog
      common/makecuts/makecuts
      common/nproc/nproc
      common/rcut/rcut
      
      data Rgg_q,Rgq_q,Rg_qq,Rg_qg,Rq_gg,Rq_gq,Rqq_g,Rqg_g,
     & Rqq_qb,Rq_qbqb,Rqb_qq,Rqbqb_q,Rgq_g,Rg_gq,Rgg_g,Rg_gg,
     & Pgg_q,Pgq_q,Pg_qq,Pg_qg,Pq_gg,Pq_gq,Pqq_g,Pqg_g,
     & Pqq_qb,Pq_qbqb,Pqb_qq,Pqbqb_q,Pgq_g,Pg_gq,Pgg_g,Pg_gg/32*0d0/

      data Rgg_g_cs,Rg_gg_cs,Pgg_g_cs,Pg_gg_cs,
     .     Rqq_qb_cs,Rq_qbqb_cs,Pqq_qb_cs,Pq_qbqb_cs,
     .     Rqbqb_q_cs,Rqb_qq_cs,Pqbqb_q_cs,Pqb_qq_cs,
     .     Rqq_g_cs,Rq_gg_cs,Pqq_g_cs,Pq_gg_cs,
     .     Rqbqb_g_cs,Rqb_gg_cs,Pqbqb_g_cs,Pqb_gg_cs,
     .     Rgg_q_cs,Rg_qq_cs,Pgg_q_cs,Pg_qq_cs,
     .     Rgg_qb_cs,Rg_qbqb_cs,Pgg_qb_cs,Pg_qbqb_cs,
     .     Rgq_g_cs,Rg_gq_cs,Rgqb_g_cs,Rg_gqb_cs,
     .     Pgq_g_cs,Pg_gq_cs,Pgqb_g_cs,Pg_gqb_cs,
     .     Rgq_qb_cs,Rq_gqb_cs,Rgqb_q_cs,Rqb_gq_cs,
     .     Pgq_qb_cs,Pq_gqb_cs,Pgqb_q_cs,Pqb_gq_cs,
     .     Pg_qg_cs,Pqg_g_cs,Rg_qg_cs,Rqg_g_cs,
     .     Rqq_q_cs,Rq_qq_cs,Rqbqb_qb_cs,Rqb_qbqb_cs,
     .     Pqq_q_cs,Pq_qq_cs,Pqbqb_qb_cs,Pqb_qbqb_cs/168*0d0/


      do cs=0,2
        Rgg_g_cs(cs)=0d0
        Rg_gg_cs(cs)=0d0
        Pgg_g_cs(cs)=0d0
        Pg_gg_cs(cs)=0d0
      enddo
      
      ntotshot=ntotshot+1
      virtint=0d0

      W=sqrts**2

      if (
     .      (case .eq. 'W_only')
     . .or. (case .eq. 'Z_only')
     . .or. (case .eq. 'Hbbbar')
     . ) then
        call gen2(r,p,pswt,*999)
        npart=2
      elseif ((case .eq. 'Z_1jet') .or. (case .eq. 'W_1jet')) then
        npart=3
        call gen3(r,p,pswt,*999)
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

c---impose experimental cuts on final state

      xlog=log(s(1,2)/scale**2)

      do j=-nf,nf
      fx1z(j)=0d0
      fx2z(j)=0d0
      enddo
            
      call fdist(pdlabel,ih1,xx(1),scale,fx1)
      call fdist(pdlabel,ih2,xx(2),scale,fx2)

      z=r(ndim)**2
      xjac=two*r(ndim)
      if (z .gt. xx(1)) then
         x1onz=xx(1)/z
         call fdist(pdlabel,ih1,x1onz,scale,fx1z)
      endif
      if (z .gt. xx(2)) then
         x2onz=xx(2)/z
         call fdist(pdlabel,ih2,x2onz,scale,fx2z)
      endif         

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
      elseif (case .eq. 'W_2jet') then         
         call qqb_w2jet(p,msq)
         call qqb_w2jet_v(p,msqv)
         call qqb_w2jet_z(p,z)
      elseif (case .eq. 'Z_only') then
         call qqb_z(p,msq)
         call qqb_z_v(p,msqv)
         call qqb_z_z(p,z)
      elseif (case .eq. 'Z_1jet') then
         call qqb_z_g(p,msq)
         call qqb_z1jet_v(p,msqv)
      elseif (case .eq. 'ZHbbar') then
         call qqb_ZH(p,msq)
         call qqb_ZH_v(p,msqv)
         call qqb_ZH_z(p,z)

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

      if ((case .eq. 'W_2jet') .or. (case .eq. 'Z_2jet')) then
c--- SUM BY COLOUR STRUCTURES: W/Z + 2 jet only
      xmsq=xmsq+fx1(j)*fx2(k)*(msq(j,k)+msqv(j,k))

      if     ((j .gt. 0) .and. (k.lt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +fx1(j)*fx2(k)*msq_cs(cs,j,k)*(-Pqq_qb_cs(cs)-Pq_qbqb_cs(cs))
     & +fx1z(j)*fx2(k)*msq_cs(cs,j,k)*(Pqq_qb_cs(cs)+Rqq_qb_cs(cs))/z
     & +fx1(j)*fx2z(k)*msq_cs(cs,j,k)*(Pq_qbqb_cs(cs)+Rq_qbqb_cs(cs))/z
      enddo
      
      elseif ((j .lt. 0) .and. (k.gt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +fx1(j)*fx2(k)*msq_cs(cs,j,k)*(-Pqbqb_q_cs(cs)-Pqb_qq_cs(cs))
     & +fx1z(j)*fx2(k)*msq_cs(cs,j,k)*(Pqbqb_q_cs(cs)+Rqbqb_q_cs(cs))/z
     & +fx1(j)*fx2z(k)*msq_cs(cs,j,k)*(Pqb_qq_cs(cs)+Rqb_qq_cs(cs))/z
      enddo
      
      elseif ((j .gt. 0) .and. (k.gt.0)) then
      do cs=0,2
      xmsq=xmsq
     & +fx1(j)*fx2(k)*msq_cs(cs,j,k)*(-Pqq_q_cs(cs)-Pq_qq_cs(cs))
     & +fx1z(j)*fx2(k)*msq_cs(cs,j,k)*(Pqq_q_cs(cs)+Rqq_q_cs(cs))/z
     & +fx1(j)*fx2z(k)*msq_cs(cs,j,k)*(Pq_qq_cs(cs)+Rq_qq_cs(cs))/z
      enddo
      
      elseif ((j .lt. 0) .and. (k.lt.0)) then
      do cs=0,2
      xmsq=xmsq
     &+fx1(j)*fx2(k)*msq_cs(cs,j,k)*(-Pqbqb_qb_cs(cs)-Pqb_qbqb_cs(cs))
     &+fx1z(j)*fx2(k)*msq_cs(cs,j,k)*(Pqbqb_qb_cs(cs)+Rqbqb_qb_cs(cs))/z
     &+fx1(j)*fx2z(k)*msq_cs(cs,j,k)*(Pqb_qbqb_cs(cs)+Rqb_qbqb_cs(cs))/z
      enddo
      
      elseif ((j .eq. 0) .and. (k.eq.0)) then
      msq_qg=msq_cs(cs,+5,k)+msq_cs(cs,+4,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+2,k)+msq_cs(cs,+1,k)
     &      +msq_cs(cs,-5,k)+msq_cs(cs,-4,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-2,k)+msq_cs(cs,-1,k)
      msq_gq=msq_cs(cs,j,+5)+msq_cs(cs,j,+4)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+2)+msq_cs(cs,j,+1)
     &      +msq_cs(cs,j,-5)+msq_cs(cs,j,-4)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-2)+msq_cs(cs,j,-1)
      do cs=0,2
      xmsq=xmsq
     & +fx1(j)*fx2(k)*(-msq_qg*Pgq_g_cs(cs)-msq_gq*Pg_gq_cs(cs))
     & +fx1z(j)*fx2(k)*(+msq_qg*(Pgq_g_cs(cs)+Rgq_g_cs(cs))/z)
     & +fx1(j)*fx2z(k)*(+msq_gq*(Pg_gq_cs(cs)+Rg_gq_cs(cs))/z)
      xmsq=xmsq
     & +fx1(j)*fx2(k)*((-Pgg_g_cs(cs)-Pg_gg_cs(cs))*msq_cs(cs,j,k))
     & +fx1z(j)*fx2(k)*((Pgg_g_cs(cs)+Rgg_g_cs(cs))/z*msq_cs(cs,j,k))
     & +fx1(j)*fx2z(k)*((Pg_gg_cs(cs)+Rg_gg_cs(cs))/z*msq_cs(cs,j,k))
      enddo
      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
      msq_qq=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      do cs=0,2
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq_cs(cs,j,k)*(-Pgg_q_cs(cs)-Pg_qq_cs(cs))
     &                   -msq_qq*Pgqb_q_cs(cs))
     &  +fx1z(j)*fx2 (k)*((Pgg_q_cs(cs)+Rgg_q_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pgqb_q_cs(cs)+Rgqb_q_cs(cs))/z*msq_qq)
     &  +fx1 (j)*fx2z(k)*((Pg_qq_cs(cs)+Rg_qq_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pg_qg_cs(cs)+Rg_qg_cs(cs))/z*msq_cs(cs,0,0))
      enddo
      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
      msq_qq=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      do cs=0,2
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq_cs(cs,j,k)*(-Pgg_q_cs(cs)-Pg_qq_cs(cs))
     &                   -msq_qq*Pgq_qb_cs(cs))
     &  +fx1z(j)*fx2 (k)*((Pgg_q_cs(cs)+Rgg_q_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pgq_qb_cs(cs)+Rgq_qb_cs(cs))/z*msq_qq)
     &  +fx1 (j)*fx2z(k)*((Pg_qq_cs(cs)+Rg_qq_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pg_qg_cs(cs)+Rg_qg_cs(cs))/z*msq_cs(cs,0,0))
      enddo
      
      elseif ((k .eq. 0) .and. (j .gt. 0)) then
      msq_qq=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
      do cs=0,2
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq_cs(cs,j,k)*(-Pq_gg_cs(cs)-Pqq_g_cs(cs))
     &                   -msq_qq*Pq_gqb_cs(cs))
     &  +fx1 (j)*fx2z(k)*((Pq_gg_cs(cs)+Rq_gg_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pq_gqb_cs(cs)+Rq_gqb_cs(cs))/z*msq_qq)
     &  +fx1z(j)*fx2 (k)*((Pqq_g_cs(cs)+Rqq_g_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pqg_g_cs(cs)+Rqg_g_cs(cs))/z*msq_cs(cs,0,0))
      enddo
      
      elseif ((k .eq. 0) .and. (j .lt. 0)) then
      msq_qq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
      do cs=0,2
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq_cs(cs,j,k)*(-Pq_gg_cs(cs)-Pqq_g_cs(cs))
     &                   -msq_qq*Pqb_gq_cs(cs))
     &  +fx1 (j)*fx2z(k)*((Pq_gg_cs(cs)+Rq_gg_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pqb_gq_cs(cs)+Rqb_gq_cs(cs))/z*msq_qq)
     &  +fx1z(j)*fx2 (k)*((Pqq_g_cs(cs)+Rqq_g_cs(cs))/z*msq_cs(cs,j,k)
     &                   +(Pqg_g_cs(cs)+Rqg_g_cs(cs))/z*msq_cs(cs,0,0))
      enddo
      endif
      
      else
c--- SUM BY TOTAL MATRIX ELEMENTS: everything else
      if     ((j .gt. 0) .and. (k.lt.0)) then
      xmsq=xmsq
     & +fx1(j)*fx2(k)*(msq(j,k)*(one-Pqq_qb-Pq_qbqb)+msqv(j,k))
     & +fx1z(j)*fx2(k)*msq(j,k)*(Pqq_qb+Rqq_qb)/z
     & +fx1(j)*fx2z(k)*msq(j,k)*(Pq_qbqb+Rq_qbqb)/z

      elseif ((j .lt. 0) .and. (k.gt.0)) then
      xmsq=xmsq
     & +fx1(j)*fx2(k)*(msq(j,k)*(one-Pqbqb_q-Pqb_qq)+msqv(j,k))
     & +fx1z(j)*fx2(k)*msq(j,k)*(Pqbqb_q+Rqbqb_q)/z
     & +fx1(j)*fx2z(k)*msq(j,k)*(Pqb_qq+Rqb_qq)/z

      elseif ((j .eq. 0) .and. (k.eq.0)) then
      msq_qg=msq(+5,k)+msq(+4,k)+msq(+3,k)+msq(+2,k)+msq(+1,k)
     &      +msq(-5,k)+msq(-4,k)+msq(-3,k)+msq(-2,k)+msq(-1,k)
      msq_gq=msq(j,+5)+msq(j,+4)+msq(j,+3)+msq(j,+2)+msq(j,+1)
     &      +msq(j,-5)+msq(j,-4)+msq(j,-3)+msq(j,-2)+msq(j,-1)
      xmsq=xmsq
     & +fx1(j)*fx2(k)*(msq(j,k)*(one)
     &                -msq_qg*Pgq_g-msq_gq*Pg_gq+msqv(j,k))
     & +fx1z(j)*fx2(k)*(+(Pgq_g+Rgq_g)/z*msq_qg)
     & +fx1(j)*fx2z(k)*(+(Pg_gq+Rg_gq)/z*msq_gq)
      if (case .eq. 'Zbbbar') then
        do cs=0,2
        xmsq=xmsq
     &   +fx1(j)*fx2(k)*((-Pgg_g_cs(cs)-Pg_gg_cs(cs))*msq_cs(cs,j,k))
     &   +fx1z(j)*fx2(k)*((Pgg_g_cs(cs)+Rgg_g_cs(cs))/z*msq_cs(cs,j,k))
     &   +fx1(j)*fx2z(k)*((Pg_gg_cs(cs)+Rg_gg_cs(cs))/z*msq_cs(cs,j,k))
        enddo
      else
      xmsq=xmsq
     & +fx1(j)*fx2(k)*((-Pgg_g-Pg_gg)*msq(j,k))
     & +fx1z(j)*fx2(k)*((Pgg_g+Rgg_g)/z*msq(j,k))
     & +fx1(j)*fx2z(k)*((Pg_gg+Rg_gg)/z*msq(j,k))
      endif
      elseif (j .eq. 0) then
       if     (k.gt.0) then
       msq_qq=msq(-1,k)+msq(-2,k)+msq(-3,k)+msq(-4,k)+msq(-5,k)
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq(j,k)*(one-Pgg_q-Pg_qq)
     &                   -msq_qq*Pgq_q+msqv(j,k))
     &  +fx1z(j)*fx2 (k)*((Pgg_q+Rgg_q)/z*msq(j,k)
     &                   +(Pgq_q+Rgq_q)/z*msq_qq)
     &  +fx1 (j)*fx2z(k)*((Pg_qq+Rg_qq)/z*msq(j,k)
     &                   +(Pg_qg+Rg_qg)/z*msq(0,0))
       elseif (k.lt.0) then
       msq_qq=msq(+1,k)+msq(+2,k)+msq(+3,k)+msq(+4,k)+msq(+5,k)
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq(j,k)*(one-Pgg_q-Pg_qq)
     &                   -msq_qq*Pgq_q+msqv(j,k))
     &  +fx1z(j)*fx2 (k)*((Pgg_q+Rgg_q)/z*msq(j,k)
     &                   +(Pgq_q+Rgq_q)/z*msq_qq)
     &  +fx1 (j)*fx2z(k)*((Pg_qq+Rg_qq)/z*msq(j,k)
     &                   +(Pg_qg+Rg_qg)/z*msq(0,0))
       endif
      elseif (k .eq. 0) then
       if     (j.gt.0) then
       msq_qq=msq(j,-1)+msq(j,-2)+msq(j,-3)+msq(j,-4)+msq(j,-5)
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq(j,k)*(one-Pq_gg-Pqq_g)
     &                   -msq_qq*Pq_gq+msqv(j,k))
     &  +fx1 (j)*fx2z(k)*((Pq_gg+Rq_gg)/z*msq(j,k)
     &                   +(Pq_gq+Rq_gq)/z*msq_qq)
     &  +fx1z(j)*fx2 (k)*((Pqq_g+Rqq_g)/z*msq(j,k)
     &                   +(Pqg_g+Rqg_g)/z*msq(0,0))
       elseif (j.lt.0) then
       msq_qq=msq(j,+1)+msq(j,+2)+msq(j,+3)+msq(j,+4)+msq(j,+5)
       xmsq=xmsq     
     &  +fx1 (j)*fx2 (k)*(msq(j,k)*(one-Pq_gg-Pqq_g)
     &                   -msq_qq*Pq_gq+msqv(j,k))
     &  +fx1 (j)*fx2z(k)*((Pq_gg+Rq_gg)/z*msq(j,k)
     &                   +(Pq_gq+Rq_gq)/z*msq_qq)
     &  +fx1z(j)*fx2 (k)*((Pqq_g+Rqq_g)/z*msq(j,k)
     &                   +(Pqg_g+Rqg_g)/z*msq(0,0))
       endif
      endif
      
      endif
      
 20   continue
      enddo
      enddo
      
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

c--- make some cuts ....
      if (makecuts) then
c--- ... either using pre-prepared ones
        if (bbproc) then
          if (cuts(pjet)) then
            ncutzero=ncutzero+1
            goto 999
          endif
        else
          if (madejetcuts(p,pjet,jets) .eqv. .false.) then
            ncutzero=ncutzero+1
            goto 999
          endif
        endif
      endif
    
      virtint=flux*xjac*pswt*xmsq/BrnRat

      if (bin) then
      val=virtint*wgt/dfloat(itmx) 
      call nplotter(r,s,pjet,val,0)
      endif

      return

 999  continue
      ntotzero=ntotzero+1
      
      return
      end


