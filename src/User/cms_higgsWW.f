
      subroutine cms_higgsWW(p,failed_cuts) 
      implicit none
      include 'constants.f' 
      include 'process.f'
      double precision p(mxpart,4)
      logical failed_cuts 
      double precision pt,etarap,m45
      double precision pt_l_hard,pt_l_soft,pttwo
      double precision pts,pth,phimax,mllmax,etmiss,etmiss_min
      double precision et_vec(4),eta_max,r2,delphi,eta_max_h
      double precision eta_max_s,eta_hard,eta_soft
      integer i 
      logical first 
      data first/.true./

      !-- cut params change to change cuts
      pts=25d0
      pth=30d0
      mllmax=50d0
      phimax=60d0
      etmiss_min=20d0
      eta_max_h=2.5d0
      eta_max_s=2.5d0
      
      if(first) then 
         first=.false.
      write(6,*)  '**************** Higgs Search cuts  ****************'
      write(6,*)  '*                                                  *'
      write(6,99)  '*     pt(lep)_max>  ', pth, '                    *'
      write(6,99)  '*     pt(lep)_min >  ', pts, '                    *'
      write(6,99)  '*     mll  <  ', mllmax, '                     *'
      write(6,99)  '*     phi(lep,lep) <  ', phimax, '               *'
      write(6,99)  '*     |eta_l | <  ', eta_max_h, '                 *'
      write(6,99) '*       ET_miss > ',etmiss_min, '                *'
      write(6,*)  '****************************************************'
      endif


!---- Subroutine for calculting CMS Higgs (m_H=160) cuts namely 
!---- pt_l(max) > 30 pt_l(min) > 25 m_ll < 50 delta_phi_ll < 60 degrees 
!---- eta_l < 2.5, ET_miss = 20
!-----Inclusive in the jet 
! '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))' 
      if(case.ne.'WWqqbr') then 
         write(6,*)'Attempted to apply Higgs search cuts to process' 
         write(6,*)'Check Runstring ' 
         stop
         return 
      endif

      
      failed_cuts = .false. 
      do i=1,4
         et_vec(i)=0d0
      enddo

      if(pt(4,p).gt.pt(5,p)) then 
         pt_l_hard=pt(4,p) 
         eta_hard=etarap(4,p)       
         pt_l_soft=pt(5,p)
         eta_soft=etarap(5,p)
      else
         pt_l_hard=pt(5,p) 
         eta_hard=etarap(5,p) 
         pt_l_soft=pt(4,p) 
         eta_soft=etarap(4,p)
      endif

   
      
!---- pt cuts
      if((pt_l_hard.lt.pth).or.(pt_l_soft.lt.pts)) then 
         failed_cuts=.true. 
         return 
      endif

!Binoth et al cuts 
!      if((pt_l_hard.gt.50d0).or.(pt_l_hard.lt.35d0)) then 
!         failed_cuts=.true. 
!         return
!      elseif((pt_l_soft.gt.25d0).or.(pt_l_soft.lt.20d0)) then 
!         failed_cuts=.true.
!         return
!      endif
     
      if(dabs(etmiss(p,et_vec)).lt.etmiss_min) then 
         failed_cuts=.true.
         return 
      endif

!---- eta_cuts
      if((dabs(eta_hard).gt.eta_max_h).or.
     &     ((dabs(eta_soft).gt.eta_max_s))) then 
         failed_cuts=.true. 
!         write(6,*) eta_hard,eta_soft
         return
      endif

     
!--- m_ll cut 
      m45=0d0
      do i=1,4
         if(i.ne.4) then 
            m45=m45-(p(4,i)+p(5,i))**2
         else
            m45=m45+(p(4,i)+p(5,i))**2 
         endif
      enddo

      if(dsqrt(max(m45,0d0)).gt.mllmax) then 
         failed_cuts=.true.
         return 
      endif

!---- phi_ll 
!--- convert phi max to radians
      phimax=phimax/(360d0)*2d0*pi

       r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     .     /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1D0
      if (r2 .lt. -0.9999999D0) r2=-1D0
      delphi=dacos(r2)
      if(delphi.gt.phimax) then 
         failed_cuts=.true. 
         return 
      endif

 99   format(1x,a29,f6.2,a17)

      return 
      end subroutine
      



        
