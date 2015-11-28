      logical function madewwcuts(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,dot,ptla,ptlb,yla,ylb,phila,
     . philb,ptmin1,ymax1,ptmin2,ymax2,mmin,missetmin,acut,pttwo
      parameter (ptmin1=10d0,ymax1=1.1d0,ptmin2=5d0,ymax2=2.5d0)
      parameter (mmin=10d0,missetmin=25d0,acut=150d0)
      
      madewwcuts=.false.
      
      ptla=pt(5,p)
      ptlb=pt(4,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(4,p))

      if ((ptla .gt. ptmin1 .and. yla .lt. ymax1 .and.
     .     ptlb .gt. ptmin2 .and. ylb .lt. ymax2)
     . .or.
     .    (ptlb .gt. ptmin1 .and. ylb .lt. ymax1 .and.
     .     ptla .gt. ptmin2 .and. yla .lt. ymax2)) then
        madewwcuts=.true.
      endif
      
      if (2d0*dot(p,5,4) .lt. mmin) then
        madewwcuts=.false.
      endif       
      
      if (pttwo(6,3,p) .lt. missetmin) then
        madewwcuts=.false.
      endif       
      
      phila=dasin(p(5,1)/ptla)
      philb=dasin(p(4,1)/ptlb)
      
      if (abs(phila-philb) .gt. pi*(acut/180d0)) then
         madewwcuts=.false.
      endif           
      
      return
      end
      
      logical function madewzcuts(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,ptla,ptlb,ptlc,yla,ylb,
     . ylc,ptmin1,ymax1,ptmin2,ptmin3,missetmin
      parameter (ptmin1=11d0,ymax1=1d0,ptmin2=7d0,ptmin3=5d0)
      parameter (missetmin=25d0)
      
      madewzcuts=.false.
      
      ptla=pt(5,p)
      ptlb=pt(6,p)
      ptlc=pt(4,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(6,p))
      ylc=abs(yrap(4,p))

      if ((ptla .gt. ptmin1 .and. yla .lt. ymax1 .and.
     .     ((ptlb .gt. ptmin2 .and. ptlc .gt. ptmin3) .or.
     .     (ptlc .gt. ptmin2 .and. ptlb .gt. ptmin3)))
     . .or.
     .    (ptlb .gt. ptmin1 .and. ylb .lt. ymax1 .and.
     .     ((ptla .gt. ptmin2 .and. ptlc .gt. ptmin3) .or.
     .     (ptlc .gt. ptmin2 .and. ptla .gt. ptmin3)))
     . .or.
     .    (ptlc .gt. ptmin1 .and. ylc .lt. ymax1 .and.
     .     ((ptla .gt. ptmin2 .and. ptlb .gt. ptmin3) .or.
     .     (ptlb .gt. ptmin2 .and. ptla .gt. ptmin3)))) then
        madewzcuts=.true.
      endif
      
      if (pt(3,p) .lt. missetmin) then
        madewzcuts=.false.
      endif       
      
      return
      end 
      
      logical function madetotcuts(p)
      implicit none
      include 'constants.f'
      integer nnproc
      double precision p(mxpart,4),pt,yrap,dot,
     . ptmin,ymax,missetmin,pttwo,mzmin,mzmax
      parameter (ptmin=20d0,ymax=2d0,missetmin=25d0)
      parameter (mzmin=75d0,mzmax=105d0)
      common/nproc/nnproc
      logical first 
      data first/.true./
      save first
      
      if (first) then
        write(*,*) 'SPECIAL CUTS OVER-RIDING THOSE ABOVE:'
        write(*,*) '   pt   > ',ptmin
        write(*,*) '  eta   < ',ymax
        write(*,*) 'miss Et > ',missetmin
        if (nnproc .ge. 70 .and. nnproc .le. 89) then
          write(*,*) mzmin,' < M_Z < ',mzmax
        endif
        first=.false.
      endif
      
      madetotcuts=.false.
      
      if     (nnproc .ge. 60 .and. nnproc .le. 69) then
        if (pt(5,p) .gt. ptmin .and. abs(yrap(5,p)) .lt. ymax .and.
     .      pt(4,p) .gt. ptmin .and. abs(yrap(4,p)) .lt. ymax .and.
     .      pttwo(6,3,p) .gt. missetmin) then
          madetotcuts=.true.
        endif
      elseif (nnproc .ge. 70 .and. nnproc .le. 74) then
        if (pt(5,p) .gt. ptmin .and. abs(yrap(5,p)) .lt. ymax .and.
     .      pt(6,p) .gt. ptmin .and. abs(yrap(6,p)) .lt. ymax .and.
     .      pt(4,p) .gt. ptmin .and. abs(yrap(4,p)) .lt. ymax .and.
     .      pt(3,p) .gt. missetmin .and. 2d0*dot(p,5,6) .gt. mzmin**2
     .      .and. 2d0*dot(p,5,6) .lt. mzmax**2) then
          madetotcuts=.true.
        endif
      elseif (nnproc .ge. 75 .and. nnproc .le. 79) then
        if (pt(5,p) .gt. ptmin .and. abs(yrap(5,p)) .lt. ymax .and.
     .      pt(6,p) .gt. ptmin .and. abs(yrap(6,p)) .lt. ymax .and.
     .      pt(3,p) .gt. ptmin .and. abs(yrap(3,p)) .lt. ymax .and.
     .      pt(4,p) .gt. missetmin .and. 2d0*dot(p,5,6) .gt. mzmin**2
     .      .and. 2d0*dot(p,5,6) .lt. mzmax**2) then
          madetotcuts=.true.
        endif
      elseif (nnproc .eq. 81) then
        if (pt(5,p) .gt. ptmin .and. abs(yrap(5,p)) .lt. ymax .and.
     .      pt(6,p) .gt. ptmin .and. abs(yrap(6,p)) .lt. ymax .and.
     .      pt(3,p) .gt. ptmin .and. abs(yrap(3,p)) .lt. ymax .and.
     .      pt(4,p) .gt. ptmin .and. abs(yrap(4,p)) .lt. ymax
     .      .and. 2d0*dot(p,5,6) .gt. mzmin**2
     .      .and. 2d0*dot(p,5,6) .lt. mzmax**2
     .      .and. 2d0*dot(p,3,4) .gt. mzmin**2
     .      .and. 2d0*dot(p,3,4) .lt. mzmax**2) then
          madetotcuts=.true.
        endif
      elseif (nnproc .eq. 82) then
        if (pt(5,p) .gt. ptmin .and. abs(yrap(5,p)) .lt. ymax .and.
     .      pt(6,p) .gt. ptmin .and. abs(yrap(6,p)) .lt. ymax .and.
     .      pttwo(3,4,p) .gt. missetmin
     .      .and. 2d0*dot(p,5,6) .gt. mzmin**2
     .      .and. 2d0*dot(p,5,6) .lt. mzmax**2
     .      .and. 2d0*dot(p,3,4) .gt. mzmin**2
     .      .and. 2d0*dot(p,3,4) .lt. mzmax**2) then
          madetotcuts=.true.
        endif
      else
        write(*,*) 'Special cuts not implemented for this case'
        stop
      endif
      
      return
      end
                 
      logical function madewzcuts2(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,ptla,ptlb,yla,ylb
      double precision ptmin,ymax
      parameter (ptmin=5d0,ymax=2.4d0)
      logical first
      data first/.true./
      save first
      
      if (first) then
        write(*,*) 'SPECIAL CUTS OVER-RIDING THOSE ABOVE:'
        write(*,*) '   pt(Z lept)   > ',ptmin
        write(*,*) '  eta(Z lept)   < ',ymax
        first=.false.
      endif
      
      madewzcuts2=.false.
      
      ptla=pt(5,p)
      ptlb=pt(6,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(6,p))

      if (ptla .gt. ptmin .and. yla .lt. ymax .and.
     .    ptlb .gt. ptmin .and. ylb .lt. ymax) then
        madewzcuts2=.true.
      endif
      
      return
      end 
      
      logical function madehiggscuts(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4)
      logical madehiggscutsww,madehiggscutswz
      logical madehiggscutszz,madehiggscutstt
      logical madehiggscutstautau
      integer nnproc
      common/nproc/nnproc
      
      madehiggscuts=.false.
      if (nnproc .eq. 61.or.nnproc .eq. 111) then
        madehiggscuts=madehiggscutsww(p)
      elseif (nnproc .eq. 70 .or. nnproc .eq. 71) then
        madehiggscuts=madehiggscutswz(p)
      elseif (nnproc .eq. 82) then
        madehiggscuts=madehiggscutszz(p)
      elseif (nnproc .eq. 151) then
        madehiggscuts=madehiggscutstt(p)
      elseif (nnproc .eq. 181) then
        madehiggscuts=madehiggscutstautau(p)
      else
      write(*,*) 'Special cuts not available for this process'
      stop
      endif
      
      return
      end      
      
      logical function madehiggscutsww(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,dot,ptla,ptlb,yla,ylb
      double precision ptemin,yemax,ptmu1min,ymu1max,ptmu2min,ymu2max
      double precision mmin,missetmin,pttwo,switch,rn
      double precision fphi,ftheta,coslpairet,mtsqlet,deltar
      double precision maxphill,maxthetall
      double precision minptll,maxcosllet,minmtlet,maxmee,maxmemu
      double precision mincosdd,maxcosdd,costhdd
      double precision ptjetveto,yjetveto,deltarmin
      integer comb
c --- parameters for basic cuts (10)
      parameter (ptemin=10d0,yemax=1.5d0,ptmu1min=10d0,ymu1max=1.5d0)
      parameter (ptmu2min=5d0,ymu2max=1.5d0,mmin=10d0,missetmin=10d0)
      parameter (deltarmin=0.4d0)
c --- parameters for extended cuts (11)-(16)
      parameter (maxphill=160d0,maxthetall=160d0)
      parameter (minptll=20d0,maxcosllet=0.5d0,minmtlet=20d0)
      parameter (maxmee=78d0,maxmemu=110d0)
      parameter (mincosdd=-0.3d0,maxcosdd=0.8d0)
      parameter (ptjetveto=95d0,yjetveto=3d0)
      
      madehiggscutsww=.false.
      
      ptla=pt(5,p)
      ptlb=pt(4,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(4,p))

c generate random number to fake ee/emu/mue/mumu with equal probability
      switch=rn(1)
      if      (switch.le.0.25d0) then
        comb=1    ! we have ee
      elseif (switch.gt.0.25d0.and.switch.le.0.5d0) then
        comb=2    ! we have emu
      elseif (switch.gt.0.5d0.and.switch.le.0.75d0) then
        comb=3    ! we have mue
      elseif (switch.gt.0.75d0.and.switch.le.1d0) then
        comb=4    ! we have mumu
      endif     
      
      if      (comb .eq. 1) then
c do e- e+
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlb .gt. ptemin .and. ylb .lt. yemax)
     .    ) then
        madehiggscutsww=.true.
        endif
       elseif (comb .eq. 2) then
c do e- mu+
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlb .gt. ptmu1min .and. ylb .lt. ymu1max)
     .    ) then
        madehiggscutsww=.true.
        endif
       elseif (comb .eq. 3) then
c do mu- e+
        if ((ptlb .gt. ptemin .and. ylb .lt. yemax .and.
     .     ptla .gt. ptmu1min .and. yla .lt. ymu1max)
     .    ) then
        madehiggscutsww=.true.
        endif
       elseif (comb .eq. 4) then
c do mu- mu+
        if ((ptla .gt. ptmu1min .and. yla .lt. ymu1max .and.
     .     ptlb .gt. ptmu2min .and. ylb .lt. ymu2max) .or.
     .      (ptlb .gt. ptmu1min .and. ylb .lt. ymu1max .and.
     .     ptla .gt. ptmu2min .and. yla .lt. ymu2max)
     .    )   then
        madehiggscutsww=.true.
        endif
      endif
      
      if (2d0*dot(p,5,4) .lt. mmin**2) then
        madehiggscutsww=.false.
      endif       
      
      if (pttwo(6,6,p) .lt. missetmin) then
        madehiggscutsww=.false.
      endif       

      if (min(deltar(5,7,p),deltar(4,7,p)) .lt. deltarmin) then
        madehiggscutsww=.false.
      endif       

c --- this completes the cuts of eq.(10)
      
      if (180d0/pi*fphi(5,4,p) .gt. maxphill) then
        madehiggscutsww=.false.
      endif       
       if (180d0/pi*ftheta(5,4,p) .gt. maxthetall) then
        madehiggscutsww=.false.
      endif       
     
c --- this completes the cuts of eq.(11)
      
      if (pttwo(5,4,p) .lt. minptll) then
        madehiggscutsww=.false.
      endif       
      if (coslpairet(5,4,6,3,p) .gt. maxcosllet) then
        madehiggscutsww=.false.
      endif       
      if (min(mtsqlet(5,6,3,p),mtsqlet(4,6,3,p)) .lt. minmtlet**2) then
        madehiggscutsww=.false.
      endif       

c --- this completes the cuts of eq.(12)

      if     ((comb .eq. 1 .or. comb .eq. 4)
     .     .and. 2d0*dot(p,5,4) .gt. maxmee**2) then
        madehiggscutsww=.false.
      elseif ((comb .eq. 2 .or. comb .eq. 3)
     .     .and. 2d0*dot(p,5,4) .gt. maxmemu**2) then
        madehiggscutsww=.false.
      endif       
            
c --- this completes the cuts of eq.(14)

      call dittdrein(p,5,4,costhdd)
      if (costhdd .lt. mincosdd .or. costhdd .gt. maxcosdd) then
        madehiggscutsww=.false.
      endif       
      
c --- this completes the cuts of eq.(15)

      if (pt(7,p) .gt. ptjetveto .and.
     .    abs(yrap(7,p)) .lt. yjetveto) then
        madehiggscutsww=.false.
      endif       

c --- this completes the cuts of eq.(16)

      return
      end
      
      logical function madehiggscutswz(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,dot,ptla,ptlb,ptlc,yla,ylb,
     . ylc,ptemin,yemax,ptmu1min,ymu1max,ptmu2min,ymu2max
      double precision mmin,missetmin,pttwo,switch,rn
      double precision fphi,ftheta,coslpairet,mtsqlet,deltar
      double precision maxphill,maxthetall
      double precision minptll,maxcosllet,minmtlet,maxmee,maxmemu
      double precision mincosdd,maxcosdd,costhdd
      double precision ptjetveto,yjetveto,deltarmin
      double precision ptdecmin,ydecmax
      integer l1,l2,m
      integer nnproc
      logical sameflav
      common/nproc/nnproc
c --- parameters for basic cuts (10)
      parameter (ptemin=10d0,yemax=1.5d0,ptmu1min=10d0,ymu1max=1.5d0)
      parameter (ptmu2min=5d0,ymu2max=1.5d0,mmin=10d0,missetmin=10d0)
      parameter (deltarmin=0.4d0)
c --- parameters for extended cuts (11)-(16)
      parameter (maxphill=160d0,maxthetall=160d0)
      parameter (minptll=20d0,maxcosllet=0.5d0,minmtlet=20d0)
      parameter (maxmee=78d0,maxmemu=110d0)
      parameter (mincosdd=-0.3d0,maxcosdd=0.8d0)
      parameter (ptjetveto=95d0,yjetveto=3d0)
c--- parameters for missing detector coverage
      parameter (ptdecmin=10d0,ydecmax=1.5d0)      
      madehiggscutswz=.false.
            
      ptla=pt(5,p)
      ptlb=pt(6,p)
      ptlc=pt(4,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(6,p))
      ylc=abs(yrap(4,p))
      
c first need to verify that all 3 leptons are not seen

      m=0
      if (ptla .lt. ptdecmin .or. yla .gt. ydecmax) m=4
      if (ptlb .lt. ptdecmin .or. ylb .gt. ydecmax) m=5
      if (ptlc .lt. ptdecmin .or. ylc .gt. ydecmax) m=7
      if (m .eq. 0) return
      
      if (nnproc .eq. 70) then
        madehiggscutswz=.true.
        return
      endif
      
      if (m .eq. 4) then
        l1=5
        l2=7
      else
        l1=5
        l2=12-m
      endif

      ptla=pt(l1,p)
      ptlb=pt(l2,p)
      yla=abs(yrap(l1,p))
      ylb=abs(yrap(l2,p))
     
      sameflav=.true.
      if (l2. eq. 4) then
c --- we might not have same flavour, decide
        switch=rn(1)
        if (switch .le. 0.5d0) sameflav=.false.
      endif
                
      switch=rn(1)
      if (sameflav) then
c --- both e's
        if     (switch. le. 0.5d0 .and.
     .      ptla .gt. ptemin .and. yla .lt. yemax .and.
     .      ptlb .gt. ptemin .and. ylb .lt. yemax) then
          madehiggscutswz=.true.
c --- both mu's
        elseif (switch. gt. 0.5d0 .and.
     .      ((ptla .gt. ptmu1min .and. yla .lt. ymu1max .and.
     .        ptlb .gt. ptmu2min .and. ylb .lt. ymu2max) .or.
     .       (ptlb .gt. ptmu1min .and. ylb .lt. ymu1max .and.
     .        ptla .gt. ptmu2min .and. yla .lt. ymu2max))) then
          madehiggscutswz=.true.
        endif
       else
c --- e mu
        if     (switch. le. 0.5d0 .and.
     .      ptla .gt. ptemin .and. yla .lt. yemax .and.
     .      ptlb .gt. ptmu1min .and. ylb .lt. ymu1max) then
          madehiggscutswz=.true.
c --- mu e
        elseif (switch. gt. 0.5d0 .and.
     .      ptlb .gt. ptemin .and. ylb .lt. yemax .and.
     .      ptla .gt. ptmu1min .and. yla .lt. ymu1max) then
          madehiggscutswz=.true.
        endif        
       endif 
            
      if (2d0*dot(p,l1,l2) .lt. mmin**2) then
        madehiggscutswz=.false.
      endif       
      
      if (pttwo(3,m,p) .lt. missetmin) then
        madehiggscutswz=.false.
      endif       

c --- this completes the cuts of eq.(10)
      
      if (180d0/pi*fphi(l1,l2,p) .gt. maxphill) then
        madehiggscutswz=.false.
      endif       
      if (180d0/pi*ftheta(l1,l2,p) .gt. maxthetall) then
        madehiggscutswz=.false.
      endif       
     
      if (min(deltar(l1,7,p),deltar(l2,7,p)) .lt. deltarmin) then
        madehiggscutswz=.false.
      endif       

c --- this completes the cuts of eq.(11)
      
      if (pttwo(l1,l2,p) .lt. minptll) then
        madehiggscutswz=.false.
      endif       
      if (coslpairet(l1,l2,3,m,p) .gt. maxcosllet) then
        madehiggscutswz=.false.
      endif       
      if (min(mtsqlet(l1,3,m,p),mtsqlet(l2,3,m,p))
     .     .lt. minmtlet**2) then
        madehiggscutswz=.false.
      endif       

c --- this completes the cuts of eq.(12)
      
      if (sameflav .and. 2d0*dot(p,l1,l2) .gt. maxmee**2) then
        madehiggscutswz=.false.
      endif       
                        
      if (sameflav .neqv. .true. .and.
     . 2d0*dot(p,l1,l2) .gt. maxmemu**2) then
        madehiggscutswz=.false.
      endif       
                        
c --- this completes the cuts of eq.(14)

      call dittdrein(p,l1,l2,costhdd)
      if (costhdd .lt. mincosdd .or. costhdd .gt. maxcosdd) then
        madehiggscutswz=.false.
      endif       
      
c --- this completes the cuts of eq.(15)

      if (pt(7,p) .gt. ptjetveto .and.
     .    abs(yrap(7,p)) .lt. yjetveto) then
        madehiggscutswz=.false.
      endif       

c --- this completes the cuts of eq.(16)

      return
      end
      
      logical function madehiggscutszz(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,dot,ptla,ptlb,ptlc,ptld
      double precision yla,ylb,ylc,yld
      double precision ptemin,yemax,ptmu1min,ymu1max,ptmu2min,ymu2max
      double precision mmin,missetmin,pttwo,switch,rn
      double precision fphi,ftheta,coslpairet,mtsqlet,deltar
      double precision maxphill,maxthetall
      double precision minptll,maxcosllet,minmtlet,maxmee,maxmemu
      double precision mincosdd,maxcosdd,costhdd
      double precision ptjetveto,yjetveto,deltarmin
c --- parameters for basic cuts (10)
      parameter (ptemin=10d0,yemax=1.5d0,ptmu1min=10d0,ymu1max=1.5d0)
      parameter (ptmu2min=5d0,ymu2max=1.5d0,mmin=10d0,missetmin=10d0)
      parameter (deltarmin=0.4d0)
c --- parameters for extended cuts (11)-(16)
      parameter (maxphill=160d0,maxthetall=160d0)
      parameter (minptll=20d0,maxcosllet=0.5d0,minmtlet=20d0)
      parameter (maxmee=78d0,maxmemu=110d0)
      parameter (mincosdd=-0.3d0,maxcosdd=0.8d0)
      parameter (ptjetveto=95d0,yjetveto=3d0)
      
      madehiggscutszz=.false.
      
      ptla=pt(5,p)
      ptlb=pt(4,p)
      ptlc=pt(6,p)
      ptld=pt(3,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(4,p))
      ylc=abs(yrap(6,p))
      yld=abs(yrap(3,p))

c generate random number to fake eevv/mumuvv with equal probability
      switch=rn(1)
      
      if      (switch.le.0.5d0) then
c do e- e+ nu nubar
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlc .gt. ptemin .and. ylc .lt. yemax)
     .    ) then
        madehiggscutszz=.true.
        endif
       elseif (switch.gt.0.5d0.and.switch.le.1d0) then
c do mu- mu+ nu nubar
        if ((ptla .gt. ptmu1min .and. yla .lt. ymu1max .and.
     .     ptlc .gt. ptmu2min .and. ylc .lt. ymu2max) .or.
     .      (ptla .gt. ptmu2min .and. yla .lt. ymu2max .and.
     .     ptlc .gt. ptmu1min .and. ylc .lt. ymu1max)
     .    ) then
        madehiggscutszz=.true.
        endif
      endif
      
      if (2d0*dot(p,5,6) .lt. mmin**2) then
        madehiggscutszz=.false.
      endif       
      
      if (pttwo(3,4,p) .lt. missetmin) then
        madehiggscutszz=.false.
      endif       
      
      if (min(deltar(5,7,p),deltar(6,7,p)) .lt. deltarmin) then
        madehiggscutszz=.false.
      endif       

c --- this completes the cuts of eq.(10)
      
      if (180d0/pi*fphi(5,6,p) .gt. maxphill) then
        madehiggscutszz=.false.
      endif       
       if (180d0/pi*ftheta(5,6,p) .gt. maxthetall) then
        madehiggscutszz=.false.
      endif       
     
c --- this completes the cuts of eq.(11)
      
      if (pttwo(5,6,p) .lt. minptll) then
        madehiggscutszz=.false.
      endif       
      if (coslpairet(5,6,3,4,p) .gt. maxcosllet) then
        madehiggscutszz=.false.
      endif       
      if (min(mtsqlet(5,3,4,p),mtsqlet(6,3,4,p)) .lt. minmtlet**2) then
        madehiggscutszz=.false.
      endif       

c --- this completes the cuts of eq.(12)

      if (2d0*dot(p,5,6) .gt. maxmee**2) then
        madehiggscutszz=.false.
      endif       
         
c --- this completes the cuts of eq.(14)

      call dittdrein(p,5,6,costhdd)
      if (costhdd .lt. mincosdd .or. costhdd .gt. maxcosdd) then
        madehiggscutszz=.false.
      endif       
      
c --- this completes the cuts of eq.(15)

      if (pt(7,p) .gt. ptjetveto .and.
     .    abs(yrap(7,p)) .lt. yjetveto) then
        madehiggscutszz=.false.
      endif       

c --- this completes the cuts of eq.(16)

      return
      end
      
      logical function madehiggscutstt(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pt,yrap,dot,ptla,ptlb
      double precision yla,ylb
      double precision ptemin,yemax,ptmu1min,ymu1max,ptmu2min,ymu2max
      double precision mmin,missetmin,pttwo,switch,rn
      double precision fphi,ftheta,coslpairet,mtsqlet,deltarpjet
      double precision maxphill,maxthetall
      double precision minptll,maxcosllet,minmtlet,maxmee,maxmemu
      double precision mincosdd,maxcosdd,costhdd
      double precision ptjetveto1,yjetveto1,deltarmin
      double precision ptjetveto2,yjetveto2
      double precision ptb,ptbb,yb,ybb,btag,bbtag,rn1,rn2
      integer comb
c --- for jet clustering
      integer njet,parts(mxpart)
      double precision pjet(mxpart,4),ptjet
c --- parameters for basic cuts (10)
      parameter (ptemin=10d0,yemax=1.5d0,ptmu1min=10d0,ymu1max=1.5d0)
      parameter (ptmu2min=5d0,ymu2max=1.5d0,mmin=10d0,missetmin=10d0)
      parameter (deltarmin=0.4d0)
c --- parameters for extended cuts (11)-(16)
      parameter (maxphill=160d0,maxthetall=160d0)
      parameter (minptll=20d0,maxcosllet=0.5d0,minmtlet=20d0)
      parameter (maxmee=78d0,maxmemu=110d0)
      parameter (mincosdd=-0.3d0,maxcosdd=0.8d0)
      parameter (ptjetveto1=95d0,yjetveto1=3d0)
      parameter (ptjetveto2=50d0,yjetveto2=3d0)
      
      madehiggscutstt=.false.

c --- cluster bb jets
      call genclust(p,5,6,0.4d0,njet,pjet,parts)
c --- with no clustering 
c      do i=1,4
c        do njet=1,2
c        pjet(njet,i)=p(njet+3,i)
c        enddo
c      enddo
c      njet=2
             
      ptb=ptjet(5,p,pjet)
      yb=dabs(yrap(5,pjet))
      if (njet .eq. 2) then
        ptbb=ptjet(6,p,pjet)
        ybb=dabs(yrap(6,pjet))
      else
        ptbb=0d0
        ybb=0d0
      endif
      
c --- now re-label the momenta so that they can be used in
c --- the existing routines and we (mostly) don't require p4,p6

c---------------??????????????????????????????????????????
c      do i=1,4
c        p(5,i)=q(7,i)
c        p(6,i)=q(8,i)
c      enddo
     
c --- we now have e-(p5),nubar(p6),nu(p3),e+(p4)     

      ptla=pt(5,p)
      ptlb=pt(4,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(4,p))

c generate random number to fake ee/emu/mue/mumu with equal probability
      switch=rn(1)
      if      (switch.le.0.25d0) then
        comb=1    ! we have ee
      elseif (switch.gt.0.25d0.and.switch.le.0.5d0) then
        comb=2    ! we have emu
      elseif (switch.gt.0.5d0.and.switch.le.0.75d0) then
        comb=3    ! we have mue
      elseif (switch.gt.0.75d0.and.switch.le.1d0) then
        comb=4    ! we have mumu
      endif     
      
      if      (comb .eq. 1) then
c do e- e+
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlb .gt. ptemin .and. ylb .lt. yemax)
     .    ) then
        madehiggscutstt=.true.
        endif
       elseif (comb .eq. 2) then
c do e- mu+
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlb .gt. ptmu1min .and. ylb .lt. ymu1max)
     .    ) then
        madehiggscutstt=.true.
        endif
       elseif (comb .eq. 3) then
c do mu- e+
        if ((ptlb .gt. ptemin .and. ylb .lt. yemax .and.
     .     ptla .gt. ptmu1min .and. yla .lt. ymu1max)
     .    ) then
        madehiggscutstt=.true.
        endif
       elseif (comb .eq. 4) then
c do mu- mu+
        if ((ptla .gt. ptmu1min .and. yla .lt. ymu1max .and.
     .     ptlb .gt. ptmu2min .and. ylb .lt. ymu2max) .or.
     .      (ptlb .gt. ptmu1min .and. ylb .lt. ymu1max .and.
     .     ptla .gt. ptmu2min .and. yla .lt. ymu2max)
     .    )   then
        madehiggscutstt=.true.
        endif
      endif
      
      if (2d0*dot(p,5,4) .lt. mmin**2) then
        madehiggscutstt=.false.
      endif       
      
      if (pttwo(6,3,p) .lt. missetmin) then
        madehiggscutstt=.false.
      endif       
      
      if (min(deltarpjet(5,1,p,pjet),deltarpjet(4,1,p,pjet))
     . .lt. deltarmin) then
        madehiggscutstt=.false.
      endif       
      if (min(deltarpjet(5,2,p,pjet),deltarpjet(4,2,p,pjet))
     . .lt. deltarmin .and. njet. eq. 2) then
        madehiggscutstt=.false.
      endif       

c --- this completes the cuts of eq.(10)
      
      if (180d0/pi*fphi(5,4,p) .gt. maxphill) then
        madehiggscutstt=.false.
      endif       
       if (180d0/pi*ftheta(5,4,p) .gt. maxthetall) then
        madehiggscutstt=.false.
      endif       
     
c --- this completes the cuts of eq.(11)
      
      if (pttwo(5,4,p) .lt. minptll) then
        madehiggscutstt=.false.
      endif       
      if (coslpairet(5,4,6,3,p) .gt. maxcosllet) then
        madehiggscutstt=.false.
      endif       
      if (min(mtsqlet(5,6,3,p),mtsqlet(4,6,3,p)) .lt. minmtlet**2) then
        madehiggscutstt=.false.
      endif       

c --- this completes the cuts of eq.(12)

      if     ((comb .eq. 1 .or. comb .eq. 4)
     .     .and. 2d0*dot(p,5,4) .gt. maxmee**2) then
        madehiggscutstt=.false.
      elseif ((comb .eq. 2 .or. comb .eq. 3)
     .     .and. 2d0*dot(p,5,4) .gt. maxmemu**2) then
        madehiggscutstt=.false.
      endif       
            
c --- this completes the cuts of eq.(14)

      call dittdrein(p,5,4,costhdd)
      if (costhdd .lt. mincosdd .or. costhdd .gt. maxcosdd) then
        madehiggscutstt=.false.
      endif       
      
c --- this completes the cuts of eq.(15)

      if (ptb .gt. ptbb) then
        if ((ptb .gt. ptjetveto1 .and. yb .lt. yjetveto1) .or.
     .      (ptbb .gt. ptjetveto2 .and. ybb .lt. yjetveto2)) then
          madehiggscutstt=.false.
        endif       
      else
        if ((ptbb .gt. ptjetveto1 .and. ybb .lt. yjetveto1) .or.
     .      (ptb .gt. ptjetveto2 .and. yb .lt. yjetveto2)) then
          madehiggscutstt=.false.
        endif       
      endif       

      btag=1.1d0*0.57d0*tanh(ptb/36.05d0)
      bbtag=1.1d0*0.57d0*tanh(ptbb/36.05d0)
     
      rn1=rn(1)
      rn2=rn(1)
      
      if ((rn1. le. btag) .or. 
     .    (rn2 .le. bbtag .and. njet .eq. 2)) then
        madehiggscutstt=.false.
      endif             
      
c --- this completes the cuts of eq.(16)

      return
      end

      logical function madehiggscutstautau(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pjet(mxpart,4),pt,yrap,dot
      double precision yla,ylb,ptla,ptlb
      double precision ptemin,yemax,ptmu1min,ymu1max,ptmu2min,ymu2max
      double precision mmin,missetmin,pttwo,switch,rn
      double precision fphi,ftheta,mtsqlet
      double precision maxphill,maxthetall
      double precision minptll,minmtlet,maxmee,maxmemu
      double precision mincosdd,maxcosdd,costhdd
      double precision ptjetveto1,yjetveto1
      double precision ptjetveto2,yjetveto2
      double precision ptb,ptbb,yb,ybb,btag,bbtag,rn1,rn2
      integer comb
c --- for jet clustering
      integer njet,parts(mxpart)
      double precision ptqfour,ptjet
c --- parameters for basic cuts (10)
      parameter (ptemin=10d0,yemax=1.5d0,ptmu1min=10d0,ymu1max=1.5d0)
      parameter (ptmu2min=5d0,ymu2max=1.5d0,mmin=10d0,missetmin=10d0)
c --- parameters for extended cuts (11)-(16)
      parameter (maxphill=160d0,maxthetall=160d0)
      parameter (minptll=20d0,minmtlet=20d0)
      parameter (maxmee=78d0,maxmemu=110d0)
      parameter (mincosdd=-0.3d0,maxcosdd=0.8d0)
      parameter (ptjetveto1=95d0,yjetveto1=3d0)
      parameter (ptjetveto2=50d0,yjetveto2=3d0)
      
c --- cluster bb jets
      call genclust(p,5,6,0.4d0,njet,pjet,parts)
c --- with no clustering 
c      do i=1,4
c        do njet=1,2
c        pjet(njet,i)=p(njet+3,i)
c        enddo
c      enddo
c      njet=2
             
      ptb=ptjet(5,p,pjet)
      yb=dabs(yrap(5,pjet))
      if (njet .eq. 2) then
        ptbb=ptjet(6,p,pjet)
        ybb=dabs(yrap(6,pjet))
      else
        ptbb=0d0
        ybb=0d0
      endif
      
c --- now re-label the momenta so that they can be used in
c --- the existing routines and we (mostly) don't require p4,p6

C----------------------------------
C?????????????????????????????????????
c      do i=1,4
c        p(5,i)=q(8,i)
c        p(6,i)=q(8,i)
c      enddo
C??????????????????????????????????????????     
c --- we now have e-(p4),nubar(p6),nu(p3),e+(p4)     

      ptla=pt(5,p)
      ptlb=pt(4,p)
      yla=abs(yrap(5,p))
      ylb=abs(yrap(4,p))

c generate random number to fake ee/emu/mue/mumu with equal probability
      switch=rn(1)
      if      (switch.le.0.25d0) then
        comb=1    ! we have ee
      elseif (switch.gt.0.25d0.and.switch.le.0.5d0) then
        comb=2    ! we have emu
      elseif (switch.gt.0.5d0.and.switch.le.0.75d0) then
        comb=3    ! we have mue
      elseif (switch.gt.0.75d0.and.switch.le.1d0) then
        comb=4    ! we have mumu
      endif     
      
      if      (comb .eq. 1) then
c do e- e+
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlb .gt. ptemin .and. ylb .lt. yemax)
     .    ) then
        madehiggscutstautau=.true.
        endif
       elseif (comb .eq. 2) then
c do e- mu+
        if ((ptla .gt. ptemin .and. yla .lt. yemax .and.
     .     ptlb .gt. ptmu1min .and. ylb .lt. ymu1max)
     .    ) then
        madehiggscutstautau=.true.
        endif
       elseif (comb .eq. 3) then
c do mu- e+
        if ((ptlb .gt. ptemin .and. ylb .lt. yemax .and.
     .     ptla .gt. ptmu1min .and. yla .lt. ymu1max)
     .    ) then
        madehiggscutstautau=.true.
        endif
       elseif (comb .eq. 4) then
c do mu- mu+
        if ((ptla .gt. ptmu1min .and. yla .lt. ymu1max .and.
     .     ptlb .gt. ptmu2min .and. ylb .lt. ymu2max) .or.
     .      (ptlb .gt. ptmu1min .and. ylb .lt. ymu1max .and.
     .     ptla .gt. ptmu2min .and. yla .lt. ymu2max)
     .    )   then
        madehiggscutstautau=.true.
        endif
      endif
      
      if (2d0*dot(p,5,4) .lt. mmin**2) then
        madehiggscutstautau=.false.
      endif       
      
      if (ptqfour(p,5,6,3,8) .lt. missetmin) then
        madehiggscutstautau=.false.
      endif       
      
c --- this completes the cuts of eq.(10)
      
      if (180d0/pi*fphi(5,4,p) .gt. maxphill) then
        madehiggscutstautau=.false.
      endif       
       if (180d0/pi*ftheta(5,4,p) .gt. maxthetall) then
        madehiggscutstautau=.false.
      endif       
     
c --- this completes the cuts of eq.(11)
      
      if (pttwo(5,4,p) .lt. minptll) then
        madehiggscutstautau=.false.
      endif       

      if (min(mtsqlet(5,5,4,p),mtsqlet(4,5,4,p)) .lt. minmtlet**2) then
        madehiggscutstautau=.false.
      endif       

c --- this completes the cuts of eq.(12)

      if     ((comb .eq. 1 .or. comb .eq. 4)
     .     .and. 2d0*dot(p,5,4) .gt. maxmee**2) then
        madehiggscutstautau=.false.
      elseif ((comb .eq. 2 .or. comb .eq. 3)
     .     .and. 2d0*dot(p,5,4) .gt. maxmemu**2) then
        madehiggscutstautau=.false.
      endif       
            
c --- this completes the cuts of eq.(14)

      call dittdrein(p,5,4,costhdd)
      if (costhdd .lt. mincosdd .or. costhdd .gt. maxcosdd) then
        madehiggscutstautau=.false.
      endif       
      
c --- this completes the cuts of eq.(15)

      if (ptb .gt. ptbb) then
        if ((ptb .gt. ptjetveto1 .and. yb .lt. yjetveto1) .or.
     .      (ptbb .gt. ptjetveto2 .and. ybb .lt. yjetveto2)) then
          madehiggscutstautau=.false.
        endif       
      else
        if ((ptbb .gt. ptjetveto1 .and. ybb .lt. yjetveto1) .or.
     .      (ptb .gt. ptjetveto2 .and. yb .lt. yjetveto2)) then
          madehiggscutstautau=.false.
        endif       
      endif       

      btag=1.1d0*0.57d0*tanh(ptb/36.05d0)
      bbtag=1.1d0*0.57d0*tanh(ptbb/36.05d0)
     
      rn1=rn(1)
      rn2=rn(1)
      
      if ((rn1. le. btag) .or. 
     .    (rn2 .le. bbtag .and. njet .eq. 2)) then
        madehiggscutstautau=.false.
      endif             
      
c --- this completes the cuts of eq.(16)

      return
      end
