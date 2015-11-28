      block data wsalam1
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      data Q(-5)/+0.333333333333333d0/
      data Q(-4)/-0.666666666666667d0/
      data Q(-3)/+0.333333333333333d0/
      data Q(-2)/-0.666666666666667d0/
      data Q(-1)/+0.333333333333333d0/
      data Q(0)/+0d0/
      data Q(+1)/-0.333333333333333d0/
      data Q(+2)/+0.666666666666667d0/
      data Q(+3)/-0.333333333333333d0/
      data Q(+4)/+0.666666666666667d0/
      data Q(+5)/-0.333333333333333d0/
      data tau/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/
      end 


      subroutine coupling
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      character*4 part
      common/part/part
      character*7 pdlabel
      integer nproc
      common/nproc/nproc
      double precision BigA,aemmz,alphas,amz,cmass,bmass
      logical verbose
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
      common/pdlabel/pdlabel
      common/verbose/verbose
      common/em/aemmz
      common/couple/amz

C---Effective field theory approach: good for all except couplings
C---to b-quarks (can be added later); Valid for scales below t-mass
C---see Georgi--NPB 363 (1991) 301 
      BigA=pi*aemmz/(sqrt(2d0)*gf)
      xw=BigA/wmass**2
c      xw=half*(one-sqrt(one-four*BigA/wmass**2*(one-BigA/wmass**2)))
c      is equal to the above


      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=sqrt(gwsq)
      call couplz(xw)

      
      


      if ((nproc .eq. 80) .and. (pdlabel .eq. 'hmrs90b')) then 
c----for check with Nason
c---ZZ
      zmass=91.18d0
      xw=0.228d0
      wmass=sqrt(zmass**2*(1d0-xw))

      aemmz=1d0/128d0      
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=sqrt(gwsq)
      scale=zmass
      musq=scale**2
      write(6,*) 'zmass for Frix/Nas',zmass
      write(6,*) 'wmass for Frix/Nas',wmass
      write(6,*) 'gwsq for Frix/Nas',gwsq
      write(6,*) 'gw for Frix/Nas',gw
      write(6,*) 'xw for Frix/Nas',xw

      elseif((pdlabel .eq. 'hmrs90b').and.(nproc .eq. 60)) then
c----for check with Frixione WW
      zmass=91.17d0
      wmass=80.0d0      
      xw=1d0-(wmass/zmass)**2
      aemmz=1d0/128d0      
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=sqrt(gwsq)
      scale=wmass
      musq=scale**2
      call couplz(xw)
      write(6,*) 'zmass for Frix/Nas',zmass
      write(6,*) 'wmass for Frix/Nas',wmass
      write(6,*) 'gwsq for Frix/Nas',gwsq
      write(6,*) 'gw for Frix/Nas',gw
      write(6,*) 'xw for Frix/Nas',xw

      elseif((pdlabel .eq. 'hmrs90b').and.
     . ((nproc .eq. 70) .or. (nproc .eq. 75))) then
c--wz
      zmass=91.17d0
      wmass=80.0d0
      xw=1d0-(wmass/zmass)**2
      aemmz=1d0/128d0      
      gwsq=fourpi*aemmz/xw
      gw=sqrt(gwsq)
      esq=gwsq*xw
      scale=0.5d0*(wmass+zmass)
      musq=scale**2
      write(6,*) 'zmass for Frix/Nas',zmass
      write(6,*) 'wmass for Frix/Nas',wmass
      write(6,*) 'gwsq for Frix/Nas',gwsq
      write(6,*) 'gw for Frix/Nas',gw
      write(6,*) 'xw for Frix/Nas',xw
      endif
      
c      write(6,*) 'zmass',zmass
c      write(6,*) 'wmass',wmass
c      write(6,*) 'gwsq',gwsq
c      write(6,*) 'gw',gw
c      write(6,*) 'xw',xw
      write(6,*) '************** Electroweak parameters **************'
      write(6,*) '*                                                  *'
      write(6,75) 'zmass',zmass,'wmass',wmass
      write(6,75) 'zwidth',zwidth,'wwidth',wwidth
      write(6,75) 'gwsq',gwsq,'gw',gw
      write(6,75) 'xw',xw,'1/alpha',1d0/aemmz
      write(6,76) 'esq',esq
      write(6,*) '****************************************************'

   75 format(' *  ',a6,f13.8,7x,a7,f13.8,'  *')
   76 format(' *  ',a6,f13.8,27x,'  *')

      call pdfwrap(pdlabel)      

      cmass=sqrt(mcsq)
      bmass=sqrt(mbsq)
      musq=scale**2
      if (part .eq. 'lord') then
        as=alphas(abs(scale),amz,1)
      else
        as=alphas(abs(scale),amz,2)
      endif
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

      call couplz(xw)

***************************************

c--- if we're doing W + jets, automatically make the CKM matrix
c--- diagonal since we're not interested in these small effects   
      if (  (nproc .eq. 10) .or. (nproc .eq. 11) .or. (nproc .eq. 16)
     . .or. (nproc .eq. 22) .or. (nproc .eq. 27)) then
        Vud=1d0
        Vus=0d0
        Vub=0d0
        Vcd=0d0
        Vcs=1d0
        Vcb=0d0
      endif

      if (verbose) then
      write(6,*)
      write(6,*) '***************** CKM mixing matrix ****************'
      write(6,*) '*                                                  *'
      write(6,47) Vud,Vus,Vub
      write(6,48) Vcd,Vcs,Vcb
      write(6,*) '****************************************************'
      if (  (nproc .eq. 10) .or. (nproc .eq. 11) .or. (nproc .eq. 16)
     . .or. (nproc .eq. 22) .or. (nproc .eq. 27)) then
      write(6,*) '* Forced to be diagonal for simplicity in W + jets *'
      write(6,*) '****************************************************'
      endif
 47   format(' *      Vud=',g10.5,'Vus=',g10.5,'Vub=',g10.5,'  *')
 48   format(' *      Vcd=',g10.5,'Vcs=',g10.5,'Vcb=',g10.5,'  *')
      endif      

      if (verbose) then
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      if (scale .gt. 0d0) then
      write(6,49) 'alpha_s (scale)',gsq/fourpi
      write(6,49) 'alpha_s (zmass)',amz
      else
      write(6,*) '*  Dynamic scale - alpha_s changed event-by-event  *'
      write(6,49) 'alpha_s (zmass)',amz
      endif
      write(6,*) '****************************************************'
 49   format(' *  ',a20,f12.8,16x,'*')
      endif
      
      return
      end




