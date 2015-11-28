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
      write(6,*) 'zmass',zmass
      write(6,*) 'wmass',wmass
      write(6,*) 'gwsq',gwsq
      write(6,*) 'gw',gw
      write(6,*) 'xw',xw


      if(scale.le.0d0) then
        if((nproc .ge. 70) .and. (nproc .le. 79)) then
          scale=0.5d0*(wmass+zmass)
        elseif((nproc .eq. 21) .or. (nproc .eq. 26)) then
          scale=0.5d0*(wmass+hmass)
        elseif((nproc .ge. 60) .and. (nproc .le. 69)) then
          scale=wmass
        elseif((nproc .ge. 80) .and. (nproc .le. 89)) then
          scale=zmass
        elseif((nproc .ge. 90) .and. (nproc .le. 96)) then
          scale=hmass
        elseif(nproc .eq. 111) then
          scale=hmass
        elseif ((nproc .ge. 150).or.(nproc .le. 152)) then
          scale=100d0
        elseif ((nproc .eq. 161)) then
          scale=100d0
        elseif ((nproc .eq. 171)) then
          scale=100d0
C appropriate scale to approximately include higher orders
        else
          write(*,*) 'Invalid input: please choose a scale'
          stop
        endif
      write(*,*) 'Scale automatically set to ',scale
      musq=scale**2
      endif 

      call pdfwrap(pdlabel)      

      cmass=sqrt(mcsq)
      bmass=sqrt(mbsq)
      musq=scale**2
      as=alphas(scale,amz,nloop)
      ason2pi=as/twopi
      gsq=fourpi*as

      if (verbose) write(6,*) 'alpha_s',gsq/fourpi,
     . ' corresponding to alpha_s(zmass)= ',amz


      call couplz(xw)

***************************************

      if (verbose) then
      write(6,*) 'Effective sin^2 theta_W xw',xw
      write(6,47) Vud,Vus,Vub
      write(6,48) Vcd,Vcs,Vcb
 47   format(' CKM:  Vud=',g10.5,'Vus=',g10.5,'Vub=',g10.5)
 48   format(' CKM   Vcd=',g10.5,'Vcs=',g10.5,'Vcb=',g10.5)
      endif      


      return
      end




