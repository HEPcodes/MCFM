      subroutine coupling
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'verbose.f'
      include 'nlooprun.f'
      include 'process.f'
      include 'ewinput.f'
      include 'nflav.f'
      include 'b0.f'
      include 'dynamicscale.f'
      character*4 part,mypart
      common/part/part
      integer i
      double precision aemmz,alphas,amz,cmass,bmass
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      double precision xsq,topwidth
      character*3 inlabel(10)
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
      common/em/aemmz
      common/couple/amz
      common/mypart/mypart

c--- blank out labels that indicate input parameters
      do i=1,10
        inlabel(i)='   '
      enddo
      inlabel(3)='(+)'
      inlabel(4)='(+)'
      inlabel(8)='(+)'

      if (ewscheme .eq. -1) then
c--- This is the MCFM default, corresponding to an effective
c--- field theory approach valid for scales below the top-mass
C--- (see Georgi, Nucl. Phys. B 363 (1991) 301).
c--- There are 4 inputs here instead of the usual 3 ...
         Gf = Gf_inp
         aemmz  = aemmz_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(8)='   '
c--- ... and as result, both xw and mtop are derived
c--- (using wmass=zmass*dsqrt(rho)*cos(theta_w)
c---    and rho=1+3d0*aemmz/16d0/pi/xw*(mt/wmass)**2 )
         xw  = fourpi*aemmz/(8d0*wmass**2*Gf/rt2)
         mt  = dsqrt(16d0*pisq/3d0/rt2/Gf*(
     .          wmass**2/zmass**2/(1d0-xw)-1d0))

      elseif (ewscheme .eq. 0) then
c------------------------------------------------------------
c     option=0 : MadEvent default (= AlpGen with iewopt=2)
c------------------------------------------------------------

c-- equal to the input values
         xw  = xw_inp
         aemmz  = aemmz_inp
         zmass  = zmass_inp
         inlabel(7)='(+)'
         inlabel(6)='(+)'
         inlabel(1)='(+)'
c-- derived
         wmass  = zmass * dsqrt( One - xw )
         Gf = aemmz * Fourpi/xw/(8d0*wmass**2/Rt2)

      elseif (ewscheme .eq. 1) then
c-----------------------------------------------------
c     option=1 : LUSIFER and AlpGen (iewopt=3) default
c-----------------------------------------------------

c-- equal to the input values
         zmass  = zmass_inp
         wmass  = wmass_inp
         Gf = Gf_inp
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(5)='(+)'
c-- derived
         xw  = One-(wmass/zmass)**2
         aemmz  = Rt2*Gf*wmass**2*xw/pi

      elseif (ewscheme .eq. 2) then
c-------------------------------------------------------------------
c     option=2 : W and Z mass are derived from couplings
c-------------------------------------------------------------------

c-- equal to the input values
         Gf = Gf_inp
         aemmz  = aemmz_inp
         xw  = xw_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(7)='(+)'
c-- derived
         wmass  = dsqrt(aemmz*pi/xw/Gf/Rt2)
         zmass  = wmass/dsqrt(One-xw)

      elseif (ewscheme .eq. 3) then
c-----------------------------------------------------------------
c     option=3 : USER choice : you should know what you're doing!!
c-----------------------------------------------------------------
         Gf = Gf_inp
         aemmz  = aemmz_inp
         xw  = xw_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(7)='(+)'
         inlabel(1)='(+)'
         inlabel(2)='(+)'

      else
         write(6,*) 'ewscheme=',ewscheme,' is not a valid input.'
         stop
      endif

c--- Now set up the other derived parameters
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=dsqrt(gwsq)
      call couplz(xw)

c--- Calculate the appropriate Higgs vacuum expectation value.
c--- This vevsq is defined so that gwsq/(4*wmass**2)=Gf*rt2=1/vevsq
c--- (ie differs from definition in ESW)
      vevsq=1d0/rt2/Gf

c--- Set-up twidth, using LO formula except when including radiation in decay
      xsq=(wmass/mt)**2
      twidth=(gw/wmass)**2*mt**3/(64d0*pi)*(1d0-xsq)**2*(1d0+2d0*xsq)
      if ((part .eq. 'todk') .or. (mypart .eq. 'todk')) then
        twidth=twidth*topwidth(mt,wmass)
      endif

      write(6,*) '************** Electroweak parameters **************'
      write(6,*) '*                                                  *'
      write(6,75) 'zmass',inlabel(1),zmass,'wmass',inlabel(2),wmass
      write(6,75) 'zwidth',inlabel(3),zwidth,'wwidth',inlabel(4),wwidth
      write(6,76) 'Gf',inlabel(5),gf,'1/aemmz',inlabel(6),1d0/aemmz
      write(6,75) 'xw',inlabel(7),xw,'mtop',inlabel(8),mt
      write(6,75) 'gwsq',inlabel(9),gwsq,'esq',inlabel(10),esq
      write(6,77) 'top width',twidth
      write(6,78) 'mb',mb,'mc',mc
      write(6,*) '*                                                  *'
      write(6,*) '* Parameters marked (+) are input, others derived  *'
      write(6,*) '****************************************************'

   75 format(' * ',a6,a3,f13.7,3x,a7,a3,f12.7,'  *')
   76 format(' * ',a6,a3,d13.6,3x,a7,a3,f12.7,'  *')
   77 format(' * ',a9,f13.7,25x,'  *')
   78 format(' * ',a5,4x,f13.7,6x,a4,2x,f13.7,'  *')

c--- set up the beta-function
      b0=(xn*11d0-2d0*nflav)/6d0

c--- initialize the pdf set
      nlooprun=0
      call pdfwrap      

      cmass=dsqrt(mcsq)
      bmass=dsqrt(mbsq)
      musq=scale**2
 
c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfwrap already
      if (nlooprun .eq. 0) then
        if (part .eq. 'lord') then
          nlooprun=1
        else
          nlooprun=2
        endif
      endif

c--- initialize alpha_s
      as=alphas(abs(scale),amz,nlooprun)

      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

***************************************

c--- if we're doing W + jets, automatically make the CKM matrix
c--- diagonal since we're not interested in these small effects   
      if ((case .eq. 'W_1jet') .or. (case .eq. 'W_2jet')) then
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
      if ((case .eq. 'W_1jet') .or. (case .eq. 'W_2jet')) then
      write(6,*) '* Forced to be diagonal for simplicity in W + jets *'
      write(6,*) '****************************************************'
      endif
 47   format(' *      Vud=',g10.5,'Vus=',g10.5,'Vub=',g10.5,'  *')
 48   format(' *      Vcd=',g10.5,'Vcs=',g10.5,'Vcb=',g10.5,'  *')
      endif      

      if ((verbose) .and. (scale .gt. 0d0)) then      
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      if (dynamicscale .eqv. .false.) then
      write(6,49) 'alpha_s (scale)',gsq/fourpi
      write(6,49) 'alpha_s (zmass)',amz
      else
      write(6,*) '*  Dynamic scale - alpha_s changed event-by-event  *'
      write(6,49) 'alpha_s (zmass)',amz
      endif
      write(6,50) ' (using ',nlooprun,'-loop running of alpha_s)'  
      write(6,*) '****************************************************'
 49   format(' *  ',a20,f12.8,16x,'*')
 50   format(' *  ',6x,a8,i1,a25,8x,'*')
      endif
      
      return
      end


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



