      subroutine scaleset(rscalestart,fscalestart,p)
c--- wrapper routine to set a dynamic scale; please refer to individual
c--- routines for exact definitions of the scales.
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'frag.f'
      include 'nlooprun.f'
      include 'qcdcouple.f'
      double precision rscalestart,fscalestart,p(mxpart,4),mu0,
     & alphas,amz
      logical first
      common/couple/amz
      data first/.true./  
      save first
      
      if     (dynstring .eq. 'm(34)') then
        call scaleset_m34(p,mu0)
      elseif (dynstring .eq. 'm(345)') then
        call scaleset_m345(p,mu0)
      elseif (dynstring .eq. 'm(3456)') then
        call scaleset_m3456(p,mu0)
      elseif (dynstring .eq. 'sqrt(M^2+pt34^2)') then
        call scaleset_Msqpt34sq(p,mu0)
      elseif (dynstring .eq. 'sqrt(M^2+pt5^2)') then
        call scaleset_Msqpt5sq(p,mu0)
      elseif (dynstring .eq. 'pt(photon)') then
        call scaleset_ptphoton(p,mu0)
      elseif (dynstring .eq. 'HT') then
        call scaleset_HT(p,mu0)
      elseif (dynstring .eq. 'DDIS') then
        call scaleset_ddis(p,mu0)
      else
        write(6,*) 'Dynamic scale choice not recognized'
	write(6,*) '   dynamicscale = ',dynstring
	stop
      endif
      
      scale=rscalestart*mu0
      facscale=fscalestart*mu0
c--- piggy-back renomalization scale for fragmentation scale 
      frag_scale=scale
	    
      if (first) then
        write(6,*)
        write(6,*)'************** Dynamic scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*                 RENORMALIZATION                  *'
        write(6,45) ' mu_ren  =',rscalestart,dynstring
        write(6,*)'*                                                  *'
        write(6,*)'*                  FACTORIZATION                   *'
        write(6,45) ' mu_fac  =',fscalestart,dynstring
	if (frag) then
        write(6,*)'*                                                  *'
        write(6,*)'*                  FRAGMENTATION                   *'
        write(6,45) ' mu_frag =',rscalestart,dynstring
	endif
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.      
      endif
  
c--- catch absurdly large and small scales      
      if  (scale .gt. 3000d0) scale=3000d0
      if  (facscale .gt. 3000d0) facscale=3000d0
      if  (frag_scale .gt. 3000d0) frag_scale=3000d0
      if  (scale .lt. 1d0) scale=1d0
      if  (facscale .lt. 1d0) facscale=1d0
      if  (frag_scale .lt. 1d0) frag_scale=1d0

c--- run alpha_s
      as=alphas(scale,amz,nlooprun)
	
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2
      
      return

 45   format(1x,'* ',a15,f6.2,' x ',a24,' *')

      end
      
