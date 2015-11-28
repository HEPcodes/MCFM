      subroutine coefsdk(s12,mtsq,ct,cv,c1)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      double precision cv,ct,c1,s12,mtsq,ddilog,rsq,omrsq,eta,Kfun,
     . wlog,rlog

      if (scheme .eq.'dred') then
C------        eta=0 4d-hel
         eta=0d0
      elseif (scheme .eq. 'tH-V') then
C------       eta=1 t'Hooft Veltman
         eta=1d0
      endif

      rsq=s12/mtsq
      omrsq=1d0-rsq
      wlog=dlog(omrsq)
      rlog=dlog(rsq)

c--- epinv here stands for (4*pi*musq/mtsq)^ep/Gamma(1-ep)/ep
      ct=epinv2*epinv
     . +epinv*(2.5d0-2d0*wlog)
     . +25d0/4d0+0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog
     . +0.5d0/omrsq+2d0*ddilog(omrsq)-5d0*pisqo6
     . -5d0*wlog+2d0*wlog**2+eta/2d0

C---- this routine has been constructed from 
C---- %\cite{Gottschalk:1980rv}
C---- \bibitem{Gottschalk:1980rv}
C---- T.~Gottschalk,
C---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
C---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
C---- %%CITATION = PHRVA,D23,56;%%
C----- Adapted from Eqs.(A8,A9)

      Kfun=1d0/rsq*wlog
      c1=2d0*Kfun      
      cv=-epinv*epinv2
     . -epinv*(2.5d0-2d0*wlog)
     . -0.5d0*(11d0+eta)-pisqo6-2d0*ddilog(rsq)
     .  +3d0*wlog-2d0*wlog**2-Kfun
      
      return
      end

