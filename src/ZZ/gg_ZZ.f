      subroutine gg_ZZ(p,msqgg)
      implicit none
c--- Author: J. M. Campbell, September 2013
c--- Matrix element squared for the process gg->ZZ
c--- Effects of massive quarks in the third generation may be included
c--- (default: included)
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'docheck.f'
      logical includegens1and2,includebottom,includetop
      integer h1,h2,h34,h56
      double precision p(mxpart,4),msqgg,fac
      double complex 
     & Mloop_uptype(2,2,2,2),Mloop_dntype(2,2,2,2),
     & Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2)

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.      
c--- set this to true to include massive bottom quark
      includebottom=.true.
c--- set this to true to include massive top quark
      includetop=.true.

c--- if set, performs check against numerical results at specific PS point
      docheck=.false.
      
c--- compute all gg->ZZ amplitudes      
      call getggZZamps(p,includegens1and2,includebottom,includetop,
     & Mloop_uptype,Mloop_dntype,Mloop_bquark,Mloop_tquark)
      
      msqgg=0d0
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
      
      msqgg=msqgg+cdabs(
     &  +2d0*Mloop_uptype(h1,h2,h34,h56)
     &  +2d0*Mloop_dntype(h1,h2,h34,h56)
     &      +Mloop_bquark(h1,h2,h34,h56)
     &      +Mloop_tquark(h1,h2,h34,h56))**2
      
      enddo
      enddo
      enddo
      enddo
      
c--- overall factor extracted (c.f. getggZZamps.f)
      fac=avegg*V*(4d0*esq*gsq/(16d0*pisq)*esq)**2
      
      msqgg=msqgg*fac
      
      return
      end
      
      


