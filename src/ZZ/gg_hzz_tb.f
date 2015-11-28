      subroutine gg_hZZ_tb(p,msq)
      implicit none
c--- Author: J. M. Campbell, September 2013
c--- Matrix element squared for gg -> H -> ZZ signal process
c--- The exact result for massive bottom and top quark loops is included
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'qlfirst.f'
      integer h1,h2,h34,h56
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      double complex ggH_bquark(2,2,2,2),ggH_tquark(2,2,2,2),Ahiggs

      if (qlfirst) then
        qlfirst=.false. 
        call qlinit
      endif
      
      msq(:,:)=0d0
      
      call getggHZZamps(p,ggH_bquark,ggH_tquark)
      
      msqgg=0d0
      do h1=1,2
      h2=h1
      do h34=1,2
      do h56=1,2
      
c--- compute total Higgs amplitude
      AHiggs=
     &  +ggH_bquark(h1,h2,h34,h56)
     &  +ggH_tquark(h1,h2,h34,h56)

      msqgg=msqgg+cdabs(AHiggs)**2

      enddo
      enddo
      enddo

c--- overall factor extracted (c.f. getggHZZamps.f)
      fac=avegg*V*(4d0*esq*gsq/(16d0*pisq)*esq)**2
      
      msq(0,0)=msqgg*fac

      return
      end
      
      
