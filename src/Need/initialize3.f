      subroutine initialize3
!******************************************************************
!     sets up masses and coupling constants for Helas3
!******************************************************************
      implicit none

      include 'coupl.inc'
      include 'qcdcouple.f'
      double precision gw_mcfm,xw_mcfm,gwsq_mcfm,esq_mcfm
      common/ewcouple/gw_mcfm,xw_mcfm,gwsq_mcfm,esq_mcfm

!-----------
! Begin Code
!-----------

      call coupsm(0)
      twidth=1.53428742d0

c      write(6,*) 'gw_mcfm',gw_mcfm
c      write(6,*) 'gwsq_mcfm',gwsq_mcfm
c      write(6,*) 'gg(1)',gg(1)
c      write(6,*) 'gg(2)',gg(2)
c      write(6,*) 'gwf(1)',gwf(1)
c      write(6,*) 'gwf(2)',gwf(2)
c      write(6,*) 'gw',gw
c      write(6,*) 'gw/sqrt(2d0)',gw/sqrt(2d0)
c      write(6,*) 'twidth',twidth
  

      gw_mcfm=0.649358543d0
      gwsq_mcfm=gw_mcfm**2
c      write(6,*) 'gwsq_mcfm',gwsq_mcfm
      gsq=gg(1)**2
c      write(6,*) 'gw_mcfm/gw',gw_mcfm/gw
c      write(6,*) 'mcfm:gsq',sqrt(gsq)/gg(1)
c      pause
      end
