      block data properties
      include 'masses.f'
      data md,mu,ms,mc,mb,mt/5d-3,5d-3,1d-1,0d0,0d0,175d0/
      data mel,mmu,mtau/0.510997d-3,0.105658389d0,1.777d0/
      data hmass,hwidth/100d0,0.0017d0/
      data wmass,wwidth/80.41d0,2.06d0/
      data zmass,zwidth/91.187d0,2.49d0/
      data twidth/1.55666215d0/
      data tauwidth/2.269d-12/
C---below are the values of the masses at about 100 GeV
c---fudged so as to get the Higgs-branching ratio right
c      data mtausq,mcsq,mbsq/3.1602d0,0.4d0,10.7d0/
      data mtausq,mcsq,mbsq/3.1602d0,0.4d0,21.344d0/
      end

      block data blckm
      implicit none
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
c---full matrix not used at present
c       data Vud,    Vus,    Vub,
c     &      Vcd,    Vcs,    Vcb
c     &    /0.9751d0,0.2220d0,0.0035d0,
c     &     0.2220d0,0.9743d0,0.0410d0/
c---Vub=Vcb=0

       data Vud,    Vus,    Vub,
     &      Vcd,    Vcs,    Vcb
     &    /0.975d0,0.22220486d0,0.0d0,
     &     0.22220486d0,0.975d0,0.0d0/

      end

      block data block8
      implicit none
      double precision aemmz
      common/em/aemmz
      data aemmz/7.7585538055706D-03/
c     value is 1/128.89
      end 







