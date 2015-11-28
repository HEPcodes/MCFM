      block data electroweak_input
************************************************************************
*     Calculation scheme for EW couplings                              *
************************************************************************
c
c     ewscheme=-1  : MCFM default 
c                    input values = Gf,alpha(m_Z),m_W,m_Z
c                    output values = sin^2(theta_W),mtop
c
c     ewscheme=0   : Old MadEvent default (= AlpGen with iewopt=2)
c                    input values = sin^2(theta_W),alpha(m_Z),m_Z
c                    output values = m_W,Gf.
c
c     ewscheme=1   : New Madevent default, "G_mu scheme"
c                    = LUSIFER and AlpGen (iewopt=3) defaults
c                    input values = G_F,m_Z,m_W
c                    output values = sin^2(theta_W),alpha(m_Z).
c
c     ewscheme=2   : input  values = G_F,sin^2(theta_W),alpha(m_Z)
c                    output values = m_W,m_Z.
c
c     ewscheme=3   : User choice. All parameters are left as they are
c                    input here. You have to know what you're doing.
c
      implicit none
      include 'ewinput.f'
      data ewscheme  / -1                  /   ! Chooses EW scheme
      data Gf_inp    / 1.16639d-5          /   ! G_F
      data aemmz_inp / 7.7585538055706d-03 /   ! alpha_EM(m_Z)=1/128.89
      data xw_inp    / 0.2312d0            /   ! sin^2(theta_W)
      data wmass_inp / 80.419d0            /   ! W mass
      data zmass_inp / 91.188d0            /   ! Z mass
      end
************************************************************************

      block data properties
      implicit none
      include 'masses.f'
      data md,mu,ms,mc,mb,mt/5d-3,5d-3,1d-1,1.5d0,4.2d0,175d0/
      data mel,mmu,mtau/0.510997d-3,0.105658389d0,1.777d0/
      data hmass,hwidth/100d0,0.0017d0/
      data wmass,wwidth/80.41d0,2.06d0/
      data zmass,zwidth/91.187d0,2.49d0/
      data twidth/1.55666215d0/
      data tauwidth/2.269d-12/
C---below are the values of the masses at about 100 GeV
c---fudged so as to get the Higgs-branching ratio right
c      data mtausq,mcsq,mbsq/3.1602d0,0.4d0,10.7d0/
      data mtausq,mcsq,mbsq/3.157729d0,2.25d0,17.64d0/
      end

      block data block_bH
      implicit none
      include 'runmb.f'
      include 'susycoup.f'
c-- whether or not to run mb (false makes no sense)      
      data runmb/.true./
c--susycoup is the deviation of Higgs coupling 
c-- from the standard model value
      data susycoup/1d0/
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

c       data Vud,    Vus,    Vub,
c     &      Vcd,    Vcs,    Vcb
c     &    /0.975d0,0.22220486d0,0.0d0,
c     &     0.22220486d0,0.975d0,0.0d0/

       data Vud,    Vus,    Vub,
     &      Vcd,    Vcs,    Vcb
     &    /1d0,0d0,0.0d0,
     &     0d0,1d0,0.0d0/

      end

c      block data block8
c      implicit none
c      double precision aemmz
c      common/em/aemmz
c--- This is the preferred value, aemmz = 1/128.89  
c      data aemmz/7.7585538055706D-03/
c--- Use this value for comparing with MADGRAPH, aemmz = 1/128
c      data aemmz/7.8125D-03/
c--- Use this value for aemmz = 1/137
c      data aemmz/7.2992901D-03/
c      end 

      block data ep
      implicit none
      include 'epinv.f'
      include 'epinv2.f'
      data epinv/ 1d3/
      data epinv2/1d3/
      end

