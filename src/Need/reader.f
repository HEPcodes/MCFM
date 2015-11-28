      subroutine reader
      implicit none  
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'zerowidth.f'
      include 'new_pspace.f'
      include 'impsample.f'
      include 'cutoff.f'
      include 'clustering.f'
      integer ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,nmin,nmax
      integer nproc
      double precision alphas
      double precision sqrts,Rcut
      double precision cmass,bmass
      character*7 pdlabel
      character*9 runstring
      character*4 part
      logical verbose,makecuts,dryrun
      logical msbar,spira
      double precision bbsqmin,bbsqmax,wsqmin,wsqmax,rtsmin
      double precision mbbmin,mbbmax,Mwmin,Mwmax,deltam
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/nproc/nproc
      common/spira/spira
      common/ranno/idum
      common/part/part
      common/runstring/runstring
      common/Rcut/Rcut


      common/limits/bbsqmin,bbsqmax,wsqmin,wsqmax
      common/verbose/verbose
      common/makecuts/makecuts
      common/dryrun/dryrun
      common/energy/sqrts
      common/qmass/cmass,bmass
      common/deltam/deltam

      common/density/ih1,ih2
      common/pdlabel/pdlabel

      common/rtsmin/rtsmin
      common/nmin/nmin
      common/nmax/nmax
      common/msbar/msbar
      
 
      data new_pspace/.false./
      data impsample/.false./
      data msbar/.false./

      open(unit=20,file='options.DAT',status='old',err=999)
      read(20,*) nproc
      write(6,*) 'nproc=',nproc
      read(20,*) part
      write(6,*) 'part=',part
      read(20,*) runstring
      write(6,*) 'runstring=',runstring
      read(20,*) verbose
      write(6,*) 'verbose=',verbose
      read(20,*) sqrts
      if (verbose) write(6,*) 'sqrts=',sqrts
      read(20,*) ih1
      if (verbose) write(6,*) 'ih1',ih1
      read(20,*) ih2
      if (verbose) write(6,*) 'ih2',ih2
      read(20,*) pdlabel
      if (verbose) write(6,*) 'pdlabel ',pdlabel
      read(20,*) hmass
      if (verbose) write(6,*) 'hmass',hmass
      read(20,*) scale
      if (verbose) write(6,*) 'scale',scale
        read(20,*) Mwmin
      if (verbose) write(6,*) 'm34min',Mwmin
      read(20,*) Mwmax 
      if (verbose) write(6,*) 'm34max',Mwmax
      read(20,*) mbbmin
      if (verbose) write(6,*) 'm56min',mbbmin
      read(20,*) mbbmax 
      if (verbose) write(6,*) 'm56max',mbbmax
      read(20,*) rtsmin
      if (verbose) write(6,*) 'rtsmin',rtsmin
      read(20,*) zerowidth
      if (verbose) write(6,*) 'zerowidth',zerowidth
      read(20,*) makecuts
      if (verbose) write(6,*) 'makecuts',makecuts
      read(20,*) Rcut
      if (verbose) write(6,*) 'Rcut',Rcut
      read(20,*) itmx1
      if (verbose) write(6,*) 'itmx1',itmx1
      read(20,*) ncall1
      if (verbose) write(6,*) 'ncall1',ncall1
      read(20,*) itmx2
      if (verbose) write(6,*) 'itmx2',itmx2
      read(20,*) ncall2
      if (verbose) write(6,*) 'ncall2',ncall2
      read(20,*) idum
      if (verbose) write(6,*) 'idum',idum
      read(20,*) realwt
      if (verbose) write(6,*) 'realwt',realwt
      read(20,*) cutoff
      if (verbose) write(6,*) 'cutoff',cutoff
      read(20,*) dryrun
      if (verbose) write(6,*) 'dryrun',dryrun
      read(20,*) debug
      if (verbose) write(6,*) 'debug',debug
      close(unit=20)

c      if (debug) then
      open(unit=21,file='debug.DAT',status='old',err=999)
      read(21,*) new_pspace
      if (verbose) write(6,*) 'new_pspace',new_pspace
      read(21,*) virtonly
      if (verbose) write(6,*) 'virtonly',virtonly
      read(21,*) realonly
      if (verbose) write(6,*) 'realonly',realonly
      read(21,*) spira
      if (verbose) write(6,*) 'spira',spira
      read(21,*) noglue
      if (verbose) write(6,*) 'noglue',noglue
      read(21,*) ggonly
      if (verbose) write(6,*) 'ggonly',ggonly
      read(21,*) ggexcl
      if (verbose) write(6,*) 'ggexcl',ggexcl
      read(21,*) nmin
      if (verbose) write(6,*) 'nmin',nmin
      read(21,*) nmax
      if (verbose) write(6,*) 'nmax',nmax
      read(21,*) clustering
      if (verbose) write(6,*) 'clustering',clustering
      close(unit=21)
c      endif

      if (ggonly .and. ggexcl) then
        write(6,*) 'ggonly and ggexcl BOTH .true. - pick one!'
        stop
      endif

      call chooser(nproc)

c-----initialize various quantities

c---initialize masses for alpha_s routine
      cmass=sqrt(mcsq)
      bmass=sqrt(mbsq)


      bbsqmin=mbbmin**2
      bbsqmax=mbbmax**2

      wsqmin=Mwmin**2
      wsqmax=Mwmax**2


c-----stange-marciano formula for resolution
      deltam=sqrt(0.64d0*hmass+0.03d0**2*hmass**2)

      if (verbose) write(6,*) 'delta m',deltam

      return
 999  continue
      write(6,*) 'Problem reading options.DAT'
      stop

      end

