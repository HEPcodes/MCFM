      subroutine reader_input
************************************************************************
*   Routine to read in the file input.DAT, which is a consolidated     *
*   form of all the input files and new to version 3.4 of MCFM         *
************************************************************************
      implicit none
      include 'constants.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'process.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'flags.f'
      include 'clustering.f'
      include 'anomcoup.f'
      include 'gridinfo.f'
      include 'verbose.f'
      include 'limits.f'
      include 'workdir.f'
      include 'jetcuts.f'
      include 'lhapdf.f'
      include 'alfacut.f'
      include 'pdlabel.f'
      include 'qcdcouple.f'
      include 'nlooprun.f'
      include 'initialscales.f'
      include 'stopscales.f'
      include 'vanillafiles.f'
      character*72 inputfile,getinput
      character*90 line
      character*4 part
      character*30 runstring
      integer nargs,lenocc,lenarg
      logical spira
      logical creatent,dswhisto,dryrun,makecuts
      integer nmin,nmax,n2,n3
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,origij
      integer NPTYPE,NGROUP,NSET
      double precision rtsmin,sqrts,factor
      double precision mbbmin,mbbmax,Mwmin,Mwmax
      double precision Rcut
      double precision leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut
      integer lbjscheme
      logical jetsopphem
      logical technicalincluded
      double precision ran2,randummy
      double precision cmass,bmass
      double precision mass2,width2,mass3,width3
      double precision amz,alphas
      
      common/couple/amz
      
      common/breit/n2,n3,mass2,width2,mass3,width3
      
      common/spira/spira
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin
 
      common/outputflags/creatent,dswhisto      

      common/nproc/nproc
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/ranno/idum
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts
      common/leptcuts/leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut,
     . lbjscheme,jetsopphem

      common/qmass/cmass,bmass

      common/origij/origij

      save /ranno/
      
      verbose=.true.

c--- work out the name of the input file and open it

      nargs=iargc()
      if (nargs .ge. 1) then
        call getarg(1,inputfile)
      else
        inputfile='input.DAT'
      endif
      
      lenarg=lenocc(inputfile)

      if ((lenarg.lt.4).or.(inputfile(lenarg-3:lenarg).ne.'.DAT')) then
        workdir=inputfile
c--- truncate if the directory / is included
        if (workdir(lenarg:lenarg) .eq. '/') then
           workdir(lenarg:lenarg)=' '
           lenarg=lenarg-1
        endif
        if (nargs .ge. 2) then
          call getarg(2,getinput)
          inputfile=workdir(1:lenarg)//'/'//getinput
        else  
          inputfile=workdir(1:lenarg)//'/input.DAT'
        endif
      else
        workdir=''
      endif
            
      write(6,*) 'Using input file named ',inputfile

      open(unit=20,file=inputfile,status='old',err=999)
      call checkversion(20,inputfile)

c--- read-in the user inputs

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- flags for the mode of MCFM   
      read(20,*) evtgen
      if (verbose) write(6,*) 'evtgen = ',evtgen
      read(20,*) creatent
      if (verbose) write(6,*) 'creatent = ',creatent
      read(20,*) skipnt
      if (verbose) write(6,*) 'skipnt = ',skipnt
      read(20,*) dswhisto
      if (verbose) write(6,*) 'dswhisto = ',dswhisto

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- general options 
      read(20,*) nproc
      if (verbose) write(6,*) 'nproc=',nproc
      read(20,*) part
      if (verbose) write(6,*) 'part=',part
      read(20,*) runstring
      if (verbose) write(6,*) 'runstring=',runstring
      read(20,*) sqrts
      if (verbose) write(6,*) 'sqrts=',sqrts
      read(20,*) ih1
      if (verbose) write(6,*) 'ih1',ih1
      read(20,*) ih2
      if (verbose) write(6,*) 'ih2',ih2
      read(20,*) hmass
      if (verbose) write(6,*) 'hmass',hmass
      read(20,*) scale
      if (verbose) write(6,*) 'scale',scale
      read(20,*) facscale
      if (verbose) write(6,*) 'facscale',facscale
      
      initrenscale_L=0d0
      initfacscale_L=0d0
      initrenscale_H=0d0
      initfacscale_H=0d0
c--- catch special scale choices for stop+b process
      if (((nproc .eq. 231) .or. (nproc .eq. 236)      
     . .or.(nproc .eq. 241) .or. (nproc .eq. 246)      
     . .or.(nproc .eq. 242) .or. (nproc .eq. 247)) .and.      
     .    (scale .eq. 0d0) .and. (facscale .eq. 0d0)) then
        read(20,*) initrenscale_L
        if (verbose) write(6,*) 'renscale_L',initrenscale_L
        read(20,*) initfacscale_L
        if (verbose) write(6,*) 'facscale_L',initfacscale_L
        read(20,*) initrenscale_H
        if (verbose) write(6,*) 'renscale_H',initrenscale_H
        read(20,*) initfacscale_H
        if (verbose) write(6,*) 'facscale_H',initfacscale_H
	scale=initrenscale_H
	facscale=initfacscale_H        
      endif

      read(20,*) dynamicscale
      if (verbose) write(6,*) 'dynamicscale',dynamicscale
      read(20,*) zerowidth
      if (verbose) write(6,*) 'zerowidth',zerowidth
      read(20,*) removebr
      if (verbose) write(6,*) 'removebr',removebr
      read(20,*) itmx1
      if (verbose) write(6,*) 'itmx1',itmx1
      read(20,*) ncall1
      if (verbose) write(6,*) 'ncall1',ncall1
      read(20,*) itmx2
      if (verbose) write(6,*) 'itmx2',itmx2
      read(20,*) ncall2
      if (verbose) write(6,*) 'ncall2',ncall2
      read(20,*) origij
      if (verbose) write(6,*) 'ij',origij
      read(20,*) dryrun
      if (verbose) write(6,*) 'dryrun',dryrun
      read(20,*) Qflag
      if (verbose) write(6,*) 'Qflag',Qflag
      read(20,*) Gflag
      if (verbose) write(6,*) 'Gflag',Gflag
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- heavy quark masses 
      read(20,*) mt
      if (verbose) write(6,*) 'mt',mt
      read(20,*) mb
      if (verbose) write(6,*) 'mb',mb
      read(20,*) mc
      if (verbose) write(6,*) 'mc',mc

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- pdf options 
      read(20,*) pdlabel
      if (verbose) write(6,*) 'pdlabel ',pdlabel
      read(20,*) NGROUP
      if (verbose) write(6,*) 'NGROUP=',NGROUP
      read(20,*) NSET
      if (verbose) write(6,*) 'NSET=',NSET
      read(20,*) PDFname
      if (verbose) write(6,*) 'PDFname=',PDFname
      read(20,*) PDFmember
      if (verbose) write(6,*) 'PDFmember=',PDFmember

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- jets and cuts options 
      read(20,*) Mwmin
      if (verbose) write(6,*) 'm34min',Mwmin
      read(20,*) Mwmax 
      if (verbose) write(6,*) 'm34max',Mwmax
      read(20,*) mbbmin
      if (verbose) write(6,*) 'm56min',mbbmin
      read(20,*) mbbmax 
      if (verbose) write(6,*) 'm56max',mbbmax
      read(20,*) inclusive
      if (verbose) write(6,*) 'inclusive',inclusive
      read(20,*) algorithm
      if (verbose) write(6,*) 'algorithm',algorithm
      read(20,*) ptjetmin
      if (verbose) write(6,*) 'ptjetmin',ptjetmin
      read(20,*) etajetmin 
      if (verbose) write(6,*) 'etajetmin',etajetmin
      read(20,*) etajetmax 
      if (verbose) write(6,*) 'etajetmax',etajetmax
      read(20,*) Rcut
      if (verbose) write(6,*) 'Rcut',Rcut
      read(20,*) makecuts
      if (verbose) write(6,*) 'makecuts',makecuts
      read(20,*) leptpt
      if (verbose) write(6,*) 'leptpt',leptpt
      read(20,*) leptrap
      if (verbose) write(6,*) 'leptrap',leptrap
      read(20,*) misspt
      if (verbose) write(6,*) 'misspt',misspt
      read(20,*) leptpt2
      if (verbose) write(6,*) 'leptpt2',leptpt2
      read(20,*) leptrap2
      if (verbose) write(6,*) 'leptrap2',leptrap2
      read(20,*) Rjlmin
      if (verbose) write(6,*) 'Rjlmin',Rjlmin
      read(20,*) Rllmin
      if (verbose) write(6,*) 'Rllmin',Rllmin
      read(20,*) delyjjmin
      if (verbose) write(6,*) 'delyjjmin',delyjjmin
      read(20,*) jetsopphem 
      if (verbose) write(6,*) 'jetsopphem',jetsopphem
      read(20,*) lbjscheme 
      if (verbose) write(6,*) 'lbjscheme',lbjscheme
      read(20,*) ptbjetmin
      if (verbose) write(6,*) 'ptbjetmin',ptbjetmin
      read(20,*) etabjetmax
      if (verbose) write(6,*) 'etabjetmax',etabjetmax
      read(20,*) gammpt
      if (verbose) write(6,*) 'gammpt',gammpt
      read(20,*) gammrap
      if (verbose) write(6,*) 'gammrap',gammrap
      read(20,*) gammcone
      if (verbose) write(6,*) 'gammcone',gammcone
      read(20,*) gammcut
      if (verbose) write(6,*) 'gammcut',gammcut

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- anomalous couplings 
      read(20,*) delg1_z
      if (verbose) write(6,*) 'delg1_z',delg1_z
      read(20,*) delk_z
      if (verbose) write(6,*) 'delk_z',delk_z
      read(20,*) delk_g
      if (verbose) write(6,*) 'delk_g',delk_g
      read(20,*) lambda_z
      if (verbose) write(6,*) 'lambda_z',lambda_z
      read(20,*) lambda_g
      if (verbose) write(6,*) 'lambda_g',lambda_g
      read(20,*) tevscale
      if (verbose) write(6,*) 'tevscale',tevscale

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) line
c--- grid information 
      read(20,*) readin
      if (verbose) write(6,*) 'readin',readin
      read(20,*) writeout
      if (verbose) write(6,*)'writeout',writeout
      read(20,*) ingridfile
      if (verbose) write(6,*) 'ingridfile',ingridfile
      read(20,*) outgridfile
      if (verbose) write(6,*) 'outgridfile',outgridfile

      if (verbose) write(6,*)

c--- check if the contents of technical.DAT are included here
c--- (the default behaviour going forward)
      technicalincluded=.false.

      read(20,99,end=88) line
      read(20,99,end=88) line
c      write(6,*) line
      if (line(2:10) .ne. 'Technical') goto 88
      technicalincluded=.true.

   88 continue
      
c--- if the contents of technical.DAT are not included, close the input
c--- file and open technical.DAT instead; otherwise continue on
      if (technicalincluded .eqv. .false.) then
        close(20)
        open(unit=20,file='technical.DAT',status='old',err=999)
        call checkversion(20,'technical.DAT')
      endif

      if (verbose) write(6,*) '[Technical parameters that'//
     .                        ' should not normally be changed]'
      if (verbose) write(6,*)

c---- read-in the technical parameters

      read(20,*) debug
      if (verbose) write(6,*) 'debug',debug
      read(20,*) verbose
      if (verbose) write(6,*) 'verbose',verbose
      read(20,*) new_pspace
      if (verbose) write(6,*) 'new_pspace',new_pspace
      read(20,*) virtonly
      if (verbose) write(6,*) 'virtonly',virtonly
      read(20,*) realonly
      if (verbose) write(6,*) 'realonly',realonly
      read(20,*) spira
      if (verbose) write(6,*) 'spira',spira
      read(20,*) noglue
      if (verbose) write(6,*) 'noglue',noglue
      read(20,*) ggonly
      if (verbose) write(6,*) 'ggonly',ggonly
      read(20,*) gqonly
      if (verbose) write(6,*) 'gqonly',gqonly
      read(20,*) vanillafiles
      if (verbose) write(6,*) 'vanillafiles',vanillafiles
      read(20,*) nmin
      if (verbose) write(6,*) 'nmin',nmin
      read(20,*) nmax
      if (verbose) write(6,*) 'nmax',nmax
      read(20,*) clustering
      if (verbose) write(6,*) 'clustering',clustering
      read(20,*) realwt
      if (verbose) write(6,*) 'realwt',realwt
      read(20,*) colourchoice
      if (verbose) write(6,*) 'colourchoice',colourchoice
      read(20,*) rtsmin
      if (verbose) write(6,*) 'rtsmin',rtsmin
      read(20,*) cutoff
      if (verbose) write(6,*) 'cutoff',cutoff
      read(20,*) aii
      if (verbose) write(6,*) 'aii',aii
      read(20,*) aif
      if (verbose) write(6,*) 'aii',aii
      read(20,*) afi
      if (verbose) write(6,*) 'aii',aii
      read(20,*) aff
      if (verbose) write(6,*) 'aii',aii
      if (verbose) write(6,*)
      close(unit=20)

      if ((etajetmin .lt. 0d0) .or. (etajetmax .lt. 0d0)) then
        write(6,*) 'etajetmin and etajetmax are absolute values,'
	write(6,*) ' please reset to a positive value.'
	stop
      endif

c--- reset values of the alpha parameters, for specific runstrings
      if (runstring(1:5) .eq. 'alpha') then
        if     (runstring(6:9) .eq. '1111') then
          aii=1d0
          aif=1d0
          afi=1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '0111') then
          aii=0.1d0
          aif=1d0
          afi=1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '1011') then
          aii=1d0
          aif=0.1d0
          afi=1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '1101') then
          aii=1d0
          aif=1d0
          afi=0.1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '1110') then
          aii=1d0
          aif=1d0
          afi=1d0
          aff=0.1d0
        elseif (runstring(6:9) .eq. '0001') then
          aii=0.1d0
          aif=0.1d0
          afi=0.1d0
          aff=1d0
        elseif (runstring(6:9) .eq. '0000') then
          aii=0.1d0
          aif=0.1d0
          afi=0.1d0
          aff=0.1d0
        else
          write(6,*) 'runstring=alpha.... reserved for internal tests'
          stop
        endif
        write(6,*) 'Values of alpha parameters reset, for testing:'
        write(6,*) '  alpha_ii=',aii
        write(6,*) '  alpha_if=',aif
        write(6,*) '  alpha_fi=',afi
        write(6,*) '  alpha_ff=',aff
        write(6,*)
      endif

      if     (index(runstring,'noglue') .gt. 0) then
        noglue=.true.
      elseif (index(runstring,'ggonly') .gt. 0) then
        ggonly=.true.
      elseif (index(runstring,'gqonly') .gt. 0) then
        gqonly=.true.
      endif
      
      if     (index(runstring,'mc1.3') .gt. 0) then
        mc=1.3d0
	mcsq=mc**2
      elseif (index(runstring,'mc1.4') .gt. 0) then
        mc=1.4d0
	mcsq=mc**2
      elseif (index(runstring,'mc1.5') .gt. 0) then
        mc=1.5d0
	mcsq=mc**2
      endif
      
      if (runstring(1:3) .eq. 'mlm') then
        write(6,*) 'WARNING: cross sections divided by Ecm**2'
	write(6,*)
      endif
      
      if (noglue) then
        write(6,*) 'WARNING: no gluon contribution included in PDF'
	write(6,*)
      endif
      if (ggonly) then
        write(6,*) 'WARNING: only gluon-gluon flux included'
	write(6,*)
      endif
      if (gqonly) then
        write(6,*) 'WARNING: only gluon-quark flux included'
	write(6,*)
      endif
      
c-----initialize various quantities

c--- save initial scale choices (that may be changed later)
      initscale=scale
      initfacscale=facscale

c--- assign squared masses for b- and c-quarks
      mbsq=mb**2
      mcsq=mc**2
      
c--- set-up mass window cuts
      bbsqmin=mbbmin**2
      bbsqmax=mbbmax**2

      wsqmin=Mwmin**2
      wsqmax=Mwmax**2

c--- set-up the variables for the process we wish to consider
      call chooser

c--- set-up the random number generator with a negative seed
      idum=-abs(origij)
      randummy=ran2()

c--- initialize masses for alpha_s routine
      cmass=dsqrt(mcsq)
      bmass=dsqrt(mbsq)

c--- E-M gauge invariance requires that delg1_g=0
      delg1_g=0d0

c--- check that we have a valid value of 'part'
      if ( (part .ne. 'lord') .and. (part .ne. 'real') .and.
     .     (part .ne. 'virt') .and. (part .ne. 'tota') ) then
        if ( (part .eq. 'todk') .and.
     .       ((case .eq. 'bq_tpq') .or. (case .eq. 't_bbar')
     .   .or. (case .eq. 'W_twdk')) ) then
c--- this is an allowed combination
        else 
          write(6,*) 'part=',part,' is not a valid option'
          write(6,*) 'for this process number.'
          stop     
        endif
      endif      

c--- set up the default choices of static scale, if required
      if (scale .lt. 0d0) then
	if     (scale .eq. -2d0) then
	  factor=0.25d0
	elseif (scale .eq. -3d0) then
	  factor=0.5d0
	elseif (scale .eq. -4d0) then
	  factor=0.75d0
	elseif (scale .eq. -5d0) then
	  factor=1d0
	elseif (scale .eq. -6d0) then
	  factor=2d0
	elseif (scale .eq. -7d0) then
	  factor=4d0
        else
	  factor=1d0
	endif	  
        if (n2+n3 .ne. 0) then
        scale=factor*(dfloat(n2)*mass2+dfloat(n3)*mass3)/dfloat(n2+n3)
c--- special cases where Higgs mass is neither mass2 nor mass3
        if (case(1:1) .eq. 'H') scale=factor*hmass	
        as=alphas(scale,amz,nlooprun)
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2
        write(6,*)
        write(6,*)'************* Strong coupling, alpha_s  ************'
        write(6,*)'*                                                  *'
        write(6,49)'alpha_s (scale)',gsq/fourpi
        write(6,49)'alpha_s (zmass)',amz
        write(6,50)' (using ',nlooprun,'-loop running of alpha_s)'  
        write(6,*)'****************************************************'
        write(6,*)
        write(6,*)'****************************************************'
        write(6,76) scale
        write(6,*)'****************************************************'
        else
        write(6,*) 'Invalid choice of renormalization scale!'
        stop
        endif
      endif
      if (facscale .lt. 0d0) then
	if     (facscale .eq. -2d0) then
	  factor=0.25d0
	elseif (facscale .eq. -3d0) then
	  factor=0.5d0
	elseif (facscale .eq. -4d0) then
	  factor=0.75d0
	elseif (facscale .eq. -5d0) then
	  factor=1d0
	elseif (facscale .eq. -6d0) then
	  factor=2d0
	elseif (facscale .eq. -7d0) then
	  factor=4d0
        else
	  factor=1d0
	endif	  
        if (n2+n3 .ne. 0) then
        facscale=factor*
     .           (dfloat(n2)*mass2+dfloat(n3)*mass3)/dfloat(n2+n3)
c--- special cases where Higgs mass is neither mass2 nor mass3
        if (case(1:1) .eq. 'H') facscale=factor*hmass	
        write(6,*)
        write(6,*)'****************************************************'
        write(6,77) facscale
        write(6,*)'****************************************************'
       else
        write(6,*) 'Invalid choice of factorization scale!'
        stop
        endif
      endif
      
      return

   49 format(' *  ',a20,f12.8,16x,'*')
   50 format(' *  ',6x,a8,i1,a25,8x,'*')
   76 format(' *      Renormalization scale =',f7.2,'              *')   
   77 format(' *        Factorization scale =',f7.2,'              *')   
   99 format(a90)

  999 continue
      write(6,*) 'Problem reading ',inputfile
      write(6,*)
      write(6,*) 'Refer to documentation for the format of input.DAT'
      write(6,*)
      stop

      end
      
