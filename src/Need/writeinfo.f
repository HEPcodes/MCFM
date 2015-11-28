      subroutine writeinfo(unitno,xsec,xsec_err)
************************************************************************
*   Routine to write out run information to a desired unit             *
************************************************************************
      implicit none
      include 'maxwt.f'
      include 'masses.f'
      include 'facscale.f'
      include 'scale.f'
      include 'zerowidth.f'
      include 'flags.f'
      include 'clustering.f'
      include 'anomcoup.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'jetcuts.f'
      include 'lhapdf.f'
      include 'pdlabel.f'
      integer unitno
      double precision xsec,xsec_err
      
      character*4 part
      character*30 runstring
      logical creatent,dswhisto,dryrun,makecuts
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,origij
      integer NPTYPE,NGROUP,NSET
      double precision sqrts
      double precision Rcut
      double precision leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut
      integer lbjscheme
      logical jetsopphem
 
      common/outputflags/creatent,dswhisto      

      common/nproc/nproc
      common/part/part
      common/runstring/runstring
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts
      common/leptcuts/leptpt,leptrap,misspt,Rjlmin,Rllmin,delyjjmin,
     . leptpt2,leptrap2,gammpt,gammrap,gammcone,gammcut,
     . lbjscheme,jetsopphem

      common/origij/origij

      write(unitno,*) '( Cross-section is: ',xsec,'+/-',xsec_err,')'
      write(unitno,*)
      write(unitno,*) '( Run corresponds to this input file)'
      
      write(unitno,*)
      write(unitno,*)
     . '( [Flags to specify the mode in which MCFM is run] )'
      write(unitno,98) evtgen,'evtgen'
      write(unitno,98) creatent,'creatent'
      write(unitno,98) skipnt,'skipnt'
      write(unitno,98) dswhisto,'dswhisto'

      write(unitno,*)
      write(unitno,*)
     . '( [General options to specify the process and execution] )'
      write(unitno,97) nproc,'nproc'
      write(unitno,96) part,'part'
      write(unitno,96) runstring,'runstring'
      write(unitno,99) sqrts,'sqrts'
      write(unitno,97) ih1,'ih1'
      write(unitno,97) ih2,'ih2'
      write(unitno,99) hmass,'hmass'
      write(unitno,99) scale,'scale'
      write(unitno,99) facscale,'facscale'
      write(unitno,98) zerowidth,'zerowidth'
      write(unitno,97) itmx1,'itmx1'
      write(unitno,97) ncall1,'ncall1'
      write(unitno,97) itmx2,'itmx2'
      write(unitno,97) ncall2,'ncall2'
      write(unitno,97) origij,'ij'
      write(unitno,98) dryrun,'dryrun'
      write(unitno,98) Qflag,'Qflag'
      write(unitno,98) Gflag,'Gflag'
      
      write(unitno,*)
      write(unitno,*) 
     . '( [Pdf selection] )'
      write(unitno,96) pdlabel,'pdlabel '
      write(unitno,97) NGROUP,'NGROUP'
      write(unitno,97) NSET,'NSET'
      write(unitno,96) PDFname,'LHAPDF group'
      write(unitno,97) PDFmember,'LHAPDF set'

      write(unitno,*)
      write(unitno,*)
     . '( [Jet definition and event cuts] )'
      write(unitno,99) dsqrt(wsqmin),'m34min'
      write(unitno,99) dsqrt(wsqmax),'m34max'
      write(unitno,99) dsqrt(bbsqmin),'m56min'
      write(unitno,99) dsqrt(bbsqmax),'m56max'
      write(unitno,98) inclusive,'inclusive'
      write(unitno,96) algorithm,'algorithm'
      write(unitno,99) ptjetmin,'ptjetmin'
      write(unitno,99) etajetmax,'etajetmax'
      write(unitno,99) Rcut,'Rcut'
      write(unitno,98) makecuts,'makecuts'
      write(unitno,99) leptpt,'leptpt'
      write(unitno,99) leptrap,'leptrap'
      write(unitno,99) misspt,'misspt'
      write(unitno,99) leptpt2,'leptpt2'
      write(unitno,99) leptrap2,'leptrap2'
      write(unitno,99) Rjlmin,'Rjlmin'
      write(unitno,99) Rllmin,'Rllmin'
      write(unitno,99) delyjjmin,'delyjjmin'
      write(unitno,98) jetsopphem,'jetsopphem'
      write(unitno,97) lbjscheme,'lbjscheme'
      write(unitno,99) gammpt,'gammpt'
      write(unitno,99) gammrap,'gammrap'
      write(unitno,99) gammcone,'gammcone'
      write(unitno,99) gammcut,'gammcut'

      write(unitno,*)
      write(unitno,*)
     . '( [Anomalous couplings of the W and Z] )'
      write(unitno,99) delg1_z,'delg1_z'
      write(unitno,99) delk_z,'delk_z'
      write(unitno,99) delk_g,'delk_g'
      write(unitno,99) lambda_z,'lambda_z'
      write(unitno,99) lambda_g,'lambda_g'
      write(unitno,99) tevscale,'tevscale'

      write(unitno,*)
      write(unitno,*) 
     . '( [How to resume/save a run] )'
      write(unitno,98) readin,'readin'
      write(unitno,98) writeout,'writeout'
      write(unitno,96) ingridfile,'ingridfile'
      write(unitno,96) outgridfile,'outgridfile'

      write(unitno,99)

      return

c--- 96 character format      
   96 format(' (',a20,12x,'[',a,']',' )')  
c--- 97 integer format      
   97 format(' (',i20,12x,'[',a,']',' )')  
c--- 98 logical format      
   98 format(' (',L20,12x,'[',a,']',' )')  
c--- 99 floating point format
   99 format(' (',f20.4,12x,'[',a,']',' )')  
      
      end
      
      
