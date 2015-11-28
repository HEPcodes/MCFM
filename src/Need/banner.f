      block data codeversion_data
      implicit none
      include 'codeversion.f'
      data codeversion/'4.1'/      
      end

      subroutine banner
************************************************************************
*  Set the version number of MCFM and write out the banner heading     *
************************************************************************
      implicit none
      include 'codeversion.f'
      character*50 line
      integer vlength,lenocc

      line='**************************************************'
      vlength=lenocc(codeversion)
      vlength=vlength+17   
      line(25-vlength/2:24+(vlength+1)/2)=
     . ' MCFM - version '//codeversion//' '

      write(6,*) line
      write(6,*) '*                                                *'
      write(6,*) '* MCFM, v'//codeversion//
     . '              January 17th, 2005  *'
      write(6,*) '*                                                *'
      write(6,*) '* Authors: John Campbell, John.Campbell@cern.ch  *'
      write(6,*) '*          Keith Ellis,   ellis@fnal.gov,        *'
      write(6,*) '*                                                *'
      write(6,*) '* For details see:                               *'
      write(6,*) '*w/Francesco Tramontano,hep-ph/0408158 (Sngl Top)*'
      write(6,*) '*  Phys.Rev.D65:113007, hep-ph/0202176 (W,Z+2j)  *'
      write(6,*) '*  Phys.Rev.D62:114012, hep-ph/0006304 (W,Z+bb)  *'
      write(6,*) '*  Phys.Rev.D60:113006, hep-ph/9905386 (diboson) *'
      write(6,*) '*                                                *'
      write(6,*) '* On the web:  http://mcfm.fnal.gov/             *'
      write(6,*) '**************************************************'
      write(6,*) 
 
      return
      end








