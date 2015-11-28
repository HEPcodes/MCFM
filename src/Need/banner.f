      subroutine banner
************************************************************************
*  Set the version number of MCFM and write out the banner heading     *
************************************************************************
      implicit none
      character*6 codeversion
      character*50 line
      integer vlength,lenocc
      common/versionnumber/codeversion
      data codeversion/'3.2'/      

      line='**************************************************'
      vlength=lenocc(codeversion)
      vlength=vlength+17   
      write(6,*) vlength
      line(25-vlength/2:24+(vlength+1)/2)=
     . ' MCFM - version '//codeversion//' '

      write(6,*) line
      write(6,*) '*                                                *'
      write(6,*) '* MCFM, v'//codeversion//
     . '              September 3rd, 2002 *'
      write(6,*) '*                                                *'
      write(6,*) '* Authors: John Campbell, johnmc@hep.anl.gov     *'
      write(6,*) '*          Keith Ellis,   ellis@fnal.gov,        *'
      write(6,*) '*                                                *'
      write(6,*) '* For details see:                               *'
      write(6,*) '*  Phys.Rev.D65:113007, hep-ph/0202176 (W,Z+2j)  *'
      write(6,*) '*  Phys.Rev.D62:114012, hep-ph/0006304 (W,Z+bb)  *'
      write(6,*) '*  Phys.Rev.D60:113006, hep-ph/9905386 (diboson) *'
      write(6,*) '*                                                *'
      write(6,*) '* On the web:  http://theory.fnal.gov/           *'
      write(6,*) '*               people/ellis/Programs/mcfm.html  *'
      write(6,*) '**************************************************'
      write(6,*) 
 
      return
      end








