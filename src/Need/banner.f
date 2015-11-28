      block data codeversion_data
      implicit none
      include 'codeversion.f'
      data codeversion/'5.4'/      
      data      prelim/.false./      ! if true, print warning message
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

c--- warning message, if necessary
      if (prelim) then
        write(6,*) '*                                                *'
        write(6,*) '*           PRELIMINARY VERSION                  *'
        write(6,*) '*                                                *'
        write(6,*) '*  NOTE: This is a private release of the MCFM   *'
        write(6,*) '*  code that has not yet been made public on the *'
        write(6,*) '*  usual website. As such:                       *'
        write(6,*) '*                                                *'
        write(6,*) '*   + Please do not redistribute without the     *'
        write(6,*) '*     knowledge of the authors;                  *'
        write(6,*) '*                                                *'
        write(6,*) '*   + Please notify the authors of any bugs      *'
        write(6,*) '*     or problems so that they can be corrected  *'
        write(6,*) '*     before the next official release.          *'
        write(6,*) '*                                                *'
        write(6,*) line
      endif
     

      write(6,*) '*                                                *'
      write(6,*) '* MCFM, v'//codeversion//
     . '                March 12th, 2008  *'
      write(6,*) '*                                                *'
      write(6,*) '* Authors: John Campbell, Keith Ellis            *'
      write(6,*) '* (J.Campbell@physics.gla.ac.uk, ellis@fnal.gov) *'
      write(6,*) '*                                                *'
      write(6,*) '* For details see:                               *'
      write(6,*) '*  N.P.B 726:109(2005), hep-ph/0506289 (W+t)     *'
      write(6,*) '*  Phys.Rev.D70:094012, hep-ph/0408158 (Sngl Top)*'
      write(6,*) '*       (with Francesco Tramontano)              *'
      write(6,*) '*                                                *'
      write(6,*) '*  Phys.Rev.D65:113007, hep-ph/0202176 (W,Z+2j)  *'
      write(6,*) '*  Phys.Rev.D62:114012, hep-ph/0006304 (W,Z+bb)  *'
      write(6,*) '*  Phys.Rev.D60:113006, hep-ph/9905386 (diboson) *'
      write(6,*) '*                                                *'
      write(6,*) '* On the web:  http://mcfm.fnal.gov/             *'
      write(6,*) '**************************************************'
      write(6,*) 
 
      return
      end








