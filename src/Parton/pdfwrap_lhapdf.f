*****************
* LHAPDF version*
*****************
      subroutine pdfwrap
      implicit none
      include 'masses.f'
      include 'lhapdf.f'
      double precision amz,alphasPDF
      logical newinput
      common/newinput/newinput
      common/lhapdf/PDFmember,PDFname

      common/couple/amz

      
      if (newinput .eqv. .false.) then
        open(unit=21,file='lhapdf.DAT',status='old',err=999)
        call checkversion(21,'lhapdf.DAT')
        read(21,*) PDFname
        read(21,*) PDFmember            
        close(21)
      endif
      
      write(6,*)
      write(6,*) '**********************************'
      write(6,*) '*     MCFM is calling LHAPDF     *'
      write(6,*) '*                                *'
      write(6,98) 'PDFname',PDFname(1:20)
      write(6,99) 'PDFmember',PDFmember
      write(6,*) '**********************************'

      call InitPDFset('PDFsets/'//PDFname(1:20))
      call InitPDF(PDFmember)
      amz=alphasPDF(zmass)

      return
 
   98 format(' *   ',a7,' ',a20,' *')
   99 format(' *  ',a10,i3,'                 *')

  999 write(6,*) 'Error reading lhapdf.DAT'
      call flush(6)
      stop

      end
 

