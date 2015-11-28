*****************
* LHAPDF version*
*****************
      subroutine pdfwrap
      implicit none
      include 'masses.f'
      include 'lhapdf.f'
      double precision amz,alphasPDF
      logical newinput,validPDF
      character*30 oldPDFname
      integer i
      common/newinput/newinput

      common/couple/amz

      
      if (newinput .eqv. .false.) then
        open(unit=21,file='lhapdf.DAT',status='old',err=999)
        call checkversion(21,'lhapdf.DAT')
        read(21,*) PDFname
        read(21,*) PDFmember            
        close(21)
      endif
      
      oldPDFname=PDFname
      validPDF=.false.
      i=0
   20 continue
      i=i+1    
      if ((oldPDFname(i:i) .eq. '.') .or.
     .    (oldPDFname(i:i) .eq. ' ') .or.
     .    (oldPDFname(i:i) .eq. '[')) then
        validPDF=.true.
        PDFname=oldPDFname(1:i-1)//'.LHpdf'
      endif  
      if ((i .lt. 20) .and. (validPDF .eqv. .false.)) goto 20
      
      if (validPDF .eqv. .false.) then
        write(6,*) 'Problem with PDFname'
        write(6,*)
        stop
      endif
      
      write(6,*)
      write(6,*) '**********************************'
      write(6,*) '*     MCFM is calling LHAPDF     *'
      write(6,*) '*                                *'
      write(6,98) 'PDFname',PDFname(1:20)
      write(6,99) 'PDFmember',PDFmember
      write(6,*) '**********************************'

      call InitPDFset('PDFsets/'//PDFname)
      call InitPDF(PDFmember)
      amz=alphasPDF(zmass)

      return
 
   98 format(' *   ',a7,' ',a20,' *')
   99 format(' *  ',a10,i3,'                 *')

  999 write(6,*) 'Error reading lhapdf.DAT'
      call flush(6)
      stop

      end
 

