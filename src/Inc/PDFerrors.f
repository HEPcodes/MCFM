      logical PDFerrors
      integer maxPDFsets,currentPDF
c--- 50 is my choice for the maximum number of PDF error sets
c--- 40 is my choice for the maximum number of dipoles
      double precision PDFxsec(0:50),PDFxsec_nd(0:50,0:40),
     . PDFwgt(0:50)
      common/PDFerrors/PDFerrors,maxPDFsets,PDFxsec,PDFwgt  
