c--- MBOOK-style routines that are supposed to create corresponding
c--- histograms for each of the error PDF's, when requested in nplotter

c--- These are the MBOOK common blocks
c      DOUBLE PRECISION HIST(100,100),XHIS(100,100),HDEL(100),HMIN(100)
c     &,HMAX(100),HAVG(100),HINT(100),HSIG(100)
c      COMMON/HISTOR/HIST,XHIS,HDEL,HMIN,HMAX,HAVG,HINT,HSIG
c      CHARACTER TITLE*100,BOOK*3
c      COMMON/HISTOC/BOOK(100),TITLE(100)                         
c      INTEGER NBIN(100),IHIS(100,100),IUSCORE(100),IOSCORE(100),
c     & IENT(100),NHIST
c      COMMON/HISTOI/NBIN,IHIS,IUSCORE,IOSCORE,IENT,NHIST


c--- This is the EBOOK common block - note that most entries are not
c--- present here, to save on storage space. The maximum number of
c--- histograms that may be calculated with errors is 4 and the
c--- maximum number of PDF error sets is 40
c      COMMON/EHISTO/HIST(4,40,100),HISTOMATCH(100),COUNTHISTO

      block data einitialize
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO 
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO                    
      data IHISTOMATCH/100*0/,ICOUNTHISTO/0/
      end

c--- sets up the histogram
      subroutine ebook(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO                    
c--- This is the MBOOK common block
      include 'histo.f'
c--- This is the EBOOK common block - note that most entries are not
c--- present here, to save on storage space. The maximum number of
c--- histograms that may be calculated with errors is 4 and the
c--- maximum number of PDF error sets is 40
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      
      ICOUNTHISTO=ICOUNTHISTO+1
      
      if (ICOUNTHISTO .GT. 4) then
        write(6,*) 'Only 4 histograms with errors allowed!'
        stop
      endif
      
      IHISTOMATCH(N)=ICOUNTHISTO
      
      DO I=1,NBIN(N)
        DO J=1,40
      EHIST(ICOUNTHISTO,J,I)=0d0
        ENDDO
      ENDDO
      
      return
      end
      
      SUBROUTINE efill(N,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'PDFerrors.f'
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO                    
c--- This is the MBOOK common block
      include 'histo.f'
c--- This is the EBOOK common block - note that most entries are not
c--- present here, to save on storage space. The maximum number of
c--- histograms that may be calculated with errors is 4 and the
c--- maximum number of PDF error sets is 40
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      I=INT((X-HMIN(N))/HDEL(N)+1D0)
      IF(I.GT.0.AND.I.LE.NBIN(N))  THEN
      NMATCH=IHISTOMATCH(N)
      
      DO J=1,maxPDFsets
        EHIST(NMATCH,J,I)=EHIST(NMATCH,J,I)
     .   +PDFwgt(J)/hdel(n)
      ENDDO
c     we are renormalising the weights by the bin width
      ENDIF
      
      return
      end
      

      SUBROUTINE etop(N,M,BTIT,LTIT,SCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'PDFerrors.f'
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO                    
      CHARACTER*(*) LTIT,BTIT,SCALE
c--- This is the MBOOK common block
      include 'histo.f'
c--- This is the EBOOK common block - note that most entries are not
c--- present here, to save on storage space. The maximum number of
c--- histograms that may be calculated with errors is 4 and the
c--- maximum number of PDF error sets is 40
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO

      IF(BOOK(N).NE.'YES') RETURN

      NMATCH=IHISTOMATCH(N)
      
c--- loop over all PDF sets
      DO K=1,maxPDFsets
      
      WRITE(97,100) TITLE(N),BTIT,LTIT,SCALE,HMIN(N),HMAX(N)
  100 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A50,'"',/1X,
     &' TITLE BOTTOM ','"',A50,'"',/1X,
     &' TITLE LEFT ','"',A50,'"',/1X,
     &' SET SCALE Y ',A5,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ')
      DO 1 J=1,NBIN(N)
      IF(EHIST(NMATCH,K,J).EQ.0.) GO TO 1
      WRITE(97,'(3X,G13.6,2(2X,G13.6))')  
     &                  XHIS(N,J),EHIST(NMATCH,K,J),0d0
    1 CONTINUE
      WRITE(97,200)
  200 FORMAT('   PLOT')
      WRITE(97,400)
  400 FORMAT('   NEW PLOT')
  
      enddo
      
      return
      end
      
      SUBROUTINE emtop(N,M,BTIT,LTIT,SCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'PDFerrors.f'
      double precision EHIST(4,40,100),maxhist(100),minhist(100)
      integer IHISTOMATCH(100),ICOUNTHISTO                    
      CHARACTER*(*) LTIT,BTIT,SCALE
c--- This is the MBOOK common block
      include 'histo.f'
c--- This is the EBOOK common block - note that most entries are not
c--- present here, to save on storage space. The maximum number of
c--- histograms that may be calculated with errors is 4 and the
c--- maximum number of PDF error sets is 40
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO

      IF(BOOK(N).NE.'YES') RETURN

c--- this bit writes out the normal plot
      WRITE(99,100) 'Errors: '//TITLE(N)(1:92),
     . BTIT,LTIT,SCALE,HMIN(N),HMAX(N)
  100 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A50,'"',/1X,
     &' TITLE BOTTOM ','"',A50,'"',/1X,
     &' TITLE LEFT ','"',A50,'"',/1X,
     &' SET SCALE Y ',A5,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ')
      WRITE(99,200)
      write(99,*) 'SET COLOR RED'
      DO J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 99
      WRITE(99,'(3X,G13.6,2(2X,G13.6))')  
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
   99 continue
      ENDDO
      
      WRITE(99,200)
  200 FORMAT('   PLOT')

c--- this bit writes out the error bounds      
      NMATCH=IHISTOMATCH(N)
      
c--- loop over all PDF sets
      DO J=1,NBIN(N)
        maxhist(J)=HIST(N,J)
        minhist(J)=HIST(N,J)
        DO K=1,maxPDFsets      
        IF(EHIST(NMATCH,K,J).ne.0.) then
          if (EHIST(NMATCH,K,J) .gt. maxhist(j))
     .      maxhist(j)=EHIST(NMATCH,K,J)
          if (EHIST(NMATCH,K,J) .lt. minhist(j))
     .      minhist(j)=EHIST(NMATCH,K,J)
        ENDIF
        ENDDO
      ENDDO

  201 FORMAT('  SET ORDER X Y DY')
  
      write(99,201)
      write(99,*) 'SET COLOR BLUE'
      DO J=1,NBIN(N)
      if (minhist(J) .ne. 0.) then
      WRITE(99,'(3X,G13.6,2(2X,G13.6))')  
     &                  XHIS(N,J),minhist(J),0d0
      endif
      ENDDO
      write(99,200)
    
      write(99,201)
      DO J=1,NBIN(N)
      if (maxhist(J) .ne. 0.) then
      WRITE(99,'(3X,G13.6,2(2X,G13.6))')  
     &                  XHIS(N,J),maxhist(J),0d0
      endif
      ENDDO
      write(99,200)
      write(99,*) 'SET COLOR WHITE'

c--- plot statistics, as normal           
      WRITE(99,300) HINT(N),HAVG(N),HSIG(N),IENT(N),IUSCORE(N)
     &   ,IOSCORE(N)
  300 FORMAT( /1x,                               
     &' BOX 7. 0.75 SIZE 9. 1.5'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' TITLE 2.8 1.2 "INTGRL =',E12.5,'   AVGE =',E12.5,
     &             '   RMS =',E12.5,'"',/1X,
     &' TITLE 2.8 0.8 "Entries =',I9,2x,'U`flow =',I9,2X
     &                                 ,'O`flow =',I9,'"',/1X,
     &' SET TITLE SIZE -2')

      WRITE(99,400)
      
  400 FORMAT('   NEW PLOT')

      return
      end
            
      
