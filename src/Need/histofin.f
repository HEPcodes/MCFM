      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/71*'lin',2*'log',19*'lin',7*'log'/
      end

      subroutine histofin(xsec,xsec_err)
      implicit none
      include 'nplot.f'
      include 'verbose.f'
      include 'PDFerrors.f'
      integer j,nlength
      character*72 runname,outfiledat,outfiletop,outfileerr
      double precision xsec,xsec_err
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO                    
      common/runname/runname
      common/nlength/nlength
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      
      write(6,*)
      write(6,*) '****************************************************'
      write(6,*) 'output files  ',runname
      write(6,*) '****************************************************'
      call flush(6)

      outfiledat=runname
      outfiletop=runname
      outfileerr=runname
      outfiledat(nlength+1:nlength+4)='.dat'
      outfiletop(nlength+1:nlength+4)='.top'
      outfileerr(nlength+1:nlength+10)='_error.top'

      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        open(unit=97,file=outfileerr,status='unknown')
      endif
      open(unit=98,file=outfiledat,status='unknown')
      open(unit=99,file=outfiletop,status='unknown')
      
c--- write out run info to top of files
      call writeinfo(98,xsec,xsec_err)      
      call writeinfo(99,xsec,xsec_err)      
      
      do j=1,nplot
      if (verbose) then
c        write(6,*) 'Finalizing plot ',j
        call flush(6)
      endif
      call mfinal(j)
      enddo

      do j=1,nplot
      if (verbose) then
c        write(6,*) 'Writing .dat for plot ',j
        call flush(6)
      endif
      call mprint(j)
      enddo
      close (unit=98)

c---generate topdrawer file
      do j=1,nplot
      if (verbose) then
c        write(6,*) 'Writing .top for plot ',j
        call flush(6)
      endif
      call mtop(j,100,'x','y',linlog(j))
      if ((PDFerrors) .and. (IHISTOMATCH(j) .ne. 0)) then
        call emtop(j,100,'x','y',linlog(j))
      endif
      enddo
      close (unit=99)

c---generate error file
      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        do j=1,nplot
          if (IHISTOMATCH(j) .ne. 0) then
            if (verbose) then
c              write(6,*) 'Writing .top for plot ',j
              call flush(6)
            endif
            call etop(j,100,'x','y',linlog(j))
          endif
        enddo
        close (unit=97)
      endif
      
      return
      end

