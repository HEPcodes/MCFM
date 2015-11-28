      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/150*'lin'/
      end

      subroutine histofin(xsec,xsec_err,itno)
c--- This outputs the final histograms for itno=0
c--- For itno>0, this is an intermediate result only
      implicit none
      include 'nplot.f'
      include 'verbose.f'
      include 'PDFerrors.f'
      include 'histo.f'
      integer j,nlength,itno,nplotmax
      character*72 runname,outfiledat,outfiletop,outfileerr
      character*4 mypart
      character*3 oldbook
      double precision xsec,xsec_err
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO                    
      common/runname/runname
      common/nlength/nlength
      common/nplotmax/nplotmax
      common/mypart/mypart
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      
      if (itno .eq. 0) then
      write(6,*)
      write(6,*) '****************************************************'
      write(6,*) 'output files  ',runname
      write(6,*) '****************************************************'
      call flush(6)
      endif

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
      call writeinfo(98,xsec,xsec_err,itno)      
      call writeinfo(99,xsec,xsec_err,itno)      

c--- calculate the errors in each plot (and store in 100+j)  - lord and virt only   
      if ((mypart .eq. 'lord') .or. (mypart .eq. 'virt')) then
        do j=1,nplotmax
        if (verbose) then
c          write(6,*) 'Finalizing plot ',j
          call flush(6)
        endif
        call mopera(j,'V',50+j,100+j,1d0,1d0)
        enddo
      endif
      
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Finalizing plot ',j
        call flush(6)
      endif
c--- ensure that MFINAL doesn't turn off booking for intermediate results
      oldbook=book(j)
      call mfinal(j)
      if (itno .gt. 0) then
      book(j)=oldbook
      endif
      enddo

      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Writing .dat for plot ',j
        call flush(6)
      endif
      call mprint(j)
      enddo
      close (unit=98)

c---generate topdrawer file
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Writing .top for plot ',j
        call flush(6)
      endif
      call mtop(j,100+j,'x','y',linlog(j))
      if ((PDFerrors) .and. (IHISTOMATCH(j) .ne. 0)) then
        call emtop(j,100,'x','y',linlog(j))
      endif
      enddo
      close (unit=99)

c---generate error file
      if ((PDFerrors) .and. (ICOUNTHISTO .gt. 0)) then
        do j=1,nplotmax
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

