      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/nplot*'lin'/
      end

      subroutine histofin(xsec,xsec_err,itno,itmx)
c--- This outputs the final histograms for itno=0
c--- For itno>0, this is an intermediate result only
      implicit none
      include 'verbose.f'
      include 'PDFerrors.f'
      include 'histo.f'
      include 'vanillafiles.f'
      integer j,nlength,itno,itmx,nplotmax,nempty
      character*255 runname,outfiledat,outfiletop,outfileerr
      character*3 oldbook
      character mop
      double precision xsec,xsec_err,scalefac,itscale
      double precision EHIST(4,40,100)   
      integer IHISTOMATCH(100),ICOUNTHISTO  
      logical scaleplots                  
      common/runname/runname
      common/nlength/nlength
      common/nplotmax/nplotmax
      COMMON/EHISTO/EHIST,IHISTOMATCH,ICOUNTHISTO
      common/scaleplots/scalefac,scaleplots
      
      
      if (itno .eq. 0) then
      write(6,*)
      write(6,*) '****************************************************'
      if (vanillafiles) then
        write(6,*) 'output file name normally  ',runname(1:nlength)
        runname='mcfm-output'
	nlength=11
	write(6,*)
        write(6,*) '  but renamed to: ',runname
      else
        write(6,*) 'output files  ',runname(1:nlength)
      endif
      write(6,*) '****************************************************'
      call flush(6)
      scaleplots=.false.
      else
      scaleplots=.true.
      scalefac=1d0/dfloat(itno)
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

c--- make sure to scale results by the maximum number of iterations
      if (itno .eq. 0) then
        itscale=1d0/dfloat(itmx)
        mop='V'
      else
        itscale=1d0/dfloat(itno)
        mop='U'
      endif
      
c--- calculate the errors in each plot (and store in 2*maxhisto+j)    
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Calculating errors for plot ',j
        call flush(6)
      endif
      call mopera(j,mop,maxhisto+j,2*maxhisto+j,itscale,1d0)
      enddo

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

      nempty=0
      do j=1,nplotmax
      if (verbose) then
c        write(6,*) 'Writing .dat for plot ',j
        call flush(6)
      endif
      call mprint(j)
      if (book(j) .ne. 'YES') nempty=nempty+1
      enddo
      close (unit=98)

c---generate topdrawer file - only for non-empty plots
      do j=1,nplotmax-nempty
      if (verbose) then
c        write(6,*) 'Writing .top for plot ',j
        call flush(6)
      endif
      call mtop(j,2*maxhisto+j,'x','y',linlog(j))
      if ((PDFerrors) .and. (IHISTOMATCH(j) .ne. 0)) then
        call emtop(j,2*maxhisto,'x','y',linlog(j))
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
            call etop(j,2*maxhisto,'x','y',linlog(j))
          endif
        enddo
        close (unit=97)
      endif
      
      return
      end

