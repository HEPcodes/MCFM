      block data linlog_data
      implicit none
      include 'nplot.f'
      data linlog/71*'lin',2*'log',19*'lin',7*'log'/
      end

      subroutine histofin(xsec,xsec_err)
      implicit none
      include 'nplot.f'
      include 'masses.f'
      include 'scale.f'
      include 'verbose.f'
      include 'process.f'
      include 'workdir.f'
      include 'pdlabel.f'
      integer j,nlength,lenocc
      character*30 runstring
      character*4 part
      character*72 outlabel1,runname,outfiledat,outfiletop,outlabeltmp
      character*3 strmh,strscale,getstr
      double precision xsec,xsec_err
      common/part/part
      common/runstring/runstring
      common/runname/runname,nlength
      
      write(6,*)
      write(6,*) '****************************************************'
      write(6,*) 'output files  ',runname
      write(6,*) '****************************************************'
      call flush(6)

      outfiledat=runname
      outfiletop=runname
      outfiledat(nlength+1:nlength+4)='.dat'
      outfiletop(nlength+1:nlength+4)='.top'

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
      enddo
      close (unit=99)

      return
      end

