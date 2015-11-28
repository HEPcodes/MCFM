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
      character*72 outlabel1,outlabel2,outfiledat,outfiletop,outlabeltmp
      character*3 strmh,strscale,getstr
      double precision xsec,xsec_err
      common/part/part
      common/runstring/runstring

      strscale=getstr(int(scale))
      if (  (case .eq. 'WHbbar')
     . .or. (case .eq. 'ZHbbar')
     . .or. (case .eq. 'qq_tth')
     . .or. (case .eq. 'tottth')
     .) then
      strmh=getstr(int(hmass))
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//'_'//strmh
      else
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale
      endif
      nlength=lenocc(outlabel1)
      outlabel2=outlabel1(1:nlength)//'_'//runstring
      nlength=lenocc(outlabel2)

c--- add working directory, if necessary 
      if (workdir .ne. '') then
        outlabeltmp=outlabel2
        outlabel2=workdir(1:lenocc(workdir))//'/'//outlabeltmp
        nlength=nlength+1+lenocc(workdir)
      endif
      
      write(6,*)
      write(6,*) '****************************************************'
      write(6,*) 'output files  ',outlabel2
      write(6,*) '****************************************************'
      call flush(6)

      outfiledat=outlabel2
      outfiletop=outlabel2
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

      character*3 function getstr(no)
c returns a string of length 3 from an integer
      integer no,i1,i2,i3,zero
      
      zero=ichar('0')

      i1=no/100
      i2=(no-i1*100)/10
      i3=no-i1*100-i2*10

      if    (i1.eq.0.and.i2.eq.0) then
        getstr=char(i3+zero)//'__'
      elseif(i1.eq.0) then
        getstr=char(i2+zero)//char(i3+zero)//'_'
      else
        getstr=char(i1+zero)//char(i2+zero)//char(i3+zero)
      endif
      
      return
      end

