      subroutine histofin
      implicit none
      include 'nplot.f'
      include 'hmass.f'
      include 'scale.f'
      integer j,nlength,lenocc
      character*9 runstring
      character*6 case
      character*4 part
      common/part/part
      common/runstring/runstring
      common/process/case
      character*72 outlabel1,outlabel2
      character*3 strmh,strscale,getstr
      character outfile28*28,outfile32*32
      character pdlabel*7
      common/pdlabel/pdlabel      
      data linlog/71*'lin',2*'log',19*'lin',7*'log'/
c      character*(*) graph

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
      if (nlength .ge. 32) then
        nlength=32
      else
        nlength=28
      endif
      outlabel2=outlabel2(1:nlength)
      write(6,*) 'output files  ',outlabel2
      if     (nlength .eq. 28) then
        outfile28=outlabel2        
        open(unit=98,file=outfile28//'.dat',status='unknown')
        open(unit=99,file=outfile28//'.top',status='unknown')
      elseif (nlength .eq. 32) then
        outfile32=outlabel2        
        open(unit=98,file=outfile32//'.dat',status='unknown')
        open(unit=99,file=outfile32//'.top',status='unknown')
      endif
      
      do j=1,nplot
      call mfinal(j)
      enddo

      do j=1,nplot
      call mprint(j)
      enddo
      close (unit=98)

c---generate topdrawer file
      do j=1,nplot
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

