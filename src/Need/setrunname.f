      subroutine setrunname(scalestart,fscalestart) 
      implicit none
      include 'flags.f'
      include 'masses.f'
      include 'process.f'
      include 'jetcuts.f'
      include 'workdir.f'
      include 'pdlabel.f'
      double precision scalestart,fscalestart
      integer nlength,lenocc
      character*30 runstring
      character*4 part
      character*72 outlabel1,runname,outlabeltmp
      character*3 strmh,strscale,getstr,strpt
      common/part/part
      common/runstring/runstring
      common/runname/runname
      common/nlength/nlength

      if (abs(scalestart-fscalestart) .lt. 1d0) then
c--- if the scales are the same, use this scale as the label
        strscale=getstr(int(scalestart))
       else
c--- .... otherwise, use the percentage of (muR/muF)
        strscale=getstr(int(scalestart/fscalestart*100d0))
      endif

      if (  (case .eq. 'WHbbar')
     . .or. (case .eq. 'ZHbbar')
     . .or. (case .eq. 'qq_tth')
     . .or. (case .eq. 'tottth')
     .) then
      strmh=getstr(int(hmass))
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//'_'//strmh
      elseif (  (case .eq. 'H_1jet') ) then
      strmh=getstr(int(hmass))
      strpt=getstr(int(ptjetmin))
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//
     . '_'//strmh//'_pt'//strpt(1:2)      
      elseif (  (case .eq. 'W_2jet')
     .     .or. (case .eq. 'Z_2jet') ) then
        if     (Gflag .eqv. .false.) then
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//'_qrk'
        elseif (Qflag .eqv. .false.) then
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale//'_glu'
        else
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale
        endif
      else
      outlabel1=case//'_'//part//'_'//pdlabel//'_'//strscale
      endif
      nlength=lenocc(outlabel1)
      runname=outlabel1(1:nlength)//'_'//runstring
      nlength=lenocc(runname)

c--- add working directory, if necessary 
      if (workdir .ne. '') then
        outlabeltmp=runname
        runname=workdir(1:lenocc(workdir))//'/'//outlabeltmp
        nlength=nlength+1+lenocc(workdir)
      endif
      
      return
      end


      character*3 function getstr(no)
c returns a string of length 3 from an integer
      integer no,i1,i2,i3,zero
      
      zero=ichar('0')

      i1=abs(no)/100
      i2=(abs(no)-i1*100)/10
      i3=abs(no)-i1*100-i2*10

      if    (i1.eq.0.and.i2.eq.0) then
        if (no .lt. 0) then
        getstr='-'//char(i3+zero)//'_'
        else
        getstr=char(i3+zero)//'__'
        endif
      elseif(i1.eq.0) then
        getstr=char(i2+zero)//char(i3+zero)//'_'
      else
        getstr=char(i1+zero)//char(i2+zero)//char(i3+zero)
      endif
      
      return
      end

