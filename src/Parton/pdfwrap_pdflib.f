      subroutine pdfwrap(pdlabel)
      implicit none
      character*7 pdlabel
      character*20 parm(20)
      double precision value(20),amz,alphas2,zmass
      parameter(zmass=91.187d0)
      common/couple/amz
      parm(1)='ngroup'
      parm(2)='nset'
      parm(3)='LO'
      value(3)=2
      if (pdlabel .eq. 'dflm160') then
        value(1)=2
        value(2)=7 
      elseif (pdlabel .eq. 'dflm260') then
        value(1)=2
        value(2)=8 
      elseif (pdlabel .eq. 'dflm360') then
        value(1)=2
        value(2)=9 
      elseif (pdlabel .eq. 'hmrs90e') then
        value(1)=3
        value(2)=14
      elseif (pdlabel .eq. 'hmrs90b') then
        value(1)=3
        value(2)=17
      elseif (pdlabel .eq. 'k_mrs_b') then
        value(1)=3
        value(2)=21
      elseif (pdlabel .eq. 'mrs92s0') then
        value(1)=3
        value(2)=26
      elseif (pdlabel .eq. 'mrs92d0') then
        value(1)=3
        value(2)=27
      elseif (pdlabel .eq. 'mrs92d-') then
        value(1)=3
        value(2)=28
      elseif (pdlabel .eq. 'mrsb135') then
        value(1)=3
        value(2)=22
      elseif (pdlabel .eq. 'mrsb235') then
        value(1)=3
        value(2)=25
      elseif (pdlabel .eq. 'mrs96r1') then
        value(1)=3
        value(2)=53
      elseif (pdlabel .eq. 'mrs96r2') then
        value(1)=3
        value(2)=54
      elseif (pdlabel .eq. 'mrs96r3') then
        value(1)=3
        value(2)=55
      elseif (pdlabel .eq. 'mrs96r4') then
        value(1)=3
        value(2)=56
      elseif (pdlabel .eq. 'mtb2_91') then
        value(1)=4
        value(2)=3
      elseif (pdlabel .eq. 'cteq4_m') then
        value(1)=4
        value(2)=34
      elseif (pdlabel .eq. 'ehlqone') then
        value(1)=1
        value(2)=8
        value(3)=1
      elseif (pdlabel .eq. 'ehlqtwo') then
        value(1)=1
        value(2)=9
        value(3)=1
      else
        write(6,*) 'Unimplemented PDF'
        write(6,*) 'Implemented are: ',
     . 'dflm160,','dflm260,','dflm360,',
     . 'hmrs90e,','hmrs90b,','k_mrs_b,',
     . 'mrs92s0,','mrs92d0,','mrs92d-,',
     . 'mrsb135,','mrsb235,','mrs96r1,',
     . 'mrs96r2,','mrs96r3,', 'mrs96r4,',
     . 'cteq4_m,','mtb2_91,',
     . 'ehlqone,','ehlqtwo'
      stop
      endif
      call pdfset(parm,value)
      amz=alphas2(zmass)
      return
      end


