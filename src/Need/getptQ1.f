      subroutine getptQ1(pt5,pt6,ptQ1)
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      double precision pt5,pt6,ptQ1,ptQ2

c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      if (jets .eq. 1) then
        if ((jetlabel(1) .eq. 'bq') .or. (jetlabel(1) .eq. 'ba')) then
          ptQ1=pt5
          return
        else
          write(6,*) 'Error in getptQ1: only 1 jet and it'
          write(6,*) ' is not a heavy quark!'
          stop
        endif
      endif

      if (jets .ne. 2) then
        write(6,*) 'Error in getptQ1: strange number of jets, ',jets,'!'
        stop
      endif

c--- now we know that we have 2 jets      
      ptQ1=-1d0
      ptQ2=-1d0

      if ((jetlabel(1) .eq. 'bq') .or. (jetlabel(1) .eq. 'ba')) ptQ1=pt5
      if ((jetlabel(2) .eq. 'bq') .or. (jetlabel(2) .eq. 'ba')) ptQ2=pt6
      
      ptQ1=max(ptQ1,ptQ2)
      
      if (ptQ1 .lt. 0d0) then
        write(6,*) 'Error in getptQ1: 2 jets, but no heavy quarks!'
        stop
      endif

      return
      end
