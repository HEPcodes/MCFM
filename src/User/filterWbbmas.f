      logical function filterWbbmas()
c--- this routine is specific to the "Wbbmas" processes 420-429;
c--- it inspects the jets to check whether an event should be included 
c--- routine returns FALSE if event does not pass the process-specific cuts
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'jetlabel.f'
      integer nproc
      common/nproc/nproc

      filterWbbmas=.true.

c---  20 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6) [massive]'
      if ((nproc .eq. 20) .or. (nproc .eq. 25)) then
        if (jets .lt. 2) then
	  filterWbbmas=.false.
	  return
c--- note: jets should already contain bq and ba because of bbproc
	endif
      endif

c--- 420 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5) [1 or 2 jets only]'
      if ((nproc .eq. 420) .or. (nproc .eq. 425)) then
        if ((jets .gt. 2) .or. (jets .lt. 1)) then
	  filterWbbmas=.false.
	  return
	endif
	if ((inclusive .eqv. .false.) .and. (jets .eq. 2)) then
	  filterWbbmas=.false.
	  return
	endif
	if (jets .eq. 1) then
	 if ((jetlabel(1) .ne. 'bq') .and. (jetlabel(1) .ne. 'ba')) then
	    filterWbbmas=.false.
	    return
	  endif
	endif
	if (jets .eq. 2) then
	  if ((jetlabel(1) .eq. 'bb') .or. (jetlabel(2) .eq. 'bb')) then
	    filterWbbmas=.false.
	    return
	  endif
	endif
      endif
      
c--- 421 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+(b+b~)(p5) [1 or 2 jets only]'
      if ((nproc .eq. 421) .or. (nproc .eq. 426)) then
        if ((jets .gt. 2) .or. (jets .lt. 1)) then
	  filterWbbmas=.false.
	  return
	endif
	if ((inclusive .eqv. .false.) .and. (jets .eq. 2)) then
	  filterWbbmas=.false.
	  return
	endif
	if (jets .eq. 1) then
	  if ((jetlabel(1) .ne. 'bb')) then
	    filterWbbmas=.false.
	    return
	  endif
	endif
	if (jets .eq. 2) then
         if ((jetlabel(1) .ne. 'bb') .and. (jetlabel(2) .ne. 'bb')) then
	    filterWbbmas=.false.
	    return
	 endif
	endif
      endif
      
c--- 422 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+jet(p6) [2 or 3 jets only]'
      if ((nproc .eq. 422) .or. (nproc .eq. 427)) then
        if (jets .lt. 2) then
	  filterWbbmas=.false.
	  return
	endif
	if ((inclusive .eqv. .false.) .and. (jets .eq. 3)) then
	  filterWbbmas=.false.
	  return
	endif
	if (jets .eq. 2) then
	  if ((jetlabel(1) .eq. 'bb') .or. (jetlabel(2) .eq. 'bb')) then
	    filterWbbmas=.false.
	    return
	  endif
	endif
      endif
      
      return
      end
      
