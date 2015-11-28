      subroutine checkorder()
c--- checks the value of nproc and part against a list of processes that
c--- are calculable only at LO. If the calculation is not possible,
c--- writes an error message and aborts
      integer nproc
      character*4 part

      common/nproc/nproc
      common/part/part

c--- special cases where there is no LO calculation
      if ((part .eq. 'lord') .and. 
     .    ( (nproc .eq. 342) .or. (nproc .eq. 352) ) ) then
        write(6,*)
        write(6,*)'This process cannot be calculated at LO, it is part'
	write(6,*)'of a NLO calculation only.'
	stop
      endif
      
c--- if we're calculating LO only, there's no problem      
      if (part .eq. 'lord') return
      
c--- otherwise, we must be performing a NLO calculation, and this list of
c--- process numbers can't be calculated beyond LO 
      if (  (nproc .eq.  14) .or. (nproc .eq.  19)
     . .or. (nproc .eq.  20) .or. (nproc .eq.  25)
     . .or. (nproc .eq.  23) .or. (nproc .eq.  28)
     . .or. (nproc .eq.  24) .or. (nproc .eq.  29)
     . .or. (nproc .eq.  45) .or. (nproc .eq.  50)
     . .or. (nproc .eq.  54) .or. (nproc .eq.  64)
     . .or. (nproc .eq.  85)
     . .or. (nproc .eq. 151)
     . .or. (nproc .eq. 152) .or. (nproc .eq. 156)
     . .or. (nproc .eq. 160)
     . .or. (nproc .eq. 190) .or. (nproc .eq. 191)
     . .or. (nproc .eq. 196) .or. (nproc .eq. 197)
     . .or. (nproc .eq. 201) .or. (nproc .eq. 202)
     . .or. (nproc .eq. 206) .or. (nproc .eq. 207)
     . .or. (nproc .eq. 216) .or. (nproc .eq. 217)
     . .or. (nproc .eq. 221) .or. (nproc .eq. 263)
     . .or. (nproc .eq. 264)
     . .or. (nproc .eq. 271) .or. (nproc .eq. 272)
     . .or. (nproc .eq. 273) .or. (nproc .eq. 274)
     . .or. (nproc .eq. 311)
     . .or. (nproc .eq. 316) .or. (nproc .eq. 321)
     . .or. (nproc .eq. 326) .or. (nproc .eq. 331)
     . .or. (nproc .eq. 336)
     . .or. (nproc .eq. 346) .or. (nproc .eq. 347)
     . .or. (nproc .eq. 356) .or. (nproc .eq. 357)
     .   ) then
        write(6,*)
        write(6,*)'This process cannot be calculated beyond LO - please'
        write(6,*)'check the values of nproc and part then try again'
        stop
      endif
      
      return
      end
  
