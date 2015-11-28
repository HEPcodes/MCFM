      subroutine checkorder(order)
c--- checks the value of nproc and part against a list of processes that
c--- are calculable only at LO. If the calculation is not possible,
c--- writes an error message and aborts
      implicit none
      include 'frag.f'
      integer nproc
      character*4 part
      character*1 order

      common/nproc/nproc
      common/part/part

c--- special cases where there is no LO calculation
      if ((part .ne. 'real') .and. (order .eq. 'R')) then 
        write(6,*)
        write(6,*)'This process can only be calculated with part=real,'
	write(6,*)'it is a subset of a NLO calculation only.'
	stop
      endif
      
c--- if we're calculating LO only, there's no problem      
      if (part .eq. 'lord') return
      
c--- otherwise, we must be performing a NLO calculation, and this list of
c--- process numbers can't be calculated beyond LO 
      if (order .eq. 'L') then
c        (nproc .eq.  14) .or. (nproc .eq.  19)
c     . .or. (nproc .eq.  23) .or. (nproc .eq.  28)
c     . .or. (nproc .eq.  24) .or. (nproc .eq.  29)
c     . .or. (nproc .eq.  34) 
c     . .or. (nproc .eq.  45) .or. (nproc .eq.  47)
c     . .or. (nproc .eq.  50)
c     . .or. (nproc .eq.  54) .or. (nproc .eq.  64).or. (nproc .eq.  66)
c     . .or. (nproc .eq.  85)
c     . .or. (nproc .eq. 121) .or. (nproc .eq. 122)
c     . .or. (nproc .eq. 156)
c     . .or. (nproc .eq. 160)
c     . .or. (nproc .eq. 190) .or. (nproc .eq. 191)
c     . .or. (nproc .eq. 196) .or. (nproc .eq. 197)
c     . .or. (nproc .eq. 201) .or. (nproc .eq. 202)
c     . .or. (nproc .eq. 206) .or. (nproc .eq. 207)
c     . .or. (nproc .eq. 216) .or. (nproc .eq. 217)
c     . .or. (nproc .eq. 221)
c     . .or. (nproc .eq. 232) .or. (nproc .eq. 237)
c     . .or. (nproc .eq. 263) .or. (nproc .eq. 264)
c     . .or. (nproc .eq. 275) .or. (nproc .eq. 276)
c     . .or. (nproc .eq. 281)
c     . .or. (nproc .eq. 292)
c     . .or. (nproc .eq. 308)
c     . .or. (nproc .eq. 311)
c     . .or. (nproc .eq. 316) .or. (nproc .eq. 321)
c     . .or. (nproc .eq. 326) .or. (nproc .eq. 331)
c     . .or. (nproc .eq. 336)
c     . .or. (nproc .eq. 346) .or. (nproc .eq. 347)
c     . .or. (nproc .eq. 356) .or. (nproc .eq. 357)
c     .   ) then
        write(6,*)
        write(6,*)'This process cannot be calculated beyond LO - please'
        write(6,*)'check the values of nproc and part then try again'
        stop
      endif

c--- check that fragmentation is not turned on for a non-fragmentation process
      if ((frag) .and. (order .ne. 'F')) then
        write(6,*)
	write(6,*) 'This process does not include photon fragmentation.'
	write(6,*) 'Please set frag=.false. in the input file.'
	stop
      endif
       
      return
      end
  
