      subroutine masscuts(s,*)
      implicit none
      include 'constants.f'
      include 'npart.f'
      logical first
      double precision bbsqmin,bbsqmax,wsqmin,wsqmax,s(mxpart,mxpart)
      integer nproc
      common/nproc/nproc
      common/limits/bbsqmin,bbsqmax,wsqmin,wsqmax      
      data first/.true./

      if (first) then
      first=.false.
      write(6,*) 'sqrt(min34)',sqrt(wsqmin)
      write(6,*) 'sqrt(max34)',sqrt(wsqmax)
      write(6,*) 'sqrt(min56)',sqrt(bbsqmin)
      write(6,*) 'sqrt(max56)',sqrt(bbsqmax)
      endif
      
      if (  (s(3,4) .lt. wsqmin) 
     . .or. (s(3,4) .gt. wsqmax))
     .  return 1
      
      if (npart .gt. 3) then
        if (  (s(5,6) .lt. bbsqmin) 
     .   .or. (s(5,6) .gt. bbsqmax))
     .    return 1
      endif
      
      return
      end

