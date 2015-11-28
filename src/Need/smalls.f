      subroutine smalls(s,npart,*)
c    cut if radiated parton too close
      implicit none
      include 'constants.f'
      include 'cutoff.f'
      integer npart
      double precision s(mxpart,mxpart)

      if (npart .eq. 2) then
      if ( 
     .      (-s(1,4) .lt. cutoff)
     . .or. (-s(2,4) .lt. cutoff)
     . .or. (+s(3,4) .lt. cutoff)
     . .or. (-s(1,3) .lt. cutoff)
     . .or. (-s(2,3) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 3) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,4) .lt. cutoff)
     . .or. (-s(2,4) .lt. cutoff)
     . .or. (+s(4,5) .lt. cutoff)
     . .or. (+s(3,4) .lt. cutoff)
     . .or. (+s(3,5) .lt. cutoff)
     . .or. (-s(1,3) .lt. cutoff)
     . .or. (-s(2,3) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 4) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . ) return 1
     
      elseif (npart .eq. 5) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . .or. (+s(5,7) .lt. cutoff)
     . .or. (+s(6,7) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 6) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . .or. (-s(1,8) .lt. cutoff)
     . .or. (-s(2,8) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . .or. (+s(5,7) .lt. cutoff)
     . .or. (+s(5,8) .lt. cutoff)
     . .or. (+s(6,7) .lt. cutoff)
     . .or. (+s(6,8) .lt. cutoff)
     . .or. (+s(7,8) .lt. cutoff)
     . ) return 1
      
      endif      

      return
      end
