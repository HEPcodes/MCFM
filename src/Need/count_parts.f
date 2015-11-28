      ! functions for counting particles 

      integer function count_photo(p)
      implicit none
      include 'constants.f'
      character*2 plabel(mxpart)
      double precision p(mxpart,4) 
      integer j
      common/plabel/plabel

      count_photo=0
      do j=1,mxpart
         if (plabel(j) .eq. 'ga') then 
            count_photo=count_photo+1
         endif
      enddo 

      return
      end

      
      integer function count_jets(p)
      implicit none
      include 'constants.f'
      include 'npart.f'
      character*2 plabel(mxpart)
      double precision p(mxpart,4) 
      integer j
      common/plabel/plabel

c---- Count final state jets only pp included initially
      count_jets=0

      do j=3,npart+2
         if (plabel(j) .eq. 'pp') then 
            count_jets=count_jets+1
         endif
      enddo 

      return 
      end
