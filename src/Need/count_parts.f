      ! functions for counting particles 

      integer function count_photo()
      implicit none
      include 'constants.f'
      include 'plabel.f'
      integer j

      count_photo=0
      do j=1,mxpart
         if (plabel(j) .eq. 'ga') then 
            count_photo=count_photo+1
         endif
      enddo 

      return
      end

      
      integer function count_jets()
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'plabel.f'
      integer j

c---- Count final state jets only pp included initially
      count_jets=0

      do j=3,npart+2
         if (plabel(j) .eq. 'pp') then 
            count_jets=count_jets+1
         endif
      enddo 

      return 
      end
