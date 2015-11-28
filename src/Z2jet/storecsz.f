      subroutine storecsz(mcs)
c-- this routine transfers the information on the colour structure
c-- for the Z2jet matrix elements into separate arrays for each
c-- incoming parton case
      implicit none
      include 'constants.f'
      include 'msq_cs.f'
      integer i,j,k
      double precision mcs(0:2,2,2)
      
      do i=0,2
        do j=1,2
          do k=1,2
        mcs(i,j,k)=mmsq_cs(i,j,k)
          enddo
        enddo
      enddo
c      write(6,*) 'storecsz:mcs',mcs      
c      write(6,*) 'storecsz:mmsq_cs',mmsq_cs
c      write(6,*) 'storecsz:msq_cs',msq_cs
      return
      end
