      subroutine zerorealhistos
c--- zero out all entries in the temporary histograms used for
c--- binning the squared weights in the real contribution
      implicit none
      include 'nplot.f'
      integer nplotmax,j
      common/nplotmax/nplotmax
      
      do j=1,nplotmax
      call mzero(3*maxhisto+j)
      enddo
      
      return
      end
      

      subroutine addrealhistos
c--- add temporay histograms to the cumulative ones
      implicit none
      include 'nplot.f'
      integer nplotmax,j
      common/nplotmax/nplotmax
      
      do j=1,nplotmax
      call mopera(maxhisto+j,'+',3*maxhisto+j,maxhisto+j,1d0,1d0)
      enddo
      
      return
      end
      
