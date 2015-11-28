      subroutine zerorealhistos
c--- zero out all entries in the temporary histograms used for
c--- binning the weights in the real contribution
      implicit none
      include 'nplot.f'
      integer nplotmax,j
      common/nplotmax/nplotmax
      
      do j=1,nplotmax
      call mzero(3*maxhisto+j)
      enddo
      
      return
      end
      

      subroutine addrealhistos(wgt)
c--- add temporay histograms to the cumulative ones
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      double precision wgt
      integer nplotmax,j
      common/nplotmax/nplotmax

c--- loop over plots      
      do I=1,nplotmax

      DO L=1,NBIN(I)

c--- add weights      
      HIST(I,L)=HIST(I,L) + HIST(3*maxhisto+I,L)
c--- add errors
      if (abs(wgt) .gt. 1d-15) then         ! for safety
      HIST(maxhisto+I,L)=HIST(maxhisto+I,L)
     . + HIST(3*maxhisto+I,L)**2/wgt*HDEL(I)
c--- we want (f**2*wgt), so we do (f*wgt)**2/wgt ; the extra
c--- factor of HDEL(I) is to account for the normalization by the
c--- bin width (c.f. MFILL in mbook.f)
      endif
      
c--- count entries
      if (IHIS(3*maxhisto+I,L) .GT. 0) then
        IHIS(I,L)=IHIS(I,L)+1
        IENT(I)=IENT(I)+1
        IHIS(maxhisto+I,L)=IHIS(maxhisto+I,L)+1
        IENT(maxhisto+I)=IENT(maxhisto+I)+1
      endif

      ENDDO
      
      enddo
      
      return
      end
      
