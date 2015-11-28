      subroutine swapjet(pjet,jetlabel,jetindex,i,j)
c--- swaps jets i..j in pjet
      implicit none
      include 'constants.f'
      integer i,j,k,itmp,jetindex(mxpart)
      double precision pjet(mxpart,4),tmp
      character jetlabel(mxpart)*2,chartmp*2
 
      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
 
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp
      
c      if ((jetlabel(i) .ne. 'pp') .and. (jetlabel(j) .ne. 'pp')) then
c      itmp=jetindex(i)
c      jetindex(i)=jetindex(j)
c      jetindex(j)=itmp
c      endif
      
      return
      end
