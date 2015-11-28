      subroutine fdist(pdlabel,ih,x,scale,xf)
      implicit none
      character*7 pdlabel,pdlbold
      double precision x,scale,xf(-5:5),dxpdf(-6:6)
      integer j,ih
      data pdlbold/'startup'/
     
      if (pdlabel .ne. pdlbold) then
         call pdfwrap(pdlabel)
         pdlbold=pdlabel
      endif
      call pftopdg(x,scale,dxpdf)
      do j=-5,5
      if (ih .eq. 1) then
        xf(j)=dxpdf(j)/x
      elseif (ih .eq. -1) then
        xf(j)=dxpdf(-j)/x
      else
        write(6,*) 'Unimplemented beam type'
        stop
      endif
      enddo
      return
      end
