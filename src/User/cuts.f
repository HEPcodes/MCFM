      logical function cuts(p)
c--- This function returns TRUE if a cut is applied
      implicit none
      include 'constants.f'
      integer j
      character*6 case
      character*5 scuts
      double precision p(mxpart,4),yrap,pt,r,pta(4),ptb(4),
     . dotpr,ptbb,costh_H,costh_max,cos1max,mtopa,mtopb
      double precision  ybmax,ptbmin,yemax,ptemin,ptnmin,rebmin,rbbmin,
     . ptbbmin,ptjetmin,yjetmax,mtopright,mtopwrong,
     . cosnew1,cosnew2,cosnew3,cosphi
      logical first,jet7,jet8,lepton7
      logical madewwcuts,madewzcuts,madewzcuts2,madezbbnncuts
      logical madetotcuts,madehiggscuts,madewbbcuts,madezbbcuts
      common/process/case

      data first/.true./

      save ybmax,ptbmin,yemax,ptemin,ptnmin,rebmin,rbbmin,ptbbmin,
     . ptjetmin,yjetmax,first,costh_max,cos1max,mtopright,mtopwrong,
     . scuts

      if (first) then
      first=.false.
      
      open(unit=55,file='cuts.DAT',status='old',err=99)
      write(6,*) 'cuts.f:Reading cuts from cuts.DAT'
      read(55,*) ptbbmin
      write(6,*) 'ptbbmin',ptbbmin
      read(55,*) rebmin
      write(6,*) 'rebmin',rebmin
      read(55,*) rbbmin
      write(6,*) 'rbbmin',rbbmin
      read(55,*) costh_max
      write(6,*) 'costh_max',costh_max
      read(55,*) cos1max
      write(6,*) 'cos1max',cos1max
      read(55,*) ybmax
      write(6,*) 'ybmax',ybmax
      read(55,*) ptbmin
      write(6,*) 'ptbmin',ptbmin
      read(55,*) yemax
      write(6,*) 'yemax',yemax
      read(55,*) ptemin
      write(6,*) 'ptemin',ptemin
      read(55,*) ptnmin
      write(6,*) 'ptnmin',ptnmin
      read(55,*) ptjetmin
      write(6,*) 'ptjetmin',ptjetmin
      read(55,*) yjetmax
      write(6,*) 'yjetmax',yjetmax
      read(55,*) mtopright
      write(6,*) 'mtopright',mtopright
      read(55,*) mtopwrong
      write(6,*) 'mtopwrong',mtopwrong
      read(55,90) scuts
      write(6,*) 'scuts ',scuts
      close(unit=55)
      endif
   90 format(a5)



      cuts=.false.

c--- Check for Special CUTS     
      if (scuts .eq. 'higgs') then
        if (madewwcuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'match') then
        if (madewzcuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'matc2') then
        if (madewzcuts2(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'total') then
        if (madetotcuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'higww') then
        if (madehiggscuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'zbbnn') then
        if (madezbbnncuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'wbbar') then
        if (madewbbcuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      elseif (scuts .eq. 'zbbar') then
        if (madezbbcuts(p)) then
          return
        else
          cuts=.true.
          return
        endif
      endif
     
      ptbb=sqrt((p(5,1)+p(6,1))**2+(p(5,2)+p(6,2))**2)
      call wconstruct(p,costh_h,cosnew1,cosnew2,cosnew3,cosphi)
       
      do j=1,4
      pta(j)=p(6,j)+p(7,j)+p(5,j)
      ptb(j)=p(6,j)+p(7,j)+p(6,j)
      enddo
      mtopa=dotpr(pta,pta)
      mtopb=dotpr(ptb,ptb)

      if (case .eq. 'tt_bbl') then
       lepton7=(abs(pt(7,p)).gt.ptemin).and.(abs(yrap(7,p)) .lt. yemax)
       if (lepton7) cuts=.true.
      elseif (case .eq. 'tt_bbh') then
       jet7=(abs(pt(7,p)).gt.ptjetmin).and.(abs(yrap(7,p)) .lt. yjetmax)
       jet8=(abs(pt(8,p)).gt.ptjetmin).and.(abs(yrap(8,p)) .lt. yjetmax)
        if ((jet7) .or. (jet8)) cuts=.true.
      else
       jet7=(abs(pt(7,p)).gt.ptjetmin).and.(abs(yrap(7,p)) .lt. yjetmax)
       if (jet7) cuts=.true.
      endif


 
      if (
c     . (abs(mtopa) .gt. mtopright**2) .or.
c     . (abs(mtopb) .gt. mtopwrong**2) .or.
     . (abs(yrap(5,p)) .gt. ybmax) .or.
     . (abs(pt(5,p)) .lt. ptbmin) .or.
     . (abs(yrap(6,p)) .gt. ybmax) .or.
     . (abs(pt(6,p)) .lt. ptbmin) .or.
     . (abs(yrap(4,p)) .gt. yemax) .or.
     . (abs(pt(4,p)) .lt. ptemin) .or. 
     . (abs(pt(3,p)) .lt. ptnmin) .or.
     . (r(p,5,6) .lt. rbbmin) .or.
     . (r(p,4,5) .lt. rebmin) .or.
     . (r(p,4,6) .lt. rebmin) .or.
     . (ptbb .lt. ptbbmin ) .or.
     . (abs(costh_h) .gt. costh_max) .or.
     . (cosnew1 .gt. cos1max) 
     . ) then
c      write(*,*) 'failed cuts'
c      write(*,*) 'yrapb',abs(yrap(5,p)),' max',ybmax
c      write(*,*) 'ptb',abs(pt(5,p)),' min',ptbmin
c      write(*,*) 'yrapbb',abs(yrap(6,p)),' max',ybmax
c      write(*,*) 'ptbb',abs(pt(6,p)),' min',ptbmin
c      write(*,*) 'yrape',abs(yrap(4,p)),' max',yemax
c      write(*,*) 'pte',abs(pt(4,p)),' min',ptemin
c      write(*,*) 'pt_nu',abs(pt(3,p)),' min',ptnmin
c      write(*,*) 'rbbbar',r(p,5,6),' min',rbbmin
c      write(*,*) 'reb',r(p,5,4),' min',rebmin
c      write(*,*) 'rebb',r(p,6,4),' min',rebmin
c      write(*,*) 'ptbb',ptbb,' min',ptbbmin
c      write(*,*) 'angle',abs(costh_h),' max',costh_max
c      write(*,*) 'angle',cosnew1,' max',cos1max 
c      pause      
        cuts=.true.
      endif
       
      return

 99   write(6,*) 'error reading cuts.DAT'
      stop
      end


      double precision function ayrap(j,p)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      if (dabs(p(j,4)) .lt. 1d-13) then
C--set to 100 if Energy vanishes
      ayrap=100d0
      else 
      ayrap=0.5d0*dabs(log((p(j,4)+p(j,3))/(p(j,4)-p(j,3))))
      endif
      return
      end

      double precision function yrap(j,p)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4)
      if (dabs(p(j,4)) .lt. 1d-13) then
C--set to 100 if Energy vanishes
      yrap=100d0
      else 
      yrap=0.5d0*log((p(j,4)+p(j,3))/(p(j,4)-p(j,3)))
      endif
      return
      end

      double precision function pt(j,p)
      implicit none
      include 'constants.f'
      integer j
      double precision p(mxpart,4),dot
      pt=sqrt(2d0*dot(p,1,j)*dot(p,2,j)/dot(p,1,2))
      return
      end

