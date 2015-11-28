      logical function madezbbcuts(p)
C--returns true if the event is to be passed
C--cuts appropriate for Z(nu nu) b bbar
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),ptminb,ymaxb,misset,etmiss,Rmin,R,
     . missetmin,ptb,ptbb,yb,ybb,pt,ayrap,ptemin,yemax,ptjetmin,yjetmax,
     . etvec(4),deltaphib,deltaphibb,azimin
      logical jetspassed,seppassed,azipassed,nolepts,noextrajets
      character plabel(mxpart)*2
      integer i
      common/plabel/plabel
      logical first 
      parameter (ptminb=15d0,ymaxb=2d0)
      parameter (missetmin=35d0,azimin=0.5d0)
      parameter (Rmin=0.7d0)
      parameter (ptemin=10d0,yemax=2d0)
      parameter (ptjetmin=15d0,yjetmax=2.5d0)
      data first/.true./
      save first
      
      if (first) then
        write(*,*) '************************************'
        write(*,*) 'SPECIAL CUTS OVER-RIDING THOSE ABOVE'
        write(*,*) '************************************'
        first=.false.
      endif
      
      madezbbcuts=.false.
      misset=etmiss(p,etvec)
      ptb=pt(5,p)
      yb=ayrap(5,p)
      ptbb=pt(6,p)
      ybb=ayrap(6,p)

c--- perform cuts on b-bbar separation
      seppassed=(R(p,5,6) .gt. Rmin)

c--- cut on azimuth between missing Et and jets
c      deltaphib =dacos((p(5,1)*etvec(1)+p(5,2)*etvec(2))/(ptb *misset))
c      deltaphibb=dacos((p(6,1)*etvec(1)+p(6,2)*etvec(2))/(ptbb*misset))
      deltaphib =(p(5,1)*etvec(1)+p(5,2)*etvec(2))/(ptb *misset)
      if (deltaphib .gt. 1d0) then
        deltaphib=0d0
        write(*,*) 'Argument of cos > 1'
      elseif (deltaphib .lt. -1d0) then
        deltaphib=pi
        write(*,*) 'Argument of cos < -1'
      else
        deltaphib=dacos(deltaphib)
      endif

      deltaphibb=(p(6,1)*etvec(1)+p(6,2)*etvec(2))/(ptbb*misset)
      if (deltaphibb .gt. 1d0) then
        deltaphibb=0d0
        write(*,*) 'Argument of cos > 1'
      elseif (deltaphibb .lt. -1d0) then
        deltaphibb=pi
        write(*,*) 'Argument of cos < -1'
      else
        deltaphibb=dacos(deltaphibb)
      endif
      azipassed=((deltaphib .gt. azimin) .and. (deltaphibb. gt. azimin))

c--- cut for invisible leptons (only 151, 152, 161 and 171)    
      if (((plabel(4) .eq. 'el') .or. (plabel(4) .eq. 'ea')) .and.
     .    ((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax))) then
        nolepts=.false.
      else
        nolepts=.true.
      endif
      
c--- additional cut for invisible leptons (151 only)    
      if (((plabel(7) .eq. 'el') .or. (plabel(7) .eq. 'ea')) .and.
     .    ((pt(7,p) .gt. ptemin) .and. (ayrap(7,p) .lt. yemax))) then
        nolepts=.false.
      endif
      
c--- add in extra Et if leptons are missed
      if (  ((plabel(4) .eq. 'el') .or. (plabel(4) .eq. 'ea')) .and.
     .      (nolepts)) then
        etvec(1)=etvec(1)+p(4,1)
        etvec(2)=etvec(2)+p(4,2)
        misset=dsqrt(etvec(1)**2+etvec(2)**2)
      endif
      if (  ((plabel(7) .eq. 'el') .or. (plabel(7) .eq. 'ea')) .and.
     .      (nolepts)) then
        etvec(1)=etvec(1)+p(7,1)
        etvec(2)=etvec(2)+p(7,2)
        misset=dsqrt(etvec(1)**2+etvec(2)**2)
      endif
           
      
c--- cut for invisible jets (only 161)    
c      if ((plabel(7) .eq. 'pj') .and. ((pt(7,p) .gt. ptjetmin)
c     .     .and. (ayrap(7,p) .lt. yjetmax))) then
c        noextrajets=.false.
c      else
        noextrajets=.true.
c      endif
      
c--- put it all together
      if (   (ptb .gt. ptminb) 
     . .and. (ptbb .gt. ptminb)
     . .and. (yb  .lt. ymaxb)
     . .and. (ybb .lt. ymaxb)
     . .and. (misset .gt. missetmin)
     . .and. (seppassed) 
     . .and. (azipassed) 
     . .and. (nolepts) 
     . .and. (noextrajets) 
     . ) madezbbcuts=.true.

c      write(*,*) ptb,yb
c      write(*,*) ptbb,ybb
c      write(*,*) misset
c      write(*,*) seppassed
c      write(*,*) azipassed
c      write(*,*) nolepts
c      write(*,*) 'PASSED? ',madezbbcuts
c      pause

      return
      end
