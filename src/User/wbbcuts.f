      logical function madewbbcuts(p)
C--returns true if the event is to be passed
C--cuts appropriate for W(nu e) H(b bbar) and backgrounds 
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),ptminb,ymaxb,misset,etmiss,
     . missetmin,ptb,ptbb,yb,ybb,pt,ayrap,ptemin,yemax,ptjetmin,yjetmax,
     . Rmin,cosmax,costh_h,cosnew1,cosnew2,cosnew3,cosphi,R,
     . ptvetoh1,ptvetoh2,ptvetol,etvec(4)
      logical noextrajets,leptspassed,seppassed,vetopassed
      character plabel(mxpart)*2
      integer i
      common/plabel/plabel
      logical first 
      parameter (ptminb=15d0,ymaxb=2d0,missetmin=20d0)
      parameter (ptemin=20d0,yemax=2.5d0)
c      These cuts are for the tau case
c      parameter (ptemin=30d0,yemax=2d0)
      parameter (ptjetmin=20d0,yjetmax=2d0)
      parameter (ptvetoh1=30d0,ptvetoh2=15d0,ptvetol=10d0)
      parameter (Rmin=0.7d0,cosmax=0.8d0)
      data first/.true./
      save first
      
      if (first) then
        write(*,*) '************************************'
        write(*,*) 'SPECIAL CUTS OVER-RIDING THOSE ABOVE'
        write(*,*) '************************************'
        first=.false.
      endif

      madewbbcuts=.false.
      misset=etmiss(p,etvec)
      ptb=pt(5,p)
      yb=ayrap(5,p)
      ptbb=pt(6,p)
      ybb=ayrap(6,p)
      
      leptspassed=.true.
c--- perform cuts on leptons in the event
      do i=4,4
        if (((plabel(i) .eq. 'el') .or. (plabel(i) .eq. 'ea')) .and.
     .      ((pt(i,p) .gt. ptemin) .and. (ayrap(i,p) .lt. yemax))) then
        continue
        else
        leptspassed=.false.
        endif        
      enddo
        if (plabel(7) .ne. 'el') goto 20
        if ((pt(7,p) .gt. ptemin) .and. (ayrap(7,p) .lt. yemax))
     .  leptspassed=.false.
 20     continue
       
      noextrajets=.true.
c--- perform cuts on extra jets in the event (eg. process 152)
c      if (nproc .eq. 152) then 
c        do i=7,8
c        if ((pt(i,p) .gt. ptjetmin) .and. (ayrap(i,p) .lt. yjetmax))
c     .    noextrajets=.false.
c        enddo
c      elseif (nproc .eq. 161) then 
c        do i=7,7
c        if ((pt(i,p) .gt. ptjetmin) .and. (ayrap(i,p) .lt. yjetmax))
c     .    noextrajets=.false.
c        enddo
c      else
c        do i=7,4+nqcdjets
c        if ((pt(i,p) .gt. ptjetmin) .and. (ayrap(i,p) .lt. yjetmax))
c     .    noextrajets=.false.
c        enddo
c      endif

c--- perform cuts on b-bbar separation
      seppassed=(R(p,5,6) .gt. Rmin)

c--- jet and lepton veto for t-tbar events (151 and 152)
c      vetopassed=.true.
c      if (plabel(7) .eq. 'el') then
c        if (((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax))
c     ..and. ((pt(7,p) .gt. ptemin) .and. (ayrap(7,p) .lt. yemax))) then
c          vetopassed=.false.
c        endif         
c        if (((pt(4,p) .gt. ptemin) .and. (ayrap(4,p) .lt. yemax))
c     ..and. ((pt(7,p) .lt. ptemin) .or. (ayrap(7,p) .gt. yemax))) then
c          leptspassed=.true.
c        endif         
c        if (((pt(7,p) .gt. ptemin) .and. (ayrap(7,p) .lt. yemax))
c     ..and. ((pt(4,p) .lt. ptemin) .or. (ayrap(4,p) .gt. yemax))) then
c          leptspassed=.true.
c        endif         
c      endif
c      if ((plabel(7) .eq. 'pp') .and. (plabel(8) .eq. 'pp')) then
c        if ((pt(7,p).gt.ptvetoh1) .or. (pt(8,p).gt.ptvetoh1)) then
c          vetopassed=.false.
c        elseif ((pt(7,p).gt.ptvetoh2) .and. (pt(8,p).gt.ptvetoh2)) then
c          vetopassed=.false.
c        endif
c      endif  
         
      call wconstruct(p,costh_h,cosnew1,cosnew2,cosnew3,cosphi)
      
C--- combine the above cuts with a cut on pt and rapidities of b-jets,
C--- and cos(theta)
      if (   (ptb .gt. ptminb) 
     . .and. (ptbb .gt. ptminb)
     . .and. (yb  .lt. ymaxb)
     . .and. (ybb .lt. ymaxb)
     . .and. (dabs(costh_h) .lt. cosmax)
     . .and. (misset .gt. missetmin)
     . .and. (leptspassed) 
     . .and. (noextrajets) 
     . .and. (seppassed) 
     . ) madewbbcuts=.true.


      return
      end
