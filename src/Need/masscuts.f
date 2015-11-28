      subroutine masscuts(s,*)
      implicit none
      include 'constants.f'
      include 'npart.f'
      logical first
      double precision bbsqmin,bbsqmax,wsqmin,wsqmax,s(mxpart,mxpart)
      integer nqcdjets,nqcdstart
      common/nqcdjets/nqcdjets,nqcdstart
      common/limits/bbsqmin,bbsqmax,wsqmin,wsqmax      
      data first/.true./

      if (first) then
      first=.false.
      write(6,*)
      write(6,*) '****************** Basic mass cuts *****************'
      write(6,*) '*                                                  *'
      write(6,99) dsqrt(wsqmin),'m34',dsqrt(wsqmax)
      if (nqcdjets .lt. 2) then
      write(6,99) dsqrt(bbsqmin),'m56',dsqrt(bbsqmax)
      else
      write(6,98) dsqrt(bbsqmin),'m(jet1,jet2)',dsqrt(bbsqmax)
      endif
      write(6,*) '****************************************************'
      endif
      
      if (  (s(3,4) .lt. wsqmin) 
     . .or. (s(3,4) .gt. wsqmax))
     .  return 1
      
      if ((npart .gt. 3) .and. (nqcdjets .lt. 2)) then
        if (  (s(5,6) .lt. bbsqmin) 
     .   .or. (s(5,6) .gt. bbsqmax))
     .    return 1
      endif
     
   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')
     
      return
      end

