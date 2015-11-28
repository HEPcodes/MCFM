!--- This is a driver function for photon isolation, the relevant photon isolation 
!--- critera are selected and applied 

      logical function iso(p,isub,phot_dip,phot_id,nd)
      implicit none
      include 'constants.f'
      include 'frag.f'
      double precision p(mxpart,4)
      integer isub,nd
      logical phot_dip,photo_iso
      integer phot_id ! refers to which photon we are isolating
      character*30 runstring
      character*2 str
      logical first 
      
      common/runstring/runstring 
      data first/.true./
      save first

      iso=.true.

! Check if no isolation required
      if ((abs(cone_ang) .lt. 1d-4).or.(abs(epsilon_h) .lt. 1d-4)) then
        if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*                                                  *'
        write(6,*)'*         No photon isolation cuts applied         *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
	endif
        return
      endif

!---- Check runstring for special isolation conditions 

!      if    (runstring(1:3).eq.'cdf') then 
!------ CDF isolation (to be written) 
!        write(6,*) 'CDF isolation routine not yet implemented.'
!	stop
!      elseif(runstring(1:2).eq.'D0') then 
!------ D0 isolation (to be written) 
!        write(6,*) 'D0 isolation routine not yet implemented.'
!	stop
!      elseif(runstring(1:5).eq.'atlas') then 
!------ Atlas isolation (to be written) 
!        write(6,*) 'ATLAS isolation routine not yet implemented.'
!	stop
!      elseif(runstring(1:3).eq.'CMS') then 
!------ CMS isolation (to be written) 
!        write(6,*) 'CMS isolation routine not yet implemented.'
!	stop
!      else
!------ standard isolation routine      
!       if(runstring(1:2).eq.'Et') then 
!------ cut on Et rather than Pt 
!         str='Et'
!       else
         str='pt'
!       endif
       iso=photo_iso(p,isub,phot_dip,phot_id,nd,str,cone_ang,epsilon_h)
c--- write out isolation parameters
       if(first.and.(epsilon_h.lt.1d0)) then 
        write(6,*)'************** Photons Isolated     ****************'
        write(6,*)'*                                                  *'
        write(6,99)'* ',str,'(had) in cone',cone_ang,' < ',epsilon_h,
     &   ' ',str,'(phot)      *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
      elseif(first.and.(epsilon_h .ge.1d0)) then
        write(6,*)'************** Photons Isolated     ****************'
        write(6,*)'*                                                  *'
        write(6,96)'* E_t (had) in cone',cone_ang,' < ',epsilon_h,
     &   'GeV    *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
      endif
      
!      endif

      return 

 99   format(1x,a2,a2,a13,f6.2,a4,f6.2,a1,a2,a16)
 96   format(1x,a19,f6.2,a4,f6.2,a17)
      end
