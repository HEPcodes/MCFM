      logical function includedipole(nd,ptrans)
c--- This function returns TRUE if the specified point ptrans,
c--- corresponding to dipole nd (nd=0 => real radiation),
c--- should be included 
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'process.f'
      include 'frag.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical gencuts,failedgencuts,photoncuts,makecuts,filterWbbmas,
     & photonfailed
      integer ij,count_photo,nphotons
      character*30 runstring
      double precision y32,y43,z3,z4,z5,z6
      double precision dphizj,pt5sq,pt6sq,pt7sq
      logical passed_frix,iso
      logical phot_dip(mxpart)
      character*2 plabel(mxpart)

      common/phot_dip/phot_dip     
      common/runstring/runstring
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag
      common/plabel/plabel

c--- default: include this contribution
      includedipole=.true.

c--- isub=1 for dipole subtractions, isub=0 for real radiation
      if (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

c--- special procedure for comparison with C. Oleari's
c---  e+e- --> QQbg calculation; run Durham algorithm
      if (runstring(1:5) .eq. 'carlo') then
        call durhamalg(ptrans,npart-isub,y32,y43,z3,z4,z5,z6)
c---     set up momentum array just in case it's needed elsewhere
        do j=1,4
          do i=1,npart+2
          ptildejet(nd,i,j)=ptrans(i,j)
          enddo
        enddo
	jets=1
        if (y43 .gt. 0d0) jets=2
	if     ((case .eq. 'qq_tbg') .or. (case .eq. 'epem3j')) then
c---     only keep events that have 3 jets when ycut=rcut
	  if ((y32 .lt. rcut) .or. (y43 .gt. rcut)) then
	    includedipole=.false.
	    jets=-1
	  endif
	  return
	elseif (case .eq. 'qqtbgg') then
c---     only keep events that have 4 jets when ycut=rcut
	  if (y43 .lt. rcut) then
	    includedipole=.false.
	    jets=-1
	  endif
	  return
	else
	  write(6,*) 'Unexpected case in includedipole.f!'
	  stop
	endif
      endif

      nphotons=count_photo(ptrans) 
      if (nphotons .gt. 0) then 
c--- Photons: Frixione isolation cuts if no fragmentation included 
        if (frag .eqv. .false.) then
          if    ((case .eq. 'dirgam') .or. (case .eq. 'gamjet'))then 
            call basic_frix(ptrans,passed_frix,isub)
          elseif((case .eq. 'Wgamma') .or. (case .eq. 'Zgamma')
     &       .or.(case .eq. 'Wgajet') .or. (case .eq. 'Zgajet')) then
            call basic_Vfrix(ptrans,passed_frix,isub)
          elseif((case .eq. 'gamgam') .or. (case .eq. 'Higaga')
     &     .or.  (case .eq. 'gmgmjt')) then 
            call basic_di_frix(ptrans,passed_frix,isub)
	  else
	    write(6,*) 'Unforeseen photon process in includedipole.f'
	    write(6,*) 'Needs a new Frixione routine.'
	    stop
          endif
          if (passed_frix .eqv. .false.) then 
            includedipole=.false.
            return
          else
            call genclustphotons(ptrans,rcut,pjet,isub)	       
          endif
        else 
c--- Photons: not Frixione, need fragmentation and isolation
          call genclust2(ptrans,rcut,pjet,isub)
!---  Isolate photon
          do j=3,mxpart
            if (plabel(j).eq.'ga') then 
              if (iso(ptrans,isub,phot_dip(nd),j,nd) .eqv. .false.) then
          	includedipole=.false.
          	return 
              endif
            endif
          enddo
        endif
      else
c--- No photons: the usual case
         call genclust2(ptrans,rcut,pjet,isub)
      endif

c--- fill ptilde array as persistent storage for the jet momenta
      do j=1,4
        do i=1,npart+2
          ptildejet(nd,i,j)=pjet(i,j)
        enddo
      enddo
     
c--- for the Wbb process, we divide up contributions in a specific way;
c--- therefore filter events using special code and skip normal jet testing   
      if (case .eq. 'Wbbmas') then
        includedipole=filterWbbmas()
	if (includedipole .eqv. .false.) return
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) includedipole=.false.
        endif
	goto 99
      endif
      

c--- if the number of jets is not correct, then do not include dipole
      if ((clustering .and. (jets .ne. nqcdjets-notag)
     &       .and. (inclusive .eqv. .false.)) .or.
     &    (clustering .and. (jets .lt. nqcdjets-notag)
     &       .and. (inclusive .eqv. .true.))) then
          includedipole=.false.

          return
      else
c--- otherwise, if it is correct, check the photon cuts if appropriate
c--- rescale photons for fragmentation processes
c	  if (rescale) call rescale_pjet(pjet)
        if (nphotons .gt. 0) then 
	  if (rescale) call rescale_pjet(pjet)
          photonfailed=photoncuts(pjet,isub,phot_dip(nd),nd)
c--- return photons to original values
	  if (rescale) call return_pjet(pjet)
	  if (photonfailed) then
	    includedipole=.false.
	    return
	  endif
	endif
c--- ... and then the lepton cuts, if necessary
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)	  
          if (failedgencuts) then
	    includedipole=.false.
	    return
	  endif
c          call writeout(ptrans)
c          pause
        endif
      endif
      
c      call Wgamcuts(ptrans,failedgencuts)
c      if (failedgencuts .eqv. .false.) includedipole=.false.
      
     
   99 continue
      
c--- added extra check here, to allow for analysis of G. Hesketh et al.
c--- that requires Z+2 jets with only one jet within cuts, to obtain
c--- prediction for Delta_phi(Z,jet) at NLO
      if (runstring(1:6) .eq. 'dphizj') then
        if    (nqcdjets .eq. 1) then
	  ij=5
	elseif (nqcdjets .eq. 2) then
	  pt5sq=pjet(5,1)**2+pjet(5,2)**2
	  pt6sq=pjet(6,1)**2+pjet(6,2)**2
	  if (pt5sq .gt. pt6sq) then
            ij=5
	  else
	    ij=6
	  endif
	elseif (nqcdjets .eq. 3) then
	  pt5sq=pjet(5,1)**2+pjet(5,2)**2
	  pt6sq=pjet(6,1)**2+pjet(6,2)**2
	  pt7sq=pjet(7,1)**2+pjet(7,2)**2
	  if     (pt5sq .gt. max(pt6sq,pt7sq)) then
            ij=5
	  elseif (pt6sq .gt. max(pt5sq,pt7sq)) then
	    ij=6
	  else
	    ij=7
	  endif
	endif
        dphizj=atan2(pjet(3,1)+pjet(4,1),pjet(3,2)+pjet(4,2))
     .        -atan2(pjet(ij,1),pjet(ij,2))
        if (dphizj .gt. pi) dphizj=twopi-dphizj
        if (dphizj .lt. -pi) dphizj=twopi+dphizj
	if (abs(dphizj) .gt. pi-1d-1) includedipole=.false.
      endif
      
      return
      end
            
