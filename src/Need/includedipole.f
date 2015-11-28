      logical function includedipole(nd,ptrans)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      include 'process.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical gencuts,failedgencuts,makecuts
      integer ij
      character*30 runstring
      double precision y32,y43,z3,z4,z5,z6
      double precision dphizj,pt5sq,pt6sq,pt7sq

      common/runstring/runstring
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag
      
      includedipole=.true.

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
      
      call genclust2(ptrans,rcut,pjet,isub)
      do j=1,4
        do i=1,npart+2
        ptildejet(nd,i,j)=pjet(i,j)
        enddo
      enddo
     
c--- if the number of jets is not correct, then do not include dipole
      if ((clustering .and. (jets .ne. nqcdjets-notag)
     .       .and. (inclusive .eqv. .false.)) .or.
     .    (clustering .and. (jets .lt. nqcdjets-notag)
     .       .and. (inclusive .eqv. .true.))) then
          includedipole=.false.
          return
      else
c--- otherwise, if it is correct, check the lepton cuts
        if (makecuts) then
          failedgencuts=gencuts(pjet,jets)
          if (failedgencuts) includedipole=.false.
        endif
      endif
      
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
            
