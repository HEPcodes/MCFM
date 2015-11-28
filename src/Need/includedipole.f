      logical function includedipole(nd,ptrans)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical gencuts,failedgencuts,makecuts

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
      
      return
      end
            
