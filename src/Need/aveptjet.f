      double precision function aveptjet(p)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'jetlabel.f'
      integer j,countjet
      character*2 plabel(mxpart)
      double precision p(mxpart,4),pjet(mxpart,4),pt,rcut
      common/plabel/plabel
      common/rcut/rcut
      
      aveptjet=0

c-- cluster jets      
      call genclust2(p,rcut,pjet,0)
      
      countjet=0
      do j=3,npart+2
        if (countjet .eq. jets) goto 99
        if (     (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'pj')
     .      .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          countjet=countjet+1
          aveptjet=aveptjet+pt(j,pjet)
        endif
      enddo
     
   99 continue  

c--- dummy value returned if countjet=0, since this process
c--- must have nqcdjets > 0 - so this point will be dumped anyway  
      if (countjet .eq. 0) then
        aveptjet=10d0
        return
      endif
     
      aveptjet=aveptjet/dfloat(countjet)
      
      return
      end
      
      
