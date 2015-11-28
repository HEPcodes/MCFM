************************************************************************
* This routine generates collinear points that satisfy the jet cuts    *
************************************************************************
      subroutine singgen(p,s,*)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'jetlabel.f'
      integer j,jlast1,jlast2,nqcdjets,nqcdstart
      double precision p(mxpart,4),pjet(mxpart,4),s(mxpart,mxpart),rcut
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      
      if     (npart .eq. 3) then
        jlast1=5
        jlast2=4
      elseif     (npart .eq. 4) then
        jlast1=6
        jlast2=5
      elseif (npart .eq. 5) then
        jlast1=7
        jlast2=5
      elseif (npart .eq. 7) then
        jlast1=7
        jlast2=6
      else
        write(6,*) 'singgen.f does not produce points for npart=',npart
        stop
      endif
      
      call genclust2(p,rcut,pjet,0)
      
c--- soft singularity
c      if    ((p(jlast1,4) .lt. 0.01d0)) then
c        write(6,*) 'Generated point with Energy = ',p(jlast1,4)
c        do j=1,npart+2
c          write(6,76) j,p(j,1),p(j,2)
c          write(6,78) p(j,3),p(j,4)
c        enddo
c        pause
c      endif     
c--- initial-final singularity
      if    ((-s(1,jlast1) .lt. 1d-1) .and. (jets .eq. nqcdjets)
     . .and. (-p(1,4) .gt. 1d0) .and. (p(jlast1,4) .gt. 1d0)) then
        write(6,*) 'Generated point with -sij = ',-s(1,jlast1)
        do j=1,npart+2
          write(6,77) j,p(j,1),p(j,2)
          write(6,78) p(j,3),p(j,4)
        enddo
        pause
      endif      
c--- final-final singularity
      if    ((+s(jlast2,jlast1) .lt. 1d-1) .and. (jets .eq. nqcdjets)
     . .and. (p(jlast2,4) .gt. 1d0) .and. (p(jlast1,4) .gt. 1d0)) then
        write(6,*) 'Generated point with +sij = ',+s(jlast2,jlast1)
        do j=1,npart+2
          write(6,79) j,p(j,1),p(j,2)
          write(6,78) p(j,3),p(j,4)
        enddo
        pause
      endif      
      
      return 1

   76 format('      data ps',i1,'/ ',f20.14,'d0,',f20.14,'d0,')
   77 format('      data pi',i1,'/ ',f20.14,'d0,',f20.14,'d0,')
   78 format('     .  ',f20.14,'d0,',f20.14,'d0/')
   79 format('      data pf',i1,'/ ',f20.14,'d0,',f20.14,'d0,')

      end
