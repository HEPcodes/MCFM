************************************************************************
* This routine checks the singularity cancellation, given the two      *
* appropriate routines - real_g, sub_gs - and the momentum array p.    *
* A given sub-process (j,k) is determined to be singular if            *
*  (npart=4) msq(j,k) > 1d3 * (smallest value of msq(j,k))             *
*  (npart=5) msq(j,k) > 1d4 * (smallest value of msq(j,k))             *
* This seems to be a reasonable criterion for this set of singular     *
* points. The singularity is considered uncancelled if                 *
*  msq(j,k)/(sum of dipoles) is more than 3% different from 1          *
************************************************************************
      subroutine singcheck(real_g,sub_gs,p)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'ptilde.f'
      integer j,k,jj,kk,jmax,kmax,nd,jets,nqcdjets,nqcdstart
      double precision p(mxpart,4),q(mxpart,4),pjet(mxpart,4),
     . msq(-nf:nf,-nf:nf),msqs(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     . s(mxpart,mxpart),rcut,debugsmall,debuglarge,xtoler
      character*32 debugmsg
      character jetlabel(mxpart)*2
      external real_g,sub_gs
      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      
      if     (npart .eq. 3) then
        jmax=5
        kmax=200
        xtoler=1d3
      elseif (npart .eq. 4) then
        jmax=5
        kmax=5
      elseif (npart .eq. 5) then
        jmax=5
        kmax=9
        xtoler=1d3
      else
        write(6,*) 'No singularity check implemented yet'
        stop
      endif
      
      do j=1,jmax
      do k=1,kmax
         if     (npart .eq. 3) then
           call coll3(p,k,j)
         elseif (npart .eq. 4) then
           call coll4a(p,k,j)
         elseif (npart .eq. 5) then
           call coll5(p,k,j)
         endif
         write(*,*) 'Point ',j,':'
         call real_g(p,msq)  
         call sub_gs(p,msqc) 
         do jj=-nf,nf
         do kk=-nf,nf
            msqs(jj,kk)=0d0
         enddo
         enddo
         call dotem(7,p,s)
         call genclust2(p,rcut,jets,pjet,jetlabel)
         if (jets .eq. -1) then
           write(6,*) 'This point does not have a final state b-jet'
           do jj=-nf,nf
           do kk=-nf,nf
              msq(jj,kk)=0d0
           enddo
           enddo
           goto 68
         endif
c         write(6,*) '0',jets,msq(2,0),(pjet(5,4)+pjet(6,4))**2
c     .                      -(pjet(5,1)+pjet(6,1))**2
c     .                      -(pjet(5,2)+pjet(6,2))**2
c     .                      -(pjet(5,3)+pjet(6,3))**2
         do nd=1,ndmax
            do jj=1,7
            do kk=1,4
               q(jj,kk)=ptilde(nd,jj,kk)
            enddo
            enddo
            call dotem(7,q,s)
            call genclust2(q,rcut,jets,pjet,jetlabel)
            if (jets .eq. nqcdjets) then
c            write(6,*) nd,jets,msqc(nd,2,0),(pjet(5,4)+pjet(6,4))**2
c     .                        -(pjet(5,1)+pjet(6,1))**2
c     .                        -(pjet(5,2)+pjet(6,2))**2
c     .                        -(pjet(5,3)+pjet(6,3))**2
               do jj=-nf,nf
               do kk=-nf,nf
                  msqs(jj,kk)=msqs(jj,kk)+msqc(nd,jj,kk)
               enddo
               enddo
            endif
         enddo

c--- find smallest value of msq
         debugsmall=1d0
         debuglarge=0d0
         do jj=-nf,nf
         do kk=-nf,nf
           if ((msq(jj,kk) .lt. debugsmall)
     .     .and. (msq(jj,kk) .gt. 0d0)) debugsmall=msq(jj,kk)        
           if ((msq(jj,kk) .gt. debuglarge)
     .     .and. (msq(jj,kk) .gt. 0d0)) debuglarge=msq(jj,kk)        
         enddo
         enddo                  
      
         if ((debuglarge/debugsmall .lt. xtoler)
     .  .or. (debuglarge .lt. 1d-10)) then
c           write(6,*) ' OK - no singular configurations'
c           goto 68
         endif
      
         if (debugsmall .gt. 1d-3) debugsmall=debugsmall/xtoler/2d0
c         if (debuglarge/debugsmall .lt. xtoler*1d2) 
c     .       debugsmall=debuglarge/xtoler
         do jj=-nf,nf
         do kk=-nf,nf
         if (msq(jj,kk) .eq. 0d0) then
            if (msqs(jj,kk) .eq. 0d0) then
               debugmsg='   OK   zero msq and subtraction'
            else
               debugmsg=' FAILED subtraction with msq=0'
            endif
            write(*,22) jj,kk,'    n/a   ',msq(jj,kk),msqs(jj,kk),
     .      debugmsg
            goto 69
            endif

         if ((msq(jj,kk)/debuglarge) .lt. 1d0/xtoler) then
            if (abs(abs(msq(jj,kk)/msqs(jj,kk))-1d0) .lt. 0.02d0)then
              debugmsg='   OK   cancelled'
            else
              debugmsg='   OK   not singular'
            endif
            write(*,21) jj,kk,msq(jj,kk)/msqs(jj,kk),
     .                    msq(jj,kk),msqs(jj,kk),debugmsg
         else
             if (msqs(jj,kk) .eq. 0d0) then
                debugmsg=' FAILED singular, no subtraction'
                write(*,22) jj,kk,'    n/a   ',msq(jj,kk),msqs(jj,kk),
     . debugmsg
                goto 69
             endif
             if (abs(abs(msq(jj,kk)/msqs(jj,kk))-1d0) .gt. 0.02d0)then
                debugmsg=' FAILED singularity uncancelled'
             else
                debugmsg='   OK   cancelled'
             endif
             write(*,21) jj,kk,msq(jj,kk)/msqs(jj,kk),
     .                  msq(jj,kk),msqs(jj,kk),debugmsg
         endif
      
   69    continue
         enddo
         enddo

   68 continue
      enddo
      enddo 
      
      write(6,*) 'Singularity check completed'
      write(6,*) 'Do ''mcfm | grep -v OK'' to look for failures'
      stop
      
      return
      
   21 format(1x,2i3,f10.6,2e14.6,1x,a32)
   22 format(1x,2i3,a10,2e14.6,1x,a32)

      end
