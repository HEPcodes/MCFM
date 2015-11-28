      subroutine scaleset(scalestart,p)
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'nwz.f'
      include 'facscale.f'
      include 'nlooprun.f'
      character*52 msg
      integer nproc,facindex
      character*4 part
      character boson
      integer nqcdjets,nqcdstart
      common/nqcdjets/nqcdjets,nqcdstart
      common/part/part
      common/nproc/nproc
      common/couple/amz
      double precision amz,scalestart,p(mxpart,4),alphas,
     . dot,pt,pttwo,ptb,facsf(8),aveptjet
      logical first  
      data facsf/0.25d0,0.33333d0,0.5d0,0.75d0,1d0,2d0,3d0,4d0/
      data first/.true./  
      save first,facsf

      msg=''

c--- factorization scale is also set here for the bg -> Hb process
      if (nproc/10 .eq. 14) then
        if     (scalestart .lt. 0d0) then
          scale=hmass
          if (scalestart .gt. -9d0) then
            facindex=int(-scalestart)
            facscale=scale*facsf(facindex)
          else
            facscale=scale
          endif
        elseif (scalestart .lt. 9d0) then
            facscale=hmass 
            facindex=int(scalestart)
            scale=facscale*facsf(facindex)
        else
          facscale=scale
        endif
        if (first) then
        write(6,*)
        write(6,*)'*** Special renormal/factor-ization scale choices **'
        write(6,*)'*                                                  *'
        write(6,76) scale
        write(6,77) facscale
        write(6,*)'****************************************************'
        first=.false.      
        endif
        goto 99
      endif

c--- hack for Wbb + jet
c      if ((nproc .eq. 24) .or. (nproc .eq. 29)) then
c        scale=dsqrt(wmass**2+(pt(5,p)**2+pt(6,p)**2+pt(7,p)**2)/3d0)
c        as=alphas(scale,amz,2)
c        goto 100
c      endif
      
c--- hack for Wbb
c      if (  (nproc .eq. 20) .or. (nproc .eq. 21)
c     . .or. (nproc .eq. 25) .or. (nproc .eq. 26)) then
c        scale=dsqrt(wmass**2+(pt(5,p)**2+pt(6,p)**2)/2d0)
c        as=alphas(scale,amz,2)
c        goto 100
c      endif
      
c--- hack for Zbb
c      if (  (nproc .eq. 50) .or. (nproc .eq. 51)) then
c        scale=dsqrt(zmass**2+(pt(5,p)**2+pt(6,p)**2)/2d0)
c        as=alphas(scale,amz,2)
c        goto 100
c      endif
      
c--- hack for Z+jet
c      if (nproc .eq. 41) then
c        scale=dsqrt(zmass**2+pt(5,p)**2)
c        as=alphas(scale,amz,2)
c        goto 100
c      endif
      
c--- hack for W+2 jet
c      if (  (nproc .eq. 22) .or. (nproc .eq. 27)) then
c        scale=dsqrt(wmass**2+(pt(5,p)**2+pt(6,p)**2)/2d0)
c        as=alphas(scale,amz,2)
c        goto 100
c      endif
      
c--- universal scale choices
      if     (scalestart .eq. -1d0) then
        if (nwz .eq. 0) then
          scale=zmass
          boson='z'
        else
          scale=wmass
          boson='w'
        endif  
        msg='*              Dynamic scale = '//
     .    boson//'mass'//'               *'
      elseif (scalestart .eq. -2d0) then
        if (nwz .eq. 0) then
          scale=dsqrt(zmass**2+pttwo(3,4,p)**2)
          boson='z'
        else
          scale=dsqrt(wmass**2+pttwo(3,4,p)**2)
          boson='w'
        endif  
        msg='*     Dynamic scale = dsqrt('//boson//'mass**2'//
     .    ' + pt_'//boson//'**2)    *'
      elseif (scalestart .eq. -3d0) then
        if (nqcdjets .eq. 0) then
          write(6,*) 'Invalid choice of scale - no jets!'
          stop
        endif
        scale=aveptjet(p)
        msg='*          Dynamic scale = < pt_jet >              *'
      elseif (scalestart .eq. -4d0) then
        if (nwz .eq. 0) then
          scale=zmass/4d0
          boson='z'
        else
          scale=wmass/4d0
          boson='w'
        endif  
        msg='*            Dynamic scale = '
     .    //boson//'mass/4'//'               *'
      elseif((nproc .eq. 18) .or. (nproc .eq. 19)) then
        scale=dsqrt(wmass**2+(pttwo(3,4,p)**2+pt(5,p)**2)/2d0)
        msg='*    Dynamic scale = dsqrt(wmass**2 + avg. pts)    *'
      elseif((nproc .ge. 70) .and. (nproc .le. 79)) then
        scale=0.5d0*(wmass+zmass)
        if (first) write(*,79) scale
      elseif((nproc .ge. 60) .and. (nproc .le. 69)) then
        if (scalestart .eq. -100d0) then
          scale=dsqrt((pttwo(3,4,p)**2+pttwo(5,6,p)**2
     .                +2d0*wmass**2)/2d0)
          msg='* Dynamic scale = dsqrt(avg. of masses and pt**2)  *'
        else
          if (first) write(*,79) scale
          scale=wmass
        endif
      elseif((nproc .ge. 80) .and. (nproc .le. 89)) then
        scale=zmass
        if (first) write(*,79) scale
      elseif((nproc .ge. 90) .and. (nproc .le. 96)) then
        scale=hmass
        if (first) write(*,79) scale
      elseif ((nproc .eq. 110) .or. (nproc .eq. 111)) then
        scale=hmass
        if (first) write(*,79) scale
      elseif ((nproc .ge. 140).and.(nproc .le. 149)) then
        if (scalestart .eq. -100d0) then
          scale=hmass
          if (first) write(*,79) scale
        else
          if (part .eq. 'real') then
            ptb=pttwo(5,6,p)
          else
            ptb=pt(5,p)
          endif
          scale=dsqrt(hmass**2+ptb**2)
          if (first) write(6,78) dabs(scalestart)
          scale=scale*dabs(scalestart) 
        endif          
      elseif ((nproc .ge. 150).and.(nproc .le. 152)) then
        scale=100d0
        if (first) write(*,79) scale
      elseif ((nproc .eq. 161)) then
        scale=100d0
        if (first) write(*,79) scale
      elseif ((nproc .eq. 171)) then
        scale=100d0
        if (first) write(*,79) scale
      elseif (nproc .eq. 190) then
        if (scalestart .eq. -1d0) then
C       scale approximating NLO corrections
        scale=0.85d0*(mt+0.5d0*hmass)
        endif
      elseif ((nproc .eq. 11)
     .  .or. (nproc .eq. 16)
     .  .or. (nproc .eq. 20)
     .  .or. (nproc .eq. 21)
     .  .or. (nproc .eq. 51)
     .  .or. (nproc .eq. 52)
     .  .or. (nproc .eq. 22)) then
        if (scalestart .eq. -1d0) then
          scale=wmass
          msg='*              Scale = Mass of the W               *'
        elseif (scalestart .eq. -2d0) then
          scale=wmass/4d0
          msg='*              Scale = wmass/4d0                   *'
        elseif (scalestart .eq. -3d0) then
          scale=dsqrt(2d0*dot(p,1,2))
          msg='*              Scale = sqrt(s-hat)                 *'
        endif
        if (scale .gt. 1000d0) scale=1000d0
        if (scale .lt. 10d0) scale=10d0
      elseif ((nproc .eq. 40)
     .  .or. (nproc .eq. 41)
     .  .or. (nproc .eq. 42)
     .  .or. (nproc .eq. 43)
     .  .or. (nproc .eq. 44)
     .  .or. (nproc .eq. 51)
     .  .or. (nproc .eq. 52)
     .  .or. (nproc .eq. 53)
     . ) then
        if (scalestart .eq. -1d0) then
          scale=zmass
          msg='*               Scale = Mass of the Z              *'
        elseif (scalestart .eq. -2d0) then
          scale=zmass/4d0
          msg='*                Scale = zmass/4d0                 *'
        endif
        if (scale .gt. 1000d0) scale=1000d0
        if (scale .lt. 10d0) scale=10d0
      else
        write(*,*) 'Invalid input: please choose a scale'
        stop
      endif

      if ((first) .and. (msg .ne. '')) then
        write(6,*)
        write(6,*)'************** Special scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*) msg
        write(6,*)'****************************************************'
        first=.false.      
      endif
      
   99 continue

c--- catch absurdly large scales      
      if  (scale .gt. 3000d0) scale=3000d0

c--- run alpha_s
      as=alphas(scale,amz,nlooprun)
      
  100 continue
  
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2
      
      return
      
   76 format(' *      Renormalization scale =',f7.2,'              *')   
   77 format(' *        Factorization scale =',f7.2,'              *')   
   78 format(' *  Dynamic scale = ',f4.2,
     .       ' x sqrt(hmass**2+pt(b)**2)  *')
   79 format(' *               Static scale =',f7.2,'              *')   
      end
