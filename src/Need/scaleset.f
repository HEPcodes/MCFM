      subroutine scaleset(rscalestart,fscalestart,p)
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'nwz.f'
      include 'facscale.f'
      include 'nlooprun.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      character*52 msg,rmsg,fmsg
      character boson
      integer nqcdjets,nqcdstart,n2,n3,j,iscale,iscalemax,
     . irscale,ifscale,iset
      double precision amz,rscalestart,fscalestart,
     . dot,pttwo,aveptjet,getEt,p(mxpart,4),alphas,
     . setscale
      logical first
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/nqcdjets/nqcdjets,nqcdstart
      common/couple/amz
      data first/.true./  
      save first

c--- convert double precision starting values of the renormalization
c--- and factorization scales to integers
      irscale=int(rscalestart+0.5d0)
      ifscale=int(fscalestart+0.5d0)
      
c--- if the dynamic scales are the same for renormalization and
c--- factorization, then we only need to do the calculation once
      if (irscale .eq. ifscale) then
        iscalemax=1
      else
        iscalemax=2
      endif
      
c--- calculate the dynamic scales
      do iscale=1,iscalemax
      if (iscale .eq. 1) iset=irscale
      if (iscale .eq. 2) iset=ifscale       

      if     (iset .eq. 1) then
c--- sqrt(boson mass + boson pt)
        if (nwz .eq. 0) then
          setscale=dsqrt(zmass**2+pttwo(3,4,p)**2)
          boson='z'
        else
          setscale=dsqrt(wmass**2+pttwo(3,4,p)**2)
          boson='w'
        endif  
        msg='*     Dynamic scale = dsqrt('//boson//'mass**2'//
     .    ' + pt_'//boson//'**2)    *'
      elseif (iset .eq. 2) then
c--- average jet pt      
        if (nqcdjets .eq. 0) then
          write(6,*) 'Invalid choice of scale - no jets!'
          stop
        endif
        setscale=aveptjet(p)
        msg='*          Dynamic scale = < pt_jet >              *'
      elseif (iset .eq. 3) then
c--- sqrt(s-hat)      
        setscale=dsqrt(2d0*dot(p,1,2))
        msg='*              Scale = sqrt(s-hat)                 *'
      elseif (iset .eq. 4) then
c--- HT      
        setscale=0d0
        do j=3,2*(1+n2+n3)+max(0,jets)
        setscale=setscale+getEt(ptildejet(0,j,4),ptildejet(0,j,1),
     .                    ptildejet(0,j,2),ptildejet(0,j,3))
        enddo
	facscale=scale
	if (first) msg='*      Dynamic scale = scalar Et sum (HT)'//
     . '          *'
      elseif ((iset .lt. 1) .or. (iset .gt. 4)) then
c--- catch invalid inputs
        write(6,*) 'Invalid dynamic scale!'
	stop
      endif
      
      if      (iscale .eq. 1) then
        scale=setscale
	rmsg=msg
      elseif ((iscale .eq. 2) .or. (iscalemax .eq. 1)) then
        facscale=setscale
        fmsg=msg
      endif
      
      enddo
      
      if (first) then
        write(6,*)
        write(6,*)'************** Special scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*                 RENORMALIZATION                  *'
        write(6,*) rmsg
        write(6,*)'*                                                  *'
        write(6,*)'*                  FACTORIZATION                   *'
        write(6,*) fmsg
        write(6,*)'****************************************************'
        first=.false.      
      endif
      

c--- catch absurdly large scales      
      if  (scale .gt. 3000d0) scale=3000d0

c--- run alpha_s
      as=alphas(scale,amz,nlooprun)
        
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2
      
      return
            
      end
