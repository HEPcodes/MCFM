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
      include 'process.f'
      include 'npart.f'
      include 'frag.f'
      character*52 msg,rmsg,fmsg
      character*4 prechar(7)
      character boson
      integer nqcdjets,nqcdstart,n2,n3,j,iscale,iscalemax,
     . irscale,ifscale,iset
      double precision amz,rscalestart,fscalestart,
     . dot,pttwo,aveptjet,getEt,p(mxpart,4),alphas,
     . setscale,pt,prescale(7),prefac
      logical first
      double precision mass2,width2,mass3,width3
      double precision b1scale,q2scale,q1scale,b2scale
      character*4 part
      common/part/part
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/nqcdjets/nqcdjets,nqcdstart
      common/couple/amz
      common/bqscale/b1scale,q2scale,q1scale,b2scale
      data prescale/0.25d0,0.5d0,0.75d0,1d0,2d0,4d0,8d0/
      data prechar/"0.25"," 0.5","0.75"," 1.0"," 2.0"," 4.0"," 8.0"/
      data first/.true./  
      save first

c--- FIRST CHECK FOR SPECIAL CASES
c---    Wgamma and Zgamma: choice of de Florian and Signer
      if       ((case .eq. 'Wgamma') .or. (case .eq. 'Zgamma'))then
         
         prefac=rscalestart
         if(prefac.gt.10d0) then 
            write(6,*) 'Invalid Dynamic scale' 
            stop
         endif
         if (rescale) then
            scale=dsqrt(mass3**2
     &           +(z_frag**2*pt(5,p)**2))	
         else
            scale=dsqrt(mass3**2+(pt(5,p)**2))	
         endif
         scale=prefac*scale
         facscale=scale
         frag_scale=scale
         if (first) then
            write(6,*)
        write(6,*)'************** Special scale choice ****************'             
        write(6,*)'*                                                  *'
        write(6,45)'*        muR = muF = ',prefac,
     &       ' sqrt(MV^2+(pt(ga)^2))  *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
            first=.false.
         endif
         goto 77
         

!----- Gamgam choice of DIPHOX = m(gamma,gamma)
      elseif((case .eq. 'gamgam')) then 
         if(rescale) then ! rescale p(4)  
            scale=z_frag*2d0*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &           -p(4,3)*p(3,3))
            scale=dsqrt(scale)
         else
             scale=2d0*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &           -p(4,3)*p(3,3))
             scale=dsqrt(scale)             
         endif
         prefac=rscalestart 
        
         if(prefac .gt. 10d0) then 
            write(6,*) 'INVALID DYNAMIC SCALE' 
            stop
         endif
	 scale=prefac*scale
         if(scale.lt.1d0) then 
            scale=1d0
         endif
         facscale=scale
         frag_scale=scale
         
        if (first) then
        write(6,*)
        write(6,*)'************** Special scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,45)'*         muR = muF = ',prefac,' M_(gam,gam)      *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
	first=.false.
	endif
 45     format(1x,a21,f6.2,a25)
	goto 77
       


!----- Dirgam choice of JETPHOX = pt(gamma)/2
      elseif((case .eq. 'dirgam').or.(case.eq.'gamjfr')) then 
         if(rescale) then 
            scale=z_frag*pt(3,p)/two
         else
            scale=pt(3,p)/two
         endif
         facscale=scale
         frag_scale=scale
         
        if (first) then
        write(6,*)
        write(6,*)'************** Special scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*              muR = muF = pt(ga)/2                *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
	first=.false.
	endif
	goto 77
     

c---    single top s-channel, dynamic scale is (pt+pb)**2
      elseif   (case .eq. 't_bbar') then
        scale=dsqrt(abs(
     .       +(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     .       -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     .       -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     .       -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2))
        facscale=scale
	if (first) then
        write(6,*)
        write(6,*)'************** Special scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*        muR = muF =  (p(top)+p(bottom))^2         *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
	first=.false.
	endif
	goto 77

      elseif    ((case .eq. 'bq_tpq') .or. (case .eq. 'qg_tbq')) then
c---    single top t-channel, dynamic scale is:
c---        Q^2=-(pt-pb)^2 on the light quark line
c---        Q^2+mt^2 on the heavy quark line
        if (case .eq. 'bq_tpq') then
          q1scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(2,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(2,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(2,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(2,3))**2)
          b2scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(2,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(2,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(2,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(2,3))**2)+mt**2
          q2scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(1,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(1,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(1,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(1,3))**2)
          b1scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(1,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(1,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(1,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(1,3))**2)+mt**2
        else
          q1scale=-(
     .         +(p(5,4)+p(1,4))**2
     .         -(p(5,1)+p(1,1))**2
     .         -(p(5,2)+p(1,2))**2
     .         -(p(5,3)+p(1,3))**2)
          b2scale=-(
     .         +(p(5,4)+p(1,4))**2
     .         -(p(5,1)+p(1,1))**2
     .         -(p(5,2)+p(1,2))**2
     .         -(p(5,3)+p(1,3))**2)+mt**2
          q2scale=-(
     .         +(p(5,4)+p(2,4))**2
     .         -(p(5,1)+p(2,1))**2
     .         -(p(5,2)+p(2,2))**2
     .         -(p(5,3)+p(2,3))**2)
          b1scale=-(
     .         +(p(5,4)+p(2,4))**2
     .         -(p(5,1)+p(2,1))**2
     .         -(p(5,2)+p(2,2))**2
     .         -(p(5,3)+p(2,3))**2)+mt**2
	endif
	q1scale=max(dsqrt(q1scale),1d0) ! min. of 1 GeV for safety
	b2scale=max(dsqrt(b2scale),1d0)
	q2scale=max(dsqrt(q2scale),1d0)
	b1scale=max(dsqrt(b1scale),1d0)
	scale=q1scale   ! for safety
	facscale=scale
	if (first) then
        write(6,*)
        write(6,*)'************** Special scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*        muR = muF = DDIS (c.f. Z. Sullivan)       *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
	first=.false.
	endif
	goto 77
      endif

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
        if (part .ne. 'lord') then
	  write(6,*) 'Error: this scale choice is not suitable'
	  write(6,*) 'for a calculation beyond leading order.'
	  stop
	endif
        setscale=dsqrt(2d0*dot(p,1,2))
        msg='*              Scale = sqrt(s-hat)                 *'
      elseif (iset .eq. 4) then
c--- HT      
        if (part .ne. 'lord') then
	  write(6,*) 'Error: this scale choice is not suitable'
	  write(6,*) 'for a calculation beyond leading order.'
	  stop
	endif
        setscale=0d0
        do j=3,2*(1+n2+n3)+max(0,jets)
        setscale=setscale+getEt(ptildejet(0,j,4),ptildejet(0,j,1),
     .                    ptildejet(0,j,2),ptildejet(0,j,3))
        enddo
        if (first) msg='*      Dynamic scale = scalar Et sum (HT)'//
     . '          *'
      elseif (iset .eq. 5) then
c--- geometric mean of pt5, pt6, pt7      
        setscale=(pt(5,p)*pt(6,p)*pt(7,p))**(1d0/3d0)
        if (first) msg='*      Dynamic scale = (pt5*pt6*pt7)**(1/3)'//
     . '        *'
      elseif (iset .eq. 6) then
c--- pt7      
        setscale=pt(7,p)
        if (first) msg='*                 Dynamic scale = pt7      '//
     . '        *'
      elseif (iset .eq. 7) then
c--- HT-hat (HT, partonic, whether or not they pass cuts)      
        setscale=0d0
        do j=3,npart+2
        setscale=setscale+pt(j,p)
        enddo
        if (first) msg='*    Dynamic scale = parton Et sum (HT-hat)'//
     . '        *'
      elseif (iset .eq. 8) then
c--- sqrt(Qsq + 4*msq) for DIS processes
        setscale=(p(1,4)+p(5,4))**2
     .   -(p(1,1)+p(5,1))**2-(p(1,2)+p(5,2))**2-(p(1,3)+p(5,3))**2
        setscale=dsqrt(-setscale+4d0*(
     .   p(3,4)**2-p(3,1)**2-p(3,2)**2-p(3,3)**2))
        if (first) msg='*        Dynamic scale = sqrt(Qsq + 4*msq) '//
     . '        *'
      elseif ((iset .gt. 10) .and. (iset .lt. 18)) then
c--- sqrt(boson mass + boson pt) multiplied by a scale factor
        prefac=prescale(iset-10)
        if (nwz .eq. 0) then
          setscale=prefac*dsqrt(zmass**2+pttwo(3,4,p)**2)
          boson='z'
        else
          setscale=prefac*dsqrt(wmass**2+pttwo(3,4,p)**2)
          boson='w'
        endif  
        msg='* Dynamic scale = '//prechar(iset-10)//' * dsqrt('
     .    //boson//'mass**2'//' + pt_'//boson//'**2) *'
      elseif ((iset .gt. 20) .and. (iset .lt. 28)) then
c--- HT-hat (HT, partonic, whether or not they pass cuts) 
c---  multiplied by a scale factor     
        setscale=0d0
        do j=3,npart+2
        setscale=setscale+pt(j,p)
        enddo
        prefac=prescale(iset-20)
        setscale=prefac*setscale
        msg='*  Dynamic scale = '//prechar(iset-20)//
     .      ' * [parton Et sum (HT-hat)] *'
      elseif ((iset .gt. 30) .and. (iset .lt. 38)) then
c--- average jet pt      
        if (nqcdjets .eq. 0) then
          write(6,*) 'Invalid choice of scale - no jets!'
          stop
        endif
        setscale=aveptjet(p)
        prefac=prescale(iset-30)
        setscale=prefac*setscale
        msg='*        Dynamic scale = '//prechar(iset-30)//
     .      ' * [< pt_jet >]       *'

        msg='*          Dynamic scale = < pt_jet >              *'

       elseif ((iset .lt. 1) .or. (iset .gt. 8)) then
c--- catch invalid inputs
        write(6,*) 'Invalid dynamic scale!'
        stop
      endif
      
c--- choose a reasonable minimum value for the scale      
      if (setscale .lt. 1d0) setscale=1d0 ! for safety

      if (iscale .eq. 1) then
        scale=setscale
        rmsg=msg
      endif
      if ((iscale .eq. 2) .or. (iscalemax .eq. 1)) then
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
  
   77 continue
      
c--- catch absurdly large scales      
      if  (scale .gt. 3000d0) scale=3000d0
      if  (facscale .gt. 3000d0) facscale=3000d0

c--- run alpha_s
      as=alphas(scale,amz,nlooprun)
	
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2
      
      return
            
      end
