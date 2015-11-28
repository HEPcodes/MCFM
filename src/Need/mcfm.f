
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'npart.f'
      include 'vegas_common.f'
      include 'realwt.f'
      include 'efficiency.f'

      double precision sig,sd,chi,sigr,sdr,chir
      double precision sqrts,xerr
      double precision virtint,realint,lowint,totint
      double precision p1ext(4),p2ext(4)
      double precision xreal,xreal2
      double precision vector(mxdim),p(mxpart,4),s(mxpart,mxpart),val
      integer itmx1,ncall1,itmx2,ncall2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      character*4 part
      common/part/part

      logical verbose,bin,dryrun
      double precision rtsmin,taumin
      common/verbose/verbose
      common/energy/sqrts
      common/bin/bin
      common/dryrun/dryrun
      common/taumin/taumin
      common/rtsmin/rtsmin
      common/xreal/xreal,xreal2
      common/pext/p1ext,p2ext
c      common/pchoice/jp,kp
      external virtint,realint,lowint,totint
      data p/mxpart*0d0,mxpart*0d0,mxpart*0d0,mxpart*0d0/

      call banner
      call reader
      
      njetzero=0
      ncutzero=0
      ntotzero=0
      ntotshot=0
      
      if (verbose) then
      write(6,*)
      write(6,*) '****************************************'
      write(6,*) '*     Cross section in femtobarns      *'
      write(6,*) '****************************************'
      write(6,*)
      endif
           
      taumin=(rtsmin/sqrts)**2

      p1ext(4)=-half*sqrts
      p1ext(1)=0d0
      p1ext(2)=0d0
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0d0
      p2ext(2)=0d0
      p2ext(3)=+half*sqrts

C---call nplotter to initialize histograms
c---we set npart as a dummy, so that all histograms are initialized
      npart=5      
      call nplotter(vector,s,p,val,1)

C---initialize all final results to zero
      sig=0d0
      sigr=0d0
      sd=0d0
      sdr=0d0
      xreal=0d0
      xerr=0d0

      if (part. eq. 'lord') then
      itmx=itmx1
      ncall=ncall1
      bin=dryrun
      call vegas(lowint,sig,sd,chi) 
      if (dryrun) goto 40
      bin=.true.
      itmx=itmx2
      ncall=ncall2
      call vegas1(lowint,sig,sd,chi) 
      go to 40

C-----
      elseif ((part .eq. 'virt') .or. (part .eq. 'tota'))  then
      ndim=ndim+1

      itmx=itmx1
      ncall=ncall1
      bin=dryrun
      call vegas(virtint,sig,sd,chi) 

      if (dryrun) then
      if (part .eq. 'virt') then 
      goto 40
      elseif (part .eq. 'tota') then
      ndim=ndim-1
      goto 30
      endif
      endif

      bin=.true.
      itmx=itmx2
      ncall=ncall2
      call vegas1(virtint,sig,sd,chi) 
      ndim=ndim-1
      endif

      if (part .eq. 'virt') goto 40
            
 30   continue
      ndim=ndim+3
      itmx=itmx1
      ncall=ncall1
      bin=dryrun

      if (realwt) then
      nprn=0
      endif
      xreal=0d0
      xreal2=0d0

      call vegas(realint,sigr,sdr,chir) 

      if (realwt) then
      sdr=sqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
      write(6,*) itmx,' iterations of ',ncall,' calls'
      write(6,*) 'Value of subtracted integral',xreal
      write(6,*) 'Error on subtracted integral',sdr
      endif
      if (dryrun) goto 40


      bin=.true.
      itmx=itmx2
      ncall=ncall2

      if (realwt) then
      nprn=0
      endif
      xreal=0d0
      xreal2=0d0

      call vegas1(realint,sigr,sdr,chir) 
      if (realwt) then
      sdr=sqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
      write(6,*) itmx,' iterations of ',ncall,' calls'
      write(6,*) 'Value of subtracted integral',xreal
      write(6,*) 'Error on subtracted integral',sdr
      endif

 40   call histofin

      if (realwt .eqv. .false.)  xreal=sigr

      xerr=sqrt(sdr**2+sd**2)
      write(6,*) 
      write(6,*)'Value of final ',part,' integral is',
     & sig+xreal,'+/-',xerr
     
      write(6,*) 
      write(6,*) 'Total number of shots       : ',ntotshot
      write(6,*) 'Total no. failing cuts      : ',ntotzero
      write(6,*) 'Number failing jet cuts     : ',njetzero
      write(6,*) 'Number failing process cuts : ',ncutzero
      write(6,*) 

c--- Calculate the actual number of shots that were passed
c--- through the jet and cut routines
      ntotshot=ntotshot-(ntotzero-njetzero-ncutzero)
      write(6,54) 'Jet efficiency : ',
     .  100d0-100d0*dfloat(njetzero)/dfloat(ntotshot)
      write(6,54) 'Cut efficiency : ',
     .  100d0-100d0*dfloat(ncutzero)/dfloat((ntotshot-njetzero))
      write(6,54) 'Total efficiency : ',
     .  100d0-100d0*dfloat((njetzero+ncutzero))/dfloat(ntotshot)
      write(6,*) 
      
   54 format(a20,f6.2,'%')
    
      stop
      end




