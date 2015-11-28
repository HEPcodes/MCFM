      subroutine mcfm_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)
************************************************************************
*                                                                      *
*  This routine should perform the sweeps of vegasnr                     *
*                                                                      *
*    Input parameters:                                                 *
*       myinit  :  the vegasnr routine entry point                       *
*       myitmx  :  the number of vegasnr sweeps                          *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error
*                                                                      *
************************************************************************
      implicit none
      include 'gridinfo.f'
      include 'realwt.f'
      include 'vegas_common.f'
      integer myitmx,myncall,myinit,i
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,chir,
     . xreal,xreal2,xinteg,xerr,adjust
      character*4 part,mypart
      common/part/part
      common/bin/bin
      common/xreal/xreal,xreal2
      double precision lowint,virtint,realint
      double precision region(2*mxdim)
      logical first
      external lowint,virtint,realint
      data first/.true./
      save first
           
c--- Initialize all integration results to zero, so that the
c--- total of virt and real may be combined at the end for 'tota'
      sig=0d0
      sigr=0d0
      sd=0d0
      sdr=0d0
      xreal=0d0
      xreal2=0d0
      
c--- Put the vegasnr parameters in the common block
      itmx=myitmx
      ncall=myncall
      bin=mybin
      
c--- Basic lowest-order integration
      if (part .eq. 'lord') then
       call boundregion(ndim,region)
       call vegasnr(region,ndim,lowint,myinit,myncall,myitmx,
     .               0,sig,sd,chi)
      endif

c--- Store value of part in mypart, which will be retained;
c--- part can be changed to make sure that the tota option works.
      mypart=part
      
c--- If we're doing the tota integration, then set up the grid info
      if (mypart .eq. 'tota') then
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_tota_virt'          
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_tota_virt'
        endif
      endif        
      
c--- Virtual integration should have one extra dimension
c--- (added and then taken away)
      if ((mypart .eq. 'virt') .or. (mypart .eq. 'tota'))  then
        part='virt'
        call boundregion(ndim+1,region)
        call vegasnr(region,ndim+1,virtint,myinit,myncall,myitmx,
     .              0,sig,sd,chi)
      endif
            
c--- If we're doing the tota integration, then set up the grid info
      if (mypart .eq. 'tota') then
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_tota_real'          
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_tota_real'
        endif
      endif        
      
c--- Real integration should have three extra dimensions
c--- 'realwt' is a special option that in general should be false
c--- ('realwt' true samples the integral according to the
c---   unsubtracted real emission weight)
      if (mypart .eq. 'real') then
        part='real'
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        call boundregion(ndim+3,region)
        call vegasnr(region,ndim+3,realint,myinit,myncall,myitmx,
     .              0,sigr,sdr,chi)
        write(6,*) 
        ncall=myncall
        if (realwt) then
          sigr=xreal
          sdr=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif
      if (mypart .eq. 'tota')  then
        part='real'
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        adjust=(dfloat(ndim+3))/(dfloat(ndim+1))
        ncall=int(dfloat(myncall)**adjust)
        write(6,*) 'Adjusting number of points for real to',ncall
        call boundregion(ndim+3,region)
        call vegasnr(region,ndim+3,realint,myinit,ncall,myitmx,
     .              0,sigr,sdr,chi)
        write(6,*) 
        ncall=myncall

        if (realwt) then
          sigr=xreal
          sdr=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif      

c--- calculate integration variables to be returned
      xinteg=sig+sigr
      xerr=dsqrt(sd**2+sdr**2)      
      
c--- return part to its real value
      part=mypart
      first=.false.
      
      return
      end
      
      
      subroutine boundregion(idim,region)
c--- Initializes integration region [0,1] for each variable
c--- in the idim-dimensional integration range
      implicit none
      include 'mxdim.f'
      integer i,idim
      double precision region(2*mxdim)
      
      do i=1,idim
      region(i)=0d0
      region(i+idim)=1d0
      enddo
      
      return
      end
      
      
