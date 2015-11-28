      subroutine mcfm_vegas(vegas_ent,myitmx,myncall,mybin,xinteg,xerr)
************************************************************************
*                                                                      *
*  This routine should perform the sweeps of Vegas                     *
*                                                                      *
*    Input parameters:                                                 *
*    vegas_ent  :  the Vegas routine (entry point)                     *
*       myitmx  :  the number of Vegas sweeps                          *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error
*                                                                      *
************************************************************************
      implicit none
      include 'realwt.f'
      include 'vegas_common.f'
      integer myitmx,myncall
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,chir,
     . xreal,xreal2,xinteg,xerr
      character*4 part
      common/part/part
      common/bin/bin
      common/xreal/xreal,xreal2
      external vegas_ent,lowint,virtint,realint

c--- Initialize all integration results to zero, so that the
c--- total of virt and real may be combined at the end for 'tota'
      sig=0d0
      sigr=0d0
      sd=0d0
      sdr=0d0
      xreal=0d0
      xreal2=0d0
      
c--- Put the Vegas parameters in the common block
      itmx=myitmx
      ncall=myncall
      bin=mybin
      
c--- Basic lowest-order integration
      if (part. eq. 'lord') then
        call vegas_ent(lowint,sig,sd,chi)
      endif
       
c--- Virtual integration should have one extra dimension
c--- (added and then taken away)
      if ((part .eq. 'virt') .or. (part .eq. 'tota'))  then
        ndim=ndim+1
        call vegas_ent(virtint,sig,sd,chi) 
        ndim=ndim-1
      endif
      
c--- Real integration should have three extra dimensions
c--- 'realwt' is a special option that in general should be false
c--- ('realwt' true samples the integral according to the
c---   unsubtracted real emission weight)
      if ((part .eq. 'real') .or. (part .eq. 'tota'))  then
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        ndim=ndim+3
        call vegas_ent(realint,sigr,sdr,chir)
        ndim=ndim-3
c--- If 'realwt' is true, xreal and xreal2 now hold the true results
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
      
      return
      end
      
