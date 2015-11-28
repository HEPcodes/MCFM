      program mcfm
************************************************************************
*                                                                      *
*  This is the main program for MCFM                                   *
*                                                                      *
*  The sequence of calls should always be:                             *
*   call mcfm_init          : basic variable initialization, print-out *
*   call mcfm_vegas(warmup) : warm-up the Vegas grid                   *
*   call mcfm_vegas(accum)  : accumulate results                       *
*   call mcfm_exit          : final processing and print-out           *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'maxwt.f'
      include 'eventbuffer.f'
      integer itmx1,ncall1,itmx2,ncall2,icall
      double precision integ,integ_err
      logical dryrun
      integer i,pflav,pbarflav
      double precision p(mxpart,4),wt
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      external vegas,vegas1,vegas2,vegas3
      

* basic variable initialization, print-out
      call mcfm_init

* in initial phases, we don't want any unweighting to take place.
* this will be set true in the first call to getevent.
      unweight = .false.

* This is the mcfm_vegas(warmup) call
* The Vegas parameters are those read from options.DAT for
* the warm-up stage (itmx1,ncall1) and binning should only take
* place if dryrun is set to true
      call mcfm_vegas(vegas,itmx1,ncall1,dryrun,integ,integ_err)
      
* This is the mcfm_vegas(accum) call
* This takes place only if dryrun is false
* The Vegas parameters are those read from options.DAT for
* the results stage (itmx2,ncall2) and binning takes place (.true.)
* wtmax may have been set during the dry run, so re-set here :
      wtmax = 0d0
      if (dryrun .eqv. .false.) then
      call mcfm_vegas(vegas1,itmx2,ncall2,.true.,integ,integ_err)
      endif
      
* So far we have not used VEGAS to generate any events.
* Make sure future calls to "getevent" are aware of this :
      numstored = 0

      if (evtgen) then
        write(6,*) 'Generate events :'
        do i=1,500
          call mcfm_getevent(p,wt,pflav,pbarflav)
          call fill_stdhep(p,0,0,wt)
c         call write_stdhep(6)
        enddo
      endif

* final processing and print-out
      call mcfm_exit(integ,integ_err)
      
      stop
      end
       
