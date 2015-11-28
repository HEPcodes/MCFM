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
      integer itmx1,ncall1,itmx2,ncall2,icall
      double precision integ,integ_err
      logical dryrun
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      external vegas,vegas1,vegas2,vegas3
      
* basic variable initialization, print-out
      call mcfm_init

* This is the mcfm_vegas(warmup) call
* The Vegas parameters are those read from options.DAT for
* the warm-up stage (itmx1,ncall1) and binning should only take
* place if dryrun is set to true
      call mcfm_vegas(vegas,itmx1,ncall1,dryrun,integ,integ_err)
      
* This is the mcfm_vegas(accum) call
* This takes place only if dryrun is false
* The Vegas parameters are those read from options.DAT for
* the results stage (itmx2,ncall2) and binning takes place (.true.)
      if (dryrun .eqv. .false.) then
      call mcfm_vegas(vegas1,itmx2,ncall2,.true.,integ,integ_err)
      endif
      
* final processing and print-out
      call mcfm_exit(integ,integ_err)
      
      stop
      end
       
