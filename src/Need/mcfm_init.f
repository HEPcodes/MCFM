      block data newinput_data
      implicit none
      logical newinput
      common/newinput/newinput
      data newinput/.true./
      end

      subroutine mcfm_init
************************************************************************
*                                                                      *
*  This routine should initialize any necessary variables and          *
*  perform the usual print-outs                                        *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'cutoff.f'
      include 'efficiency.f'
      include 'limits.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'scale.f'
      include 'verbose.f'
      include 'workdir.f'
      logical creatent,dswhisto
      double precision rtsmin,sqrts,p1ext(4),p2ext(4),
     . vector(mxdim),p(mxpart,4),s(mxpart,mxpart),val
      integer j,k
      logical newinput
      common/newinput/newinput
      common/rtsmin/rtsmin
      common/energy/sqrts
      common/pext/p1ext,p2ext
      data p/mxpart*3d0,mxpart*4d0,mxpart*0d0,mxpart*5d0/
      data vector/mxdim*0d0/

* Welcome banner
      call banner
* Read-in the options.DAT file
      if (newinput) then
        call reader_input
      else
        call reader
        workdir=''
      endif

* Initialize efficiency variables      
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

* Set-up incoming beams and PS integration cut-offs
      rtsmin=min(rtsmin,dsqrt(wsqmin+cutoff))
      rtsmin=min(rtsmin,dsqrt(bbsqmin+cutoff))
      taumin=(rtsmin/sqrts)**2
      xmin=taumin

      p1ext(4)=-half*sqrts
      p1ext(1)=0d0
      p1ext(2)=0d0
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0d0
      p2ext(2)=0d0
      p2ext(3)=+half*sqrts

* Initialize all histograms
* npart=6 is a dummy value, to ensure that all histograms are included
      npart=6
      call dotem(8,p,s)
      val=1d-15   
      call nplotter(vector,s,p,val,1)
       
      do j=1,mxpart
      do k=1,4
      p(j,k)=0d0
      enddo
      enddo 

* Set-up run name
      call setrunname(scale)
           
      return
      end
            
