c --- Dave Waters, 28.05.2001
c --- =======================
c --- Provide routines that are called from mcfm and fill
c --- histograms in an hbook file. It is more convenient 
c --- and flexible to obtain hbook histograms
c --- in this way rather than parse the standard mcfm
c --- output files.

      subroutine dswhbook(n,titlex,dx,xmin,xmax)
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      integer n
      character titlex*8
      real*8 dx,xmin,xmax
      logical first
      data first /.true./
      save first
      
      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer NHISTOMAX
      parameter(NHISTOMAX=200)
      integer BOOKED(NHISTOMAX)
      common/BOOKPATTERN/BOOKED

      integer ISTAT,i
      character*100 outfile

      character*72 runname
      integer nlength, lenocc
      common/runname/runname,nlength
c     ------------------------------------------------------------------
      
      if (first) then
c ---   Open the file :
	outfile=runname(1:lenocc(runname))//'.rz'
        call hlimit(NWPAWC)
c ---   A record length of 8192 words will produce a warning but is 
c ---   OK for most systems (see HBOOK manual, p.21) : 
        call hropen(30,'HISTOS',outfile,'N',8192,ISTAT)
        if (ISTAT.NE.0) then
          write(6,*) 'ERROR ReadInOut : hropen (HISTOS), ISTAT =',
     +                                                   ISTAT
          stop
        endif
        do i=1,NHISTOMAX
          BOOKED(i)=0
        enddo
        first = .false.
      endif

c --- Book the histogram :
      call hcdir('//HISTOS',' ')
      call hbook1(n,titlex,int(sngl((xmax-xmin)/dx)),
     +            sngl(xmin),sngl(xmax),0.)
      BOOKED(n)=1

      return
      end
c

      subroutine dswhfill(n,var,wgt)
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      integer n
      real*8 var,wgt

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR
c     ------------------------------------------------------------------

c --- Fill the histogram :
      call hcdir('//HISTOS',' ')
      call hf1(n,sngl(var),sngl(wgt))

      return
      end

c

      subroutine dswhrout
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer NHISTOMAX
      parameter(NHISTOMAX=200)
      integer BOOKED(NHISTOMAX)
      common/BOOKPATTERN/BOOKED

c --- Common block to control output
c --- (histograms or ntuple)
      logical creatent
      common /outputflags/creatent

      integer ICYCLE,i
c     ------------------------------------------------------------------

      if (.NOT.creatent) then
c ---   Read out the booked histograms :
        do i=1,NHISTOMAX
          if (BOOKED(i).EQ.1) then
            call hcdir('//HISTOS',' ')
c           write(6,*) 'Calling hrout'
            call hrout(i,ICYCLE,' ')
            if (ICYCLE.NE.1) then
              write(6,*) 'ERROR dswhrout : ICYCLE =',ICYCLE
              write(6,*) 'ERROR dswhrout : Dupicated HISTO ID'
            endif
          endif
        enddo

      else

c ---   Read out the ntuple :
c ---   =====================
c       write(6,*) 'Calling hrout'
c ---   For some reason, this call to hrout produces an ICYCLE of 2.
c ---   Possibly because the ntuple is explicitly a disk ntuple.
        call hrout(300,ICYCLE,' ')
c       if (ICYCLE.NE.1) then
c         write(6,*) 'ERROR dswhrout : ICYCLE =',ICYCLE
c         write(6,*) 'ERROR dswhrout : Dupicated NTUPLE ID'
c       endif

      endif

      return
      end

c

      subroutine dswclose
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer ICYCLE
c     ------------------------------------------------------------------

      call hrend('HISTOS')
      close(30)
      
      write(6,*) '<----- Completed a batch of n-tuples ----->' 
      call flush(6)
      
      return
      end

c

      subroutine bookfill(tag,p,wt)
      implicit none
      include 'constants.f'
      include 'maxwt.f'

      character tag*4
      double precision p(mxpart,4)
      double precision wt 

      if (.not.skipnt) then
        if (tag.eq.'book') then
          call dswntuplebook
        elseif (tag .eq. 'plot') then
          call dswntuplefill(p,wt)
        endif
      endif

      return
      end

c
      subroutine dswntuplebook
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
c--- Included for the extra code below
      include 'npart.f'
      include 'mxdim.f'
      include 'scale.f'
      include 'facscale.f'

      double precision scale_store,facscale_store

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR
      integer IQUEST
      common/QUEST/IQUEST(100)

c--- Added to keep track of number of momenta entries to be filled
      character*4 part
      common/part/part
c--- Extra definitions to facilitate dummy call to lowint
      double precision dummy,wgt,r(mxdim),lowint
      integer ifill
      integer imaxmom
      common/imaxmom/imaxmom
            
      character*3 CHTAGS4(13)
      character*3 CHTAGS5(17)
      character*3 CHTAGS6(21)
      character*3 CHTAGS7(25)
      character*3 CHTAGS8(29)
      data CHTAGS4/
     +     'px3','py3','pz3','E3 ',
     +     'px4','py4','pz4','E4 ',
     +     'wt ','gg ','gq ','qq ','qqb'/
      data CHTAGS5/
     +     'px3','py3','pz3','E3 ',
     +     'px4','py4','pz4','E4 ',
     +     'px5','py5','pz5','E5 ',
     +     'wt ','gg ','gq ','qq ','qqb'/
      data CHTAGS6/
     +     'px3','py3','pz3','E3 ',
     +     'px4','py4','pz4','E4 ',
     +     'px5','py5','pz5','E5 ',
     +     'px6','py6','pz6','E6 ',
     +     'wt ','gg ','gq ','qq ','qqb'/
      data CHTAGS7/
     +     'px3','py3','pz3','E3 ',
     +     'px4','py4','pz4','E4 ',
     +     'px5','py5','pz5','E5 ',
     +     'px6','py6','pz6','E6 ',
     +     'px7','py7','pz7','E7 ',
     +     'wt ','gg ','gq ','qq ','qqb'/
      data CHTAGS8/
     +     'px3','py3','pz3','E3 ',
     +     'px4','py4','pz4','E4 ',
     +     'px5','py5','pz5','E5 ',
     +     'px6','py6','pz6','E6 ',
     +     'px7','py7','pz7','E7 ',
     +     'px8','py8','pz8','E8 ',
     +     'wt ','gg ','gq ','qq ','qqb'/

      integer ISTAT
      character*100 outfile

      character*72 runname
      integer nlength, lenocc
      common/runname/runname,nlength

      logical first
      integer batchno
      character*3 batchstr,getstr
      data first/.true./
      save first,batchno
c     ------------------------------------------------------------------

c--- Need to ascertain the correct size for momenta n-tuples when this routine
c--- is called for the first time, achieved via a dummy call to lowint
      if (first) then      
        do ifill=1,mxdim
          r(ifill)=0.5d0
        enddo
c--- Be careful that dynamic scale choices aren't ruined
c--- (in versions 5.1 and before, this occured when calling lowint)
        scale_store=scale
	facscale_store=facscale
        dummy=lowint(r,wgt)
        scale=scale_store
	facscale=scale_store
	
        imaxmom=npart
        if ((part.eq.'real').or.(part.eq.'tota').or.(part.eq.'todk'))
     .    imaxmom=imaxmom+1
	  batchno=-1
	  first=.false.
      endif
      
c---- Increment the batch number counter and convert to a string    
      batchno=batchno+1
      batchstr=getstr(batchno)
      
c --- Create the output file :
      outfile=runname(1:lenocc(runname))//batchstr//'.rz'
      write(6,*) '<----- Creating batch ',batchno,' of n-tuples ----->' 
      call flush(6)

      call hlimit(NWPAWC)
      IQUEST(10) = 65000
c --- A record length of 8192 words will produce a warning but is 
c --- OK for most systems (see HBOOK manual, p.21) : 
      call hropen(30,'HISTOS',outfile,'N',8192,ISTAT)
      if (ISTAT.NE.0) then
        write(6,*) 'ERROR ReadInOut : hropen (HISTOS), ISTAT =',
     +        ISTAT
      endif

c --- Book an extremely simple row-wise ntuple. Make it explicitly
c --- a disk resident ntuple by specifying the top directory name
c --- of the previously opened RZ file in the 4th argument 
c --- (see HBOOK manual, p.19) :
      if     (imaxmom .eq. 2) then
        call hbookn(300,'MCFM',13,'//HISTOS',4096,CHTAGS4)
      elseif (imaxmom .eq. 3) then
        call hbookn(300,'MCFM',17,'//HISTOS',4096,CHTAGS5)
      elseif (imaxmom .eq. 4) then
        call hbookn(300,'MCFM',21,'//HISTOS',4096,CHTAGS6)
      elseif (imaxmom .eq. 5) then
        call hbookn(300,'MCFM',25,'//HISTOS',4096,CHTAGS7)
      elseif (imaxmom .eq. 6) then
        call hbookn(300,'MCFM',29,'//HISTOS',4096,CHTAGS8)
      else
        write(6,*) 'Problem in dswntuplebook - value npart=',npart
        write(6,*) 'not anticipated. Program halted.'
        stop
      endif
      
      return
      end

c
      subroutine dswntuplefill(p,wt)
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      include 'constants.f'
      include 'wts_bypart.f'
      
      double precision p(mxpart,4)
      double precision wt 

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

c--- Extra common block to carry the information about maximum momenta entries
      integer imaxmom
      common/imaxmom/imaxmom

      integer i
      real pfill(imaxmom*4+5)

c--- Variables to count the number of ntuples filled so far and set the
c--- maximum number filled per batch; the current value is somewhat arbitrary
c--- and may be modified by the user
      integer icount,batchlimit
      data icount/0/
      save icount
      
      parameter(batchlimit=1000000)

c     ------------------------------------------------------------------

c--- if the weight is zero, don't bother to add the n-tuple
      if (wt .eq. 0d0) then
        return
      endif

c--- The n-tuple should be filled, but first check if the current batch
c--- is already full; if so, close and then open a new one
      icount=icount+1
      if (mod(icount,batchlimit) .eq. 0) then
        call dswhrout
	call dswclose
	call dswntuplebook
      endif
      
c --- Fill the ntuple :
      do i=1,imaxmom
        pfill(4*(i-1)+1)=sngl(p(i+2,1))
        pfill(4*(i-1)+2)=sngl(p(i+2,2))
        pfill(4*(i-1)+3)=sngl(p(i+2,3))
        pfill(4*(i-1)+4)=sngl(p(i+2,4))
      enddo
      pfill(imaxmom*4+1)=sngl(wt)
      pfill(imaxmom*4+2)=sngl(wt_gg)
      pfill(imaxmom*4+3)=sngl(wt_gq)
      pfill(imaxmom*4+4)=sngl(wt_qq)
      pfill(imaxmom*4+5)=sngl(wt_qqb)
            
c     write(6,*) 'Filling ntuple with weight',wt

      call hfn(300,pfill)
      
      return
      end

c
