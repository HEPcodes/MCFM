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
c     ------------------------------------------------------------------
      
      if (first) then
c ---   Open the file :
        outfile='mcfm.rz'
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

      return
      end

c

      subroutine bookfill(tag,p,wt)
      implicit none
      include 'constants.f'

      character tag*4
      double precision p(mxpart,4)
      double precision wt 

      if (tag.eq.'book') then
        call dswntuplebook
      elseif (tag .eq. 'plot') then
        call dswntuplefill(p,wt)
      endif

      return
      end

c
      subroutine dswntuplebook
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR
      integer IQUEST
      common/QUEST/IQUEST(100)

      character*3 CHTAGS(17)
      data CHTAGS/
     +     'px3','py3','pz3','E3',
     +     'px4','py4','pz4','E4',
     +     'px5','py5','pz5','E5',
     +     'px6','py6','pz6','E6',
     +     'wt' /

      integer ISTAT
      character*100 outfile

c     ------------------------------------------------------------------

c --- Create the output file :
      outfile='mcfm.rz'
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
      call hbookn(300,"MCFM",17,'//HISTOS',4096,CHTAGS)

      return
      end

c
      subroutine dswntuplefill(p,wt)
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      include 'constants.f'

      double precision p(mxpart,4)
      double precision wt 

      integer NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer i
      real pfill(17)

c     ------------------------------------------------------------------

c --- Fill the ntuple :
      do i=1,4
        pfill(4*(i-1)+1)=sngl(p(i+2,1))
        pfill(4*(i-1)+2)=sngl(p(i+2,2))
        pfill(4*(i-1)+3)=sngl(p(i+2,3))
        pfill(4*(i-1)+4)=sngl(p(i+2,4))
      enddo
      pfill(17)=sngl(wt)

      call hfn(300,pfill)

      return
      end

c

