      subroutine mcfm_exit(xinteg,xinteg_err)
************************************************************************
*                                                                      *
*  This routine should perform the final processing and print-outs     *
*                                                                      *
************************************************************************
      implicit none
      include 'efficiency.f'
      integer j,k
      double precision xinteg,xinteg_err
      double precision lord_bypart(-1:1,-1:1),lordnorm
      character*4 part
      logical creatent,dswhisto
      common/outputflags/creatent,dswhisto
      common/part/part
      common/bypart/lord_bypart
      
c--- Print-out the value of the integral and its error
      write(6,*) 
      write(6,*)'Value of final ',part,' integral is',
     . xinteg,' +/-',xinteg_err, ' fb'
     
c--- Print-out a summary of the effects of jets and cuts
      write(6,*) 
      write(6,*) 'Total number of shots       : ',ntotshot
      write(6,*) 'Total no. failing cuts      : ',ntotzero
      write(6,*) 'Number failing jet cuts     : ',njetzero
      write(6,*) 'Number failing process cuts : ',ncutzero
      write(6,*) 
      call flush(6)

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
      
      if (part .eq. 'lord') then
        lordnorm=0d0
        do j=-1,1
        do k=-1,1
          lordnorm=lordnorm+lord_bypart(j,k)
        enddo
        enddo
        write(6,*) 'Contribution from parton sub-processes:'
        write(6,*) '---------------------------------------'
        write(6,55) '   GG    ',
     .  lord_bypart( 0, 0)/lordnorm*xinteg,
     .  lord_bypart( 0, 0)/lordnorm*100d0
        write(6,55) 'GQ + GQB ',
     .  (lord_bypart( 0,+1)+lord_bypart( 0,-1))/lordnorm*xinteg,
     .  (lord_bypart( 0,+1)+lord_bypart( 0,-1))/lordnorm*100d0
        write(6,55) 'QG + QBG ',
     .  (lord_bypart(+1, 0)+lord_bypart(-1, 0))/lordnorm*xinteg,
     .  (lord_bypart(+1, 0)+lord_bypart(-1, 0))/lordnorm*100d0
        write(6,55) 'QQ + QBQB',
     .  (lord_bypart(+1,+1)+lord_bypart(-1,-1))/lordnorm*xinteg,
     .  (lord_bypart(+1,+1)+lord_bypart(-1,-1))/lordnorm*100d0
        write(6,55) '   QQB   ',
     .  (lord_bypart(+1,-1)+lord_bypart(-1,+1))/lordnorm*xinteg,
     .  (lord_bypart(+1,-1)+lord_bypart(-1,+1))/lordnorm*100d0
        write(6,*) '---------------------------------------'
      endif
      call flush(6)

   54 format(a20,f6.2,'%')
   55 format(4x,a9,'  |',f14.5,f8.2,'%')

c--- Finalize the histograms, if we're not filling ntuples instead
      if (creatent .eqv. .false.) then
        if (dswhisto .eqv. .false.) then
c--- Traditional MCFM histograms
          call histofin(xinteg,xinteg_err)
        else
c--- DSW histograms - store the information
c          call dswhbook(200,'Sigma',1.0d0,0.0d0,10.0d0)
c          call dswhfill(200,0.5d0,xinteg)
c          call dswhfill(200,1.5d0,xinteg_err)
c--- DSW histograms - output and close file
          call dswhrout
          call dswclose
        endif
      else
c        call dswhrout
c        call dswclose
      endif
      
      return
      
      end
