      subroutine compare_2jet_mad(p)
      implicit none
      
      include 'constants.f'
      include 'zerowidth.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'flags.f'
      double precision msqx(-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     . p(mxpart,4),msq(-nf:nf,-nf:nf),msqc(-nf:nf,-nf:nf)
      integer i,j,k,l,m
      double precision aemmz,ans,Vfac,ee,xmsqLL,xmsqLR,xmsq1LR,xmsqud,
     . xmsquu,xmsqRL,
     . P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                 
      integer nwz,n1,n2
      double precision mqq(0:2,fn:nf,fn:nf)
      common/nwz/nwz
      common/msqx/msqx

c--- implement the momentum exchange      
      do i=1,4
        if (i.lt.4) then
          j=i
        else
          j=0
        endif 
        p1(j)=-p(1,i)
        p2(j)=-p(2,i)
        p3(j)=p(3,i)
        p4(j)=p(4,i)
      enddo

c--- make call to MCFM routines
      aemmz=1d0/128d0      
      ee=sqrt(fourpi*aemmz)
      write(6,*) 'ee',ee
      gwsq=fourpi*aemmz/xw
      gw=sqrt(gwsq)
      gsq=1d0
      call ckmfill(nwz)
      
      write(*,*) 'zmass, wmass',zmass,wmass
      write(*,*) 'new xw,gw',xw,gw

      write(*,*)

      call initialize

      call qqb_2jet(p,msq)
      call qqb_2jet_mad(p,msqc)


      write(*,*) '        MCFM             MADG               ratio'

      
c      do j=-nf,nf
c      do k=-nf,nf
c      write(6,98) j,k,msq(j,k),msqc(j,k)
c      enddo
c      enddo

c      pause


c      do j=-nf,nf
c      do k=-nf,nf
      do j=-2,2
      do k=-2,2
c      if ((j.lt.0 .and. k.lt.0) .or. (j.gt.0 .and. k.gt.0)) then
      if (abs(msq(j,k)/msqc(j,k)-1d0) .gt.1d-6)
c      if ((msq(j,k).gt.0d0) .or. (msqc(j,k).gt.0d0) )
     .write(6,99) j,k,msq(j,k),msqc(j,k),msq(j,k)/msqc(j,k)
c      endif
   20 enddo
      enddo
      pause


   98 format(2i3,g17.6,g17.6)
   99 format(2i3,2e17.9,f15.9)

      
      return
      end
