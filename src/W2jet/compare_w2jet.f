      subroutine compare(p)
      implicit none
      
      include 'constants.f'
      include 'zerowidth.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'flags.f'
      double precision msqx(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqxa(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      common/msqx/msqx
    
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf),
     . msqa(-nf:nf,-nf:nf)
      integer i,j,k,l,m
      double precision aemmz,ans,Vfac,ee,xmsqLL,xmsqLR,xmsq1LR,xmsqud,
     . xmsquu,xmsqRL,
     . P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                 
      integer nwz,n1,n2
      common/nwz/nwz

      Qflag=.true.
      Gflag=.false.
      write(6,*) 'Gflag',Gflag
      write(6,*) 'Qflag',Qflag

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
        p5(j)=p(5,i)
        p6(j)=p(6,i)
        p7(j)=p(7,i)
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

c      call qqb_w2jet_g_mad(p,msqa)

c      call qqb_w2jet_g(p,msq)     

c      call testem(1,5,2,6,7,3,4,xmsqLR,xmsq1LR,xmsqRL,xmsqLL,
c     . xmsqud,xmsquu)
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uL dR -> dL dR',0.5d0*xmsqLR
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uL dR -> dR dL',0.5d0*xmsq1LR
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uL sL -> dL sL',xmsqLL
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uL sR -> dL sR',xmsqLR
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uL dL -> dL dL',0.5d0*xmsqud
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uL uL -> dL uL',xmsquu
c      write(*,*) 'RKE: 5 1 6 2 7 4 3 : uR uL -> dL uR',xmsqRL
c      write(*,*)
      
      call qqb_wp2jet_mad2(p,msqa)
      do j=-nf,nf
      do k=-nf,nf
      do l=-nf,nf
      do m=-nf,nf
        msqxa(j,k,l,m)=msqx(j,k,l,m)
      enddo
      enddo
      enddo
      enddo

      call qqb_w2jetx(p,msq)
      
      write(*,*) '        MCFM             MADG               ratio'

      
      do j=-nf,nf
      do k=-nf,nf
      do l=-nf,nf
      do m=-nf,nf
c      if (abs(msqx(j,k,l,m)/msqxa(j,k,l,m)-1d0) .gt.1d-6)
      if ((msqx(j,k,l,m).gt.0d0) .or. (msqxa(j,k,l,m).gt.0d0) )
     . write(6,98) j,k,l,m,msqx(j,k,l,m),msqxa(j,k,l,m),
     .  msqx(j,k,l,m)/msqxa(j,k,l,m)
      enddo
      enddo
      enddo
      enddo

      pause


c      do j=-nf,nf
c      do k=-nf,nf
      do j=-5,5
      do k=-5,5
c      if ((j.lt.0 .and. k.lt.0) .or. (j.gt.0 .and. k.gt.0)) then
c      if (abs(msq(j,k)/msqa(j,k)-1d0) .gt.1d-6)
      if ((msq(j,k).gt.0d0) .or. (msqa(j,k).gt.0d0) )
     .write(6,99) j,k,msq(j,k),msqa(j,k),msq(j,k)/msqa(j,k)
c      endif
   20 enddo
      enddo
      pause

c      call qqb_z2jet_g(p,msq)
c      call initialize
c      call qqb_z2jet_g_mad(p,msqa)

c      write(*,*) '        MCFM new        MADG            ratio'
c      do j=-nf,nf
c      do k=-nf,nf
c      if (abs(msq(j,k)/msqa(j,k)-1d0) .gt. 1d-6)
c     . write(6,*) j,k,msq(j,k),msqa(j,k),msq(j,k)/msqa(j,k)
c      enddo
c      enddo

c      pause


   98 format(4i3,2e17.9,f15.9)
   99 format(2i3,2e17.9,f15.9)

      
      return
      end
