      subroutine qqb_w2jet_g_old(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W^+ + g(p7)
c                           |     |
c                           |     --> nu(p3)+e^+(p4)
c                           |
c                           ---> f(p5)+f(p6)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'prods.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'hardscale.f'
      include 'flags.f'
      include 'lc.f'
      integer j,k,nu,hq,Qh,hg,lh,n1,n2,f3,f4,jj,kk
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision mmsq_gg,mmsq_qqb,mmsq_qbq,mmsq_qg
      double precision mmsq_gq,mmsq_gqb,mmsq_qbg
      double precision fac,LRq(2),LRb(2),lr1(2),Vfac,stat
      double precision propsq,ofac,
     .                 msq_qqbs(5,5,5,5),msq_qbqs(5,5,5,5),
     .                 msq_qqb(5,5,5,5),msq_qbq(5,5,5,5),
     .                 msq_qg(5,5,5,5),msq_qbg(5,5,5,5),
     .                 msq_gqb(5,5,5,5),msq_gq(5,5,5,5),
     .                 msq_qq(5,5,5,5),msq_qbqb(5,5,5,5)
      double complex czq,czb
      logical first
      data first/.true./
      save first

      if (first) then
      first=.false.
        if (Gflag) then
          write(*,*) 'Using QQGG+G (REAL) matrix elements'
          write(*,*) '[LC is     N   ]'
          write(*,*) '[SLC is   1/N  ]'
          write(*,*) '[SSLC is 1/N**3]'
        endif
        if (Qflag) then
          write(*,*) 'Using QQBQQB+G (REAL) matrix elements'
          write(*,*) '[LC is   1 ]'
          write(*,*) '[SLC is 1/N]'
        endif
        if     (colourchoice .eq. 1) then
          write(*,*) 'Leading colour only in REAL'
        elseif (colourchoice .eq. 2) then
          write(*,*) 'Sub-leading colour only in REAL'
        elseif (colourchoice .eq. 3) then
          write(*,*) 'Sub-sub-leading colour only in REAL'
        elseif (colourchoice .eq. 0) then
          write(*,*) 'Total of all colour structures in REAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
 
      call spinoru(7,p,za,zb)
      propsq=s(3,4)**2/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      if (Gflag) then
************************************************************************
*     Calculate contributions from the QQGG matrix elements            *
************************************************************************

      call spinoru(7,p,za,zb)
      call xwqqggg(5,1,2,7,6,3,4,mmsq_gg)
      call xwqqggg(1,5,6,7,2,3,4,mmsq_qqb)
      call xwqqggg(2,5,6,7,1,3,4,mmsq_qbq)
      call xwqqggg(1,2,6,7,5,3,4,mmsq_qg)
      call xwqqggg(2,1,6,7,5,3,4,mmsq_gq)
      call xwqqggg(5,1,6,7,2,3,4,mmsq_gqb)
      call xwqqggg(5,2,6,7,1,3,4,mmsq_qbg)
      
      endif
      
      if (Qflag) then
************************************************************************
*     Calculate contributions from the QQBQQB matrix elements          *
************************************************************************

c      call msq_WqqQQg(5,1,2,6,7,4,3,msq_qqbs)
c      call msq_WqqQQg(6,1,5,2,7,4,3,msq_qqbs)
c      call msq_WqqQQg(1,6,5,2,7,4,3,msq_qbqs)

c      call msq_WqqQQg(2,1,5,6,7,4,3,msq_qqb)
c      call msq_WqqQQg(1,2,5,6,7,4,3,msq_qbq)

c      call msq_WqqQQg(1,5,2,6,7,4,3,msq_qbqb)
      ofac=8d0*gsq**3*gwsq**2*aveqq
      call msq_WqqQQg(5,1,6,2,7,4,3,msq_qq,ofac)
      
c      call msq_WqqQQg(7,1,5,6,2,4,3,msq_qg)
c      call msq_WqqQQg(2,7,5,6,1,4,3,msq_gqb)
c      call msq_WqqQQg(7,2,5,6,1,4,3,msq_gq)
c      call msq_WqqQQg(1,7,5,6,2,4,3,msq_qbg)

      endif

C---call spinor routine and load common block twopij

c--- note the factor of 4d0*xw**2 relative to wbb
      fac=4d0*gsq**3*(gwsq/2d0)**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8d0
       
      LRb(1)=L(1)
      LRb(2)=R(1)

      do j=-nf,nf
      do k=-nf,nf
c      do j=2,2
c      do k=2,2
      
      msq(j,k)=0d0

      if (Gflag) then
************************************************************************
*     Sum the contributions from the QQGG matrix elements              *
************************************************************************

c--- note the identical particle factor of 1/6 for the
c--- q-qb initial states, due to 3 gluons in the final state     
      if     ((j .eq. 0) .and. (k .eq. 0)) then
        Vfac=0d0
        do n1=1,nf
          do n2=-nf,-1
            Vfac=Vfac+Vsq(n1,n2)
          enddo
        
        enddo
        msq(j,k)=propsq*mmsq_gg*Vfac*(gwsq**2/4d0/esq**2)
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
        msq(j,k)=propsq*mmsq_qqb*Vsq(j,k)
     .            *(aveqq/avegg)*(gwsq**2/4d0/esq**2)/6d0
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        msq(j,k)=propsq*mmsq_qbq*Vsq(j,k)
     .            *(aveqq/avegg)*(gwsq**2/4d0/esq**2)/6d0
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
        msq(j,k)=half*propsq*mmsq_qg
     .            *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
     .            *(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
        msq(j,k)=half*propsq*mmsq_qbg
     .            *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
     .            *(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        msq(j,k)=half*propsq*mmsq_gq
     .            *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
     .            *(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
        msq(j,k)=half*propsq*mmsq_gqb
     .            *(aveqg/avegg)*(gwsq**2/4d0/esq**2) 
     .            *(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))
      endif
      endif
      
      if (Qflag) then
************************************************************************
*     Sum the contributions from the QQBQQB matrix elements            *
************************************************************************

c--- note the factor of 4d0*xw**2 relative to wbb
      fac=gsq**3*gwsq**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8d0

      do f3=1,5
      do f4=f3,5
      if (f3.eq.f4) then
        stat=0.5d0
      else
        stat=1d0
      endif
      if     ((j .eq. 0) .and. (k .eq. 0)) then
c--- no glue-glue contribution here
      elseif ((j .gt. 0) .and. (k .gt. 0)) then
      if (j .le. k) then
      msq(j,k)=msq(j,k)+stat*fac*aveqq*msq_qq(j,f3,k,f4)
      endif
      if (j .gt. k) then
      msq(j,k)=msq(j,k)+stat*fac*aveqq*msq_qq(j,f4,k,f3)
      endif
      if (msq_qq(j,f3,k,f4) .ne. 0d0) then
c        write(*,*) j,k,': contribution from ',f3,f4
      endif
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=msq(j,k)+fac*aveqq*msq_qqbs(j,f3,-k,f4)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=msq(j,k)+fac*aveqq*msq_qbq(-j,k,f3,f4)
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
      msq(j,k)=msq(j,k)+stat*fac*aveqq*msq_qbqb(f3,-j,f4,-k)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
      do kk=1,5
      msq(j,k)=msq(j,k)+fac*aveqg*msq_qg(j,kk,f3,f4)
      enddo
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
      do kk=1,5
      msq(j,k)=msq(j,k)+fac*aveqg*msq_qbg(-j,kk,f3,f4)
      enddo
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
      do jj=1,5
      msq(j,k)=msq(j,k)+fac*aveqg*msq_qg(jj,k,f3,f4)
      enddo
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
      do jj=1,5
      msq(j,k)=msq(j,k)+fac*aveqg*msq_qbg(jj,-k,f3,f4)
      enddo
      endif

      enddo
      enddo

      endif

      enddo
      enddo

      return
      end
