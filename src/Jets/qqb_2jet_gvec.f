      subroutine qqb_twojet_gvec(p,n,in,msq)
C*********************************************************************** 
c     Author: R.K. Ellis                                               *
c     October, 2002.                                                   *
c     Matrix element for twojet production                             *
c     averaged over initial colours and spins                          *
c     with line "in" contracted with the vector n(mu)                  *
c     (orthogonal to p(in))                                            * 
c     q(-p1)+qbar(-p2)--> g(p3)+ g(p4)                                 *
C*********************************************************************** 
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k,in
C--in is the label of the parton dotted with n
      double precision msq(-nf:nf,-nf:nf),n(4),p(mxpart,4),fac,
     .  aq_gg,qa_gg,gg_gg,gg_qa,
     .  gq_gq,gq_qg,qg_qg,qg_gq,
     .  ga_ga,ga_ag,ag_ag,ag_ga,
     .  ss,tt,uu,smallc,smalld,nDn,nDp4 
      call dotem(4,p,s)
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      nDn=n(4)**2-n(1)**2-n(2)**2-n(3)**2
      fac=-0.5d0*nDn*gsq**2
      nDp4=n(4)*p(in,4)-n(3)*p(in,3)-n(2)*p(in,2)-n(1)*p(in,1)

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp4).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error for in=',in
         write(*,*) 'cutoff',1d-3*abs(p(in,4))
         write(6,*) 'nDp4',nDp4
         call flush(6)
         stop
      endif


      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      qg_gq=-fac*aveqg*smallc(uu,ss,tt)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      ag_ga=-fac*aveqg*smallc(uu,tt,ss)

      gq_gq=-fac*aveqg*smallc(tt,uu,ss)
      gq_qg=-fac*aveqg*smallc(uu,tt,ss)
      ga_ga=-fac*aveqg*smallc(tt,ss,uu)
      ga_ag=-fac*aveqg*smallc(uu,ss,tt)

      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)

      gg_gg=fac*avegg*smalld(ss,tt,uu)*half

      write(6,*) 'pause in qqb_2jet_gvec'
      write(6,*) 
     .  aq_gg,qa_gg,gg_gg,gg_qa,
     .  gq_gq,gq_qg,qg_qg,qg_gq,
     .  ga_ga,ga_ag,ag_ag,ag_ga
      pause
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      if (in.eq.1) then
         if (j.eq.0) then
                if (k.eq.0) msq(j,k)=gg_gg+dfloat(nf)*gg_qa
                if (k.gt.0) msq(j,k)=gq_qg
                if (k.lt.0) msq(j,k)=ga_ag
         endif
      elseif (in.eq.2) then
         if (k.eq.0) then
                if (k.eq.0) msq(j,k)=gg_gg+dfloat(nf)*gg_qa
                if (k.gt.0) msq(j,k)=qg_qg
                if (k.lt.0) msq(j,k)=ag_ag
         endif
      elseif (in.eq.3) then
C--qa      
         if ((j .gt. 0) .and. (k .eq. -j)) then
                msq(j,k)=qa_gg
C--aq      
         elseif ((j .lt. 0) .and. (k .eq. -j)) then
                msq(j,k)=aq_gg
C--qg_gq      
         elseif ((j .gt. 0) .and. (k .eq. 0)) then
                msq(j,k)=qg_gq
C--ag_ga      
         elseif ((j .lt. 0) .and. (k .eq. 0)) then
                msq(j,k)=ag_ga
         elseif ((j .eq. 0) .and. (k .eq. 0)) then
                msq(j,k)=gg_gg
         endif
      elseif (in.eq.4) then
C--qa      
         if ((j .gt. 0) .and. (k .eq. -j)) then
                msq(j,k)=qa_gg
C--aq      
         elseif ((j .lt. 0) .and. (k .eq. -j)) then
                msq(j,k)=aq_gg
C--qg_qg      
         elseif ((j .gt. 0) .and. (k .eq. 0)) then
                msq(j,k)=qg_qg
C--ag_ag      
         elseif ((j .lt. 0) .and. (k .eq. 0)) then
                msq(j,k)=ag_ag
         elseif ((j .eq. 0) .and. (k .eq. 0)) then
                msq(j,k)=gg_gg
         endif
      endif

      enddo
      enddo
      return
      end
