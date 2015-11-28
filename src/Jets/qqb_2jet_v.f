      subroutine qqb_2jet_v(p,msq)
      implicit none
c----Matrix element for W + jet production
c----in order alpha_s^2
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->f(p3)+f(p4)
c---
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qqij_ij,aaij_ij,aqij_ij,qaij_ij,
     .  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     .  aqii_jj,qaii_jj,
     .  gq_gq,ga_ga,qg_qg,ag_ag,gg_gg,gg_qa,qa_gg,aq_gg,ss,tt,uu,
     .  virta,virtb,virtc,virtd 
c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      scheme='dred'

      call dotem(4,p,s)
      fac=gsq**2*ason2pi
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=aveqq*virta(ss,tt,uu)
      aaij_ij=qqij_ij
      aqij_ij=aveqq*virta(uu,tt,ss)
      qaij_ij=aqij_ij
      qaii_jj=aveqq*virta(tt,ss,uu)
      aqii_jj=qaii_jj

      qqii_ii=fac*half*aveqq*virtb(ss,tt,uu)
      aaii_ii=qqii_ii
      qaii_ii=fac*aveqq*virtb(uu,tt,ss)
      aqii_ii=qaii_ii

c      gg_qa=+fac*avegg*virtc(ss,tt,uu)
      gg_qa=+fac*avegg*virtc(ss,uu,tt)
      qa_gg=+half*gg_qa
      aq_gg=+half*gg_qa

      qg_qg=-fac*aveqg*virtc(tt,ss,uu)
      gq_gq=+qg_qg

      ag_ag=-fac*aveqg*virtc(tt,uu,ss)
      ga_ga=+ag_ag

      gg_gg=fac*avegg*virtd(ss,tt,uu)

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=qaii_ii+dfloat(nf-1)*qaii_jj+qa_gg
          else
            msq(j,k)=qaij_ij
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aqii_ii+dfloat(nf-1)*aqii_jj+half*aq_gg
          else
            msq(j,k)=aqij_ij
          endif

C--qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_gg+dfloat(nf)*gg_qa
      endif

      enddo
      enddo

      return
      end
