      subroutine qqb_3jet(p,msq)
      implicit none
c----Lowest order matrix element for 3 jet production
C----averaged over initial colours and spins and with statistical
C--- factor for final states
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer j,k,m(-5:5)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  Biga,Bigb,Bigc,Bigd,
     . gg_ggg,gg_aqg,aq_ggg,
     . qa_ggg,ga_agg,ag_agg,qg_qgg,gq_qgg,qq_qqg,qa_qag,aq_aqg,
     . aa_aag,qr_qrg,qa_rbg,qb_qbg,bq_bqg,
     . qg_qrb,gq_qrb,ag_rba,ga_rba,
     . qg_qqa,gq_qqa,ag_qaa,ga_qaa,ab_abg,aq_brg

      data m/-1,-2,-1,-2,-1,0,1,2,1,2,1/ 

      call dotem(5,p,s)

      gg_ggg=+Bigd()/6d0

      qa_ggg=+Bigc(2,1,3,4,5)/6d0
      aq_ggg=+Bigc(1,2,3,4,5)/6d0
      ga_agg=-Bigc(3,2,1,4,5)/2d0
      ag_agg=-Bigc(3,1,2,4,5)/2d0
      qg_qgg=-Bigc(1,3,2,4,5)/2d0
      gq_qgg=-Bigc(2,3,1,4,5)/2d0
      gg_aqg=+Bigc(4,3,1,2,5)

      qa_qag=+Bigb(1,4,3,2,5)
c      aq_qag=+Bigb(2,4,3,1,5)
      aq_aqg=+Bigb(2,3,4,1,5)
      qq_qqg=+Bigb(1,2,3,4,5)/2d0
      aa_aag=+Bigb(3,4,1,2,5)/2d0
      qg_qqa=-Bigb(1,5,3,4,2)/2d0
      ag_qaa=-Bigb(5,4,1,3,2)/2d0
      gq_qqa=-Bigb(2,5,3,4,1)/2d0
      ga_qaa=-Bigb(5,4,2,3,1)/2d0


c      aq_rbg=+Biga(2,4,1,3,5)
      aq_brg=+Biga(2,3,1,4,5)
      qa_rbg=+Biga(1,4,2,3,5)
      qb_qbg=+Biga(1,4,3,2,5)
c      bq_qbg=+Biga(2,4,3,1,5)
      bq_bqg=+Biga(2,3,4,1,5)

      qr_qrg=+Biga(1,2,3,4,5)
      ab_abg=+Biga(3,4,1,2,5)

      qg_qrb=-Biga(1,5,3,4,2)
      gq_qrb=-Biga(2,5,3,4,1)
      ag_rba=-Biga(5,4,1,3,2)
      ga_rba=-Biga(5,4,2,3,1)


      do j=-2,2
      do k=-2,2
c--set msq=0 to initalize
      msq(j,k)=0d0
C--gg      
      if     ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=avegg*(gg_ggg+dfloat(nf)*gg_aqg)
C--qq      
      elseif ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aveqq*qq_qqg
          else
            msq(j,k)=aveqq*qr_qrg
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aveqq*(qa_qag+dfloat(nf-1)*qa_rbg+qa_ggg)
          else
            msq(j,k)=aveqq*qb_qbg
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aveqq*aa_aag
          else
            msq(j,k)=aveqq*ab_abg
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aveqq*(aq_aqg+dfloat(nf-1)*aq_brg+aq_ggg)
          else
            msq(j,k)=aveqq*bq_bqg
          endif

C--qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=aveqg*(qg_qgg+dfloat(nf-1)*qg_qrb+qg_qqa)
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=aveqg*(ag_agg+dfloat(nf-1)*ag_rba+ag_qaa)
C--gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=aveqg*(gq_qgg+dfloat(nf-1)*gq_qrb+gq_qqa)
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=aveqg*(ga_agg+dfloat(nf-1)*ga_rba+ga_qaa)
      endif

      enddo
      enddo
      
      
c--- the other flavour combinations are easy now
      do j=-nf,nf
      do k=-nf,nf
        if ((abs(j) .gt. 2) .or. (abs(k) .gt. 2)) then
          msq(j,k)=msq(m(j),m(k))
        endif
      enddo
      enddo



      return
      end
