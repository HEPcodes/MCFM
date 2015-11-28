************************************************************************
*     This is the W+ routine                                           *
************************************************************************
      subroutine qqb_wp2jet_g(p,msq)
************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2001.                                                      *
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
      include 'hardscale.f'
      include 'flags.f'
      include 'lc.f'
      integer j,k,n1,n2
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision mmsq_gg,mmsq_qqb,mmsq_qbq,mmsq_qg
      double precision mmsq_gq,mmsq_gqb,mmsq_qbg
      double precision fac,LRb(2),Vfac
      double precision propsq
      double precision 
     . QQ_ud_dd,QQ_us_ds,QQ_uu_du,QQ_du_dd,QQ_su_sd,QQ_uu_ud,
     .QbQb_ds_us,QbQb_du_uu,QbQb_dd_ud,QbQb_sd_su,QbQb_ud_uu,QbQb_dd_du,
     . RRb_dd_du,RRb_ss_du,RRb_uu_du,RRb_ud_dd,RRb_ud_ss,RRb_ud_uu,
     . QQb_ud_dd,QQb_us_ds,QQb_ud_uu,QQb_du_dd,QQb_su_sd,QQb_uu_ud,
     . RbR_dd_du,RbR_ss_du,RbR_uu_du,RbR_ud_dd,RbR_ud_ss,RbR_ud_uu,
     . QbQ_ud_dd,QbQ_us_ds,QbQ_ud_uu,QbQ_du_dd,QbQ_su_sd,QbQ_uu_ud,
     . QG_d_ddu,QG_d_dsc,QG_u_ddd,QG_u_dcc,QG_u_duu,
     . GQ_d_ddu,GQ_d_dsc,GQ_u_ddd,GQ_u_dcc,GQ_u_duu,GQ_u_udu,
     . QbG_d_ddu,QbG_u_ucs,QbG_u_uud,QbG_d_uuu,QbG_d_ucc,QbG_d_udd,
     . GQb_d_ddu,GQb_u_ucs,GQb_u_uud,GQb_d_udd,GQb_d_ucc,GQb_d_uuu
      logical first
      integer jj(-nf:nf),kk(-nf:nf)
      data jj/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      data kk/-1,-2,-1,-2,-1,0,1,2,1,2,1/
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

c--- basic Q-Q amplitudes 
      call addhel(1,5,2,6,7,3,4,QQ_ud_dd,QQ_us_ds,QQ_uu_du)
      call addhel(2,6,1,5,7,3,4,QQ_du_dd,QQ_su_sd,QQ_uu_ud)
c--- basic Qb-Qb amplitudes 
      call addhel(5,1,6,2,7,3,4,QbQb_dd_ud,QbQb_ds_us,QbQb_du_uu)
      call addhel(6,2,5,1,7,3,4,QbQb_dd_du,QbQb_sd_su,QbQb_ud_uu)
c--- basic Q-Qb amplitudes
c--- annihilation 
      call addhel(6,5,1,2,7,3,4,RRb_dd_du,RRb_ss_du,RRb_uu_du)
      call addhel(1,2,6,5,7,3,4,RRb_ud_dd,RRb_ud_ss,RRb_ud_uu)
c--- scattering 
      call addhel(1,5,6,2,7,3,4,QQb_ud_dd,QQb_us_ds,QQb_ud_uu)
      call addhel(6,2,1,5,7,3,4,QQb_du_dd,QQb_su_sd,QQb_uu_ud)
c--- basic Qb-Q amplitudes
c--- annihilation 
      call addhel(6,5,2,1,7,3,4,RbR_dd_du,RbR_ss_du,RbR_uu_du)
      call addhel(2,1,6,5,7,3,4,RbR_ud_dd,RbR_ud_ss,RbR_ud_uu)
c--- scattering 
      call addhel(2,5,6,1,7,3,4,QbQ_ud_dd,QbQ_us_ds,QbQ_ud_uu)
      call addhel(6,1,2,5,7,3,4,QbQ_du_dd,QbQ_su_sd,QbQ_uu_ud)
c--- basic Q-G amplitudes
      call addhel(7,6,1,5,2,3,4,QG_d_ddu,QG_d_dsc,QG_u_ddd)
      call addhel(1,5,7,6,2,3,4,QG_u_ddd,QG_u_dcc,QG_u_duu)
c--- basic QB-G amplitudes
      call addhel(6,7,5,1,2,3,4,QbG_d_ddu,QbG_u_ucs,QbG_u_uud)
      call addhel(5,1,6,7,2,3,4,QbG_d_udd,QbG_d_ucc,QbG_d_uuu)
c--- basic G-Q amplitudes
      call addhel(7,5,2,6,1,3,4,GQ_d_ddu,GQ_d_dsc,GQ_u_udu)
      call addhel(2,6,7,5,1,3,4,GQ_u_ddd,GQ_u_dcc,GQ_u_duu)
c--- basic G-QB amplitudes
      call addhel(5,7,6,2,1,3,4,GQb_d_udd,GQb_u_ucs,GQb_u_uud)
      call addhel(6,2,5,7,1,3,4,GQb_d_ddu,GQb_d_ucc,GQb_d_uuu)

      endif

c--- note the factor of 4d0*xw**2 relative to wbb
      fac=4d0*gsq**3*(gwsq/2d0)**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8d0
       
      LRb(1)=L(1)
      LRb(2)=R(1)

      do j=-nf,nf
      do k=-nf,nf
      
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

      if     ((j .eq. 0) .and. (k .eq. 0)) then
c--- no glue-glue contribution here
      elseif ((j .gt. 0) .and. (k .gt. 0)) then
c--- Q Q --> Q Q
        if ((jj(j) .eq. 2) .and. (kk(k) .eq. 1)) then
          msq(j,k)=msq(j,k)+Vsq(j,-k)*QQ_ud_dd*0.5d0
     .           +(Vsum(j)-Vsq(j,-k))*QQ_us_ds
        elseif ((jj(j) .eq. 1) .and. (kk(k) .eq. 2)) then
          msq(j,k)=msq(j,k)+Vsq(k,-j)*QQ_du_dd*0.5d0
     .           +(Vsum(k)-Vsq(k,-j))*QQ_su_sd
        elseif ((jj(j) .eq. 2) .and. (kk(k) .eq. 2)) then
          if (j .eq. k) msq(j,k)=msq(j,k)+Vsum(j)*QQ_uu_du
          if (j .ne. k) msq(j,k)=msq(j,k)+Vsum(j)*QQ_us_ds
     .                                   +Vsum(k)*QQ_su_sd
        endif
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
c--- Qb Qb --> Qb Qb
        if ((jj(j) .eq. -2) .and. (kk(k) .eq. -1)) then
          msQ(j,k)=msq(j,k)+Vsq(k,-j)*QbQb_ud_uu*0.5d0
     .           +(Vsum(k)-Vsq(k,-j))*QbQb_sd_su
        elseif ((jj(j) .eq. -1) .and. (kk(k) .eq. -2)) then
          msq(j,k)=msq(j,k)+Vsq(j,-k)*QbQb_du_uu*0.5d0
     .           +(Vsum(j)-Vsq(j,-k))*QbQb_ds_us
        elseif ((jj(j) .eq. -1) .and. (kk(k) .eq. -1)) then
          if (j .eq. k) msq(j,k)=msq(j,k)+Vsum(j)*QbQb_dd_ud
          if (j .ne. k) msq(j,k)=msq(j,k)+Vsum(j)*QbQb_ds_us
     .                                   +Vsum(k)*QbQb_sd_su
        endif
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
c--- Q Qb --> Q Qb
        if ((jj(j).eq.2) .and. (kk(k).eq.-1)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*(RRb_ud_dd+RRb_ud_uu
     .                               +dfloat(nf-2)*RRb_ud_ss)
     .           +(Vsum(j)-Vsq(j,k))*QQb_us_ds
     .           +(Vsum(k)-Vsq(j,k))*QQb_su_sd
        elseif ((jj(j).eq.2) .and. (kk(k).eq.-2)) then
          if (j.eq.-k) then
            Vfac=0d0
            do n1=1,nf
            do n2=-nf,-1
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(j)*RRb_uu_du
     .                 +Vfac*RRb_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(j)*QQb_us_ds
          endif
        elseif ((jj(j).eq.1) .and. (kk(k).eq.-1)) then
          if (j.eq.-k) then
            Vfac=0d0
            do n1=1,nf
            do n2=-nf,-1
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(k)*RRb_dd_du
     .                 +Vfac*RRb_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(k)*QQb_su_sd
          endif
        endif
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
c--- Qb Q --> Q Qb
        if ((jj(j).eq.-1) .and. (kk(k).eq.2)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*(RbR_ud_dd+RbR_ud_uu
     .                               +dfloat(nf-2)*RbR_ud_ss)
     .           +(Vsum(k)-Vsq(j,k))*QbQ_us_ds
     .           +(Vsum(j)-Vsq(j,k))*QbQ_su_sd
        elseif ((jj(j).eq.-2) .and. (kk(k).eq.2)) then
          if (j.eq.-k) then
            Vfac=0d0
            do n1=-nf,-1
            do n2=1,nf
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(k)*RbR_uu_du
     .                 +Vfac*RbR_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(k)*QbQ_us_ds
          endif
        elseif ((jj(j).eq.-1) .and. (kk(k).eq.1)) then
          if (j.eq.-k) then
            Vfac=0d0
            do n1=-nf,-1
            do n2=1,nf
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(j)*RbR_dd_du
     .                 +Vfac*RbR_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(j)*QbQ_su_sd
          endif
        endif
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
c--- Q G --> Q q qb
        if( jj(j) .eq. 2) then
              msq(j,k)=msq(j,k)+Vsq(j,-j+1)*(
     .                   (aveqg/aveqq)*QG_u_duu
     .                  +dfloat(nf-2)*(aveqg/aveqq)*QG_u_dcc
     .                  +(aveqg/aveqq)*QG_u_ddd*0.5d0)
     .                         +(2d0-Vsq(j,-j+1))*
     .                   (aveqg/aveqq)*QG_d_dsc
        elseif( jj(j) .eq. 1) then
         msq(j,k)=msq(j,k)+Vsq(-j,j+1)*(aveqg/aveqq)*QG_d_ddu*0.5d0
     .                    +(2d0-Vsq(-j,j+1))*(aveqg/aveqq)*QG_d_dsc
        endif
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
c--- QB G --> QB qb q
        if( jj(j) .eq. -1) then
              msq(j,k)=msq(j,k)+Vsq(j,-j+1)*(
     .                   (aveqg/aveqq)*QbG_d_udd
     .                  +(aveqg/aveqq)*QbG_d_uuu*0.5d0
     .                  +dfloat(nf-2)*(aveqg/aveqq)*QbG_d_ucc)
     .                         +(2d0-Vsq(j,-j+1))*
     .                   (aveqg/aveqq)*QbG_u_ucs
        elseif( jj(j) .eq. -2) then
              msq(j,k)=msq(j,k)+(aveqg/aveqq)*QbG_u_uud*0.5d0
     .                         +(aveqg/aveqq)*QbG_u_ucs
        endif
c--- G Q --> Q q qb
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
        if( kk(k) .eq. 2) then
         msq(j,k)=msq(j,k)+Vsum(k)*(
     .              (aveqg/aveqq)*GQ_u_udu
     .             +dfloat(nf-2)*(aveqg/aveqq)*GQ_u_dcc
     .             +(aveqg/aveqq)*GQ_u_ddd*0.5d0
     .             +(aveqg/aveqq)*GQ_d_dsc)
        elseif( kk(k) .eq. 1) then
         msq(j,k)=msq(j,k)+Vsq(-k,k+1)*(aveqg/aveqq)*GQ_d_ddu*0.5d0
     .                    +(2d0-Vsq(-k,k+1))*(aveqg/aveqq)*GQ_d_dsc
        endif
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
c--- G QB --> QB qb q
        if( kk(k) .eq. -1) then
              msq(j,k)=msq(j,k)+Vsq(k,-k+1)*(
     .                   (aveqg/aveqq)*GQb_d_udd
     .                  +(aveqg/aveqq)*GQb_d_uuu*0.5d0
     .                  +dfloat(nf-2)*(aveqg/aveqq)*GQb_d_ucc)
     .                         +(2d0-Vsq(k,-k+1))*
     .                   (aveqg/aveqq)*GQb_u_ucs
        elseif( kk(k) .eq. -2) then
         msq(j,k)=msq(j,k)+Vsq(-k,k+1)*(aveqg/aveqq)*GQb_u_uud*0.5d0
     .                    +(2d0-Vsq(-k,k+1))*(aveqg/aveqq)*GQb_u_ucs
        endif
      endif

      endif

      enddo
      enddo

      return
      end


      subroutine addhel(i1,i2,i3,i4,i5,i6,i7,xmsq_ud_dd,
     . xmsq_us_ds,xmsq_uu_du)
      implicit none
      include 'constants.f'
      integer i1,i2,i3,i4,i5,i6,i7
      double precision xmsq_ud_dd,xmsq_us_ds,xmsq_uu_du,
     .uLdR_dLdRp,uLdR_dRdLp,uRuL_dLuRp,uLsL_dLsLp,uLdL_dLdLp,uLuL_dLuLp,
     .uLdR_dLdRm,uLdR_dRdLm,uRuL_dLuRm,uLsL_dLsLm,uLdL_dLdLm,uLuL_dLuLm

c--- basic set of +ve gluon helicity amplitudes for QQ -> QQ       
      call testem(i1,i2,i3,i4,i5,i6,i7,uLdR_dLdRp,uLdR_dRdLp,
     .           uRuL_dLuRp,uLsL_dLsLp,uLdL_dLdLp,uLuL_dLuLp)

c--- generate -ve gluon helicity amplitudes this way for now
c--- under this symmetry, uL dR -> dR dL  <===>  uR uL -> dL uR
c---                and   uL dL -> dL dL  <===>  uL uL -> dL uL
      call testem(i2,i1,i4,i3,i5,i7,i6,uLdR_dLdRm,uRuL_dLuRm,
     .           uLdR_dRdLm,uLsL_dLsLm,uLuL_dLuLm,uLdL_dLdLm)

      xmsq_ud_dd=(uLdR_dLdRp+uLdR_dRdLp+uLdL_dLdLp
     .           +uLdR_dLdRm+uLdR_dRdLm+uLdL_dLdLm)

      xmsq_us_ds=(uLdR_dLdRp+uLsL_dLsLp
     .           +uLdR_dLdRm+uLsL_dLsLm)

c--- note that uLuR_uLuR = uLdR_dLdR
      xmsq_uu_du=(uLdR_dLdRp+uRuL_dLuRp+uLuL_dLuLp
     .           +uLdR_dLdRm+uRuL_dLuRm+uLuL_dLuLm)

      return
      end
      
