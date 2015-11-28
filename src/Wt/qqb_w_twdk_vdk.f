      subroutine qqb_w_twdk_vdk(p,msqv)
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     Virtual corrections in the decay, averaged over initial          * 
*      colours and spins                                               *
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
      implicit none
      integer i,j,k
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'scheme.f'
      double precision p(mxpart,4),msqv(-nf:nf,-nf:nf),
     . virtgqdk,gq,qg,t(4),fac
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba

      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

c---fill matrices of spinor products

      do i=1,4
        t(i)=p(5,i)+p(6,i)+p(7,i)
      enddo

      call spinoru(7,p,za,zb)
      call spinork(7,p,zab,zba,t)

      fac=ason2pi*cf
      fac=aveqg*gwsq**4*gsq*V*fac
      gq=virtgqdk(1,2,3,4,5,6,7)
      qg=virtgqdk(2,1,3,4,5,6,7)

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msqv(j,k)=0d0
      if     ((j .eq. +5) .and. (k .eq. 0)) then
          msqv(j,k)=fac*qg
      elseif ((j .eq. 0) .and. (k .eq. +5)) then
          msqv(j,k)=fac*gq
      endif
      enddo
      enddo
      
      return
      end


      double precision function virtgqdk(ig,is,ie,in,je,jn,jc)
      implicit none

      integer ig,is,ie,in,je,jn,jc
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double precision snec,taugt,prop,mtsq,cv,ct,c1
      double complex amp(2),amp1(2),ampho(2)
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba

      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)
      taugt=s(ig,je)+s(ig,jn)+s(ig,jc)

      call coefsdk(s(jn,je),mtsq,ct,cv,c1)

      prop=(s(ie,in)-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((s(je,jn)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      
c--- tree level amplitudes for the decay t->Wb, "c0 contribution"
      amp(1) =
     &  - za(ie,in)*za(jc,jn)*zb(is,in)**2/zb(ig,is)*zab(ig,je)*
     & taugt**(-1)
     &  + za(ig,ie)*za(jc,jn)*zb(is,in)*zb(is,je)/zb(ig,is)*mtsq*
     & taugt**(-1)
     &

      amp(2) =
     &  - za(is,ie)*za(jc,jn)*zb(is,in)*zb(ig,je)/za(ig,is)*mtsq*
     & taugt**(-1)
     &  + za(jc,jn)*zb(is,in)/za(ig,is)*zab(is,je)*zab(ie,ig)*
     & taugt**(-1)
     &  + za(jc,jn)*zb(ig,in)/za(ig,is)*zab(ie,je)
     &

c--- virtual amplitudes for the decay t->Wb, "c1 contribution"
      amp1(1) =
     &  - za(ie,in)*za(ig,jc)*za(jc,jn)*zb(is,in)**2*zb(jc,je)/zb(ig,
     & is)*taugt**(-1)
     &  + za(ig,ie)*za(jc,jn)*zb(is,in)*zb(jc,je)/zb(ig,is)*zab(jc,is)
     & *taugt**(-1)
     &

      amp1(2) =
     &  - za(is,ie)*za(jc,jn)*zb(is,in)*zb(jc,je)/za(ig,is)*zab(jc,ig)
     & *taugt**(-1)
     &  + za(is,jc)*za(jc,jn)*zb(is,in)*zb(jc,je)/za(ig,is)*zab(ie,ig)
     & *taugt**(-1)
     &  + za(ie,jc)*za(jc,jn)*zb(ig,in)*zb(jc,je)/za(ig,is)
     &

      ampho(1)=dcmplx(ct+cv)*amp(1)+dcmplx(0.5d0*c1)*amp1(1)
      ampho(2)=dcmplx(ct+cv)*amp(2)+dcmplx(0.5d0*c1)*amp1(2)

      virtgqdk=dble(amp(1)*dconjg(ampho(1)))/prop
     .        +dble(amp(2)*dconjg(ampho(2)))/prop
      return
      end


