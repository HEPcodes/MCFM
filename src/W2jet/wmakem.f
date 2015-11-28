      subroutine wmakem(i1,i2,i3,i4,i5,i6,i7,
     . m1_1234,m2_1234,m3_3412,m4_3412)
C     Author: R.K. Ellis, March 2001
C     A subroutine calculating Nagy and Trocsanyi, PRD59 014020 (1999) 
C     Eq. A.23 with a factor of 2*i*e^2*g^3/s removed
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'ckm1.f'
      include 'ewcharge.f'
      include 'prods.f'
      include 'hardscale.f'
      integer fa,fb,fc,fd,ha,hb,hq,Qh,hg,lh,i1,i2,i3,i4,i5,i6,i7
      double precision vl(2),gq(2,2)
      double complex m1_1234(5,5,5,5,2,2,2),m2_1234(5,5,5,5,2,2,2),
     .               m3_3412(5,5,5,5,2,2,2),m4_3412(5,5,5,5,2,2,2)
      double complex a1(2,2,2,2),a2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2),
     . prop

      call nagy1(i1,i2,i3,i4,i5,i6,i7,a1,a2)
      call nagy2(i3,i4,i1,i2,i5,i6,i7,b3,b4)
      
      prop=s(i6,i7)/dcmplx((s(i6,i7)-wmass**2),wmass*wwidth)
      
      lh=1
      do fa=1,5
      do fb=1,5
      do fc=1,5
      do fd=1,5
        do ha=1,2
        do hb=1,2
        do hg=1,2
        if ((VV(fa,-fb) .ne. 0d0) .and. (fc.eq.fd)
     .    .and. (ha .eq. 1)) then
        m1_1234(fa,fb,fc,fd,ha,hb,hg)=
     .   VV(fa,-fb)*prop*a1(ha,hb,hg,lh)
        m2_1234(fa,fb,fc,fd,ha,hb,hg)=
     .   VV(fa,-fb)*prop*a2(ha,hb,hg,lh)
        else
        m1_1234(fa,fb,fc,fd,ha,hb,hg)=czip
        m2_1234(fa,fb,fc,fd,ha,hb,hg)=czip
        endif
        if ((VV(fc,-fd) .ne. 0d0) .and. (fa.eq.fb)
     .    .and. (hb .eq. 1)) then
        m3_3412(fc,fd,fa,fb,hb,ha,hg)=
     .   VV(fc,-fd)*prop*b3(hb,ha,hg,lh)
        m4_3412(fc,fd,fa,fb,hb,ha,hg)=
     .   VV(fc,-fd)*prop*b4(hb,ha,hg,lh)
        else
        m3_3412(fc,fd,fa,fb,hb,ha,hg)=czip
        m4_3412(fc,fd,fa,fb,hb,ha,hg)=czip
        endif

        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      enddo
      
      return
      end

