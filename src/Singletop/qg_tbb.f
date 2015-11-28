      subroutine qg_tbb(p,msq)
C----q(p1)+g(p2)--> t(nu(p3)+e^+(p4))+b(p5))+bb(p6)+q'(p7)
      implicit none
      include 'constants.f'
      include 'hardscale.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     . q_g,qbar_g,g_q,g_qbar
      double precision p1(4),p2(4),p3(4),q1(4),q2(4),q3(4),k0(4),xmatsq
      double precision dot,pt,ptbsq,ptbbsq,s56
      include 'ckm.f'
      integer j,k

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
c----put basic cuts to avoid mass singularity for zero mass
      ptbsq=pt(5,p)**2
      ptbbsq=pt(6,p)**2
      s56=two*dot(p,5,6)
      if (
     .     (s56 .lt. four*hscalesq) 
     . .or. (ptbsq.lt. hscalesq) 
     . .or. (ptbbsq .lt. hscalesq) ) return


      do j=1,4
      p1(j)=p(1,j)
      p2(j)=p(7,j)
      p3(j)=p(6,j)
      q1(j)=p(3,j) 
      q2(j)=p(4,j) 
      q3(j)=p(5,j) 
      k0(j)=p(2,j)
      enddo

      q_g=xmatsq(p1,p2,p3,q1,q2,q3,k0)
      g_q=xmatsq(k0,p2,p3,q1,q2,q3,p1)

      qbar_g=xmatsq(p2,p1,p3,q1,q2,q3,k0)
      g_qbar=xmatsq(p2,k0,p3,q1,q2,q3,p1)


      do j=-nf,nf
      do k=-nf,nf
      if (k.eq.0) then
      if (j .gt. 0) msq(j,k)=Vsum(j)*q_g
      if (j .lt. 0) msq(j,k)=Vsum(j)*qbar_g
      endif
      if (j.eq.0) then
      if (k .gt. 0) msq(j,k)=Vsum(k)*g_q
      if (k .lt. 0) msq(j,k)=Vsum(k)*g_qbar
      endif
c      if (msq(j,k) .ne. 0d0) write(6,*) j,k,msq(j,k)

      enddo
      enddo

      return
      end

      double precision function xmatsq(p1,p2,p3,q1,q2,q3,k0)
c     matrix element with all momenta outgoing
c     ie the process
c     u(-p1)+g(-k0) --> d(p2)+c(p3)+n(q1)+e^+(q2)+b(q3)
c     averaged(summed) over all initial(final) spins and colors
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j
      double precision dot2,sym,det3,det4s
      double precision a11,a22,a33,a12,a13,a23
      double precision p1(4),p2(4),p3(4),q1(4),q2(4),q3(4),k0(4)
      double precision ps(4),qs(4)
      double precision pssq,qssq,propsq,proqsq,dd1,dd2,dd3

      double precision p1Dp2,p1Dp3,p1Dq3,p1Dk0
      double precision p2Dps,p2Dq2,p2Dqs
      double precision p3Dq1,p3Dq3,p3Dk0
      double precision psDps,psDq2
      double precision q1Dq2,q1Dq3,q1Dk0
      double precision q2Dqs
      double precision q3Dk0
      double precision qsDqs

      double precision mtsq
      do j=1,4
      ps(j)=p1(j)+p2(j)+p3(j)
      qs(j)=q1(j)+q2(j)+q3(j)
      enddo

      p1Dp2 = dot2(p1,p2)
      p1Dp3 = dot2(p1,p3)
      p1Dq3 = dot2(p1,q3)
      p1Dk0 = dot2(p1,k0)
      p2Dps = dot2(p2,ps)
      p2Dq2 = dot2(p2,q2)
      p2Dqs = dot2(p2,qs)
      p3Dq1 = dot2(p3,q1)
      p3Dq3 = dot2(p3,q3)
      p3Dk0 = dot2(p3,k0)
      psDps = dot2(ps,ps)
      psDq2 = dot2(ps,q2)
      q1Dq2 = dot2(q1,q2)
      q1Dq3 = dot2(q1,q3)
      q1Dk0 = dot2(q1,k0)
      q2Dqs = dot2(q2,qs)
      q3Dk0 = dot2(q3,k0)
      qsDqs = dot2(qs,qs)
     
      mtsq=mt**2
      pssq=psDps-mtsq
      qssq=qsDqs-mtsq
      propsq=(2*p1Dp2-wmass**2)**2
      proqsq=(2*q1Dq2-wmass**2)**2
c      write(6,*) pssq,qssq,mt
c      write(6,*) sqrt(propsq),sqrt(proqsq),wmass
c      pause
      a11=q1Dq3*p1Dp3/(pssq*qssq)**2*(
     1 +2*(qssq+pssq)*p2Dps*q2Dqs-qssq*pssq*p2Dq2
     1 +2*qsDqs*(p2Dps*psDq2-half*mtsq*p2Dq2)
     1 +2*psDps*(q2Dqs*p2Dqs-half*mtsq*p2Dq2))

      a22=half*q1Dq3/(p3Dk0*qssq)**2
     1 *SYM(p2,qs,q2,qs)*(mb**2*(p1Dk0+p1Dp3)-p1Dk0*p3Dk0)

      a33=half*p1Dp3/(q3Dk0*pssq)**2
     1 *SYM(q2,ps,p2,ps)*(mb**2*(q1Dk0+q1Dq3)-q1Dk0*q3Dk0)

      dd1=DET3(p2,q2,qs,p1,p3,k0)
      dd2=DET3(k0,p2,ps,k0,p1,p3)
      dd3=DET3(p2,ps,q2,p1,p3,k0)
      a12=half*q1Dq3/(p3Dk0*pssq*qssq**2)*(
     1 (2*p1Dp3+p1Dk0)*(2*q2Dqs*SYM(qs,p3,ps,p2)-qsDqs*SYM(q2,p3,ps,p2)
     1 -mtsq*SYM(p3,p2,qs,q2))
     1 +p1Dp3*(2*q2Dqs*SYM(qs,k0,ps,p2)-qsDqs*SYM(q2,k0,ps,p2)
     1 -mtsq*SYM(k0,p2,qs,q2))
     1 +p3Dk0*(-2*q2Dqs*SYM(qs,p1,ps,p2)+qsDqs*SYM(q2,p1,ps,p2)
     1 +mtsq*SYM(p1,p2,qs,q2))
     1 -mtsq*dd1+2*q2Dqs*dd2+qsDqs*dd3)

      dd1=DET3(q2,p2,ps,q1,q3,k0)
      dd2=DET3(q2,qs,k0,q1,q3,k0)
      dd3=DET3(q2,p2,qs,q1,q3,k0)
      a13=half*p1Dp3/(q3Dk0*qssq*pssq**2)*(
     1 (2*q1Dq3+q1Dk0)*(2*p2Dps*SYM(ps,q3,qs,q2)-psDps*SYM(p2,q3,qs,q2)
     1 -mtsq*SYM(q3,q2,ps,p2))
     1 +q1Dq3*(2*p2Dps*SYM(ps,k0,qs,q2)-psDps*SYM(p2,k0,qs,q2)
     1 -mtsq*SYM(k0,q2,ps,p2))
     1 +q3Dk0*(-2*p2Dps*SYM(ps,q1,qs,q2)+psDps*SYM(p2,q1,qs,q2)
     1 +mtsq*SYM(q1,q2,ps,p2))
     1 -mtsq*dd1+2*p2Dps*dd2-psDps*dd3)

      dd1=DET4S(k0,p2,q2,qs,k0,q1,q3,p3)
      dd2=DET4S(k0,q2,p2,ps,k0,p1,p3,q3)
      a23=half/(p3Dk0*q3Dk0*pssq*qssq)*(SYM(q2,qs,p2,ps)
     1 *(-q1Dq3*(p1Dq3*p3Dk0-p1Dk0*p3Dq3-p1Dp3*(p3Dq3+p3Dk0))
     1   -p1Dp3*(p3Dq1*q3Dk0-q1Dk0*p3Dq3-q1Dq3*(p3Dq3+q3Dk0)))
     1   +p1Dp3*dd1+q1Dq3*dd2)
      xmatsq=8d0*pi*as*gwsq**4*(a11+a22+a33+a12+a13+a23)/(propsq*proqsq)
      xmatsq=xmatsq*qssq**2/(qssq**2+(mt*twidth)**2)
      xmatsq=xmatsq*proqsq/(proqsq+(wmass*twidth)**2)
      return
      end

      double precision function SYM(p1,p2,p3,p4)
      implicit none
      double precision p1(4),p2(4),p3(4),p4(4),dot2
      SYM=dot2(p1,p2)*dot2(p3,p4)
     1 +dot2(p2,p3)*dot2(p1,p4)
     1 -dot2(p1,p3)*dot2(p2,p4)
      return
      end

      double precision function det4S(p1,p2,p3,p4,q1,q2,q3,q4)
      implicit none
      double precision p1(4),p2(4),p3(4),p4(4)
      double precision q1(4),q2(4),q3(4),q4(4)
      double precision det3,dot2
      det4S=
     1 -dot2(p2,q1)*det3(p1,p3,p4,q2,q3,q4)
     1 +dot2(p3,q1)*det3(p1,p2,p4,q2,q3,q4)
     1 -dot2(p4,q1)*det3(p1,p2,p3,q2,q3,q4)
      return
      end

c      double precision function det4(p1,p2,p3,p4,q1,q2,q3,q4)
c      implicit none
c      double precision p1(4),p2(4),p3(4),p4(4)
c      double precision q1(4),q2(4),q3(4),q4(4)
c      double precision det3,dot2
c      det4=
c     1 +dot2(p1,q1)*det3(p2,p3,p4,q2,q3,q4)
c     1 -dot2(p2,q1)*det3(p1,p3,p4,q2,q3,q4)
c     1 +dot2(p3,q1)*det3(p1,p2,p4,q2,q3,q4)
c     1 -dot2(p4,q1)*det3(p1,p2,p3,q2,q3,q4)
c      return
c      end

      double precision function det3(p1,p2,p3,q1,q2,q3)
      implicit none
      double precision p1(4),p2(4),p3(4),q1(4),q2(4),q3(4)
      double precision dot2,det2
      det3=
     1 +dot2(p1,q1)*det2(p2,p3,q2,q3)
     1 -dot2(p1,q2)*det2(p2,p3,q1,q3)
     1 +dot2(p1,q3)*det2(p2,p3,q1,q2)
      return
      end

      double precision function det2(p1,p2,q1,q2)
      double precision p1(4),p2(4),q1(4),q2(4),dot2
      det2=dot2(p1,q1)*dot2(p2,q2)-dot2(p1,q2)*dot2(p2,q1)
      return
      end

      DOUBLE PRECISION FUNCTION dot2(P1,P2)
      IMPLICIT NONE
      DOUBLE PRECISION P1(4),P2(4)
      dot2=P1(4)*P2(4)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3)
      RETURN
      END



