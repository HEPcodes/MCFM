      subroutine h4q(p1,p2,p3,p4,ampsq,ampsqid)
      implicit none
C     Taken from Kauffman,Desai,Risal
C     PRD 55 1997 (4009)
c     q(-p1)+qp(-p2)--> h --> q(p3)+qp(p4)
C     returns overall matrix element squared
C     summed over colors and spins with factor of g^4*A^2 removed 
      include 'constants.f'
      include 'sprods_com.f'
      double precision ampsq,ampsqid,ammsq,ammsqi
      integer p1,p2,p3,p4
c      include 'zprods_com.f'
c      integer j1,j2
c      double complex a(2,2),b(2,2),amm
C  Eq 29 (suitably modified for my momentum configuration)
C====statement function
C--left-left amplitude
c      amm(p1,p2,p3,p4)=
c     . +za(p3,p4)**2/(za(p1,p3)*za(p2,p4))
c     . +zb(p1,p2)**2/(zb(p1,p3)*zb(p2,p4))
C--The above amplitude squared + color factor
      ammsq(p1,p2,p3,p4)=V/4d0
     . *((s(p1,p2)-s(p3,p4))**2+
     . (s(p1,p3)*s(p2,p4)+s(p3,p4)*s(p1,p2)-s(p1,p4)*s(p2,p3))**2
     .  /(s(p1,p3)*s(p2,p4)))/(s(p1,p3)*s(p2,p4))
C--The interference with color factor Eq.(A20)
      ammsqi(p1,p2,p3,p4)=cf/2d0
     . *((s(p1,p2)-s(p3,p4))**2
     . *(s(p1,p3)*s(p2,p4)+s(p1,p4)*s(p2,p3)-s(p1,p2)*s(p3,p4))
     . -2d0*(s(p1,p2)*s(p3,p4)+s(p2,p3)*s(p1,p4)-s(p1,p3)*s(p2,p4))
     .     *(s(p3,p4)*s(p1,p2)+s(p1,p3)*s(p2,p4)-s(p1,p4)*s(p2,p3))
     .  )/(s(p1,p3)*s(p2,p3)*s(p1,p4)*s(p2,p4))

c      a(1,1)=+amm(p1,p2,p3,p4)
c      a(2,1)=+amm(p3,p2,p1,p4)
c      a(1,2)=+amm(p1,p4,p3,p2)
c      a(2,2)=+amm(p3,p4,p1,p2)

C-----sign for FD statistics
c      b(1,1)=-amm(p1,p2,p4,p3)
c      b(2,1)=-amm(p4,p2,p1,p3)
c      b(1,2)=-amm(p1,p3,p4,p2)
c      b(2,2)=-amm(p4,p3,p1,p2)
c      ampsq=0d0
c      ampsqid=0d0
c      do j1=1,2
c      do j2=1,2
c      ampsq=ampsq+V/4d0*(cdabs(a(j1,j2)))**2
c      ampsqid=ampsqid+V/4d0*(cdabs(a(j1,j2))**2+cdabs(b(j1,j2))**2)
c      if (j1.eq.j2) ampsqid=ampsqid-cf*dble(a(j1,j2)*dconjg(b(j1,j2)))
c      enddo
c      enddo

      ampsq=ammsq(p1,p2,p3,p4)+ammsq(p3,p2,p1,p4)
     .     +ammsq(p1,p4,p3,p2)+ammsq(p3,p4,p1,p2)
      ampsqid=ampsq
     .     +ammsq(p1,p2,p4,p3)+ammsq(p4,p2,p1,p3)
     .     +ammsq(p1,p3,p4,p2)+ammsq(p4,p3,p1,p2)
     .     +ammsqi(p1,p2,p3,p4)+ammsqi(p3,p4,p1,p2)

      return
      end
