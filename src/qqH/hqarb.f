      subroutine hqarb(p1,p2,p3,p4,ampsq,ampsqid)
      implicit none
C     Taken from Kauffman
C     PRD 55 1997 (4009)
      h --> q(p1)+a(-p2)+r(p3)+b(p4)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,j1,j2
      double complex a(2,2),b(2,2),app
      double complex xa(mxpart,mxpart),xb(mxpart,mxpart)
      double precision ampsq
C  Eq 29
C====statement function
      app(p1,p2,p3,p4)=
     . +za(p2,p4)**2/(za(p1,p2)*za(p3,p4))
     . +zb(p1,p3)**2/(zb(p1,p2)*zb(p3,p4))

      ap(2,2)=app(p1,p2,p3,p4)
      ap(1,2)=app(p2,p1,p3,p4)
      ap(2,1)=app(p1,p2,p4,p3)
      ap(1,1)=app(p2,p1,p4,p3)

      bp(2,2)=app(p1,p4,p3,p2)
      bp(1,2)=app(p4,p1,p3,p2)
      bp(2,1)=app(p1,p4,p2,p3)
      bp(1,1)=app(p4,p1,p2,p3)



      ampsq=0d0
      do j1=1,2
      do j2=1,2
      ampsq=ampsq
     .  +cf**2*xn*(cdabs(a(j1,j2))**2+cdabs(ba(j1,j2,j3))**2)
     .  -cf*cdabs(ab(j1,j2,j3)*Dconjg(ba(j1,j2,j3)))
      enddo
      enddo
      enddo

      return

      end
