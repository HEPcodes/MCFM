      subroutine hqqgg(p1,p2,p3,p4,ampsq)
      implicit none
C     Taken from Kauffman,Desai,Risal
C     PRD 55 1997 (4009)
C     and checked with hep-ph/9903330 
      include 'constants.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,j1,j2,j3
      double complex ab(2,2,2),ba(2,2,2),
     . abppp,abmmm,abppm,abmmp,bappm,bammp
      double precision ampsq
C====statement functions
C  Eq 25
      abppp(p1,p2,p3,p4)=
     . +(za(p2,p1)*zb(p1,p3)+za(p2,p4)*zb(p4,p3))**2*zb(p1,p4)
     . /((s(p1,p2)+s(p1,p4)+s(p2,p4))*za(p2,p4))
     . *(1d0/s(p1,p2)+1d0/s(p1,p4))
     . -(za(p2,p1)*zb(p1,p4)+za(p2,p3)*zb(p3,p4))**2
     . *zb(p1,p3)/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*za(p2,p3))
     . +(+za(p2,p3)*zb(p3,p1)+za(p2,p4)*zb(p4,p1))**2
     . /(zb(p1,p2)*za(p2,p3)*za(p2,p4)*za(p3,p4))

      abmmm(p1,p2,p3,p4)=
     . +(zb(p2,p1)*za(p1,p3)+zb(p2,p4)*za(p4,p3))**2
     . *za(p1,p4)/((s(p1,p2)+s(p1,p4)+s(p2,p4))
     . *zb(p2,p4))*(1d0/s(p1,p2)+1d0/s(p1,p4))
     . -(+zb(p2,p1)*za(p1,p4)+zb(p2,p3)*za(p3,p4))**2
     . *za(p1,p3)/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*zb(p2,p3))
     . +(+zb(p2,p3)*za(p3,p1)+zb(p2,p4)*za(p4,p1))**2
     . /(za(p1,p2)*zb(p2,p3)*zb(p2,p4)*zb(p3,p4))
C--Eq 26
      abppm(p1,p2,p3,p4)=
     . -za(p2,p4)**3/(za(p1,p2)*za(p2,p3)*za(p3,p4))
     . +zb(p1,p3)**3/(zb(p1,p2)*zb(p1,p4)*zb(p3,p4))
      abmmp(p1,p2,p3,p4)=
     . -zb(p2,p4)**3/(zb(p1,p2)*zb(p2,p3)*zb(p3,p4))
     . +za(p1,p3)**3/(za(p1,p2)*za(p1,p4)*za(p3,p4))

C--Eq 27
      bappm(p1,p2,p3,p4)=
     . -zb(p1,p3)**2*zb(p2,p3)/(zb(p1,p2)*zb(p2,p4)*zb(p3,p4))
     . +za(p1,p4)*za(p2,p4)**2/(za(p1,p2)*za(p1,p3)*za(p3,p4))
      
      bammp(p1,p2,p3,p4)=
     . -za(p1,p3)**2*za(p2,p3)/(za(p1,p2)*za(p2,p4)*za(p3,p4))
     . +zb(p1,p4)*zb(p2,p4)**2/(zb(p1,p2)*zb(p1,p3)*zb(p3,p4))


C====end statement functions
      
      ab(2,2,2)=abppp(p1,p2,p3,p4)
      ab(1,2,2)=abppp(p2,p1,p3,p4)
      ba(2,2,2)=abppp(p1,p2,p4,p3)
      ba(1,2,2)=abppp(p2,p1,p4,p3)

      ab(1,1,1)=abmmm(p1,p2,p3,p4)
      ab(2,1,1)=abmmm(p2,p1,p3,p4)
      ba(1,1,1)=abmmm(p1,p2,p4,p3)
      ba(2,1,1)=abmmm(p2,p1,p4,p3)

      ab(2,2,1)=abppm(p1,p2,p3,p4)
      ba(2,1,2)=abppm(p1,p2,p4,p3)
      ab(1,1,2)=abmmp(p1,p2,p3,p4)
      ba(1,2,1)=abmmp(p1,p2,p4,p3)

      ba(2,2,1)=bappm(p1,p2,p3,p4)
      ab(2,1,2)=bappm(p1,p2,p4,p3)
      ba(1,1,2)=bammp(p1,p2,p3,p4)
      ab(1,2,1)=bammp(p1,p2,p4,p3)

      ampsq=0d0
      do j1=1,2
      do j2=1,2
      do j3=1,2
      ampsq=ampsq
     .  +cf**2*xn*(cdabs(ab(j1,j2,j3))**2+cdabs(ba(j1,j2,j3))**2)
     .  -cf*dble(ab(j1,j2,j3)*Dconjg(ba(j1,j2,j3)))
      enddo
      enddo
      enddo

      return
      end
