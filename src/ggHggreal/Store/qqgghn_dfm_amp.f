c--- NOTE: THIS ROUTINE FAILS FOR SOME CROSSINGS. FOR NOW, WE JUST
c---  USE THE SQUARED MATRIX ELEMENTS FOR THIS PIECE
      subroutine qqgghn_dfm_amp(p1,p2,p3,p4,p,n,
     . qqgghn_ab,qqgghn_ba,qqgghn_sym)
      implicit none
c---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> H((p5+p6)+g(p3)+g(p4)
c   with momentum 4 contracted with the vector n(mu),
c   separated into colour orderings AB, BA and the sum
c     calculated by the program qqgghn.frm
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      double precision p(mxpart,4),n(4),nDn,nDp1,nDp2,nDp3,nDp4
      double precision s123,s124,qqgghn_ab,qqgghn_ba,qqgghn_sym
      integer p1,p2,p3,p4,j12,j3,j4,h1,h2,h3
      double complex ab(2,2,2),ba(2,2,2),abppp,abppm,abpmp
      double complex abmmm,bammp,abmmp
      double complex ndem,ndep,vecm(mxpart,mxpart)
      double precision ampsq

C====Taken from DFM hep-ph/0404013
C====statement functions
C  Eq A.9 #1
      abppp(p1,p2,p3,p4)=
     . +(za(p2,p1)*zb(p1,p4)+za(p2,p3)*zb(p3,p4))**2*zb(p1,p3)
     . /((s(p1,p2)+s(p1,p3)+s(p2,p3))*za(p2,p3))
     . *(1d0/s(p1,p2)+1d0/s(p1,p3))
     . -(za(p2,p1)*zb(p1,p3)+za(p2,p4)*zb(p4,p3))**2
     . *zb(p1,p4)/((s(p1,p2)+s(p1,p4)+s(p2,p4))*s(p1,p2)*za(p2,p4))
     . -(+za(p2,p4)*zb(p4,p1)+za(p2,p3)*zb(p3,p1))**2
     . /(zb(p1,p2)*za(p2,p4)*za(p2,p3)*za(p3,p4))

C  Eq A.9 #2
      abpmp(p1,p2,p3,p4)=
     . -(za(p2,p3)**3/(za(p1,p2)*za(p2,p4)*za(p3,p4))
     .  -zb(p1,p4)**3/(zb(p1,p2)*zb(p1,p3)*zb(p3,p4)))

C  Eq A.9 #3
      abppm(p1,p2,p3,p4)=
     . -(-zb(p1,p3)**2*zb(p2,p3)/(zb(p1,p2)*zb(p2,p4)*zb(p3,p4))
     .   +za(p1,p4)*za(p2,p4)**2/(za(p1,p2)*za(p1,p3)*za(p3,p4)))


C====end statement functions
      

      nDp4=n(4)*p(p4,4)-n(3)*p(p4,3)-n(2)*p(p4,2)-n(1)*p(p4,1)

c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      if (abs(nDp4).gt.1d-3*abs(p(1,4))) then 
         write(*,*) 'Error in qqgghn for :',p1,p2,p3,p4
         write(*,*) 'cutoff',1d-3*abs(p(1,4))
         write(6,*) 'nDp4',nDp4
         call flush(6)
         stop
      endif


      ab(2,2,2)=abppp(p1,p2,p3,p4)
      ab(1,1,1)=-dconjg(ab(2,2,2))

      ab(2,2,1)=abppm(p1,p2,p3,p4)
      ab(1,1,2)=-dconjg(ab(2,2,1))

      ab(2,1,2)=abpmp(p1,p2,p3,p4)
      ab(1,2,1)=-dconjg(ab(2,1,2))

      ab(1,2,2)=-abppp(p2,p1,p4,p3)
      ab(2,1,1)=-dconjg(ab(1,2,2))

      ba(2,2,2)=abppp(p1,p2,p4,p3)
      ba(1,1,1)=-dconjg(ba(2,2,2))

      ba(2,2,1)=abpmp(p1,p2,p4,p3)
      ba(1,1,2)=-dconjg(ba(2,2,1))

      ba(2,1,2)=abppm(p1,p2,p4,p3)
      ba(1,2,1)=-dconjg(ba(2,1,2))

      ba(1,2,2)=-abppp(p2,p1,p3,p4)
      ba(2,1,1)=-dconjg(ba(1,2,2))


      call ndveccur(p1,p4,n,p,vecm)
      ndep=vecm(p1,p4)/rt2/za(p1,p4)
      ndem=vecm(p4,p1)/rt2/zb(p4,p1)

      qqgghn_ab=0d0
      qqgghn_ba=0d0
      qqgghn_sym=0d0

c      write(6,*) 'qqgghn_dfm_amp:ab(2,2,2)',ab(2,2,2)
c      write(6,*) 'qqgghn_dfm_amp:ab(2,2,1)',ab(2,2,1)
c      write(6,*) 'qqgghn_dfm_amp:ab(2,1,2)',ab(2,1,2)
c      write(6,*) 'qqgghn_dfm_amp:ab(2,1,1)',ab(2,1,1)
c      write(6,*) 'qqgghn_dfm_amp:ab(1,2,2)',ab(1,2,2)
c      write(6,*) 'qqgghn_dfm_amp:ab(1,2,1)',ab(1,2,1)
c      write(6,*) 'qqgghn_dfm_amp:ab(1,1,2)',ab(1,1,2)
c      write(6,*) 'qqgghn_dfm_amp:ab(1,1,1)',ab(1,1,1)

c      write(6,*)
c      write(6,*) 'qqgghn_dfm_amp:ba(2,2,2)',ba(2,2,2)
c      write(6,*) 'qqgghn_dfm_amp:ba(2,2,1)',ba(2,2,1)
c      write(6,*) 'qqgghn_dfm_amp:ba(2,1,2)',ba(2,1,2)
c      write(6,*) 'qqgghn_dfm_amp:ba(2,1,1)',ba(2,1,1)
c      write(6,*) 'qqgghn_dfm_amp:ba(1,2,2)',ba(1,2,2)
c      write(6,*) 'qqgghn_dfm_amp:ba(1,2,1)',ba(1,2,1)
c      write(6,*) 'qqgghn_dfm_amp:ba(1,1,2)',ba(1,1,2)
c      write(6,*) 'qqgghn_dfm_amp:ba(1,1,1)',ba(1,1,1)
c      write(6,*)
c      call qqgghamp(p1,p2,p3,p4,p)
c      pause
      do j12=1,2
      do j3=1,2

      qqgghn_ab=qqgghn_ab
     .  +2d0*xn*cdabs(ndem*ab(j12,j3,2)+ndep*ab(j12,j3,1))**2
      qqgghn_ba=qqgghn_ba
     .  +2d0*xn*cdabs(ndem*ba(j12,j3,2)+ndep*ba(j12,j3,1))**2
      qqgghn_sym=qqgghn_sym
     .  -2d0/xn*cdabs(
     . +ndem*(ab(j12,j3,2)+ba(j12,j3,2))
     . +ndep*(ab(j12,j3,1)+ba(j12,j3,1)))**2

      enddo
      enddo

      return
      end
