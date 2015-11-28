      subroutine ggtth(p,dens,dent,denx,deny,wtgg)
      implicit none
C***********************************************************************
C     Author: R.K. Ellis                                               *
C     December, 1999.                                                  *
C     calculate the element squared and subtraction terms              *
C     for the process                                                  *
c     My notation                                                      *
C     g(-p1) +g(-p2)=bbar(p6)+e-(p7)+nubar(p8)+b(p5)+nu(p3)+e+(p4)  *
C     +b(p9)+bbar(p10)                                                 *
C                                                                      *
C   There are eight diagrams
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprodx.f'
      integer q4,r7,s1,s2,t1,t2,x1,x2,j,pol1,pol2
      parameter(q4=3,r7=5,s1=6,s2=8,t1=9,t2=10,x1=1,x2=2)
C  for the moment only!!!!!!!!!!!!!!!!!!!!
      double precision wtgg,p(mxpart,4),dens,dent,denx,deny,nab(8)
      double complex dia(8,2,2),sump(2,2),summ(2,2),N1(2),N2(2),
     . bita1,bita2,bitb1,bitb2
      data nab/+1d0,+1d0,+1d0,-1d0,-1d0,-1d0,+2d0,+2d0/
      save nab
       call spinoru(10,p,za,zb)

      N1(1)=-1d0/zb(1,2)
      N1(2)=+1d0/za(1,2)

      N2(1)=-1d0/zb(2,1)
      N2(2)=+1d0/za(2,1)
C   L=1
C   R=2
C--diagram 1
      bita1=mt*(zb(1,t1)*za(t1,7)+zb(1,r7)*za(r7,7))
      bita2=mt*(zb(2,t2)*za(t2,7)+zb(2,r7)*za(r7,7))
      bitb1=za(1,t1)*zb(t1,r7)*za(r7,7)+mt**2*za(1,r7)
      bitb2=za(2,t2)*zb(t2,r7)*za(r7,7)+mt**2*za(2,r7)
      dia(1,1,1)=N1(1)*N2(1)/(dent*denx)*(
     . +zb(4,q4)*za(q4,1)*zb(2,x2)*za(x2,2)*bita1
     . +zb(4,q4)*za(q4,1)*mt*zb(2,1)*bitb2
     . +mt*zb(4,2)*mt*za(1,2)*bita1
     . +zb(4,q4)*za(1,x1)*zb(x1,1)*bitb2)
      dia(1,2,1)=N1(2)*N2(1)/(dent*denx)*(
     . +mt*zb(4,1)*za(2,x1)*zb(x1,2)*bita1
     . +zb(4,q4)*za(q4,2)*zb(1,x1)*zb(s1,2)*bitb1)
      dia(1,2,2)=N1(2)*N2(2)/(dent*denx)*(
     . +mt*zb(4,1)*za(2,x1)*zb(x1,2)*bitb1
     . +mt*zb(4,1)*mt*za(2,1)*bita2
     . +zb(4,q4)*za(q4,2)*mt*zb(1,2)*bitb1
     . +zb(4,q4)*za(q4,2)*zb(1,x1)*za(x1,1)*bita2)

C--diagram 2
C--diagram 3
C--diagram 4
C--diagram 5
C--diagram 6
C--diagram 7
      dia(7,1,2)=czip
      dia(7,2,1)=czip
      dia(7,1,1)=2d0*za(2,1)*zb(2,1)*mt/dens*(
     . +zb(4,q4)*za(q4,s2)*zb(s2,2)*za(2,7)
     . +zb(4,q4)*za(q4,2)*zb(2,r7)*za(r7,7)
     . +zb(4,s2)*za(s2,2)*zb(2,r7)*za(r7,7)
     . +mt**2*zb(4,2)*za(2,7))
      dia(7,2,2)=dia(7,1,1)
C--diagram 8
      dia(8,1,2)=czip
      dia(8,2,1)=czip
      dia(8,1,1)=2d0*za(2,1)*zb(2,1)*mt/dent*(
     . +zb(4,q4)*za(q4,2)*zb(2,t2)*za(t2,7)
     . +zb(4,q4)*za(q4,2)*zb(2,r7)*za(r7,7)
     . +zb(4,2)*za(2,t2)*zb(t2,r7)*za(r7,7)
     . +mt**2*zb(4,2)*za(2,7))
      dia(8,2,2)=dia(8,1,1)

      wtgg=0d0
      do pol1=1,2
      do pol2=1,2
      sump(pol1,pol2)=czip
      summ(pol1,pol2)=czip
      enddo
      enddo

      do pol1=1,2
      do pol2=1,2
      do j=1,8
      summ(pol1,pol2)=summ(pol1,pol2)+nab(j)*dia(j,pol1,pol2)
      if (j .lt. 7) sump(pol1,pol2)=sump(pol1,pol2)+dia(j,pol1,pol2)
      enddo
C Overall factor of V/4=2 taken outside nonab factor=N/4, ab factor=(V-1)/2/N
      wtgg=wtgg
     . +0.75d0*abs(summ(pol1,pol2))**2+7d0/6d0*abs(sump(pol1,pol2))**2

      enddo
      enddo
      return
      end
