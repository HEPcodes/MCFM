      double precision function gggampsq() 
      implicit none
      include 'constants.f'
      integer i1(6),i2(6),i3(6),j,h1,h2,h3,hr1,hr2,hr3,h(3)
      double complex a(6,2,2,2)
      double complex ttbgggppp,ttbgggmpp,ttbgggpmp,ttbgggppm,
     .               ttbgggmmm,ttbgggpmm,ttbgggmpm,ttbgggmmp
      data i1/1,2,3,3,1,2/
      data i2/2,3,1,2,3,1/
      data i3/3,1,2,1,2,3/

      do j=1,6
      do h1=1,2
      do h2=1,2
      do h3=1,2
C     Consider for example ppm. Helicity of the gluon
C     in the last place is always negative (ie,1). 
C     Consider j=4. 
C     In that case gluon in last position is number 1
c     Hence h(1)=h3 
      h(i1(j))=h1
      h(i2(j))=h2
      h(i3(j))=h3

      if ((h1.eq.1).and.(h2.eq.1).and.(h3.eq.1))then
      a(j,h(1),h(2),h(3))=ttbgggmmm(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.2).and.(h2.eq.1).and.(h3.eq.1))then
      a(j,h(1),h(2),h(3))=ttbgggpmm(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.1).and.(h2.eq.2).and.(h3.eq.1))then
      a(j,h(1),h(2),h(3))=ttbgggmpm(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.1).and.(h2.eq.1).and.(h3.eq.2))then
      a(j,h(1),h(2),h(3))=ttbgggmmp(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.2).and.(h2.eq.2).and.(h3.eq.2))then
      a(j,h(1),h(2),h(3))=ttbgggppp(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.1).and.(h2.eq.2).and.(h3.eq.2))then
      a(j,h(1),h(2),h(3))=ttbgggmpp(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.2).and.(h2.eq.1).and.(h3.eq.2))then
      a(j,h(1),h(2),h(3))=ttbgggpmp(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      if ((h1.eq.2).and.(h2.eq.2).and.(h3.eq.1))then
      a(j,h(1),h(2),h(3))=ttbgggppm(i1(j),i2(j),i3(j),4,5,6,7)
      endif

      enddo
      enddo
      enddo
      enddo

      gggampsq=0d0
      do hr1=1,2
      do hr2=1,2
      do hr3=1,2
      gggampsq=gggampsq
     . +(xn**2-2d0)*(
     . cdabs(a(1,hr1,hr2,hr3))**2+cdabs(a(2,hr1,hr2,hr3))**2
     .+cdabs(a(3,hr1,hr2,hr3))**2+cdabs(a(4,hr1,hr2,hr3))**2
     .+cdabs(a(5,hr1,hr2,hr3))**2+cdabs(a(6,hr1,hr2,hr3))**2)

      gggampsq=gggampsq
     .  +cdabs(a(1,hr1,hr2,hr3)+a(2,hr1,hr2,hr3)
     .        +a(3,hr1,hr2,hr3)+a(4,hr1,hr2,hr3)
     .        +a(5,hr1,hr2,hr3)+a(6,hr1,hr2,hr3))**2/xn**2
      gggampsq=gggampsq
     . +(a(1,hr1,hr2,hr3)
     . *dconjg(a(4,hr1,hr2,hr3)-a(5,hr1,hr2,hr3)-a(6,hr1,hr2,hr3))
     . +a(2,hr1,hr2,hr3)
     . *dconjg(a(5,hr1,hr2,hr3)-a(6,hr1,hr2,hr3)-a(4,hr1,hr2,hr3))
     . +a(3,hr1,hr2,hr3)
     . *dconjg(a(6,hr1,hr2,hr3)-a(4,hr1,hr2,hr3)-a(5,hr1,hr2,hr3))
     . +a(4,hr1,hr2,hr3)
     . *dconjg(a(1,hr1,hr2,hr3)-a(2,hr1,hr2,hr3)-a(3,hr1,hr2,hr3))
     . +a(5,hr1,hr2,hr3)
     . *dconjg(a(2,hr1,hr2,hr3)-a(3,hr1,hr2,hr3)-a(1,hr1,hr2,hr3))
     . +a(6,hr1,hr2,hr3)
     . *dconjg(a(3,hr1,hr2,hr3)-a(1,hr1,hr2,hr3)-a(2,hr1,hr2,hr3)))
      
      enddo
      enddo
      enddo

C====V/8d0 (color factor) x (factor of 8 removed along with gsq in 
C====form program)
      gggampsq=V*gggampsq

      return
      end
