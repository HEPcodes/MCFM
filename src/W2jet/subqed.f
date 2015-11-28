*** NOTE: This subroutine isn't actually used, but can be used to check the
***  subleading amplitude generated in qqb_w2jet by two perms of subqcd

      subroutine subqed(i1,i2,i3,i4,i5,i6,za,zb,m)
c     the matrix element for the helicity amplitudes for the QED process
c     q(-p1)+a(-p4) --> q(p2)+l(p3)+g(p5)+g(p6)
c     multipled by ((a+l)^2-M**2)/(a+l)**2
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      real*8 pssq,qssq,alsq
      double complex m(-1:1,-1:1)
      double complex p1,p2,p3
      integer i1,i2,i3,i4,i5,i6

      pssq=s(i1,i5)+s(i1,i6)+s(i5,i6)
      qssq=s(i2,i5)+s(i2,i6)+s(i5,i6)
      alsq=s(i3,i4)

C--- qed results
      m(+1,+1)=four*za(i3,i2)**2/za(i4,i3)
     & *za(i2,i1)/(za(i5,i2)*za(i6,i2)*za(i5,i1)*za(i6,i1))

      p1=four*(
     & +za(i3,i2)*zb(i5,i1)/(za(i5,i1)*zb(i6,i1)*pssq)
     & *(zb(i4,i1)*za(i6,i1)+zb(i4,i5)*za(i6,i5))
     &  )/alsq

      p2=four*(
     & -(zb(i1,i2)*za(i3,i2)+zb(i1,i6)*za(i3,i6))
     & *(zb(i4,i1)*za(i2,i1)+zb(i4,i5)*za(i2,i5))
     & /(zb(i6,i2)*zb(i6,i1)*za(i5,i2)*za(i5,i1))
     &  )/alsq

      p3=four*(
     & +zb(i4,i1)*za(i6,i2)/(zb(i6,i2)*za(i5,i2)*qssq)
     & *(zb(i5,i2)*za(i3,i2)+zb(i5,i6)*za(i3,i6))
     &  )/alsq
 
      m(+1,-1)=p1-p2+p3

       p1=
     & four*(
     & +za(i3,i2)*zb(i6,i1)/(za(i6,i1)*zb(i5,i1)*pssq)
     &  *(zb(i4,i1)*za(i5,i1)+zb(i4,i6)*za(i5,i6)))/alsq
       p2=
     & four*(-(zb(i1,i2)*za(i3,i2)+zb(i1,i5)*za(i3,i5))
     & *(zb(i4,i1)*za(i2,i1)+zb(i4,i6)*za(i2,i6))
     & /(za(i6,i2)*za(i6,i1)*zb(i5,i2)*zb(i5,i1)))/alsq
       p3=
     & four*(+zb(i4,i1)*za(i5,i2)/(zb(i5,i2)*za(i6,i2)*qssq)
     & *(zb(i6,i2)*za(i3,i2)+zb(i6,i5)*za(i3,i5))
     & )/alsq
      m(-1,+1)=p1-p2+p3

      m(-1,-1)=four*zb(i4,i1)**2/zb(i4,i3)
     & *zb(i2,i1)/(zb(i5,i2)*zb(i6,i2)*zb(i5,i1)*zb(i6,i1))

      return
      end

