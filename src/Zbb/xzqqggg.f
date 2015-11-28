      subroutine xzqqggg(j1,j2,j3,j4,j5,j6,j7,mqqb)
      implicit none
C     Author J.M.Campbell, February 2000
C     Returns the amplitudes squared for the process
C     0 ---> q(p1)+g(p2)+g(p3)+g(p4)+qbar(p5)+l(p6)+a(p7)
C     mqqb(2,2) has two indices;the first for the helicity quark line;
C     the second for helicity of lepton line.
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'lc.f'
      integer i2(6),i3(6),i4(6),j,lh,h2,h3,h4,hq,h(7)
      integer j1,j2,j3,j4,j5,j6,j7
      double precision mqqb(2,2),m1,m2,m0,fac
      double complex tempm0,m(6),amp_qqggg
      double complex mppppm,mpmmmm,mpppmm,mppmpm,
     .    mpmppm,mppmmm,mpmpmm,mpmmpm
c        parameter(x=0.5d0*xn/cf,omx=1d0-x)
c      data i2/2,2,4,3,3,4/
c      data i3/3,4,2,4,2,3/
c      data i4/4,3,3,2,4,2/
C first argument is quark line helicity
C second argument is lepton line helicity
      
C ---final matrix element squared is needed as function of quark line helicity
C----and lepton line helicity
      
      fac=avegg*gsq**3*esq**2*xn**3*cf*8d0
c--- extra factor of 8 due to colour matrix normalization (rt2**6)

      i2(1)=j2
      i3(1)=j3
      i4(1)=j4
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3
      i2(3)=j4
      i3(3)=j2
      i4(3)=j3
      i2(4)=j3
      i3(4)=j4
      i4(4)=j2
      i2(5)=j3
      i3(5)=j2
      i4(5)=j4
      i2(6)=j4
      i3(6)=j3
      i4(6)=j2

      do hq=1,2
      do lh=1,2
C initialize loop sums to zero
      mqqb(hq,lh)=0d0

      do h2=1,2
      do h3=1,2
      do h4=1,2

        h(j2)=h2
        h(j3)=h3
        h(j4)=h4
        tempm0=czip
        m2=zip

        do j=1,6
          m(j)=amp_qqggg(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     .                                       i4(j),h(i4(j)),j5,lh,j6,j7)
          tempm0=tempm0+m(j)
          m2=m2+abs(m(j))**2
        enddo
c        write(6,*) 'm2',m2  
        
        m0=cdabs(tempm0)**2
c--- check sign of the last three terms: original ver. had + vs. (B33)
        m1=-2d0*m2
     .   -2d0*Dble(Dconjg(m(1))*(m(2)+m(5)-m(6)))
     .   -2d0*Dble(Dconjg(m(4))*(m(5)+m(6)-m(2)))
     .   -2d0*Dble(Dconjg(m(3))*(m(6)+m(2)-m(5)))

c--- here we have (2,3,4)+(2,4,3)+(4,2,3) [4 is photon-like]
c---         plus (3,4,2)+(3,2,4)+(2,3,4) [2 is photon-like]
c---         plus (4,2,3)+(4,3,2)+(3,4,2) [3 is photon-like]
c--- (plus perms)
        m1=cdabs(m(1)+m(2)+m(3))**2
     .    +cdabs(m(4)+m(5)+m(1))**2
     .    +cdabs(m(3)+m(6)+m(4))**2
     .    +cdabs(m(5)+m(4)+m(6))**2
     .    +cdabs(m(2)+m(1)+m(5))**2
     .    +cdabs(m(6)+m(3)+m(2))**2

c--- note that the x defined in this routine differs by a factor
c--- of 2 from the x defined in Nagy-Trocs, cf. (39), (B30)
c        mqqb(hq,lh)=mqqb(hq,lh)+fac*(omx**2*m0-x*omx*m1+x**2*m2)
c--- re-written to make colour hierarchy explicit
      if (LConly) then
      mqqb(hq,lh)=mqqb(hq,lh)+fac*m2
      else        
      mqqb(hq,lh)=mqqb(hq,lh)+fac*(m2-m1/xnsq+(xnsq+1d0)*m0/xnsq**2)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      
      return
      end
