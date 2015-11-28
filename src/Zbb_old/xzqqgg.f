      subroutine xzqqgg(mqqb)
      implicit none
C     Author R.K.Ellis, January 2000
C     Returns the amplitudes squared for the process
C     0---> q(p1)+g(p2)+g(p3)+qbar(p4)+l(p5)+a(p6)
C     mqqb(2,2) has two indices;the first for the helicity quark line;
C     the second for helicity of lepton line.
      include 'constants.f'
      include 'prods.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer i2(2),i3(2),j,lh,h2,h3,hq,h(2:3)
      double precision mqqb(2,2),m1,m0,x,fac
      double complex tempm0,m(2),amp_qqgg
      parameter(x=xn/cf)
      data i2/2,3/
      data i3/3,2/

      
C ---final matrix element squared is needed as function of quark line helicity
C----and lepton line helicity
C----first argument is quark line helicity
C----second argument is lepton line helicity
      mqqb(1,1)=0d0
      mqqb(1,2)=0d0
      mqqb(2,1)=0d0
      mqqb(2,2)=0d0

      fac=avegg*2d0*gsq**2*esq**2*cf**2*xn
      do hq=1,2
      do lh=1,2

      do h2=1,2
      do h3=1,2
        h(2)=h2
        h(3)=h3
C initialize loop sum to zero
        tempm0=czip
         do j=1,2
         m(j)=amp_qqgg(1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),4,lh)
         tempm0=tempm0+m(j)
         enddo
      m1=Dble(Dconjg(m(1))*m(2))
      m0=abs(tempm0)**2
      mqqb(hq,lh)=mqqb(hq,lh)+fac*(m0-x*m1)
      enddo
      enddo

      enddo
      enddo

      return
      end
