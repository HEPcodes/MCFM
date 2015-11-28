      subroutine xzqqggg(mqqb)
      implicit none
      include 'constants.f'
      include 'sprodx.f'
      include 'dprodx.f'
      integer i2(6),i3(6),i4(6),i6(2),i7(2),j,lh,h2,h3,h4,hq,h(2:4)
      double precision mqqb(2,2),m1,m2,m0,x,omx
      double complex tempm0(2),m(6),amp_qqggg
      parameter(x=0.5d0*xn/cf,omx=1d0-x)
      data i2/2,3,4,2,3,4/
      data i3/3,4,2,4,2,3/
      data i4/4,2,3,3,4,2/
      data i6/7,6/
      data i7/6,7/
C first argument is quark line helicity
C second argument is lepton line helicity
      
C ---final matrix element squared is needed as function of quark line helicity
C----and lepton line helicity
      
C initialize loop sums to zero
      do lh=1,2
      tempm0=czip
      m1=zip
      m2=zip
      do hq=1,2
      mqqb(hq)=0d0
      mqqb(hq)=0d0
      enddo
      enddo

      do hq=1,2
      do lh=1,2

      do h2=1,2
      do h3=1,2
      do h4=1,2
      h(2)=h2      
      h(3)=h3      
      h(4)=h4      
         do j=1,6
         h(i2(j))=h2
         h(i3(j))=h3
         h(i4(j))=h4
         m(j)=amp_qqggg(i1,hq,i2(j),h(i2(j),i3(j),h(i3(j),i4(j),i5,hq)
         tempm0=tempm0+m(j)
         m2=m2+abs(m(j))**2
         enddo
      m0=abs(tempm0)**2
      m1=-2d0*m2
     . +2d0*Dble(Dconjg(m(1))*(m(4)+m(5)-m(6)))
     . +2d0*Dble(Dconjg(m(2))*(m(5)+m(6)-m(4)))
     . +2d0*Dble(Dconjg(m(3))*(m(6)+m(4)-m(5)))

      mqqb(hq,lh)=mqqb(hq,lh)+cf**3*xn*(omx**2*m0-x*omx*m1+x**2*m2)
      enddo
      enddo
      enddo
      enddo

      return
      end
