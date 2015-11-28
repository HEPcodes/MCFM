      subroutine diag1(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T1)
      implicit none
C---returns the amplitude for the left-left process diagram1
C---only the two of the seven propagators are included.
c--- first index of T1 refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T1(2,2),fac,pr,pr1687,pr4351
      integer bpass,i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)

      fac=16d0*za(i5,i3)*zb(i8,i6)
      do bpass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it
      else
      j9=it
      jt=i9
      endif
      pr1687=pr(i1,i6,i8,i7)
      pr4351=pr(i4,i3,i5,i1)
      T1(bpass,1)=-fac*za(i2,j9)
     . *(+zb(jt,i2)*(pr(i4,i3,i5,i2)*pr1687-mt**2*zb(i4,i1)*za(i2,i7))
     .   +zb(jt,j9)*(pr(i4,i3,i5,j9)*pr1687-mt**2*zb(i4,i1)*za(j9,i7)))
      T1(bpass,2)=-fac*zb(i2,jt)
     . *(+za(j9,i2)*(pr4351*pr(i2,i6,i8,i7)-mt**2*zb(i4,i2)*za(i1,i7))
     .   +za(j9,jt)*(pr4351*pr(jt,i6,i8,i7)-mt**2*zb(i4,jt)*za(i1,i7)))
      enddo

      return
      end
       

      subroutine diag2(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T2)
C---returns the amplitude for the left-left process-diagram2
C---only the two of the seven propagators are included.
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T2(2,2),fac,pr,pr2687,pr4352
      integer bpass,i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)
      fac=+16d0*za(i5,i3)*zb(i8,i6)
      do bpass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it
      else
      j9=it
      jt=i9
      endif
      pr4352=pr(i4,i3,i5,i2)
      pr2687=pr(i2,i6,i8,i7)
      T2(bpass,1)=fac*zb(jt,i1)
     . *(+za(i1,j9)*(pr4352*pr(i1,i6,i8,i7)-mt**2*zb(i4,i1)*za(i2,i7))
     .   +za(jt,j9)*(pr4352*pr(jt,i6,i8,i7)-mt**2*zb(i4,jt)*za(i2,i7)))
      T2(bpass,2)=fac*za(j9,i1)
     . *(+zb(i1,jt)*(pr(i4,i3,i5,i1)*pr2687-mt**2*zb(i4,i2)*za(i1,i7))
     .   +zb(j9,jt)*(pr(i4,i3,i5,j9)*pr2687-mt**2*zb(i4,i2)*za(j9,i7)))
      enddo
      return
      end
       

      subroutine diag3(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T3)
      implicit none
C---returns the amplitude for the left-left and left right process 
C   diagram3
C---only the two of the seven propagators are included.
c--- first index refers to b-line(9,t)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T3(2,2),fac,part1,part2,part3,part4
      integer bpass,ipass,i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j1,j2,j9,jt
      do bpass=1,2
      do ipass=1,2

      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      else
      j9=it
      jt=i9 
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif

      part1=zb(jt,i6)*za(i6,i7)+zb(jt,i8)*za(i8,i7)
      part2=zb(i4,i3)*za(i3,j2)+zb(i4,i5)*za(i5,j2)
      part3=zb(j1,i6)*za(i6,j9)+zb(j1,i7)*za(i7,j9)+zb(j1,i8)*za(i8,j9)
     . +zb(j1,jt)*za(jt,j9)
      part4=za(j2,i6)*zb(i6,jt)+za(j2,i7)*zb(i7,jt)+za(j2,i8)*zb(i8,jt)
     . +za(j2,j9)*zb(j9,jt)
      fac=+16d0*za(i5,i3)*zb(i8,i6)

      T3(bpass,ipass)=fac
     . *(part2*part3*part1
     . +mt**2*(part2*zb(j1,jt)*za(j9,i7)-part1*zb(i4,j1)*za(j2,j9)
     . -part4*zb(i4,j1)*za(j9,i7)))

      enddo
      enddo
      return
      end

      subroutine diag4(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T4)
      implicit none
C---returns the amplitude for the left-left
C   diagram4
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T4(2,2),fac,part1,part2,part3,part4
      integer bpass,ipass,i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j1,j2,j9,jt

      do bpass=1,2
      do ipass=1,2
      
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      else
      j9=it
      jt=i9 
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif


      part1=zb(j1,i6)*za(i6,i7)+zb(j1,i8)*za(i8,i7)
      part2=zb(i4,i3)*za(i3,j9)+zb(i4,i5)*za(i5,j9)
      part3=-zb(jt,i3)*za(i3,j2)-zb(jt,i4)*za(i4,j2)-zb(jt,i5)*za(i5,j2)
     . -zb(jt,j9)*za(j9,j2)
      part4=-za(j9,i3)*zb(i3,j1)-za(j9,i4)*zb(i4,j1)-za(j9,i5)*zb(i5,j1)
     . -za(j9,jt)*zb(jt,j1)
      fac=+16d0*za(i5,i3)*zb(i8,i6)

      T4(bpass,ipass)=fac
     . *(part2*part3*part1
     . +mt**2*(part2*zb(jt,j1)*za(j2,i7)-part1*zb(i4,jt)*za(j9,j2)
     . -part4*zb(i4,jt)*za(j2,i7)))
      enddo
      enddo

      return
      end

      subroutine diag5(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T5)
      implicit none
C---returns the amplitude for the left-left
C   diagram5
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T5(2,2),fac,part1,part2,part3
      integer bpass,ipass,i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j1,j2,j9,jt
      
      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      else
      j9=it
      jt=i9 
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif

      part1=zb(i8,i6)*za(i6,j9)+zb(i8,jt)*za(jt,j9)
      part2=zb(i4,i3)*za(i3,j2)+zb(i4,i5)*za(i5,j2)
      part3=zb(j1,i6)*za(i6,i7)+zb(j1,i8)*za(i8,i7)+zb(j1,j9)*za(j9,i7)
c     . +zb(j1,jt)*za(jt,i7)
      fac=+16d0*za(i5,i3)*zb(jt,i6)*part1
      T5(bpass,ipass)=fac*(part3*part2-mt**2*zb(i4,j1)*za(j2,i7))
      enddo
      enddo
      return
      end
       

      subroutine diag6(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T6)
      implicit none
C---returns the amplitude for the left-left
C   diagram6
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T6(2,2),fac,part1,part2,part3
      integer bpass,ipass,i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j1,j2,j9,jt
      
      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      else
      j9=it
      jt=i9 
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif

      part1=zb(jt,i5)*za(i5,i3)+zb(jt,j9)*za(j9,i3)
      part2=zb(j1,i6)*za(i6,i7)+zb(j1,i8)*za(i8,i7)
      part3=zb(i4,i3)*za(i3,j2)+zb(i4,i5)*za(i5,j2)
     .     +zb(i4,j9)*za(j9,j2)+zb(i4,jt)*za(jt,j2)

      fac=16d0*zb(i8,i6)*za(i5,j9)*part1
      T6(bpass,ipass)=fac*(-part3*part2+mt**2*zb(i4,j1)*za(j2,i7))
      enddo
      enddo
      return
      end
       

      subroutine diag7(i1,i2,i3,i4,i5,i6,i7,i8,i9,it,T7)
      implicit none
C---returns the amplitude for the left-left process diagram1
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      double complex T7(2,2),fac,p,part1,part2,part3
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,it,j1,j2,j9,jt,
     . jw,jx,jy,jz,ipass,bpass
C--statement function
      p(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)
      fac=+16d0*za(i5,i3)*zb(i8,i6)
       
      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      else
      j9=it
      jt=i9 
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif


      part1=-0.5d0*(
     . +p(i4,i3,i5,j1)*p(j1,i6,i8,i7)+p(i4,i3,i5,j2)*p(j2,i6,i8,i7)
     . -p(i4,i3,i5,j9)*p(j9,i6,i8,i7)-p(i4,i3,i5,jt)*p(jt,i6,i8,i7)
     . -mt**2*(p(i4,j1,j2,i7)-p(i4,j9,jt,i7)))

      part2=-(p(i4,i3,i5,j9)*p(jt,i6,i8,i7)-mt**2*zb(i4,jt)*za(j9,i7))
     . *(za(j2,j9)*zb(j9,j1)+za(j2,jt)*zb(jt,j1))

      part3=+(p(i4,i3,i5,j2)*p(j1,i6,i8,i7)-mt**2*zb(i4,j1)*za(j2,i7))
     . *(za(j9,j1)*zb(j2,jt)+za(j9,j2)*zb(j2,jt))

      T7(bpass,ipass)=fac*(part1+part2+part3)

      enddo
      enddo

      return
      end




