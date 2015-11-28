      subroutine diagg1(i1,i2,i9,it,T1)
      implicit none
C---returns the amplitude for the left-left process diagram1
C---only the two of the seven propagators are included.
c--- first index of T1 refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T1(2,2),fac,pr,pr1687,pr4351,den,polnorm
      integer bpass,i1,i2,i9,it,j9,jt,jw,jx,jy,jz
      double precision qDq,rDr,aDa,bDb,mtsq
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)

      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))

      den=den*s(i2,i9)
     . *(s(i1,i2)+s(i1,i9)+s(i2,i9))
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(qDq-mtsq)*(aDa-mtsq)


      do bpass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      fac=16d0*za(5,3)*zb(8,6)/(den*polnorm)

      pr1687=pr(i1,6,8,7)
      pr4351=pr(4,3,5,i1)
      T1(bpass,1)=fac*za(i2,j9)
     . *(+zb(jt,i2)*(-pr(4,3,5,i2)*pr1687+mt**2*zb(4,i1)*za(i2,7))
     .   +zb(jt,i9)*(-pr(4,3,5,i9)*pr1687+mt**2*zb(4,i1)*za(i9,7)))

      T1(bpass,2)=fac*za(i2,j9)
     . *(+zb(jt,i2)*(-pr4351*pr(i2,6,8,7)+mt**2*zb(4,i2)*za(i1,7))
     .   +zb(jt,i9)*(-pr4351*pr(i9,6,8,7)+mt**2*zb(4,i9)*za(i1,7)))
      enddo

      return
      end
       

      subroutine diagg2(i1,i2,i9,it,T2)
C---returns the amplitude for the left-left process-diagram2
C---only the two of the seven propagators are included.
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T2(2,2),fac,pr,pr2687,pr4352,den,polnorm
      double precision qDq,rDr,aDa,bDb,mtsq
      integer bpass,i1,i2,i9,it,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)

      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))


      den=den*s(i1,i9)*(s(i1,i2)+s(i1,i9)+s(i2,i9))
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(qDq-mt**2)*(aDa-mt**2)

      do bpass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      fac=16d0*za(5,3)*zb(8,6)/(den*polnorm)

      pr2687=pr(i2,6,8,7)
      pr4352=pr(4,3,5,i2)
      T2(bpass,1)=fac*zb(jt,i1)
     . *(+za(i1,j9)*(pr4352*pr(i1,6,8,7)-mt**2*zb(4,i1)*za(i2,7))
     .   +za(i9,j9)*(pr4352*pr(i9,6,8,7)-mt**2*zb(4,i9)*za(i2,7)))
      T2(bpass,2)=fac*za(j9,i1)
     . *(+zb(i1,jt)*(pr(4,3,5,i1)*pr2687-mt**2*zb(4,i2)*za(i1,7))
     .   +zb(i9,jt)*(pr(4,3,5,i9)*pr2687-mt**2*zb(4,i2)*za(i9,7)))
      enddo
      return
      end
       

      subroutine diagg3(i1,i2,i9,it,T3)
      implicit none
C---returns the amplitude for the left-left and left right process 
C   diagram3
c--- first index refers to gluon-line 9
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T3(2,2),pr,fac,pg687,p4352,p167899,m26789g,
     . den,polnorm
      double precision qDq,rDr,aDa,bDb,mtsq
      integer bpass,ipass,i1,i2,i9,it,j1,j2,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)

      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))

      den=den*s(i1,i2)
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(s(6,7)+s(6,8)+s(7,8)+s(6,i9)+s(7,i9)+s(8,i9)-mt**2)
     . *(qDq-mt**2)*(aDa-mt**2)


      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif


      pg687=pr(jt,6,8,7)
      p4352=pr(4,3,5,j2)
      p167899=pr(j1,6,7,j9)+pr(j1,8,i9,j9)
      m26789g=za(j2,6)*zb(6,jt)+za(j2,7)*zb(7,jt)
     .       +za(j2,8)*zb(8,jt)+za(j2,i9)*zb(i9,jt)
      fac=+16d0*za(5,3)*zb(8,6)/(den*polnorm)

      T3(bpass,ipass)=(p4352*(p167899*pg687+mt**2*zb(j1,jt)*za(j9,7))
     . -mt**2*zb(4,j1)*(m26789g*za(j9,7)+pg687*za(j2,j9)))*fac

      enddo
      enddo
      return
      end

      subroutine diagg4(i1,i2,i9,it,T4)
      implicit none
C---returns the amplitude for the left-left
C   diaggram4
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T4(2,2),fac,pr,p1687,p4359,m934591,pg34592,
     . den,polnorm
      double precision qDq,rDr,aDa,bDb,mtsq
      integer bpass,ipass,i1,i2,i9,it,j1,j2,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)

      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))

      den=den*s(i1,i2)
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(rDr-mt**2)*(qDq-mt**2)*(aDa-mt**2)

      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif



      p1687=pr(j1,6,8,7)
      p4359=pr(4,3,5,j9)
      pg34592=pr(jt,3,4,j2)+pr(jt,5,i9,j2)
      m934591=+za(j9,3)*zb(3,j1)+za(j9,4)*zb(4,j1)
     .       +za(j9,5)*zb(5,j1)+za(j9,i9)*zb(i9,j1)
      fac=+16d0*za(5,3)*zb(8,6)/(den*polnorm)


      T4(bpass,ipass)=fac
     . *(p4359*(-pg34592*p1687+mt**2*zb(jt,j1)*za(j2,7))
     . +mt**2*zb(4,jt)*(-za(j9,j2)*p1687+m934591*za(j2,7)))
      enddo
      enddo

      return
      end

      subroutine diagg5(i1,i2,i9,it,T5)
      implicit none
C---returns the amplitude for the left-left
C   diagram5
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T5(2,2),fac,pr,p4352,p16897,den,polnorm
      double precision qDq,rDr,aDa,bDb,mtsq
      integer bpass,ipass,i1,i2,i9,it,j1,j2,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)
      
      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))

      den=den*s(i1,i2)*s(6,i9)
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(qDq-mt**2)*(bDb-mt**2)

      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif

      p4352=pr(4,3,5,j2)
      p16897=zb(j1,6)*za(6,7)+zb(j1,8)*za(8,7)+zb(j1,i9)*za(i9,7)
      fac=-16d0*za(5,3)*zb(jt,6)*pr(8,6,i9,j9)/(den*polnorm)
      T5(bpass,ipass)=fac*(-p4352*p16897+mt**2*zb(4,j1)*za(j2,7))

      enddo
      enddo
      return
      end
       

      subroutine diagg6(i1,i2,i9,it,T6)
      implicit none
C---returns the amplitude for the left-left
C   diagram6
c--- first index refers to gluon-line(9)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T6(2,2),pr,fac,pg953,p1687,p43592,den,polnorm
      double precision qDq,rDr,aDa,bDb,mtsq
      integer bpass,ipass,i1,i2,i9,it,j1,j2,j9,jt,jw,jx,jy,jz
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)

      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))

      den=den*s(i1,i2)*s(5,i9)
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(rDr-mt**2)*(aDa-mt**2)
      
      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif

      pg953=pr(jt,5,i9,3)
      p1687=pr(j1,6,8,7)
      p43592=zb(4,3)*za(3,j2)+zb(4,5)*za(5,j2)+zb(4,i9)*za(i9,j2)
      fac=16d0*zb(8,6)*za(5,j9)*pg953/(den*polnorm)
      T6(bpass,ipass)=fac*(-p43592*p1687+mt**2*zb(4,j1)*za(j2,7))
      enddo
      enddo
      return
      end
       

      subroutine diagg7(i1,i2,i9,it,T7)
      implicit none
C---returns the amplitude for the left-left process diagram1
c--- first index refers to b-line(9,10)
c--- second index refers to initial-line(1,2)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      double complex T7(2,2),den,fac,pr,part1,part2,part3,polnorm
      double precision qDq,rDr,aDa,bDb,mtsq
      integer i1,i2,i9,it,j1,j2,j9,jt,jw,jx,jy,jz,ipass,bpass
C--statement function
      pr(jw,jx,jy,jz)=zb(jw,jx)*za(jx,jz)+zb(jw,jy)*za(jy,jz)
       
      qDq=s(3,4)+s(3,5)+s(4,5)
      rDr=qDq+s(3,i9)+s(4,i9)+s(5,i9)
      aDa=s(6,7)+s(6,8)+s(7,8)
      bDb=aDa+s(6,i9)+s(7,i9)+s(8,i9)
      mtsq=mt**2
      den=dcmplx(one,mt*twidth/(qDq-mtsq))
     .   *dcmplx(one,mt*twidth/(rDr-mtsq))
     .   *dcmplx(one,mt*twidth/(aDa-mtsq))
     .   *dcmplx(one,mt*twidth/(bDb-mtsq))


      den=den*s(i1,i2)
     . *(s(i1,i2)+s(i1,i9)+s(i2,i9))
     . *dcmplx((s(3,4)-wmass**2),wmass*wwidth)
     . *dcmplx((s(7,8)-wmass**2),wmass*wwidth)
     . *(qDq-mt**2)*(aDa-mt**2)
      

      do bpass=1,2
      do ipass=1,2
      if (bpass .eq. 1) then
      j9=i9
      jt=it 
      polnorm=sqrt(2d0)*zb(j9,jt)
      else
      j9=it
      jt=i9 
      polnorm=sqrt(2d0)*za(j9,jt)
      endif

      if (ipass .eq. 1) then
      j1=i1
      j2=i2 
      else
      j1=i2
      j2=i1 
      endif
      fac=16d0*za(5,3)*zb(8,6)/(den*polnorm)

      part1=-(-pr(4,3,5,i9)*pr(i9,6,8,7)+mt**2*zb(4,i9)*za(i9,7))
     . *za(j2,j9)*zb(jt,j1)

      part2=+(-pr(4,3,5,j9)*pr(jt,6,8,7)+mt**2*zb(4,jt)*za(j9,7))
     . *za(j2,i9)*zb(i9,j1)

      part3=-(-pr(4,3,5,j2)*pr(j1,6,8,7)+mt**2*zb(4,j1)*za(j2,7))
     . *(0.5d0*zb(jt,i9)*za(i9,j9)+pr(jt,i1,i2,j9))

      T7(bpass,ipass)=fac*(part1+part2+part3)
     
      enddo
      enddo

      return
      end




