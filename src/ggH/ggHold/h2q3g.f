      subroutine h2q3g(p1,p2,p3,p4,p5,Hqaggg)
      implicit none
C-----calculates the matrix element squared
C-----for q(p1)+q~(p2)-->g(p3)+g(p4)+g(p5)+H
C-----has to be preceded by a call to spinoru to set up the za,zb
      include 'constants.f'
      include 'zprods_com.f'
      double precision Hqaggg
      double complex qamp(6,2,2,2,2),temp(6,2,2,2,2), Appp
!     . Apmp,Appm,Ampp ! no longer used
      double complex A2q3g_mmmpp, A2q3g_mmpmp, A2q3g_mpmmp
      integer i3(6),i4(6),i5(6),j,k,h1,h3,h4,h5,p1,p2,p3,p4,p5,n(5)
      data i3/3,3,4,4,5,5/
      data i4/4,5,3,5,3,4/
      data i5/5,4,5,3,4,3/

      double precision DJK(6,6),xa,xb,xc,xd
c--- Full result
      parameter(xa=xn*cf**2,xb=-cf/2d0,xc=0.25d0/xn,xd=(xn**2+1d0)*xc)
c--- Leading in colour
c      parameter(xa=xn**3/4d0,xb=0d0,xc=0d0,xd=0d0)
c--- 1/N^2 suppressed
c      parameter(xa=-xn/2d0,xb=-xn/4d0,xc=0d0,xd=xn/4d0)
c--- 1/N^4 suppressed
c      parameter(xa=0.25d0/xn,xb=0.25d0/xn,xc=0.25d0/xn,xd=0.25d0/xn)
c--- Note that the definition of DJK here differs from the one in (B.22) of
c--- the paper by an overall factor of Cf which is restored in the sum below
      DATA  (DJK(J,1),J=1,6)/xa,xb,xb,xc,xc,xd/
      DATA  (DJK(J,2),J=1,6)/xb,xa,xc,xd,xb,xc/
      DATA  (DJK(J,3),J=1,6)/xb,xc,xa,xb,xd,xc/
      DATA  (DJK(J,4),J=1,6)/xc,xd,xb,xa,xc,xb/
      DATA  (DJK(J,5),J=1,6)/xc,xb,xd,xc,xa,xb/
      DATA  (DJK(J,6),J=1,6)/xd,xc,xc,xb,xb,xa/

      n(1)=p1
      n(2)=p2
      n(3)=p3
      n(4)=p4
      n(5)=p5

C----definition of helicities is for outgoing lines
C----labelling is as follows
C     temp(j,h1,h3,h4,h5) since h2 can be obtained from h1

      do j=1,6

      temp(j,2,2,2,2)=Appp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
     . za,zb)
      temp(j,2,1,1,1)=Appp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
     . zb,za)

      temp(j,2,2,2,1) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     . n(2),zb,za)
      temp(j,2,2,1,1) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     . n(1),za,zb)

      temp(j,2,2,1,2) = A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j))
     . ,n(2),zb,za)
      temp(j,2,1,2,1) = A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     . n(1),za,zb)

      temp(j,2,1,2,2) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     . n(2),zb,za)
      temp(j,2,1,1,2) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     .n(1),za,zb)


      temp(j,1,1,1,1)=Appp(n(1),n(2),n(i3(j)),n(i4(j)),n(i5(j)),
     . zb,za) 
      temp(j,1,2,2,2)=Appp(n(2),n(1),n(i5(j)),n(i4(j)),n(i3(j)),
     . za,zb)

      temp(j,1,2,1,1) = A2q3g_mpmmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     . n(2),za,zb)
      temp(j,1,2,2,1) = A2q3g_mpmmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     . n(1),zb,za)

      temp(j,1,1,1,2) = A2q3g_mmmpp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     . n(2),za,zb)
      temp(j,1,1,2,2) = A2q3g_mmmpp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     . n(1),zb,za)
      
      temp(j,1,1,2,1)= A2q3g_mmpmp(n(1),n(i3(j)),n(i4(j)),n(i5(j)),
     . n(2),za,zb)
      temp(j,1,2,1,2)= A2q3g_mmpmp(n(2),n(i5(j)),n(i4(j)),n(i3(j)),
     . n(1),zb,za)
            
      enddo

C----At this stage we have setup the amplitudes but failed 
C----to assign the helicities properly. So we now reshuffle 
C----to get these right.     
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      do j=1,6
      if (j.eq. 1) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h4,h5)
      if (j.eq. 2) qamp(j,h1,h3,h4,h5)=temp(j,h1,h3,h5,h4)
      if (j.eq. 3) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h3,h5)
      if (j.eq. 4) qamp(j,h1,h3,h4,h5)=temp(j,h1,h4,h5,h3)
      if (j.eq. 5) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h3,h4)
      if (j.eq. 6) qamp(j,h1,h3,h4,h5)=temp(j,h1,h5,h4,h3)
      enddo
      enddo
      enddo
      enddo
      enddo
 
c--- now perform the sum with the appropriate weights
      Hqaggg=0d0
      do j=1,6
      do k=1,6
      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      Hqaggg=Hqaggg
     .+Cf*djk(j,k)*dble(qamp(j,h1,h3,h4,h5)*dconjg(qamp(k,h1,h3,h4,h5)))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
