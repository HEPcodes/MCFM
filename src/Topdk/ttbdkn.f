      double precision function ttbdkn(j1,j2,ps)
      implicit none
C     p1 incoming gluon contracted with n
C     p2 incoming gluon
C     p3 top momentum
C     p4 positron momentum
C     p5 anti-top momentum
C     p6 vector n
C     p7 electron momentum
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer j1,j2
      double precision ps(mxpart,4),s12,nDn,mtsq,t1,t2,ro
c--- appropriate scale is approx 1d-3*energy(incoming)
c--- so of order(1) for the Tevatron
      call dotem(7,ps,s)
      s12=s(j1,j2)
      s12=s(j1,j2)
      ro=4d0*mt**2/s12
      t2=-s(j2,3)/s12
      t1=-s(j2,5)/s12
      nDn=ps(6,4)**2-ps(6,1)**2-ps(6,2)**2-ps(6,3)**2
      mtsq=mt**2
      if (abs(s(1,6)).gt.1d-3*abs(ps(j1,4))) then
         write(6,*) 'nDp1',0.5d0*s(1,6)
         call flush(6)
         stop
      endif
      ttbdkn =  + s12**(-1) * (  - s(3,4)*s(3,7)*s(6,3)**2*t1**(-1)*ro
     &     + s(3,4)*s(3,7)*s(6,3)**2*ro + 2.D0*s(3,4)*s(3,7)*s(6,3)*s(6
     &    ,5)*ro - s(3,4)*s(3,7)*s(6,5)**2*t2**(-1)*ro + s(3,4)*s(3,7)*
     &    s(6,5)**2*ro - s(3,4)*s(6,3)**2*s(j1,7)*t1**(-1)*ro + 1.D0/2.D
     &    0*s(3,4)*s(6,3)**2*s(j1,7)*ro - 1.D0/2.D0*s(3,4)*s(6,3)**2*s(
     &    j2,7)*t1**(-1)*ro + 1.D0/2.D0*s(3,4)*s(6,3)**2*s(j2,7)*ro + 1.
     &    D0/2.D0*s(3,4)*s(6,3)*s(6,5)*s(j1,7)*t2**(-1)*ro + s(3,4)*s(6
     &    ,3)*s(6,5)*s(j1,7)*ro + s(3,4)*s(6,3)*s(6,5)*s(j2,7)*ro - 1.D0
     &    /2.D0*s(3,4)*s(6,5)**2*s(j1,7)*t2**(-1)*ro + 1.D0/2.D0*s(3,4)
     &    *s(6,5)**2*s(j1,7)*ro - 1.D0/2.D0*s(3,4)*s(6,5)**2*s(j2,7)*
     &    t2**(-1)*ro + 1.D0/2.D0*s(3,4)*s(6,5)**2*s(j2,7)*ro - 1.D0/2.D
     &    0*s(3,7)*s(6,3)**2*s(j1,4)*t1**(-1)*ro + 1.D0/2.D0*s(3,7)*s(6
     &    ,3)**2*s(j1,4)*ro - 1.D0/2.D0*s(3,7)*s(6,3)**2*s(j2,4)*
     &    t1**(-1)*ro + 1.D0/2.D0*s(3,7)*s(6,3)**2*s(j2,4)*ro - 1.D0/2.D
     &    0*s(3,7)*s(6,3)*s(6,5)*s(j1,4)*t1**(-1)*ro + s(3,7)*s(6,3)*s(
     &    6,5)*s(j1,4)*ro )
      ttbdkn = ttbdkn + s12**(-1) * ( s(3,7)*s(6,3)*s(6,5)*s(j2,4)*ro
     &     + 1.D0/2.D0*s(3,7)*s(6,5)**2*s(j1,4)*ro - 1.D0/2.D0*s(3,7)*
     &    s(6,5)**2*s(j2,4)*t2**(-1)*ro + 1.D0/2.D0*s(3,7)*s(6,5)**2*s(
     &    j2,4)*ro - 1.D0/2.D0*s(6,3)**2*s(j1,4)*s(j1,7)*t1**(-1)*ro + 
     &    1.D0/2.D0*s(6,3)**2*s(j1,4)*s(j1,7)*ro - 1.D0/2.D0*s(6,3)**2*
     &    s(j1,4)*s(j2,7)*t1**(-1)*ro + 1.D0/2.D0*s(6,3)**2*s(j1,4)*s(
     &    j2,7)*ro - 1.D0/2.D0*s(6,3)**2*s(j1,7)*s(j2,4)*t1**(-1)*ro + 
     &    1.D0/2.D0*s(6,3)**2*s(j1,7)*s(j2,4)*ro - 1.D0/2.D0*s(6,3)**2*
     &    s(j2,4)*s(j2,7)*t1**(-1)*ro + 1.D0/2.D0*s(6,3)**2*s(j2,4)*s(
     &    j2,7)*ro - 1.D0/2.D0*s(6,3)*s(6,5)*s(j1,4)*s(j1,7)*t1**(-1)*
     &    ro + s(6,3)*s(6,5)*s(j1,4)*s(j1,7)*ro - 1.D0/2.D0*s(6,3)*s(6,
     &    5)*s(j1,4)*s(j2,7)*t1**(-1)*ro + s(6,3)*s(6,5)*s(j1,4)*s(j2,7
     &    )*ro + s(6,3)*s(6,5)*s(j1,7)*s(j2,4)*ro + s(6,3)*s(6,5)*s(j2,
     &    4)*s(j2,7)*ro + 1.D0/2.D0*s(6,5)**2*s(j1,4)*s(j1,7)*ro + 1.D0/
     &    2.D0*s(6,5)**2*s(j1,4)*s(j2,7)*ro - 1.D0/2.D0*s(6,5)**2*s(j1,
     &    7)*s(j2,4)*t2**(-1)*ro )
      ttbdkn = ttbdkn + s12**(-1) * ( 1.D0/2.D0*s(6,5)**2*s(j1,7)*s(j2,
     &    4)*ro - 1.D0/2.D0*s(6,5)**2*s(j2,4)*s(j2,7)*t2**(-1)*ro + 1.D0
     &    /2.D0*s(6,5)**2*s(j2,4)*s(j2,7)*ro )
      ttbdkn = ttbdkn + nDn * ( 1.D0/2.D0*s(j1,4)*s(j1,7)*t1**(-1)*ro
     &     - 1.D0/4.D0*s(j1,4)*s(j1,7)*t1**(-1)*ro**2 - 1.D0/4.D0*s(j1,
     &    4)*s(j1,7)*t2**(-1)*ro**2 + 1.D0/2.D0*s(j1,4)*s(j2,7)*
     &    t1**(-1)*ro + s(j1,4)*s(j2,7)*t1*ro - s(j1,4)*s(j2,7)*ro - s(
     &    j1,7)*s(j2,4)*t1*ro + s(1,2)*s(4,7)*t1*ro - s(1,2)*s(4,7)*
     &    t1**2*ro + 2.D0*s(3,4)*s(j1,7) - 4.D0*s(3,4)*s(j1,7)*t1 + 4.D0
     &    *s(3,4)*s(j1,7)*t1**2 - 1.D0/2.D0*s(3,4)*s(j1,7)*t2**(-1)*ro
     &     + 2.D0*s(3,4)*s(j2,7) - 4.D0*s(3,4)*s(j2,7)*t1 + 4.D0*s(3,4)
     &    *s(j2,7)*t1**2 + 2.D0*s(3,4)*s(3,7) - 4.D0*s(3,4)*s(3,7)*t1
     &     + 4.D0*s(3,4)*s(3,7)*t1**2 + 1.D0/2.D0*s(3,7)*s(j1,4)*
     &    t1**(-1)*ro )
      ttbdkn = ttbdkn - 1.D0/2.D0*s(3,4)*s(6,3)*s(6,7)*ro + 1.D0/2.D0*
     &    s(3,4)*s(6,5)*s(6,7)*t2**(-1)*ro - 1.D0/2.D0*s(3,4)*s(6,5)*s(
     &    6,7)*ro - 1.D0/2.D0*s(3,7)*s(6,3)*s(6,4)*t1**(-1)*ro + 1.D0/2.
     &    D0*s(3,7)*s(6,3)*s(6,4)*ro + 1.D0/2.D0*s(3,7)*s(6,4)*s(6,5)*
     &    ro + 1.D0/4.D0*s(4,7)*s(6,3)**2*t1**(-1)*ro**2 - 1.D0/4.D0*s(
     &    4,7)*s(6,3)**2*ro**2 - 1.D0/2.D0*s(4,7)*s(6,3)*s(6,5)*ro**2
     &     + 1.D0/4.D0*s(4,7)*s(6,5)**2*t2**(-1)*ro**2 - 1.D0/4.D0*s(4,
     &    7)*s(6,5)**2*ro**2 - 1.D0/2.D0*s(6,3)*s(6,4)*s(j1,7)*t1**(-1)
     &    *ro + 1.D0/4.D0*s(6,3)*s(6,4)*s(j1,7)*t1**(-1)*ro**2 + 1.D0/2.
     &    D0*s(6,3)*s(6,4)*s(j1,7)*ro - 1.D0/2.D0*s(6,3)*s(6,4)*s(j2,7)
     &    *t1**(-1)*ro + 1.D0/2.D0*s(6,3)*s(6,4)*s(j2,7)*ro - 1.D0/4.D0
     &    *s(6,3)*s(6,7)*s(j1,4)*t1**(-1)*ro**2 - 1.D0/4.D0*s(6,4)*s(6,
     &    5)*s(j1,7)*t2**(-1)*ro**2 + 1.D0/2.D0*s(6,4)*s(6,5)*s(j1,7)*
     &    ro + 1.D0/2.D0*s(6,4)*s(6,5)*s(j2,7)*ro + 1.D0/4.D0*s(6,5)*s(
     &    6,7)*s(j1,4)*t2**(-1)*ro**2

      ttbdkn = (V/(t1*t2)-2d0*xn**2)*ttbdkn
      return
      end
