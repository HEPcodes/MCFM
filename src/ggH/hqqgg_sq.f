      subroutine hqqgg_sq(p1,p2,p3,p4,ampsq)
      implicit none
C     Taken from Kauffman
C     PRD 55 1997 (4009)
C     and checked with hep-ph/9903330 
      include 'constants.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,j,j1,j2,j3
      integer i1,i2,i3,i4,i5,i6,k1,k2
      double precision ampsq
      double precision sq(2,2,2),Trace4,Trace6,mssq,masq,m34sq,m43sq
      double precision n11,n12,n22,n33,n13,n23,d1,d2,d3,d4,d5,c1,c2
      double precision m34int
      parameter(c1=cf**2*xn,c2=-0.5d0*cf)

      Trace4(i1,i2,i3,i4)=
     . +s(i1,i2)*s(i3,i4)+s(i2,i3)*s(i1,i4)-s(i1,i3)*s(i2,i4)
      Trace6(i1,i2,i3,i4,i5,i6)=
     . +0.5d0*s(i1,i2)*Trace4(i3,i4,i5,i6)
     . -0.5d0*s(i1,i3)*Trace4(i2,i4,i5,i6)
     . +0.5d0*s(i1,i4)*Trace4(i2,i3,i5,i6)
     . -0.5d0*s(i1,i5)*Trace4(i2,i3,i4,i6)
     . +0.5d0*s(i1,i6)*Trace4(i2,i3,i4,i5)
      n11(i1,i2,i3,i4)=
     . 0.25d0*s(i1,i4)*s(i2,i4)
     . *(+Trace4(i2,i1,i3,i1)+Trace4(i2,i1,i3,i4)
     .   +Trace4(i2,i4,i3,i1)+Trace4(i2,i4,i3,i4))**2
      n22(i1,i2,i3,i4)=n11(i1,i2,i4,i3)
      n33(i1,i2,i3,i4)=
     . 0.25d0*s(i1,i2)*s(i2,i3)*s(i2,i4)*s(i3,i4)
     . *(+Trace4(i1,i3,i2,i3)+Trace4(i1,i3,i2,i4)
     .   +Trace4(i1,i4,i2,i3)+Trace4(i1,i4,i2,i4))**2
      n12(i1,i2,i3,i4)=
     . (+Trace6(i2,i1,i3,i1,i4,i1)+Trace6(i2,i1,i3,i1,i4,i3)
     .  +Trace6(i2,i4,i3,i1,i4,i1)+Trace6(i2,i4,i3,i1,i4,i3))
     .*(+Trace6(i2,i1,i3,i2,i4,i1)+Trace6(i2,i1,i3,i2,i4,i3)
     .  +Trace6(i2,i4,i3,i2,i4,i1)+Trace6(i2,i4,i3,i2,i4,i3))
     .  -Trace4(i1,i3,i2,i4)
     . *(s(i1,i2)*s(i1,i3)+s(i2,i4)*s(i3,i4)+Trace4(i1,i2,i4,i3))
     . *(s(i1,i2)*s(i1,i4)+s(i2,i3)*s(i3,i4)+Trace4(i1,i2,i3,i4))
      n13(i1,i2,i3,i4)=-s(i2,i4)*n12(i4,i2,i3,i1)
      n23(i1,i2,i3,i4)=n13(i1,i2,i4,i3)
      do j=1,2
      if (j.eq.1) then
         k1=p1
         k2=p2
      elseif (j.eq.2) then
         k1=p2
         k2=p1
      endif

      d1=s(k1,p4)*s(k2,p4)*(s(k1,p4)+s(k2,p4)+s(k1,k2))
      d2=s(k1,p3)*s(k2,p3)*(s(k1,p3)+s(k2,p3)+s(k1,k2))
      d3=s(k2,p4)*(s(k1,p4)+s(k2,p4)+s(k1,k2))
     . /(1d0/s(k1,k2)+0.5d0/s(k1,p4))
      d4=s(k2,p3)*(s(k1,p3)+s(k2,p3)+s(k1,k2))
     . /(1d0/s(k1,k2)+0.5d0/s(k1,p3))
      d5=s(k1,k2)*s(k2,p3)*s(k2,p4)*s(p3,p4)

      mssq=0.25d0*(n11(k1,k2,p3,p4)/d1**2+n22(k1,k2,p3,p4)/d2**2
     . + n12(k1,k2,p3,p4)/(d1*d2))

      masq=
     . +n11(k1,k2,p3,p4)/d3**2
     . +n22(k1,k2,p3,p4)/d4**2
     . +n33(k1,k2,p3,p4)/d5**2
     . -n12(k1,k2,p3,p4)/(d3*d4)
     . +n13(k1,k2,p3,p4)/(d3*d5)
     . +n23(k1,k2,p3,p4)/(d4*d5)


      m34sq=c1*(
     . +s(k1,p3)**3/(s(k1,k2)*s(k1,p4)*s(p3,p4))
     . +s(k2,p4)**3/(s(k1,k2)*s(k2,p3)*s(p3,p4))
     . +1d0/(s(k1,p4)*s(k2,p3)*(s(k1,k2)*s(p3,p4))**2)
     . *(-Trace4(k1,k2,p4,p3)**2*Trace4(k1,p3,k2,p4)
     . -s(k1,p3)*s(k2,p4)*Trace4(k1,k2,p3,p4)*Trace4(k1,k2,p4,p3)
     . +s(k1,k2)*s(k1,p3)*s(k2,p4)*s(p3,p4)*Trace4(k1,p3,k2,p4)))

      m43sq=c1*(
     . +s(k1,p3)**2*s(k2,p3)/(s(k1,k2)*s(k2,p4)*s(p3,p4))
     . +s(k1,p4)*s(k2,p4)**2/(s(k1,k2)*s(k1,p3)*s(p3,p4))
     . +(Trace4(k1,k2,p4,p3)*Trace4(k1,k2,p3,p4)
     .  +s(k1,k2)*s(p3,p4)*Trace4(k1,p3,k2,p4))
     . /(s(k1,k2)*s(p3,p4))**2)

      m34int=c2*(-Trace4(k1,p3,k2,p4)/(s(k1,k2)*s(p3,p4))
     . *(s(k1,p3)**2/s(k1,p4)/s(k2,p4)+s(k2,p4)**2/s(k1,p3)/s(k2,p3))
     . +2d0/(s(k1,k2)*s(p3,p4))**2
     .*(Trace4(k1,k2,p4,p3)**2-2d0*s(k1,k2)*s(k1,p3)*s(k2,p4)*s(p3,p4)))

c      m34inta=c2*(+2d0
c     . +s(k1,p3)**2/s(k1,p4)/s(k2,p4)+s(k2,p4)**2/s(k1,p3)/s(k2,p3)
c     .  + 2d0/s(k1,k2)**2/s(p3,p4)**2
c     . *(s(k1,p3)*s(k2,p4)-s(k1,p4)*s(k2,p3))**2
c     .  -((s(k2,p3)*s(k1,p4)+s(k1,p3)*s(k2,p4))
c     . *(s(k1,p3)**2/s(k1,p4)/s(k2,p4)+s(k2,p4)**2/s(k1,p3)/s(k2,p3))
c     .  +4d0*s(k1,p4)*s(k2,p3))/s(k1,k2)/s(p3,p4))

      if (j .eq. 1) then
      sq(2,2,2) = 2d0*(C1*(mssq+masq)+C2*(mssq-masq))
      sq(2,2,1) = m34sq+m43sq+m34int
      sq(1,1,1)=sq(2,2,2)
      sq(1,1,2)=sq(2,2,1)
      elseif (j .eq. 2) then
      sq(1,2,2) = 2d0*(C1*(mssq+masq)+C2*(mssq-masq))
      sq(1,2,1) = m34sq+m43sq+m34int
      sq(2,1,1)=sq(1,2,2)
      sq(2,1,2)=sq(1,2,1)
      endif
      enddo
      
      ampsq=0d0
      do j1=1,2
      do j2=1,2
      do j3=1,2
      ampsq=ampsq+sq(j1,j2,j3)
      enddo
      enddo
      enddo

      return
      end
