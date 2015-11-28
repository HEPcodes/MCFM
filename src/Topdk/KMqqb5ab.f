      subroutine KMqqb5ab(mt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm,Ma,Mb)
      implicit none
      include 'constants.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (4.12)


      double complex zp1(4),zp2(4),zp3(4),zp4(4),Ma,Mb
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4)
      double complex z1(4),z2(4),z3(4),z4(4)
      double complex string0,string1,string2,string3
      double complex qbqQbQ(0:6),sum(0:6)
      double precision mt
      data z1/(1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0)/
      data z2/(0d0,0d0),(1d0,0d0),(0d0,0d0),(0d0,0d0)/
      data z3/(0d0,0d0),(0d0,0d0),(1d0,0d0),(0d0,0d0)/
      data z4/(0d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)/
      integer j

      qbqQbQ(0)=
     & +string1(Vbl,z1,Ul)*string1(Ubm,z1,Vm)
     & -string1(Vbl,z2,Ul)*string1(Ubm,z2,Vm)
     & -string1(Vbl,z3,Ul)*string1(Ubm,z3,Vm)
     & -string1(Vbl,z4,Ul)*string1(Ubm,z4,Vm)

      qbqQbQ(1)=string1(Vbl,zp3,Ul)*string1(Ubm,zp1,Vm)

      qbqQbQ(2)=
     & +string3(Vbl,z1,zp3,z1,Ul)*string3(Ubm,z1,zp1,z1,Vm)
     & -string3(Vbl,z1,zp3,z2,Ul)*string3(Ubm,z2,zp1,z1,Vm)
     & -string3(Vbl,z1,zp3,z3,Ul)*string3(Ubm,z3,zp1,z1,Vm)
     & -string3(Vbl,z1,zp3,z4,Ul)*string3(Ubm,z4,zp1,z1,Vm)

     & -string3(Vbl,z2,zp3,z1,Ul)*string3(Ubm,z1,zp1,z2,Vm)
     & +string3(Vbl,z2,zp3,z2,Ul)*string3(Ubm,z2,zp1,z2,Vm)
     & +string3(Vbl,z2,zp3,z3,Ul)*string3(Ubm,z3,zp1,z2,Vm)
     & +string3(Vbl,z2,zp3,z4,Ul)*string3(Ubm,z4,zp1,z2,Vm)

     & -string3(Vbl,z3,zp3,z1,Ul)*string3(Ubm,z1,zp1,z3,Vm)
     & +string3(Vbl,z3,zp3,z2,Ul)*string3(Ubm,z2,zp1,z3,Vm)
     & +string3(Vbl,z3,zp3,z3,Ul)*string3(Ubm,z3,zp1,z3,Vm)
     & +string3(Vbl,z3,zp3,z4,Ul)*string3(Ubm,z4,zp1,z3,Vm)

     & -string3(Vbl,z4,zp3,z1,Ul)*string3(Ubm,z1,zp1,z4,Vm)
     & +string3(Vbl,z4,zp3,z2,Ul)*string3(Ubm,z2,zp1,z4,Vm)
     & +string3(Vbl,z4,zp3,z3,Ul)*string3(Ubm,z3,zp1,z4,Vm)
     & +string3(Vbl,z4,zp3,z4,Ul)*string3(Ubm,z4,zp1,z4,Vm)

      qbqQbQ(3)=
     &  + string3(Vbl,z1,z1,z1,Ul)*string3(Ubm,z1,z1,z1,Vm)
     &  - string3(Vbl,z1,z1,z2,Ul)*string3(Ubm,z2,z1,z1,Vm)
     &  - string3(Vbl,z1,z1,z3,Ul)*string3(Ubm,z3,z1,z1,Vm)
     &  - string3(Vbl,z1,z1,z4,Ul)*string3(Ubm,z4,z1,z1,Vm)
     &  - string3(Vbl,z1,z2,z1,Ul)*string3(Ubm,z1,z2,z1,Vm)
     &  + string3(Vbl,z1,z2,z2,Ul)*string3(Ubm,z2,z2,z1,Vm)
     &  + string3(Vbl,z1,z2,z3,Ul)*string3(Ubm,z3,z2,z1,Vm)
     &  + string3(Vbl,z1,z2,z4,Ul)*string3(Ubm,z4,z2,z1,Vm)
     &  - string3(Vbl,z1,z3,z1,Ul)*string3(Ubm,z1,z3,z1,Vm)
     &  + string3(Vbl,z1,z3,z2,Ul)*string3(Ubm,z2,z3,z1,Vm)
     &  + string3(Vbl,z1,z3,z3,Ul)*string3(Ubm,z3,z3,z1,Vm)
     &  + string3(Vbl,z1,z3,z4,Ul)*string3(Ubm,z4,z3,z1,Vm)
     &  - string3(Vbl,z1,z4,z1,Ul)*string3(Ubm,z1,z4,z1,Vm)
     &  + string3(Vbl,z1,z4,z2,Ul)*string3(Ubm,z2,z4,z1,Vm)
      qbqQbQ(3)= qbqQbQ(3)
     &  + string3(Vbl,z1,z4,z3,Ul)*string3(Ubm,z3,z4,z1,Vm)
     &  + string3(Vbl,z1,z4,z4,Ul)*string3(Ubm,z4,z4,z1,Vm)
     &  - string3(Vbl,z2,z1,z1,Ul)*string3(Ubm,z1,z1,z2,Vm)
     &  + string3(Vbl,z2,z1,z2,Ul)*string3(Ubm,z2,z1,z2,Vm)
     &  + string3(Vbl,z2,z1,z3,Ul)*string3(Ubm,z3,z1,z2,Vm)
     &  + string3(Vbl,z2,z1,z4,Ul)*string3(Ubm,z4,z1,z2,Vm)
     &  + string3(Vbl,z2,z2,z1,Ul)*string3(Ubm,z1,z2,z2,Vm)
     &  - string3(Vbl,z2,z2,z2,Ul)*string3(Ubm,z2,z2,z2,Vm)
     &  - string3(Vbl,z2,z2,z3,Ul)*string3(Ubm,z3,z2,z2,Vm)
     &  - string3(Vbl,z2,z2,z4,Ul)*string3(Ubm,z4,z2,z2,Vm)
     &  + string3(Vbl,z2,z3,z1,Ul)*string3(Ubm,z1,z3,z2,Vm)
     &  - string3(Vbl,z2,z3,z2,Ul)*string3(Ubm,z2,z3,z2,Vm)
     &  - string3(Vbl,z2,z3,z3,Ul)*string3(Ubm,z3,z3,z2,Vm)
     &  - string3(Vbl,z2,z3,z4,Ul)*string3(Ubm,z4,z3,z2,Vm)
     &  + string3(Vbl,z2,z4,z1,Ul)*string3(Ubm,z1,z4,z2,Vm)
      qbqQbQ(3)= qbqQbQ(3)
     &  - string3(Vbl,z2,z4,z2,Ul)*string3(Ubm,z2,z4,z2,Vm)
     &  - string3(Vbl,z2,z4,z3,Ul)*string3(Ubm,z3,z4,z2,Vm)
     &  - string3(Vbl,z2,z4,z4,Ul)*string3(Ubm,z4,z4,z2,Vm)
     &  - string3(Vbl,z3,z1,z1,Ul)*string3(Ubm,z1,z1,z3,Vm)
     &  + string3(Vbl,z3,z1,z2,Ul)*string3(Ubm,z2,z1,z3,Vm)
     &  + string3(Vbl,z3,z1,z3,Ul)*string3(Ubm,z3,z1,z3,Vm)
     &  + string3(Vbl,z3,z1,z4,Ul)*string3(Ubm,z4,z1,z3,Vm)
     &  + string3(Vbl,z3,z2,z1,Ul)*string3(Ubm,z1,z2,z3,Vm)
     &  - string3(Vbl,z3,z2,z2,Ul)*string3(Ubm,z2,z2,z3,Vm)
     &  - string3(Vbl,z3,z2,z3,Ul)*string3(Ubm,z3,z2,z3,Vm)
     &  - string3(Vbl,z3,z2,z4,Ul)*string3(Ubm,z4,z2,z3,Vm)
     &  + string3(Vbl,z3,z3,z1,Ul)*string3(Ubm,z1,z3,z3,Vm)
     &  - string3(Vbl,z3,z3,z2,Ul)*string3(Ubm,z2,z3,z3,Vm)
     &  - string3(Vbl,z3,z3,z3,Ul)*string3(Ubm,z3,z3,z3,Vm)
     &  - string3(Vbl,z3,z3,z4,Ul)*string3(Ubm,z4,z3,z3,Vm)
      qbqQbQ(3)= qbqQbQ(3)
     &  + string3(Vbl,z3,z4,z1,Ul)*string3(Ubm,z1,z4,z3,Vm)
     &  - string3(Vbl,z3,z4,z2,Ul)*string3(Ubm,z2,z4,z3,Vm)
     &  - string3(Vbl,z3,z4,z3,Ul)*string3(Ubm,z3,z4,z3,Vm)
     &  - string3(Vbl,z3,z4,z4,Ul)*string3(Ubm,z4,z4,z3,Vm)
     &  - string3(Vbl,z4,z1,z1,Ul)*string3(Ubm,z1,z1,z4,Vm)
     &  + string3(Vbl,z4,z1,z2,Ul)*string3(Ubm,z2,z1,z4,Vm)
     &  + string3(Vbl,z4,z1,z3,Ul)*string3(Ubm,z3,z1,z4,Vm)
     &  + string3(Vbl,z4,z1,z4,Ul)*string3(Ubm,z4,z1,z4,Vm)
     &  + string3(Vbl,z4,z2,z1,Ul)*string3(Ubm,z1,z2,z4,Vm)
     &  - string3(Vbl,z4,z2,z2,Ul)*string3(Ubm,z2,z2,z4,Vm)
     &  - string3(Vbl,z4,z2,z3,Ul)*string3(Ubm,z3,z2,z4,Vm)
     &  - string3(Vbl,z4,z2,z4,Ul)*string3(Ubm,z4,z2,z4,Vm)
     &  + string3(Vbl,z4,z3,z1,Ul)*string3(Ubm,z1,z3,z4,Vm)
     &  - string3(Vbl,z4,z3,z2,Ul)*string3(Ubm,z2,z3,z4,Vm)
     &  - string3(Vbl,z4,z3,z3,Ul)*string3(Ubm,z3,z3,z4,Vm)
      qbqQbQ(3)= qbqQbQ(3)
     &  - string3(Vbl,z4,z3,z4,Ul)*string3(Ubm,z4,z3,z4,Vm)
     &  + string3(Vbl,z4,z4,z1,Ul)*string3(Ubm,z1,z4,z4,Vm)
     &  - string3(Vbl,z4,z4,z2,Ul)*string3(Ubm,z2,z4,z4,Vm)
     &  - string3(Vbl,z4,z4,z3,Ul)*string3(Ubm,z3,z4,z4,Vm)
     &  - string3(Vbl,z4,z4,z4,Ul)*string3(Ubm,z4,z4,z4,Vm)


      qbqQbQ(4)=mt*string1(Vbl,zp3,Ul)*string0(Ubm,Vm)

      qbqQbQ(5)=mt*(
     & +string1(Vbl,z1,Ul)*string2(Ubm,z1,zp1,Vm)
     & -string1(Vbl,z2,Ul)*string2(Ubm,z2,zp1,Vm)
     & -string1(Vbl,z3,Ul)*string2(Ubm,z3,zp1,Vm)
     & -string1(Vbl,z4,Ul)*string2(Ubm,z4,zp1,Vm))

      qbqQbQ(6)=mt*(
     & +string3(Vbl,z1,zp3,z1,Ul)*string2(Ubm,z1,z1,Vm)
     & -string3(Vbl,z2,zp3,z1,Ul)*string2(Ubm,z1,z2,Vm)
     & -string3(Vbl,z3,zp3,z1,Ul)*string2(Ubm,z1,z3,Vm)
     & -string3(Vbl,z4,zp3,z1,Ul)*string2(Ubm,z1,z4,Vm)

     & -string3(Vbl,z1,zp3,z2,Ul)*string2(Ubm,z2,z1,Vm)
     & +string3(Vbl,z2,zp3,z2,Ul)*string2(Ubm,z2,z2,Vm)
     & +string3(Vbl,z3,zp3,z2,Ul)*string2(Ubm,z2,z3,Vm)
     & +string3(Vbl,z4,zp3,z2,Ul)*string2(Ubm,z2,z4,Vm)

     & -string3(Vbl,z1,zp3,z3,Ul)*string2(Ubm,z3,z1,Vm)
     & +string3(Vbl,z2,zp3,z3,Ul)*string2(Ubm,z3,z2,Vm)
     & +string3(Vbl,z3,zp3,z3,Ul)*string2(Ubm,z3,z3,Vm)
     & +string3(Vbl,z4,zp3,z3,Ul)*string2(Ubm,z3,z4,Vm)

     & -string3(Vbl,z1,zp3,z4,Ul)*string2(Ubm,z4,z1,Vm)
     & +string3(Vbl,z2,zp3,z4,Ul)*string2(Ubm,z4,z2,Vm)
     & +string3(Vbl,z3,zp3,z4,Ul)*string2(Ubm,z4,z3,Vm)
     & +string3(Vbl,z4,zp3,z4,Ul)*string2(Ubm,z4,z4,Vm))

c--- Switch signs of contributions (1), (2) and (5) (which all
c---  contain p1-slash) due to difference in incoming/outgoing momenta
      qbqQbQ(1)=-qbqQbQ(1)
      qbqQbQ(2)=-qbqQbQ(2)
      qbqQbQ(5)=-qbqQbQ(5)
 
      call coeffboxa(mt,zp1,zp2,zp3,zp4,sum)
      Ma=im*sum(0)*qbqQbQ(0)
      do j=1,6
      Ma=Ma+im*sum(j)*qbqQbQ(j)
      enddo      

      call coeffboxb(mt,zp1,zp2,zp3,zp4,sum)
      Mb=im*sum(0)*qbqQbQ(0)
      do j=1,6
      Mb=Mb+im*sum(j)*qbqQbQ(j)
      enddo    

c--- now multiply by overall color factors
c---  Eq. (4.13) for Boxa:
c---   (Ta Tb) (Tb Ta) -> (Cf-1/2/Nc)(Ta)(Ta) + delta term that vanishes
      Ma=Ma*(Cf-1d0/2d0/xn)
c---  Eq. (4.15) for Boxb:
c---   (Ta Tb) (Ta Tb) -> (-1/Nc)(Ta)(Ta) + delta term that vanishes
      Mb=Mb*(-1d0/xn)
        
      return
      end
