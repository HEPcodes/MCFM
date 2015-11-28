      double complex function Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      implicit none

c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (4.1)

      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex Ul(4),Vbl(4),Ubm(4),Vm(4)
      double complex z1(4),z2(4),z3(4),z4(4),cdot,im,s,string1
      parameter(im=(0d0,1d0))
      data z1/(1d0,0d0),(0d0,0d0),(0d0,0d0),(0d0,0d0)/
      data z2/(0d0,0d0),(1d0,0d0),(0d0,0d0),(0d0,0d0)/
      data z3/(0d0,0d0),(0d0,0d0),(1d0,0d0),(0d0,0d0)/
      data z4/(0d0,0d0),(0d0,0d0),(0d0,0d0),(1d0,0d0)/
      s=2d0*cdot(zp1,zp2)
      Bqq=im*(
     & +string1(Vbl,z1,Ul)*string1(Ubm,z1,Vm)
     & -string1(Vbl,z2,Ul)*string1(Ubm,z2,Vm)
     & -string1(Vbl,z3,Ul)*string1(Ubm,z3,Vm)
     & -string1(Vbl,z4,Ul)*string1(Ubm,z4,Vm))/s

c      write(6,*) 'string0(Vbl,Ul)',string0(Vbl,Ul)
c      write(6,*) 'string1(Vbl,zp1,Ul)',string1(Vbl,zp1,Ul)
c      write(6,*) 'string1(Vbl,zp2,Ul)',string1(Vbl,zp2,Ul)
c      write(6,*) 'string1(Vbl,z1,Ul)',string1(Vbl,z1,Ul)
c      write(6,*) 'string1(Vbl,z2,Ul)',string1(Vbl,z2,Ul)
c      write(6,*) 'string1(Vbl,z3,Ul)',string1(Vbl,z3,Ul)
c      write(6,*) 'string1(Vbl,z4,Ul)',string1(Vbl,z4,Ul)
      return
      end
