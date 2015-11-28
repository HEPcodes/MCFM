      subroutine KMampsqqqbQQb(q1,q2,q3,q4,xmt,res0,res1)
      implicit none
c----notation                                                       *
C     0--> qbar(q1)+q(q2)+Q(p3)+Qb(p4)
c---
c--- Computed using the helicity amplitudes from:
c---
c--- \bibitem{Korner:2002hy}
c--- J.~G.~Korner and Z.~Merebashvili,
c--- %``One-loop corrections to four-point functions with two external massive
c--- %fermions and two external massless partons,''
c--- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c--- [arXiv:hep-ph/0207054].
c--- %%CITATION = PHRVA,D66,054023;%%
      include 'constants.f'
      include 'scale.f'
      include 'epinv.f'
      include 'spinorsw.f'
      integer j,nu,j1,j2,j3,j4
      double precision q1(4),q2(4),q3(4),q4(4),xmt,res0,res1,nlf
      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex Ubm(4),Vm(4),Vbl(4),Ul(4)
      double complex lo1,Da,Db,Dc,Dd,De,Df,Dg
      double complex Bqq,KMqqb5c,KMqqb5d,KMqqb5e,KMqqb5f,KMqqb5g
      double complex ren,KMqqb
      parameter(nlf=5d0)
            
c--- set up common variable to indicate KM routines should be used
c--- to calculate spinor strings
      spinorsw='KM'      

c--- these routines expect energy in the 1st component and, as usual,
c--- we modifiy this solution by exchanging (px <--> pz) and changing 
C--- the sign of py, (to remain in right handed system) [c.f. spinor_w.f]
      do j=1,4
      if (j .eq. 1) nu=4
      if (j .eq. 2) nu=3
      if (j .eq. 3) nu=2
      if (j .eq. 4) nu=1
      zp1(j)=dcmplx(q1(nu))
      zp2(j)=dcmplx(q2(nu))
      zp3(j)=dcmplx(q3(nu))
      zp4(j)=dcmplx(q4(nu))
      if (j .eq. 3) then
      zp1(j)=-zp1(j)
      zp2(j)=-zp2(j)
      zp3(j)=-zp3(j)
      zp4(j)=-zp4(j)
      endif
      enddo

c--- initialize to zero
      res0=0d0
      res1=0d0

      do j1=-1,1,2
      call vspinor(zp4,xmt,j1,Vm)

      do j2=-1,1,2
      call ubarspinor(zp3,xmt,j2,Ubm)

      do j3=-1,1,2
      call vbar0spinor(zp2,j3,Vbl)

      j4=-j3            
      call u0spinor(zp1,j4,Ul)

      lo1=Bqq(zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      call KMqqb5ab(xmt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm,Da,Db)
      Dc=KMqqb5c(xmt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      Dd=KMqqb5d(xmt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      De=KMqqb5e(xmt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      Df=KMqqb5f(xmt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      Dg=KMqqb5g(xmt,zp1,zp2,zp3,zp4,Ul,Vbl,Ubm,Vm)
      
c--- put in overall wave function and charge renormalization by hand
      ren=lo1
      
      ren=ren*(
     & +2d0*((nlf/3d0-11d0/6d0*xn)*epinv+xn/6d0)
     & +2d0/3d0*(epinv+log(musq/xmt**2))
     & -(xn**2-1d0)/2d0/xn*(3d0*(epinv+log(musq/xmt**2))+5d0))
     
c--- sum all diagrams (ab)
      KMqqb=Da+Db+Dc+Dd+De+Df+Dg+ren
      
      res0=res0+cdabs(lo1)**2
      res1=res1+dreal(dconjg(lo1)*KMqqb)

      enddo
      enddo
      enddo

      res0=V/4d0*res0
      res1=V/4d0*res1
      
      return
      end
