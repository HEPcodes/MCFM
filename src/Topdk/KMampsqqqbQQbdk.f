      subroutine KMampsqqqbQQbdk(q1,q2,q3,q4,q5,q6,q7,q8,
     & xmt,xmw,res0,res1)
      implicit none
C     0--> qb(q1)+q(q2)+nu(q3)+e+(q4)+b(q5)+bbar(q6)+e-(q7)+nubar(q8)
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
      include 'epinv.f'
      include 'scale.f'
      include 'spinorsw.f'
      integer j,nu,j3,j4
      double precision q1(4),q2(4),q3(4),q4(4),q5(4),q6(4),q7(4),q8(4)
      double precision k1(4),k2(4),k3(4),k4(4),k5(4),k6(4),k7(4),k8(4),
     & xmt,xmw,res0,res1,nlf
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
      k1(j)=q1(nu)
      k2(j)=q2(nu)
      k3(j)=q3(nu)
      k4(j)=q4(nu)
      k5(j)=q5(nu)
      k6(j)=q6(nu)
      k7(j)=q7(nu)
      k8(j)=q8(nu)
      if (j .eq. 3) then
      k1(j)=-k1(j)
      k2(j)=-k2(j)
      k3(j)=-k3(j)
      k4(j)=-k4(j)
      k5(j)=-k5(j)
      k6(j)=-k6(j)
      k7(j)=-k7(j)
      k8(j)=-k8(j)
      endif
      enddo

C-----remember that diagrams are Qbar(k1)+Q(k2)+q(k3)+qbar(k4)
C     0--> qbar(q1)+q(q2)+Q(q345)+Qb(q678)
C     0--> qbar(k4)+q(k3)+Q(k2)+Qb(k1)
      do nu=1,4
      zp1(nu)=dcmplx(k1(nu))
      zp2(nu)=dcmplx(k2(nu))
      zp3(nu)=dcmplx(k3(nu)+k4(nu)+k5(nu))
      zp4(nu)=dcmplx(k6(nu)+k7(nu)+k8(nu))
      enddo

      call ubarWspinor(k3,k4,k5,xmt,xmw,Ubm)
      call uWspinor(k6,k7,k8,xmt,xmw,Vm)

c--- initialize to zero
      res0=0d0
      res1=0d0

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

      res0=V/4d0*res0
      res1=V/4d0*res1

      return
      end
