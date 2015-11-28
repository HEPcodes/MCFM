      subroutine wconstruct(p,costh_h,cosnew1,cosnew2,cosnew3,cosphi)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'nn.f'
      logical exact_w_mom
      double precision ptn,pte,ple,pln,Ee,En,ptwsq,mtsq,a,b,cosnew1,
     & cosnew2
      double precision cosnew3,nh_cmw1(4),nh_cmw2(4),nw_cmh1(4),
     & nw_cmh2(4)
      double precision discr,plntrue,p(mxpart,4),pw(4),pbb(4),pwbb(4),
     & p1(4),p2(4),
     & p1cm(4),pbbcm(4),pcscm(4),pwcm(4),plnp,plnm,costh_h,
     & dotpr,lambda,e1,e2,s2,shat,
     & p4cm(4),p5cm(4),ncm(4),p3(4),p4(4),p5(4),p6(4),xx,
     & pdiff(4),pdiffcm(4),cosphi
      integer j

      exact_w_mom=.true.

      do j=1,4
      pbb(j)=p(5,j)+p(6,j)
      pw(j)=p(3,j)+p(4,j)
      pwbb(j)=pw(j)+pbb(j)
      p1(j)=-p(1,j)
      p2(j)=-p(2,j)
      p3(j)=+p(3,j)
      p4(j)=+p(4,j)
      p5(j)=+p(5,j)
      p6(j)=+p(6,j)
      if (p(6,j) .eq. 0d0) then
      write(*,*) 1,j,p(1,j)
      write(*,*) 2,j,p(2,j)
      write(*,*) 3,j,p(3,j)
      write(*,*) 4,j,p(4,j)
      write(*,*) 5,j,p(5,j)
      write(*,*) 6,j,p(6,j)
      write(*,*) 7,j,p(7,j)
      pause
      endif
      enddo

      if (exact_w_mom) go to 20
      pte=sqrt(p(4,1)**2+p(4,2)**2)
      ptn=sqrt(p(3,1)**2+p(3,2)**2)
      ptwsq=(pw(1)**2+pw(2)**2)
      Ee=p(4,4)
      ple=p(4,3)
      plntrue=p(3,3)
      mtsq=(pte+ptn)**2-ptwsq

c  + 1/2*Mw^2*MT^2 - 1/4*Mw^4 - 1/4*MT^4 - pte*ptn*Mw^2 + pte*ptn*MT^2
c   + ptn^2*ple^2
c 
c + pln
c  * ( - 2*pte*ptn*ple - ple*Mw^2 + ple*MT^2 )
c 
c + pln^2
c  * ( pte^2 ) + 0.

c      write(6,*) 'wmass**2,mtsq',wmass**2,mtsq
      if (wmass**2 .gt. mtsq) then
      discr=Ee*sqrt((wmass**2-mtsq)*(wmass**2-mtsq+4*pte*ptn))
      a=pte**2
      b=-(2*pte*ptn+wmass**2-mtsq)*ple
c      c=-0.25d0*(wmass**2-mtsq+2*pte*ptn)**2
c     # +ptn**2*Ee**2
c      write(6,*) 'discr',sqrt(b**2-4*a*c),discr
      else
      discr=0
      a=pte**2
      b=-2*pte*ptn*ple
      endif

c      wsq=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2
c      write(6,*) 'cal:mw',wsq,plntrue
      plnp=(-b+discr)/(two*a)
      plnm=(-b-discr)/(two*a)

c---substitute calculated value of pln  
      pln=plnm
      En=sqrt(pln**2+ptn**2)
      pw(3)=pln+p(4,3)
      pw(4)=En+p(4,4)
c      wsq=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2
c      write(6,*) 'meas:mwm',wsq,pln

      pln=plnp
      En=sqrt(pln**2+ptn**2)
      pw(4)=En+p(4,4)
      pw(3)=pln+p(4,3)
c      wsq=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2

c      write(6,*) 'meas:mwp',wsq,pln
c      pause
      if (abs(plntrue-plnp) .lt. abs(plntrue-plnm)) then
      ntrue=ntrue+1
      else
      nfalse=nfalse+1
      endif
      ncount=nfalse+ntrue
      if (mod(ncount,10000) .eq. 0) 
     . write(6,*) dfloat(ntrue)/dfloat(ncount)

      pwbb(4)=pw(4)+pbb(4)
      pwbb(3)=pw(3)+pbb(3)
 
 20   continue
      

c---boost to center of mass
      s2=dotpr(pwbb,pwbb)
      shat=2d0*dotpr(p1,p2)
c energy of particle one in the wbb centre of mass
      E1=dotpr(p1,pwbb)/sqrt(s2)
c energy of particle two in the wbb centre of mass
      E2=dotpr(p2,pwbb)/sqrt(s2)

      do j=1,4
      pcscm(j)=E2*p1(j)-E1*p2(j)
      enddo
 
c----calculate w momentum in centre of mass
      call boosta(pwbb,pw,pwcm)
c----calculate bb momentum in centre of mass
      call boosta(pwbb,pbb,pbbcm)
      call boosta(pwbb,p1,p1cm)
c---momentum of higgs in cm
      lambda=sqrt(pbbcm(1)**2+pbbcm(2)**2+pbbcm(3)**2)
c----beta of W in cm
c      betaw=sqrt(pwcm(1)**2+pwcm(2)**2+pwcm(3)**2)/pwcm(4)
c---Collins-Soper angle
      costh_h=-dotpr(pcscm,pbb)/(lambda*sqrt(E1*E2*shat))

c---construct basis vectors
      call angle(p1,p2,pw,pbb,costh_h,
     . nh_cmw1,nh_cmw2,nw_cmh1,nw_cmh2)


c      write(6,*) 'nw_cmh',nw_cmh(4),nw_cmh(1),nw_cmh(2),nw_cmh(3)
c      write(6,*) 'nhcm in wc',nhcm(4),nhcm(1),nhcm(2),nhcm(3)
c      write(6,*) 'p7cm in wc',p7cm(4),p7cm(1),p7cm(2),p7cm(3)
c      pause
      call boosta(pw,nh_cmw1,ncm)
      call boosta(pw,p4,p4cm)
      xx=sqrt(p4cm(1)**2+p4cm(2)**2+p4cm(3)**2)
      cosnew1=-dotpr(p4cm,ncm)/xx

      call boosta(pw,nh_cmw2,ncm)
      xx=sqrt(p4cm(1)**2+p4cm(2)**2+p4cm(3)**2)
      cosnew2=-dotpr(p4cm,ncm)/xx

      call boosta(pbb,nw_cmh1,ncm)
      call boosta(pbb,p5,p5cm)
      xx=sqrt(p5cm(1)**2+p5cm(2)**2+p5cm(3)**2)
      cosnew3=-dotpr(p5cm,ncm)/xx

c----Now a new attempt
      do j=1,4
      pdiff(j)=p5(j)-p6(j)
      enddo
c express pdiff in bb centre of mass, just to check
      call boosta(pbb,pdiff,pdiffcm)
      xx=sqrt(pdiffcm(1)**2+pdiffcm(2)**2+pdiffcm(3)**2)
      do j=1,4
      pdiffcm(j)=pdiffcm(j)/xx
      enddo
      call boosta(pbb,p5,p5cm)
      cosnew1=dotpr(p5cm,pdiffcm)/sqrt(p5cm(1)**2+p5cm(2)**2+p5cm(3)**2)
c      write(6,*) 'cosnew1',cosnew1
c express pdiff in w centre of mass.
      call boosta(pw,pdiff,pdiffcm)
      xx=sqrt(pdiffcm(1)**2+pdiffcm(2)**2+pdiffcm(3)**2)
      pdiffcm(4)=0d0
      do j=1,3
      pdiffcm(j)=pdiffcm(j)/xx
      enddo
      call boosta(pw,p4,p4cm)

      cosnew2=dotpr(p4cm,pdiffcm)/sqrt(p4cm(1)**2+p4cm(2)**2+p4cm(3)**2)

      call anglephi(p1,p2,p5,p6,p3,p4,cosphi)

c      costh_h=(pbbcm(1)*p1cm(1)+pbbcm(2)*p1cm(2)+pbbcm(3)*p1cm(3))
c     .   /lambda/sqrt(p1cm(1)**2+p1cm(2)**2+p1cm(3)**2)

      return
      end


 


