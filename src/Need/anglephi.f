      subroutine anglephi(p1,p2,p4,p5,p6,p7,cosphi)
      implicit none
      integer j
      double precision cosphi

      double precision p1(4),p2(4),p4(4),p5(4),p6(4),p7(4),pwh(4)
      double precision p1r(4),p4r(4),p5r(4),p6r(4),p7r(4),pwr(4),phr(4),
     . p4n(4)
      double precision p4cm(4),p5cm(4),p6cm(4),p7cm(4),pw(4),ph(4),
     . pwcm(4),phcm(4)
      double precision pwnew(4),px(4),p1cm(4),p2cm(4)
      do j=1,4
      pwh(j)=p4(j)+p5(j)+p6(j)+p7(j)
      pw(j)=p6(j)+p7(j)
      ph(j)=p4(j)+p5(j)
      enddo


c----calculate l momentum in centre of mass
c      call boosta(pwh,p1,p1cm)
c      call boosta(pwh,p2,p2cm)
c      call boosta(pwh,p4,p4cm)
      call boosta(pwh,p5,p5cm)
c      call boosta(pwh,p6,p6cm)
      call boosta(pwh,p7,p7cm)
      call boosta(pwh,pw,pwcm)
c      call boosta(pwh,ph,phcm)


c      call rotate(pwcm,p1cm,p1r)
c      call rotate(pwcm,p4cm,p4r)
c      call rotate(pwcm,p4cm,p4r)
      call rotate(pwcm,p5cm,p5r)
      call rotate(pwcm,p7cm,p7r)
c      call rotate(pwcm,p6cm,p6r)

      cosphi=(p5r(1)*p7r(1)+p5r(2)*p7r(2))
     . /(sqrt(p5r(1)**2+p5r(2)**2)*sqrt(p7r(1)**2+p7r(2)**2))

      return
      end

      subroutine rotate(pw,pl,pr)
      implicit none
      double precision costh,sinth,temp3
c ---Rotate pl to frame in which pw is along the 3 axis (with unchanged sign)
c----output is 4-vector pr
      double precision pw(4),pl(4),pr(4),pt(3)

      pr(4)=pl(4)

      costh=pw(3)/sqrt(pw(2)**2+pw(3)**2)
      sinth=pw(2)/sqrt(pw(2)**2+pw(3)**2)

      pt(1)=pl(1)
      pr(2)=+costh*pl(2)-sinth*pl(3)
      pt(3)=+sinth*pl(2)+costh*pl(3)
      temp3=sinth*pw(2)+costh*pw(3)
      costh=temp3/sqrt(pw(1)**2+temp3**2)
      sinth=pw(1)/sqrt(pw(1)**2+temp3**2)

      pr(1)=costh*pt(1)-sinth*pt(3)
      pr(3)=sinth*pt(1)+costh*pt(3)

      return
      end

