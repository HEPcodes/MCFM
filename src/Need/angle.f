      subroutine angle(p1,p2,pw,ph,costh,
     . nh_cmw1,nh_cmw2,nw_cmh1,nw_cmh2)
      implicit none
      double precision cosxi,s,MW,MH,betaw,costh
     .  ,nh_cmw1(4),nh_cmw2(4),nw_cmh1(4),nw_cmh2(4)
      double precision p1(4),p2(4),pw(4),ph(4),dotpr
      double precision pwDph,pwDpw,phDph


      
      pwDph=dotpr(pw,ph)
      pwDpw=dotpr(pw,pw)
      phDph=dotpr(ph,ph)
c---calculation of betaw
      s=phDph+pwDpw+2d0*pwDph
      MH=dsqrt(phDph)
      MW=dsqrt(pwDpw)
      
      betaw=dsqrt((s-(MW+MH)**2)*(s-(MW-MH)**2))/(s-MH**2+MW**2)
      cosxi=costh*dsqrt(1d0-betaw**2)/dsqrt(1d0-(betaw*costh)**2)

c      write(6,*) 'betaw in angle',betaw
c      write(6,*) 'cosxi',cosxi

      call vec(pw,ph,p1,cosxi,nh_cmw1,nh_cmw2)
      call vec(ph,pw,p1,cosxi,nw_cmh1,nw_cmh2)
 

c      write(6,*) 'nw_cmh angle',nw_cmh(4),nw_cmh(1),nw_cmh(2),nw_cmh(3)
c      call boosta(pw,n2,n2cm)
c      write(6,*) 'dotpr(n1cm,n1cm)',dotpr(n1cm,n1cm)
c      write(6,*) 'n1:',n1(4),n1(1),n1(2),n1(3)
c      write(6,*) 'n1cm:',n1cm(4),n1cm(1),n1cm(2),n1cm(3)
c      write(6,*) 'dotpr(n2cm,n2cm)',dotpr(n2cm,n2cm)
c      write(6,*) 'n2:',n2(4),n2(1),n2(2),n2(3)
c      write(6,*) 'n2cm:',n2cm(4),n2cm(1),n2cm(2),n2cm(3)
c      write(6,*) 'n1Dn2',dotpr(n1,n2)

      return
      end

      subroutine vec(pw,ph,p1,cosxi,nh_cmw1,nh_cmw2)
      implicit none
c---given pw,ph,p1 and cosxi contstruct nh_cmw
      double precision cosxi,sinxi,nh_cmw1(4),nh_cmw2(4)
      double precision dotpr,n1(4),n2(4),a,b,norm1,norm2,pw(4),
     . ph(4),p1(4)
      double precision pwDph,pwDpw,pwDp1,phDp1,phDph
      integer j

      sinxi=dsqrt(1d0-cosxi**2)

      pwDp1=dotpr(pw,p1)
      phDp1=dotpr(ph,p1)
      pwDph=dotpr(pw,ph)
      pwDpw=dotpr(pw,pw)
      phDph=dotpr(ph,ph)


      a=-(pwDph*phDp1-pwDp1*phDph)/(pwDph**2-phDph*pwDpw)
      b=-pwDp1/pwDph-pwDpw*a/pwDph

      do j=1,4
      n1(j)=ph(j)-pwDph/pwDpw*pw(j)
      n2(j)=a*pw(j)+b*ph(j)+p1(j)
      enddo

      norm1=dsqrt(-dotpr(n1,n1))
      norm2=dsqrt(-dotpr(n2,n2))
      do j=1,4
      n1(j)=n1(j)/norm1
      n2(j)=n2(j)/norm2
      nh_cmw1(j)=cosxi*n1(j)+sinxi*n2(j)
      nh_cmw2(j)=cosxi*n1(j)-sinxi*n2(j)
      enddo
      return
      end
