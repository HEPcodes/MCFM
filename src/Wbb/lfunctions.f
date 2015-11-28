************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
      double complex function L0(x,y)
      implicit none
      include 'constants.f'
      double complex Lnrat
      double precision x,y
      L0=Lnrat(x,y)/(one-x/y)
      return
      end

      double complex function L1(x,y)
      implicit none
      include 'constants.f'
      double precision x,y,r
      double complex l0
      r=x/y
      L1=(L0(x,y)+one)/(one-r)
      return
      end

      double complex function Ls0(x1,y1,x2,y2)
      implicit none
      include 'constants.f'
      double precision x1,x2,y1,y2,r1,r2
      double complex Lsm1
      r1=x1/y1
      r2=x2/y2
      Ls0=Lsm1(x1,y1,x2,y2)/(one-r1-r2)
      return
      end

      double complex function Ls1(x1,y1,x2,y2)
      implicit none
      include 'constants.f'
      double precision x1,x2,y1,y2,r1,r2
      double complex Ls0,L0
      r1=x1/y1
      r2=x2/y2
      Ls1=(Ls0(x1,y1,x2,y2)+L0(x1,y1)+L0(x2,y2))/(one-r1-r2)
      return
      end



      double complex function Lsm1(x1,y1,x2,y2)
      implicit none
      include 'constants.f'
      double precision x1,x2,y1,y2,r1,r2,omr1,omr2,ddilog
      double complex dilog1,dilog2,Lnrat
      r1=x1/y1
      r2=x2/y2
      omr1=one-r1
      omr2=one-r2
      if (omr1 .gt. one) then 
          dilog1=-ddilog(r1)+pisqo6-Lnrat(x1,y1)*log(omr1)
      else
          dilog1=dcmplx(ddilog(omr1))
      endif
      if (omr2 .gt. one) then 
          dilog2=-ddilog(r2)+pisqo6-Lnrat(x2,y2)*log(omr2)
      else
          dilog2=dcmplx(ddilog(omr2))
      endif
      lsm1=dilog1+dilog2+Lnrat(x1,y1)*Lnrat(x2,y2)-pisqo6
      return
      end

      double complex function Lsm1_2mh(s,t,m1sq,m2sq)
      implicit none
      include 'constants.f'
      double precision s,t,m1sq,m2sq
      double complex lsm1_2mht,I3m
      Lsm1_2mh=Lsm1_2mht(s,t,m1sq,m2sq)
     & +(half*(s-m1sq-m2sq)+m1sq*m2sq/t)*I3m(s,m1sq,m2sq)
      return
      end

      double complex function Lsm1_2mht(s,t,m1sq,m2sq)
      implicit none
      include 'constants.f'
      double precision s,t,m1sq,m2sq,ddilog,r1,r2,omr1,omr2
      double complex Lnrat,dilog1,dilog2
      r1=m1sq/t
      r2=m2sq/t
      omr1=one-r1
      omr2=one-r2

      if (omr1 .gt. one) then 
      dilog1=-ddilog(r1)+pisqo6-Lnrat(-m1sq,-t)*log(omr1)
      else
      dilog1=dcmplx(ddilog(omr1))
      endif
      if (omr2 .gt. one) then 
      dilog2=-ddilog(r2)+pisqo6-Lnrat(-m2sq,-t)*log(omr2)
      else
      dilog2=dcmplx(ddilog(omr2))
      endif
      lsm1_2mht=-dilog1-dilog2
     & +half*(Lnrat(-s,-m1sq)*Lnrat(-s,-m2sq)-Lnrat(-s,-t)**2)
      return
      end


      double complex function Lsm1_2me(s,t,m1sq,m3sq)
      implicit none
      include 'constants.f'
      integer j
      double precision s,t,m1sq,m3sq,ddilog,arg(5),omarg(5)
      double complex Lnrat,Li2(5),wlog(5)
      arg(1)=m1sq/s
      wlog(1)=Lnrat(-m1sq,-s)

      arg(2)=m1sq/t
      wlog(2)=Lnrat(-m1sq,-t)

      arg(3)=m3sq/s
      wlog(3)=Lnrat(-m3sq,-s)

      arg(4)=m3sq/t
      wlog(4)=Lnrat(-m3sq,-t)

      arg(5)=arg(1)*arg(4)
      wlog(5)=Lnrat(-m1sq,-s)+Lnrat(-m3sq,-t)

      do j=1,5
         omarg(j)=one-arg(j)
         if (omarg(j) .gt. one) then 
             Li2(j)=-ddilog(arg(j))+pisqo6-wlog(j)*log(omarg(j))
          else
             Li2(j)=dcmplx(ddilog(omarg(j)))
         endif
      enddo
      Lsm1_2me=Li2(5)-Li2(1)-Li2(2)-Li2(3)-Li2(4)-half*Lnrat(-s,-t)**2
      return
      end

