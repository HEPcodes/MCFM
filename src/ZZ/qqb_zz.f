      subroutine qqb_zz(p,msq)
      implicit none

C----Author R.K.Ellis December 1998
C----modified by JMC to include supplementary diagrams February 1999
C----modified by RKE (following suggestion of GZ) 
c----to correct supplementary diagrams April 2011
c----Matrix element for ZZ production
c----NB: we also include virtual photons
c    in the notation of DKS
C    averaged over initial colours and spins
c    u(-p1)+dbar(-p2)-->e^-(p3) + e^+(p4)   + \mu^-(p5)+ \mu^+(p6)
c    q(-p1)+qbar(-p2)-->l'(p3)  + lbar'(p4) + l(p5)    + lbar(p6)
c    with Z-leptons couplings l1 for (5,6) and l2 for(3,4)
c          and lepton charges q2 for (5,6) and q1 for (3,4)
c----No statistical factor of 1/2 included.

      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'zerowidth.f'
      include 'srdiags.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4),
     & ave,v34(2),v56(2),q34,q56     
      double complex qqb(2,2,2),qbq(2,2,2),q_qb,qb_q
      double complex aq12(2,2,2),qa12(2,2,2),
     & aq34(2,2,2),qa34(2,2,2),aq56(2,2,2),qa56(2,2,2)
      double complex propz1,propz2,props,a6trees,cprop
      double complex prop12,prop34,prop56
      double precision FAC
      integer j,k,polq,pol34,pol56
      integer h1,h2,h3
      parameter(ave=0.25d0/xn)

      fac=-4D0*esq**2

      v34(1)=l1
      v34(2)=r1
      q34=q1
      v56(1)=l2
      v56(2)=r2
      q56=q2
c--set msq=0 to initalize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--   s returned from sprod (common block) is 2*dot product
      call spinoru(6,p,za,zb)

c--   calculate propagators
      if (zerowidth) then
      prop12=s(1,2)/dcmplx(s(1,2)-zmass**2,zmass*zwidth)
      prop34=s(3,4)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      cprop=dcmplx(1d0)
      else
      prop12=dcmplx(s(1,2)/(s(1,2)-zmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-zmass**2))
      prop56=dcmplx(s(5,6)/(s(5,6)-zmass**2))
      propz1=(s(3,4)-zmass**2)/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      propz2=(s(5,6)-zmass**2)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      props=(s(1,2)-zmass**2)/dcmplx(s(1,2)-zmass**2,zmass*zwidth)
      cprop=propz1*propz2*props
      endif
            
c-- here the labels correspond to the polarizations of the
c-- quark, outgoing lepton 3 and outgoing lepton 5 respectively
      call zzamp(1,2,3,4,5,6,za,zb,aq12,aq34,aq56)
      call zzamp(2,1,3,4,5,6,za,zb,qa12,qa34,qa56)
      
c      aq12(1,1,1)=A6trees(1,2,6,5,4,3,za,zb) 
c      aq12(1,1,2)=A6trees(1,2,6,5,3,4,za,zb) 
c      aq12(1,2,1)=A6trees(1,2,5,6,4,3,za,zb) 
c      aq12(1,2,2)=A6trees(1,2,5,6,3,4,za,zb) 

c      qa12(1,1,1)=A6trees(2,1,6,5,4,3,za,zb)
c      qa12(1,1,2)=A6trees(2,1,6,5,3,4,za,zb) 
c      qa12(1,2,1)=A6trees(2,1,5,6,4,3,za,zb) 
c      qa12(1,2,2)=A6trees(2,1,5,6,3,4,za,zb)
      
c      if (srdiags) then
cc---for supplementary diagrams.
c      aq56(1,1,1)=+A6trees(3,4,1,2,5,6,za,zb)
c      aq34(1,1,1)=+A6trees(6,5,1,2,4,3,za,zb)
c      aq56(1,1,2)=-A6trees(4,3,1,2,5,6,za,zb)
c      aq34(1,1,2)=+A6trees(6,5,1,2,3,4,za,zb)      
c      aq56(1,2,1)=+A6trees(3,4,1,2,6,5,za,zb)
c      aq34(1,2,1)=-A6trees(5,6,1,2,4,3,za,zb)
c      aq56(1,2,2)=-A6trees(4,3,1,2,6,5,za,zb)
c      aq34(1,2,2)=-A6trees(5,6,1,2,3,4,za,zb)

c      qa56(1,1,1)=+A6trees(3,4,2,1,5,6,za,zb)
c      qa34(1,1,1)=+A6trees(6,5,2,1,4,3,za,zb)
c      qa56(1,1,2)=-A6trees(4,3,2,1,5,6,za,zb)
c      qa34(1,1,2)=+A6trees(6,5,2,1,3,4,za,zb)      
c      qa56(1,2,1)=+A6trees(3,4,2,1,6,5,za,zb)
c      qa34(1,2,1)=-A6trees(5,6,2,1,4,3,za,zb)
c      qa56(1,2,2)=-A6trees(4,3,2,1,6,5,za,zb)
c      qa34(1,2,2)=-A6trees(5,6,2,1,3,4,za,zb)
c      endif

c      do j=1,2
c      do k=1,2
c      aq12(2,j,k)=-qa12(1,j,k)
c      qa12(2,j,k)=-aq12(1,j,k)
c      aq56(2,j,k)=qa56(1,j,k)
c      qa56(2,j,k)=aq56(1,j,k)
c      aq34(2,j,k)=qa34(1,j,k)
c      qa34(2,j,k)=aq34(1,j,k)
c      enddo
c      enddo

      do j=-nf,nf
      k=-j
      msq(j,k)=0d0

      if (j.eq.0) go to 20
      if ((j .gt. 0).and.(k .lt. 0)) then
      do polq=1,2
      do pol56=1,2
      do pol34=1,2
      if     (polq .eq. 1) then
       q_qb=(prop56*v56(pol56)*l(j)+q56*q(j))
     &     *(prop34*v34(pol34)*l(j)+q34*q(j))*qa12(polq,pol34,pol56)
        if (srdiags) then
          q_qb=q_qb
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*l(j)+q34*q(j))*qa56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*l(j)+q56*q(j))*qa34(polq,pol34,pol56)
        endif
      elseif (polq .eq. 2) then
       q_qb=(prop56*v56(pol56)*r(j)+q56*q(j))
     &     *(prop34*v34(pol34)*r(j)+q34*q(j))*qa12(polq,pol34,pol56)
        if (srdiags) then
         q_qb=q_qb
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*r(j)+q34*q(j))*qa56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*r(j)+q56*q(j))*qa34(polq,pol34,pol56)
        endif
      endif
C-- Inclusion of width a la Baur and Zeppenfeld
      q_qb=fac*q_qb*cprop
      msq(j,k)=msq(j,k)+ave*abs(q_qb)**2
      enddo
      enddo
      enddo

      elseif ((j .lt. 0).and.(k .gt. 0)) then

      do polq=1,2
      do pol34=1,2
      do pol56=1,2
      if     (polq .eq. 1) then
       qb_q=(prop56*v56(pol56)*l(k)+q56*q(k))
     &     *(prop34*v34(pol34)*l(k)+q34*q(k))*aq12(polq,pol34,pol56)
        if (srdiags) then
         qb_q=qb_q
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*l(k)+q34*q(k))*aq56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*l(k)+q56*q(k))*aq34(polq,pol34,pol56)
        endif
      elseif (polq .eq. 2) then
       qb_q=(prop56*v56(pol56)*r(k)+q56*q(k))
     &     *(prop34*v34(pol34)*r(k)+q34*q(k))*aq12(polq,pol34,pol56)
        if (srdiags) then
         qb_q=qb_q
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*r(k)+q34*q(k))*aq56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*r(k)+q56*q(k))*aq34(polq,pol34,pol56)
        endif
      endif
C-- Inclusion of width a la Baur and Zeppenfeld
      qb_q=fac*qb_q*cprop
      msq(j,k)=msq(j,k)+ave*abs(qb_q)**2

      enddo
      enddo
      enddo

      endif


 20   continue
      enddo

      return
      end
