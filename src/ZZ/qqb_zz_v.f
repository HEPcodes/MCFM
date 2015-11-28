      subroutine qqb_zz_v(p,msqv)
      implicit none

C----Author R.K.Ellis December 1998
C----modified by JMC to include supplementary diagrams February 1999
c--- calculate the virtual matrix element squared
c----and subtraction terms for ZZ production
C----modified by RKE (following suggestion of GZ) 
c----to correct supplementary diagrams April 2011
c----NB: we also include virtual photons
C    averaged over initial colours and spins
c    u(-p1)+dbar(-p2)-->\mu^-(p4)+\mu^+(p5)+e^-(p6)+e^+(p7)
c    Notation to allow room for p3 --- gluon emission.
c----No statistical factor of 1/2 included.

      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scheme.f'
      include 'ewcharge.f'
      include 'zerowidth.f'
      include 'srdiags.f'
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & p(mxpart,4),qdks(mxpart,4),v56(2),v34(2),q56,q34,virt,
     & fac,facnlo,ave     
      double complex qa12(2,2,2),aq12(2,2,2),lqa12(2,2,2),laq12(2,2,2)
      double complex qa56(2,2,2),aq56(2,2,2),qa34(2,2,2),aq34(2,2,2)
      double complex propz1,propz2,props,a6trees,a6loops,cprop
      double complex aqqb,aqbq,bqqb,bqbq,Vpole,Vpole12,suppl
      double complex prop12,prop34,prop56
      integer j,k,polq,pol34,pol56
      parameter(ave=0.25d0/xn)

      scheme='dred'
      
      fac=-4D0*esq**2
      facnlo=ason2pi*cf

      v34(1)=l1
      v34(2)=r1
      q34=q1
      v56(1)=l2
      v56(2)=r2
      q56=q2

c--set msq=0 to initalize
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

c---calculate the lowest order matrix element
      call qqb_zz(p,msq)

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
      

      laq12(1,1,1)=A6loops(1,2,5,6,4,3,za,zb) 
      laq12(1,2,1)=A6loops(1,2,5,6,3,4,za,zb) 
      laq12(1,1,2)=A6loops(1,2,6,5,4,3,za,zb) 
      laq12(1,2,2)=A6loops(1,2,6,5,3,4,za,zb) 

      lqa12(1,1,1)=A6loops(2,1,5,6,4,3,za,zb)
      lqa12(1,2,1)=A6loops(2,1,5,6,3,4,za,zb) 
      lqa12(1,1,2)=A6loops(2,1,6,5,4,3,za,zb) 
      lqa12(1,2,2)=A6loops(2,1,6,5,3,4,za,zb)

      do j=1,2
      do k=1,2
      laq12(2,j,k)=-lqa12(1,j,k)
      lqa12(2,j,k)=-laq12(1,j,k)
      enddo
      enddo

      if (srdiags) then
c---for supplementary diagrams, loops just tree*Vpole since they're all triangle-type
      Vpole12=Vpole(s(1,2))
      endif

      do j=-nf,nf
      k=-j
      virt=0d0
      if (j.eq.0) go to 20

      if ((j .gt. 0).and.(k .lt. 0)) then
      do polq=1,2
      do pol34=1,2
      do pol56=1,2
      if     (polq .eq. 1) then
       aqqb=(prop56*v56(pol56)*l(j)+q56*q(j))
     &     *(prop34*v34(pol34)*l(j)+q34*q(j))* qa12(polq,pol34,pol56)
       bqqb=(prop56*v56(pol56)*l(j)+q56*q(j))
     &     *(prop34*v34(pol34)*l(j)+q34*q(j))*lqa12(polq,pol34,pol56)
         if (srdiags) then
         suppl=
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*l(j)+q34*q(j))*qa56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*l(j)+q56*q(j))*qa34(polq,pol34,pol56)
         aqqb=aqqb+suppl
         bqqb=bqqb+suppl*Vpole12
         endif
      elseif (polq .eq. 2) then
       aqqb=(prop56*v56(pol56)*r(j)+q56*q(j))
     &     *(prop34*v34(pol34)*r(j)+q34*q(j))* qa12(polq,pol34,pol56)
       bqqb=(prop56*v56(pol56)*r(j)+q56*q(j))
     &     *(prop34*v34(pol34)*r(j)+q34*q(j))*lqa12(polq,pol34,pol56)
         if (srdiags) then
         suppl=
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*r(j)+q34*q(j))*qa56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*r(j)+q56*q(j))*qa34(polq,pol34,pol56)
         aqqb=aqqb+suppl
         bqqb=bqqb+suppl*Vpole12
         endif
      endif
C-- Inclusion of width a la Baur and Zeppenfeld
      aqqb=FAC*aqqb*cprop
      bqqb=FAC*bqqb*cprop
      virt=virt+facnlo*ave*two*dble(dconjg(aqqb)*bqqb)
      enddo
      enddo
      enddo

      elseif ((j .lt. 0).and.(k .gt. 0)) then

      do polq=1,2
      do pol56=1,2
      do pol34=1,2
      if     (polq .eq. 1) then
       aqbq=(prop56*v56(pol56)*l(k)+q56*q(k))
     &     *(prop34*v34(pol34)*l(k)+q34*q(k))* aq12(polq,pol34,pol56)
       bqbq=(prop56*v56(pol56)*l(k)+q56*q(k))
     &     *(prop34*v34(pol34)*l(k)+q34*q(k))*laq12(polq,pol34,pol56)
         if (srdiags) then
         suppl=
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*l(k)+q34*q(k))*aq56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*l(k)+q56*q(k))*aq34(polq,pol34,pol56)
         aqbq=aqbq+suppl
         bqbq=bqbq+suppl*Vpole12
         endif
      elseif (polq .eq. 2) then
       aqbq=(prop56*v56(pol56)*r(k)+q56*q(k))
     &     *(prop34*v34(pol34)*r(k)+q34*q(k))* aq12(polq,pol34,pol56)
       bqbq=(prop56*v56(pol56)*r(k)+q56*q(k))
     &     *(prop34*v34(pol34)*r(k)+q34*q(k))*laq12(polq,pol34,pol56)
         if (srdiags) then
         suppl=
     &       +(prop56*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v34(pol34)*r(k)+q34*q(k))*aq56(polq,pol34,pol56)
     &       +(prop34*v34(pol34)*v56(pol56)+q34*q56)
     &       *(prop12*v56(pol56)*r(k)+q56*q(k))*aq34(polq,pol34,pol56)
         aqbq=aqbq+suppl
         bqbq=bqbq+suppl*Vpole12
         endif
      endif
C-- Inclusion of width a la Baur and Zeppenfeld
      aqbq=FAC*aqbq*cprop
      bqbq=FAC*bqbq*cprop
      virt=virt+facnlo*ave*two*dble(dconjg(aqbq)*bqbq)
      enddo
      enddo
      enddo

      endif

      msqv(j,k)=virt

 20   continue
      enddo

      return
      end
