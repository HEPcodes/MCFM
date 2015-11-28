      subroutine qqb_trigam_v(p,msq)
c--- Virtual matrix element for the process
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + gam(p5
c---
c--- J. M. Campbell and C. Williams, March 2013
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'scheme.f'
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      double complex virt_trigam_MHV,test 
      double complex qqb3gam_lo(2,2,2,2),qbq3gam_lo(2,2,2,2)
      double complex qqb3gam_v(2,2,2,2),qbq3gam_v(2,2,2,2)
      double complex ALO1,ALO2,ANLO1,ANLO2
      double precision qqbsum(nf),qbqsum(nf)
      integer h1,h2,h3,h4
      integer j,k
      double precision fac,statfac
      parameter(statfac=one/6d0)
      

      scheme='tH-V'
      msq(:,:)=zip

      fac=esq**3*xn*8d0*aveqq*statfac

      qqbsum(:)=0d0 
      qbqsum(:)=0d0 
      do h1=1,2 
         do h2=1,2 
            do h3=1,2 
               do h4=1,2
                  qqb3gam_lo(h1,h2,h3,h4)=czip
                  qbq3gam_lo(h1,h2,h3,h4)=czip
                  qqb3gam_v(h1,h2,h3,h4)=czip
                  qbq3gam_v(h1,h2,h3,h4)=czip
               enddo
            enddo
         enddo
      enddo
      
      call spinoru(5,p,za,zb)
     

!------ fill qqb and qbq => 3 gam lo amplitudes
      call amp_lo_3gam(1,2,3,4,5,za,zb,qqb3gam_lo) 
      call amp_lo_3gam(2,1,3,4,5,za,zb,qbq3gam_lo)
!------ fill qqb and qbq => 3 gam vitual  amplitudes
      call amp_virt_3gam(1,2,3,4,5,za,zb,qqb3gam_v) 
      call amp_virt_3gam(2,1,3,4,5,za,zb,qbq3gam_v) 


      do j=1,nf
         do h1=1,2 
            do h2=1,2 
               do h3=1,2 
                  do h4=1,2
                     ALO1=qqb3gam_lo(h1,h2,h3,h4)
                     ANLO1=qqb3gam_v(h1,h2,h3,h4)
                     qqbsum(j)=qqbsum(j)+Q(j)**6*
     &    cf*ason2pi*(Dble(Dconjg(ALO1)*ANLO1)+Dble(Dconjg(ANLO1)*ALO1))
                     ALO2=qbq3gam_lo(h1,h2,h3,h4)
                     ANLO2=qbq3gam_v(h1,h2,h3,h4)
                     qbqsum(j)=qbqsum(j)+Q(j)**6*
     &    cf*ason2pi*(Dble(Dconjg(ALO2)*ANLO2)+Dble(Dconjg(ANLO2)*ALO2))
                  enddo
            enddo
         enddo     
      enddo
      enddo

      
      do j=-nf,nf 
         k=-j 
         if(j.lt.0) then 
            msq(j,k)=fac*qbqsum(-j)
         elseif(j.eq.0) then 
            msq(j,k)=0d0 
         elseif(j.gt.0) then 
            msq(j,k)=fac*qqbsum(j)
         endif
      enddo
      
      return 
      end              
      
      

      subroutine amp_virt_3gam(p1,p2,p3,p4,p5,za,zb,amp) 
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5
      double complex amp(2,2,2,2),virt_trigam_MHV
      integer h1,h2,h3
!----- sum over q-qbar helicities 

!====== zero amplitudes
      amp(2,1,1,1)=czip 
      amp(2,2,2,2)=czip
      amp(1,2,2,2)=czip
      amp(1,1,1,1)=czip

!======= 3,4,5 symmetry 
      amp(1,2,2,1)=virt_trigam_MHV(p1,p2,p3,p4,p5,za,zb)
      amp(1,2,1,2)=virt_trigam_MHV(p1,p2,p3,p5,p4,za,zb)
      amp(1,1,2,2)=virt_trigam_MHV(p1,p2,p4,p5,p3,za,zb)
      
!======= line reversal 
      amp(2,2,2,1)=virt_trigam_MHV(p2,p1,p3,p4,p5,za,zb)
      amp(2,2,1,2)=virt_trigam_MHV(p2,p1,p3,p5,p4,za,zb)
      amp(2,1,2,2)=virt_trigam_MHV(p2,p1,p4,p5,p3,za,zb)

!======= conjugation 
      amp(2,1,1,2)=-virt_trigam_MHV(p1,p2,p3,p4,p5,zb,za)
      amp(2,1,2,1)=-virt_trigam_MHV(p1,p2,p3,p5,p4,zb,za)      
      amp(2,2,1,1)=-virt_trigam_MHV(p1,p2,p4,p5,p3,zb,za)
!=====conjugation and line revseral 
      amp(1,1,1,2)=-virt_trigam_MHV(p2,p1,p3,p4,p5,zb,za)
      amp(1,1,2,1)=-virt_trigam_MHV(p2,p1,p3,p5,p4,zb,za)
      amp(1,2,1,1)=-virt_trigam_MHV(p2,p1,p4,p5,p3,zb,za)

      return 
      end
