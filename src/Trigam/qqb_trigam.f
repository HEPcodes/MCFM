

      subroutine qqb_trigam(p,msq) 
      implicit none 
!----- fill msq for pp->3gam 
      include 'constants.f' 
      include 'zprods_decl.f' 
      include 'ewcouple.f' 
      include 'ewcharge.f'

      double precision p(mxpart,4),msq(-nf:nf,-nf:nf) 
      integer j,k 
      integer h1,h2,h3,h4
      double complex qqb3gam(2,2,2,2),qbq3gam(2,2,2,2)  
      double precision qqbsum(nf),qbqsum(nf)
      double precision fac,statfac
      parameter(statfac=one/6d0)
      complex Kphase

      fac=esq**3*xn*8d0*aveqq*statfac 

      do j=-nf,nf 
         do k=-nf,nf 
            msq(j,k)=0d0 
         enddo
      enddo
      qqbsum(:)=0d0 
      qbqsum(:)=0d0 
      do h1=1,2 
         do h2=1,2 
            do h3=1,2 
               do h4=1,2
                  qqb3gam(h1,h2,h3,h4)=czip
                  qbq3gam(h1,h2,h3,h4)=czip
               enddo
            enddo
         enddo
      enddo
      
      call spinoru(5,p,za,zb)
      
!------ fill qqb and qbq => 3 gam helicity amplitudes
      call amp_lo_3gam(1,2,3,4,5,za,zb,qqb3gam) 
      call amp_lo_3gam(2,1,3,4,5,za,zb,qbq3gam) 

      do j=1,nf
      do h1=1,2 
         do h2=1,2 
            do h3=1,2 
               do h4=1,2
              qqbsum(j)=qqbsum(j)+Q(j)**6*cdabs(qqb3gam(h1,h2,h3,h4))**2
              qbqsum(j)=qbqsum(j)+Q(j)**6*cdabs(qbq3gam(h1,h2,h3,h4))**2  
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
      
      subroutine amp_lo_3gam(p1,p2,p3,p4,p5,za,zb,amp) 
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5
      double complex amp(2,2,2,2),trigam
      integer h1,h2,h3
!----- sum over q-qbar helicities 

!====== zero amplitudes
      amp(2,1,1,1)=czip 
      amp(2,2,2,2)=czip
      amp(1,2,2,2)=czip
      amp(1,1,1,1)=czip

!======= 3,4,5 symmetry 
      amp(1,2,2,1)=trigam(p1,p2,p3,p4,p5,za,zb)
      amp(1,2,1,2)=trigam(p1,p2,p3,p5,p4,za,zb)
      amp(1,1,2,2)=trigam(p1,p2,p4,p5,p3,za,zb)
       
!======= line reversal 
      amp(2,2,2,1)=trigam(p2,p1,p3,p4,p5,za,zb)
      amp(2,2,1,2)=trigam(p2,p1,p3,p5,p4,za,zb)
      amp(2,1,2,2)=trigam(p2,p1,p4,p5,p3,za,zb)

!======= conjugation 
      amp(2,1,1,2)=-trigam(p1,p2,p3,p4,p5,zb,za)
      amp(2,1,2,1)=-trigam(p1,p2,p3,p5,p4,zb,za)      
      amp(2,2,1,1)=-trigam(p1,p2,p4,p5,p3,zb,za)
!=====conjugation and line revseral 
      amp(1,1,1,2)=-trigam(p2,p1,p3,p4,p5,zb,za)
      amp(1,1,2,1)=-trigam(p2,p1,p3,p5,p4,zb,za)
      amp(1,2,1,1)=-trigam(p2,p1,p4,p5,p3,zb,za)

      return 
      end 



      double complex function trigam(p1,p2,p3,p4,p5,za,zb) 
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5 
!----- amplitude for q(-,-p1),qb(+,-p2),gam(3,+),gam(4,+),gam(5,-) 
!----- all momentum outgoing 
      
      trigam=za(p2,p1)*za(p1,p5)**2/
     &     (za(p1,p3)*za(p1,p4)*za(p2,p3)*za(p2,p4))
      
      return 
      end 
