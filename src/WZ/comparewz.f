      subroutine comparewz(p)
      implicit none
      
      include 'constants.f'
      include 'zerowidth.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
 
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      integer i,j
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                     
      real*8 aemmz,ans,Vfac
      real*8 sudb_vmmuemep,sudb_veepvmvm,sudb_veepssb,
     .       sdub_muvmemep,sdub_emvevmvm,sdub_emvessb
      integer n1,n2

c--- implement the momentum exchange      
      do i=1,4
        if (i.lt.4) then
          j=i
        else
          j=0
        endif 
        p1(j)=-p(1,i)
        p2(j)=-p(2,i)
        p3(j)=p(3,i)
        p4(j)=p(4,i)
        p5(j)=p(5,i)
        p6(j)=p(6,i)
        p7(j)=p(7,i)
      enddo

c--- make call to MCFM routines
      aemmz=1d0/128d0      
      gwsq=fourpi*aemmz/xw
      gw=sqrt(gwsq)
      gsq=1d0
      call ckmfill(nwz)
      
      write(*,*) 'zmass, wmass',zmass,wmass
      write(*,*) 'new xw,gw',xw,gw

      write(*,*)

      call qqb_wz(p,msq)
      
      write(*,*) 'MCFM (u db)',msq(2,-1)
      write(*,*) 'MCFM (db u)',msq(-1,2)
      write(*,*) 'MCFM (d ub)',msq(1,-2)
      write(*,*) 'MCFM (ub d)',msq(-2,1)
c      write(*,*) 'MCFM ( g g)',msq(0,0)
c      write(*,*) 'MCFM ( u g)',msq(2,0)
c      write(*,*) 'MCFM (db g)',msq(-1,0)
c      write(*,*) 'MCFM ( g u)',msq(0,2)
c      write(*,*) 'MCFM (g db)',msq(0,-1)
      
      call initialize

      ans=sudb_vmmuemep(p1,p2,p3,p4,p5,p6)
      write(*,*) 'Madgraph u db -> vm, mu, em, ep',ans
      ans=sudb_vmmuemep(p2,p1,p3,p4,p5,p6)
      write(*,*) 'Madgraph db u -> vm, mu, em, ep',ans

      ans=sudb_veepvmvm(p1,p2,p3,p4,p5,p6)
      write(*,*) 'Madgraph u db -> ve, ep, vm, vm',ans
      ans=sudb_veepvmvm(p2,p1,p3,p4,p5,p6)
      write(*,*) 'Madgraph db u -> ve, ep, vm, vm',ans

      ans=sudb_veepssb(p1,p2,p3,p4,p5,p6)
      write(*,*) 'Madgraph u db -> ve, ep, s , sb',ans
      ans=sudb_veepssb(p2,p1,p3,p4,p5,p6)
      write(*,*) 'Madgraph db u -> ve, ep, s , sb',ans
            
      ans=sdub_muvmemep(p1,p2,p3,p4,p5,p6)
      write(*,*) 'Madgraph d ub -> mu, vm, em, ep',ans
      ans=sdub_muvmemep(p2,p1,p3,p4,p5,p6)
      write(*,*) 'Madgraph db u -> mu, vm, em, ep',ans
        
      ans=sdub_emvevmvm(p1,p2,p3,p4,p5,p6)
      write(*,*) 'Madgraph d ub -> em, ve, vm, vm',ans
      ans=sdub_emvevmvm(p2,p1,p3,p4,p5,p6)
      write(*,*) 'Madgraph db u -> em, ve, vm, vm',ans

      ans=sdub_emvessb(p1,p2,p3,p4,p5,p6)
      write(*,*) 'Madgraph d ub -> em, ve, b , bb',ans
      ans=sdub_emvessb(p2,p1,p3,p4,p5,p6)
      write(*,*) 'Madgraph db u -> em, ve, b , bb',ans
      pause
        
      return
      end
