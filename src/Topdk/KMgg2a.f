      subroutine KMgg2a(mt,zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm,
     & Ma1ba,Ma1dd,Ma2ba,Ma2dd,Ma34ba,Ma34dd)
      implicit none
      include 'constants.f'
c----- \bibitem{Korner:2002hy}
c----- J.~G.~Korner and Z.~Merebashvili,
c----- %``One-loop corrections to four-point functions with two external massive
c----- %fermions and two external massless partons,''
c----- Phys.\ Rev.\  D {\bf 66}, 054023 (2002)
c----- [arXiv:hep-ph/0207054].
c----- %%CITATION = PHRVA,D66,054023;%%
C-----This is an implementation of KM, Eq. (3.1)


      double complex zp1(4),zp2(4),zp3(4),zp4(4)
      double complex e1(4),e2(4),Ubm(4),Vm(4)
      double complex cdot,t
      double complex Bt,Ma1ba,Ma1dd,Ma2ba,Ma2dd,Ma34ba,Ma34dd
      double complex string0,string1,string2
      double complex ggQbQ(0:7),ggQbQL(1:7,1:5)
      double complex b0,bp(5),bm(2),bn(2),bc,bd(2),be(2),bg(5)
      double precision mt
      integer j

      t=2d0*cdot(zp1,zp3)

c--- Eq. (3.3) (note additional factor wrt. Eq. (2.6))
      ggQbQ(0)=im*t*Bt(zp1,zp2,zp3,zp4,e1,e2,Ubm,Vm)

      ggQbQ(1)=-string1(Ubm,zp1,Vm)
      ggQbqL(1,1)=cdot(e1,e2)
      ggQbqL(1,2)=cdot(e1,zp3)*cdot(e2,zp3)
      ggQbqL(1,3)=cdot(e1,zp3)*cdot(e2,zp4)
      ggQbqL(1,4)=cdot(e1,zp4)*cdot(e2,zp3)
      ggQbqL(1,5)=cdot(e1,zp4)*cdot(e2,zp4)

      ggQbQ(2)=string1(Ubm,e1,Vm)
      ggQbqL(2,1)=cdot(e2,zp3)
      ggQbqL(2,2)=cdot(e2,zp4)

      ggQbQ(3)=string1(Ubm,e2,Vm)
      ggQbqL(3,1)=cdot(e1,zp3)
      ggQbqL(3,2)=cdot(e1,zp4)
      
      ggQbQ(4)=mt*string2(Ubm,e1,e2,Vm)
      
      ggQbQ(5)=-mt*string2(Ubm,e1,zp1,Vm)
      ggQbqL(5,1)=cdot(e2,zp3)
      ggQbqL(5,2)=cdot(e2,zp4)

      ggQbQ(6)=-mt*string2(Ubm,e2,zp1,Vm)
      ggQbqL(6,1)=cdot(e1,zp3)
      ggQbqL(6,2)=cdot(e1,zp4)

      ggQbQ(7)=mt*string0(Ubm,Vm)
      ggQbqL(7,1)=cdot(e1,e2)
      ggQbqL(7,2)=cdot(e1,zp3)*cdot(e2,zp3)
      ggQbqL(7,3)=cdot(e1,zp3)*cdot(e2,zp4)
      ggQbqL(7,4)=cdot(e1,zp4)*cdot(e2,zp3)
      ggQbqL(7,5)=cdot(e1,zp4)*cdot(e2,zp4)

c---  Box 2a1
      call coeffboxa1(mt,zp1,zp2,zp3,zp4,b0,bp,bm,bn,bc,bd,be,bg)
      
c      write(6,*) 'b0',b0
c      write(6,*) 'bp',bp
c      write(6,*) 'bm',bm
c      write(6,*) 'bn',bn
c      write(6,*) 'bc',bc
c      write(6,*) 'bd',bd
c      write(6,*) 'be',be
c      write(6,*) 'bg',bg
            
      Ma1ba=b0*ggQbQ(0)+bc*ggQbQ(4)
      do j=1,2
      Ma1ba=Ma1ba+ggQbQ(2)*ggQbQL(2,j)*bm(j)+ggQbQ(3)*ggQbQL(3,j)*bn(j)
     &           +ggQbQ(5)*ggQbQL(5,j)*bd(j)+ggQbQ(6)*ggQbQL(6,j)*be(j)
      enddo   
      do j=1,5
      Ma1ba=Ma1ba+ggQbQ(1)*ggQbQL(1,j)*bp(j)+ggQbQ(7)*ggQbQL(7,j)*bg(j)
      enddo   
      Ma1ba=Ma1ba*im   

c--- now multiply by overall color factors
c---  Eq. (3.6) for Boxa1:
      Ma1dd=Ma1ba/4d0
      Ma1ba=Ma1ba*(Cf-xn/2d0)

c---  Box 2a2
      call coeffboxa2(mt,zp1,zp2,zp3,zp4,b0,bp,bm,bn,bc,bd,be,bg)
      
c      write(6,*) 'b0',b0
c      write(6,*) 'bp',bp
c      write(6,*) 'bm',bm
c      write(6,*) 'bn',bn
c      write(6,*) 'bc',bc
c      write(6,*) 'bd',bd
c      write(6,*) 'be',be
c      write(6,*) 'bg',bg
            
      Ma2ba=b0*ggQbQ(0)+bc*ggQbQ(4)
      do j=1,2
      Ma2ba=Ma2ba+ggQbQ(2)*ggQbQL(2,j)*bm(j)+ggQbQ(3)*ggQbQL(3,j)*bn(j)
     &           +ggQbQ(5)*ggQbQL(5,j)*bd(j)+ggQbQ(6)*ggQbQL(6,j)*be(j)
      enddo   
      do j=1,5
      Ma2ba=Ma2ba+ggQbQ(1)*ggQbQL(1,j)*bp(j)+ggQbQ(7)*ggQbQL(7,j)*bg(j)
      enddo   
      Ma2ba=Ma2ba*im   

c--- now multiply by overall color factors
c---  Eq. (3.6) for Boxa2:
      Ma2dd=Ma2ba/4d0
      Ma2ba=Ma2ba*(xn/2d0)

c---  Box 2a4
      call coeffboxa4(mt,zp1,zp2,zp3,zp4,b0,bp,bm,bn,bc,bd,be,bg)
      
c      write(6,*) 'b0',b0
c      write(6,*) 'bp',bp
c      write(6,*) 'bm',bm
c      write(6,*) 'bn',bn
c      write(6,*) 'bc',bc
c      write(6,*) 'bd',bd
c      write(6,*) 'be',be
c      write(6,*) 'bg',bg
            
      Ma34ba=b0*ggQbQ(0)+bc*ggQbQ(4)
      do j=1,2
      Ma34ba=Ma34ba
     &           +ggQbQ(2)*ggQbQL(2,j)*bm(j)+ggQbQ(3)*ggQbQL(3,j)*bn(j)
     &           +ggQbQ(5)*ggQbQL(5,j)*bd(j)+ggQbQ(6)*ggQbQL(6,j)*be(j)
      enddo   
      do j=1,5
      Ma34ba=Ma34ba
     &           +ggQbQ(1)*ggQbQL(1,j)*bp(j)+ggQbQ(7)*ggQbQL(7,j)*bg(j)
      enddo   
      Ma34ba=Ma34ba*im   

c--- now multiply by overall color factors
c---  Eq. (3.6) for Boxa34:
      Ma34dd=Ma34ba/4d0
      Ma34ba=czip
	
      return
      end
