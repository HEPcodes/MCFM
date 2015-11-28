      subroutine preclus(p,npar)
      implicit double precision (a-h,o-z)
      include 'constants.f'
      double precision p(mxpart,4)
      integer npar
      common /jetdef/ etminj,etmaxj,delrjj,rapmaxj,rapminj 
      common /clusdef/ rsep,jalg1,jalg2
c      common /jetcom/ icol,ji,jj,jk
      common /parmom/ ppar(4,10)
      common /phypar/ w,ipp1,ipp2,rmw,rgw,rmz,rgz,sw2,qcdl
      common/energy/sqrts
      logical first
      data first/.true./
      save first
      w=sqrts

      if (first) then
      first=.false.
      
      open(unit=55,file='jetcuts.dat',status='old',err=99)
      write(6,*) 'Reading cuts from jetcuts.dat'
      read(55,*) etminj
      write(6,*) 'etminj',etminj
      read(55,*) etmaxj
      write(6,*) 'etmaxj',etmaxj
      read(55,*) delrjj
      write(6,*) 'delrjj',delrjj
      read(55,*) rapmaxj
      write(6,*) 'rapmaxj',rapmaxj
      read(55,*) rapminj
      write(6,*) 'rapminj',rapminj
      read(55,*) jalg1
      write(6,*) 'jalg1',jalg1
      read(55,*) jalg2
      write(6,*) 'jalg2',jalg2
      close(unit=55)
      endif
      rsep=1.3d0*delrjj
      npar=3

      do j=1,4
      ppar(j,1) =p(1,j)
      ppar(j,2) =p(2,j)
      ppar(j,3) =p(6,j)
      ppar(j,4) =p(7,j)
      ppar(j,5) =p(4,j)
      ppar(j,6) =p(5,j)
      ppar(j,7) =p(3,j)
      ppar(j,8) =0d0
      ppar(j,9) =0d0
      ppar(j,10)=0d0
      enddo
      
      
      return
 99   continue
      write(6,*) 'Error reading jetscuts.dat'
      stop
      end
c      etminj=10d0
c      etmaxj=500d0
c      rapmaxj=3.5d0
c      rapminj=0d0
*
* clustering criterion
*       jalg1 = 1 ; deltaR(i,j)   < delrjj
*       jalg1 = 2 ; deltaR(i,jet) < delrjj and deltaR(j,jet) < delrjj
*       jalg1 = 3 ; kt algorithm; R = delrjj
*       jalg1 = 4 ; deltaR(i,jet) < delrjj and deltaR(j,jet) < delrjj
*                      but deltaR(i,j) < Rsep
*
c      jalg1=4
* recombination scheme
*       jalg2 = 1 is D0 eta/phi
*       jalg2 = 2 is Snowmass
*       jalg2 = 3 is 4 momentum - ET = sqrt(px**2+py**2)
*       jalg2 = 4 is 4 momentum - ET = E sin(theta)
*
c      jalg2=1
* conesize
c      delrjj=0.7d0      

*
*     if jalg1 = 4, must set rsep
*
c      rsep=1.3d0*delrjj
*
* experimental lepton cuts
*
c      etminl=25d0
c      etmis =25d0
c      delrjl=0.4d0
c      rapmaxl=1.1d0
c      rlepmin =60d0
c      rlepmax =100d0
*
* hadron rapidity coverage (missing E_t reconstruction)
*
c      raphad=4d0
*
