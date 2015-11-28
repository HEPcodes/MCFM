      subroutine qq_H_mad(p1,p2,p3,p4,p5,msqm)
      implicit none
      
      include 'constants.f'
 
      double precision msqm(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      integer j,k,m,n
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                     
      real*8 sdbd_dbdh,sdbd_ddbh,sdbd_ubuh,sdbd_uubh,sdbdb_dbdbh,
     .sdbu_dbuh,sdbu_udbh,sdbub_dbubh,sdbub_ubdbh,sdd_ddh,sddb_dbdh,
     .sddb_ddbh,sddb_ubuh,sddb_uubh,sdu_duh,sdu_udh,sdub_dubh,sdub_ubdh,
     .subd_dubh,subd_ubdh,subdb_dbubh,subdb_ubdbh,subu_dbdh,subu_ddbh,
     .subu_ubuh,subu_uubh,subub_ububh,sud_duh,sud_udh,sudb_dbuh,
     .sudb_udbh,suu_uuh,suub_dbdh,suub_ddbh,suub_ubuh,suub_uubh  

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
        msqm(j,k,m,n)=0d0

      if     ((j.eq.-2).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq.-2)) then
        msqm( -2, -2, -2, -2)=subub_ububh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq.-1).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( -2, -1, -2, -1)=subdb_ubdbh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( -2, -1, -1, -2)=subdb_dbubh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq. 1).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( -2, 1, -2, 1)=subd_ubdh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( -2, 1, 1, -2)=subd_dubh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( -2, 2, -2, 2)=subu_ubuh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( -2, 2, -1, 1)=subu_dbdh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( -2, 2, 1, -1)=subu_ddbh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( -2, 2, 2, -2)=subu_uubh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( -1, -2, -2, -1)=sdbub_ubdbh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq.-2).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( -1, -2, -1, -2)=sdbub_dbubh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq.-1)) then
        msqm( -1, -1, -1, -1)=sdbdb_dbdbh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( -1, 1, -2, 2)=sdbd_ubuh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( -1, 1, -1, 1)=sdbd_dbdh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( -1, 1, 1, -1)=sdbd_ddbh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( -1, 1, 2, -2)=sdbd_uubh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq. 2).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( -1, 2, -1, 2)=sdbu_dbuh(p1,p2,p3,p4,p5)
      elseif ((j.eq.-1).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( -1, 2, 2, -1)=sdbu_udbh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( 1, -2, -2, 1)=sdub_ubdh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq.-2).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( 1, -2, 1, -2)=sdub_dubh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 1, -1, -2, 2)=sddb_ubuh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 1, -1, -1, 1)=sddb_dbdh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 1, -1, 1, -1)=sddb_ddbh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 1, -1, 2, -2)=sddb_uubh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq. 1)) then
        msqm( 1, 1, 1, 1)=sdd_ddh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq. 2).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 1, 2, 1, 2)=sdu_duh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 1).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 1, 2, 2, 1)=sdu_udh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 2, -2, -2, 2)=suub_ubuh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 2, -2, -1, 1)=suub_dbdh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 2, -2, 1, -1)=suub_ddbh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 2, -2, 2, -2)=suub_uubh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( 2, -1, -1, 2)=sudb_dbuh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq.-1).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( 2, -1, 2, -1)=sudb_udbh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 2, 1, 1, 2)=sud_duh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq. 1).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 2, 1, 2, 1)=sud_udh(p1,p2,p3,p4,p5)
      elseif ((j.eq. 2).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq. 2)) then
        msqm( 2, 2, 2, 2)=suu_uuh(p1,p2,p3,p4,p5)
      endif
     
      enddo
      enddo
      enddo
      enddo
      
      return
      end
      
