      subroutine qq_H_mad_g(p1,p2,p3,p4,p5,p6,msqm)
      implicit none
      
      include 'constants.f'
 
      double precision msqm(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      integer j,k,m,n
      REAL*8 P1(0:3),P2(0:3),P3(0:3),P4(0:3),P5(0:3),P6(0:3),P7(0:3)                     
      real*8 
     . sdbd_dbdhg,sdbd_ddbhg,sdbd_ubuhg,sdbd_uubhg,sdbdb_dbdbhg,
     . sdbg_dbdbhd,sdbg_dbdhdb,sdbg_dbubhu,sdbg_dbuhub,sdbg_ddbhdb,
     . sdbg_ubdbhu,sdbg_ubuhdb,sdbg_udbhub,sdbg_uubhdb,sdbu_dbuhg,
     . sdbu_udbhg,sdbub_dbubhg,sdbub_ubdbhg,sdd_ddhg,sddb_dbdhg,
     . sddb_ddbhg,sddb_ubuhg,sddb_uubhg,sdg_dbdhd,sdg_ddbhd,sdg_ddhdb,
     . sdg_dubhu,sdg_duhub,sdg_ubdhu,sdg_ubuhd,sdg_udhub,sdg_uubhd,
     . sdu_duhg,sdu_udhg,sdub_dubhg,sdub_ubdhg,sgd_dbdhd,sgd_ddbhd,
     . sgd_ddhdb,sgd_dubhu,sgd_duhub,sgd_ubdhu,sgd_ubuhd,sgd_udhub,
     . sgd_uubhd,sgdb_dbdbhd,sgdb_dbdhdb,sgdb_dbubhu,sgdb_dbuhub,
     . sgdb_ddbhdb,sgdb_ubdbhu,sgdb_ubuhdb,sgdb_udbhub,sgdb_uubhdb,
     . sgu_dbdhu,sgu_dbuhd,sgu_ddbhu,sgu_duhdb,sgu_ubuhu,sgu_udbhd,
     . sgu_udhdb,sgu_uubhu,sgu_uuhub,sgub_dbdhub,sgub_dbubhd,
     . sgub_ddbhub,sgub_dubhdb,sgub_ubdbhd,sgub_ubdhdb,sgub_ububhu,
     . sgub_ubuhub,sgub_uubhub,subd_dubhg,subd_ubdhg,subdb_dbubhg,
     . subdb_ubdbhg,subg_dbdhub,subg_dbubhd,subg_ddbhub,subg_dubhdb,
     . subg_ubdbhd,subg_ubdhdb,subg_ububhu,subg_ubuhub,subg_uubhub,
     . subu_dbdhg,subu_ddbhg,subu_ubuhg,subu_uubhg,subub_ububhg,
     . sud_duhg,sud_udhg,sudb_dbuhg,sudb_udbhg,sug_dbdhu,sug_dbuhd,
     . sug_ddbhu,sug_duhdb,sug_ubuhu,sug_udbhd,sug_udhdb,sug_uubhu,
     . sug_uuhub,suu_uuhg,suub_dbdhg,suub_ddbhg,suub_ubuhg,suub_uubhg
                                                                          
      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
        msqm(j,k,m,n)=0d0

      if     ((j.eq.-2).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq.-2)) then
        msqm( -2, -2, -2, -2)=subub_ububhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq.-1).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( -2, -1, -2, -1)=subdb_ubdbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( -2, -1, -1, -2)=subdb_dbubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq.-2)) then
        msqm( -2, 0, -2, -2)=subg_ububhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( -2, 0, -2, -1)=subg_ubdbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( -2, 0, -2, 1)=subg_ubdhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( -2, 0, -2, 2)=subg_ubuhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( -2, 0, -1, -2)=subg_dbubhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( -2, 0, -1, 1)=subg_dbdhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( -2, 0, 1, -2)=subg_dubhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( -2, 0, 1, -1)=subg_ddbhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( -2, 0, 2, -2)=subg_uubhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 1).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( -2, 1, -2, 1)=subd_ubdhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( -2, 1, 1, -2)=subd_dubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( -2, 2, -2, 2)=subu_ubuhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( -2, 2, -1, 1)=subu_dbdhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( -2, 2, 1, -1)=subu_ddbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-2).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( -2, 2, 2, -2)=subu_uubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( -1, -2, -2, -1)=sdbub_ubdbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq.-2).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( -1, -2, -1, -2)=sdbub_dbubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq.-1)) then
        msqm( -1, -1, -1, -1)=sdbdb_dbdbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( -1, 0, -2, -1)=sdbg_ubdbhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( -1, 0, -2, 2)=sdbg_ubuhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( -1, 0, -1, -2)=sdbg_dbubhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq.-1)) then
        msqm( -1, 0, -1, -1)=sdbg_dbdbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( -1, 0, -1, 1)=sdbg_dbdhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( -1, 0, -1, 2)=sdbg_dbuhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( -1, 0, 1, -1)=sdbg_ddbhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( -1, 0, 2, -2)=sdbg_uubhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( -1, 0, 2, -1)=sdbg_udbhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( -1, 1, -2, 2)=sdbd_ubuhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( -1, 1, -1, 1)=sdbd_dbdhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( -1, 1, 1, -1)=sdbd_ddbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 1).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( -1, 1, 2, -2)=sdbd_uubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 2).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( -1, 2, -1, 2)=sdbu_dbuhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq.-1).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( -1, 2, 2, -1)=sdbu_udbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq.-2)) then
        msqm( 0, -2, -2, -2)=sgub_ububhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( 0, -2, -2, -1)=sgub_ubdbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( 0, -2, -2, 1)=sgub_ubdhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 0, -2, -2, 2)=sgub_ubuhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( 0, -2, -1, -2)=sgub_dbubhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 0, -2, -1, 1)=sgub_dbdhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( 0, -2, 1, -2)=sgub_dubhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 0, -2, 1, -1)=sgub_ddbhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-2).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 0, -2, 2, -2)=sgub_uubhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq.-2).and.(n.eq.-1)) then
        msqm( 0, -1, -2, -1)=sgdb_ubdbhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 0, -1, -2, 2)=sgdb_ubuhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq.-2)) then
        msqm( 0, -1, -1, -2)=sgdb_dbubhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq.-1)) then
        msqm( 0, -1, -1, -1)=sgdb_dbdbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 0, -1, -1, 1)=sgdb_dbdhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( 0, -1, -1, 2)=sgdb_dbuhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 0, -1, 1, -1)=sgdb_ddbhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 0, -1, 2, -2)=sgdb_uubhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq.-1).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( 0, -1, 2, -1)=sgdb_udbhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( 0, 1, -2, 1)=sgd_ubdhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 0, 1, -2, 2)=sgd_ubuhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 0, 1, -1, 1)=sgd_dbdhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( 0, 1, 1, -2)=sgd_dubhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 0, 1, 1, -1)=sgd_ddbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq. 1)) then
        msqm( 0, 1, 1, 1)=sgd_ddhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 0, 1, 1, 2)=sgd_duhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 0, 1, 2, -2)=sgd_uubhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 1).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 0, 1, 2, 1)=sgd_udhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 0, 2, -2, 2)=sgu_ubuhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 0, 2, -1, 1)=sgu_dbdhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( 0, 2, -1, 2)=sgu_dbuhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 0, 2, 1, -1)=sgu_ddbhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 0, 2, 1, 2)=sgu_duhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 0, 2, 2, -2)=sgu_uubhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( 0, 2, 2, -1)=sgu_udbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 0, 2, 2, 1)=sgu_udhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 0).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq. 2)) then
        msqm( 0, 2, 2, 2)=sgu_uuhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( 1, -2, -2, 1)=sdub_ubdhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq.-2).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( 1, -2, 1, -2)=sdub_dubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 1, -1, -2, 2)=sddb_ubuhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 1, -1, -1, 1)=sddb_dbdhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 1, -1, 1, -1)=sddb_ddbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq.-1).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 1, -1, 2, -2)=sddb_uubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq. 1)) then
        msqm( 1, 0, -2, 1)=sdg_ubdhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 1, 0, -2, 2)=sdg_ubuhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 1, 0, -1, 1)=sdg_dbdhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq.-2)) then
        msqm( 1, 0, 1, -2)=sdg_dubhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 1, 0, 1, -1)=sdg_ddbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq. 1)) then
        msqm( 1, 0, 1, 1)=sdg_ddhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 1, 0, 1, 2)=sdg_duhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 1, 0, 2, -2)=sdg_uubhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 1, 0, 2, 1)=sdg_udhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq. 1)) then
        msqm( 1, 1, 1, 1)=sdd_ddhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 2).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 1, 2, 1, 2)=sdu_duhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 1).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 1, 2, 2, 1)=sdu_udhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 2, -2, -2, 2)=suub_ubuhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 2, -2, -1, 1)=suub_dbdhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 2, -2, 1, -1)=suub_ddbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq.-2).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 2, -2, 2, -2)=suub_uubhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq.-1).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( 2, -1, -1, 2)=sudb_dbuhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq.-1).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( 2, -1, 2, -1)=sudb_udbhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq.-2).and.(n.eq. 2)) then
        msqm( 2, 0, -2, 2)=sug_ubuhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq. 1)) then
        msqm( 2, 0, -1, 1)=sug_dbdhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq.-1).and.(n.eq. 2)) then
        msqm( 2, 0, -1, 2)=sug_dbuhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq.-1)) then
        msqm( 2, 0, 1, -1)=sug_ddbhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 2, 0, 1, 2)=sug_duhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq.-2)) then
        msqm( 2, 0, 2, -2)=sug_uubhu(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq.-1)) then
        msqm( 2, 0, 2, -1)=sug_udbhd(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 2, 0, 2, 1)=sug_udhdb(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 0).and.(m.eq. 2).and.(n.eq. 2)) then
        msqm( 2, 0, 2, 2)=sug_uuhub(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 1).and.(m.eq. 1).and.(n.eq. 2)) then
        msqm( 2, 1, 1, 2)=sud_duhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 1).and.(m.eq. 2).and.(n.eq. 1)) then
        msqm( 2, 1, 2, 1)=sud_udhg(p1,p2,p3,p4,p5,p6)
      elseif ((j.eq. 2).and.(k.eq. 2).and.(m.eq. 2).and.(n.eq. 2)) then
        msqm( 2, 2, 2, 2)=suu_uuhg(p1,p2,p3,p4,p5,p6)
      endif
     
      enddo
      enddo
      enddo
      enddo
      
      return
      end
      
