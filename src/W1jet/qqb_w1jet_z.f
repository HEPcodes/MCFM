      subroutine qqb_w1jet_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR.f'
      include 'flags.f'
      double precision z,xl12,xl15,xl25,p(mxpart,4),dot
      double precision ii_qg,ii_gq,ii_gg,if_qg,if_gg,fi_qg,fi_gg,fi_gq

      xl12=log( two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl25=log(-two*dot(p,2,5)/musq)

      Rqq_qb=0d0
      Rq_qbqb=0d0
      Rqbqb_q=0d0
      Rqb_qq=0d0
      Rgg_q=0d0
      Rg_qq=0d0
      Rgq_q=0d0
      Rq_gg=0d0
      Rqq_g=0d0
      Rq_gq=0d0
      Rgq_g=0d0
      Rg_gq=0d0

      Pqq_qb=0d0
      Pq_qbqb=0d0
      Pqbqb_q=0d0
      Pqb_qq=0d0
      Pgg_q=0d0
      Pg_qq=0d0
      Pgq_q=0d0
      Pq_gg=0d0
      Pqq_g=0d0
      Pq_gq=0d0
      Pgq_g=0d0
      Pg_gq=0d0

      Rgg_g=0d0
      Rg_gg=0d0
      Pgg_g=0d0
      Pg_gg=0d0

      if (Gflag) then
c--- Regular terms

c--- (q,qb) terms
      Rqq_qb =ason2pi*0.5d0*xn*(if_qg(z,xl15,2)+0.5d0*fi_gg(z,xl15,2)
     &                            -ii_qg(z,xl12,2)/xnsq)
      Rq_qbqb=ason2pi*0.5d0*xn*(if_qg(z,xl25,2)+0.5d0*fi_gg(z,xl25,2)
     &                            -ii_qg(z,xl12,2)/xnsq)
      Rqbqb_q=Rq_qbqb
      Rqb_qq =Rqq_qb

c--- (q,g) and (g,q)
      Rgg_q=0.5d0*xn*(ii_gg(z,xl12,2)+if_gg(z,xl15,2)+fi_qg(z,xl15,2))
      Rg_qq=0.5d0*xn*(ii_qg(z,xl12,2))
     &     -0.5d0/xn*(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
      Rgq_q=tr*ii_gq(z,xl12,2)

      Rgg_q=Rgg_q*ason2pi
      Rg_qq=Rg_qq*ason2pi
      Rgq_q=Rgq_q*ason2pi
      
      Rq_gg=0.5d0*xn*(ii_gg(z,xl12,2)+if_gg(z,xl25,2)+fi_qg(z,xl25,2))
      Rqq_g=0.5d0*xn*(ii_qg(z,xl12,2))
     &     -0.5d0/xn*(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
      Rq_gq=tr*ii_gq(z,xl12,2)

      Rq_gg=Rq_gg*ason2pi
      Rqq_g=Rqq_g*ason2pi
      Rq_gq=Rq_gq*ason2pi

c--- (g,g)

      Rgq_g=0.5d0*tr*ii_gq(z,xl12,2)
      Rg_gq=0.5d0*tr*ii_gq(z,xl12,2)
      
      Rgq_g=Rgq_g*ason2pi
      Rg_gq=Rg_gq*ason2pi

c--- Plus terms

c--- (q,qb) terms
      Pqq_qb =ason2pi*0.5d0*xn*(if_qg(z,xl15,3)+0.5d0*fi_gg(z,xl15,3)
     &                            -ii_qg(z,xl12,3)/xnsq)
      Pq_qbqb=ason2pi*0.5d0*xn*(if_qg(z,xl25,3)+0.5d0*fi_gg(z,xl25,3)
     &                            -ii_qg(z,xl12,3)/xnsq)
      Pqbqb_q=Pq_qbqb
      Pqb_qq =Pqq_qb

c--- (q,g) and (g,q)
      Pgg_q=0.5d0*xn*(ii_gg(z,xl12,3)+if_gg(z,xl15,3)+fi_qg(z,xl15,3))
      Pg_qq=0.5d0*xn*(ii_qg(z,xl12,3))
     &     -0.5d0/xn*(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
      Pgq_q=tr*ii_gq(z,xl12,3)

      Pgg_q=Pgg_q*ason2pi
      Pg_qq=Pg_qq*ason2pi
      Pgq_q=Pgq_q*ason2pi
      
      Pq_gg=0.5d0*xn*(ii_gg(z,xl12,3)+if_gg(z,xl25,3)+fi_qg(z,xl25,3))
      Pqq_g=0.5d0*xn*(ii_qg(z,xl12,3))
     &     -0.5d0/xn*(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
      Pq_gq=tr*ii_gq(z,xl12,3)

      Pq_gg=Pq_gg*ason2pi
      Pqq_g=Pqq_g*ason2pi
      Pq_gq=Pq_gq*ason2pi

c--- (g,g)

      Pgq_g=0.5d0*tr*ii_gq(z,xl12,3)
      Pg_gq=0.5d0*tr*ii_gq(z,xl12,3)

      Pgq_g=Pgq_g*ason2pi
      Pg_gq=Pg_gq*ason2pi
      endif
      
      if (Qflag) then
      Rq_qbqb=Rq_qbqb+ason2pi*tr*fi_gq(z,xl25,2)
      Rq_gg=Rq_gg+ason2pi*0.5d0*(xn-1d0/xn)*
     .      (ii_gq(z,xl12,2)+ii_gq(z,xl25,2))
      Rgg_q=Rgg_q+ason2pi*0.5d0*(xn-1d0/xn)*
     .      (ii_gq(z,xl12,2)+ii_gq(z,xl15,2))
      Rqb_qq=Rqb_qq+ason2pi*tr*fi_gq(z,xl25,2)

      Pq_qbqb=Pq_qbqb+ason2pi*tr*fi_gq(z,xl25,3)
      Pq_gg=Pq_gg+ason2pi*0.5d0*(xn-1d0/xn)*
     .      (ii_gq(z,xl12,3)+ii_gq(z,xl25,3))
      Pgg_q=Pgg_q+ason2pi*0.5d0*(xn-1d0/xn)*
     .      (ii_gq(z,xl12,3)+ii_gq(z,xl15,3))
      Pqb_qq=Pqb_qq+ason2pi*tr*fi_gq(z,xl25,3)
      endif
      
      return
      end
