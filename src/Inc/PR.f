c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Rab_c) and leg 2 (Rc_ab)
c--- In each case the partion labelling is:
c---       emitter   a
c---       emitted   b
c---     spectator   c
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      double precision Rgg_q,Rgq_q,Rg_qq,Rg_qg
      double precision Rq_gg,Rq_gq,Rqq_g,Rqg_g
      double precision Rqq_qb,Rq_qbqb,Rqb_qq,Rqbqb_q
      double precision Rgq_g,Rg_gq,Rgg_g,Rg_gg
      double precision Pgg_q,Pgq_q,Pg_qq,Pg_qg
      double precision Pq_gg,Pq_gq,Pqq_g,Pqg_g
      double precision Pqq_qb,Pq_qbqb,Pqb_qq,Pqbqb_q
      double precision Pgq_g,Pg_gq,Pgg_g,Pg_gg
      common/PR/Rgg_q,Rgq_q,Rg_qq,Rg_qg,Rq_gg,Rq_gq,Rqq_g,Rqg_g,
     & Rqq_qb,Rq_qbqb,Rqb_qq,Rqbqb_q,Rgq_g,Rg_gq,Rgg_g,Rg_gg,
     & Pgg_q,Pgq_q,Pg_qq,Pg_qg,Pq_gg,Pq_gq,Pqq_g,Pqg_g,
     & Pqq_qb,Pq_qbqb,Pqb_qq,Pqbqb_q,Pgq_g,Pg_gq,Pgg_g,Pg_gg
