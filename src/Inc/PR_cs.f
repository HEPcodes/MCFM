c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Rab_c) and leg 2 (Rc_ab)
c--- In each case the partion labelling is:
c---       emitter   a
c---       emitted   b
c---     spectator   c
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      double precision Rgg_g_cs(0:2),Rg_gg_cs(0:2),
     .                 Pgg_g_cs(0:2),Pg_gg_cs(0:2),
     .                 Rqq_qb_cs(0:2),Rq_qbqb_cs(0:2),
     .                 Pqq_qb_cs(0:2),Pq_qbqb_cs(0:2),
     .                 Rqbqb_q_cs(0:2),Rqb_qq_cs(0:2),
     .                 Pqbqb_q_cs(0:2),Pqb_qq_cs(0:2),
     .                 Rqq_g_cs(0:2),Rq_gg_cs(0:2),
     .                 Pqq_g_cs(0:2),Pq_gg_cs(0:2),
     .                 Rqbqb_g_cs(0:2),Rqb_gg_cs(0:2),
     .                 Pqbqb_g_cs(0:2),Pqb_gg_cs(0:2),
     .                 Rgg_q_cs(0:2),Rg_qq_cs(0:2),
     .                 Pgg_q_cs(0:2),Pg_qq_cs(0:2),
     .                 Rgg_qb_cs(0:2),Rg_qbqb_cs(0:2),
     .                 Pgg_qb_cs(0:2),Pg_qbqb_cs(0:2),
     .                 Rgq_g_cs(0:2),Rg_gq_cs(0:2),
     .                 Rgqb_g_cs(0:2),Rg_gqb_cs(0:2),
     .                 Pgq_g_cs(0:2),Pg_gq_cs(0:2),
     .                 Pgqb_g_cs(0:2),Pg_gqb_cs(0:2),
     .                 Rgq_qb_cs(0:2),Rq_gqb_cs(0:2),
     .                 Rgqb_q_cs(0:2),Rqb_gq_cs(0:2),
     .                 Pgq_qb_cs(0:2),Pq_gqb_cs(0:2),
     .                 Pgqb_q_cs(0:2),Pqb_gq_cs(0:2),
     .                 Pg_qg_cs(0:2),Pqg_g_cs(0:2),
     .                 Rg_qg_cs(0:2),Rqg_g_cs(0:2),
     .                 Rqq_q_cs(0:2),Rq_qq_cs(0:2),
     .                 Rqbqb_qb_cs(0:2),Rqb_qbqb_cs(0:2),
     .                 Pqq_q_cs(0:2),Pq_qq_cs(0:2),
     .                 Pqbqb_qb_cs(0:2),Pqb_qbqb_cs(0:2)

      common/RP_cols/Rgg_g_cs,Rg_gg_cs,
     .               Pgg_g_cs,Pg_gg_cs,
     .               Rqq_qb_cs,Rq_qbqb_cs,
     .               Pqq_qb_cs,Pq_qbqb_cs,
     .               Rqbqb_q_cs,Rqb_qq_cs,
     .               Pqbqb_q_cs,Pqb_qq_cs,
     .               Rqq_g_cs,Rq_gg_cs,
     .               Pqq_g_cs,Pq_gg_cs,
     .               Rqbqb_g_cs,Rqb_gg_cs,
     .               Pqbqb_g_cs,Pqb_gg_cs,
     .               Rgg_q_cs,Rg_qq_cs,
     .               Pgg_q_cs,Pg_qq_cs,
     .               Rgg_qb_cs,Rg_qbqb_cs,
     .               Pgg_qb_cs,Pg_qbqb_cs,
     .               Rgq_g_cs,Rg_gq_cs,
     .               Rgqb_g_cs,Rg_gqb_cs,
     .               Pgq_g_cs,Pg_gq_cs,
     .               Pgqb_g_cs,Pg_gqb_cs,
     .               Rgq_qb_cs,Rq_gqb_cs,
     .               Rgqb_q_cs,Rqb_gq_cs,
     .               Pgq_qb_cs,Pq_gqb_cs,
     .               Pgqb_q_cs,Pqb_gq_cs,
     .               Pg_qg_cs,Pqg_g_cs,
     .               Rg_qg_cs,Rqg_g_cs,
     .               Rqq_q_cs,Rq_qq_cs,
     .               Rqbqb_qb_cs,Rqb_qbqb_cs,
     .               Pqq_q_cs,Pq_qq_cs,
     .               Pqbqb_qb_cs,Pqb_qbqb_cs

