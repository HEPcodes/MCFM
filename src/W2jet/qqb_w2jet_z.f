      subroutine qqb_w2jet_z(p,z)
************************************************************************
*     Author: J.M. Campbell                                            *
*     January, 2001.                                                   *
*     Additions Aug. 2001, for 4Q piece                                *
************************************************************************
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the terms for the QQGG piece                     *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR.f'
      include 'PR_cs.f'
      include 'lc.f'
      include 'flags.f'
      double precision z,p(mxpart,4),dot
      double precision xl12,xl15,xl16,xl25,xl26,xl56
      double precision ii_qg,ii_gq,if_qg,fi_qg,if_gq,
     .                 ff_qg,ii_gg,if_gg,fi_gg,
     .                 ff_gq,ff_gg
      double precision dip_gq
      integer cs
      
      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(+two*dot(p,5,6)/musq)

      do cs=0,2
        Rqq_qb_cs(cs)=0d0
        Rq_qbqb_cs(cs)=0d0
        Pqq_qb_cs(cs)=0d0
        Pq_qbqb_cs(cs)=0d0
        Rqbqb_q_cs(cs)=0d0
        Rqb_qq_cs(cs)=0d0
        Pqbqb_q_cs(cs)=0d0
        Pqb_qq_cs(cs)=0d0
        Rgg_g_cs(cs)=0d0
        Rg_gg_cs(cs)=0d0
        Pgg_g_cs(cs)=0d0
        Pg_gg_cs(cs)=0d0
        Rqq_g_cs(cs)=0d0
        Rq_gg_cs(cs)=0d0
        Pqq_g_cs(cs)=0d0
        Pq_gg_cs(cs)=0d0
        Rgg_q_cs(cs)=0d0
        Rg_qq_cs(cs)=0d0
        Pgg_q_cs(cs)=0d0
        Pg_qq_cs(cs)=0d0
        Rg_gq_cs(cs)=0d0
        Rgq_g_cs(cs)=0d0
        Pg_gq_cs(cs)=0d0
        Pgq_g_cs(cs)=0d0
        Rg_gqb_cs(cs)=0d0
        Rgqb_g_cs(cs)=0d0
        Pg_gqb_cs(cs)=0d0
        Pgqb_g_cs(cs)=0d0
        Rq_gqb_cs(cs)=0d0
        Rgq_qb_cs(cs)=0d0
        Pq_gqb_cs(cs)=0d0
        Pgq_qb_cs(cs)=0d0
        Rqb_gq_cs(cs)=0d0
        Rgqb_q_cs(cs)=0d0
        Pqb_gq_cs(cs)=0d0
        Pgqb_q_cs(cs)=0d0
      enddo

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************            
      if (Gflag) then      
c--- QUARK-ANTIQUARK contributions
c--- additional final-final pieces that are 0:
c--- ason2pi/2d0*xn*ff_gg(z,xl56,2)*( (1) + (2) )
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rqq_qb_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl15,2)+fi_gg(z,xl15,2)/2d0)
      Rqq_qb_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl16,2)+fi_gg(z,xl16,2)/2d0)
      Rq_qbqb_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl26,2)+fi_gg(z,xl26,2)/2d0)
      Rq_qbqb_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl25,2)+fi_gg(z,xl25,2)/2d0)

      Pqq_qb_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl15,3)+fi_gg(z,xl15,3)/2d0)
      Pqq_qb_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl16,3)+fi_gg(z,xl16,3)/2d0)
      Pq_qbqb_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl26,3)+fi_gg(z,xl26,3)/2d0)
      Pq_qbqb_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl25,3)+fi_gg(z,xl25,3)/2d0)
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      Rqq_qb_cs(0)=Rqq_qb_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl15,2)+fi_gg(z,xl15,2)/2d0
     .  +if_qg(z,xl16,2)+fi_gg(z,xl16,2)/2d0)
      Rqq_qb_cs(1)=Rqq_qb_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,2)
      Rqq_qb_cs(2)=Rqq_qb_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,2)
      Rq_qbqb_cs(0)=Rq_qbqb_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl25,2)+fi_gg(z,xl25,2)/2d0
     .  +if_qg(z,xl26,2)+fi_gg(z,xl26,2)/2d0)
      Rq_qbqb_cs(1)=Rq_qbqb_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,2)
      Rq_qbqb_cs(2)=Rq_qbqb_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,2)

      Pqq_qb_cs(0)=Pqq_qb_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl15,3)+fi_gg(z,xl15,3)/2d0
     .  +if_qg(z,xl16,3)+fi_gg(z,xl16,3)/2d0)
      Pqq_qb_cs(1)=Pqq_qb_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      Pqq_qb_cs(2)=Pqq_qb_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      Pq_qbqb_cs(0)=Pq_qbqb_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl25,3)+fi_gg(z,xl25,3)/2d0
     .  +if_qg(z,xl26,3)+fi_gg(z,xl26,3)/2d0)
      Pq_qbqb_cs(1)=Pq_qbqb_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      Pq_qbqb_cs(2)=Pq_qbqb_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      Rqq_qb_cs(0)=Rqq_qb_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,2)
      Rq_qbqb_cs(0)=Rq_qbqb_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,2)

      Pqq_qb_cs(0)=Pqq_qb_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,3)
      Pq_qbqb_cs(0)=Pq_qbqb_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,3)
      endif

c--- ANTIQUARK-QUARK contributions
c--- additional final-final pieces that are 0:
c--- ason2pi/2d0*xn*ff_gg(z,xl56,2)*( (1) + (2) )
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rqbqb_q_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl16,2)+fi_gg(z,xl16,2)/2d0)
      Rqbqb_q_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl15,2)+fi_gg(z,xl15,2)/2d0)
      Rqb_qq_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl25,2)+fi_gg(z,xl25,2)/2d0)
      Rqb_qq_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl26,2)+fi_gg(z,xl26,2)/2d0)

      Pqbqb_q_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl16,3)+fi_gg(z,xl16,3)/2d0)
      Pqbqb_q_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl15,3)+fi_gg(z,xl15,3)/2d0)
      Pqb_qq_cs(1)=ason2pi/2d0*xn*(if_qg(z,xl25,3)+fi_gg(z,xl25,3)/2d0)
      Pqb_qq_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl26,3)+fi_gg(z,xl26,3)/2d0)
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      Rqbqb_q_cs(0)=Rqbqb_q_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl15,2)+fi_gg(z,xl15,2)/2d0
     .  +if_qg(z,xl16,2)+fi_gg(z,xl16,2)/2d0)
      Rqbqb_q_cs(1)=Rqbqb_q_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,2)
      Rqbqb_q_cs(2)=Rqbqb_q_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,2)
      Rqb_qq_cs(0)=Rqb_qq_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl25,2)+fi_gg(z,xl25,2)/2d0
     .  +if_qg(z,xl26,2)+fi_gg(z,xl26,2)/2d0)
      Rqb_qq_cs(1)=Rqb_qq_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,2)
      Rqb_qq_cs(2)=Rqb_qq_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,2)

      Pqbqb_q_cs(0)=Pqbqb_q_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl15,3)+fi_gg(z,xl15,3)/2d0
     .  +if_qg(z,xl16,3)+fi_gg(z,xl16,3)/2d0)
      Pqbqb_q_cs(1)=Pqbqb_q_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      Pqbqb_q_cs(2)=Pqbqb_q_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      Pqb_qq_cs(0)=Pqb_qq_cs(0)+ason2pi/2d0*xn*
     .  (if_qg(z,xl25,3)+fi_gg(z,xl25,3)/2d0
     .  +if_qg(z,xl26,3)+fi_gg(z,xl26,3)/2d0)
      Pqb_qq_cs(1)=Pqb_qq_cs(1)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      Pqb_qq_cs(2)=Pqb_qq_cs(2)-ason2pi/2d0/xn*ii_qg(z,xl12,3)
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      Rqbqb_q_cs(0)=Rqbqb_q_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,2)
      Rqb_qq_cs(0)=Rqb_qq_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,2)

      Pqbqb_q_cs(0)=Pqbqb_q_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,3)
      Pqb_qq_cs(0)=Pqb_qq_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .  *ii_qg(z,xl12,3)
      endif

c--- GLUON-GLUON contributions
c--- no additional final-final pieces
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rgg_g_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl15,2)+fi_qg(z,xl15,2)
     .                           +ii_gg(z,xl12,2))
      Rgg_g_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl16,2)+fi_qg(z,xl16,2)
     .                           +ii_gg(z,xl12,2))
      Rg_gg_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl26,2)+fi_qg(z,xl26,2)
     .                           +ii_gg(z,xl12,2))
      Rg_gg_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl25,2)+fi_qg(z,xl25,2)
     .                           +ii_gg(z,xl12,2))

      Pgg_g_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl15,3)+fi_qg(z,xl15,3)
     .                           +ii_gg(z,xl12,3))
      Pgg_g_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl16,3)+fi_qg(z,xl16,3)
     .                           +ii_gg(z,xl12,3))
      Pg_gg_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl26,3)+fi_qg(z,xl26,3)
     .                           +ii_gg(z,xl12,3))
      Pg_gg_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl25,3)+fi_qg(z,xl25,3)
     .                           +ii_gg(z,xl12,3))
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      Rgg_g_cs(0)=Rgg_g_cs(0)+ason2pi/2d0*xn*
     .  (if_gg(z,xl15,2)+fi_qg(z,xl15,2)
     .  +if_gg(z,xl16,2)+fi_qg(z,xl16,2))
      Rg_gg_cs(0)=Rg_gg_cs(0)+ason2pi/2d0*xn*
     .  (if_gg(z,xl25,2)+fi_qg(z,xl25,2)
     .  +if_gg(z,xl26,2)+fi_qg(z,xl26,2))
      Rgg_g_cs(1)=Rgg_g_cs(1)-ason2pi/2d0/xn*ff_qg(z,xl56,2)
      Rg_gg_cs(1)=Rg_gg_cs(1)-ason2pi/2d0/xn*ff_qg(z,xl56,2)
      Rgg_g_cs(2)=Rgg_g_cs(2)-ason2pi/2d0/xn*ff_qg(z,xl56,2)
      Rg_gg_cs(2)=Rg_gg_cs(2)-ason2pi/2d0/xn*ff_qg(z,xl56,2)

      Pgg_g_cs(0)=Pgg_g_cs(0)+ason2pi/2d0*xn*
     .  (if_gg(z,xl15,3)+fi_qg(z,xl15,3)
     .  +if_gg(z,xl16,3)+fi_qg(z,xl16,3))
      Pg_gg_cs(0)=Pg_gg_cs(0)+ason2pi/2d0*xn*
     .  (if_gg(z,xl25,3)+fi_qg(z,xl25,3)
     .  +if_gg(z,xl26,3)+fi_qg(z,xl26,3))
      Pgg_g_cs(1)=Pgg_g_cs(1)-ason2pi/2d0/xn*ff_qg(z,xl56,3)
      Pg_gg_cs(1)=Pg_gg_cs(1)-ason2pi/2d0/xn*ff_qg(z,xl56,3)
      Pgg_g_cs(2)=Pgg_g_cs(2)-ason2pi/2d0/xn*ff_qg(z,xl56,3)
      Pg_gg_cs(2)=Pg_gg_cs(2)-ason2pi/2d0/xn*ff_qg(z,xl56,3)
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      Rgg_g_cs(0)=Rgg_g_cs(0)-ason2pi/2d0*(xn+1d0/xn)*ff_qg(z,xl56,2)
      Rg_gg_cs(0)=Rg_gg_cs(0)-ason2pi/2d0*(xn+1d0/xn)*ff_qg(z,xl56,2)

      Pgg_g_cs(0)=Pgg_g_cs(0)-ason2pi/2d0*(xn+1d0/xn)*ff_qg(z,xl56,3)
      Pg_gg_cs(0)=Pg_gg_cs(0)-ason2pi/2d0*(xn+1d0/xn)*ff_qg(z,xl56,3)
      endif

c--- QUARK-GLUON contributions
c--- additional final-final pieces that are 0:
c--- ason2pi/2d0*xn*(ff_qg(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rqq_g_cs(1)=ason2pi/4d0*xn*(ii_qg(z,xl12,2)+ii_gg(z,xl12,2))
      Rqq_g_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl16,2)+fi_gg(z,xl16,2)/2d0)
      Rq_gg_cs(1)=ason2pi/4d0*xn*(2d0*if_gg(z,xl26,2)+fi_gg(z,xl26,2)
     .                           +ii_qg(z,xl12,2)+ii_gg(z,xl12,2))
      Rq_gg_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl26,2)+fi_gg(z,xl26,2)/2d0
     .                           +if_gg(z,xl25,2)+fi_qg(z,xl25,2))

      Pqq_g_cs(1)=ason2pi/4d0*xn*(ii_qg(z,xl12,3)+ii_gg(z,xl12,3))
      Pqq_g_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl16,3)+fi_gg(z,xl16,3)/2d0)
      Pq_gg_cs(1)=ason2pi/4d0*xn*(2d0*if_gg(z,xl26,3)+fi_gg(z,xl26,3)
     .                           +ii_qg(z,xl12,3)+ii_gg(z,xl12,3))
      Pq_gg_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl26,3)+fi_gg(z,xl26,3)/2d0
     .                           +if_gg(z,xl25,3)+fi_qg(z,xl25,3))
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      Rqq_g_cs(1)=Rqq_g_cs(1)
     .           -ason2pi/2d0/xn*(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
      Rqq_g_cs(2)=Rqq_g_cs(2)
     .           -ason2pi/2d0/xn*(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
      Rqq_g_cs(0)=Rqq_g_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_qg(z,xl16,2)+fi_gg(z,xl16,2)
     .            +ii_qg(z,xl12,2)+ii_gg(z,xl12,2)
     .            +ff_qg(z,xl56,2)+ff_gg(z,xl56,2)/2d0)
      Rq_gg_cs(0)=Rq_gg_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_gg(z,xl25,2)+2d0*fi_qg(z,xl25,2)
     .            +ii_qg(z,xl12,2)+ii_gg(z,xl12,2)
     .            +ff_qg(z,xl56,2)+ff_gg(z,xl56,2)/2d0)

      Pqq_g_cs(1)=Pqq_g_cs(1)
     .           -ason2pi/2d0/xn*(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
      Pqq_g_cs(2)=Pqq_g_cs(2)
     .           -ason2pi/2d0/xn*(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
      Pqq_g_cs(0)=Pqq_g_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_qg(z,xl16,3)+fi_gg(z,xl16,3)
     .            +ii_qg(z,xl12,3)+ii_gg(z,xl12,3)
     .            +ff_qg(z,xl56,3)+ff_gg(z,xl56,3)/2d0)
      Pq_gg_cs(0)=Pq_gg_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_gg(z,xl25,3)+2d0*fi_qg(z,xl25,3)
     .            +ii_qg(z,xl12,3)+ii_gg(z,xl12,3)
     .            +ff_qg(z,xl56,3)+ff_gg(z,xl56,3)/2d0)
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      Rqq_g_cs(0)=Rqq_g_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .                       *(if_qg(z,xl15,2)+fi_qg(z,xl15,2))

      Pqq_g_cs(0)=Pqq_g_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .                       *(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
      endif

c--- GLUON-QUARK contributions
c--- additional final-final pieces that are 0:
c--- ason2pi/2d0*xn*(ff_qg(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rgg_q_cs(1)=ason2pi/4d0*xn*(ii_qg(z,xl12,2)+ii_gg(z,xl12,2)
     .                           +2d0*if_gg(z,xl16,2)+fi_gg(z,xl16,2))
      Rgg_q_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl16,2)+fi_gg(z,xl16,2)/2d0
     .                           +if_gg(z,xl15,2)+fi_qg(z,xl15,2))
      Rg_qq_cs(1)=ason2pi/4d0*xn*(ii_qg(z,xl12,2)+ii_gg(z,xl12,2))
      Rg_qq_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl26,2)+fi_gg(z,xl26,2)/2d0)

      Pgg_q_cs(1)=ason2pi/4d0*xn*(ii_qg(z,xl12,3)+ii_gg(z,xl12,3)
     .                           +2d0*if_gg(z,xl16,3)+fi_gg(z,xl16,3))
      Pgg_q_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl16,3)+fi_gg(z,xl16,3)/2d0
     .                           +if_gg(z,xl15,3)+fi_qg(z,xl15,3))
      Pg_qq_cs(1)=ason2pi/4d0*xn*(ii_qg(z,xl12,3)+ii_gg(z,xl12,3))
      Pg_qq_cs(2)=ason2pi/2d0*xn*(if_qg(z,xl26,3)+fi_gg(z,xl26,3)/2d0)
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      Rg_qq_cs(1)=Rg_qq_cs(1)
     .           -ason2pi/2d0/xn*(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
      Rg_qq_cs(2)=Rg_qq_cs(2)
     .           -ason2pi/2d0/xn*(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
      Rgg_q_cs(0)=Rgg_q_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_gg(z,xl15,2)+2d0*fi_qg(z,xl15,2)
     .            +ii_qg(z,xl12,2)+ii_gg(z,xl12,2)
     .            +ff_qg(z,xl56,2)+ff_gg(z,xl56,2)/2d0)
      Rg_qq_cs(0)=Rg_qq_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_qg(z,xl26,2)+fi_gg(z,xl26,2)
     .            +ii_qg(z,xl12,2)+ii_gg(z,xl12,2)
     .            +ff_qg(z,xl56,2)+ff_gg(z,xl56,2)/2d0)

      Pg_qq_cs(1)=Pg_qq_cs(1)
     .           -ason2pi/2d0/xn*(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
      Pg_qq_cs(2)=Pg_qq_cs(2)
     .           -ason2pi/2d0/xn*(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
      Pgg_q_cs(0)=Pgg_q_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_gg(z,xl15,3)+2d0*fi_qg(z,xl15,3)
     .            +ii_qg(z,xl12,3)+ii_gg(z,xl12,3)
     .            +ff_qg(z,xl56,3)+ff_gg(z,xl56,3)/2d0)
      Pg_qq_cs(0)=Pg_qq_cs(0)+ason2pi/4d0*xn*
     .            (2d0*if_qg(z,xl26,3)+fi_gg(z,xl26,3)
     .            +ii_qg(z,xl12,3)+ii_gg(z,xl12,3)
     .            +ff_qg(z,xl56,3)+ff_gg(z,xl56,3)/2d0)
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      Rg_qq_cs(0)=Rg_qq_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .                       *(if_qg(z,xl25,2)+fi_qg(z,xl25,2))

      Pg_qq_cs(0)=Pg_qq_cs(0)-ason2pi/2d0*(xn+1d0/xn)
     .                       *(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
      endif

c--- GLUON-ANTIQUARK contributions
c--- additional final-final pieces that are 0:
c--- ason2pi/2d0*xn*(ff_qg(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      Rgg_qb_cs(1) = Rgg_q_cs(2)
      Rgg_qb_cs(2) = Rgg_q_cs(1)
      Rgg_qb_cs(0) = Rgg_q_cs(0)
      Rg_qbqb_cs(1)= Rg_qq_cs(2)    
      Rg_qbqb_cs(2)= Rg_qq_cs(1)
      Rg_qbqb_cs(0)= Rg_qq_cs(0)

      Pgg_qb_cs(1) = Pgg_q_cs(2)
      Pgg_qb_cs(2) = Pgg_q_cs(1)
      Pgg_qb_cs(0) = Pgg_q_cs(0)
      Pg_qbqb_cs(1)= Pg_qq_cs(2)    
      Pg_qbqb_cs(2)= Pg_qq_cs(1)
      Pg_qbqb_cs(0)= Pg_qq_cs(0)

c--- ANTIQUARK-GLUON contributions
c--- additional final-final pieces that are 0:
c--- ason2pi/2d0*xn*(ff_qg(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      Rqbqb_g_cs(1)= Rqq_g_cs(2)
      Rqbqb_g_cs(2)= Rqq_g_cs(1)
      Rqbqb_g_cs(0)= Rqq_g_cs(0)
      Rqb_gg_cs(1) = Rq_gg_cs(2)
      Rqb_gg_cs(2) = Rq_gg_cs(1)
      Rqb_gg_cs(0) = Rq_gg_cs(0)
      
      Pqbqb_g_cs(1)= Pqq_g_cs(2)
      Pqbqb_g_cs(2)= Pqq_g_cs(1)
      Pqbqb_g_cs(0)= Pqq_g_cs(0)
      Pqb_gg_cs(1) = Pq_gg_cs(2)
      Pqb_gg_cs(2) = Pq_gg_cs(1)      
      Pqb_gg_cs(0) = Pq_gg_cs(0)      
      
c--- (gq) dipoles from (g,g) contribution      
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rgq_g_cs(1)=ason2pi/2d0*xn*(avegg/aveqg)*ii_gq(z,xl12,2)
      Rgq_g_cs(2)=Rgq_g_cs(1)
      Rgqb_g_cs(1)=Rgq_g_cs(1)
      Rgqb_g_cs(2)=Rgq_g_cs(1)
      Rg_gq_cs(1)=Rgq_g_cs(1)
      Rg_gq_cs(2)=Rgq_g_cs(1)
      Rg_gqb_cs(1)=Rgq_g_cs(1)
      Rg_gqb_cs(2)=Rgq_g_cs(1)

      Pgq_g_cs(1)=ason2pi/2d0*xn*(avegg/aveqg)*ii_gq(z,xl12,3)
      Pgq_g_cs(2)=Pgq_g_cs(1)
      Pgqb_g_cs(1)=Pgq_g_cs(1)
      Pgqb_g_cs(2)=Pgq_g_cs(1)
      Pg_gq_cs(1)=Pgq_g_cs(1)
      Pg_gq_cs(2)=Pgq_g_cs(1)
      Pg_gqb_cs(1)=Pgq_g_cs(1)
      Pg_gqb_cs(2)=Pgq_g_cs(1)
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      dip_gq=ason2pi/2d0*(avegg/aveqg)*ii_gq(z,xl12,2)
      Rgq_g_cs(1)=Rgq_g_cs(1)-dip_gq/xn
      Rgq_g_cs(2)=Rgq_g_cs(2)-dip_gq/xn
      Rgq_g_cs(0)=Rgq_g_cs(0)+dip_gq*2d0*xn
      Rgqb_g_cs(1)=Rgqb_g_cs(1)-dip_gq/xn
      Rgqb_g_cs(2)=Rgqb_g_cs(2)-dip_gq/xn
      Rgqb_g_cs(0)=Rgqb_g_cs(0)+dip_gq*2d0*xn
      Rg_gq_cs(1)=Rg_gq_cs(1)-dip_gq/xn
      Rg_gq_cs(2)=Rg_gq_cs(2)-dip_gq/xn
      Rg_gq_cs(0)=Rg_gq_cs(0)+dip_gq*2d0*xn
      Rg_gqb_cs(1)=Rg_gqb_cs(1)-dip_gq/xn
      Rg_gqb_cs(2)=Rg_gqb_cs(2)-dip_gq/xn
      Rg_gqb_cs(0)=Rg_gqb_cs(0)+dip_gq*2d0*xn
      
      dip_gq=ason2pi/2d0*(avegg/aveqg)*ii_gq(z,xl12,3)
      Pgq_g_cs(1)=Pgq_g_cs(1)-dip_gq/xn
      Pgq_g_cs(2)=Pgq_g_cs(2)-dip_gq/xn
      Pgq_g_cs(0)=Pgq_g_cs(0)+dip_gq*2d0*xn
      Pgqb_g_cs(1)=Pgqb_g_cs(1)-dip_gq/xn
      Pgqb_g_cs(2)=Pgqb_g_cs(2)-dip_gq/xn
      Pgqb_g_cs(0)=Pgqb_g_cs(0)+dip_gq*2d0*xn
      Pg_gq_cs(1)=Pg_gq_cs(1)-dip_gq/xn
      Pg_gq_cs(2)=Pg_gq_cs(2)-dip_gq/xn
      Pg_gq_cs(0)=Pg_gq_cs(0)+dip_gq*2d0*xn
      Pg_gqb_cs(1)=Pg_gqb_cs(1)-dip_gq/xn
      Pg_gqb_cs(2)=Pg_gqb_cs(2)-dip_gq/xn
      Pg_gqb_cs(0)=Pg_gqb_cs(0)+dip_gq*2d0*xn
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      dip_gq=ason2pi/2d0*(avegg/aveqg)*ii_gq(z,xl12,2)
      Rgq_g_cs(0)=Rgq_g_cs(0)-dip_gq*(xn+1d0/xn)
      Rgqb_g_cs(0)=Rgqb_g_cs(0)-dip_gq*(xn+1d0/xn)
      Rg_gq_cs(0)=Rgq_g_cs(0)-dip_gq*(xn+1d0/xn)
      Rg_gqb_cs(0)=Rgqb_g_cs(0)-dip_gq*(xn+1d0/xn)
      
      dip_gq=ason2pi/2d0*(avegg/aveqg)*ii_gq(z,xl12,3)
      Pgq_g_cs(0)=Pgq_g_cs(0)-dip_gq*(xn+1d0/xn)
      Pgqb_g_cs(0)=Pgqb_g_cs(0)-dip_gq*(xn+1d0/xn)
      Pg_gq_cs(0)=Pgq_g_cs(0)-dip_gq*(xn+1d0/xn)
      Pg_gqb_cs(0)=Pgqb_g_cs(0)-dip_gq*(xn+1d0/xn)      
      endif
      
c--- (gq) dipoles from (q,g) contribution      
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      Rq_gqb_cs(1)=ason2pi/2d0*xn*(aveqg/aveqq)*ii_gq(z,xl12,2)
      Rq_gqb_cs(2)=Rq_gqb_cs(1)
      Pq_gqb_cs(1)=ason2pi/2d0*xn*(aveqg/aveqq)*ii_gq(z,xl12,3)
      Pq_gqb_cs(2)=Pq_gqb_cs(1)
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      dip_gq=ason2pi/2d0*(aveqg/aveqq)*ii_gq(z,xl12,2)
      Rq_gqb_cs(1)=Rq_gqb_cs(1)-dip_gq/xn
      Rq_gqb_cs(2)=Rq_gqb_cs(2)-dip_gq/xn
      Rq_gqb_cs(0)=Rq_gqb_cs(0)+dip_gq*2d0*xn

      dip_gq=ason2pi/2d0*(aveqg/aveqq)*ii_gq(z,xl12,3)
      Pq_gqb_cs(1)=Pq_gqb_cs(1)-dip_gq/xn
      Pq_gqb_cs(2)=Pq_gqb_cs(2)-dip_gq/xn
      Pq_gqb_cs(0)=Pq_gqb_cs(0)+dip_gq*2d0*xn
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      dip_gq=ason2pi/2d0*(aveqg/aveqq)*ii_gq(z,xl12,2)
      Rq_gqb_cs(0)=Rq_gqb_cs(0)-dip_gq*(xn+1d0/xn)
      dip_gq=ason2pi/2d0*(aveqg/aveqq)*ii_gq(z,xl12,3)
      Pq_gqb_cs(0)=Pq_gqb_cs(0)-dip_gq*(xn+1d0/xn)
      endif      
      
c--- (gq) dipoles from (g,q) contribution      
      Rgqb_q_cs(1)=Rq_gqb_cs(1)
      Rgqb_q_cs(2)=Rq_gqb_cs(1)
      Rgqb_q_cs(0)=Rq_gqb_cs(0)
      Pgqb_q_cs(1)=Pq_gqb_cs(1)
      Pgqb_q_cs(2)=Pq_gqb_cs(1)
      Pgqb_q_cs(0)=Pq_gqb_cs(0)
c--- (gq) dipoles from (g,qb) contribution      
      Rgq_qb_cs(1)=Rq_gqb_cs(1)
      Rgq_qb_cs(2)=Rq_gqb_cs(1)
      Rgq_qb_cs(0)=Rq_gqb_cs(0)
      Pgq_qb_cs(1)=Pq_gqb_cs(1)
      Pgq_qb_cs(2)=Pq_gqb_cs(1)
      Pgq_qb_cs(0)=Pq_gqb_cs(0)
c--- (gq) dipoles from (qb,g) contribution      
      Rqb_gq_cs(1)=Rq_gqb_cs(1)
      Rqb_gq_cs(2)=Rq_gqb_cs(1)
      Rqb_gq_cs(0)=Rq_gqb_cs(0)
      Pqb_gq_cs(1)=Pq_gqb_cs(1)
      Pqb_gq_cs(2)=Pq_gqb_cs(1)
      Pqb_gq_cs(0)=Pq_gqb_cs(0)
      endif
      
************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************            
      if (Qflag) then
c--- QUARK-QUARK contributions
c--- additional final-final: ason2pi/2d0*ff_qg(z,xl56,2)*2d0*
c---  ( (0) * (xn+one/xn) +  (1) * (two/xn) + (2) * (two/xn) ) 
      Rqq_q_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*one/xn
     . -(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*one/xn
     . +ii_qg(z,xl12,2)*(xn+one/xn))
      Rqq_q_cs(1)=ason2pi/2d0*(
     . -(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*one/xn
     . +(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*(xn-two/xn)
     . +ii_qg(z,xl12,2)*two/xn)
      Rqq_q_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*(xn-two/xn)
     . -(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*one/xn
     . +ii_qg(z,xl12,2)*two/xn)
      Rq_qq_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*one/xn
     . -(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*one/xn
     . +ii_qg(z,xl12,2)*(xn+one/xn))
      Rq_qq_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*(xn-two/xn)
     . -(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*one/xn
     . +ii_qg(z,xl12,2)*two/xn)
      Rq_qq_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*one/xn
     . +(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*(xn-two/xn)
     . +ii_qg(z,xl12,2)*two/xn)

      Pqq_q_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*one/xn
     . -(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*one/xn
     . +ii_qg(z,xl12,3)*(xn+one/xn))
      Pqq_q_cs(1)=ason2pi/2d0*(
     . -(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*one/xn
     . +(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*(xn-two/xn)
     . +ii_qg(z,xl12,3)*two/xn)
      Pqq_q_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*(xn-two/xn)
     . -(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*one/xn
     . +ii_qg(z,xl12,3)*two/xn)
      Pq_qq_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*one/xn
     . -(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*one/xn
     . +ii_qg(z,xl12,3)*(xn+one/xn))
      Pq_qq_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*(xn-two/xn)
     . -(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*one/xn
     . +ii_qg(z,xl12,3)*two/xn)
      Pq_qq_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*one/xn
     . +(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*(xn-two/xn)
     . +ii_qg(z,xl12,3)*two/xn)
     
c--- ANTIQUARK-ANTIQUARK contributions
c--- additional final-final: ason2pi/2d0*ff_qg(z,xl56,2)*2d0*
c---  ( (0) * (xn+one/xn) +  (1) * (two/xn) + (2) * (two/xn) ) 
      Rqbqb_qb_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*one/xn
     . -(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*one/xn
     . +ii_qg(z,xl12,2)*(xn+one/xn))
      Rqbqb_qb_cs(1)=ason2pi/2d0*(
     . -(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*one/xn
     . +(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*(xn-two/xn)
     . +ii_qg(z,xl12,2)*two/xn)
      Rqbqb_qb_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*(xn-two/xn)
     . -(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*one/xn
     . +ii_qg(z,xl12,2)*two/xn)
      Rqb_qbqb_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*one/xn
     . -(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*one/xn
     . +ii_qg(z,xl12,2)*(xn+one/xn))
      Rqb_qbqb_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*(xn-two/xn)
     . -(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*one/xn
     . +ii_qg(z,xl12,2)*two/xn)
      Rqb_qbqb_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*one/xn
     . +(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*(xn-two/xn)
     . +ii_qg(z,xl12,2)*two/xn)

      Pqbqb_qb_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*one/xn
     . -(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*one/xn
     . +ii_qg(z,xl12,3)*(xn+one/xn))
      Pqbqb_qb_cs(1)=ason2pi/2d0*(
     . -(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*one/xn
     . +(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*(xn-two/xn)
     . +ii_qg(z,xl12,3)*two/xn)
      Pqbqb_qb_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*(xn-two/xn)
     . -(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*one/xn
     . +ii_qg(z,xl12,3)*two/xn)
      Pqb_qbqb_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*one/xn
     . -(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*one/xn
     . +ii_qg(z,xl12,3)*(xn+one/xn))
      Pqb_qbqb_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*(xn-two/xn)
     . -(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*one/xn
     . +ii_qg(z,xl12,3)*two/xn)
      Pqb_qbqb_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*one/xn
     . +(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*(xn-two/xn)
     . +ii_qg(z,xl12,3)*two/xn)

c--- QUARK-ANTIQUARK contributions
c--- additional final-final: ason2pi/2d0*ff_qg(z,xl56,2)*2d0*
c---  ( (0) * (-one/xn) +  (1) * (-one/xn) + (2) * (xn-two/xn) ) 
      Rqq_qb_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*one/xn
     . +(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*(xn+one/xn)
     . -ii_qg(z,xl12,2)*one/xn)
      Rqq_qb_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*(xn-two/xn)
     . +(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*two/xn
     . -ii_qg(z,xl12,2)*one/xn)
      Rqq_qb_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*one/xn
     . +(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*two/xn
     . +ii_qg(z,xl12,2)*(xn-two/xn))
      Rq_qbqb_cs(0)=ason2pi/2d0*(
     . +(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*(xn+one/xn)
     . -(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*one/xn
     . -ii_qg(z,xl12,2)*one/xn)
      Rq_qbqb_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*two/xn
     . +(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*(xn-two/xn)
     . -ii_qg(z,xl12,2)*one/xn)
      Rq_qbqb_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*two/xn
     . -(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*one/xn
     . +ii_qg(z,xl12,2)*(xn-two/xn))
      
      Pqq_qb_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*one/xn
     . +(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*(xn+one/xn)
     . -ii_qg(z,xl12,3)*one/xn)
      Pqq_qb_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*(xn-two/xn)
     . +(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*two/xn
     . -ii_qg(z,xl12,3)*one/xn)
      Pqq_qb_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*one/xn
     . +(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*two/xn
     . +ii_qg(z,xl12,3)*(xn-two/xn))
      Pq_qbqb_cs(0)=ason2pi/2d0*(
     . +(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*(xn+one/xn)
     . -(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*one/xn
     . -ii_qg(z,xl12,3)*one/xn)
      Pq_qbqb_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*two/xn
     . +(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*(xn-two/xn)
     . -ii_qg(z,xl12,3)*one/xn)
      Pq_qbqb_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*two/xn
     . -(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*one/xn
     . +ii_qg(z,xl12,3)*(xn-two/xn))
      
c--- ANTIQUARK-QUARK contributions
c--- additional final-final: ason2pi/2d0*ff_qg(z,xl56,2)*2d0*
c---  ( (0) * (-one/xn) +  (1) * (-one/xn) + (2) * (xn-two/xn) ) 
      Rqbqb_q_cs(0)=ason2pi/2d0*(
     . +(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*(xn+one/xn)
     . -(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*one/xn
     . -ii_qg(z,xl12,2)*one/xn)
      Rqbqb_q_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*two/xn
     . +(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*(xn-two/xn)
     . -ii_qg(z,xl12,2)*one/xn)
      Rqbqb_q_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl15,2)+fi_qg(z,xl15,2))*two/xn
     . -(if_qg(z,xl16,2)+fi_qg(z,xl16,2))*one/xn
     . +ii_qg(z,xl12,2)*(xn-two/xn))
      Rqb_qq_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*one/xn
     . +(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*(xn+one/xn)
     . -ii_qg(z,xl12,2)*one/xn)
      Rqb_qq_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*(xn-two/xn)
     . +(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*two/xn
     . -ii_qg(z,xl12,2)*one/xn)
      Rqb_qq_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl25,2)+fi_qg(z,xl25,2))*one/xn
     . +(if_qg(z,xl26,2)+fi_qg(z,xl26,2))*two/xn
     . +ii_qg(z,xl12,2)*(xn-two/xn))

      Pqbqb_q_cs(0)=ason2pi/2d0*(
     . +(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*(xn+one/xn)
     . -(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*one/xn
     . -ii_qg(z,xl12,3)*one/xn)
      Pqbqb_q_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*two/xn
     . +(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*(xn-two/xn)
     . -ii_qg(z,xl12,3)*one/xn)
      Pqbqb_q_cs(2)=ason2pi/2d0*(
     . +(if_qg(z,xl15,3)+fi_qg(z,xl15,3))*two/xn
     . -(if_qg(z,xl16,3)+fi_qg(z,xl16,3))*one/xn
     . +ii_qg(z,xl12,3)*(xn-two/xn))
      Pqb_qq_cs(0)=ason2pi/2d0*(
     . -(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*one/xn
     . +(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*(xn+one/xn)
     . -ii_qg(z,xl12,3)*one/xn)
      Pqb_qq_cs(1)=ason2pi/2d0*(
     . +(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*(xn-two/xn)
     . +(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*two/xn
     . -ii_qg(z,xl12,3)*one/xn)
      Pqb_qq_cs(2)=ason2pi/2d0*(
     . -(if_qg(z,xl25,3)+fi_qg(z,xl25,3))*one/xn
     . +(if_qg(z,xl26,3)+fi_qg(z,xl26,3))*two/xn
     . +ii_qg(z,xl12,3)*(xn-two/xn))
      endif

      return
      end

