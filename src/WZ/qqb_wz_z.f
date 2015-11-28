      subroutine qqb_wz_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR.f'
      double precision z,xl12,p(mxpart,4),dot,ii_qg,ii_gq

      xl12=log(two*dot(p,1,2)/musq)
c----contributions for one leg

      Rqq_qb=ason2pi*cf*ii_qg(z,xl12,2)
      Rq_qbqb=Rqq_qb
      Rqbqb_q=Rqq_qb
      Rqb_qq=Rqq_qb

      Pqq_qb=ason2pi*cf*ii_qg(z,xl12,3)
      Pq_qbqb=Pqq_qb
      Pqbqb_q=Pqq_qb
      Pqb_qq=Pqq_qb
      
      Rgq_q=ason2pi*tr*ii_gq(z,xl12,2)
      Rq_gq=Rgq_q

      return
      end
