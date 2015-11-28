      subroutine qqb_wbb_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR.f'
      double precision z,p(mxpart,4),dot
      double precision xl12,xl15,xl16,xl25,xl26
      double precision ii_qg,ii_gq,if_qg,fi_qg

      xl12=log(two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)

      Rqq_qb=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
     .                       +two/xn *(if_qg(z,xl16,2)+fi_qg(z,xl16,2))
     .                       -one/xn * ii_qg(z,xl12,2))
      Pqq_qb=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
     .                       +two/xn *(if_qg(z,xl16,3)+fi_qg(z,xl16,3))
     .                       -one/xn * ii_qg(z,xl12,3))
      Rqbqb_q=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
     .                        +two/xn *(if_qg(z,xl26,2)+fi_qg(z,xl26,2))
     .                        -one/xn * ii_qg(z,xl12,2))
      Pqbqb_q=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
     .                        +two/xn *(if_qg(z,xl26,3)+fi_qg(z,xl26,3))
     .                        -one/xn * ii_qg(z,xl12,3))
      Rq_qbqb=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl26,2)+fi_qg(z,xl26,2))
     .                        +two/xn *(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
     .                        -one/xn * ii_qg(z,xl12,2))
      Pq_qbqb=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl26,3)+fi_qg(z,xl26,3))
     .                        +two/xn *(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
     .                        -one/xn * ii_qg(z,xl12,3))
      Rqb_qq=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl16,2)+fi_qg(z,xl16,2))
     .                       +two/xn *(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
     .                       -one/xn * ii_qg(z,xl12,2))
      Pqb_qq=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl16,3)+fi_qg(z,xl16,3))
     .                       +two/xn *(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
     .                       -one/xn * ii_qg(z,xl12,3))

      Rgq_q=ason2pi*tr*ii_gq(z,xl12,2)
      Rq_gq=Rgq_q

      return
      end

