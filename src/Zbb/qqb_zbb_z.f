      subroutine qqb_zbb_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR.f'
      include 'PR_cs.f'
      double precision z,p(mxpart,4),dot
      double precision xl12,xl15,xl16,xl25,xl26,xl56
      double precision ii_qg,ii_gq,if_qg,fi_qg,if_gq,
     .                 ff_qg,ii_gg,if_gg,fi_gg

      xl12=log(two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(two*dot(p,5,6)/musq)

      Rqq_qb =ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
     .                        +two/xn *(if_qg(z,xl16,2)+fi_qg(z,xl16,2))
     .                       -one/xn *(ii_qg(z,xl12,2)+ff_qg(z,xl56,2)))
      Pqq_qb =ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
     .                        +two/xn *(if_qg(z,xl16,3)+fi_qg(z,xl16,3))
     .                       -one/xn *(ii_qg(z,xl12,3)+ff_qg(z,xl56,3)))
      Rqbqb_q=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
     .                        +two/xn *(if_qg(z,xl26,2)+fi_qg(z,xl26,2))
     .                       -one/xn *(ii_qg(z,xl12,2)+ff_qg(z,xl56,2)))
      Pqbqb_q=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
     .                        +two/xn *(if_qg(z,xl26,3)+fi_qg(z,xl26,3))
     .                       -one/xn *(ii_qg(z,xl12,3)+ff_qg(z,xl56,3)))
      Rq_qbqb=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl26,2)+fi_qg(z,xl26,2))
     .                        +two/xn *(if_qg(z,xl25,2)+fi_qg(z,xl25,2))
     .                       -one/xn *(ii_qg(z,xl12,2)+ff_qg(z,xl56,2)))
      Pq_qbqb=ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl26,3)+fi_qg(z,xl26,3))
     .                        +two/xn *(if_qg(z,xl25,3)+fi_qg(z,xl25,3))
     .                       -one/xn *(ii_qg(z,xl12,3)+ff_qg(z,xl56,3)))
      Rqb_qq =ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl16,2)+fi_qg(z,xl16,2))
     .                        +two/xn *(if_qg(z,xl15,2)+fi_qg(z,xl15,2))
     .                       -one/xn *(ii_qg(z,xl12,2)+ff_qg(z,xl56,2)))
      Pqb_qq =ason2pi/2d0*((xn-two/xn)*(if_qg(z,xl16,3)+fi_qg(z,xl16,3))
     .                        +two/xn *(if_qg(z,xl15,3)+fi_qg(z,xl15,3))
     .                       -one/xn *(ii_qg(z,xl12,3)+ff_qg(z,xl56,3)))

      Rgq_q=ason2pi*tr*ii_gq(z,xl12,2)
      Rq_gq=Rgq_q

c--- need to check these factors
      Rgg_g=ason2pi/2d0*(xn*ii_gq(z,xl12,2)
     .             +one/xn*(if_gq(z,xl15,2)+if_gq(z,xl16,2)))
      Rg_gg=ason2pi/2d0*(xn*ii_gq(z,xl12,2)
     .            +one/xn*(if_gq(z,xl25,2)+if_gq(z,xl26,2)))

c--- new
      Rgg_g_cs(0)=-ason2pi/2d0/xn*(ff_qg(z,xl56,2))
     .            +ason2pi/2d0*xn*(if_gg(z,xl15,2)+if_gg(z,xl16,2)
     .                            +fi_qg(z,xl15,2)+fi_qg(z,xl16,2)
     .                            -ff_qg(z,xl56,2))                 
      Rgg_g_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl15,2)+fi_qg(z,xl15,2)
     .                           +ii_gg(z,xl12,2))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,2))
      Rgg_g_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl16,2)+fi_qg(z,xl16,2)
     .                           +ii_gg(z,xl12,2))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,2))
      Rg_gg_cs(0)=-ason2pi/2d0/xn*(ff_qg(z,xl56,2))
     .            +ason2pi/2d0*xn*(if_gg(z,xl25,2)+if_gg(z,xl26,2)
     .                            +fi_qg(z,xl25,2)+fi_qg(z,xl26,2)                 
     .                            -ff_qg(z,xl56,2))                 
      Rg_gg_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl26,2)+fi_qg(z,xl26,2)
     .                           +ii_gg(z,xl12,2))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,2))
      Rg_gg_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl25,2)+fi_qg(z,xl25,2)
     .                           +ii_gg(z,xl12,2))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,2))

      Pgg_g_cs(0)=-ason2pi/2d0/xn*(ff_qg(z,xl56,3))
     .            +ason2pi/2d0*xn*(if_gg(z,xl15,3)+if_gg(z,xl16,3)
     .                            +fi_qg(z,xl15,3)+fi_qg(z,xl16,3)
     .                            -ff_qg(z,xl56,3))                 
      Pgg_g_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl15,3)+fi_qg(z,xl15,3)
     .                           +ii_gg(z,xl12,3))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,3))
      Pgg_g_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl16,3)+fi_qg(z,xl16,3)
     .                           +ii_gg(z,xl12,3))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,3))
      Pg_gg_cs(0)=-ason2pi/2d0/xn*(ff_qg(z,xl56,3))
     .            +ason2pi/2d0*xn*(if_gg(z,xl25,3)+if_gg(z,xl26,3)
     .                            +fi_qg(z,xl25,3)+fi_qg(z,xl26,3)                 
     .                            -ff_qg(z,xl56,3))                 
      Pg_gg_cs(1)=ason2pi/2d0*xn*(if_gg(z,xl26,3)+fi_qg(z,xl26,3)
     .                           +ii_gg(z,xl12,3))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,3))
      Pg_gg_cs(2)=ason2pi/2d0*xn*(if_gg(z,xl25,3)+fi_qg(z,xl25,3)
     .                           +ii_gg(z,xl12,3))
     .           -ason2pi/2d0/xn*(ff_qg(z,xl56,3))

      return
      end

