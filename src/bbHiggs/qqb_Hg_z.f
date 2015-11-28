      subroutine qqb_Hg_z(p,z)
************************************************************************
*     Author: John M. Campbell                                         *
*     February, 2002                                                   *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'flags.f'
      include 'facscale.f'
      integer is
      double precision z,xl12,xl15,xl25,p(mxpart,4),dot
      double precision ii_qq,ii_qg,ii_gq,ii_gg,
     .                 if_qq,if_gg,
     .                 fi_qq,fi_gg,
     .                 ii_gg_fac,if_gg_fac,ii_qq_fac,if_qq_fac
      double precision xl12_ren,xl15_ren,xl25_ren,
     .                 xl12_fac,xl15_fac,xl25_fac

      xl12_ren=dlog(+two*dot(p,1,2)/musq)
      xl15_ren=dlog(-two*dot(p,1,5)/musq)
      xl25_ren=dlog(-two*dot(p,2,5)/musq)
      xl12_fac=dlog(+two*dot(p,1,2)/facscale**2)
      xl15_fac=dlog(-two*dot(p,1,5)/facscale**2)
      xl25_fac=dlog(-two*dot(p,2,5)/facscale**2)

c--- 2-quark terms
c--- sum over regular and plus terms
      do is=1,3
c--- No (q,qb) terms here

c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_gg_fac(z,xl12_ren,xl12_fac,is)
     &  +if_gg_fac(z,xl25_ren,xl25_fac,is)+fi_qq(z,xl25_ren,is))
      Q1(q,q,g,is)=ason4pi*xn
     & *(ii_qq_fac(z,xl12_ren,xl12_fac,is)
     & -(if_qq_fac(z,xl15_ren,xl15_fac,is)+fi_qq(z,xl15_ren,is))/xnsq)
c--- (qb,g)
      Q2(g,g,a,is)=Q2(g,g,q,is)
      Q1(a,a,g,is)=Q1(q,q,g,is)
c      Q2(q,g,a,is)=Q2(a,g,q,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_gg_fac(z,xl12_ren,xl12_fac,is)
     &  +if_gg_fac(z,xl15_ren,xl15_fac,is)+fi_qq(z,xl15_ren,is))
      Q2(q,q,g,is)=ason4pi*xn
     & *(ii_qq_fac(z,xl12_ren,xl12_fac,is)
     & -(if_qq_fac(z,xl25_ren,xl25_fac,is)+fi_qq(z,xl25_ren,is))/xnsq)
c--- (g,qb)
      Q1(g,g,a,is)=Q1(g,g,q,is)
      Q2(a,a,g,is)=Q2(q,q,g,is)
c      Q1(q,g,a,is)=Q1(a,g,q,is)

c--- (g,g)
      Q1(q,g,g,is)=ason4pi*2d0*tr*ii_qg(z,xl12_fac,is)
      Q1(a,g,g,is)=Q1(q,g,g,is)
      Q2(q,g,g,is)=Q1(q,g,g,is)
      Q2(a,g,g,is)=Q1(q,g,g,is)
      
      enddo

c--- 4-quark terms
      do is=1,3
      Q1(g,q,q,is)=ason4pi*(xn-1d0/xn)*ii_gq(z,xl12_fac,is)
      Q2(g,q,q,is)=ason4pi*(xn-1d0/xn)*ii_gq(z,xl12_fac,is)
      Q1(g,a,a,is)=Q1(g,q,q,is)
      Q2(g,a,a,is)=Q2(g,q,q,is)
      Q1(g,a,q,is)=Q1(g,q,q,is)
      Q2(g,a,q,is)=Q2(g,q,q,is)
      Q1(g,q,a,is)=Q1(g,q,q,is)
      Q2(g,q,a,is)=Q2(g,q,q,is)

      enddo
      
      return
      end
