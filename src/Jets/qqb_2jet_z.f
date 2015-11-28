      subroutine qqb_2jet_z(p,z)
************************************************************************
*     Author: J.M. Campbell                                            *
*     November, 2002.                                                  *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'agq.f'
      include 'PR_twojet.f'
      double precision z,p(mxpart,4),dot
      double precision xl12,xl13,xl14,xl23,xl24,xl34
      double precision 
     .                 ii_qq,ii_qg,ii_gq,ii_gg,
     .                 if_qq,if_gg,
     .                 fi_qq,fi_gg,
     .                 ff_qq,ff_gg
      double precision tempgq,tempqg
      integer is,cs
      
      xl12=dlog(+two*dot(p,1,2)/musq)
      xl13=dlog(-two*dot(p,1,3)/musq)
      xl14=dlog(-two*dot(p,1,4)/musq)
      xl23=dlog(-two*dot(p,2,3)/musq)
      xl24=dlog(-two*dot(p,2,4)/musq)
      xl34=dlog(+two*dot(p,3,4)/musq)

c---- gg piece
      do is=1,3
c      S1(g,g,g,gf_gf,0,is)=ason4pi*xn*(
c     . +1.5d0*if_gg(z,xl13,is)+1.5d0*fi_gg(z,xl13,is)
c     . +1.5d0*if_gg(z,xl14,is)+1.5d0*fi_gg(z,xl14,is)
c     . +3d0*ii_gg(z,xl12,is)+1.5d0*ff_gg(z,xl34,is))
c      S1(g,g,g,gf_gf,1,is)=ason4pi*xn*(
c     . +2.5d0*if_gg(z,xl13,is)+2.5d0*fi_gg(z,xl13,is)
c     . +2d0*if_gg(z,xl14,is)+2d0*fi_gg(z,xl14,is)
c     . +1.5d0*ii_gg(z,xl12,is)+0.75d0*ff_gg(z,xl34,is))
c      S1(g,g,g,gf_gf,2,is)=ason4pi*xn*(
c     . +2.5d0*if_gg(z,xl14,is)+2.5d0*fi_gg(z,xl14,is)
c     . +2d0*if_gg(z,xl13,is)+2d0*fi_gg(z,xl13,is)
c     . +1.5d0*ii_gg(z,xl12,is)+0.75d0*ff_gg(z,xl34,is))
c      S2(g,g,g,gf_gf,0,is)=ason4pi*xn*(
c     . +1.5d0*if_gg(z,xl23,is)+1.5d0*fi_gg(z,xl23,is)
c     . +1.5d0*if_gg(z,xl24,is)+1.5d0*fi_gg(z,xl24,is)
c     . +3d0*ii_gg(z,xl12,is)+1.5d0*ff_gg(z,xl34,is))
c      S2(g,g,g,gf_gf,1,is)=ason4pi*xn*(
c     . +2.5d0*if_gg(z,xl24,is)+2.5d0*fi_gg(z,xl24,is)
c     . +2d0*if_gg(z,xl23,is)+2d0*fi_gg(z,xl23,is)
c     . +1.5d0*ii_gg(z,xl12,is)+0.75d0*ff_gg(z,xl34,is))
c      S2(g,g,g,gf_gf,2,is)=ason4pi*xn*(
c     . +2.5d0*if_gg(z,xl23,is)+2.5d0*fi_gg(z,xl23,is)
c     . +2d0*if_gg(z,xl24,is)+2d0*fi_gg(z,xl24,is)
c     . +1.5d0*ii_gg(z,xl12,is)+0.75d0*ff_gg(z,xl34,is))

      S1(g,g,g,qf_af,0,is)=ason4pi*xn*(
     . +if_gg(z,xl13,is)+fi_qq(z,xl13,is)
     . +if_gg(z,xl14,is)+fi_qq(z,xl14,is)
     . -ff_qq(z,xl34,is)*(1d0+1d0/xnsq))
      S1(g,g,g,qf_af,1,is)=ason4pi*xn*(
     . +if_gg(z,xl13,is)+fi_qq(z,xl13,is)
     . +ii_gg(z,xl12,is)
     . -ff_qq(z,xl34,is)/xnsq)
      S1(g,g,g,qf_af,2,is)=ason4pi*xn*(
     . +if_gg(z,xl14,is)+fi_qq(z,xl14,is)
     . +ii_gg(z,xl12,is)
     . -ff_qq(z,xl34,is)/xnsq)
      S2(g,g,g,qf_af,0,is)=ason4pi*xn*(
     . +if_gg(z,xl23,is)+fi_qq(z,xl23,is)
     . +if_gg(z,xl24,is)+fi_qq(z,xl24,is)
     . -ff_qq(z,xl34,is)*(1d0+1d0/xnsq))
      S2(g,g,g,qf_af,1,is)=ason4pi*xn*(
     . +if_gg(z,xl24,is)+fi_qq(z,xl24,is)
     . +ii_gg(z,xl12,is)
     . -ff_qq(z,xl34,is)/xnsq)
      S2(g,g,g,qf_af,2,is)=ason4pi*xn*(
     . +if_gg(z,xl23,is)+fi_qq(z,xl23,is)
     . +ii_gg(z,xl12,is)
     . -ff_qq(z,xl34,is)/xnsq)
     
      do cs=0,2
      S1(q,g,g,qf_gf,cs,is)=
     . ason4pi*2d0*cf*(avegg/aveqg)*ii_qg(z,xl12,is)
      S2(q,g,g,qf_gf,cs,is)=
     . ason4pi*2d0*cf*(avegg/aveqg)*ii_qg(z,xl12,is)
      S1(a,g,g,af_gf,cs,is)=
     . ason4pi*2d0*cf*(avegg/aveqg)*ii_qg(z,xl12,is)
      S2(a,g,g,af_gf,cs,is)=
     . ason4pi*2d0*cf*(avegg/aveqg)*ii_qg(z,xl12,is)
      enddo
      
      enddo

      return
      end
      
