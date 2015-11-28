# Makefile routine.

MCFMHOME        = /scratch1/johnmc/MCFM_download
MYHOME          = /scratch1/johnmc/MCFM_download
SOURCEDIR	= $(MCFMHOME)/src
VPATH		= $(DIRS)
BIN		= $(MYHOME)/Bin
INCPATH  	= $(SOURCEDIR)/Inc

DIRS	=	$(SOURCEDIR)/User:\
		$(SOURCEDIR)/Spinor:$(SOURCEDIR)/Vol:\
		$(SOURCEDIR)/Need:$(SOURCEDIR)/Lib:$(SOURCEDIR)/Phase:\
		$(SOURCEDIR)/Parton:$(SOURCEDIR)/Integrate:\
		$(SOURCEDIR)/Wbb:$(SOURCEDIR)/Zbb:\
		$(SOURCEDIR)/WHbbar:$(SOURCEDIR)/ZHbbar:\
		$(SOURCEDIR)/WW:$(SOURCEDIR)/WZ:$(SOURCEDIR)/ZZ:\
		$(SOURCEDIR)/Top:$(SOURCEDIR)/Singletop:\
		$(SOURCEDIR)/TopH:$(SOURCEDIR)/TopZ:$(SOURCEDIR)/Topg:\
		$(SOURCEDIR)/HWW:$(SOURCEDIR)/HZZ:$(SOURCEDIR)/Tau:\
		$(SOURCEDIR)/Hbbbar:$(SOURCEDIR)/Httbar:\
		$(SOURCEDIR)/W:$(SOURCEDIR)/Z:\
		$(SOURCEDIR)/W1jet:$(SOURCEDIR)/Z1jet:\
		$(SOURCEDIR)/W2jet:$(SOURCEDIR)/W2jetvirt:\
		$(SOURCEDIR)/Z2jet:


FC = g77
FFLAGS 	=	-c -I$(INCPATH) 
#CRNLIB = /home/johnmc/Fortran/Cernlib/

# -----------------------------------------------------------------------------
# Specify the object files. 

HBBBARFILES = \
qqb_hbbbar.o \
qqb_hbbbar_g.o \
qqb_hbbbar_gs.o \
qqb_hbbbar_gvec.o \
qqb_hbbbar_v.o \
qqb_hbbbar_z.o

HWWFILES = \
qqb_hww.o \
qqb_hww_g.o \
qqb_hww_gvec.o \
qqb_hww_gs.o \
qqb_hww_z.o \
qqb_hww_v.o

HTAUTAUFILES = \
qqb_higgs.o \
ehsv.o \
li2.o \
Htautau.o

HZZFILES = \
qqb_hzz.o \
qqb_hzz_g.o \
qqb_hzz_gvec.o \
qqb_hzz_gs.o \
qqb_hzz_z.o \
qqb_hzz_v.o

INTEGRATEFILES = \
dvegas.o \
lenocc.o \
mbook.o \
ran0.o \
ran1.o \
rn.o

LIBFILES = \
dclaus.o \
ddilog.o

NEEDFILES = \
angle.o \
anglephi.o \
banner.o \
boost.o \
boosta.o \
branch.o \
breitw.o \
breitw1.o \
chooser.o \
ckmfill.o \
clust.o \
conserve.o \
coupling.o \
couplz.o \
cuts.o \
dittdrein.o \
dipoles.o \
dipolesub.o \
dot.o \
dotem.o \
dotpr.o \
etmiss.o \
gasdev.o \
storeptilde.o \
gtperp.o \
higgsp.o \
higgsw.o \
histofin.o \
itransform.o \
lowint.o \
masscuts.o \
mcfm.o \
preclus.o \
r.o \
reader.o \
realint.o \
smalls.o \
totint.o \
transform.o \
virtint.o \
wconstruct.o \
writeout.o \
dsigdy.o \
ddvdif.o \
gen2a.o \
gen3a.o \
gen3b.o \
phase2.o \
gen3from2.o \
wtgen.o \
wt2gen.o 

PARTONFILES = \
Ctq4Fn.o \
Ctq5Pdf.o \
alfamz.o \
fdist.o \
mrs96.o \
mrs98.o \
mrs99.o \
mrsebh.o \
mrsg.o \
mt.o \
newton1.o \
pdfwrap_linux.o

PHASEFILES = \
gen2.o \
gen3.o \
gen4.o \
gen4a.o \
gen3m.o \
gen4from3.o \
gen5.o \
gen5a.o \
gen5from4.o \
gen6.o \
gen7.o \
gencol.o \
genff.o \
genif.o \
genii.o \
genrad.o \
genrff.o \
genrif.o \
genrii.o \
kingen.o \
phase3.o \
phase3m.o \
phase4.o \
phase41.o \
phase4m.o \
phase5.o \
phase51.o \
phase5_bgd.o \
phase6.o \
phase7.o \
phase7m.o \
phi1_2.o \
phi1_2m.o \
phi1_2m_nobw.o \
phi3.o \
phi3m.o \
phi3m0.o \
rangen.o \
wt4gen.o

SINGLETOPFILES = \
qg_tbb.o \
qqb_tbb.o

TAUTAUFILES = \
qqb_tautau.o \
tautauww.o

TOPFILES = \
dotks.o \
ggttww.o \
qqb_ttb.o \
std.o \
ttbbww.o

TOPHFILES = \
ggtth.o \
qqbtth.o \
qqb_tth.o \
gen8.o \
qqb_tottth.o \
phase8.o 

TOPZFILES = \
ttbbd1.o \
qqb_ttz.o

TOPGFILES = \
qqb_ttb_g.o \
xzqqgg.o \
ttgdiags.o

USERFILES = \
genclust.o \
genclust2.o \
mdata.o \
nplotter.o \
jetcuts.o \
speccuts.o \
rkecuts.o \
wbbcuts.o \
zbbcuts.o

VOLFILES = \
qqb_vol.o \
vol.o

WFILES = \
qqb_w.o \
qqb_w_g.o \
donothing_gvec.o \
qqb_w_gs.o \
qqb_w_v.o \
qqb_w_z.o

W2JETFILES_jmc = \
aqqb_wgg.o \
qqb_w2jet.o\
qqb_w2jet_g.o\
qqb_w2jet_gvec.o\
itwojet.o \
twojet.o \
initqqqq.o \
initqqgg.o \
xmqqqq.o \
xmqqgg.o \
xwqqqq.o \
xwqqbqqb.o \
xzqqqq.o \
couplant.o \
subqcd.o \
jetcuts.o \
subqcd_new.o \
subqed.o

W2JETFILES = \
qqb_w2jet.o \
qqb_w2jet_g.o \
qqb_w2jet_gvec.o \
qqb_w2jet_gs.o \
xwqqqq.o \
xwqqbqqb.o \
twojet.o \
itwojet.o \
initqqqq.o \
initqqgg.o \
xmqqqq.o \
xzqqqq.o \
xmqqgg.o \
subqcd.o \
subqcdn.o \
prod4.o \
couplant.o \
aqqb_wgg.o

W1JETFILES = \
qqb_w1jet_gs.o \
qqb_w1jet_soft.o \
qqb_w1jet_z.o \
a5nlo.o \
virt5.o \
qqb_w1jet_v.o \
qqb_w_gvec.o 

Z1JETFILES = \
qqb_z1jet_v.o \
qqb_z1jet_z.o \
qqb_z2jet_gvec.o \
qqb_z_gvec.o \
qqb_z1jet_gs.o \
qqb_z1jet_soft.o

Z2JETFILES = \
qqb_z2jet.o \
z2jetsq.o

WHBBARFILES = \
qqb_wh.o \
qqb_wh_g.o \
qqb_wh_gs.o \
qqb_wh_v.o \
qqb_wh_z.o 

SPINORFILES = \
spinork.o \
spinoru.o

WWFILES = \
BigT.o \
L34_12.o \
a6loop.o \
a7tree.o \
b7tree.o \
fa.o \
qqb_ww.o \
qqb_ww_g.o \
qqb_ww_gs.o \
qqb_ww_v.o \
qqb_ww_z.o \
wwamps.o

WZFILES = \
qqb_wz.o \
qqb_wz_g.o \
qqb_wz_gs.o \
qqb_wz_v.o \
qqb_wz_z.o \
wzamps.o

WBBFILES = \
a6.o \
a61.o \
a61LLL.o \
a61LRL.o \
a62.o \
a6f.o \
a6routine.o \
a6tree.o \
aqqb_wbb.o \
atrLLL.o \
atrLRL.o \
atree.o \
atreepm.o \
atreepp.o \
atreesl.o \
d1.o \
d10.o \
d11.o \
d12.o \
d2.o \
d3.o \
d4.o \
d5.o \
d6.o \
d78.o \
d9.o \
fpm.o \
fpp.o \
fsl.o \
i3m.o \
lfunctions.o \
lnrat.o \
msqwbb.o \
qqb_wbb.o \
qqb_wbb_g.o \
qqb_wbb_gs.o \
qqb_wbb_soft.o \
qqb_wbb_v.o \
qqb_wbb_z.o \
sumd.o \
sumqcda.o \
sumqcdb.o \
t.o \
t1234.o \
vv.o \
wbb.o

ZFILES = \
qqb_z.o \
qqb_z_g.o \
qqb_z_gs.o \
qqb_z_v.o \
qqb_z_z.o

ZHBBARFILES = \
qqb_zh.o \
qqb_zh_g.o \
qqb_zh_gs.o \
qqb_zh_v.o \
qqb_zh_z.o

ZZFILES = \
qqb_zz.o \
qqb_zz_g.o \
qqb_zz_gs.o \
qqb_zz_v.o \
qqb_zz_z.o \
zzamps.o

ZBBFILES = \
aqqb_zbb.o \
atreez.o \
qqb_zbb.o \
qqb_zbb_g.o \
qqb_zbb_gs.o \
qqb_zbb_soft.o \
qqb_zbb_v.o \
qqb_zbb_z.o \
a6treeg.o \
xzqqgg_v.o \
xzqqggg.o \
amp_qqgg.o \
amp_qqggg.o \
a61g.o \
a6g.o \
fcc.o \
fsc.o \
fsc1.o \
fsc2.o \
fsc3.o \
fsc4.o \
fsc5.o \
fsc6.o \
fsc7.o \
fsc8.o \
fax.o \
faxsl.o \
fvf.o \
fvs.o \
vvg.o \
Amptree2q3g_debr.o \
qqb_zbb_gvec.o \
spassign.o \
nagyqqqqg.o

# Master program.
# extra lines: -L$(CRNLIB) -L/home/johnmc/madgraph/lib/ 
#              -ldhelas

mcfm:		\
            $(INTEGRATEFILES) $(LIBFILES) $(NEEDFILES) \
            $(PARTONFILES) $(PHASEFILES) $(SINGLETOPFILES) $(TOPFILES) \
            $(TOPHFILES) $(TOPZFILES) $(TOPGFILES) \
            $(USERFILES) $(VOLFILES) $(WFILES) $(W2JETFILES) $(W2JETVIRT) $(WHBBARFILES) \
            $(WWFILES) $(WZFILES) $(WBBFILES) $(ZFILES) $(ZHBBARFILES) \
            $(ZZFILES) $(ZBBFILES) $(ZGFILES) $(W1JETFILES) $(Z2JETFILES) \
		$(Z1JETFILES) $(HBBBARFILES) $(HWWFILES) $(HZZFILES) \
            $(TAUTAUFILES) $(HTAUTAUFILES) $(SPINORFILES) 
		$(FC) -o $@ \
            $(INTEGRATEFILES) $(LIBFILES) $(NEEDFILES) \
            $(PARTONFILES) $(PHASEFILES) $(SINGLETOPFILES) $(TOPFILES) \
            $(TOPHFILES) $(TOPZFILES) $(TOPGFILES) \
            $(USERFILES) $(VOLFILES) $(WFILES) $(W2JETFILES) $(W2JETVIRT) $(WHBBARFILES) \
            $(WWFILES) $(WZFILES) $(WBBFILES) $(ZFILES) $(ZHBBARFILES) \
            $(ZZFILES) $(ZBBFILES) $(ZGFILES) $(W1JETFILES) $(Z2JETFILES) \
		$(Z1JETFILES) $(HBBBARFILES) $(HWWFILES) $(HZZFILES) \
            $(TAUTAUFILES) $(HTAUTAUFILES) $(SPINORFILES)

	mv mcfm Bin/
# -----------------------------------------------------------------------------
# Specify the dependencies of the .o files and the rules to make them.

FTNCHEKPATH = /home/ellis/Fortran/Ftnchek/ftnchek-2.10.1

.SUFFIXES: .prj

# improved so that .prj files are moved out of src directory and
# into base, only for .f files that don't already exist there
.f.prj: 
		$(FTNCHEKPATH)/ftnchek -project -noextern -library -quiet \
            -include=$(INCPATH) $< ; \
            if ! [ -e $(MYHOME)/$(notdir $<) ] ; then \
            mv $(basename $<).prj $(MYHOME) ; fi
            
PRJS =      $(INTEGRATEFILES:.o=.prj) $(LIBFILES:.o=.prj) \
            $(NEEDFILES:.o=.prj) $(PARTONFILES:.o=.prj) \
            $(PHASEFILES:.o=.prj) $(SINGLETOPFILES:.o=.prj) \
            $(TOPFILES:.o=.prj) $(TOPHFILES:.o=.prj) \
            $(TOPZFILES:.o=.prj) $(TOPGFILES:.o=.prj) \
            $(USERFILES:.o=.prj) $(VOLFILES:.o=.prj) \
            $(WFILES:.o=.prj) $(W2JETFILES:.o=.prj) \
            $(W2JETVIRT:.o=.prj) $(WHBBARFILES:.o=.prj) \
            $(WWFILES:.o=.prj) $(WZFILES:.o=.prj) \
            $(WBBFILES:.o=.prj) $(ZFILES:.o=.prj) \
            $(ZHBBARFILES:.o=.prj) $(ZZFILES:.o=.prj) \
            $(ZBBFILES:.o=.prj) $(ZGFILES:.o=.prj) \
            $(W1JETFILES:.o=.prj) $(Z2JETFILES:.o=.prj) \
		$(Z1JETFILES:.o=.prj) $(HBBBARFILES:.o=.prj) \
            $(HWWFILES:.o=.prj) $(HZZFILES:.o=.prj) \
            $(TAUTAUFILES:.o=.prj) $(HTAUTAUFILES:.o=.prj) \
            $(SPINORFILES:.o=.prj)

check:      $(PRJS) 
		$(FTNCHEKPATH)/ftnchek $(PRJS)

# -----------------------------------------------------------------------------
# Specify other options.

cleanprj:
	- rm -f *.prj $(SOURCEDIR)/*/*.prj

clean:
	- rm -f *.o *.s *.prj *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE
