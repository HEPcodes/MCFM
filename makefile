# Makefile routine.

# Replace this with the location of Cernlib on your system (if desired)
CERNLIB     = 
# Replace this with the location of LHAPDF on your system (if desired)
LHAPDFLIB   = 

MCFMHOME        = /Users/johnmc/MCFM
SOURCEDIR       = /Users/johnmc/MCFM/src
VPATH		= $(DIRS)
BIN		= $(MCFMHOME)/Bin
INCPATH  	= $(SOURCEDIR)/Inc
OUTPUT_OPTION	= -o $(MCFMHOME)/obj/$@

# Set this to NATIVE/PDFLIB/LHAPDF
#   NATIVE -- internal routines
#   PDFLIB -- PDFLIB v8.04
#   LHAPDF -- Les Houches library
PDFROUTINES = NATIVE

# Set this to NO/YES/FROOT
#   NO  -- no n-tuple output or unweighting is possible
#   YES -- n-tuples, unweighting available - but only if CERNLIB exists
#   FROOT -- n-tuples, unweighting available - but only if CERNLIB
#            exists; output uses the FROOT package of P. Nadolsky.
NTUPLES = NO

FC = gfortran
FFLAGS 	= -fno-automatic -fno-f2c -O0 -g -I$(INCPATH)

# If using FROOT package for ROOT ntuples, first specify C++ compiler:
CXX = g++
CXXFLAGS=$(CXXFLAGS0) $(DROOT) 
# ROOTLIBS and ROOTINCLUDE are locations of ROOT libraries and header files.
# Find the ROOT directory automatically by running root-config 
# or specify it manually by editing the variable ROOTDIR 

#ROOTDIR= $(shell root-config --prefix)
#DMYROOT= -DMYROOT
#ROOTLIBS     := $(shell root-config --prefix=$(ROOTDIR)  --libs)
#ROOTINCLUDE     := -I $(shell root-config --prefix=$(ROOTDIR) --incdir)


DIRS	=	$(MCFMHOME):\
		$(MCFMHOME)/obj:\
		$(SOURCEDIR)/User:$(SOURCEDIR)/Vol:\
		$(SOURCEDIR)/Need:$(SOURCEDIR)/Lib:$(SOURCEDIR)/Phase:\
		$(SOURCEDIR)/Parton:$(SOURCEDIR)/Integrate:\
		$(SOURCEDIR)/Wbb:$(SOURCEDIR)/Zbb:\
		$(SOURCEDIR)/WHbbar:$(SOURCEDIR)/ZHbbar:\
		$(SOURCEDIR)/WW:$(SOURCEDIR)/WZ:$(SOURCEDIR)/ZZ:\
		$(SOURCEDIR)/Top:$(SOURCEDIR)/Topdk:$(SOURCEDIR)/Singletop:\
		$(SOURCEDIR)/TopH:$(SOURCEDIR)/TopZ:\
		$(SOURCEDIR)/HWW:$(SOURCEDIR)/HZZ:$(SOURCEDIR)/Tau:\
		$(SOURCEDIR)/Httbar:\
		$(SOURCEDIR)/W:$(SOURCEDIR)/Z:\
		$(SOURCEDIR)/W1jet:$(SOURCEDIR)/Z1jet:\
		$(SOURCEDIR)/Wcjet:$(SOURCEDIR)/Wbjet:\
		$(SOURCEDIR)/W2jet:$(SOURCEDIR)/W2jetvirt:\
		$(SOURCEDIR)/WZbbm:\
		$(SOURCEDIR)/Wgam:$(SOURCEDIR)/Zgam:\
                $(SOURCEDIR)/Z2jet:$(SOURCEDIR)/Zb:\
		$(SOURCEDIR)/bbHiggs:$(SOURCEDIR)/Wt:\
                $(SOURCEDIR)/qqH:$(SOURCEDIR)/qqHWW:$(SOURCEDIR)/ZQjet:\
                $(SOURCEDIR)/ggH:$(SOURCEDIR)/ggHg:\
                $(SOURCEDIR)/ggHggreal:$(SOURCEDIR)/ggHggvirt:\
                $(SOURCEDIR)/H4pCode:\
                $(SOURCEDIR)/WHWW:$(SOURCEDIR)/ZHWW:\
		$(SOURCEDIR)/Stopb:$(SOURCEDIR)/Stopjet:\
		$(SOURCEDIR)/epem3j

# -----------------------------------------------------------------------------
# Specify the object files. 

STOPBFILES = \
qg_tbq.o \
qg_tbq_g.o \
qg_tbq_gs.o \
qg_tbq_gvec.o \
qg_tbq_v.o \
qg_tbq_z.o \
inter.o \
inter_gg.o \
inter_qq.o \
reals1.o \
reals2.o \
reals3.o \
reals_qq.o \
stopf0.o \
stopf1.o \
stop_def.o \
Ammm.o \
Ammp.o \
Ampm.o \
Ampp.o \
Apmm.o \
Apmp.o \
Appm.o \
Appp.o \
Bmmm.o \
Bmmp.o \
Bmpm.o \
Bmpp.o \
Bpmm.o \
Bpmp.o \
Bppm.o \
Bppp.o

STOPJETFILES = \
inters.o \
inters_qq.o \
qq_tbg.o \
qq_tbg_g.o \
qq_tbg_gs.o \
qq_tbg_gvec.o \
qq_tbg_v.o \
qq_tbg_z.o

EPEM3JFILES = \
epem3j.o \
epem3j_g.o \
epem3j_gs.o \
epem3j_gvec.o \
epem3j_v.o \
epem3j_z.o

HWWJETFILES = \
gg_hZZg.o \
gg_hZZgg.o \
gg_hZZg_gs.o \
gg_hZZg_gvec.o \
gg_hZZg_v.o \
gg_hZZg_z.o \
gg_hWWg.o \
gg_hWWgg.o \
gg_hWWg_gs.o \
gg_hWWg_gvec.o \
gg_hWWg_v.o \
gg_hWWg_z.o \
gg_hWWgg_v.o \
gg_hWWgg_gvec.o \
gg_hWWggg.o \
gg_hWWgg_gs.o \
gg_hWWgg_z.o \
gg_hZZgg_v.o \
gg_hZZgg_gvec.o \
gg_hZZggg.o \
gg_hZZgg_gs.o \
gg_hZZgg_z.o \

BBHIGGSFILES = \
bbaqh.o \
bbbbh.o \
bbggh.o \
bbghvirt.o \
qqb_H_gvec.o \
qqb_Hg.o \
qqb_Hg_g.o \
qqb_Hg_gs.o \
qqb_Hg_v.o \
qqb_Hg_z.o

GGHFILES = \
finitemtcorr.o \
gg_h.o \
gg_hg.o \
gg_h_gs.o \
gg_h_gvec.o \
gg_h_v.o \
gg_h_z.o \
hqqgg.o \

GGHGFILES = \
amplonumer.o \
Hqarbsq.o \
Hqaqasq.o \
qqgghn.o \
gg_hgg.o \
gg_hgg_v.o \
gg_hgg_z.o \
gg_hgg_gs.o \
gg_hgg_gvec.o \
gg_hg_gs.o \
gg_hg_gvec.o \
gg_hg_v.o \
gg_hg_z.o \
H4plo.o \
HAQggnew.o \
h4gnew.o \
Amplo_AQgg.o \
Ampsq_AQaq.o \
h4q.o \
h4g.o \
hjetfill.o 

GGHGGrealFILES = \
gg_hggg.o \
gg_ggg.o \
amp_h5g.o \
h5g.o \
Ampppp.o \
Appppp.o \
Ammppp.o \
iperm.o \
h2q3g.o \
a2q3g_mmmpp.o  \
a2q3g_mmpmp.o \
a2q3g_mpppp.o \
a2q3g_mpmmp.o \
h4qg.o \
q4ghppp1.o \
q4ghppp3.o \
q4ghpmp1.o \
q4ghpmp3.o \
extract.o \
gggghn_amp.o \
hqqggdfm.o \
ndveccur.o

GGHGGvirtFILES = \
checkegzres.o \
checkscheme.o \
H4prenorm.o \
Ampvirt_AQgg.o \
Ampvirt_gggg.o \
Ampvirtsq_AQaq.o \
Hggggvsqanal.o \
HAQggvsqanal.o \
Hqarbvsqanal.o \
Hqaqavsqanal.o \
Hqaggvsq.o \
Hggggvsq.o \
Hqarbvsq.o \
Hqaqavsq.o \
GZHggggvsqPoles.o \
GZHqaggvsqPoles.o

H4PCODEFILES = \
A0phiggggpppp.o \
A0phiggggpmmm.o \
A0phiggggmmpp.o \
A0phiggggmpmp.o \
A0phiggggmmmm.o \
A0Hggggpppp.o \
A0Hggggpmmm.o \
A0Hggggmmpp.o \
A0Hggggmpmp.o \
A1phiggggmmmm.o \
A1phiggggmmpp.o \
A1phiggggmpmp.o \
A1Hggggpppp.o \
A1Hggggpmmm.o \
A1Hggggmmpp.o \
A1Hggggmpmp.o \
A0phiAQgg.o \
A0phiAgQg.o \
A0HAQgg.o \
A1phiAQggmpmm.o \
A1phiAQggmpmp.o \
A1phiAQggmppm.o \
A1phiAQggmppp.o \
A1phiAgQgmmpp.o \
A1phiAgQgmppm.o \
A1phiAgQgmppp.o \
A41phiAQgg.o \
A41HAQgg.o \
A43phiAQggmpmm.o \
A43phiAQggmpmp.o \
A43phiAQggmppm.o \
A43HAQgg.o \
A0Hqarb.o \
A0phiqarb.o \
A41Hqarb.o \
A41phiqarb.o \
A42Hqarb.o \
A42phiqarb.o \
Afphiqarbmpmp.o \
Afphiqarbmppm.o \
Alcphiqarbmpmp.o \
Alcphiqarbmppm.o \
Aslcphiqarbmpmp.o \
Aslcphiqarbmppm.o \
BGRL.o \
F31m.o \
F33m.o \
F41m.o \
F41mF.o \
F42me.o \
F42meF.o \
F42mhF.o

HWWFILES = \
qqb_hww.o \
qqb_hww_g.o \
qqb_hww_gvec.o \
qqb_hww_gs.o \
qqb_hww_z.o \
qqb_hww_v.o

HZZFILES = \
qqb_hzz.o \
qqb_hzz_g.o \
qqb_hzz_gvec.o \
qqb_hzz_gs.o \
qqb_hzz_z.o \
qqb_hzz_v.o

HTTBARFILES = \
qqb_higgs.o \
qqb_higgs_odd.o \
ehsv.o \
ehsv_odd.o \
li2.o \

INTEGRATEFILES = \
vegas.o \
ebook.o \
mbook.o \
ran0.o \
ran1.o \
rn.o

LIBFILES = \
ddilog.o \
lenocc.o

NEEDFILES = \
aveptjet.o \
banner.o \
boost.o \
boostx.o \
branch.o \
checkndotp.o \
checkversion.o \
checkjets.o \
checkorder.o \
chooser.o \
ckmfill.o \
coupling.o \
coupling2.o \
couplz.o \
dclaus.o \
dipoles.o \
dipoles_fac.o \
dipoles_mass.o \
dipolesub.o \
dipolesubx.o \
dipolesubxx.o \
donothing_gvec.o \
dot.o \
dotem.o \
etmiss.o \
getbs.o \
getptildejet.o \
getptQ1.o \
gtperp.o \
higgsp.o \
higgsw.o \
histofin.o \
itransform.o \
includedipole.o \
lowint_incldip.o \
masscuts.o \
mcfm.o \
mcfm_exit.o \
mcfm_init.o \
mcfm_vegas.o \
mfrun.o \
ptyrap.o \
r.o \
read_jetcuts.o \
readcoup.o \
reader_input.o \
realhistos.o \
realint.o \
scaleset.o \
sethparams.o \
setmb_msbar.o \
setrunname.o \
smalls.o \
spinork.o \
spinoru.o \
storedip.o \
storeptilde.o \
swapjet.o \
transform.o \
virtint_incldip.o \
writeinfo.o \
writeout.o \
dips_mass.o \
storedip_mass.o \
transform_mass.o \
zeromsq.o

PARTONFILES = \
checkpath.o \
newton1.o

PHASEFILES = \
breitw.o \
breitw1.o \
gen2.o \
gen2a.o \
gen2jet.o \
gen2m.o \
gen3.o \
gen3a.o \
gen3b.o \
gen3jet.o \
gen3m.o \
gen3m_rap.o \
gen3from2.o \
gen4.o \
gen4a.o \
gen4from3.o \
gen4h.o \
gen5.o \
gen5a.o \
gen5from4.o \
gen6.o \
gen6_rap.o \
gen7.o \
gen7m.o \
phase7m_alt.o \
gen7_rap.o \
gen8.o \
gen_njets.o \
gen_soft.o \
genff.o \
genii.o \
genif.o \
genrad.o \
genrff.o \
genrif.o \
genrii.o \
gen_stop.o \
phase2.o \
phase3.o \
phase3m.o \
phase4.o \
phase41.o \
phase4m.o \
phase5.o \
phase5a.o \
phase51.o \
phase6.o \
phase6a.o \
phase6b.o \
phase7.o \
phase7a.o \
phase7b.o \
phase7m.o \
phase8.o \
phi1_2.o \
phi1_2m.o \
phi1_2m_bw.o \
phi1_2m_nobw.o \
phi3m.o \
phi3m0.o \
wt2gen.o \
wt4gen.o \
wtgen.o

QQHFILES = \
qq_Hqq_g.o \
VV_Hqq.o \
VV_Hqq_g.o \
VV_Hqq_gs.o \
VV_Hqq_v.o \
VV_Hqq_z.o \
WW_Hqq.o \
WW_Hqq_g.o \
WW_Hqq_gs.o \
ZZ_Hqq.o \
ZZ_Hqq_g.o \
ZZ_Hqq_gs.o

QQHWWFILES = \
VV_HWW.o \
WW_HWW.o \
ZZ_HWW.o \
VV_HWW_g.o \
WW_HWW_g.o \
ZZ_HWW_g.o \
VV_HWW_gs.o \
WW_HWW_gs.o \
ZZ_HWW_gs.o \
VV_HWW_v.o \
VV_HWW_z.o

SINGLETOPFILES = \
topwidth.o \
extend_trans.o \
qqb_tbb.o \
qqb_tbb_g.o \
qqb_tbb_v.o \
qqb_tbb_z.o \
qqb_tbb_gs.o \
qqb_tbb_gdk.o \
qqb_tbb_gsdk.o \
qqb_tbb_vdk.o \
coefsdk.o \
bq_tpq.o \
bq_tpq_v.o \
bq_tpq_z.o \
bq_tpq_gdk.o \
bq_tpq_gsdk.o \
bq_tpq_vdk.o

TAUTAUFILES = \
qqb_tautau.o \
tautauww.o \
dotks.o \
std.o

TOPFILES = \
qqb_QQb.o \
qqb_QQb_g.o \
qqb_QQb_gs.o \
qqb_QQb_gvec.o \
qqb_QQb_v.o \
qqb_QQb_z.o

TOPDKFILES = \
ggttww1.o \
qqb_QQbdk.o \
qqb_QQbdk_g.o \
qqb_QQbdk_gs.o \
qqb_QQbdk_gvec.o \
qqb_QQbdk_v.o \
qqb_QQbdk_z.o \
KMampsqggQQb.o \
KMampsqqqbQQb.o \
KMampsqggQQbdk.o \
KMampsqqqbQQbdk.o \
KMampsqggQQbdkn.o \
KMgg3h.o \
KMgg2b.o \
KMgg2c1.o \
KMgg2c2.o \
KMgg2c3.o \
KMgg2c4.o \
KMgg2d1.o \
KMgg3f1.o \
KMgg3f2.o \
KMgg2a.o \
KMggtrig.o \
KMggtriq.o \
KMggtriHQ.o \
KMqqb5ab.o \
KMqqb5c.o \
KMqqb5d.o \
KMqqb5e.o \
KMqqb5f.o \
KMqqb5g.o \
KMqqb5h.o \
KMqqb5i.o \
coeffboxa1.o \
coeffboxa2.o \
coeffboxa4.o \
coeffboxa.o \
coeffboxb.o \
Bqq.o \
Bs.o \
Bt.o \
Bu.o \
cdot.o \
pol_mless.o \
psp.o \
spb.o \
bsp.o \
u0spinor.o \
ubar0spinor.o \
ubarspinor.o \
v0spinor.o \
vbar0spinor.o \
vspinor.o \
string.o \
string_KM.o \
ubarWspinor.o \
uWspinor.o \
extend_trans_ttb.o \
polarization.o \
ttbgggppp.o \
ttbgggppm.o \
ttbgggpmm.o \
ttbgggmmm.o \
ttbgggmpp.o \
ttbgggmmp.o \
ttbgggmpm.o \
ttbgggpmp.o \
ttbqqbsqpp.o \
ttbqqbsqpm.o \
ttbqqbsqmp.o \
ttbqqbsqmm.o \
ttbqqbtqpp.o \
ttbqqbtqpm.o \
ttbqqbtqmp.o \
ttbqqbtqmm.o \
ttbqqbqqpp.o \
ttbqqbqqpm.o \
ttbqqbqqmp.o \
ttbqqbqqmm.o \
ttbqqbrqpp.o \
ttbqqbrqpm.o \
ttbqqbrqmp.o \
ttbqqbrqmm.o

TOPHFILES = \
ggtth.o \
qqbtth.o \
qqb_tth.o \
qqb_tottth.o

TOPZFILES = \
ttbbd1.o \
qqb_ttz.o

USERFILES = \
cdfhwwcuts.o \
deltarj.o \
durhamalg.o \
eventhandler.o \
etdoublebin.o \
fill_stdhep.o \
genclust2.o \
genclust_kt.o \
genclust_cone.o \
gencuts.o \
getet.o \
hwwcuts.o \
hwwjetplots.o \
jetlabel_to_stdhep.o \
jetreorder.o \
mdata.o \
miscclust.o \
nplotter.o \
plots_stop_cfmt.o \
stopcuts.o \
wbfcuts.o \
wbfcuts_jeppe.o \
wconstruct.o 

VOLFILES = \
qqb_vol.o \
vol.o \
vol3_mass.o \
vol_mass.o

WFILES = \
qqb_w.o \
qqb_w_g.o \
qqb_w_gs.o \
qqb_w_v.o \
qqb_w_z.o

W1JETFILES = \
A5NLO.o \
A51.o \
A52.o \
qqb_w1jet_gs.o \
qqb_w1jet_v.o \
qqb_w1jet_z.o \
qqb_w_gvec.o \
virt5.o

WBJETFILES = \
qqb_wbjet.o 
 
WCJETFILES = \
qqb_w_cjet.o \
qqb_w_cjet_g.o \
qqb_w_cjet_gs.o \
qqb_w_cjet_gvec.o \
qqb_w_cjet_v.o \
qqb_w_cjet_z.o \
qqb_w_cjet_massless.o \
qqb_w_cjet_massless_g.o \
wqq_sc.o \
w2jetsq_mass.o \
subqcdm.o

W2JETFILES = \
ampwqqqqg.o \
qqb_w2jet.o \
qqb_w2jet_g.o \
qqb_wp2jet_g.o \
qqb_wm2jet_g.o \
qqb_w2jet_gs.o \
qqb_wp2jet_gs.o \
qqb_wm2jet_gs.o \
qqb_w2jet_gvec.o \
qqb_w2jet_gvecx.o \
qqb_w2jet_v.o \
qqb_w2jet_z.o \
qqb_w2jetx.o \
qqbw2j_loop.o \
storecsv_px.o \
storecsv_qx.o \
subqcd.o \
subqcdn.o \
w2jetnx.o \
w2jetsq.o \
xwqqgg_v.o \
xwqqggg.o 


W2JETVIRTFILES = \
a6g.o \
a61g.o \
a6treeg.o \
fax.o \
faxsl.o \
fcc.o \
fcc_qpgmgpqm.o \
fcc_qpgpgmqm.o \
fcc_qpgpgpqm.o \
fcc_qpgpqmgm.o \
fcc_qpgpqmgp.o \
fcc_qpqmgmgp.o \
fcc_qpqmgpgm.o \
fcc_qpqmgpgp.o \
fsc.o \
fsc1.o \
fsc2.o \
fsc3.o \
fsc4.o \
fsc5.o \
fsc6.o \
fsc7.o \
fsc8.o \
fvf.o \
fvs.o \
vvg.o

WHBBARFILES = \
qqb_wh.o \
qqb_wh_g.o \
qqb_wh_gs.o \
qqb_wh_v.o \
qqb_wh_z.o

WHWWFILES = \
qqb_wh_ww.o \
qqb_wh_ww_g.o \
qqb_wh_ww_gs.o \
qqb_wh_ww_v.o

WWFILES = \
BigT.o \
L34_12.o \
a6loop.o \
a7tree.o \
amps_anom.o \
b7tree.o \
fa.o \
qqb_ww.o \
qqb_ww_unpol.o \
susana.o \
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

WZBBMFILES = \
gamps0.o \
gampsabc.o \
gampsdef.o \
gampsgh.o \
mamps.o \
qqb_wbbm.o \
qqb_zbbm.o \
qqb_wbbm_g.o \
wbbgmassa.o \
wbbgmassb.o \
wbbgmassi.o \
wbbgmassf.o

WBBFILES = \
a6.o \
a61.o \
a61LLL.o \
a61LRL.o \
a62.o \
a6routine.o \
a6tree.o \
aqqb_wbb.o \
atrLLL.o \
atrLRL.o \
atree.o \
fpm.o \
fpp.o \
fsl.o \
i3m.o \
lfunctions.o \
lnrat.o \
msqwbb.o \
nagyqqqqg.o \
qqb_wbb.o \
qqb_wbb_g.o \
wbbgamp.o \
qqb_wbb_gs.o \
qqb_wbb_v.o \
qqb_wbb_z.o \
t.o \
t1234.o \
vv.o

WGAMFILES = \
qqb_wgam.o \
qqb_wgam_g.o \
qqb_wgam_gs.o \
qqb_wgam_v.o \
qqb_wgam_z.o \
fagamma.o \
fbgamma.o \
vpole.o

WTFILES = \
BBamps.o \
BBamps_nores.o \
C0fa2m.o \
C0fb2m.o \
I3me.o \
Lsm1_2m.o \
Lsm2_2m.o \
extend_trans_wt.o \
functions.o \
functions1.o \
gs_wc_dg.o \
gs_wt_prog.o \
gs_wt_prog_nores.o \
qb_wtq.o \
qqb_w_tndk.o \
qqb_w_tndk_g.o \
qqb_w_tndk_gs.o \
qqb_w_tndk_gvec.o \
qqb_w_tndk_v.o \
qqb_w_tndk_z.o \
qqb_w_twdk.o \
qqb_w_twdk_g.o \
qqb_w_twdk_gdk.o \
qqb_w_twdk_gs.o \
qqb_w_twdk_gsdk.o \
qqb_w_twdk_gvec.o \
qqb_w_twdk_v.o \
qqb_w_twdk_vdk.o \
qqb_w_twdk_z.o \
qqb_wtbndk.o \
qqb_wtbwdk.o \
tree.o \
virt_mm.o \
virt_mp.o \
virt_pm.o \
virt_pp.o \
vol_wt.o \
wamp.o \
wampd.o \
wcjetn.o \
wtransform_wt.o

ZGAMFILES = \
qqb_zgam.o \
qqb_zgam_g.o \
qqb_zgam_gs.o \
qqb_zgam_v.o \
qqb_zgam_z.o

ZFILES = \
qqb_z.o \
qqb_z1jet.o \
qqb_z_gs.o \
qqb_z_v.o \
qqb_z_z.o

Z1JETFILES = \
qqb_z1jet_gs.o \
qqb_z1jet_v.o \
qqb_z1jet_z.o \
qqb_z_gvec.o

Z2JETFILES = \
a61z.o \
a62z.o \
a63.o \
a63z.o \
a6ax.o \
atreez.o \
fmt.o \
fzip.o \
makem.o \
makemb.o \
msq_ZqqQQg.o \
msq_z2jetx.o \
qqb_z2jet.o \
qqb_z2jet_g.o \
qqb_z2jet_gs.o \
qqb_z2jet_gvec.o \
qqb_z2jet_gvecx.o \
qqb_z2jet_v.o \
qqb_z2jet_z.o \
qqb_z2jetx.o \
storecsz.o \
xzqqgg_v_sym.o \
z2jetsq.o \
z2jetsqn.o

ZHBBARFILES = \
qqb_zh.o \
qqb_zh_g.o \
qqb_zh_gs.o \
qqb_zh_v.o \
qqb_zh_z.o

ZHWWFILES = \
qqb_zh_ww.o \
qqb_zh_ww_g.o \
qqb_zh_ww_gs.o \
qqb_zh_ww_v.o

ZZFILES = \
qqb_zz.o \
qqb_zz_g.o \
qqb_zz_gs.o \
qqb_zz_v.o \
qqb_zz_z.o \
zzamps.o

ZBBFILES = \
amp_qqggg.o \
ampqqb_qqb.o \
aqqb_zbb.o \
msq_qqQQg.o \
qqb_zbb.o \
qqb_zbb_g.o \
qqb_zbb_gs.o \
qqb_zbb_gvec.o \
qqb_zbb_v.o \
qqb_zbb_z.o \
xzqqgg.o \
xzqqgg_v.o \
xzqqggg.o

ZQFILES = \
gQ_zQ.o \
gQ_zQ_g.o \
gQ_zQ_gs.o \
gQ_zQ_v.o \
gQ_zQ_z.o \
qqb_zccm.o

ZQJETFILES = \
genclust_hqrk.o \
msq_ZqqQQg_noid.o \
qqb_zbjet.o \
qqb_zbjet_g.o \
qqb_zbjet_gvec.o \
qqb_zbjet_gs.o \
qqb_zbjet_v.o \
qqb_zbjet_z.o

LIBDIR=.
LIBFLAGS=

# Check NTUPLES flag
ifeq ($(NTUPLES),FROOT)
  ifeq ($(CERNLIB),)
    ERRORMSG=Please specify the path to CERNLIB to use FROOT n-tuples
    $(error $(ERRORMSG))
  endif
  USERFILES += mcfm_froot.o froot.co
  LIBDIR=$(CERNLIB)
  LIBFLAGS = -lmathlib -lpacklib -lkernlib
  LIBFLAGS += $(ROOTLIBS)
  NTUPMSG='   ----> MCFM compiled with FROOT n-tuple output <----'
 else
  ifeq ($(NTUPLES),YES)
    ifeq ($(CERNLIB),)
      ERRORMSG=Please specify the path to CERNLIB to use n-tuples
      $(error $(ERRORMSG))
    endif
    USERFILES += dswhbook.o
    LIBDIR=$(CERNLIB)
# Note: some versions of the CERN libraries (e.g. Redhat v7.2)
#       require that "-lpacklib" be changed to "-lpacklib_noshift"
    LIBFLAGS = -lmathlib -lpacklib -lkernlib
    NTUPMSG='   ----> MCFM compiled with optional n-tuple output <----'
  else
    ifeq ($(NTUPLES),NO)
      USERFILES += dsw_dummy.o
    else
      ERRORMSG=Please set NTUPLES equal to NO/YES/FROOT
      $(error $(ERRORMSG))
    endif
    NTUPMSG = '   ----> MCFM compiled with histogram output only <----'
  endif
endif

# Check PDFROUTINES flag and add appropriate files
ifeq ($(PDFROUTINES),PDFLIB)
   PARTONFILES += \
   alfamz.o \
   fdist_pdflib.o \
   pdfwrap_pdflib.o
   LIBDIR=$(CERNLIB)
   LIBFLAGS := -lpdflib804 $(LIBFLAGS)
   ifeq (,$(findstring packlib,$(LIBFLAGS)))
     LIBFLAGS += -lpacklib -lmathlib
   endif
   PDFMSG='   ----> MCFM compiled with PDFLIB routines <----'
else
ifeq ($(PDFROUTINES),LHAPDF)
   PARTONFILES += \
   alfamz_lhapdf.o \
   fdist_lhapdf.o \
   pdfwrap_lhapdf.o
   LIBDIR += -L$(LHAPDFLIB)
   LIBFLAGS += -lLHAPDF
   PDFMSG='   ----> MCFM compiled with LHAPDF routines <----'
else
ifeq ($(PDFROUTINES),NATIVE)
   PARTONFILES += \
   alfamz.o \
   Ctq4Fn.o \
   Ctq5Par.o \
   Ctq5Pdf.o \
   Cteq6Pdf-2008.o \
   cteq3.o \
   mrs96.o \
   mrs98.o \
   mrs98lo.o \
   mrs98ht.o \
   mrs99.o \
   mrsebh.o \
   mrsg.o \
   mrst2001lo.o \
   jeppelo.o \
   mrst2001.o \
   mrst2002.o \
   mrst2004.o \
   mrst2004f3.o \
   mrst2004f4.o \
   mstwpdf.o \
   mt.o \
   eks98r.o \
   fdist_linux.o \
   pdfwrap_linux.o
   PDFMSG='   ----> MCFM compiled with its own PDFs only <----'
else
   ERRORMSG=Please set PDFROUTINES equal to NATIVE/PDFLIB/LHAPDF
   $(error $(ERRORMSG))
endif
endif
endif

# Master program.
# extra lines: -L$(CRNLIB) -L/home/johnmc/madgraph/lib/ 
#              -ldhelas

OURCODE = $(LIBFILES) $(NEEDFILES) \
          $(PHASEFILES) $(SINGLETOPFILES) \
          $(TOPHFILES) $(TOPZFILES) $(TOPDKFILES) \
          $(USERFILES) $(VOLFILES) $(WFILES) $(W2JETFILES) \
          $(WCJETFILES) $(WBJETFILES) \
	  $(W2JETVIRTFILES) $(WHBBARFILES) $(WGAMFILES) $(ZGAMFILES) \
          $(WWFILES) $(WZFILES) $(ZFILES) $(ZHBBARFILES) \
          $(ZZFILES) $(ZGFILES) $(W1JETFILES) $(Z2JETFILES) \
	    $(Z1JETFILES) $(HWWFILES) $(HZZFILES) \
          $(TAUTAUFILES) $(HTTBARFILES) \
          $(BBHIGGSFILES) $(WBBFILES) $(ZBBFILES) $(WZBBMFILES) \
          $(QQHFILES) $(QQHWWFILES) $(GGHFILES) $(GGHGFILES) \
          $(GGHGGrealFILES) $(GGHGGvirtFILES) $(H4PCODEFILES) \
	  $(TOPFILES) $(ZQFILES) $(ZQJETFILES) $(WTFILES) $(HWWJETFILES) \
	  $(WHWWFILES) $(ZHWWFILES) \
	  $(STOPBFILES) $(STOPJETFILES) $(EPEM3JFILES) \
	  $(CHECKINGFILES)
          
OTHER = $(INTEGRATEFILES) $(PARTONFILES) 
ALLMCFM = $(OTHER) $(OURCODE)

          
# CERNLIB libraries for PDFLIB: -lpdflib804 -lmathlib -lpacklib 

mcfm: $(ALLMCFM)
	$(FC) $(FFLAGS) -L$(LIBDIR) -o $@ \
	$(patsubst %,obj/%,$(ALLMCFM)) $(LIBFLAGS) 
	mv mcfm Bin/
	@echo $(PDFMSG)
	@echo $(NTUPMSG)

# for FROOT package
%.co: %.c
	$(CXX) -c $(CXXFLAGS) $(DMYROOT) $(ROOTINCLUDE) -o obj/$@ $<

# -----------------------------------------------------------------------------
# Specify other options.

FTNCHEKPATH = /home/ellis/Fortran/Ftnchek/ftnchek-3.1.2
FORCHKPATH = /home/ellis/bin/

# Specify the dependencies of the .o files and the rules to make them.

FOROPTS = -include=$(INCPATH) -nonovice -nopretty -quiet

.SUFFIXES: .prj

# improved so that .prj files are moved out of src directory and
# into base, only for .f files that don't already exist there
.f.prj: 
		$(FTNCHEKPATH)/ftnchek -project -noextern\
            $(FOROPTS) $< ; \
            if ! [ -e $(MCFMHOME)/$(notdir $<) ] ; then \
            mv $(basename $<).prj $(MCFMHOME) ; fi
            
PRJS =      $(OURCODE:.o=.prj) 

check:      $(PRJS) 
		$(FTNCHEKPATH)/ftnchek $(FOROPTS) $(PRJS)



PRJSF =      $(OURCODE:.o=.f) 

checkf:      
		$(FORCHKPATH)/forchk -I $(INCPATH) $(PRJSF)

clean:
	- rm -f *.o obj/*.o *.s *.prj *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE

