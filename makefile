# Makefile routine.

# Replace this with the location of Cernlib on your system (if desired)
CERNLIB     = 
# Replace this with the location of LHAPDF on your system (if desired)
LHAPDFLIB   = 

MCFMHOME        = /home/johnmc/MCFM
SOURCEDIR       = /home/johnmc/MCFM/src
VPATH		= $(DIRS)
BIN		= $(MCFMHOME)/Bin
INCPATH  	= $(SOURCEDIR)/Inc
OUTPUT_OPTION	= -o $(MCFMHOME)/obj/$@

# Set this to NATIVE/PDFLIB/LHAPDF
#   NATIVE -- internal routines
#   PDFLIB -- PDFLIB v8.04
#   LHAPDF -- Les Houches library
PDFROUTINES = NATIVE

DIRS	=	$(MCFMHOME):\
		$(MCFMHOME)/obj:\
		$(SOURCEDIR)/User:\
		$(SOURCEDIR)/Vol:\
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
		$(SOURCEDIR)/WZbbm:$(SOURCEDIR)/Wgam:\
		$(SOURCEDIR)/W2jet/Mad:\
		$(SOURCEDIR)/Z2jet:$(SOURCEDIR)/Z2jet_mad:\
		$(SOURCEDIR)/Z2jet/Mad:$(SOURCEDIR)/madtest:$(MCFMHOME)/madwm2:\
		$(SOURCEDIR)/bbHiggs:$(SOURCEDIR)/bbHiggs_mad:\
		$(SOURCEDIR)/Topdk:\
		$(SOURCEDIR)/madtest/w2jet:\
		$(SOURCEDIR)/madtest/z2jet

FC = g77
FFLAGS 	= -c -fno-f2c -malign-double -O2 -I$(INCPATH) 

# -----------------------------------------------------------------------------
# Specify the object files. 

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

HBBBARFILES = \
qqb_hbbbar.o \
qqb_hbbbar_g.o \
qqb_hbbbar_gs.o \
qqb_hbbbar_gvec.o \
qqb_hbbbar_v.o \
qqb_hbbbar_z.o

HTTBARFILES = \
qqb_higgs.o \
qqb_higgs_odd.o \
ehsv.o \
ehsv_odd.o \
li2.o \
Htautau.o

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
checkversion.o \
chooser.o \
ckmfill.o \
clust.o \
coll.o \
coll3.o \
coll3m.o \
coll4.o \
coll4a.o \
coll5.o \
conserve.o \
coupling.o \
couplz.o \
ddvdif.o \
dipoles.o \
dipoles_fac.o \
dipolesub.o \
dipolesubx.o \
dips_mass.o \
dittdrein.o \
donothing_gvec.o \
dot.o \
dotem.o \
dotpr.o \
dsigdy.o \
etmiss.o \
gasdev.o \
gtperp.o \
higgsp.o \
higgsw.o \
histofin.o \
itransform.o \
lowint.o \
masscuts.o \
mcfm.o \
mcfm_exit.o \
mcfm_init.o \
mcfm_vegas.o \
mfrun.o \
preclus.o \
prod4.o \
ptyrap.o \
r.o \
readcoup.o \
reader.o \
realint.o \
scaleset.o \
singcheck.o \
singgen.o \
smalls.o \
spinork.o \
spinoru.o \
storedip.o \
storedip_mass.o \
storeptilde.o \
swapjet.o \
transform.o \
transform_mass.o \
virtint.o \
wconstruct.o \
writeout.o

PARTONFILES = \
alfamz.o \
newton1.o

PHASEFILES = \
breitw.o \
breitw1.o \
gen2.o \
gen2a.o \
gen2m.o \
gen3.o \
gen3a.o \
gen3b.o \
gen3m.o \
gen3m_rap.o \
gen3from2.o \
gen4.o \
gen4a.o \
gen4from3.o \
gen5.o \
gen5a.o \
gen5from4.o \
gen6.o \
gen6_rap.o \
gen7.o \
gen7_rap.o \
gen_njets.o \
gencol.o \
genff.o \
genif.o \
genii.o \
genrad.o \
genrff.o \
genrif.o \
genrii.o \
kingen.o \
phase2.o \
phase3.o \
phase3m.o \
phase4.o \
phase41.o \
phase4m.o \
phase5.o \
phase51.o \
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
wt2gen.o \
wt4gen.o \
wtgen.o

SINGLETOPFILES = \
qg_tbb.o \
qqb_tbb.o

TAUTAUFILES = \
qqb_tautau.o \
tautauww.o

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
dotks.o \
std.o \
ttgdiags.o

USERFILES = \
deltarj.o \
etdoublebin.o \
genclust.o \
genclust2.o \
genclust_kt.o \
genclust_cone.o \
gencuts.o \
mdata.o \
nplotter.o \
threebee.o

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
a5nlo.o \
qqb_w1jet_gs.o \
qqb_w1jet_soft.o \
qqb_w1jet_v.o \
qqb_w1jet_z.o \
qqb_w_gvec.o \
virt5.o

W2JETFILES = \
aqqb_wgg.o \
bit.o \
couplant.o \
initqqgg.o \
initqqqq.o \
msq_WqqQQg.o \
nagy1.o \
nagy2.o \
qqb_w2jet.o \
qqb_w2jet_g.o \
qqb_wp2jet_g.o \
qqb_wm2jet_g.o \
qqb_w2jet_gs.o \
qqb_wp2jet_gs.o \
qqb_wm2jet_gs.o \
qqb_w2jet_gvec.o \
qqb_w2jet_gvecx.o \
qqb_w2jet_soft.o \
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
wmakem.o \
wmakemb.o \
work.o \
xmqqgg.o \
xmqqqq.o \
xwqqbqqb.o \
xwqqgg_v.o \
xwqqggg.o \
xwqqqq.o \
xzqqqq.o

W2JETVIRTFILES = \
a61g.o \
a6g.o \
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

WWFILES = \
BigT.o \
L34_12.o \
a6loop.o \
a7tree.o \
amps_anom.o \
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

WZBBMFILES = \
gamps0.o \
gampsabc_old.o \
gampsdef_old.o \
gampsgh_old.o \
mamps.o \
qqb_wbbm.o \
qqb_zbbm.o

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
nagyqqqqg.o \
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

WGAMFILES = \
qqb_wgam.o \
qqb_wgam_g.o \
qqb_wgam_gs.o \
qqb_wgam_v.o \
qqb_wgam_z.o \
fagamma.o \
fbgamma.o \
vpole.o

ZFILES = \
qqb_z.o \
qqb_z1jet.o \
qqb_z_g.o \
qqb_z_gs.o \
qqb_z_v.o \
qqb_z_z.o

Z1JETFILES = \
qqb_z1jet_gs.o \
qqb_z1jet_soft.o \
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

ZZFILES = \
qqb_zz.o \
qqb_zz_g.o \
qqb_zz_gs.o \
qqb_zz_v.o \
qqb_zz_z.o \
zzamps.o

ZBBFILES = \
amp_qqgg.o \
amp_qqggg.o \
ampqqb_qqb.o \
aqqb_zbb.o \
msq_qqQQg.o \
qqb_zbb.o \
qqb_zbb_g.o \
qqb_zbb_gs.o \
qqb_zbb_gvec.o \
qqb_zbb_soft.o \
qqb_zbb_v.o \
qqb_zbb_z.o \
spassign.o \
xzqqgg.o \
xzqqgg_v.o \
xzqqggg.o

BBHIGGSFILES = \
bbaqh.o \
bbbbh.o \
bbggh.o \
bbghvirt.o \
gq_qqqb.o \
qqb_H_gvec.o \
qqb_Hg.o \
qqb_Hg_g.o \
qqb_Hg_gs.o \
qqb_Hg_v.o \
qqb_Hg_z.o

# Check PDFLIB flag and add appropriate files
ifeq ($(PDFROUTINES),PDFLIB)
   PARTONFILES += \
   fdist_pdflib.o \
   pdfwrap_pdflib.o
   LIBDIR=$(CERNLIB)
   LIBFLAGS= -lpdflib804 -lmathlib -lpacklib
   MSG='   ----> MCFM compiled with PDFLIB routines <----'
else
ifeq ($(PDFROUTINES),LHAPDF)
   PARTONFILES += \
   fdist_lhapdf.o \
   pdfwrap_lhapdf.o
   LIBDIR=$(LHAPDFLIB)
   LIBFLAGS=-lLHAPDF
   MSG='   ----> MCFM compiled with LHAPDF routines <----'
else
ifeq ($(PDFROUTINES),NATIVE)
   PARTONFILES += \
   Ctq4Fn.o \
   Ctq5Pdf.o \
   Ctq6Pdf.o \
   cteq3.o \
   mrs96.o \
   mrs98.o \
   mrs98ht.o \
   mrs99.o \
   mrsebh.o \
   mrsg.o \
   mrst2001.o \
   mt.o \
   fdist_linux.o \
   pdfwrap_linux.o
   LIBDIR=.
   LIBFLAGS=
   MSG='   ----> MCFM compiled with its own PDFs only <----'
else
   ERRORMSG=Please set PDFROUTINES equal to NATIVE/PDFLIB/LHAPDF
   $(error $(ERRORMSG))
endif
endif
endif

# Master program.

ALLMCFM = $(INTEGRATEFILES) $(LIBFILES) $(NEEDFILES) \
          $(PARTONFILES) $(PHASEFILES) $(SINGLETOPFILES) \
          $(TOPHFILES) $(TOPZFILES) $(TOPGFILES) \
          $(USERFILES) $(VOLFILES) $(WFILES) $(W2JETFILES) \
	  $(W2JETVIRTFILES) $(WHBBARFILES) $(WGAMFILES) \
          $(WWFILES) $(WZFILES) $(ZFILES) $(ZHBBARFILES) \
          $(ZZFILES) $(ZGFILES) $(W1JETFILES) $(Z2JETFILES) \
	  $(Z1JETFILES) $(HBBBARFILES) $(HWWFILES) $(HZZFILES) \
          $(TAUTAUFILES) $(HTTBARFILES) \
          $(BBHIGGSFILES) $(WBBFILES) $(ZBBFILES) $(WZBBMFILES)

mcfm: $(ALLMCFM)
	$(FC) -L$(LIBDIR) -o $@ \
	$(patsubst %,obj/%,$(ALLMCFM)) $(LIBFLAGS)
	mv mcfm Bin/mcfm
	@echo $(MSG)

# -----------------------------------------------------------------------------
# Specify other options.

clean:
	- rm -f *.o obj/*.o *.s *.prj *~ core

# -----------------------------------------------------------------------------

# DO NOT DELETE

