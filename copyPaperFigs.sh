#!/bin/sh

# /bin/cp -f figures/paperplots.m3body.pdf               /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_01a.pdf
/bin/cp -f figures/paperplots.unblinded.m3body.logy.pdf  /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_01a.pdf
# /bin/cp -f figures/paperplots.score.pdf                /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_01b.pdf
/bin/cp -f figures/paperplots.score.logy.pdf             /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_01b.pdf

/bin/cp -f figures/paperplots.calo_mt.pdf          /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02a.pdf
/bin/cp -f figures/paperplots.trk_met.pdf          /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02b.pdf
/bin/cp -f figures/paperplots.isolation020.pdf     /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02c.pdf
/bin/cp -f figures/paperplots.ht.pdf               /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02d.pdf
/bin/cp -f figures/paperplots.trk_mt.pdf           /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02e.pdf
/bin/cp -f figures/paperplots.calo_trk_dphi.pdf    /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02f.pdf
/bin/cp -f figures/paperplots.calo_met.pdf         /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02g.pdf
/bin/cp -f figures/paperplots.dptreltrk.pdf        /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_02h.pdf

/bin/cp -f figures/paperplots.calo_dphi3mu.pdf     /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03a.pdf
/bin/cp -f figures/paperplots.pvalue.pdf           /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03b.pdf
/bin/cp -f figures/paperplots.Sa0xy.pdf            /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03c.pdf
/bin/cp -f figures/paperplots.trksfitprob.pdf      /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03d.pdf
/bin/cp -f figures/paperplots.pT3body.pdf          /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03e.pdf
/bin/cp -f figures/paperplots.PVNtrk.pdf           /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03f.pdf
/bin/cp -f figures/paperplots.SLxy.pdf             /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03g.pdf
/bin/cp -f figures/paperplots.dptrelcal.pdf        /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_03h.pdf

/bin/cp -f figures/paperplots.ht_dphimet_calo.pdf  /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_04a.pdf
/bin/cp -f figures/paperplots.ht_dphimet_trk.pdf   /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_04b.pdf
/bin/cp -f figures/paperplots.isolation030.pdf     /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_04c.pdf
/bin/cp -f figures/paperplots.mSS.pdf              /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_04d.pdf
/bin/cp -f figures/paperplots.mOS1.pdf             /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_04e.pdf
/bin/cp -f figures/paperplots.mOS2.pdf             /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_04f.pdf

/bin/cp -f figures/BkgEstimateWithSignal.BDT.1450-1690-1870-2110.neg090.pos1000.pdf  /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_05.pdf
/bin/cp -f figures/BkgEstimateWithSignal.SB0.1450-1690-1870-2110.neg090.pos1000.pdf  /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/fig_06.pdf
/bin/cp -f figures/BkgEstimateWithSignal.*.1450-1690-1870-2110.neg090.pos1000.*      /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/
cd /Users/hod/Google\ Drive/Papers/LFV/paper/paper/figures/
source convert.sh . BkgEstimateWithSignal.BDT.1450-1690-1870-2110.neg090.pos1000
source convert.sh . BkgEstimateWithSignal.SB0.1450-1690-1870-2110.neg090.pos1000
cd -
