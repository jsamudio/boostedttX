"""ZH candidate variables of the Run 3 analysis (candidate = leading fatjet in GloParT XbbVsQCD)."""
import awkward as ak

deltaR = (lambda obj1_, obj2_: ak.flatten(obj1_.metric_table(obj2_)))


def zh_helper(events):
    ZHCand = events.FatJetSorted[:, 0]
    ak4 = events.JetGood
    bjet = events.bJetGood
    qjet = events.qJetGood

    events['ZH_pt'] = ZHCand.pt
    events['ZH_M'] = ZHCand.mass                      # GloParT 2-prong regressed mass (set in the workflow)
    events['MET_pt'] = events.PuppiMET.pt

    # Taggers of the candidate
    events['ZH_bbvLscore'] = ZHCand.btagBB            # GloParT XbbVsQCD
    events['ZH_GloParT3_Xbb'] = ZHCand.globalParT3_Xbb
    events['ZH_GloParT3_QCD'] = ZHCand.globalParT3_QCD
    events['ZH_PNet_XbbVsQCD'] = ZHCand.particleNet_XbbVsQCD
    events['ZH_PNetLegacy_XbbVsQCD'] = ZHCand.particleNetLegacy_Xbb / (ZHCand.particleNetLegacy_Xbb + ZHCand.particleNetLegacy_QCD)

    # AK4 jets inside / outside the candidate cone
    ZH_b_dr = deltaR(ZHCand, bjet)
    ZH_q_dr = deltaR(ZHCand, qjet)
    events['n_b_inZH'] = ak.sum(ZH_b_dr <= 0.8, axis=1)
    events['n_b_outZH'] = ak.sum(ZH_b_dr > 0.8, axis=1)
    events['n_q_inZH'] = ak.sum(ZH_q_dr <= 0.8, axis=1)
    events['n_q_outZH'] = ak.sum(ZH_q_dr > 0.8, axis=1)
    events['n_ak4_outZH'] = ak.sum(deltaR(ZHCand, ak4) > 0.8, axis=1)
    events['n_ak4jets'] = ak.count(ak4.pt, axis=1)
