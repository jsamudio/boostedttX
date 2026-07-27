import math
import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

def get_controlvars(events):
    ak4 = events.JetGood
    ak8 = events.FatJetGood
    bjet = events.bJetGood
    qjet = events.qJetGood
    lep = events.LeptonGood
    ele = events.ElectronGood
    muon = events.MuonGood
    met = events.PuppiMET

    def sortPt(obj):
        # 1. Get the indices that would sort the jets by pT in descending order
        sort_indices = ak.argsort(obj.pt, ascending=False)

        # 2. Apply those indices to the Jet collection
        sorted_obj = obj[sort_indices]

        return sorted_obj

    # leading is 1 and subleading is 2

    events['nPV'] = events.PV.npvs
    # print("nPV:", events['nPV'])
    
    events['nPVGood'] = events.PV.npvsGood
    # print("nPVGood:", events['nPVGood'])
    
    events['MET_pt'] = met.pt
    # print("MET_pt:", events['MET_pt'])
    
    events['MET_phi'] = met.phi
    # print("MET_phi:", events['MET_phi'])
    
    events['lep_pt'] = ak.firsts(lep.pt)
    # print("lep_pt:", events['lep_pt'])
    
    events['ele_pt'] = ak.firsts(ele.pt)
    # print("ele_pt:", events['ele_pt'])
    
    events['muon_pt'] = ak.firsts(muon.pt)
    # print("muon_pt:", events['muon_pt'])
    
    events['lep_eta'] = ak.firsts(lep.eta)
    # print("lep_eta:", events['lep_eta'])
    
    events['ele_eta'] = ak.firsts(ele.eta)
    # print("ele_eta:", events['ele_eta'])
    
    events['muon_eta'] = ak.firsts(muon.eta)
    # print("muon_eta:", events['muon_eta'])
    
    events['n_ak4'] = ak.count(ak4.pt, axis=1)
    # print("n_ak4:", events['n_ak4'])
    
    events['n_bjet'] = ak.count(bjet.pt, axis=1)
    # print("n_bjet:", events['n_bjet'])
    
    events['n_ak8'] = ak.count(ak8.pt, axis=1)
    # print("n_ak8:", events['n_ak8'])
    
    events['jet1_pt'] = sortPt(ak4)[:,0].pt
    # print("jet1_pt:", events['jet1_pt'])
    
    events['jet2_pt'] = sortPt(ak4)[:,1].pt
    # print("jet2_pt:", events['jet2_pt'])
    
    events['bjet1_pt'] = sortPt(bjet)[:,0].pt
    # print("bjet1_pt:", events['bjet1_pt'])
    
    #events['bjet2_pt'] = sortPt(bjet)[:,1].pt
    # print("bjet2_pt:", events['bjet2_pt'])
    
    events['jet1_eta'] = sortPt(ak4)[:,0].eta
    # print("jet1_eta:", events['jet1_eta'])
    
    events['jet2_eta'] = sortPt(ak4)[:,1].eta
    # print("jet2_eta:", events['jet2_eta'])
    
    events['bjet1_eta'] = sortPt(bjet)[:,0].eta
    # print("bjet1_eta:", events['bjet1_eta'])
    
    #events['bjet2_eta'] = sortPt(bjet)[:,1].eta
    # print("bjet2_eta:", events['bjet2_eta'])
    
    events['jet1_btag'] = sortPt(ak4)[:,0].btagB
    # print("jet1_btag:", events['jet1_btag'])
    
    events['jet2_btag'] = sortPt(ak4)[:,1].btagB
    # print("jet2_btag:", events['jet2_btag'])
    
    events['bjet1_btag'] = sortPt(bjet)[:,0].btagB
    # print("bjet1_btag:", events['bjet1_btag'])
    
    # FIXME
    #events['bjet2_btag'] = sortPt(bjet)[:,1].btagB 
    # print("bjet2_btag:", events['bjet2_btag'])
    
    events['fatjet1_pt'] = sortPt(ak8)[:,0].pt
    # print("fatjet1_pt:", events['fatjet1_pt'])
    
    events['fatjet1_eta'] = sortPt(ak8)[:,0].eta
    # print("fatjet1_eta:", events['fatjet1_eta'])
    
    events['fatjet1_mass'] = (sortPt(ak8)[:,0].globalParT3_massCorrX2p * sortPt(ak8)[:,0].mass * (1-sortPt(ak8)[:,0].rawFactor))
    # print("fatjet1_mass:", events['fatjet1_mass'])

    