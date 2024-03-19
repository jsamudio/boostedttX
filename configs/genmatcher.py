import math
import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

def match_gen_lep(events):
    lep_eta = events.LeptonGood.eta
    lep_phi = events.LeptonGood.phi
    #
    gen_id = events.GenPart.pdgId
    print(gen_id)
