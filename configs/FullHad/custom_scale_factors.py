import numpy as np
import awkward as ak
import correctionlib


def sf_btag(params, jets, year, njets, variations=["central"]):
    '''
    DeepJet (or other taggers) AK4 btagging SF.
    See https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/BTV_2018_UL_btagging.html
    The scale factors have 8 default uncertainty
    sources (hf,lf,hfstats1/2,lfstats1/2,cferr1/2) (all of this up_*var*, and down_*var*).
    All except the cferr1/2 uncertainties are to be
    applied to light and b jets. The cferr1/2 uncertainties are to be applied to c jets.
    hf/lfstats1/2 uncertainties are to be decorrelated between years, the others correlated.
    Additional jes-varied scale factors are supplied to be applied for the jes variations.

    if variation is not one of the jes ones both the up and down sf is returned.
    If variation is a jet variation the argument must be up_jes* or down_jes* since it is applied on the specified
    Jes variation jets.
    '''
    bbtagSF = params.jet_scale_factors.bbtagSF[year]
    bbtag_discriminator = 'xbbVsQCD'
    cset = correctionlib.CorrectionSet.from_file(bbtagSF.file)
    corr = cset[btagSF.name]

    flavour = ak.to_numpy(ak.flatten(jets.hadronFlavour))
    abseta = np.abs(ak.to_numpy(ak.flatten(jets.eta)))
    pt = ak.to_numpy(ak.flatten(jets.pt))
    discr = ak.to_numpy(ak.flatten(jets[btag_discriminator]))

    central_SF_byjet = corr.evaluate("central", flavour, abseta, pt, discr)

    def _get_sf_variation_with_mask(variation, mask):
        index = (np.indices(discr.shape)).flatten()[mask]
        # Copying the central SF
        sf = np.copy(central_SF_byjet)
        w = corr.evaluate(variation, flavour[mask], abseta[mask], pt[mask], discr[mask])
        sf[index] = w
        sf_out = ak.prod(ak.unflatten(sf, njets), axis=1)
        return sf_out

    output = {}
    for variation in variations:
        if variation == "central":
            output[variation] = [ak.prod(ak.unflatten(central_SF_byjet, njets), axis=1)]
        else:
            # Nominal sf==1
            nominal = np.ones(ak.num(njets, axis=0))
            # Systematic variations
            if "cferr" in variation:
                # Computing the scale factor only on c-flavour jets
                c_mask = flavour == 4
                output[variation] = [
                    nominal,
                    _get_sf_variation_with_mask(f"up_{variation}", c_mask),
                    _get_sf_variation_with_mask(f"down_{variation}", c_mask),
                ]

            elif variation.startswith("JES") and "AK4" in variation:
                # We need to convert the name of the variation
                # from JES_VariationUp to  up_jesVariation
                if variation.startswith("JES_Total") and variation[-2:] == "Up":
                    btag_jes_var = "up_jes"
                elif variation.startswith("JES_Total") and variation[-4:] == "Down":
                    btag_jes_var = "down_jes"
                else:
                    # we need to remove the possible jet type
                    variation = variation.replace("_AK4PFchs", "")
                    variation = variation.replace("_AK4PFPuppi", "")
                    if variation[-2:] == "Up":
                        btag_jes_var = f"up_jes{variation[4:-2]}"
                    elif variation[-4:] == "Down":
                        btag_jes_var = f"down_jes{variation[4:-4]}"
                # This is a special case where a dedicate btagSF is computed for up and down Jes shape variations.
                # This is not an up/down variation, but a single modified SF.
                # N.B: It is a central SF
                notc_mask = flavour != 4
                output["central"] = [_get_sf_variation_with_mask(btag_jes_var, notc_mask)]
            else:
                # Computing the scale factor only NON c-flavour jets
                notc_mask = flavour != 4
                output[variation] = [
                    nominal,
                    _get_sf_variation_with_mask(f"up_{variation}", notc_mask),
                    _get_sf_variation_with_mask(f"down_{variation}", notc_mask),
                ]

    return output