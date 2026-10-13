from coffea.util import load
import numpy as np

out = load('/cms/data/store/user/moanwar/ttx_output_hadronic/out_onnx_v17/output_all.coffea')
sum_signOf_gw = dict(out['sum_signOf_genweights'])
cols = out['columns']

QCD_CUT  = 0.0004
CATEGORY = 'btag_mask'

samples = {
    'ttH->bb': ['ttHTobb__genMatch', 'ttHTobb__non_genMatch'],
    'TTZ->bb': ['TTZToBB__genMatch', 'TTZToBB__non_genMatch'],
    'tt+B':    ['TTbb_Hadronic__tt+B',  'TTToHadronic__tt+B'],
    'tt+C':    ['TTToHadronic__tt+C',   'TTbb_Hadronic__tt+C'],
    'tt+LF':   ['TTToHadronic__tt+LF',  'TTbb_Hadronic__tt+LF'],
    'QCD':     ['QCD_HT'],
}

def get_val(col):
    return col.value if hasattr(col, 'value') else np.array(col)

header = (
    f"{'Sample':<12} | {'Raw_total':>10} | {'Raw_afterQCD':>12} | "
    f"{'Wgt_total':>12} | {'Wgt_afterQCD':>12} | {'MeanWgt':>12} | {'QCD_eff':>8}"
)
print(header)
print('-' * 90)

for label, skeys in samples.items():
    raw_total = 0
    raw_qcd   = 0
    wgt_total = 0.0
    wgt_qcd   = 0.0
    mean_wgts = []

    for skey in skeys:
        if skey not in cols:
            continue
        for dataset in cols[skey]:
            if CATEGORY not in cols[skey][dataset]:
                continue
            if dataset not in sum_signOf_gw:
                continue
            cat = cols[skey][dataset][CATEGORY]
            c   = cat['nominal'] if 'nominal' in cat else cat
            if 'spanet_output_qcd' not in c:
                continue

            qcd  = get_val(c['spanet_output_qcd'])
            norm = get_val(c['events_norm_weight'])
            genw = get_val(c['events_genWeight'])
            wgts = (norm / sum_signOf_gw[dataset]) * np.sign(genw)

            mask = qcd < QCD_CUT

            raw_total += len(qcd)
            raw_qcd   += int(mask.sum())
            wgt_total += float(wgts.sum())
            wgt_qcd   += float(wgts[mask].sum())
            if mask.sum() > 0:
                mean_wgts.extend(wgts[mask].tolist())

    mean_w  = float(np.mean(mean_wgts)) if mean_wgts else 0.0
    qcd_eff = raw_qcd / raw_total if raw_total > 0 else 0.0

    row = (
        f"{label:<12} | {raw_total:>10d} | {raw_qcd:>12d} | "
        f"{wgt_total:>12.2f} | {wgt_qcd:>12.2f} | {mean_w:>12.6f} | {qcd_eff:>8.4f}"
    )
    print(row)

print()
print("QCD_eff = fraction of raw events passing QCD cut")
print("Wgt_afterQCD = sum of event weights after QCD cut (physics yield)")
