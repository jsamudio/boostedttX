from base_conf import cfg, run_options

# Minimal test: one signal + one background, 2017 only
# Use this to verify the workflow runs before launching full jobs
cfg.datasets["filter"]["year"]    = ["2017"]
cfg.datasets["filter"]["samples"] = ["ttHTobb", "TTToHadronic", "QCD_HT"]
