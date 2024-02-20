# boostedttX
Contains pocket-coffea version of the boosted ttX analysis

# Run analysis on `lxplus`
Copy repo preferably onto EOS area for storage concerns

No need to install anything, just run in an apptainer with:

```bash
apptainer shell --bind /afs -B /cvmfs/cms.cern.ch \
                --bind /tmp  --bind /eos/cms/ \
    --env KRB5CCNAME=$KRB5CCNAME --bind /etc/sysconfig/ngbauth-submit  \
    /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-stable
```
### Build datasets using European sources:

```bash
build_datasets.py --cfg datasets/datasets_definitions.json -o -rs 'T[123]_(FR|IT|DE|BE|CH|UK)_\w+'
```

### Run local test of processor:

```bash
runner.py --cfg testconf.py --test -o .
```

### Make test plots:
```bash
make_plots.py -i output_all.coffea --cfg parameters_dump.yaml -o plots --overwrite
```
