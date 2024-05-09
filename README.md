# boostedttX
Contains pocket-coffea version of the boosted ttX analysis

# Run analysis on `lxplus`
Copy repo preferably onto EOS area for storage concerns. Then activate proxy to access datasets.

```bash
voms-proxy-init --voms cms
```

No need to install anything, just run in an apptainer with:

```bash
apptainer shell --bind /afs -B /cvmfs/cms.cern.ch \
                --bind /tmp  --bind /eos/cms/ \
    --env KRB5CCNAME=$KRB5CCNAME --bind /etc/sysconfig/ngbauth-submit  \
    /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-stable
```
### Build datasets using European sources:

```bash
build_datasets.py --cfg datasets/datasets_definitions.json -o -rs 'T[123]_(FR|IT|DE|BE|UK)_\w+'
```

### Considerations for DNN

To accommodate the DNN processing, the user needs to make a local virtual environment inside the singularity container:

```bash
python -m venv --system-site-packages myenv
source myenv/bin/activate
python -m pip install tensorflow keras
```

In addition, for the processor `runner.py` to function with the custom packages, I found it best to also have a local copy of PocketCoffea.


### Run local test of processor using cloned repo:

```bash
python PocketCoffea/pocket_coffea/scripts/runner.py --cfg testconf.py --test -o test
```

For a full run:

```bash
python PocketCoffea/pocket_coffea/scripts/runner.py --cfg testconf.py --executor futures --full -s 10 -o test
```

### Using output columns, train DNN:

```bash
python dnn_datasets -i test/output_all.coffea

python dnn_model.py
```

The analysis will need to be rerun to apply the trained DNN model.
There might be an issue with missing `.h5` model file on inital run, if so I would suggest modifying the workflow for the first run to remove the DNN application.

### Make test plots:
```bash
make_plots.py -i output_all.coffea --cfg parameters_dump.yaml -o plots --overwrite -op params/plotting_style.yaml
```
