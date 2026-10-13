#BASE_CFG=/home/moanwar/ttx_hadronic/boostedttX/configs/FullHad
#BASE_OUT=/cms/data/store/user/moanwar/ttx_output_hadronic/full_run2_v36

# ── 2017 ──
#pocket-coffea run --cfg $BASE_CFG/run_2017_MC.py   -e futures -s 8 -o $BASE_OUT/2017/mc
#pocket-coffea run --cfg $BASE_CFG/run_2017_data.py -e futures -s 8 -o $BASE_OUT/2017/data

# ── 2018 ──
#pocket-coffea run --cfg $BASE_CFG/run_2018_MC.py   -e futures -s 8 -o $BASE_OUT/2018/mc
#pocket-coffea run --cfg $BASE_CFG/run_2018_data.py -e futures -s 8 -o $BASE_OUT/2018/data

# ── 2016 postVFP ──
#pocket-coffea run --cfg $BASE_CFG/run_2016_MC.py   -e futures -s 8 -o $BASE_OUT/2016/mc
#pocket-coffea run --cfg $BASE_CFG/run_2016_data.py -e futures -s 8 -o $BASE_OUT/2016/data

# ── 2016 preVFP (wait for rucio transfers first) ──
#pocket-coffea run --cfg $BASE_CFG/run_2016APV_MC.py  -e futures -s 8 -o $BASE_OUT/2016apv/mc
#pocket-coffea run --cfg $BASE_CFG/run_2016_data.py   -e futures -s 8 -o $BASE_OUT/2016apv/data

# ── Full Run2 ──
#pocket-coffea run --cfg $BASE_CFG/run_all.py -e futures -s 8 -o $BASE_OUT/run2


#!/bin/bash
# run_inference.sh
# Runs pocket-coffea per-sample to avoid OOM during merge phase.
# Each sample runs independently then all outputs are merged at the end.

#!/bin/bash
# run_inference.sh
# Runs pocket-coffea per-sample to avoid OOM during merge phase.
# Requires base_conf.py and merge_coffea.py in the same directory.

BASE_CFG=/home/moanwar/ttx_hadronic/boostedttX/configs/FullHad
BASE_OUT=/cms/data/store/user/moanwar/ttx_output_hadronic/full_run2_v37
SCALEOUT=8

# ── Sample lists ──────────────────────────────────────────────────────────────
#"TTTT"
MC_SAMPLES=(
    "ttHTobb" "ttHToNonbb"
    "TTZToBB" "TTZToQQ"
    "TTWJetsToLNu" "TTWJetsToQQ"
    "TTGJets"
    "TTbb_Hadronic" "TTbb_SemiLeptonic" "TTbb_2L2Nu"
    "TTToHadronic" "TTToSemiLeptonic" "TTTo2L2Nu"
    "ST_tW_top" "ST_tW_antitop" "ST_s_channel"
    "ST_t_channel_top" "ST_t_channel_antitop"
    "QCD_HT"
    "WJetsToLNu_HT" "WJetsToQQ_HT"
    "DYJetsToLL_HT" "ZJetsToQQ_HT"
)



DATA_SAMPLES=("JetHT" "BTagCSV")

FAILED=()

# ── Helper: run one sample ────────────────────────────────────────────────────
run_sample() {
    local SAMPLE=$1
    local YEAR=$2
    local OUTDIR=$3

    local SAMPLE_OUTDIR="${OUTDIR}/_${SAMPLE}"
    local COFFEA_FILE="${SAMPLE_OUTDIR}/output_all.coffea"
    local TMP_CFG="${BASE_CFG}/_tmp_${SAMPLE}_${YEAR}.py"

    # Skip if already done
    if [ -f "$COFFEA_FILE" ]; then
        echo "[SKIP] ${SAMPLE} (${YEAR}) — already done"
        return 0
    fi

    mkdir -p "$SAMPLE_OUTDIR"

    # Write single-sample config — uses base_conf.make_cfg
    cat > "$TMP_CFG" <<PYEOF
from base_conf import make_cfg, run_options
cfg = make_cfg(samples=["${SAMPLE}"], years=["${YEAR}"])
PYEOF

    echo ""
    echo "============================================================"
    echo "Running: ${SAMPLE} (${YEAR})"
    echo "Output:  ${SAMPLE_OUTDIR}"
    echo "============================================================"

    pocket-coffea run \
        --cfg "$TMP_CFG" \
        -e futures \
        -s $SCALEOUT \
        -o "$SAMPLE_OUTDIR"

    local EXIT_CODE=$?
    rm -f "$TMP_CFG"

    if [ $EXIT_CODE -ne 0 ] || [ ! -f "$COFFEA_FILE" ]; then
        echo "[ERROR] ${SAMPLE} (${YEAR}) failed"
        return 1
    fi

    echo "[OK] ${SAMPLE} (${YEAR}) done"
    return 0
}

# ── Helper: merge using external python script ────────────────────────────────
merge_outputs() {
    local OUTDIR=$1
    shift
    local SAMPLES=("$@")

    echo ""
    echo "============================================================"
    echo "Merging ${#SAMPLES[@]} samples in: ${OUTDIR}"
    echo "============================================================"

    # Build comma-separated sample list for Python
    local SAMPLE_LIST=""
    for S in "${SAMPLES[@]}"; do
        SAMPLE_LIST="${SAMPLE_LIST}\"${S}\","
    done

    python3 "$BASE_CFG/merge_coffea.py" \
        --outdir "$OUTDIR" \
        --samples "${SAMPLES[@]}"
}

# ── Run a year+type block ─────────────────────────────────────────────────────
run_block() {
    local YEAR=$1
    local TYPE=$2    # "mc" or "data"
    local OUTDIR=$3

    echo ""
    echo "################################################################"
    echo "# ${YEAR} ${TYPE^^}"
    echo "################################################################"

    mkdir -p "$OUTDIR"

    if [ "$TYPE" = "mc" ]; then
        local SAMPLES=("${MC_SAMPLES[@]}")
    else
        local SAMPLES=("${DATA_SAMPLES[@]}")
    fi

    for SAMPLE in "${SAMPLES[@]}"; do
        run_sample "$SAMPLE" "$YEAR" "$OUTDIR" || FAILED+=("${TYPE}_${YEAR}_${SAMPLE}")
    done

    merge_outputs "$OUTDIR" "${SAMPLES[@]}"
}

# ── Execute all blocks ────────────────────────────────────────────────────────
#run_block "2017" "mc"   "${BASE_OUT}/2017/mc"
#run_block "2017" "data" "${BASE_OUT}/2017/data"

#run_block "2018" "mc"   "${BASE_OUT}/2018/mc"
#run_block "2018" "data" "${BASE_OUT}/2018/data"

#run_block "2016_PostVFP" "mc"   "${BASE_OUT}/2016Post/mc"
#run_block "2016_PostVFP" "data"   "${BASE_OUT}/2016Post/data"

run_block "2016_PreVFP"  "mc"   "${BASE_OUT}/2016Pre/mc"
#run_block "2016_PreVFP"  "data"   "${BASE_OUT}/2016Pre/data"


# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "################################################################"
echo "# SUMMARY"
echo "################################################################"
if [ ${#FAILED[@]} -eq 0 ]; then
    echo "[SUCCESS] All samples completed!"
else
    echo "[WARNING] Failed samples:"
    for F in "${FAILED[@]}"; do
        echo "  $F"
    done
    echo ""
    echo "Rerun this script — completed samples will be skipped."
fi
