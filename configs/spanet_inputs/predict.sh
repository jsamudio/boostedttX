MODEL=/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/spanet_output/version_0
PREDICTION=/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/spanet_output/version_0/predictions_tthbb_testing.h5
VALIDATION=/cms/data/jsamudio/boosted/boostedttX/configs/spanet_inputs/uscms_test_630869.h5

# Create output directory if it does not exist
mkdir -p $(dirname $PREDICTION)
python -m spanet.predict $MODEL $PREDICTION -tf $VALIDATION --gpu

