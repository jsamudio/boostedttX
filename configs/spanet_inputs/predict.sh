MODEL=/cms/data/jsamudio/boosted/boostedttX/configs/spanet_output/version_25
PREDICTION=/cms/data/jsamudio/boosted/boostedttX/configs/spanet_output/version_25/predictions_tthbb_testing.h5
VALIDATION=/cms/data/jsamudio/boosted/boostedttX/configs/spanet_inputs/jan28v2_test_603083.h5

# Create output directory if it does not exist
mkdir -p $(dirname $PREDICTION)
python -m spanet.predict $MODEL $PREDICTION -tf $VALIDATION --gpu

