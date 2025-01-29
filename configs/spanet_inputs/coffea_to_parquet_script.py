import argparse

from parquet import ParquetDataset

# Read arguments from command line: input file and output directory. Description: script to convert ntuples from coffea file to parquet file.
parser = argparse.ArgumentParser(description='Convert awkward ntuples in coffea files to parquet files.')
parser.add_argument('-c', '--cfg', type=str, required=True, help='YAML configuration file with input features and features to pad')
parser.add_argument('-i', '--input', type=str, required=True, help='Input coffea file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output parquet file')
parser.add_argument('--cat', type=str, default="semilep_LHE", required=False, help='Event category')
parser.add_argument('-n', '--ntuples', default=None, type=str, required=False, help='Additional input parquet file with ntuples')

args = parser.parse_args()

dataset = ParquetDataset(
    args.input,
    args.output,
    args.cfg,
    args.cat,
    args.ntuples
)
dataset.save_parquet()
