import os
import argparse
import pyarrow.parquet as pq

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input folder with parquet files", required=True)
parser.add_argument("-o", "--output", help="Output file", required=True)
args = parser.parse_args()

# Check if the input folder exists
if not os.path.exists(args.input):
    raise Exception(f"Input folder {args.input} does not exist")
# Create output folder if it does not exist
os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".parquet")]
schema = pq.ParquetFile(files[0]).schema_arrow

with pq.ParquetWriter(args.output, schema=schema) as writer:
    for i, file in enumerate(files):
        print(f"[{i+1}/{len(files)}] Reading {file}")
        writer.write_table(pq.read_table(file, schema=schema))

print(f"Output file {args.output} created with {len(files)} files")