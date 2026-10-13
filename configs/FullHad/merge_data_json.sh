# Create a single merged file by combining the existing JSON files
python3 << 'EOF'
import json
import glob

# Load all the JSON files created by build_dataset.py
merged = {}

# Load JetHT files
for era in ['B', 'C', 'D', 'E', 'F']:
    try:
        with open(f'datasets/JetHT_2017{era}.json', 'r') as f:
            data = json.load(f)
            merged.update(data)
    except FileNotFoundError:
        print(f"JetHT_2017{era}.json not found")

# Load BTagCSV files
for era in ['B', 'C', 'D', 'E', 'F']:
    try:
        with open(f'datasets/BTagCSV_2017{era}.json', 'r') as f:
            data = json.load(f)
            merged.update(data)
    except FileNotFoundError:
        print(f"BTagCSV_2017{era}.json not found")

# Save merged file
with open('datasets/data_2017.json', 'w') as f:
    json.dump(merged, f, indent=2)

print(f"Merged {len(merged)} datasets into datasets/data_2017.json")
EOF
