import uproot
import glob

# 1. Define your path
directory_path = "/cms/data/store/user/ttxeft/NanoAODv15/2024/TTH/*.root"
file_list = glob.glob(directory_path)

total_events = 0
corrupted_files = []

print(f"Scanning {len(file_list)} files...")

# 2. Iterate manually to avoid the uproot batching bug
for i, file_path in enumerate(file_list):
    # Print progress every 500 files so you know it hasn't frozen
    if i % 500 == 0 and i > 0:
        print(f"Processed {i}/{len(file_list)} files...")
        
    try:
        # Open each file explicitly and grab the event count
        with uproot.open(file_path) as root_file:
            total_events += root_file["Events"].num_entries
            
    except Exception:
        # If a file is completely broken/corrupted, it skips it and logs it here
        corrupted_files.append(file_path)

print(f"\n✅ Grand Total Events: {total_events}")

# 3. Report any broken files
if corrupted_files:
    print(f"⚠️ Warning: Found {len(corrupted_files)} corrupted or unreadable files.")
    # print("Corrupted files:", corrupted_files) # Uncomment to see the bad files
