import json
import glob
import os

for filepath in glob.glob("datasets/*.json"):
    with open(filepath) as f:
        content = f.read()
    
    # Replace the xrootd redirector path with local path
    content = content.replace(
        "root://cms-xrdr.sdfarm.kr:1094//xrd/store/mc",
        "file:/cms/data/store/mc"
    )
    # Also catch any other redirector variations
    content = content.replace(
        "root://cmsxrootd.fnal.gov//store/mc",
        "file:/cms/data/store/mc"
    )
    content = content.replace(
        "root://xrootd.unl.edu//store/mc",
        "file:/cms/data/store/mc"
    )
    
    with open(filepath, "w") as f:
        f.write(content)

print("Done! All JSON files updated.")
