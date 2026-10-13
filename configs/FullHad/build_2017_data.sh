cat > datasets/data_2017_definitions.json << 'EOF'
{
    "JetHT_2017B": {
        "sample": "JetHT",
        "json_output": "datasets/JetHT_2017B.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/JetHT/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "B",
                    "isMC": false
                }
            }
        ]
    },
    "JetHT_2017C": {
        "sample": "JetHT",
        "json_output": "datasets/JetHT_2017C.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/JetHT/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "C",
                    "isMC": false
                }
            }
        ]
    },
    "JetHT_2017D": {
        "sample": "JetHT",
        "json_output": "datasets/JetHT_2017D.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/JetHT/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "D",
                    "isMC": false
                }
            }
        ]
    },
    "JetHT_2017E": {
        "sample": "JetHT",
        "json_output": "datasets/JetHT_2017E.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/JetHT/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "E",
                    "isMC": false
                }
            }
        ]
    },
    "JetHT_2017F": {
        "sample": "JetHT",
        "json_output": "datasets/JetHT_2017F.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/JetHT/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "F",
                    "isMC": false
                }
            }
        ]
    },
    "BTagCSV_2017B": {
        "sample": "BTagCSV",
        "json_output": "datasets/BTagCSV_2017B.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/BTagCSV/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "B",
                    "isMC": false
                }
            }
        ]
    },
    "BTagCSV_2017C": {
        "sample": "BTagCSV",
        "json_output": "datasets/BTagCSV_2017C.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/BTagCSV/Run2017C-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "C",
                    "isMC": false
                }
            }
        ]
    },
    "BTagCSV_2017D": {
        "sample": "BTagCSV",
        "json_output": "datasets/BTagCSV_2017D.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/BTagCSV/Run2017D-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "D",
                    "isMC": false
                }
            }
        ]
    },
    "BTagCSV_2017E": {
        "sample": "BTagCSV",
        "json_output": "datasets/BTagCSV_2017E.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/BTagCSV/Run2017E-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "E",
                    "isMC": false
                }
            }
        ]
    },
    "BTagCSV_2017F": {
        "sample": "BTagCSV",
        "json_output": "datasets/BTagCSV_2017F.json",
        "year": "2017",
        "isMC": false,
        "files": [
            {
                "das_names": ["/BTagCSV/Run2017F-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD"],
                "metadata": {
                    "year": "2017",
                    "era": "F",
                    "isMC": false
                }
            }
        ]
    }
}
EOF


pocket-coffea build-datasets \
  --cfg datasets/data_2017_definitions.json \
  -o \
  -rs 'T[123]_(FR|IT|DE|BE|CH|UK)_\w+'
