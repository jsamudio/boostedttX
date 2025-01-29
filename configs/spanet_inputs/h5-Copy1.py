import os
from enum import Enum
from collections import defaultdict

import numba
import vector
vector.register_awkward()
vector.register_numba()

import numpy as np
import awkward as ak
import h5py

from omegaconf import OmegaConf

# Inherited from https://github.com/Alexanders101/SPANet/blob/master/spanet/dataset/types.py
class SpecialKey(str, Enum):
    Mask = "MASK"
    Event = "EVENT"
    Inputs = "INPUTS"
    Targets = "TARGETS"
    Particle = "PARTICLE"
    Regressions = "REGRESSIONS"
    Permutations = "PERMUTATIONS"
    Classifications = "CLASSIFICATIONS"
    Embeddings = "EMBEDDINGS"
    Weights = "WEIGHTS"

class H5Dataset:
    def __init__(self, input_file, output_file, cfg, fully_matched=False, shuffle=True, reweigh=False, entrystop=None, has_data=False):
        # Load several input files into a list
        if type(input_file) == str:
            self.input_files = [input_file]
        elif type(input_file) == list:
            self.input_files = input_file
        else:
            raise ValueError(f"Input file {input_file} should be a string or a list of strings.")
        self.output_file = output_file
        self.cfg = cfg
        self.fully_matched = fully_matched
        self.shuffle = shuffle
        self.reweigh = reweigh
        self.entrystop = entrystop
        self.has_data = has_data

        self.sample_dict = defaultdict(dict)

        self.load_config()
        self.check_output()
        self.dataset = Dataset(self)
        if self.fully_matched:
            self.select_fully_matched()

    def get_sample_name(self, input_file):
        '''Get the sample name from the input file name.'''
        sample_list = []
        for sample in self.mapping_sample.keys():
            if sample in input_file:
                sample_list.append(sample)
        if len(sample_list) == 0:
            raise ValueError(f"Sample name not found in the input file name: {input_file}.\nAvailable samples: {self.mapping_sample.keys()}")
        elif len(sample_list) > 1:
            raise ValueError(f"""Multiple sample names found in the input file name: {input_file}.
                                A single sample name should be specified in the file name.\n
                                Available samples: {self.mapping_sample.keys()}""")
        return sample_list[0]

    def build_labels(self, df, sample):
        '''Build labels for the classification nodes.'''
        for s, label in self.mapping_sample.items():
            if s == sample:
                df[label] = ak.values_astype(np.ones(len(df), dtype=int), int)
            else:
                df[label] = ak.values_astype(np.zeros(len(df), dtype=int), int)

        # Define one hot encoded label for multiclassifier
        if self.one_hot_encoding:
            for label, mapping in self.mapping_encoding.items():
                df[label] = ak.values_astype(np.ones(len(df), dtype=int) * mapping[sample], int)
        return df

    def scale_weights(self, df, sample):
        '''Scale the event weights by a factor as specified in the configuration file.'''
        if sample not in self.weights_scale.keys():
            raise ValueError(f"Sample {sample} not found in the weights_scale dictionary.")
        for s, factor in self.weights_scale.items():
            if s == sample:
                # Overwrite the event weights with the scaled weights
                df["event"] = ak.with_field(df.event, factor * df.event.weight, "weight")
        return df

    def load_input(self):
        '''Load the input file.'''
        # Check if the parquet files exist
        for input_file in self.input_files:
            if not os.path.exists(input_file):
                raise ValueError(f"Input file {input_file} does not exist.")
            if not input_file.endswith(".parquet"):
                raise ValueError(f"Input file {input_file} should have the `.parquet` extension.")
        # Read the parquet files
        print(f"Reading {len(self.input_files)} parquet files: ", self.input_files)
        dfs = []
        for input_file in self.input_files:
            print("Reading file: ", input_file)
            df = ak.from_parquet(input_file)
            # For data, save weights of 1 for all events
            if self.has_data & (not "event" in df.fields):
                df["event"] = ak.zip({"weight": np.ones(len(df), dtype=np.float)})
            # Reweigh the events by a factor
            if self.reweigh:
                df = self.scale_weights(df, self.get_sample_name(input_file))
            # Get sample name from the input file name
            df = self.build_labels(df, self.get_sample_name(input_file))
            dfs.append(df)
        # Return the concatenated dataframe
        # If shuffle is True, the events are randomly shuffled
        df_concat = ak.concatenate(dfs)
        if self.shuffle:
            df_concat = df_concat[np.random.permutation(len(df_concat))]
        if self.entrystop:
            df_concat = df_concat[:self.entrystop]
        return df_concat

    def load_config(self):
        '''Load the config file with OmegaConf and read the input features.'''
        self.cfg = OmegaConf.load(self.cfg)
        print("Reading configuration file: ")
        print(OmegaConf.to_yaml(self.cfg))
        self.input_features = self.cfg["input"]
        self.collection = self.cfg["collection"]
        self.targets = self.cfg["particles"]
        self.classification_targets = self.cfg["classification"]
        self.frac_train = self.cfg["frac_train"]
        self.mapping_sample = self.cfg["mapping_sample"]
        self.one_hot_encoding = True if "mapping_encoding" in self.cfg else False
        if self.one_hot_encoding:
            self.mapping_encoding = self.cfg["mapping_encoding"]
        if "weights_scale" in self.cfg:
            self.weights_scale = self.cfg["weights_scale"]

    def check_output(self):
        '''Check the output file extension and if it already exists.'''
        # Check the output file extension
        filename, file_extension = os.path.splitext(self.output_file)
        if not file_extension == ".h5":
            raise ValueError(f"Output file {self.output_file} should be in .h5 format.")
        # Check if output file exists
        if os.path.exists(self.output_file):
            raise ValueError(f"Output file {self.output_file} already exists.")
        os.makedirs(os.path.abspath(os.path.dirname(self.output_file)), exist_ok=True)

    def select_fully_matched(self):
        '''Select only fully matched events.'''
        mask_fullymatched = ak.sum(self.dataset.df[self.collection["Jet"]].matched == True, axis=1) >= 5
        mask_fullymatched = mask_fullymatched & (ak.sum(self.dataset.df[self.collection["FatJet"]].matched == True, axis = 1) == 1)
        #mask_fullymatched = ak.sum(self.dataset.df[self.collection["Jet"]].matched == True, axis=1) >= 6
        df = self.dataset.df[mask_fullymatched]
        jets = df[self.collection["Jet"]]

        # We require exactly 2 jets from the Higgs, 3 jets from the W or hadronic top, and 1 lepton from the leptonic top
        #jets_higgs = jets[jets.prov == 1]
        #mask_match = ak.num(jets_higgs) == 2

        jets_w_thadr = jets[(jets.prov == 5) | (jets.prov == 2)]
        mask_match = mask_match & (ak.num(jets_w_thadr) == 3)

        jets_tlep = jets[jets.prov == 3]
        mask_match = mask_match & (ak.num(jets_tlep) == 1)

        df = df[mask_match]
        print(f"Selected {len(df)} fully matched events")

        self.dataset.df = df

    def create_groups(self):
        '''Create the groups in the h5 file.'''
        for target in self.targets:
            self.file.create_group(f"{SpecialKey.Targets}/{target}")
        for object in self.input_features.keys():
            self.file.create_group(f"{SpecialKey.Inputs}/{object}")
        for group in self.classification_targets:
            self.file.create_group(f"{SpecialKey.Classifications}/{group}")

    def create_targets(self, df):
        jets = df[self.collection["Jet"]]
        fatjets = df[self.collection["FatJet"]]
        indices = ak.local_index(jets)
        print("Jets:", ak.num(jets))
        print("FatJets:", ak.local_index(fatjets))
        print("indices:", indices)

        for target in self.targets:
            # no jets to attribute to the H, so we need to define it as just the FatJet
            '''
            if target == "h":
                mask = jets.prov == 1 # H->b1b2
                # We select the local indices of jets matched with the Higgs
                # The indices are padded with None such that there are 2 entries per event
                # The None values are filled with -1 (a nan value).
                indices_prov = ak.fill_none(ak.pad_none(indices[mask], 2), -1)

                index_b1 = indices_prov[:,0]
                print("index_b1", index_b1)
                index_b2 = indices_prov[:,1]
                print("index_b2", index_b2)

                self.file.create_dataset(f"{SpecialKey.Targets}/h/b1", np.shape(index_b1), dtype='int64', data=index_b1)
                self.file.create_dataset(f"{SpecialKey.Targets}/h/b2", np.shape(index_b2), dtype='int64', data=index_b2)
            '''

            if target == "h":
                mask = fatjets.prov == 1 # H->b1b2
                # We select the local indices of jets matched with the Higgs
                # The indices are padded with None such that there are 2 entries per event
                # The None values are filled with -1 (a nan value).
                indices_prov = ak.fill_none(ak.pad_none(ak.local_index(fatjets), 1), -1)

                index_h = indices_prov[:,0]
                print("index_h", index_h)

                self.file.create_dataset(f"{SpecialKey.Targets}/h/h", np.shape(index_h), dtype='int64', data=index_h)
                
            elif target == "t1":
                mask = jets.prov == 5 # W->q1q2 from t1
                indices_prov = ak.fill_none(ak.pad_none(indices[mask], 2), -1)

                index_q1 = indices_prov[:,0]
                index_q2 = indices_prov[:,1]
                print("index_q1", index_q1)

                mask = jets.prov == 2 # t1->Wb
                index_b_hadr = ak.fill_none(ak.pad_none(indices[mask], 1), -1)[:,0]

                self.file.create_dataset(f"{SpecialKey.Targets}/t1/q1", np.shape(index_q1), dtype='int64', data=index_q1)
                self.file.create_dataset(f"{SpecialKey.Targets}/t1/q2", np.shape(index_q2), dtype='int64', data=index_q2)
                self.file.create_dataset(f"{SpecialKey.Targets}/t1/b", np.shape(index_b_hadr), dtype='int64', data=index_b_hadr)

            elif target == "t2":
                mask = jets.prov == 3 # t2->b
                index_b_lep = ak.fill_none(ak.pad_none(indices[mask], 1), -1)[:,0]

                self.file.create_dataset(f"{SpecialKey.Targets}/t2/b", np.shape(index_b_lep), dtype='int64', data=index_b_lep)
            else:
                raise NotImplementedError

    def create_classifications(self, df):
        '''Create the classification targets in the h5 file.'''
        for group, targets in self.classification_targets.items():
            for target in targets:
                if group == SpecialKey.Event:
                    try:
                        values = df[target]
                    except:
                        raise Exception(f"Target {target} not found in the dataframe.")
                    self.file.create_dataset(f"{SpecialKey.Classifications}/{group}/{target}", np.shape(values), dtype='int64', data=values)
                    print(group, target, values)
                else:
                    raise NotImplementedError

    def create_weights(self, df):
        '''Create the weights in the h5 file.'''
        weights = df.event.weight
        print("Creating dataset: ", f"{SpecialKey.Weights}/weight")
        self.file.create_dataset(f"{SpecialKey.Weights}/weight", np.shape(weights), dtype='float32', data=weights)

    def create_inputs(self, df):
        '''Create the input arrays in the h5 file.'''
        features = self.get_object_features(df)

        for obj, feats in features.items():
            for feat, val in feats.items():
                if feat == "MASK":
                    dtype = 'bool'
                else:
                    dtype = 'float32'
                dataset_name = f"{SpecialKey.Inputs}/{obj}/{feat}"
                print("Creating dataset: ", dataset_name)
                print("Shape is:", np.shape(val))
                ds = self.file.create_dataset(dataset_name, np.shape(val), dtype=dtype, data=val)

    def get_object_features(self, df):

        df_features = defaultdict(dict)
        print(df_features)
        for obj, features in self.input_features.items():

            features_dict = {}
            collection = self.collection[obj]

            if (collection in df.fields):
                objects = df[collection]
                for feat in features:
                    if feat == "MASK":
                        continue
                    if feat in ["sin_phi", "cos_phi"]:
                        phi = objects["phi"]
                        if feat == "sin_phi":
                            values = np.sin(phi)
                        elif feat == "cos_phi":
                            values = np.cos(phi)
                    elif feat == "is_electron":
                        values = ak.values_astype(abs(objects["pdgId"]) == 11, int)
                    else:
                        values = objects[feat]
                    if objects.ndim == 1:
                        features_dict[feat] = ak.to_numpy(values)
                    elif objects.ndim == 2:
                        features_dict[feat] = ak.to_numpy(ak.fill_none(ak.pad_none(values, 16, clip=True), 0))
                    else:
                        raise NotImplementedError

                if "MASK" in features:
                    if not "pt" in features:
                        raise NotImplementedError
                    features_dict["MASK"] = ~(features_dict["pt"] == 0)
            else:
                raise ValueError(f"Collection {collection} not found in the parquet file.")
            df_features[obj] = features_dict
        return df_features

    def h5_tree(self, val, pre=''):
        items = len(val)
        for key, val in val.items():
            items -= 1
            if items == 0:
                # the last item
                if type(val) == h5py._hl.group.Group:
                    print(pre + '└── ' + key)
                    self.h5_tree(val, pre+'    ')
                else:
                    print(pre + '└── ' + key + ' (%d)' % len(val))
            else:
                if type(val) == h5py._hl.group.Group:
                    print(pre + '├── ' + key)
                    self.h5_tree(val, pre+'│   ')
                else:
                    print(pre + '├── ' + key + ' (%d)' % len(val))

    def print(self):
        '''Print the h5 file tree.'''
        self.h5_tree(self.file)

    def save_h5(self, dataset_type):
        '''Save the h5 file.'''
        assert dataset_type in ["train", "test"], f"Dataset type `{dataset_type}` not recognized."

        df = getattr(self.dataset, dataset_type)
        output_file = self.output_file.replace(".h5", f"_{dataset_type}_{len(df)}.h5")
        
        print("Creating output file: ", output_file)
        self.file = h5py.File(output_file, "w")
        self.create_groups()
        if not self.has_data:
            self.create_targets(df)
        self.create_classifications(df)
        self.create_inputs(df)
        self.create_weights(df)
        print(self.file)
        self.print()
        self.file.close()

    def save_h5_all(self):
        '''Save the h5 file for both the training and testing datasets.'''
        for dataset in ["train", "test"]:
            print(f"Processing dataset: {dataset} ({len(getattr(self.dataset, dataset))} events)")
            self.save_h5(dataset)

class Dataset:
    def __init__(self, h5_dataset : H5Dataset):
        self.df = h5_dataset.load_input()
        self.frac_train = h5_dataset.frac_train

    @property
    def train(self):
        '''Return the training dataset according to the fraction `frac_train`.'''
        return self.df[:int(self.frac_train*len(self.df))]

    @property
    def test(self):
        '''Return the testing dataset according to the fraction `frac_train`.'''
        return self.df[int(self.frac_train*len(self.df)):]
