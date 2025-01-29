import os
from collections import defaultdict

import numba
import vector
vector.register_awkward()
vector.register_numba()

import numpy as np
import awkward as ak

from omegaconf import OmegaConf

from coffea.util import load
from coffea.processor.accumulator import column_accumulator
from coffea.processor import accumulate

class ParquetDataset:
    def __init__(self, input_file, output_file, cfg, cat="btag_mask", input_ntuples=None):
        self.input_file = input_file
        self.output_file = output_file
        self.cfg = cfg
        self.cat = cat
        self.input_ntuples = input_ntuples

        self.load_config()
        self.load_input()
        self.check_output()
        self.normalize_genweights()
        self.load_features()
        if self.schema == "coffea":
            self.build_arrays_from_accumulator()
        elif self.schema == "parquet":
            self.build_arrays_from_ntuples()
        else:
            raise ValueError(f"Schema {self.schema} not recognized.")

    @property
    def nevents(self):
        if self.schema == "coffea":
            pass
        elif self.schema == "parquet":
            return len(self.events)

    def get_datasets(self, sample=None):
        if self.schema == "coffea":
            datasets = list(self.columns[sample].keys())
        elif self.schema == "parquet":
            datasets_all = [dataset for dataset, d in self.cutflow.items() if d.get(sample, None) != None]
            datasets = []
            for dataset in datasets_all:
                for sample, nevt in self.cutflow[dataset].items():
                    if nevt == self.nevents:
                        datasets.append(dataset)
        else:
            raise ValueError(f"Schema {self.schema} not recognized.")

        if len(datasets) == 0:
            raise ValueError(f"No dataset found for sample {sample}.")
        if len(datasets) > 1:
            raise ValueError(f"Multiple datasets found for sample {sample}.\nDatasets: {datasets}")

        return datasets

    @property
    def samples(self):
        if self.schema == "coffea":
            return list(self.columns.keys())
        elif self.schema == "parquet":
            samples = []
            datasets = list(self.cutflow.keys())
            for dataset in datasets:
                for sample, nevt in self.cutflow[dataset].items():
                    if nevt == self.nevents:
                        samples.append(sample)
            if len(samples) == 0:
                raise ValueError(f"No sample found with {self.nevents} events.")
            return samples

    def is_mc(self, sample):
        if sample.startswith("DATA"):
            return False
        else:
            return True

    def load_input(self):
        '''Load the input file and check if the event category `cat` is present.'''
        if not os.path.exists(self.input_file):
            raise ValueError(f"Input file {self.input_file} does not exist.")
        else:
            print(f"Reading input file: {self.input_file}")
            self.df = load(self.input_file)
            self.columns = self.df["columns"]
            if not self.cat in self.df["cutflow"].keys():
                raise ValueError(f"Event category `{self.cat}` not found in the input file.")
            self.sum_genweights = self.df["sum_genweights"]
            self.cutflow = self.df["cutflow"][self.cat]

        if self.input_ntuples:
            self.schema = "parquet"
            print("Reading input ntuples: ", self.input_ntuples)
            self.events = ak.from_parquet(self.input_ntuples)
            nevents = self.nevents
            print(f"Number of events: {nevents}")
        else:
            self.schema = "coffea"

        # Build the array dictionary and the zipped dictionary for each sample
        self.array_dict = {s : defaultdict(dict) for s in self.samples}
        self.zipped_dict = {s : defaultdict(dict) for s in self.samples}
        self.features_dict = {s : defaultdict(dict) for s in self.samples}
        self.features_pad_dict = {s : defaultdict(dict) for s in self.samples}

    def load_config(self):
        '''Load the features, features_pad, awkward_collections and matched_collections_dict
        dictionaries with OmegaConf.'''
        print(f"Reading configuration file: {self.cfg}")
        self.cfg = OmegaConf.load(self.cfg)
        self.features = self.cfg["features"]
        self.features_pad = self.cfg["features_pad"]
        self.awkward_collections = self.cfg["awkward_collections"]
        self.matched_collections_dict = self.cfg["matched_collections"]

    def check_output(self):
        # Check if the output file is already existing
        if os.path.exists(self.output_file):
            raise ValueError(f"Output file {self.output_file} already exists.")
        # Check that the format of the output file is `.parquet`
        if not self.output_file.endswith(".parquet"):
            raise ValueError(f"Output file {self.output_file} should have the `.parquet` extension.")
        # Create the output directory if it does not exist
        os.makedirs(os.path.dirname(os.path.abspath(self.output_file)), exist_ok=True)

    def get_weight(self, sample, dataset):
        if self.schema == "coffea":
            weight = self.columns[sample][dataset][self.cat]["weight"].value
        elif self.schema == "parquet":
            weight = self.events.weight
        return weight

    @staticmethod
    def get_collection_dim(collection_dict):
        '''Get the dimension of the collection from the dictionary of arrays.'''
        ndims = []
        for key, values in collection_dict.items():
            ndims.append(values.ndim)
        ndims = set(ndims)
        if len(ndims) > 1:
            raise ValueError(f"The collection dimension cannot be determined as there are columns of different dimensions: {ndims}")
        return list(ndims)[0]

    def normalize_genweights(self):
        '''Since the array `weight` is filled on the fly with the weight associated with the event, 
        it does not take into account the overall scaling by the sum of genweights (`sum_genweights`).
        In order to correct for this, we have to scale by hand the `weight` array dividing by the sum of genweights.'''

        print("Normalizing genweights...")
        for sample in self.samples:
            if self.is_mc(sample):
                for dataset in self.get_datasets(sample):
                    weight = self.get_weight(sample, dataset)
                    if self.schema == "coffea":
                        weight_new = column_accumulator(weight / self.sum_genweights[dataset])
                        self.columns[sample][dataset][self.cat]["weight"] = weight_new
                    elif self.schema == "parquet":
                        self.events = ak.with_field(self.events, weight / self.sum_genweights[dataset], "weight")

    def load_features(self):
        '''Load the features dictionary with common features and sample-specific features.'''
        
        print("Loading features...")
        # Compose the features dictionary with common features and sample-specific features
        for sample in self.samples:
            self.features_dict[sample] = self.features["common"].copy()
            if sample in self.features["by_sample"].keys():
                self.features_dict[sample].update(self.features["by_sample"][sample])

            # Compose the dictionary of features to pad
            self.features_pad_dict[sample] = self.features_pad["common"].copy()
            if sample in self.features_pad["by_sample"].keys():
                self.features_pad_dict[sample].update(self.features_pad["by_sample"][sample])

            # Create a default dictionary of dictionaries to store the arrays
            self.array_dict[sample] = {k : defaultdict(dict) for k in self.features_dict[sample].keys()}

    def build_arrays_from_ntuples(self):
        '''Build the Momentum4D arrays for the jets, partons, leptons, met and higgs, reading from unflattened ntuples.'''

        print("Samples: ", self.samples)

        for sample in self.samples:

            # Build the array_dict from the unflattened ntuples read from parquet
            for collection, variable in self.features_dict[sample].items():
                for key_feature, key_coffea in variable.items():
                    if (collection == "JetGoodMatched") & (key_coffea == "provenance"):
                        try:
                            self.array_dict[sample][collection][key_feature] = self.events[f"JetGoodMatched_provenance"]
                        except:
                            self.array_dict[sample][collection][key_feature] = self.events[f"PartonMatched_provenance"]
                    else:
                        values = self.events[f"{collection}_{key_coffea}"]
                        if (collection not in self.awkward_collections) & (values.ndim > 1):
                            values = ak.flatten(values)
                        self.array_dict[sample][collection][key_feature] = values

                # Add padded features to the array, according to the features dictionary
                if collection in self.features_pad_dict[sample].keys():
                    for key_feature, value in self.features_pad_dict[sample][collection].items():
                        self.array_dict[sample][collection][key_feature] = value * ak.ones_like(self.array_dict[sample][collection]["pt"])

            # The awkward arrays are zipped together to form the Momentum4D arrays.
            # If the collection is not a Momentum4D array, it is zipped as it is.
            # Here we assume that the ntuples are unflattened, therefore we don't need to unflatten them before zipping.
            for collection in self.array_dict[sample].keys():
                ndim = self.get_collection_dim(self.array_dict[sample][collection])
                if (collection not in self.awkward_collections) & (ndim > 1):
                    # We flatten the collections that should not be saved as awkward collections after zipping
                    self.zipped_dict[sample][collection] = ak.flatten(ak.zip(self.array_dict[sample][collection], with_name='Momentum4D'))
                else:
                    self.zipped_dict[sample][collection] = ak.zip(self.array_dict[sample][collection], with_name='Momentum4D')
                print(f"Collection: {collection}")
                print("Fields: ", self.zipped_dict[sample][collection].fields)

            for collection in self.zipped_dict[sample].keys():
                # Pad the matched collections with None if there is no matching
                if collection in self.matched_collections_dict.keys():
                    matched_collection = self.matched_collections_dict[collection]
                    masked_arrays = ak.mask(self.zipped_dict[sample][matched_collection], self.zipped_dict[sample][matched_collection].pt==-999, None)
                    self.zipped_dict[sample][matched_collection] = masked_arrays
                    # Add the matched flag and the provenance to the matched jets
                    if collection == "JetGood":
                        is_matched = ~ak.is_none(masked_arrays, axis=1)
                        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], is_matched, "matched")
                        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], ak.fill_none(masked_arrays.prov, -1), "prov")

            # Add the remaining keys to the zipped dictionary
            if self.is_mc(sample):
                self.zipped_dict[sample]["event"] = ak.zip({"weight" : self.events.weight})

    def build_arrays_from_accumulator(self):
        '''Build the Momentum4D arrays for the jets, partons, leptons, met and higgs.'''
        
        print("Samples: ", self.samples)

        for sample in self.samples:
            datasets = self.get_datasets(sample)
            print("Datasets: ", datasets)

            # Accumulate ntuples from different data-taking eras
            # FIXME I think I just need to change what is done here since columns are already in the accumulator
            cs = accumulate([self.columns[sample][dataset][self.cat] for dataset in datasets])
            #cs = [self.columns[sample][dataset][self.cat] for dataset in datasets]

           # print(cs[0])

            # In order to get the numpy array from the column_accumulator, we have to access the `value` attribute.
            for collection, variables in self.features_dict[sample].items():
                for key_feature, key_coffea in variables.items():
                    print(collection, key_coffea)
                    if (collection == "JetGoodMatched") & (key_coffea == "provenance"):
                        self.array_dict[sample][collection][key_feature] = cs[f"PartonMatched_{key_coffea}"].value
                    elif (collection == "FatJetMatched") & (key_coffea == "provenance"):
                        self.array_dict[sample][collection][key_feature] = cs[f"HiggsMatched_{key_coffea}"].value
                    else:
                        self.array_dict[sample][collection][key_feature] = cs[f"{collection}_{key_coffea}"].value

                # Add padded features to the array, according to the features dictionary
                if collection in self.features_pad_dict[sample].keys():
                    for key_feature, value in self.features_pad_dict[sample][collection].items():
                        self.array_dict[sample][collection][key_feature] = value * np.ones_like(cs[f"{collection}_pt"].value)

            # The awkward arrays are zipped together to form the Momentum4D arrays.
            # If the collection is not a Momentum4D array, it is zipped as it is,
            # otherwise the Momentum4D arrays are zipped together and unflattened depending on the number of objects in the collection.
            for collection in self.array_dict[sample].keys():
                if collection in self.awkward_collections:
                    self.zipped_dict[sample][collection] = ak.unflatten(ak.zip(self.array_dict[sample][collection], with_name='Momentum4D'), cs[f"{collection}_N"].value)
                    print("N:", cs[f"{collection}_N"].value)
                #elif collection == "FatJetSorted":
                #    self.zipped_dict[sample][collection] = ak.unflatten(ak.zip(self.array_dict[sample][collection], with_name='Momentum4D'), np.ones_like(cs["FatJetSorted_phi"].value, dtype=int))
                #    if sample == "ttHTobb":
                #        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], ak.fill_none(masked_arrays.prov, -1), "prov")
                #    print("N?", np.ones_like(cs["FatJetSorted_phi"].value))
                else:
                    self.zipped_dict[sample][collection] = ak.zip(self.array_dict[sample][collection], with_name='Momentum4D')
                print(f"Collection: {collection}")
                print("Fields: ", self.zipped_dict[sample][collection].fields)

            for collection in self.zipped_dict[sample].keys():
                # Pad the matched collections with None if there is no matching
                if collection in self.matched_collections_dict.keys():
                    matched_collection = self.matched_collections_dict[collection]
                    masked_arrays = ak.mask(self.zipped_dict[sample][matched_collection], self.zipped_dict[sample][matched_collection].pt==-999, None)
                    self.zipped_dict[sample][matched_collection] = masked_arrays
                    # Add the matched flag and the provenance to the matched jets
                    if collection == "JetGood":
                        is_matched = ~ak.is_none(masked_arrays, axis=1)
                        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], is_matched, "matched")
                        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], ak.fill_none(masked_arrays.prov, -1), "prov")
                    if collection == "FatJetSorted":
                        is_matched = ~ak.is_none(masked_arrays, axis=1)
                        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], is_matched, "matched")
                        self.zipped_dict[sample][collection] = ak.with_field(self.zipped_dict[sample][collection], ak.fill_none(masked_arrays.prov, -1), "prov")
                        

            # Add the remaining keys to the zipped dictionary
            if self.is_mc(sample):
                self.zipped_dict[sample]["event"] = ak.zip({"weight" : cs["weight"].value})

    def save_parquet(self):
        '''Create the parquet files with the zipped dictionary.'''

        for sample in self.samples:
            # The Momentum4D arrays are zipped together to form the final dictionary of arrays.
            print("Zipping the collections into a single dictionary...")
            df_out = ak.zip(self.zipped_dict[sample], depth_limit=1)

            if len(self.samples) > 1:
                output_file = self.output_file.replace(".parquet", f"_{sample}.parquet")
            else:
                output_file = self.output_file
            print(f"Saving the output dataset to file: {output_file}")
            ak.to_parquet(df_out, output_file)
