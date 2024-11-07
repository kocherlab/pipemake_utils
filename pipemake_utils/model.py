import os
import json
import yaml
import logging
import copy

import numpy as np

from collections import defaultdict


class ModelFile(dict):
    def __init__(self, *arg, **kw):
        super(ModelFile, self).__init__(*arg, **kw)
        self.inds = []
        self.ind_file = ""
        self.exclude_file = ""

        if arg and self.confirm_model_instance(arg[1]):
            self.update_inds(arg[1])

    def __setitem__(self, *arg, **kw):
        super(ModelFile, self).__setitem__(*arg, **kw)

        if arg and self.confirm_model_instance(arg[1]):
            self.update_inds(model=arg[1])

    def __delitem__(self, key):
        super(ModelFile, self).__delitem__(key)
        self.update_inds()

    @classmethod
    def read(cls, filename):
        # Create ModelFile object
        models_to_return = cls()

        # Read the model file
        models_dict = _readFile(filename, _sniffFormat(filename))

        # Loop the parsed models
        for model_dict in models_dict:
            # Check if model already exists
            if str(model_dict["name"]) in models_to_return:
                raise Exception(f'Model {model_dict["name"]} already exists')

            # Create the model
            model = Model.from_dict(model_dict)

            # Save the model
            models_to_return[str(model.name)] = model

        logging.info(f"Finished reading model file ({filename})")

        # Return the models
        return models_to_return

    def write(self, filename, format="YAML", overwrite=True):
        # Check if the file is to be overwritten
        if not overwrite:
            # Check if the file exists
            if os.path.exists(filename):
                raise Exception(f"{filename} already exists")

        # Open the output file
        output_file = open(filename, "w")

        # Create the model data to be written
        models_data = [_d.to_dict() for _d in super(ModelFile, self).values()]

        # Write the json-formmated data to the output file
        if format.upper() == "JSON":
            output_file.write(json.dumps(models_data, indent=4))
        elif format.upper() == "YAML":
            output_file.write(yaml.dump(models_data, default_flow_style=False))
        else:
            raise Exception(
                f"Unknown format for {format}. Supported model file formats are JSON and YAML."
            )

        # Close the output file
        output_file.close()

        logging.info(f"Finished writing model file: {filename}")

    def copy_model(self, src_model_name, new_model_name):
        src_model = super(ModelFile, self).__getitem__(src_model_name)

        src_model_copy = copy.deepcopy(src_model)

        src_model_copy.name = new_model_name

        super(ModelFile, self).__setitem__(new_model_name, src_model_copy)

    def rename_model(self, src_model_name, new_model_name):
        src_model = super(ModelFile, self).pop(src_model_name)

        src_model.name = new_model_name

        super(ModelFile, self).__setitem__(new_model_name, src_model)

    def update_inds(self, model=None):
        if self.confirm_model_instance(model):
            # Return error if inds is empty
            if not model.inds:
                raise IOError("No individuals found in %s." % model.name)

            # Create a list of the unique individuals
            unique_inds = list(set(self.inds + model.inds))

        else:
            # Create an empty list for the unique individuals
            unique_inds = []

            # Loop the models in the file
            for model_in_file in super(ModelFile, self).values():
                # Create a list of the unique individuals
                unique_inds = list(set(unique_inds + model_in_file.inds))

        # Store the individuals
        self.inds = unique_inds

    def create_ind_file(self, file_ext="", file_path="", overwrite=False):
        # Assign the filename for the population file
        ind_filename = "unique_individuals" + file_ext

        # If a path is assigned, create the file at the specified location
        if file_path:
            ind_filename = os.path.join(file_path, ind_filename)

        # Check if previous files should be overwriten
        if not overwrite:
            # Check if the file already exists
            if os.path.isfile(ind_filename):
                raise IOError("Individuals file exists. Use --overwrite to ignore")

        # Create the population file
        ind_file = open(ind_filename, "w")
        ind_file.write("%s\n" % "\n".join(self.inds))
        ind_file.close()

        # Save the individuals filename
        self.ind_file = ind_filename

    def delete_ind_file(self):
        # Check if an individuals file was created
        if self.ind_file:
            # Delete the individuals file
            os.remove(self.ind_file)

            # Remove the filename
            self.ind_file = ""

    def create_exclude_ind_file(
        self, inds_to_include=[], file_ext="", file_path="", overwrite=False
    ):
        # Assign the filename for the population file
        ind_filename = "exclude_individuals" + file_ext

        # If a path is assigned, create the file at the specified location
        if file_path:
            ind_filename = os.path.join(file_path, ind_filename)

        # Check if previous files should be overwriten
        if not overwrite:
            # Check if the file already exists
            if os.path.isfile(ind_filename):
                raise IOError("Individuals file exists.")

        # Create exclude list by removing included individuals
        exclude_inds = list(set(self.inds) - set(inds_to_include))

        # Create the population file
        ind_file = open(ind_filename, "w")
        ind_file.write("%s\n" % "\n".join(exclude_inds))
        ind_file.close()

        # Save the individuals filename
        self.exclude_file = ind_filename

    def delete_exclude_ind_file(self):
        # Check if an individuals file was created
        if self.exclude_file:
            # Delete the individuals file
            os.remove(self.exclude_file)

            # Remove the filename
            self.exclude_file = ""

    @staticmethod
    def confirm_model_instance(unknown):
        if isinstance(unknown, Model):
            return True

        else:
            return False


class Model:
    def __init__(self, name, pops={}, **kwargs):
        self.name = name
        self.ind_dict = defaultdict(list)
        self._pop_kwargs = {}
        self._kwargs = kwargs
        self.pop_files = []
        self.ind_file = ""

        # Loop the populations in the pops dict
        for pop, pop_kwargs in pops.items():
            # Check if the pop already exists
            if pop in self.ind_dict:
                raise Exception(f"{pop} already exists")

            # Check if the inds are provided
            if "inds" not in pop_kwargs:
                raise Exception(f"No individuals found for {pop}")

            # Assign the inds
            inds = pop_kwargs["inds"]

            # Create a set of individuals that are already assigned to a population
            ind_set = set(inds).intersection(self.inds)

            if ind_set:
                raise Exception(
                    f'Individuals already assigned to another population: {", ".join(ind_set)}'
                )

            # Assign the individuals to the population
            self.ind_dict[pop] = [str(ind) for ind in inds]

            # Save the population kwargs
            del pop_kwargs["inds"]
            if pop_kwargs:
                self._pop_kwargs[pop] = pop_kwargs

            logging.info(f"Population {pop} added to {self.name}")

    def __str__(self):
        return self.name

    @property
    def inds(self):
        return [ind for pop in self.ind_dict for ind in self.ind_dict[pop]]

    @property
    def npop(self):
        return len(self.ind_dict)

    @property
    def nind(self):
        return {pop: len(inds) for pop, inds in self.ind_dict.items()}

    @classmethod
    def from_dict(cls, model_dict):
        return cls(**model_dict)

    def assign_pop(self, pop, inds=[]):
        if not inds:
            raise Exception("No individuals found for %s" % pop)

        if pop in self.ind_dict:
            raise Exception(f"{pop} already exists")

        # Assign the dict
        self.ind_dict[pop] = [str(ind) for ind in inds]

    def remove_pop(self, pop):
        # Confirm the pop is in the model
        if str(pop) not in self.ind_dict:
            # Raise error if pop not found
            raise Exception(f"{pop} not found")

        del self.ind_dict[pop]

    def update_pop(self, pop, inds=[], rm_inds=[]):
        # Confirm the pop is in the model
        if str(pop) not in self.ind_dict:
            # Raise error if pop not found
            raise Exception(f"{pop} not found")

        if inds:
            self.ind_dict[pop].extend([str(ind) for ind in inds])

        if rm_inds:
            for rm_ind in rm_inds:
                if str(rm_ind) in self.ind_dict[pop]:
                    self.ind_dict[pop].remove(str(rm_ind))

    def sample_pop(self, pop, sample_size, with_replacements=False):
        # Confirm the pop is in the model
        if str(pop) not in self.ind_dict:
            # Raise error if pop not found
            raise Exception("%s not found" % pop)

        # Confirm the sample size is an int
        try:
            sample_size = int(sample_size)

        except Exception as e:
            raise e

        # Check if the sample size is larger than the pop
        if int(sample_size) > self.nind[pop]:
            # Raise error if sample_size is larger
            raise Exception(f"The samples size ({sample_size}) is larger than {pop}")

        # Use numpy choice to randomly sample the pop
        sampled_inds = np.random.choice(
            self.ind_dict[pop], sample_size, replace=with_replacements
        )

        # Save the sampled inds as a list
        self.ind_dict[pop] = list(sampled_inds)

    def sample_pops(self, sample_size, with_replacements=False):
        # Confirm the sample size is an int
        try:
            sample_size = int(sample_size)

        except Exception as e:
            raise e

        # Loop each pop in the pop list
        for pop in self.ind_dict:
            # Check if the sample size is larger than the pop
            if int(sample_size) > self.nind[pop]:
                # Raise error if sample_size is larger
                raise Exception(
                    f"The samples size ({sample_size}) is larger than {pop}"
                )

        # Loop each pop in the pop list, if no error raised
        for pop in self.ind_dict:
            # Use numpy choice to randomly sample the pop
            sampled_inds = np.random.choice(
                self.ind_dict[pop], sample_size, replace=with_replacements
            )

            # Save the sampled inds as a list
            self.ind_dict[pop] = list(sampled_inds)

    def create_pop_files(self, file_ext="", file_path="", overwrite=False):
        for pop in self.ind_dict:
            # Assign the filename for the population file
            pop_filename = f"{pop}.{file_ext}"

            # If a path is assigned, create the file at the specified location
            if file_path:
                pop_filename = os.path.join(file_path, pop_filename)

            # Check if previous files should be overwriten
            if not overwrite:
                # Check if the file already exists
                if os.path.isfile(pop_filename):
                    raise IOError(f"Population file exists: {pop_filename}")

            # Create the population file
            pop_file = open(pop_filename, "w")
            pop_file.write("%s\n" % "\n".join(self.ind_dict[pop]))
            pop_file.close()

            # Save the population filename
            self.pop_files.append(pop_filename)

        logging.info(f"Population files created for {self.name}")

    def delete_pop_files(self):
        # Check if pop files were created
        if len(self.pop_files) != 0:
            # Loop the created pop files
            for pop_file in self.pop_files:
                # Delete the pop file
                os.remove(pop_file)

            # Remove the filenames
            self.pop_files = []

            logging.info(f"Population files deleted for {self.name}")

    def create_ind_file(self, file_ext="", file_path="", overwrite=False):
        # Assign the filename for the population file
        ind_filename = "individual.keep" + file_ext

        # If a path is assigned, create the file at the specified location
        if file_path:
            ind_filename = os.path.join(file_path, ind_filename)

        # Check if previous files should be overwriten
        if not overwrite:
            # Check if the file already exists
            if os.path.isfile(ind_filename):
                raise IOError("Individuals file exists.")

        # Create the population file
        ind_file = open(ind_filename, "w")
        ind_file.write("%s\n" % "\n".join(self.inds))
        ind_file.close()

        # Save the individuals filename
        self.ind_file = ind_filename

        logging.info(f"Individuals file created for {self.name}")

    def delete_ind_file(self):
        # Check if an individuals file was created
        if self.ind_file:
            # Delete the individuals file
            os.remove(self.ind_file)

            # Remove the filename
            self.ind_file = ""

            logging.info(f"Individuals file deleted for {self.name}")

    def return_pop(self, ind):
        # Loop the pops within the model
        for pop, inds in self.ind_dict.items():
            # Check if the ind belongs to the pop
            if ind in inds:
                return pop

        # Raise error if ind not found
        raise Exception(f"Individual ({ind}) not found within {self.name}")

    def to_dict(self):
        # Create a dictionary to store the model data
        model_dict = {}

        # Save the model name
        model_dict["name"] = self.name

        # Save the model tree if it exists
        for arg, value in self._kwargs.items():
            model_dict[arg] = value

        # Create a dictionary to store the populations
        pop_dict = {}

        # Loop the populations in the model
        for pop in self.ind_dict:
            # Store the individuals in the population
            pop_dict[pop] = {}
            pop_dict[pop]["inds"] = self.ind_dict[pop]

        # Save the populations in the model
        model_dict["pops"] = pop_dict

        return model_dict


def _sniffFormat(filename):
    # Try to open the model file in JSON format
    try:
        with open(filename, "r") as model_file:
            json.load(model_file)
        return "JSON"
    except FileNotFoundError as e:
        raise e
    except Exception:
        pass

    # Try to open the model file in YAML format
    try:
        with open(filename, "r") as model_file:
            yaml.load(model_file, Loader=yaml.FullLoader)
        return "YAML"
    except FileNotFoundError as e:
        raise e
    except Exception:
        pass

    raise Exception(
        f"Unknown format for {filename}. Supported model file formats are JSON and YAML."
    )


def _readFile(filename, format):
    class UniqueKeyLoader(yaml.SafeLoader):
        def construct_mapping(self, node, deep=False):
            mapping = set()
            for key_node, value_node in node.value:
                if ":merge" in key_node.tag:
                    continue
                key = self.construct_object(key_node, deep=deep)
                if key in mapping:
                    raise ValueError(f"Duplicate {key!r} key found in YAML.")
                mapping.add(key)
            return super().construct_mapping(node, deep)

    with open(filename, "r") as model_file:
        if format == "JSON":
            model_dict = json.load(model_file)
        elif format == "YAML":
            model_dict = yaml.load(model_file, Loader=UniqueKeyLoader)

    return model_dict
