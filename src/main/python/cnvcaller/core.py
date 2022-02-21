#!/usr/bin/env python3

"""
Created:      28/09/2021
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2021 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import enum
import io
import os
import pickle
import re
import sys
import argparse
import yaml

import numpy as np
import pandas as pd
import pyranges
import sklearn.decomposition
import IlluminaBeadArrayFiles

# Metadata
__program__ = "CNV-caller"
__author__ = "C.A. (Robert) Warmerdam"
__email__ = "c.a.warmerdam@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


# Constants
DEFAULT_FINAL_REPORT_COLS = ["Sample ID", "SNP Name", "GType", "SNP",
                             "X", "Y", "B Allele Freq", "Log R Ratio"]
AUTOSOMES_CHR = ["{}".format(chrom) for chrom in range(1, 23)]


# Classes
class ArgumentParser:
    def __init__(self):
        self.sub_commands = list()
        self.parser = self.create_argument_parser()
        self.add_command_choice_argument(self.parser)
        self.add_bead_pool_manifest_argument(self.parser)
        self.add_corrective_variants_argument(self.parser)
        self.add_sample_sheet_argument()
        self.add_bed_path_parameter()
        self.add_debug_parameter()
        self.add_config_parameter()
        self.add_out_argument(self.parser)
    class SubCommand(enum.Enum):
        VARIANTS = "variants"
        DATA = "data"
        FIT = "fit"
        CALL = "call"
        @classmethod
        def list(cls):
            return list(map(lambda c: c.value, cls))
        def __str__(self):
            return self.name.lower()
    def add_subparsers(self):
        subparsers = self.parser.add_subparsers(help='procedure to run')
        parser_for_input_preparation = subparsers.add_parser(
            'stage-data', help="Process final report files to pickled panda DataFrames.")
        self.add_final_report_path_argument(parser_for_input_preparation)
        self.add_corrective_variants_argument(parser_for_input_preparation)
        self.add_staged_data_output_argument(parser_for_input_preparation)
        parser_for_correction = subparsers.add_parser(
            'correction', help='Perform decomposition for adjustment of raw intensities.')
        self.add_staged_data_argument(parser_for_correction)
        self.add_out_argument(parser_for_correction)
        parser_for_fit = subparsers.add_parser('fit', help='Perform decomposition in locus of interest"')
        self.add_staged_data_argument(parser_for_fit)
        parser_for_calling = subparsers.add_parser('call', help="Call CNVs using correction and calling parameters")
        self.add_staged_data_argument(parser_for_calling)
        self.add_calling_parameter_argument(parser_for_calling)
    def parse_input(self, argv):
        """
        Parse command line input.
        :param argv: given arguments
        :return: parsed arguments
        """
        command_parser = self.create_command_parser()
        args_partial = command_parser.parse_known_args(argv)
        self.sub_commands = [self.SubCommand[command.upper()] for command in args_partial[0].command]
        if "data" in self.sub_commands and len(self.sub_commands) > 1:
            raise argparse.ArgumentError(command_parser._actions[1],
                                         "The command data cannot be used in combination with other commands")
        self.extend_argument_parser()
        args_remainder = self.parser.parse_args(argv)
        return args_remainder
    def is_action_requested(self, sub_command):
        return sub_command in self.sub_commands
    @staticmethod
    def create_command_parser():
        parser = argparse.ArgumentParser(
            description="CNV-calling algorithm",
            usage=(os.linesep.join([
                "core.py <command> [<args>]",
                "",
                "The available commands are",
                "  variants       Samples variants in proportion to chromosome",
                "                 length.",
                "  data           Prepares arrays of raw intensity data from",
                "                 final report files.",
                "  fit            Perform decomposition in locus of interest",
                "  call           Call CNVs using correction and calling parameters"])))
        ArgumentParser.add_command_choice_argument(parser)
        return parser
    @classmethod
    def add_command_choice_argument(cls, parser):
        parser.add_argument('command',
                            help='Command(s) to run', nargs="+",
                            choices=cls.SubCommand.list())
    @staticmethod
    def create_argument_parser():
        """
        Method creating an argument parser
        :param command: 
        :return: parser
        """
        parser = argparse.ArgumentParser(description="CNV-calling algorithm",
                                         formatter_class=argparse.RawDescriptionHelpFormatter)
        return parser
    def add_sample_sheet_argument(self):
        self.parser.add_argument('-s', '--sample-sheet', type=self.is_readable_file,
                                 required=True,
                                 default=None,
                                 help="Samplesheet")
    def add_final_report_path_argument(self, parser):
        parser.add_argument('-g', '--final-report-file-path', type=self.is_readable_file, required=True, default=None,
                            help="Path to where final report files are located")
    def add_out_argument(self, parser):
        parser.add_argument('-o', '--out', type=self.can_write_to_file_path,
                            required=True, default=None,
                            help="File prefix the output can be written to. ")
    def add_bed_path_parameter(self):
        self.parser.add_argument('-b', '--bed-file', type=self.is_readable_file,
                                 required=True,
                                 default=None,
                                 help="Bed file detailing a locus of interest."
                                      "This is excluded in corrections, and exclusively"
                                      "assessed in the fitting and calling steps")
    @classmethod
    def can_write_to_file_path(cls, file):
        """
        Checks whether the given directory is readable
        :param file: a path to a directory in string format
        :return: file
        :raises: Exception: if the dirname of the given path is invalid
        :raises: Exception: if the dirname of the given directory is not writable
        """
        directory = os.path.dirname(file)
        if not os.path.isdir(directory):
            raise argparse.ArgumentTypeError("directory: {0} is not a valid path".format(directory))
        if os.access(directory, os.R_OK):
            return file
        else:
            raise argparse.ArgumentTypeError("directory: {0} is not a readable dir".format(directory))
    @classmethod
    def is_readable_dir(cls, directory):
        """
        Checks whether the given directory is readable
        :param directory: a path to a directory in string format
        :return: train_directory
        :raises: Exception: if the given path is invalid
        :raises: Exception: if the given directory is not accessible
        """
        if not os.path.isdir(directory):
            raise argparse.ArgumentTypeError("directory: {0} is not a valid path".format(directory))
        if os.access(directory, os.R_OK):
            return directory
        else:
            raise argparse.ArgumentTypeError("directory: {0} is not a readable dir".format(directory))
    @classmethod
    def is_readable_file(cls, file_path):
        """
        Checks whether the given directory is readable
        :param file_path: a path to a file in string format
        :return: file_path
        :raises: Exception: if the given path is invalid
        :raises: Exception: if the given directory is not accessible
        """
        if not os.path.isfile(file_path):
            raise argparse.ArgumentTypeError("file path:{0} is not a valid file path".format(file_path))
        if os.access(file_path, os.R_OK):
            return file_path
        else:
            raise argparse.ArgumentTypeError("file path:{0} is not a readable file".format(file_path))
    @classmethod
    def is_data(cls, path, extension_expected):
        """
        Parses a string that specifies the path(s) of the data.
        :param path:
        :param extension_expected:
        :return:
        """
        extension_actual = os.path.splitext(path)[1]
        if extension_actual == extension_expected:
            file = argparse.FileType('r')(path)
            return path
        else:
            raise argparse.ArgumentTypeError(
                "file path:{0} is not of type .'{1}'. (.'{2}' received)".format(
                    path, extension_actual, extension_expected))
    @classmethod
    def is_writable_location(cls, path):
        if os.access(path, os.W_OK):
            return path
        else:
            raise argparse.ArgumentTypeError("directory: {0} is not a writable path".format(path))
    def add_corrective_variants_argument(self, parser):
        parser.add_argument(
            '-v', '--corrective-variants', type=self.is_readable_file,
            help="filters out all variants that are not listed here"
        )
    def extend_argument_parser(self):
        sub_command_mapping = {
            self.SubCommand.VARIANTS:
                {},
            self.SubCommand.DATA:
                {self.add_final_report_path_argument},
            self.SubCommand.FIT:
                {self.add_staged_data_argument},
            self.SubCommand.CALL:
                {self.add_staged_data_argument,
                 self.add_batch_weights_argument,
                 self.add_calling_cluster_weight_argument}}
        methods_to_run = set.union(*[sub_command_mapping.get(sub_command) for sub_command in self.sub_commands])
        if self.add_staged_data_argument in methods_to_run and self.add_staged_data_output_argument in methods_to_run:
            methods_to_run.discard(self.add_staged_data_argument)
        for method in methods_to_run:
            method(self.parser)
    def add_batch_weights_argument(self, parser):
        parser.add_argument('-C', '--correction', type=self.is_readable_dir,
                            required=True, nargs='+', default=None,
                            help="path where batch correction weights, and corrected data are stored."
                                 "output of the batch weighting step")
    def add_calling_cluster_weight_argument(self, parser):
        parser.add_argument('-F', '--cluster-file', type=self.is_readable_dir,
                            required=True, nargs='+', default=None,
                            help="path where cluster weights are stored."
                                 "output of the fitting step")
    def add_bead_pool_manifest_argument(self, parser):
        parser.add_argument('-bpm', '--bead-pool-manifest', type=self.is_readable_file,
                                 required=True, default=None,
                                 help="path to a .bpm file corresponding to the genotyping array")
    def add_staged_data_output_argument(self, parser):
        parser.add_argument(
            '--out', type=self.can_write_to_file_path,
            metavar="PATH_TO_NEW_PICKLE_FILE",
            help="path to a pickle file")
    def add_staged_data_argument(self, parser):
        parser.add_argument(
            '--input', '-i', type=lambda arg: self.is_data(arg, '.pkl'), nargs="+",
            metavar="PATHS_TO_PICKLE_INTENSITY FILES",
            help=os.linesep.join([
                "paths directing to .pkl files with intensity data.",
                "Columns should correspond to samples and rows should correspond to variants."]))
    def add_debug_parameter(self):
        self.parser.add_argument(
            '--debug', '-d', action="store_true", default=False,
            help="write files useful for debugging"
        )
    def add_config_parameter(self):
        self.parser.add_argument(
            '--config', '-c', type=self.config, required=True,
            help="config file")
    @classmethod
    def config(cls, path):
        return yaml.safe_load(open(path))


class IntensityDataReader:
    def __init__(self, sample_list):
        self._sample_list = sample_list
    def load(self, data):
        data_frame_list = list()
        for file in data:
            data_frame = pd.read_pickle(file, compression=None)
            # if data_frame.index != self._variant_list:
            #     raise Exception("variants don't match")
            data_frame_list.append(data_frame)
        intensity_data = pd.concat(data_frame_list)
        sample_list = intensity_data["Sample ID"].unique()
        missing_samples = np.setdiff1d(self._sample_list, sample_list)
        if len(missing_samples) > 0:
            print(missing_samples)
            print("warning: {}".format(missing_samples))
        excess_samples = np.setdiff1d(sample_list, self._sample_list)
        if len(excess_samples) > 0:
            print(excess_samples)
            print("warning excess: {}".format(excess_samples))
        return intensity_data


class FinalReportReaderException(Exception):
    """
    Exception raised for errors in the input final reports
    """
    def __init__(self, message, file, line_index):
        self.file = file
        self.message = message
        self.line_index = line_index
        super().__init__(self.message)
    def __str__(self):
        return os.linesep.join(["Exception encountered in file:",
                                "'{1}', on line {2}: {3}"]).format(
            self.file, self.line_index, self.message)


class FinalReportGenotypeDataReader:
    """
    Read final report files
    """
    new_part_pattern = re.compile(r"^\[\w+]$")
    def __init__(self, path, sample_list, variant_list):
        self._part_key = None
        self.sep = "\t"
        self._path = path
        self._sample_list = sample_list
        self._variants_to_include = variant_list
        self._variants_to_include_indices = None
        self._line_counter = 0
    def read_intensity_data(self):
        """
        Method that reads the intensity values from a final report file.
        :return: Data frame.
        """
        with open(self._path, self.get_reading_mode()) as buffer:
            part_buffer = list()
            for line in buffer:
                self._line_counter += 1
                # The header is in key value format
                key_match = self.new_part_pattern.match(line)
                if key_match:
                    print(key_match.group(0))
                    if self._part_key is None:
                        pass
                    elif self._part_key == "[Header]":
                        self.parse_header(part_buffer)
                    else:
                        raise NotImplementedError(
                            "Parts other than the '[Header]' and '[Data]' part not supported"
                        )
                    del part_buffer[:]
                    self._part_key = key_match.group(0)
                    if self._part_key == "[Data]":
                        data_frame = self._read_data(buffer)
                else:
                    part_buffer.append(line)
        return data_frame
    def get_reading_mode(self):
        reading_mode = "r"
        if self._path.endswith(".gz"):
            reading_mode = "rb"
        return reading_mode
    def parse_header(self, part_buffer):
        pass
    def _read_data(self, buffer):
        data_array_list = list()
        sample_list = list()
        columns = pd.read_csv(io.StringIO(buffer.readline()), nrows=0, sep=self.sep).columns.to_list()
        sample_id_index = columns.index("Sample ID")
        sample_buffer = io.StringIO()
        current_sample = None
        sample_counter = 0
        for line in buffer:
            self._line_counter += 1
            # key_match = self.new_part_pattern.match(line)
            # if key_match:
            #     self._part_key = key_match.group(0)
            #     break
            splitted = line.split(self.sep, sample_id_index + 1)
            sample_id = splitted[sample_id_index]
            if sample_id == current_sample or current_sample is None:
                sample_buffer.write(line)
                if current_sample is None:
                    current_sample = sample_id
                    print(sample_counter, current_sample)
            else:
                if current_sample in self._sample_list:
                    data_array_list.append(self._read_sample_intensities(sample_buffer, columns))
                    sample_list.append(current_sample)
                    sample_counter += 1
                else:
                    print("Skipping sample since it is not in the sample list: {}".format(current_sample))
                # Reset buffer
                sample_buffer.truncate(0)
                sample_buffer.seek(0)
                sample_buffer.write(line)
                current_sample = sample_id
                print(sample_counter, current_sample)

        if current_sample is not None and current_sample in self._sample_list:
            data_array_list.append(self._read_sample_intensities(sample_buffer, columns))
            sample_list.append(current_sample)
            sample_counter += 1
        else:
            print("Skipping sample since it is not in the sample list: {}".format(current_sample))

        # if len(np.intersect1d(self._sample_list,
        #                       sample_list)) != len(self._sample_list):
        #     raise FinalReportReaderException(
        #         "Samples in manifest do not perfectly intersect with samples in sample sheet",
        #         self._path,
        #         self._line_counter)
        return pd.concat(data_array_list)
    def _read_sample_intensities(self, buffer, columns):
        buffer.seek(0)
        sample_data_frame = pd.read_csv(buffer, names=columns,
                                        sep=self.sep,
                                        usecols=DEFAULT_FINAL_REPORT_COLS,
                                        dtype={"Sample ID": str, "SNP Name": str},
                                        index_col="SNP Name")
        sample_data_frame = sample_data_frame.loc[
            self._variants_to_include.Name.to_numpy(),:]
        sample_data_frame["R"] = sample_data_frame[["X", "Y"]].sum(axis=1).values
        if not np.all(self._variants_to_include.Name.to_numpy() == sample_data_frame.index.to_numpy()):
            raise FinalReportReaderException(
                ("Variants in variants file do not perfectly intersect with variants of sample {}"
                    .format(sample_data_frame["Sample ID"][0])), self._path, self._line_counter)
        return sample_data_frame


class IntensityCorrection:
    def __init__(self, variant_list_for_locus, pca_n_components=None,
                 pca_scaling=True, regression_fit_intercept=False):
        self._variant_list_for_locus_of_interest = variant_list_for_locus
        self._scale = pca_scaling
        self._pca = sklearn.decomposition.PCA(
            n_components=pca_n_components)
        self._correction_model = sklearn.linear_model.LinearRegression(
            fit_intercept=regression_fit_intercept)
        self._standardize_scaler = sklearn.preprocessing.StandardScaler()
        self._corrected = None
    def fit(self, intensity_data):
        if self._scale:
            print("Scaling")
            intensity_data_preprocessed = pd.DataFrame(
                self._standardize_scaler.fit_transform(intensity_data.loc[:, self.variant_indices(intensity_data)]),
                columns=intensity_data.columns[self.variant_indices(intensity_data)],
                index=intensity_data.index)
            print("Scaling performed")
        else:
            intensity_data_preprocessed = intensity_data.loc[:,
                                          self.indices_not_in_locus_of_interest(intensity_data)]
        # Calculate the eigenvectors used to correct the correction variants.
        # These eigenvectors represent how to scale each variant in a sample
        # so that the result explaines covariability among the variants.
        # I.e., the projected principal components (PCs) explain batch effects.
        # We want to regress out these batch effects in the locus of interest.
        self._principal_components = pd.DataFrame(
            self._pca.fit_transform(intensity_data_preprocessed),
            index=intensity_data.index)
        # The projected principal components explain batch effects.
        # We try to explain as much of the locus of interest using the PCs
        # The residuals can be used in further analyses.
        self._correction_model.fit(
            self._principal_components, intensity_data.loc[:, self._variant_list_for_locus_of_interest])
        # Write intensities of locus of interest corrected for batch effects.
        self._corrected = self._correct_batch_effects(intensity_data, self._principal_components)
        return
    def variant_indices(self, intensity_data):
        return np.logical_and(
            self.indices_not_in_locus_of_interest(intensity_data),
            ~intensity_data.isnull().any(axis=0))
    def correct_intensities(self, intensity_data):
        if self._scale:
            intensity_data_preprocessed = pd.DataFrame(
                self._standardize_scaler.transform(intensity_data.loc[:,
                                                   self.indices_not_in_locus_of_interest(intensity_data)]),
                columns=intensity_data.columns[self.indices_not_in_locus_of_interest(intensity_data)],
                index=intensity_data.index)
        else:
            intensity_data_preprocessed = intensity_data.loc[:,
                                          self.indices_not_in_locus_of_interest(intensity_data)]
        # Get batch effects by calculating principal components
        principal_components = pd.DataFrame(
            self._pca.transform(intensity_data_preprocessed),
            index=intensity_data.index)
        # The principal components depict batch effects.
        # Here, we predict the batch effects on the locus of interest.
        # Using the predicted batch effects, we can correct the locus of interest for the
        # expected batch effects.
        corrected_intensities = self._correct_batch_effects(intensity_data, principal_components)
        return corrected_intensities
    def indices_not_in_locus_of_interest(self, intensity_data):
        return ~intensity_data.columns.isin(
            self._variant_list_for_locus_of_interest)
    def _correct_batch_effects(self, intensity_data, principal_components):
        # The principal components depict batch effects.
        # Here, we predict the batch effects on the locus of interest.
        # Using the predicted batch effects, we can correct the locus of interest for the
        # expected batch effects.
        predicted_batch_effects_on_locus_of_interest = self._correction_model.predict(
            principal_components)
        # We can correct the locus of interest by subtracting the predicted batch effects
        # from the raw intensity data.
        corrected_intensities = intensity_data[self._variant_list_for_locus_of_interest] - \
                                predicted_batch_effects_on_locus_of_interest
        return corrected_intensities
    def write_output(self, path):
        pickle.dump(self, open(
            ".".join([path, "intensity_correction", "mod", "pkl"]), "wb"))
        self._principal_components.to_csv(
            ".".join([path, "intensity_correction", "pcs", "csv", "gz"]))
        self._pca.explained_variance_.to_csv(
            ".".join([path, "intensity_correction", "eigenvalues", "csv", "gz"]))
        self._corrected.to_csv(
            ".".join([path, "intensity_correction", "corrected", "csv", "gz"]))
    @classmethod
    def load_instance(cls, path):
        return pickle.load(open(
            ".".join([path, "intensity_correction", "mod", "pkl"]), "rb"))


# Functions
def calculate_downsampling_factor(grouped_data_frame, N):
    grouped_data_frame['proportionsObserved'] = (
        grouped_data_frame.shape[0] / N)
    grouped_data_frame['downsamplingFactor'] = (
        grouped_data_frame.proportionsExpected / grouped_data_frame.proportionsObserved)
    return grouped_data_frame


def draw_variants_proportionate(grouped_data_frame, max_downsampling_factor):
    print(grouped_data_frame)
    # Get the downsampling factor for the specific chromosome.
    downsampling_factor = grouped_data_frame.downsamplingFactor.unique() / max_downsampling_factor
    print(grouped_data_frame.downsamplingFactor.unique())
    print(max_downsampling_factor)
    print(downsampling_factor)
    # Perform sampling
    return (grouped_data_frame
            .sample(frac=float(downsampling_factor), replace=False))


def sample_corrective_variants_proportionally(corrective_variant_path, manifest_ranges):
    # Read the names of those variants that adhere to a number of criteria.
    # This list can be obtained by running PLINK for instance.
    corrective_variant_names = pd.read_csv(corrective_variant_path, header=None)[0].to_list()
    # Now get details for these variants from the manifest file.
    corrective_variants = manifest_ranges[manifest_ranges.Name.isin(corrective_variant_names)]
    print(len(corrective_variants))
    # We need to sample a number of variants so that a number of variants are
    # drawn that are proportionate to the total length of each autosome.
    # We therefore get the autosome sizes here.
    chromosome_sizes = pyranges.data.chromsizes().as_df()
    chromosome_sizes.Chromosome = chromosome_sizes.Chromosome.apply(lambda x: x.lstrip("chr"))
    filtered_chromosome_sizes = chromosome_sizes.loc[chromosome_sizes.Chromosome.isin(AUTOSOMES_CHR)]
    # Now, we calculate the proportion of the total autosomes each autosome represents.
    filtered_chromosome_sizes.loc[:,'proportionsExpected'] = \
        filtered_chromosome_sizes.End / np.sum(filtered_chromosome_sizes.End)
    # We rename the data frame for easy merging.
    filtered_chromosome_sizes = filtered_chromosome_sizes.rename(
        columns={"Start": "ChromSizeStart", "End": "ChromSizeEnd"})
    # We need to select variants so that these are equally distributed across chromosomes.
    # We do this by sampling n variants in each chromosome,
    # where n denotes the proportional length of every chromosome, multiplied by the number of variants to
    # sample
    # Calculate what proportion of variants in each chromosome should be
    # discarded.
    corrective_variants_extended = (corrective_variants.df
        .merge(filtered_chromosome_sizes, on='Chromosome')
        .groupby('Chromosome')
        .apply(calculate_downsampling_factor, len(corrective_variants)))
    # Calculate the maximum downsampling factor
    corrective_variants_grouped = corrective_variants_extended.groupby('Chromosome')
    max_downsampling_factor = np.max(corrective_variants_extended.downsamplingFactor)
    sampled_corrective_variants_list = list()
    for group in corrective_variants_grouped.groups.keys():
        grouped_data_frame = corrective_variants_grouped.get_group(group)
        downsampling_factor = grouped_data_frame.downsamplingFactor.unique() / max_downsampling_factor
        sampled_corrective_variants_list.append(grouped_data_frame
                                                .sample(frac=float(downsampling_factor), replace=False))
    return pd.concat(sampled_corrective_variants_list).loc[:,["Chromosome", "Start", "End", "Name"]]


# Main

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = ArgumentParser()
    args = parser.parse_input(argv[1:])

    # Read locus of interest
    locus_of_interest = pd.read_csv(
        args.bed_file, index_col=False,
        names=("Chromosome", "Start", "End", "Name"),
        dtype={"Chromosome": str},
        sep="\t")

    # Convert the locus of interest to a pyranges object
    locus_ranges = pyranges.PyRanges(locus_of_interest)

    # Read the bead pool manifest file
    manifest = IlluminaBeadArrayFiles.BeadPoolManifest(args.bead_pool_manifest)

    # Process manifest variant list
    variant_list = list()
    for variant_index in range(manifest.num_loci):
        variant_list.append((
            manifest.chroms[variant_index],
            manifest.map_infos[variant_index],
            manifest.map_infos[variant_index],
            manifest.names[variant_index],
            manifest.snps[variant_index]))

    # Get table for the manifest
    manifest_data_frame = pd.DataFrame(
        variant_list, columns=("Chromosome", "Start", "End", "Name", "Alleles"))

    manifest_data_frame[['Ref', 'Alt']] = manifest_data_frame['Alleles'].str.split('\[(\w)/(\w)\]', expand=True).iloc[:,[1,2]]
    manifest_ranges = pyranges.PyRanges(manifest_data_frame)

    # Get the intersect between the variants that are in the manifest and the locus of interest
    variants_in_locus = manifest_ranges.intersect(locus_ranges)

    # Read the sample sheet
    sample_sheet = pd.read_csv(args.sample_sheet, sep=",")

    if parser.is_action_requested(ArgumentParser.SubCommand.VARIANTS):

        # Sample corrective variants
        sampled_corrective_variants = sample_corrective_variants_proportionally(
            args.corrective_variants, manifest_ranges)

        corrective_variants_dataframe = (pd
            .merge(sampled_corrective_variants, variants_in_locus,
                   how="outer", indicator=True)
            .mask["merge"] == "left_only")

        # variants as pyranges object
        corrective_variants_dataframe.to_csv(
            "{}.corrective.bed".format(args.out), index=False, sep="\t", header=False)

        variants_in_locus.as_df().to_csv(
            "{}.locus.bed".format(args.out), index=False, sep="\t", header=False)

    else:

        # Read locus of interest
        sampled_corrective_variants = pd.read_csv(
            args.corrective_variants, index_col=False,
            names=("Chromosome", "Start", "End", "Name"),
            dtype={"Chromosome": str},
            sep="\t")

    # Convert the locus of interest to a pyranges object
    variants_to_read = pyranges.concat([
        pyranges.PyRanges(sampled_corrective_variants),
        variants_in_locus])

    if parser.is_action_requested(ArgumentParser.SubCommand.DATA):

        # Get intensity data
        intensity_data_reader = FinalReportGenotypeDataReader(
            args.final_report_file_path,
            sample_sheet["Sample_ID"].values,
            variants_to_read.as_df())

        intensity_data = intensity_data_reader.read_intensity_data()

        #intensity_data.columns.to_csv(args.out)

        intensity_data.to_pickle(args.out)

    if parser.is_action_requested(ArgumentParser.SubCommand.FIT):

        # Get batch correction configuration
        value_to_use = args.config['base']['value']
        intensity_correction_parameters = args.config['batch correction']

        # Load intensity data
        intensity_data_frame_reader = IntensityDataReader(sample_sheet["Sample_ID"])
        intensity_data_frame = intensity_data_frame_reader.load(args.input)

        intensity_data_frame[variants_in_locus.Name].to_csv(
            ".".join([args.out, "intensity_data_frame", "csv.gz"]),
            sep="\t", index_label='variant')

        # Intensity matrix
        intensity_matrix = intensity_data_frame.pivot(columns = "Sample ID", values = value_to_use)
        intensity_matrix.to_csv(
            ".".join([args.out, "mat", value_to_use.replace(" ", "_"), "csv"]),
            sep="\t", index_label='variant')

        value_to_use = "B Allele Freq"

        # Intensity matrix
        intensity_matrix = intensity_data_frame.pivot(columns = "Sample ID", values = value_to_use)
        intensity_matrix.to_csv(
            ".".join([args.out, "mat", value_to_use.replace(" ", "_"), "csv"]),
            sep="\t", index_label='variant')

        value_to_use = "Log R Ratio"

        intensity_matrix = intensity_data_frame.pivot(columns = "Sample ID", values = value_to_use)
        intensity_matrix.to_csv(
            ".".join([args.out, "mat", value_to_use.replace(" ", "_"), "csv"]),
            sep="\t", index_label='variant')

        # Do correction of intensities
        intensity_correction = IntensityCorrection(variants_in_locus.Name, **intensity_correction_parameters)
        intensity_correction.fit(intensity_matrix.T)

        # Write output for intensity correction
        intensity_correction.write_output(args.out)

    # if parser.is_action_requested(ArgumentParser.SubCommand.FIT):
    #
    #     if intensity_correction is None:
    #         intensity_correction = IntensityCorrection.load_instance(args.out)
    #
    #     corrected_intensities = intensity_correction.correct_intensities(intensity_data)



    if parser.is_action_requested(ArgumentParser.SubCommand.CALL):
        pass
    # args.

    # Perform method
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
