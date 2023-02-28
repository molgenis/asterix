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

import argparse
# Standard imports.
import enum
import io
import os
import pickle
import re
import sys

import IlluminaBeadArrayFiles
import numpy as np
import pandas as pd
import pyranges
import scipy.stats
import sklearn.decomposition
import sklearn.discriminant_analysis
import sklearn.mixture
import yaml

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
from scipy.special import logsumexp
from scipy.special import factorial
from sklearn.cluster import MeanShift
from sklearn.mixture._gaussian_mixture import _estimate_gaussian_covariances_full, _estimate_log_gaussian_prob, \
    _compute_precision_cholesky, GaussianMixture
from sklearn.neighbors import KNeighborsClassifier

# Define a dictionary with all the columns types that are written to the final report files

DEFAULT_FINAL_REPORT_COLS = {"Sample ID": 'str', "SNP Name": 'str', "GType": 'str', "SNP": 'str',
                             "X": 'float', "Y": 'float', "B Allele Freq": 'float', "Log R Ratio": 'float'}
# Define a list of the autosomes to consider in the analysis
AUTOSOMES_CHR = ["{}".format(chrom) for chrom in range(1, 22)]


# Classes
class GenomicWindow:
    """
    Class that is responsible for handling genomic windows, specifically for argument parsing of parameters related
    to genomic windows.

    The use case is based on automatic conversion of a string to a formal genomic window class instance
    with correct data types. (an integer, and the associated units).

    Units can be bp (for basepairs), kb (for kilobases), mb (for megabases) and variants (number of variants).

    Makes use of pyranges package to be able to select the variants within the defined window.
    """
    BP_UNITS = {"bp": 1, "kb": 1000, "mb": 1000000}
    UNITS = tuple(BP_UNITS.keys()) + ("variants",)
    def __init__(self, window, unit="variants"):
        self.unit = unit
        self.window = window
    def get_variants(self, locus_ranges, manifest_ranges):
        if self.unit == "variants":
            return pyranges.concat([
                self.get_n_variants(locus_ranges, manifest_ranges),
                manifest_ranges.intersect(locus_ranges)[["Chromosome", "Start", "End", "Name"]]])
        else:
            return manifest_ranges.intersect(
                locus_ranges.extend(int(self.get_window("bp"))))[["Chromosome", "Start", "End", "Name"]]
    def get_n_variants(self, locus_ranges, manifest_ranges):
        return (pyranges.concat([
            locus_ranges[[]].k_nearest(
                manifest_ranges, how="downstream", overlap=False,
                k=self.window),
            locus_ranges[[]].k_nearest(
                manifest_ranges, how="upstream", overlap=False,
                k=self.window)])
            .new_position("swap")[["Chromosome", "Start", "End", "Name"]])
    @classmethod
    def from_string(cls, window_as_string):
        regex_match = re.fullmatch(r"(\d+)(\w+)?", window_as_string)
        if regex_match is None:
            raise ValueError("window declaration {} not valid".format(window_as_string))
        window, unit = regex_match.group(1, 2)
        cls._unit_check(unit)
        if unit is None:
            unit = "variants"
        return GenomicWindow(int(window), unit)
    @classmethod
    def _unit_check(cls, unit):
        if unit is not None and unit not in cls.UNITS and unit is not None:
            raise ValueError("window unit not valid. {} not in {}".format(unit, cls.UNITS))
    def __str__(self):
        return " ".join([self.window, self.unit])
    def get_window(self, unit):
        self._unit_check(unit)
        return self.window * (self.BP_UNITS[self.unit] / self.BP_UNITS[unit])


class ArgumentParser:
    """
    Class that is responsible for parsing command line arguments.

    Defines multiple subcommand-like features that each perform a separate step in the CNV calling procedure.
    """
    def __init__(self):
        self.sub_commands = list()
        self.parser = self.create_argument_parser()
        self.add_command_choice_argument(self.parser)
        self.add_bead_pool_manifest_argument(self.parser)
        self.add_sample_list_argument()
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
    def add_sample_list_argument(self):
        self.parser.add_argument('-s', '--sample-list', type=self.is_readable_file,
                                 required=True,
                                 default=None,
                                 help="List of samples to include. This should be the list of samples that"
                                      "passed quality control.")
    def add_final_report_path_argument(self, parser):
        parser.add_argument('-g', '--final-report-file-path', type=self.is_readable_file, required=True, default=None,
                            help="Path to where final report files are located")
    def add_out_argument(self, parser):
        parser.add_argument('-o', '--out', type=self.can_write_to_file_path,
                            required=True, default=None,
                            help="File prefix the output can be written to. ")
    def add_bed_path_parameter(self, parser):
        parser.add_argument('-b', '--bed-file', type=self.is_readable_file,
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
    def is_prefix_pointing_to_readables(cls, prefix, suffixes):
        """
        Checks whether the given directory is readable
        :param suffixes: a list of suffixes to check
        :param prefix: a prefix to a file in string format
        :return: file_path
        :raises: Exception: if the given path is invalid
        :raises: Exception: if the given directory is not accessible
        """
        for suffix in suffixes:
            file_path = "{}.{}".format(prefix, suffix)
            assert file_path == cls.is_readable_file(file_path)
        return prefix
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
    def add_variant_prefix_argument(self, parser):
        parser.add_argument(
            '-V', '--variants-prefix',
            type=lambda x: self.is_prefix_pointing_to_readables(x, ["corrective.bed", "locus.bed"]),
            help="matches .locus.bed and .corrective.bed files."
        )
    def extend_argument_parser(self):
        sub_command_mapping = {
            self.SubCommand.VARIANTS:
                {self.add_window_argument,
                 self.add_bed_path_parameter,
                 self.add_corrective_variants_argument},
            self.SubCommand.DATA:
                {self.add_final_report_path_argument,
                 self.add_variant_prefix_argument},
            self.SubCommand.FIT:
                {self.add_staged_data_argument,
                 self.add_variant_prefix_argument},
            self.SubCommand.CALL:
                {self.add_staged_data_argument,
                 self.add_variant_prefix_argument,
                 self.add_batch_weights_argument,
                 self.add_calling_cluster_weight_argument}}
        methods_to_run = set.union(*[sub_command_mapping.get(sub_command) for sub_command in self.sub_commands])
        if self.add_staged_data_argument in methods_to_run and self.add_staged_data_output_argument in methods_to_run:
            methods_to_run.discard(self.add_staged_data_argument)
        for method in methods_to_run:
            method(self.parser)
    def add_window_argument(self, parser):
        parser.add_argument('-w', '--window', type=GenomicWindow.from_string,
                            required=False, default=0,
                            help="number of variants or kb to extend the locus of interest"
                                 "for variants")
    def add_batch_weights_argument(self, parser):
        parser.add_argument('-C', '--correction', type=self.is_readable_dir,
                            required=True, default=None,
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
    """
    Loads intensity data, stored in pickle files, which has been prepared from final report files.
    """
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
            print("warning: {}".format(missing_samples))
        excess_samples = np.setdiff1d(sample_list, self._sample_list)
        if len(excess_samples) > 0:
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
        data_frame = self._empty_dataframe()
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
                    data_array_list.append(self._empty_dataframe())
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
    def _empty_dataframe(self):
        final_report_columns = DEFAULT_FINAL_REPORT_COLS.copy()
        final_report_columns['R'] = 'float'
        columns = {key: pd.Series(dtype=value) for key, value in final_report_columns.items()}
        empty_sample_data_frame = pd.DataFrame(columns)
        empty_sample_data_frame.set_index('SNP Name', inplace=True)
        return empty_sample_data_frame
    def _read_sample_intensities(self, buffer, columns):
        buffer.seek(0)
        sample_data_frame = pd.read_csv(buffer, names=columns,
                                        sep=self.sep,
                                        usecols=list(DEFAULT_FINAL_REPORT_COLS.keys()),
                                        dtype=DEFAULT_FINAL_REPORT_COLS,
                                        index_col="SNP Name")
        sample_data_frame = sample_data_frame.loc[
                            self._variants_to_include.Name.to_numpy(), :]
        sample_data_frame["R"] = sample_data_frame[["X", "Y"]].sum(axis=1).values
        if not np.all(self._variants_to_include.Name.to_numpy() == sample_data_frame.index.to_numpy()):
            raise FinalReportReaderException(
                ("Variants in variants file do not perfectly intersect with variants of sample {}"
                 .format(sample_data_frame["Sample ID"][0])), self._path, self._line_counter)
        return sample_data_frame


class IntensityCorrection:
    """
    Class that is responsible for intensity correction of raw intensities.
    """
    def __init__(self, variant_list_for_locus, pca_n_components=None,
                 pca_over_samples=True,
                 pca_scaling=True, regression_fit_intercept=True):
        self._regression_fit_intercept = regression_fit_intercept
        self._reference_sample_means = None
        self._pca_over_samples = pca_over_samples
        self._target_variants = variant_list_for_locus
        self._scale = pca_scaling
        self.pca_n_components = pca_n_components
        self._correction_models = [None] * len(self._target_variants)
        self._standardize_scaler = sklearn.preprocessing.StandardScaler()
        self._fitted_reference_variants = None
        self._pca_fit = None
        self._pca_explained_variance_dataframe = None
    def fit(self, batch_effects, target_intensity_data):
        # The projected principal components explain batch effects.
        # We try to explain as much of the locus of interest using the PCs
        # The residuals can be used in further analyses.
        target_intensity_data_sliced = target_intensity_data.loc[
                                       :, self._target_variants]
        target_intensity_data_preprocessed = pd.DataFrame(
            sklearn.preprocessing.StandardScaler(with_std=False)
                .fit_transform(target_intensity_data_sliced),
            columns=target_intensity_data_sliced.columns,
            index=target_intensity_data_sliced.index
        )
        # Fit the correction model
        self._correction_models = target_intensity_data_preprocessed.apply(
            lambda y: self._fit_correction_model(batch_effects, y))
        # Write intensities of locus of interest corrected for batch effects.
        corrected_intensities = self._correct_batch_effects(
            target_intensity_data_sliced, batch_effects)
        return corrected_intensities
    def _fit_correction_model(self, batch_effects, y):
        regression_model = sklearn.linear_model.LinearRegression(
            fit_intercept=self._regression_fit_intercept)
        nan_mask = y.notna()
        return regression_model.fit(batch_effects[nan_mask], y[nan_mask])
    def batch_effects_train(self, reference_intensity_data):
        # First prepare the reference intensity data.
        # On this reference intensity data, we base the batch effects.
        # Within the reference intensity data, we confine ourselves to the
        # variants that are not in the locus of interest
        reference_intensity_data_sliced = (
            reference_intensity_data.loc[:,
            self.variant_indices_outside_locus_of_interest(reference_intensity_data)])
        # Set the variant names we create a fit for.
        self._fitted_reference_variants = reference_intensity_data_sliced.columns.values
        # Now, if requested, we must scale the intensities of variants.
        if self._scale:
            # The intensity data matrix (transposed or not) we can center and scale.
            print("Scaling reference intensity data...")
            intensity_data_preprocessed = self._scale_fit_transform(
                reference_intensity_data_sliced)
        else:
            intensity_data_preprocessed = reference_intensity_data_sliced
        # Calculate the eigenvectors used to correct the correction variants.
        # These eigenvectors represent how to scale each variant in a sample
        # so that the result explaines covariability among the variants.
        # I.e., the projected principal components (PCs) explain batch effects.
        # We want to regress out these batch effects in the locus of interest.
        batch_effects = self._pca_fit_transform(intensity_data_preprocessed)
        return batch_effects
    def correct_intensities(self, batch_effects, target_intensity_data):
        # The principal components depict batch effects.
        # Here, we predict the batch effects on the locus of interest.
        # Using the predicted batch effects, we can correct the locus of interest for the
        # expected batch effects.
        target_intensity_data_sliced = (
            target_intensity_data.loc[
            :, self._target_variants[self._target_variants.isin(target_intensity_data.columns)]])
        residual_intensities = self._correct_batch_effects(
            target_intensity_data_sliced, batch_effects)
        return residual_intensities
    def batch_effects(self, reference_intensity_data):
        # First prepare the reference intensity data.
        # On this reference intensity data, we base the batch effects.
        # Within the reference intensity data, we confine ourselves to the
        # variants that are not in the locus of interest
        reference_intensity_data_sliced = (
            reference_intensity_data.loc[:,
            self.in_fitted_reference_variants(reference_intensity_data)])
        # Now, if requested, we must scale the intensities of variants.
        if self._scale:
            intensity_data_preprocessed = self._scale_transform(
                reference_intensity_data_sliced)
        else:
            intensity_data_preprocessed = reference_intensity_data_sliced
        # Get batch effects by calculating principal components
        batch_effects = self._pca_transform(
            intensity_data_preprocessed)
        return batch_effects
    def _scale_fit_transform(self, reference_intensity_data):
        return pd.DataFrame(
            self._standardize_scaler.fit_transform(
                reference_intensity_data),
            columns=reference_intensity_data.columns,
            index=reference_intensity_data.index)
    def _scale_transform(self, reference_intensity_data):
        return pd.DataFrame(
            self._standardize_scaler.transform(
                reference_intensity_data),
            columns=reference_intensity_data.columns,
            index=reference_intensity_data.index)
    def _pca_fit_transform(self, intensity_data_preprocessed):
        pca = sklearn.decomposition.PCA(
            n_components=self.pca_n_components)
        if self._pca_over_samples:
            print("Calculating principal components over samples")
            self._reference_sample_means = intensity_data_preprocessed.mean(axis=1)
            intensity_data_centered = (
                    intensity_data_preprocessed.T - self._reference_sample_means)
            # Now, if we should do a pca over the samples instead of
            # over the columns (variants), we transpose the
            # the intensity data frame, and fit the PCA on this matrix
            self._pca_fit = pd.DataFrame(
                pca.fit_transform(intensity_data_centered),
                index=self._fitted_reference_variants)
            # The components_ attribute (n_components, n_samples),
            # represent the batch effects.
            batch_effects = pd.DataFrame(
                pca.components_.T,
                index=intensity_data_preprocessed.index)
            print(batch_effects)
        else:
            # Now, if we should not do a pca over the samples,
            # but instead over the columns (variants), we fit
            # the pca object on the intensity data (not transposing)
            print("Calculating principal components")
            batch_effects = pd.DataFrame(
                pca.fit_transform(intensity_data_preprocessed),
                index=intensity_data_preprocessed.index)
            # We now assign the eigenvectors to the _pca_fit attribute
            self._pca_fit = pd.DataFrame(
                pca.components_.T,
                index=self._fitted_reference_variants)
        # Store the explained variance data in a dataframe
        self._pca_explained_variance_dataframe = (
            pd.DataFrame({'explained_variance': pca.explained_variance_,
                          'explained_variance_ratio': pca.explained_variance_ratio_}))
        return batch_effects
    def _pca_transform(self, intensity_data_preprocessed):
        if self._pca_over_samples:
            intensity_data_centered = (
                    intensity_data_preprocessed.T
                    - intensity_data_preprocessed.mean(axis=1))
            # We need to get the projections of the fit over samples now
            # We can calculate the batch effects for our samples using the following code,
            # wherein the x matrix represents the least squares solution to the
            # linear matrix equation.
            # (a @ x = b,
            # wherein a represents the intensity data,
            # x represents the batch effects (projected data of over sample pca),
            # and b depicts the pca fit (eigenvectors of the over sample pca))
            x, residuals, rank, s = np.linalg.lstsq(
                intensity_data_centered, self._pca_fit, rcond=None)
            # Assign x to be the batch effects
            batch_effects = pd.DataFrame(
                x, index=intensity_data_preprocessed.index)
        else:
            # replicate the transform method of sklearn.decomposition.PCA
            intensity_data_centered = (intensity_data_preprocessed
                                       - intensity_data_preprocessed.mean(axis=0))
            batch_effects = pd.DataFrame(
                np.dot(intensity_data_centered, self._pca_fit),
                index=intensity_data_preprocessed.index)
        return batch_effects
    def variant_indices_outside_locus_of_interest(self, intensity_data):
        return np.logical_and(
            self.indices_not_in_locus_of_interest(intensity_data),
            ~intensity_data.isnull().any(axis=0))
    def indices_not_in_locus_of_interest(self, intensity_data):
        return ~intensity_data.columns.isin(
            self._target_variants)
    def _correct_batch_effects(self, target_intensity_data_sliced, principal_components):
        # The principal components depict batch effects.
        # Here, we predict the batch effects on the locus of interest.
        # Using the predicted batch effects, we can correct the locus of interest for the
        # expected batch effects.
        print(target_intensity_data_sliced)
        predicted_batch_effects_in_locus_of_interest = np.array(
            [self._correction_models[index].predict(principal_components)
             for index, column in enumerate(target_intensity_data_sliced)]).T
        print(predicted_batch_effects_in_locus_of_interest.shape)
        print(np.mean(predicted_batch_effects_in_locus_of_interest, axis=None))
        # We can correct the locus of interest by subtracting the predicted batch effects
        # from the raw intensity data.
        residual_intensities = (
                target_intensity_data_sliced
                - predicted_batch_effects_in_locus_of_interest)
        return residual_intensities
    def write_output(self, path, corrected_intensities, batch_effects):
        corrected_intensities.to_csv(
            ".".join([path, "intensity_correction", "corrected", "csv", "gz"]))
        batch_effects.to_csv(
            ".".join([path, "intensity_correction", "batcheffects", "csv", "gz"]))
    def load_corrected_intensities(self, path):
        return pd.read_csv(
            ".".join([path, "intensity_correction", "corrected", "csv", "gz"]),
            index_col="Sample ID").rename_axis('SNP Name', axis=1)
    def write_fit(self, path):
        pickle.dump(self, open(
            ".".join([path, "intensity_correction", "mod", "pkl"]), "wb"))
        self._pca_fit.to_csv(
            ".".join([path, "intensity_correction", "pca", "csv", "gz"]))
        self._pca_explained_variance_dataframe.to_csv(
            ".".join([path, "intensity_correction", "eigenvalues", "csv", "gz"]))
    @classmethod
    def load_instance(cls, path):
        return pickle.load(open(
            ".".join([path, "intensity_correction", "mod", "pkl"]), "rb"))
    def outliers(self):
        raise NotImplementedError()
        # Loop through all requested principal components,
    # marking each sample that is outside of the mean -/+ x*sd
    def in_fitted_reference_variants(self, intensity_data):
        return np.logical_and(
            intensity_data.columns.isin(self._fitted_reference_variants),
            ~intensity_data.isnull().any(axis=0))


# class MultiDimensionalHweCalculator:
#     """
#     Class that is responsible for optimizing HWE for both a genotype dimension, a deletion and duplication dimension
#     """
#     def __init__(self):
#         pass
#
#     def calculate_expected_frequencies(self):
#         # Calculate p
#         # Calculate q
#
#         # Calculate -1
#         # Calculate +1
#
#         # Calculate +1:genotype effect
#         # Calculate -1:genotype effect


def calculate_theta(dataframe):
    theta = np.rad2deg(np.arctan(
        dataframe.xs("A", level="subspace", axis=1)/dataframe.xs("B", level="subspace", axis=1)))
    theta[theta < -45] = 180 + theta[theta < -45]
    return theta


class NaiveGenotypeClustering:
    """
    Class that performs a clustering of samples within a given variant to get a first approximation of genotypes
    present in the variant.

    The naive clustering attempts to separate genotypes within a given copy number
    """
    def __init__(self, copy_number_labels=None, ref_copy_number=2):
        if copy_number_labels is None:
            copy_number_labels = [-2, -1, 0, 1]
        self.ref_copy_number = ref_copy_number
        self.copy_number_labels = copy_number_labels
    def cluster_genotypes(self, intensities, resp):
        """
        Method that clusters samples to form genotypes.
        :param intensities:
        :param resp:
        :return: cluster to genotype assignment, genotype responsibilities
        """
        # Declare clustering table
        naive_genotype_clustering = {copy_number: None for copy_number in self.copy_number_labels}
        # Find median in -2 copy number assignment
        # Declare median table list
        medians = list()
        for index, assay_data_frame in intensities.groupby(level="SNP Name", axis=1):
            medians.append(assay_data_frame.iloc[(resp[:, 0] >= 0.95).astype(bool)].median())
        median_data_frame = pd.concat(medians)
        intensities_reanchored = intensities.subtract(median_data_frame)
        theta = (intensities_reanchored.groupby(level="SNP Name", axis=1, group_keys=False)
                 .apply(calculate_theta))
        genotype_assignment_list = list()
        cluster_to_genotype_list = list()
        print(np.unique(resp_to_assignments(resp), return_counts=True))
        # Loop through copy numbers.
        for copy_number_index, copy_number in enumerate(naive_genotype_clustering.keys()):
            print("cnv", copy_number_index, copy_number)
            if copy_number_index == 0:
                genotype_assignment_list.append(resp[:, copy_number_index][:,np.newaxis])
                cluster_to_genotype_list.append(np.array([[0], [0]]))
                continue
            # Use only individuals that fit to this copy number
            # Therefore, we create a mask.
            selection_mask = (resp[:, copy_number_index] >= 0.95).astype(bool)
            # Based on the mask, we select the proper theta
            theta_selected = theta.iloc[selection_mask]
            cluster_assignments = self.genotype_clustering(theta_selected)
            cluster_resp = assignments_to_resp(cluster_assignments)
            genotype_assignments = np.full((resp.shape[0], cluster_resp.shape[1]), 0)
            genotype_assignments[selection_mask, :] = cluster_resp
            genotype_assignment_list.append(
                genotype_assignments)
            cluster_to_genotype_list.append(
                self.assign_genotypes(theta_selected, cluster_assignments, copy_number))
        print("Cluster to genotype list")
        cluster_mapping = np.concatenate(cluster_to_genotype_list, axis=1)
        print(cluster_mapping)
        cluster_resp = np.concatenate(genotype_assignment_list, axis=1)
        print(cluster_resp)
        return (cluster_mapping,
                cluster_resp)
    def genotype_clustering(self, theta_selected):
        ms = MeanShift(bin_seeding=True, min_bin_freq=2,
                       bandwidth=12)
        ms.fit(X=theta_selected)
        # Get the unique labels and counts
        return ms.labels_
    def assign_genotypes(self, theta, cluster_assignments, copy_number):
        # Which genotypes are possible given the copy number
        # 0 => 00
        # 1 => A,B
        # 2 => 2A, AB, 2B
        # 3 => 3A, 2AB, A2B, 3B
        copy_number_count = copy_number + self.ref_copy_number
        candidate_genotypes = np.array(range(0, copy_number_count + 1))
        # Make an assignment of clusters to the various genotypes
        # average_theta = theta_selected[cluster_assignments]
        # Order the clusters according to the average theta,
        # and assign the highest theta to the most A genotype.
        labels = np.unique(cluster_assignments)
        n_genotype_clusters: int = len(labels)
        resp = assignments_to_resp(cluster_assignments)
        nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
        print("frequencies:", nk)
        theta_means = np.mean((np.dot(resp.T, theta) / nk[:, np.newaxis]), axis=1)
        print("Theta means:", theta_means)
        mean_indices = np.argsort(theta_means)
        print("Mean indices:", mean_indices)
        # If there is only one cluster assignment
        # we cannot deduce much from the allele frequencies.
        # Thus, we have to do assignment based on the theta.
        # Is the theta > 60, we assign the A genotype,
        # if the theta < 30, we assign the B genotype,
        # if else: raise exception.
        if n_genotype_clusters == 1:
            if theta_means[0] < 30:
                # B genotype
                best_insertion_pos = 0
            elif theta_means[0] > 60:
                # A genotype
                best_insertion_pos = len(candidate_genotypes) - 1
            else:
                raise ValueError("Uncertain cluster to genotype assignment")
        else:
            best_insertion_pos = self.hwe_informed_cluster_assignment(
                candidate_genotypes, copy_number_count, n_genotype_clusters, nk[mean_indices])
        # print("Best:", candidate_genotypes[
        #                best_insertion_pos:
        #                best_insertion_pos+n_genotype_clusters])
        # Mean indices gives the indices on how to sort the resp;
        # [0, 1] would mean that 0 should be first, 1 should be last
        # [1, 0] would mean that 1 should be first, 0 should be last
        reverting_array = np.empty((n_genotype_clusters,))
        reverting_array[mean_indices] = np.arange(n_genotype_clusters)
        # summing this with the best_insertion_pos results in the
        # genotypes that each resp part represents
        a_dosages = candidate_genotypes[best_insertion_pos:best_insertion_pos+n_genotype_clusters]
        b_dosages = (len(candidate_genotypes) - 1 - candidate_genotypes)
        b_dosages = b_dosages[best_insertion_pos:best_insertion_pos+n_genotype_clusters]
        genotypes = np.stack((a_dosages, b_dosages))[:, reverting_array]
        # print("Genotypes:")
        # print(genotypes)
        return genotypes
    def hwe_informed_cluster_assignment(
            self,
            candidate_genotypes,
            copy_number_count,
            n_genotype_clusters,
            nk):
        # Determine which means match with the possible genotypes,
        # so that hwe is matched.
        # Loop through all possibilities of fitting the theta_means in the possible genotypes
        chi_squared_min = np.Inf
        best_insertion_pos = 0
        print("copy number count:")
        print(copy_number_count)
        print("frequencies:", nk)
        print("looping through insertion positions to find best fitting position")
        for insertion_pos in range(0, len(candidate_genotypes) - n_genotype_clusters + 1):
            print("insertion position:")
            print(insertion_pos)
            # Make an array that lists the frequencies
            # for the samples currently considered.
            # First make an empty array, than fill
            # it with the appropriate frequencies
            frequencies_observed = np.array([0.0] * len(candidate_genotypes))
            nk_ = (nk / np.sum(nk))
            frequencies_observed[insertion_pos:insertion_pos + n_genotype_clusters] = nk_
            print(nk_)
            # Print the results
            print(frequencies_observed)
            print("This genotype according to this insertion position:")
            print(candidate_genotypes[
                  insertion_pos:
                  insertion_pos + n_genotype_clusters])
            print("Testing this setup...")
            # Now determine which matches HWE best
            chi_squared = self.hardy_chi_squared_extended(
                candidate_genotypes, copy_number_count, frequencies_observed)
            # Which insertion pos gives best HWE match
            if chi_squared < chi_squared_min:
                chi_squared_min = chi_squared
                best_insertion_pos = insertion_pos
        return best_insertion_pos
    def hardy_chi_squared_extended(self, candidate_genotypes, copy_number_count, frequencies_observed):
        # print("Starting Hardy Weinberg calculation")
        # print(candidate_genotypes)
        # print(copy_number_count)
        # print(frequencies_observed)
        freq_a = np.sum(candidate_genotypes * frequencies_observed) / len(candidate_genotypes)
        dosage_b = len(candidate_genotypes) - 1 - candidate_genotypes
        freq_b = np.sum((dosage_b) * frequencies_observed) / len(candidate_genotypes)
        # print(freq_a, freq_b)
        estimated = (
                np.power(freq_a, candidate_genotypes) *
                np.power(freq_b, dosage_b))
        # print(estimated)
        # print(dosage_b)
        permutation_possibilities = (
                factorial(copy_number_count) / (
                factorial(candidate_genotypes) *
                factorial(dosage_b)))
        # print(permutation_possibilities)
        frequencies_estimated = (estimated * permutation_possibilities)
        # print(frequencies_estimated)
        chi_squared = np.nansum((
                np.power(frequencies_observed - frequencies_estimated, 2) /
                frequencies_estimated))
        # print(chi_squared)
        return chi_squared
    def genotype(self, intensity_dataset, resp):
        for variant in intensity_dataset.variants():
            self.cluster_genotypes(variant.data_frame(), resp)


class NaiveHweInformedClustering:
    """
    Class that performs a clustering of samples within a given variant to get a first approximation of CNV status.

    The naive clustering makes use of expected frequencies of duplications and deletions,
    and assigns a proportional number of samples with the
    highest and lowest intensities to duplicated and deleted CNV status respectively.
    """
    COLUMN_INDEX = 1
    def __init__(
            self,
            copy_number_allele_frequencies = None,
            copy_number_allele_markers = None,
            genotype_labels = None,
            ref_genotype = 2
    ):
        self.ref_genotype = ref_genotype
        if (copy_number_allele_frequencies is None
                or copy_number_allele_markers is None):
            copy_number_allele_frequencies = [0.0295, 0.0275]
            copy_number_allele_markers = [-1, 1]
        self.allele_table = pd.DataFrame({
            "freq": copy_number_allele_frequencies,
            "marker": copy_number_allele_markers
        })
        self.add_reference_allel()
        self.genotypes = np.array(genotype_labels)
        self.genotype_label = "marker_genotype"
    def get_copy_number_genotype_frequencies(self):
        allele_table = self.allele_table
        cartesian_allele_table = allele_table.join(
            allele_table, how='cross',
            lsuffix="_first", rsuffix="_second")
        cartesian_allele_table['freq_genotype'] = (
                cartesian_allele_table['freq_first'] * cartesian_allele_table['freq_second'])
        cartesian_allele_table[self.genotype_label] = (
                cartesian_allele_table['marker_first'] + cartesian_allele_table['marker_second'])
        genotype_table = cartesian_allele_table.groupby(
            [self.genotype_label]
        ).agg(
            {
                'freq_genotype': 'sum',
            }).reset_index()
        genotype_table = genotype_table.loc[np.logical_or(
            self.genotypes is None,
            genotype_table[self.genotype_label].isin(np.array(self.genotypes))),]
        genotype_table['freq_genotype']/=genotype_table['freq_genotype'].sum()
        return genotype_table
    def add_reference_allel(self):
        ref_frequency = 1 - self.allele_table['freq'].sum()
        self.allele_table = pd.concat([
            self.allele_table,
            pd.DataFrame({
                "freq": [ref_frequency],
                "marker": [0]
            })
        ])
    def get_copy_number_assignment(self, intensities, genotype_frequencies=None):
        if genotype_frequencies is None:
            genotype_frequencies = self.get_copy_number_genotype_frequencies()
        genotype_frequencies_ordered = genotype_frequencies.sort_values(by=[self.genotype_label])
        cumulative_frequencies = np.cumsum(
            genotype_frequencies_ordered['freq_genotype'])
        ordered_genotype_markers = genotype_frequencies_ordered[self.genotype_label]
        genotype_marker_matrix = self.get_votes(
            cumulative_frequencies, intensities, ordered_genotype_markers)
        # Do at least half of the intensity variables agree?
        # Get the frequency of the most frequent item in the list
        voting_result = np.apply_along_axis(
            lambda x: self.vote(x),
            axis=self.COLUMN_INDEX,
            arr=genotype_marker_matrix)
        return pd.Series(voting_result, index=intensities.index, name=self.genotype_label)
    def get_votes(self, cumulative_frequencies, intensities, ordered_genotype_markers):
        break_points = np.quantile(intensities, cumulative_frequencies, axis=0)
        genotype_marker_matrix = np.empty(intensities.shape, dtype=int)
        for column_index in range(intensities.shape[self.COLUMN_INDEX]):
            genotype_marker_matrix[:, column_index] = (
                ordered_genotype_markers[np.digitize(
                    intensities.iloc[:, column_index],
                    break_points[:, column_index],
                    right=True)]
            )
        return genotype_marker_matrix
    def vote(self, votes):
        mode_result = scipy.stats.mode(votes, nan_policy='propagate')
        if mode_result.count[0] > (0.5 * len(votes)):
            return float(mode_result.mode[0])
        print(mode_result.count[0] > (0.5 * len(votes)))
        return np.NaN
    def get_centroids(self, intensities, copy_number_assignment=None, copy_number_frequencies=None):
        if copy_number_assignment is None:
            copy_number_assignment = self.get_copy_number_assignment(intensities).values
        if copy_number_frequencies is None:
            copy_number_frequencies = self.get_copy_number_genotype_frequencies()
        copy_number_frequencies.index = copy_number_frequencies[self.genotype_label]
        copy_number_centroids = copy_number_frequencies.apply(
            lambda x:
                intensities[copy_number_assignment == x[self.genotype_label]].mean(axis=0),
            axis=1, result_type='expand')
        print(copy_number_centroids)
        return copy_number_centroids
    def maximum_haplotype_counts(self):
        return dict(zip(self.genotypes, self.genotypes + self.ref_genotype + 1))
    def write_output(self, path, naive_copy_number_assignment):
        naive_copy_number_assignment.to_csv(
            ".".join([path, "naive_clustering", "assignments", "csv", "gz"]))
    def update_assignments_gaussian(self, intensities, genotype_frequencies=None):
        if genotype_frequencies is None:
            genotype_frequencies = self.get_copy_number_genotype_frequencies()
        assignments = self.get_copy_number_assignment(intensities, genotype_frequencies)
        # Resp
        # resp = np.equal.outer(assignments, self.genotypes).astype(float)
        resp = np.array(
            [(assignments == genotype) for genotype in self.genotypes]).astype(float).T
        # nk
        nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
        # means
        means = self.get_centroids(
            intensities,
            copy_number_assignment=assignments,
            copy_number_frequencies=genotype_frequencies).values
        # Covariances
        covariances = _estimate_gaussian_covariances_full(
            resp, intensities.values, nk, means, reg_covar=1e-6)
        # Precision
        precisions_cholesky = _compute_precision_cholesky(
            covariances, "full"
        )
        probabilities = (
            _estimate_log_gaussian_prob(intensities.values, means, precisions_cholesky, "full")
            + (np.log(genotype_frequencies['freq_genotype'].values)[np.newaxis, :]))
        labels_updated = probabilities.argmax(axis=1)
        assignments[...] = (
            labels_updated + min(genotype_frequencies[self.genotype_label]))
        return pd.DataFrame(probabilities, assignments.index, columns=self.genotypes), assignments


class SnpIntensityCnvCaller:
    """
    Class responsible for doing cnv calling.
    The class expects initial probabilities of CNV status per sample.
    The class attempts to do CNV calling in a number of steps.

    Initially, the class performs an initialization step.
    In the initialization step, a mean-shift clustering algorithm is performed. This aims to find clustering taking
    both X and Y intensities into account. This requires a bandwidth parameter that tells the clustering algorithm
    Which distance to consider between points for cluster association. Only clusters with more or equal points than a
    min_bin_freq are created.

    In the second step, the overlap between mean-shift clustering and the naive cnv clusters are calculated.
    This step aims to optimize the cluster assignment given the contraints of the expected clusters. For homozygous
    deletions for instance, only 1 cluster is expected at most, while at most 2 are expected for heterozygous deletions.
    """
    def __init__(self, counts, init_probabilities, frequency_min=5):
        self._components_map = dict()
        self._probabilities = init_probabilities
        self._weight = 0
        self._fitted_models = dict()
        self._frequency_min = frequency_min
        self._counts = counts
        self._genotype_clustering = NaiveGenotypeClustering()
    def fit(self, variant):
        # Expand probabilities to cover clusters
        # in this variant.
        mask = variant.mask
        masked_variant_dataframe = variant.data_frame().iloc[mask, :]
        # Get hard-call-assignments for each sample
        probabilities = self._probabilities[mask, :]
        assignments = resp_to_assignments(probabilities)
        # Based on neighbouring samples, assess if the CNV status for each sample is in-line with the
        # Naive clustering.
        sample_concordances = self.nearest_neighbour_concordance(masked_variant_dataframe, assignments)
        # Limiting ourselves to the samples that have a corresponding assignment, we
        # attempt to find the genotype clusters per CNV status.
        # We need to record both responsibility matrix,
        # as well as how each column translates to genotypes
        self._components_map[variant], resp_partial = self._genotype_clustering.cluster_genotypes(
            masked_variant_dataframe.iloc[sample_concordances, :],
            probabilities[sample_concordances, :])
        # Add missing rows in map
        print(resp_partial.shape)
        resp = np.full((probabilities.shape[0], resp_partial.shape[1]), 0).astype(float)
        resp[sample_concordances] = resp_partial
        print(resp.shape)
        print(resp.dtype)
        # Using the newly identified clusters,
        # optimize fit using a gaussian mixture model.
        mixtures = self._e_m(masked_variant_dataframe.values[sample_concordances, :], resp_partial)
        em_resp = mixtures.predict_proba(masked_variant_dataframe.values[sample_concordances, :])
        print(mixtures.converged_)
        print(mixtures.n_iter_)
        # print(em_resp)
        # print(resp_partial)
        em_assignments = resp_to_assignments(em_resp)[np.amax(em_resp, axis=1) > 0.99]
        partial_assignments = resp_to_assignments(resp_partial)[[np.amax(em_resp, axis=1) > 0.99]]
        matched_cluster = (em_assignments == partial_assignments)
        print(np.unique(partial_assignments, return_counts=True))
        print(np.unique(em_assignments, return_counts=True))
        self._fitted_models[variant] = mixtures
        # Somehow
    def fit_over_variants(self, dataset):
        for variant in dataset.variants():
            print("")
            print(variant.identifier)
            print("----------------")
            self.fit(variant)
    def predict(self, dataset):
        predicted = dict()
        for variant in dataset.variants():
            print("")
            print(variant.identifier)
            print("----------------")
            predicted[variant.identifier] = self.predict_variant(variant)
    def predict_variant(self, variant):
        # Expand probabilities to cover clusters
        # in this variant.
        mask = variant.mask
        masked_variant_dataframe = variant.data_frame().iloc[mask, :]
        em_resp = self._fitted_models[variant].predict_proba(
            masked_variant_dataframe.values)
        print(variant)
        pd.DataFrame(em_resp)
        self._components_map[variant]
        # Somehow
    def _e_m(self, values, resp):
        # Fit mixture model in this variant
        mixture_model = (
            IterativeGaussianMixture(
                n_components=resp.shape[1],
                resp_init=resp)
            .fit(values))
        return mixture_model
    def _estimate_centroids(self, values, resp):
        nk = (resp.sum(axis=0, keepdims=False)
              + 10 * np.finfo(resp.dtype).eps)
        centroids = np.dot(resp.T, values) / nk[:, np.newaxis]
        return centroids
    def write_fit(self, dataset, path):
        probabilities_out = list()
        probabilities_cnv = list()
        variant_identifiers = list()
        for variant in dataset.variants():
            mask = variant.mask
            samples = variant.data_frame().index[mask]
            # Write init responsibilities per variant
            variant_identifiers.append(variant._identifiers[0])
            # Get predicted probabilities
            predicted = (
                self._fitted_models[variant]
                .predict_proba(variant.values()[mask]))
            components_map = pd.DataFrame(self._components_map[variant].T, columns=["A", "B"])
            # Generate dataframe combining init and fitted probabilities
            probabilities_data_frame = (
                pd.melt(pd.DataFrame(predicted, index=samples).reset_index(),
                        id_vars=["Sample ID"],
                        var_name="Cluster_ID",
                        value_name='Probability'))
            out_data_frame = probabilities_data_frame.merge(
                components_map, how='inner',
                left_on='Cluster_ID', right_index=True,
                validate='m:1')
            print(out_data_frame)
            probabilities_out.append(out_data_frame)
            # Append probabilities mapped to CNVs
            # probabilities_cnv.append(pd.DataFrame(
            #     self._map_probabilities(predicted, components_map),
            #     index=samples, columns=self._counts))
        # Concatenate dataframes
        pd.concat(probabilities_out, keys=variant_identifiers,
                  names=['Variant', 'Row_ID'], ignore_index=False).to_csv(
            ".".join([path, "cnv_probabilities", "fitted", "csv", "gz"]))
        # pd.concat(probabilities_cnv, keys=variant_identifiers,
        #           names=['Sample_ID', 'variant']).to_csv(
        #     ".".join([path, "cnv_probabilities", "mapped", "csv", "gz"]))
    def nearest_neighbour_concordance(self, intensities, assignments):
        """
        Identifies the samples that have concordant cnv status according
        to a nearest neighbours approach.

        :param intensities: intensity dataframe for a variant.
        :return: boolean array indicating concordance.
        """
        # Perform nearest neighbours.
        neigh = KNeighborsClassifier(n_neighbors=5, weights='uniform')
        neigh.fit(intensities.values, y=assignments)
        proba = neigh.predict_proba(intensities.values)
        sample_concordance = resp_to_assignments(proba, threshold=0.8) == assignments
        return sample_concordance
    @classmethod
    def load_instance(cls, path):
        return pickle.load(open(
            ".".join([path, "cnv_calling", "mod", "pkl"]), "rb"))


class IterativeGaussianMixture(GaussianMixture):
    """
    Gaussian mixture model that can be initialized with initial responsibilities instead of initial means nd weights.
    This method should refrain from swapping mixture identities.

    Currently, there is no mechanism in place to prevent mixture identities.

    TODO: test if mixture identities change.
    TODO: if so, either fix the weights according to HWE informed weights,
    TODO: or recalculate weights, with the constraint of HWE.
    """
    def __init__(
            self,
            n_components=1,
            *,
            covariance_type="full",
            tol=1e-3,
            reg_covar=1e-6,
            max_iter=100,
            n_init=1,
            init_params="kmeans",
            resp_init=None,
            weights_init=None,
            means_init=None,
            precisions_init=None,
            random_state=None,
            warm_start=False,
            verbose=0,
            verbose_interval=10
    ):
        super().__init__(
            n_components=n_components,
            tol=tol,
            reg_covar=reg_covar,
            max_iter=max_iter,
            n_init=n_init,
            init_params=init_params,
            covariance_type=covariance_type,
            weights_init=weights_init,
            means_init=means_init,
            precisions_init=precisions_init,
            random_state=random_state,
            warm_start=warm_start,
            verbose=verbose,
            verbose_interval=verbose_interval,
        )
        self.resp_init=resp_init
    def _initialize_parameters(self, X, random_state):
        """
        Initializes parameters for
        :param X:
        """
        self._initialize(X, self.resp_init)


class ComplexIntensityDataset:
    """
    Variant intensity dataset class that is used for convenient slicing of specific variant intensities
    """
    def __init__(self, X, variant_map, feature_level="SNP Name"):
        self._feature_level = feature_level
        # Set X, with levels reordered to have the feature level first.
        self._X = pd.DataFrame.copy(X)
        levels = list(self._X.columns.names)
        print(levels)
        levels.remove(self._feature_level)
        levels.insert(0, self._feature_level)
        print(levels)
        self._X.columns.reorder_levels(levels)
        # Get complex features.
        self._variant_map = variant_map
        self._variants = self._setup_variants(self._variant_map)
    def get_intensities(self, identifiers):
        return self._X.loc[:,identifiers]
    def feature_labels(self):
        return self._X.columns.get_level_values(self._feature_level).unique()
    def _setup_variants(self, map):
        variants = dict()
        for key, feature_identifiers in enumerate(map):
            variants[key] = Variant(
                np.array(feature_identifiers), self)
        return variants
    def variants(self):
        for _, variant in self._variants.items():
            yield variant
    def write_dataset(self, path):
        for variant in self.variants():
            variant.data_frame().stack().to_csv(
                ".".join([path, "ab_intensity_dataset", variant.identifier, "csv", "gz"]))
    @property
    def n_samples(self):
        return self._X.shape[0]
    @property
    def n_variants(self):
        return len(self._variants)


class Variant:
    """
    Variant class that is used with the complex intensity dataset
    """
    def __init__(self, identifiers, intensity_provider):
        self._intensity_provider = intensity_provider
        self._identifiers = identifiers
        self._mask = self.mask
    def data_frame(self):
        return self._intensity_provider.get_intensities(self._identifiers)
    def values(self):
        return self.data_frame().values
    @property
    def identifier(self):
        return self._identifiers[0]
    @property
    def mask(self):
        return ~np.isnan(self.values()).any(axis=1)
    def __hash__(self):
        return hash(self._identifiers.tobytes())
    def __eq__(self, other):
        return (
            self.__class__ == other.__class__ and
            np.all(self._identifiers == other._identifiers)
        )


# Functions
def ranges_from_file(file_path):
    """
    Function that loads a pyranges object from a bedfile.
    :param file_path: File path to load the file from.
    :return: pyranges object representing the genomic ranges from within the given file.
    """
    return pyranges.PyRanges(pd.read_csv(file_path, index_col=False,
        names=("Chromosome", "Start", "End", "Name"),
        dtype={"Chromosome": str},
        sep="\t"))


def calculate_downsampling_factor(grouped_data_frame, N):
    """
    Function that calculates for a data frame how much each group should be downsampled to obtain a required number.
    :param grouped_data_frame: The dataframe group to calculate downsampling for.
    :param N: The target number of rows.
    :return: The grouped data frame, with downsamplingFactor column, representing the downsampling factor.
    """
    grouped_data_frame['proportionsObserved'] = (
            grouped_data_frame.shape[0] / N)
    grouped_data_frame['downsamplingFactor'] = (
            grouped_data_frame.proportionsExpected / grouped_data_frame.proportionsObserved)
    return grouped_data_frame


def load_ranges_from_config(path_or_list_of_dicts):
    """
    Load the genomic ranges as pyranges object, given either a path or a list of dicts.
    :param path_or_list_of_dicts: Path pointing to a bed file, or a list of dicts, representing bed file contents.
    :return: Pyranges.
    """
    ranges_data_frame = None
    if path_or_list_of_dicts is None:
        return None
    elif isinstance(path_or_list_of_dicts, str):
        ranges_data_frame = pd.read_csv(
            path_or_list_of_dicts, index_col=False, header=None, dtype={0: str}, sep="\t")
    else:
        ranges_data_frame = pd.DataFrame(path_or_list_of_dicts)
    ranges_data_frame.columns = np.array(["Chromosome", "Start", "End", "Name"])[0:ranges_data_frame.shape[1]]
    return pyranges.PyRanges(
        ranges_data_frame)


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
    filtered_chromosome_sizes = filtered_chromosome_sizes.assign(
        proportionsExpected=lambda d: d["End"] / np.sum(d["End"]))
    # We rename the data frame for easy merging.
    filtered_chromosome_sizes = filtered_chromosome_sizes.rename(
        columns={"Start": "ChromSizeStart", "End": "ChromSizeEnd"})
    # We need to select variants so that these are equally distributed across chromosomes.
    # We do this by sampling n variants in each chromosome,
    # where n denotes the proportional length of every chromosome, multiplied by the number of variants to
    # sample
    # Calculate what proportion of variants in each chromosome should be
    # discarded.
    corrective_variants_extended = (corrective_variants
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
    return pyranges.PyRanges(pd.concat(sampled_corrective_variants_list).loc[:, ["Chromosome", "Start", "End", "Name"]])


def assignments_to_resp(assignments, labels=None, astype=float):
    """
    Converts assignments (per sample assignment of label)
    to resp (per sample float assignment per component)

    :param assignments: 1D array of assignments to components
    :param labels: labels to create a column for

    :return: array-like of shape (n_samples, n_components)
    """
    if labels is None:
        labels = np.unique(assignments)
    resp = np.equal.outer(assignments, labels).astype(astype)
    return resp


def resp_to_assignments(resp, threshold=None):
    """
    Converts resp/probabilities (per sample float assignment per component)
    to assignments (per sample assignment of label)

    :param resp:
    :param threshold:
    :return: array-like of shape (n_samples)
    """
    assignments = resp.argmax(axis=1)
    if threshold is not None:
        assignments[np.amax(resp, axis=1) < threshold] = -1
    return assignments


# Main

def extract_intensity_matrices(intensity_data_frame, value_to_use, variants_in_locus):
    # Intensity matrix
    intensity_data_frame["Sample ID"] = pd.Categorical(
        intensity_data_frame["Sample ID"])
    intensity_matrix = intensity_data_frame.pivot(
        columns="Sample ID", values=value_to_use)
    a_matrix = intensity_data_frame.pivot(
        columns="Sample ID", values="X").T[variants_in_locus.Name]
    b_matrix = intensity_data_frame.pivot(
        columns="Sample ID", values="Y").T[variants_in_locus.Name]
    return a_matrix, b_matrix, intensity_matrix


def load_intensity_data_frame(in_directory, out_directory, sample_list, value_to_use):
    print("Loading intensity data...")
    intensity_data_frame_reader = IntensityDataReader(sample_list["Sample_ID"])
    intensity_data_frame = intensity_data_frame_reader.load(in_directory)
    intensity_data_frame_file_path = ".".join(
        [out_directory, "intensity_data_frame", "csv.gz"])
    intensity_matrix_file_path = (
        ".".join([out_directory, "mat", value_to_use.replace(" ", "_"), "csv"]))
    print("Intensity data loaded of shape: ".format(
        intensity_data_frame.shape), end=os.linesep * 2)
    print("Writing intensity data to: {}    {}".format(
        os.linesep, intensity_data_frame_file_path))
    print(
        "Writing matrix to: {}    {}"
        .format(os.linesep, intensity_matrix_file_path),
        end=os.linesep * 2)
    return intensity_data_frame


def get_corrected_complex_dataset(
        a_matrix, b_matrix, corrected_intensities, intensity_matrix, variants_in_locus, variant_selection):
    correction_factor = (
            corrected_intensities /
            intensity_matrix.T[variants_in_locus.Name])
    feature_matrix_a = a_matrix * correction_factor
    feature_matrix_b = b_matrix * correction_factor
    feature_matrix_ab = pd.concat(
        [feature_matrix_a, feature_matrix_b], keys=["A", "B"],
        names=["subspace", corrected_intensities.columns.name],
        axis=1).swaplevel(axis=1)
    variant_map = list(variant_selection.df.groupby("Start")["Name"].apply(tuple))
    return ComplexIntensityDataset(feature_matrix_ab, variant_map)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = ArgumentParser()
    args = parser.parse_input(argv[1:])

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

    manifest_data_frame[['Ref', 'Alt']] = (
        manifest_data_frame['Alleles']
            .str.split('\[(\w)/(\w)\]', expand=True)
            .iloc[:, [1, 2]])
    manifest_ranges = pyranges.PyRanges(manifest_data_frame)

    # Read the sample sheet
    sample_list = pd.read_csv(
        args.sample_list, sep=",", header=None, names=["Sample_ID"],
        dtype={"Sample_ID": str})

    if parser.is_action_requested(ArgumentParser.SubCommand.VARIANTS):

        # Read locus of interest
        locus_of_interest = pd.read_csv(
            args.bed_file, index_col=False,
            names=("Chromosome", "Start", "End", "Name"),
            dtype={"Chromosome": str},
            sep="\t")

        # Convert the locus of interest to a pyranges object
        locus_ranges = pyranges.PyRanges(locus_of_interest)

        # Get the intersect between the variants that are in the manifest and the locus of interest
        variants_in_locus = manifest_ranges.intersect(locus_ranges)

        window_from_locus = args.window

        variants_in_locus_extended = window_from_locus.get_variants(locus_ranges, manifest_ranges)

        manifest_ranges_locus_excluded = (
            pd.merge(manifest_ranges.as_df(), variants_in_locus.as_df(),
                     how="outer", indicator=True)
                .query('_merge == "left_only"'))

        # Sample corrective variants
        sampled_corrective_variants = sample_corrective_variants_proportionally(
            args.corrective_variants, manifest_ranges_locus_excluded)

        # variants as pyranges object
        sampled_corrective_variants.as_df().to_csv(
            "{}.corrective.bed".format(args.out), index=False, sep="\t", header=False)

        variants_in_locus_extended.as_df().to_csv(
            "{}.locus.bed".format(args.out), index=False, sep="\t", header=False)

    else:

        # Read locus of interest
        sampled_corrective_variants = ranges_from_file(
            "{}.corrective.bed".format(args.variants_prefix))

        # Read locus of interest
        variants_in_locus = ranges_from_file(
            "{}.locus.bed".format(args.variants_prefix))

    print(variants_in_locus)

    # Convert the locus of interest to a pyranges object
    variants_to_read = pyranges.concat([
        variants_in_locus,
        sampled_corrective_variants])

    variants_to_read = variants_to_read[~variants_to_read.Name.duplicated()]

    if parser.is_action_requested(ArgumentParser.SubCommand.DATA):
        # Get intensity data
        intensity_data_reader = FinalReportGenotypeDataReader(
            args.final_report_file_path,
            sample_list["Sample_ID"].values,
            variants_to_read.as_df())

        intensity_data = intensity_data_reader.read_intensity_data()

        intensity_data.to_pickle(args.out)

    if parser.is_action_requested(ArgumentParser.SubCommand.FIT):
        # Get batch correction configuration
        value_to_use = args.config['base']['value']
        intensity_correction_parameters = args.config['batch correction']
        cnv_calling_parameters = args.config['cnv calling']

        # Load intensity data
        print("Loading intensity data...")
        intensity_data_frame = load_intensity_data_frame(
            args.input, args.out, sample_list, value_to_use)

        a_matrix, b_matrix, intensity_matrix = extract_intensity_matrices(
            intensity_data_frame, value_to_use, variants_in_locus)

        print("Starting intensity correction...")

        # Do correction of intensities
        intensity_correction = IntensityCorrection(
            variants_in_locus.Name,
            **intensity_correction_parameters)
        print("Calculating PCA loadings for genome-wide batch effects...")
        batch_effects = intensity_correction.batch_effects_train(
            reference_intensity_data=intensity_matrix.T)
        print("Calculating batch effect corrections...")
        corrected_intensities = intensity_correction.fit(
            batch_effects=batch_effects,
            target_intensity_data=intensity_matrix.T)
        print("Intensity correction complete.", end=os.linesep*2)

        # Write output for intensity correction
        intensity_correction.write_output(
            args.out,
            corrected_intensities=corrected_intensities,
            batch_effects=batch_effects)
        intensity_correction.write_fit(args.out)

        # Get variants to use for initial clustering
        print("Starting naive clustering...")
        ranges_for_naive_clustering = load_ranges_from_config(
            args.config['naive clustering'])
        variant_selection = load_ranges_from_config(
            args.config['variant selection'])

        variants_for_naive_clustering = (
            variants_in_locus.overlap(ranges_for_naive_clustering))

        naive_clustering = NaiveHweInformedClustering(
            genotype_labels=[-2, -1, 0, 1], ref_genotype = 2)
        copy_number_frequencies = naive_clustering.get_copy_number_genotype_frequencies()
        naive_copy_number_assignment = naive_clustering.get_copy_number_assignment(
            corrected_intensities[variants_for_naive_clustering.Name.values],
            copy_number_frequencies)
        naive_clustering.write_output(
            args.out, naive_copy_number_assignment)
        print("Naive clustering complete", end=os.linesep*2)
        updated_naive_cnv_probabilities, updated_naive_cnv_assignment = (
            naive_clustering.update_assignments_gaussian(
            corrected_intensities[variants_for_naive_clustering.Name.values]))
        updated_naive_cnv_probabilities.to_csv(
            ".".join([args.out, "naive_clustering", "updated_probabilities", "csv", "gz"]))

        resp = assignments_to_resp(updated_naive_cnv_assignment.values)

        intensity_dataset = get_corrected_complex_dataset(
            a_matrix, b_matrix, corrected_intensities, intensity_matrix,
            variants_in_locus, variant_selection)

        intensity_dataset.write_dataset(args.out)
        #
        # naive_genotyping = NaiveGenotypeClustering(
        #     genotype_labels=[-2, -1, 0, 1], ref_genotype = 2)
        # naive_genotyping = naive_genotyping.genotype(intensity_dataset, resp)

        cnv_caller = SnpIntensityCnvCaller(
            (np.array([-2, -1, 0, 1]) + 2),
            resp)
        cnv_caller.fit_over_variants(intensity_dataset)
        SnpIntensityCnvCaller.write_fit(cnv_caller, intensity_dataset, args.out)
        #cnv_caller.write_fit(intensity_dataset, args.out)

    if parser.is_action_requested(ArgumentParser.SubCommand.CALL):
        # Get batch correction configuration
        value_to_use = args.config['base']['value']

        # Load intensity data
        intensity_data_frame = load_intensity_data_frame(
            args.input, args.out, sample_list, value_to_use)

        a_matrix, b_matrix, intensity_matrix = extract_intensity_matrices(
            intensity_data_frame, value_to_use, variants_in_locus)

        # Do correction of intensities
        intensity_correction = IntensityCorrection.load_instance(args.correction)
        batch_effects = intensity_correction.batch_effects(
            reference_intensity_data=intensity_matrix.T)
        corrected_intensities = intensity_correction.correct_intensities(
            batch_effects=batch_effects,
            target_intensity_data=intensity_matrix.T)
        print("Intensity correction complete.", end=os.linesep*2)

        # Write output for intensity correction
        intensity_correction.write_output(
            args.out,
            corrected_intensities=corrected_intensities,
            batch_effects=batch_effects)
        intensity_correction.write_fit(args.out)

        intensity_dataset = get_corrected_complex_dataset(a_matrix, b_matrix, corrected_intensities, intensity_matrix, variants_in_locus)

        cnv_caller = SnpIntensityCnvCaller.load_instance(args.fit)
        cnv_caller.predict(intensity_dataset)

        # TODO: Write output of cnv prediction.
    # Write output for intensity correction
        # intensity_correction.write_output(args.out,,

    # args.

    # Perform method
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
