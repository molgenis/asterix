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
import re
import copy

import numpy as np
import pandas as pd
import pyranges
import scipy.stats
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
DEFAULT_FINAL_REPORT_COLS = {"Sample ID": 'str', "SNP Name": 'str', "GType": 'str', "SNP": 'str',
                             "X": 'float', "Y": 'float', "B Allele Freq": 'float', "Log R Ratio": 'float'}
AUTOSOMES_CHR = ["{}".format(chrom) for chrom in range(1, 23)]


# Classes
class GenomicWindow:
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
    def __init__(self, variant_list_for_locus, pca_n_components=None,
                 pca_over_samples=True,
                 pca_scaling=True, regression_fit_intercept=True):
        self._reference_sample_means = None
        self._pca_over_samples = pca_over_samples
        self._target_variants = variant_list_for_locus
        self._scale = pca_scaling
        self.pca_n_components = pca_n_components
        self._correction_model = sklearn.linear_model.LinearRegression(
            fit_intercept=regression_fit_intercept)
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
        self._correction_model.fit(
            batch_effects, target_intensity_data_preprocessed)
        # Write intensities of locus of interest corrected for batch effects.
        corrected_intensities = self._correct_batch_effects(
            target_intensity_data_sliced, batch_effects)
        return corrected_intensities
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
        predicted_batch_effects_in_locus_of_interest = self._correction_model.predict(
            principal_components)
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

class SnpIntensityCnvCaller():
    """
    Class for calling copy numbers
    The class combines a series of strategies to perform cnv calling.
    - For each position, we perform a pca-like procedure on the intensities
    - Using the principal components, we apply gaussian mixture models on
    predefined ranges of the gene of interest.
    - For each sample we determine which (combination),
    of the predefined ranges performs
    _
    """
    def __init__(self, variants, clustering_ranges, ranges_for_dimensionality_reduction):

        self._variants = variants
        self._ranges_for_dimensionality_reduction = self.insert_indices(
            ranges_for_dimensionality_reduction)
        self._clustering_ranges = (
            clustering_ranges.insert(np.arange(len(clustering_ranges))))

        self._dimensionality_reduction_models = [None] * len(self._ranges_for_dimensionality_reduction)
        self._mixture_models = [None] * len(self._clustering_ranges)

        # For each dimensionality reduction range, get which ranges for the dimensionality reduction are involved
        self.dimensionality_reduction_range_mapping = self._ranges_for_dimensionality_reduction.join(
            self._variants,
            how='left', slack=1).as_df()[['Index', 'Index_b']]

        # For each clustering range, get which ranges for the dimensionality reduction are involved
        self._mixture_model_range_mapping = self._clustering_ranges.join(
            self._ranges_for_dimensionality_reduction,
            how='left', slack=1).as_df()[['Index', 'Index_b']]

    def fit_transform_dimensionality_reduction(self, intensities, cluster_assignment):

        intensities_projected = np.empty(
            shape=(intensities.shape[0], len(self._ranges_for_dimensionality_reduction)))

        for index, range in self._ranges_for_dimensionality_reduction.as_df().iterrows():

            range_index = range.Index
            indices_from_input = self._mixture_model_range_mapping[
                self._mixture_model_range_mapping['Index'] == range_index, 'Index_b']

            discriminant_analysis = sklearn.discriminant_analysis.LinearDiscriminantAnalysis(n_components=1)
            intensities_projected[:,index] = discriminant_analysis.fit_transform(intensities[:,indices_from_input], cluster_assignment)

            self._dimensionality_reduction_models[index] = discriminant_analysis

        return intensities_projected

    def fit_mixture_models(self, intensities, weights, centroids):
        # Range is a pyranges object.
        self._clustering_ranges.as_df().apply(
            self._fit_mixture_model,
            args=(intensities,), axis=1, centroids=centroids, weights=weights)

    def _fit_mixture_model(self, range, intensities, centroids, weights):
        # Slice the x_matrix columns to keep those that are in this range
        # x_matrix columns should match with the ranges in dimensionality reduction

        range_index = range.Index
        indices_from_input = self._mixture_model_range_mapping[
            self._mixture_model_range_mapping['Index'] == range_index, 'Index_b']

        centroids_per_feature = centroids[indices_from_input]

        mixture_model = sklearn.mixture.GaussianMixture(
            n_components=weights.shape[0], weights_init=weights[range_index], means_init=centroids_per_feature)
        mixture_model.fit(intensities[:,indices_from_input])

        self._mixture_models[intensities] = mixture_model
        return

    def mixture_model_probabilities(self, intensities):

        probabilities = self._clustering_ranges.as_df().apply(
            self._predict_mixture_model,
            args=(intensities,), axis=1)

        return probabilities

    def copy_numbers(self, intensities):
        probabilities = self.mixture_model_probabilities(intensities)

        # What to do when the mixture models are performed?
        # Per sample, per range, per class we have a probability.
        # Per sample we must then ascertain what cnv per range
        # is most probable.

        # We can do this by ascertaining for each sample what
        # combination of cnv calls is most probable

        pass

    def _predict_mixture_model(self, index, intensities):
        range_index = range.Index
        indices_from_input = self._mixture_model_range_mapping[self._mixture_model_range_mapping['Index'] == range_index, 'Index_b']

        probabilities = self._mixture_models[index].predict_proba(intensities[:,indices_from_input])

        return probabilities

    def insert_indices(self, ranges):
        return ranges.insert(pd.Series(
            np.arange(len(ranges)),
            "Index"))

    def write_output(self, path, projected_intensities, copy_number_probabilities):
        copy_number_probabilities.to_csv(
            ".".join([path, "copy_number_assignment", "probabilities", "csv", "gz"]))
        projected_intensities.to_csv(
            ".".join([path, "copy_number_assignment", "projected_intensities", "csv", "gz"]))

    def write_fit(self, path):
        pickle.dump(self, open(
            ".".join([path, "copy_number_assignment", "mod", "pkl"]), "wb"))

class NaiveHweInformedClustering:
    def __init__(self,
                 copy_number_allele_frequencies = None,
                 copy_number_allele_markers = None):

        if (copy_number_allele_frequencies is None
                or copy_number_allele_markers is None):
            copy_number_allele_frequencies = [0.0295, 0.0275]
            copy_number_allele_markers = [-1, 1]

        self.allele_table = pd.DataFrame({
            "freq": copy_number_allele_frequencies,
            "marker": copy_number_allele_markers
        })
        
        self.add_reference_allel()

    def get_copy_number_genotype_frequencies(self):
        allele_table = self.allele_table
        cartesian_allele_table = allele_table.join(
            allele_table, how='cross',
            lsuffix="_first", rsuffix="_second")

        cartesian_allele_table['freq_genotype'] = (
                cartesian_allele_table['freq_first'] * cartesian_allele_table['freq_second'])

        cartesian_allele_table['marker_genotype'] = (
                cartesian_allele_table['marker_first'] + cartesian_allele_table['marker_second'])

        genotype_table = cartesian_allele_table.groupby(
            ['marker_genotype']
        ).agg(
            {
                'freq_genotype': 'sum',
            }).reset_index()

        return(genotype_table)

    def add_reference_allel(self):
        ref_frequency = 1 - self.allele_table['freq'].sum()

        self.allele_table = pd.concat([
            self.allele_table,
            pd.DataFrame({
                "freq": [ref_frequency],
                "marker": [0]
            })
        ])

    def get_copy_number_assignment(self, x_matrix):
        column_axis = 1

        genotype_frequencies_ordered = (
            self.get_copy_number_genotype_frequencies()
                .sort_values(by=['marker_genotype']))

        cumulative_frequencies = np.cumsum(
            genotype_frequencies_ordered['freq_genotype'])

        ordered_genotype_markers = genotype_frequencies_ordered['marker_genotype']

        break_points = np.quantile(x_matrix, cumulative_frequencies, axis=0)

        genotype_marker_matrix = np.empty(x_matrix.shape, dtype=int)

        for column_index in range(x_matrix.shape[column_axis]):
            genotype_marker_matrix[:,column_index] = (
                ordered_genotype_markers[np.digitize(
                    x_matrix[:, column_index],
                    break_points[:,column_index],
                    right=True)]
            )

        sample_ranking = scipy.stats.rankdata(
            genotype_marker_matrix.sum(axis=column_axis),
            method='ordinal')

        # Get for every sample, get the appropriate genotype marker according to their rank
        copy_number_assignment = (
            ordered_genotype_markers[np.digitize(
                sample_ranking,
                np.quantile(sample_ranking, cumulative_frequencies),
                right=True)])

        return copy_number_assignment

    def get_centroids(self, x_matrix, copy_number_assignment=None, copy_number_frequencies=None):
        if copy_number_assignment is None:
            copy_number_assignment = self.get_copy_number_assignment(x_matrix).values

        if copy_number_frequencies is None:
            copy_number_frequencies = self.get_copy_number_genotype_frequencies()

        copy_number_frequencies.index = copy_number_frequencies['marker_genotype']

        copy_number_centroids = copy_number_frequencies.apply(
            lambda x:
                x_matrix[copy_number_assignment == x['marker_genotype']].mean(axis=0),
            axis=1, result_type='expand')

        return copy_number_centroids

    def write_output(self, path, naive_copy_number_assignment):
        naive_copy_number_assignment.to_csv(
            ".".join([path, "naive_clustering", "assignments", "csv", "gz"]))

# Functions
def ranges_from_file(file_path):
    return pyranges.PyRanges(pd.read_csv(file_path, index_col=False,
        names=("Chromosome", "Start", "End", "Name"),
        dtype={"Chromosome": str},
        sep="\t"))

def calculate_downsampling_factor(grouped_data_frame, N):
    grouped_data_frame['proportionsObserved'] = (
            grouped_data_frame.shape[0] / N)
    grouped_data_frame['downsamplingFactor'] = (
            grouped_data_frame.proportionsExpected / grouped_data_frame.proportionsObserved)
    return grouped_data_frame

def load_ranges_from_config(path_or_list_of_dicts):
    if isinstance(path_or_list_of_dicts, str):
        return pyranges.PyRanges(pd.read_csv(
            path_or_list_of_dicts, index_col=False,
            names=("Chromosome", "Start", "End"),
            dtype={"Chromosome": str},
            sep="\t"))
    return pyranges.PyRanges(
        pd.DataFrame(path_or_list_of_dicts, columns=["Chromosome", "Start", "End"]))

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

# Main

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
        sampled_corrective_variants = ranges_from_file("{}.corrective.bed".format(args.variants_prefix))

        # Read locus of interest
        variants_in_locus = ranges_from_file("{}.locus.bed".format(args.variants_prefix))

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

    intensity_data_frame_file_path = ".".join([args.out, "intensity_data_frame", "csv.gz"])

    if parser.is_action_requested(ArgumentParser.SubCommand.FIT):
        # Get batch correction configuration
        value_to_use = args.config['base']['value']
        intensity_correction_parameters = args.config['batch correction']

        # Load intensity data
        print("Loading intensity data...")
        intensity_data_frame_reader = IntensityDataReader(sample_list["Sample_ID"])
        intensity_data_frame = intensity_data_frame_reader.load(args.input)

        print("Intensity data loaded of shape: ".format(intensity_data_frame.shape), end=os.linesep*2)
        print("Writing intensity data to: {}    {}".format(os.linesep, intensity_data_frame_file_path))

        intensity_data_frame.loc[variants_in_locus.Name].to_csv(
            intensity_data_frame_file_path,
            sep="\t", index_label='variant')

        intensity_matrix_file_path = ".".join([args.out, "mat", value_to_use.replace(" ", "_"), "csv"])
        print("Writing matrix to: {}    {}".format(os.linesep, intensity_matrix_file_path), end=os.linesep*2)

        # Intensity matrix
        intensity_data_frame["Sample ID"] = pd.Categorical(intensity_data_frame["Sample ID"])
        intensity_matrix = intensity_data_frame.pivot(
            columns="Sample ID", values=value_to_use)
        intensity_matrix.to_csv(
            intensity_matrix_file_path,
            sep="\t", index_label='variant')

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
        ranges_for_dimensionality_reduction = load_ranges_from_config(
            args.config['dimensionality reduction'])
        clustering_ranges = load_ranges_from_config(
            args.config['clustering ranges'])

        variants_for_naive_clustering = variants_in_locus.overlap(ranges_for_naive_clustering)

        naive_clustering = NaiveHweInformedClustering()
        naive_copy_number_assignment = naive_clustering.get_copy_number_assignment(
            corrected_intensities[:,variants_for_naive_clustering.Name.values].values)
        naive_clustering.write_output(
            args.out, naive_copy_number_assignment)
        print("Naive clustering complete", end=os.linesep*2)

        print("Starting intensity CNV calling...")
        cyp2d6_intensity_cnv_caller = SnpIntensityCnvCaller(clustering_ranges=clustering_ranges,
                                                            ranges_for_dimensionality_reduction=ranges_for_dimensionality_reduction)

        print("Performing dimensionality reduction...")
        projected_intensities = cyp2d6_intensity_cnv_caller.fit_transform_dimensionality_reduction(
            corrected_intensities, naive_copy_number_assignment)
        print("Disentangling mixtures...")
        weights = np.tile(
            naive_clustering.get_copy_number_genotype_frequencies()['freq_genotypes'].values,
            (1, len(clustering_ranges)))
        cyp2d6_intensity_cnv_caller.fit_mixture_models(
            projected_intensities,
            weights,
            naive_clustering.get_centroids(projected_intensities).values)
        cyp2d6_cnv_probabilities = cyp2d6_intensity_cnv_caller.mixture_model_probabilities(
            projected_intensities)
        print("Mixtures calculated, samples scored.", end=os.linesep*2)
        print("Writing output to {}".format(args.out))
        cyp2d6_intensity_cnv_caller.write_output(
            args.out,
            projected_intensities=projected_intensities,
            copy_number_probabilities=cyp2d6_cnv_probabilities)
        cyp2d6_intensity_cnv_caller.write_fit(
            args.out
        )

    if parser.is_action_requested(ArgumentParser.SubCommand.CALL):
        # Get batch correction configuration
        value_to_use = args.config['base']['value']

        # Load intensity data
        intensity_data_frame_reader = IntensityDataReader(sample_list["Sample_ID"])
        intensity_data_frame = intensity_data_frame_reader.load(args.input)

        # Write intensity data
        intensity_data_frame.loc[variants_in_locus.Name].to_csv(
            intensity_data_frame_file_path,
            sep="\t", index_label='variant')

        # Intensity matrix
        intensity_matrix = intensity_data_frame.pivot(
            columns="Sample ID", values=value_to_use)

        # Do correction of intensities
        intensity_correction = IntensityCorrection.load_instance(args.correction)
        batch_effects = intensity_correction.batch_effects(
            reference_intensity_data=intensity_matrix.T)
        corrected_intensities = intensity_correction.correct_intensities(
            batch_effects=batch_effects,
            target_intensity_data=intensity_matrix.T)

        # Write output for intensity correction
        intensity_correction.write_output(
            args.out,
            corrected_intensities=corrected_intensities,
            batch_effects=batch_effects)

    # args.

    # Perform method
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
