"""
Created:      13/07/2022
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2022 C.A. Warmerdam

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
import numbers
import sys

import numpy as np
from scipy.special import logsumexp
from sklearn.cluster import MeanShift
from sklearn.mixture import GaussianMixture
from sklearn.mixture._gaussian_mixture import _compute_precision_cholesky, _estimate_log_gaussian_prob, \
    _estimate_gaussian_covariances_full

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

# Classes
from sklearn.utils import check_scalar, check_random_state


class ComplexFeatureDataset:
    def __init__(self, X, feature_level="SNP Name"):
        self._feature_level = feature_level
        self._X = X
        self._complex_features = self._setup_complex_features()
        self.shape = self._shape()
        print(self.shape)
    def _shape(self):
        feature_labels = self.feature_labels()
        return (self._X.shape[0], ) + feature_labels.shape
    def feature_labels(self):
        return self._X.columns.get_level_values(self._feature_level).unique()
    def _setup_complex_features(self):
        feature_label = self.feature_labels()
        return [self._setup_complex_feature(label) for label in feature_label]
    def _setup_complex_feature(self, label):
        feature = ComplexFeature(label, feature_level=self._feature_level)
        feature.missing(feature.X_(self._X).isna().any(axis=1))
        return feature
    def featurize_over_components(self, labels, max_features_across_components):
        for complex_feature in self._complex_features:
            complex_feature.initialize(
                complex_feature.X_(self._X), labels, max_features_across_components)
    def iter_nests(self):
        for complex_feature in self._complex_features:
            yield complex_feature.X_(self._X), complex_feature
    def featurization(self):
        identifiers = [complex_feature.identifier for complex_feature in self._complex_features]
        return pd.concat(
            [complex_feature.featurization(complex_feature.X_(self._X), self._X.index)
             for complex_feature in self._complex_features],
            axis=1, keys=identifiers, names=["feature", "component", "inner"])


class ComplexFeature:
    def __init__(self, identifier, feature_level="SNP Name"):
        self.identifier = identifier
        self._feature_level = feature_level
        self._inner_features = list()
        self._missing = None
        self._inner_model = None
        self._resp = None
    def initialize(self, X, sample_labels, subdivisions_max_per_component):
        if self._missing is None:
            self._missing = X.isna().any(axis=1)
        X_non_missing = X[~self._missing]
        sample_labels_non_missing = sample_labels[~self._missing]
        print(self.identifier)
        resp_per_component_list = list()
        for label, subdivisions_max in subdivisions_max_per_component.items():
            samples_pass = sample_labels_non_missing == label
            # filtered subdivision X
            component_filtered_x = X_non_missing[samples_pass]
            # theta
            theta = np.arctan(
                component_filtered_x.iloc[:, 0] /
                component_filtered_x.iloc[:, 1])
            # Resp
            resp_non_null = np.array([[1.0]*component_filtered_x.shape[0]]).T
            if subdivisions_max != 1:
                # Get the responsibilities per subdivision
                resp_non_null = self._resp_init(theta.values.reshape(-1, 1), subdivisions_max)
                # We need to update resp?
            samples_null = (sample_labels_non_missing == label)
            resp = np.zeros((
                X_non_missing.shape[0], # Number of samples with data for this complex feature
                resp_non_null.shape[1])) # Number of components in inner feature
            resp[samples_null] = resp_non_null
            # array-like of shape (n_samples_not_null, n_components*n_features)
            self._inner_features.append(tuple(range(resp_non_null.shape[1])))
            resp_per_component_list.append(resp)
        self._resp = np.concatenate(resp_per_component_list, axis=1)
        self._inner_model = InternalFeature(n_components=self._resp.shape[1], resp_init=self._resp)
        self._inner_model.fit(X_non_missing)
    def _resp_init(self, X, subdivisions_max, frequency_min=5):
        ms = MeanShift(bin_seeding=True, bandwidth=0.12, min_bin_freq = frequency_min)
        ms.fit(X = X)
        values, counts = np.unique(ms.labels_, return_counts=True)
        print(values, counts)
        # Now, select those values with a count lower than
        # The 0 - subdivisions_max value should be good enough.
        # However, if there are values equal to that, both should be removed
        partitioning_index = 0 - subdivisions_max
        frequency_threshold = frequency_min
        if len(counts) > subdivisions_max:
            partitioned = np.partition(counts, partitioning_index)[partitioning_index]
            frequency_threshold = partitioned if np.count_nonzero(counts == partitioned) == 1 else partitioned + 1
        values_filtered = values[counts >= np.max([frequency_threshold, frequency_min])]
        print(values_filtered)
        responsibilities = np.equal.outer(ms.labels_, values_filtered).astype(float)
        return responsibilities
    def X_(self, X_):
        return X_.xs(self.identifier, level=self._feature_level, axis=1)
    def e_m_round(self, X_, outer):
        X_non_missing = X_[~self._missing]
        self.update(outer[~self._missing])
        # means = self.estimate_means(X_non_missing)
        # diffs = self.diffs(X_non_missing, means)
        # covariances = _estimate_gaussian_covariances_full(
        #     self._resp, nk_inner, means, diffs, reg_covar)
        # precisions_cholesky = _compute_precision_cholesky(
        #     covariances, "full"
        # )
        # log_prob = (
        #         _estimate_log_gaussian_prob(X_, means, precisions_cholesky, "full"))
        # log_prob
    def resp(self):
        return self._resp
    def update(self, outer):
        # Resp for each internal component should be multiplied with the corresponding
        # Outer component, and the result should be normalised
        outer_repeated_log = np.repeat(
            np.log(outer[~self._missing]),
            repeats=[len(inner_comps) for inner_comps in self._inner_features],
            axis=0)
        resp_total_log = outer_repeated_log + np.log(self._resp)
        log_prob_norm = logsumexp(resp_total_log, axis=1)
        with np.errstate(under="ignore"):
            # ignore underflow
            self._resp = np.exp(resp_total_log - log_prob_norm[:, np.newaxis])
    def estimate_means_old(self, X_, resp):
        """
        Estimate means for each subfeature
        :param X_: intensities_matrix (n_samples, n_components)
        :param resp: array-like of shape (n_samples, n_components)
            The responsibilities for each data sample in X. 0-1
        :return: array-like of shape (n_components, n_subdivisions, n_subspace)
        """
        subdivision_means = list()
        # List of array-like. Shape (n_components, n_features)
        # Now loop through components
        resp_non_missing = resp[~self._missing].T
        for inner_component, outer_resp in zip(self._inner_features, resp_non_missing):
            log_resp = np.log(inner_component.resp)
            # Resp for inner features. [0-1]
            # Could have all 0s for a sample.
            # We should multiply the inner features, and normalise over all components*inner features
            resp_total = outer_resp[:,np.newaxis]
            if log_resp.shape[1] > 1:
                # Now, optimize centroids to fit responsibilities
                # Calculate common resp
                log_resp_sum = np.log(resp_total) + log_resp
                log_prob_norm = logsumexp(log_resp_sum, axis=1)
                with np.errstate(under="ignore"):
                    # ignore underflow
                    resp_total = np.exp(log_resp_sum - log_prob_norm[:, np.newaxis])
                    resp_total[np.isnan(resp_total)] = 0
            nk = (resp_total.sum(axis=0, keepdims=False)
                  + 10 * np.finfo(resp_total.dtype).eps)
            component_means = np.dot(resp_total.T, X_[~self._missing]) / nk[:, np.newaxis]
            print(component_means)
            # array-like of shape (n_subdivisions, n_subspace)
            subdivision_means.append(component_means)
        # Now we have calculated means for this component
        return subdivision_means
    def gaussian_parameters_featurized(self, X_):
        # Estimate means for each cluster
        means = self.estimate_means(X_)
        # Array-like of shape (n_components, n_subspace)
        # Estimate diffs for each cluster
        diffs = self.diffs(X_, means)
        # Array-like of shape (n_samples, n_components, n_subspace)
        means_list = (
            [means[np.array(inner_comps)+outer_comp,:]
             for outer_comp, inner_comps in enumerate(self._inner_features)])
        # Lists (features, components) of arrays (n_subspace)
        diffs_list = (
            [diffs[:,outer_comp+inner_comps,:]
             for outer_comp, inner_comps in enumerate(self._inner_features)])
        # Lists (features, components) of arrays (n_samples, n_subspace)
    def estimate_means(self, X_):
        nk = (self._resp.sum(axis=0, keepdims=False)
              + 10 * np.finfo(self._resp.dtype).eps)
        means = np.dot(self._resp.T, X_[~self._missing]) / nk[:, np.newaxis]
        return means
    def diffs(self, X_, means_over_inner_features):
        """
        Calculate the difference between X_ and means.

        :param X_: intensities matrix of shape (n_samples, n_subspace)
        :param means_over_inner_features: List with array-like of shape (n_inner_features, n_components, n_subspace)
        :return: difference between samples and means
        """
        # For the means
        diffs = X_[~self._missing][np.newaxis,...] - means
            # array-like of shape (n_components, n_samples, n_subspace)
        return diffs
    def featurization(self, X, index):
        return pd.concat(
            [pd.DataFrame(
                np.exp(inner_component._estimate_weighted_log_prob(X[~self._missing])),
                index=index[~self._missing])
             for inner_component in self._inner_features],
            names=["component", "inner"], axis=1,
            keys=range(len(self._inner_features)))
    def missing(self, missing):
        self._missing=missing


class InternalFeature(GaussianMixture):
    def __init__(
        self,
        n_components=1,
        *,
        covariance_type="full",
        tol=1e-3,
        reg_covar=1e-6,
        max_iter=0,
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
    def _check_initial_parameters(self, X):
        """Check values of the basic parameters.
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
        """
        check_scalar(
            self.n_components,
            name="n_components",
            target_type=numbers.Integral,
            min_val=1,
        )
        check_scalar(self.tol, name="tol", target_type=numbers.Real, min_val=0.0)
        check_scalar(
            self.n_init, name="n_init", target_type=numbers.Integral, min_val=1
        )
        check_scalar(
            self.max_iter, name="max_iter", target_type=numbers.Integral, min_val=0
        )
        check_scalar(
            self.reg_covar, name="reg_covar", target_type=numbers.Real, min_val=0.0
        )
        # Check all the parameters values of the derived class
        self._check_parameters(X)
    def fit_predict(self, X, y=None):
        """Estimate model parameters using X and predict the labels for X.
        The method fits the model n_init times and sets the parameters with
        which the model has the largest likelihood or lower bound. Within each
        trial, the method iterates between E-step and M-step for `max_iter`
        times until the change of likelihood or lower bound is less than
        `tol`, otherwise, a :class:`~sklearn.exceptions.ConvergenceWarning` is
        raised. After fitting, it predicts the most probable label for the
        input data points.
        .. versionadded:: 0.20
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.
        y : Ignored
            Not used, present for API consistency by convention.
        Returns
        -------
        labels : array, shape (n_samples,)
            Component labels.
        """
        if X.shape[0] < self.n_components:
            raise ValueError(
                "Expected n_samples >= n_components "
                f"but got n_components = {self.n_components}, "
                f"n_samples = {X.shape[0]}"
            )
        self._check_initial_parameters(X)
        # if we enable warm_start, we will have a unique initialisation
        do_init = not (self.warm_start and hasattr(self, "converged_"))
        n_init = self.n_init if do_init else 1
        max_lower_bound = -np.inf
        self.converged_ = False
        random_state = check_random_state(self.random_state)
        n_samples, _ = X.shape
        for init in range(n_init):
            self._print_verbose_msg_init_beg(init)
            if do_init:
                self._initialize_parameters(X, random_state)
            lower_bound = -np.inf if do_init else self.lower_bound_
            n_iter = 0
            for n_iter in range(1, self.max_iter + 1):
                prev_lower_bound = lower_bound
                log_prob_norm, log_resp = self._e_step(X)
                self._m_step(X, log_resp)
                lower_bound = self._compute_lower_bound(log_resp, log_prob_norm)
                change = lower_bound - prev_lower_bound
                self._print_verbose_msg_iter_end(n_iter, change)
                if abs(change) < self.tol:
                    self.converged_ = True
                    break
            self._print_verbose_msg_init_end(lower_bound)
            if lower_bound > max_lower_bound or max_lower_bound == -np.inf:
                max_lower_bound = lower_bound
                best_params = self._get_parameters()
                best_n_iter = n_iter
        self._set_parameters(best_params)
        self.n_iter_ = best_n_iter
        self.lower_bound_ = max_lower_bound
        # Always do a final e-step to guarantee that the labels returned by
        # fit_predict(X) are always consistent with fit(X).predict(X)
        # for any value of max_iter and tol (and any random_state).
        _, log_resp = self._e_step(X)
        return log_resp.argmax(axis=1)


class NestedGaussianMixture(GaussianMixture):
    def __init__(
            self,
            n_components=1,
            *,
            covariance_type="full",
            tol=1e-3,
            reg_covar=1e-6,
            max_iter=0,
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
    def _validate_data(
            self,
            X="no_validation",
            *args, **kwargs):
        pass
        # super()._validate_data(X="no_validation", *args, **kwargs)
    def _initialize(self, X, resp):
        """Initialization of the Gaussian mixture parameters.
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
        resp : array-like of shape (n_samples, n_components)
        """
        n_samples, _ = X.shape
        weights, means, covariances = self._estimate_gaussian_parameters(
            X, resp, self.reg_covar, self.covariance_type
        )
        weights /= n_samples
        self.weights_ = weights if self.weights_init is None else self.weights_init
        self.means_ = means if self.means_init is None else self.means_init
        if self.precisions_init is None:
            self.covariances_ = covariances
            self.precisions_cholesky_ = _compute_precision_cholesky(
                covariances, self.covariance_type
            )
        else:
            raise NotImplementedError()
    def _estimate_log_prob(self, X):
        log_prob = _estimate_log_gaussian_prob(X, self.means_, self.precisions_cholesky_, self.covariance_type)
        # Now we perhaps should make sure that the probabilities of features in the same subspace are
        # not overly used.
        return log_prob
    def _estimate_gaussian_parameters(self, X, resp, reg_covar, covariance_type):
        """Estimate the Gaussian distribution parameters.
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input data array.
        resp : array-like of shape (n_samples, n_components)
            The responsibilities for each data sample in X. 0-1
        reg_covar : float
            The regularization added to the diagonal of the covariance matrices.
        covariance_type : {'full', 'tied', 'diag', 'spherical'}
            The type of precision matrices.
        Returns
        -------
        nk : array-like of shape (n_components,)
            The numbers of data samples in the current components.
        means : array-like of shape (n_components, n_features)
            The centers of the current components.
        covariances : array-like
            The covariance matrix of the current components.
            The shape depends of the covariance_type.
        """
        if covariance_type != "full":
            raise NotImplementedError()
        nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
        means_list = list()
        diffs_list = list()
        for X_, complex_feature in X.iter_nests():
            complex_feature.e_m_round(X_)
            means, diffs = complex_feature.gaussian_parameters_featurized(X_)
            means_list.extend(means)
            diffs_list.extend(diffs)
        # List containing array-like of shape (n_samples, n_components, n_inner_features)
        means = np.concatenate([means for means in means_list], axis=1)
        diffs = np.concatenate([diffs for diffs in diffs_list], axis=2)
        # array-like of shape (n_components, n_features_expanded)
        covariances = self._estimate_gaussian_covariances_full(resp, nk, means, diffs, reg_covar)
        return nk, means, covariances
    def _estimate_gaussian_covariances_full(self, resp, nk, means, diffs_over_components, reg_covar):
        """Estimate the full covariance matrices.
        Parameters
        ----------
        resp : array-like of shape (n_samples, n_components)
        X : array-like of shape (n_samples, n_features)
        nk : array-like of shape (n_components,)
        means : array-like of shape (n_components, n_features)
        reg_covar : float
        Returns
        -------
        covariances : array, shape (n_components, n_features, n_features)
            The covariance matrix of the current components.
        """
        n_components, n_features = means.shape
        covariances = np.empty((n_components, n_features, n_features))
        for k in range(n_components):
            diffs = diffs_over_components[k]
            # array-like of shape (n_samples, n_features_expanded)
            covariances[k] = np.dot(resp[:, k] * diffs.T, diffs) / nk[k]
            covariances[k].flat[:: n_features + 1] += reg_covar
        return covariances
    def fit_predict(self, X, y=None):
        """Estimate model parameters using X and predict the labels for X.

        The method fits the model n_init times and sets the parameters with
        which the model has the largest likelihood or lower bound. Within each
        trial, the method iterates between E-step and M-step for `max_iter`
        times until the change of likelihood or lower bound is less than
        `tol`, otherwise, a :class:`~sklearn.exceptions.ConvergenceWarning` is
        raised. After fitting, it predicts the most probable label for the
        input data points.

        .. versionadded:: 0.20

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        labels : array, shape (n_samples,)
            Component labels.
        """
        if X.shape[0] < self.n_components:
            raise ValueError(
                "Expected n_samples >= n_components "
                f"but got n_components = {self.n_components}, "
                f"n_samples = {X.shape[0]}"
            )
        self._check_initial_parameters(X)
        # if we enable warm_start, we will have a unique initialisation
        do_init = not (self.warm_start and hasattr(self, "converged_"))
        n_init = self.n_init if do_init else 1
        max_lower_bound = -np.inf
        self.converged_ = False
        random_state = check_random_state(self.random_state)
        n_samples, _ = X.shape
        for init in range(n_init):
            self._print_verbose_msg_init_beg(init)
            if do_init:
                self._initialize_parameters(X, random_state)
            lower_bound = -np.inf if do_init else self.lower_bound_
            for n_iter in range(1, self.max_iter + 1):
                prev_lower_bound = lower_bound
                log_prob_norm, log_resp = self._e_step(X)
                self._m_step(X, log_resp)
                lower_bound = self._compute_lower_bound(log_resp, log_prob_norm)
                change = lower_bound - prev_lower_bound
                self._print_verbose_msg_iter_end(n_iter, change)
                if abs(change) < self.tol:
                    self.converged_ = True
                    break
            self._print_verbose_msg_init_end(lower_bound)
            if lower_bound > max_lower_bound or max_lower_bound == -np.inf:
                max_lower_bound = lower_bound
                best_params = self._get_parameters()
                best_n_iter = n_iter
        self._set_parameters(best_params)
        self.n_iter_ = best_n_iter
        self.lower_bound_ = max_lower_bound
        # Always do a final e-step to guarantee that the labels returned by
        # fit_predict(X) are always consistent with fit(X).predict(X)
        # for any value of max_iter and tol (and any random_state).
        _, log_resp = self._e_step(X)
        return log_resp.argmax(axis=1)



def assignments_to_resp(assignments, labels=None):
    """
    Converts assignments (per sample assignment of label)
    to resp (per sample float assignment per component)

    :param assignments: 1D array of assignments to components
    :param labels: labels to create a column for

    :return: array-like of shape (n_samples, n_components)
    """
    if labels is None:
        labels = np.unique(assignments)
    resp = np.equal.outer(assignments, labels).astype(float)
    return resp

# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    # Perform method
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
