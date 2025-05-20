"""
Author: Shujia Huang
Date  : 2014-05-20 08:50:06
"""
import sys
import numpy as np
from sklearn.mixture import GaussianMixture as GMM
from scipy.special import logsumexp

from .variant_data_manager import VariantDatum
from .variant_data_manager import VariantRecalibratorArgumentCollection as VRAC


class VariantRecalibratorEngine(object):
    def __init__(self, vrac=None):
        self.vrac = vrac if vrac else VRAC()
        self.MIN_PROB_CONVERGENCE     = 2e-3
        self.MIN_ACCEPTABLE_LOD_SCORE = -2000.0

    def ClassifyData(self, data_size):
        """
        Classify the data into TrainingSet, Cross-ValidationSet and TestSet. 
        Reture the data indexes

        Call in GenerateModel
        """
        train_set_size = int(np.round(self.vrac.TRAIN_SIZE_RATE * data_size))
        test_set_size  = int(np.round(self.vrac.TEST_SIZE_RATE * data_size))
        cv_set_size    = int(np.round(self.vrac.CV_SIZE_RATE * data_size))

        # The index array of training data
        train_set_idx = range(train_set_size)

        # The index array of cross-validation data 
        cv_set_idx = range(train_set_size, cv_set_size + train_set_size)
        
        # The index array of Test data
        test_set_idx = range(cv_set_size + test_set_size, data_size)

        return train_set_idx, cv_set_idx, test_set_idx

    def GenerateModel(self, data, max_gaussians):
        if len(data) == 0:
            raise ValueError('[ERROR] No data found. The size is %d\n' % len(data))

        if not isinstance(data[0], VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in GenerateModel() of class VariantRecalibrato-'
                             'rEngine(), but found %s\n' % str(type(data[0])))

        if max_gaussians <= 0:
            raise ValueError(f'[ERROR] maxGaussians must be a positive integer '
                             f'but found: {max_gaussians}\n')

        gmms = [GMM(n_components=n + 1,
                    covariance_type='full',
                    tol=self.MIN_PROB_CONVERGENCE,
                    max_iter=self.vrac.NITER,
                    n_init=self.vrac.NINIT) for n in range(max_gaussians)]
        training_data = np.array([d.annotations for d in data])

        # find a best components for GMM model
        min_bic, bics = np.inf, []
        for n, g in enumerate(gmms, start=1):
            sys.stderr.write(f'[INFO] Trying {n} gaussian in GMM process training ...\n')
            g.fit(training_data)
            bic = g.bic(training_data)
            bics.append(bic)
            if bic == float('inf') or (bic < min_bic and g.converged_):
                bestgmm, min_bic = g, bic
                
            sys.stderr.write(f'  -- Converge infomation of training process: {g.converged_}')

        sys.stderr.write(f'[INFO] All the BIC: {bics}\n')
        sys.stderr.write('[INFO] Model Training Done. And take the model '
                         'with %d gaussiones which is the best with BIC %f.\n' %
                         (len(bestgmm.means_), min_bic))
        
        return bestgmm

    def EvaluateData(self, data, gmm, evaluate_contrastively=False):
        if not isinstance(data[0], VariantDatum):
            raise ValueError('[ERROR] The data type should be "VariantDatum" '
                             'in EvaluateData() of class VariantRecalibrator-'
                             'Engine(), but found %s\n' % str(type(data[0])))

        sys.stderr.write('[INFO] Evaluating full set of %d variants ...\n' % len(data))
        for i, _ in enumerate(data):
            # log likelihood and the base is 10
            this_lod = gmm.score(data[i].annotations[np.newaxis,:]) / np.log(10)
            if np.math.isnan(this_lod):
                gmm.converged_ = False
                return

            if evaluate_contrastively:
                # data[i].lod must has been assigned by good model or something like that.
                # contrastive evaluation: (prior + positive model - negative model)
                data[i].lod = data[i].prior + data[i].lod - this_lod
                if this_lod == float('inf'):
                    data[i].lod = self.MIN_ACCEPTABLE_LOD_SCORE * (1.0 + np.random.rand(1)[0])
            else:
                # positive model only so set the lod and return 
                data[i].lod = this_lod

        return self

    def CalculateWorstPerformingAnnotation(self, data, good_model, bad_model):
        for i, d in enumerate(data):
            prob_diff = [self.EvaluateDatumInOneDimension(good_model, d, k) -
                         self.EvaluateDatumInOneDimension(bad_model, d, k)
                        for k in range(len(d.annotations))]

            # Get the index of the worst annotations
            data[i].worst_annotation = np.argsort(prob_diff)[0]

        return self

    def EvaluateDatumInOneDimension(self, gmm, datum, iii):
        p_in_gaussian_logE = [
            np.log(w) + NormalDistributionLoge(gmm.means_[k][iii], gmm.covars_[k][iii][iii], datum.annotations[iii]) 
            for k, w in enumerate(gmm.weights_)
        ]
        # logsumexp is a numerically stable way to compute log(sum(exp(x)))
        return logsumexp(np.array(p_in_gaussian_logE))/np.log(10) # log10(Sum(pi_k * p(v|n,k)))


def NormalDistributionLoge(mu, sigma, x):
    if sigma <= 0:
        raise ValueError(f'[ERROR] sd: Standard deviation of normal must be > 0 but found: {sigma}\n')
    
    if (mu == float('inf') or mu == float('-inf') or
        sigma == float('inf') or sigma == float('-inf') or
        x == float('inf') or x == float('-inf')):
        raise ValueError('[ERROR] mean, sd, or, x: Normal parameters must '
                         'be well formatted (non-INF, non-NAN)')

    a = -1.0 * (np.log(sigma) + 0.5 * np.log(2 * np.pi))
    b = -0.5 * ((x - mu) / sigma) ** 2

    return a + b  # The Natural log
