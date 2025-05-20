"""
=========================================
Variant quality score recalibrator (VQSR)
=========================================
Author: Shujia Huang
Date  : 2014-05-23 11:21:53
"""
import sys

from . import variant_data_manager as VDM
from . import variant_recalibrator_engine as VRE


class VariantRecalibrator(object):
    def __init__ (self, vrac=None):
        self.vrac = vrac if vrac else VDM.VariantRecalibratorArgumentCollection()
        self.data_manger      = VDM.VariantDataManager()
        self.engine           = VRE.VariantRecalibratorEngine(self.vrac)
        self.bad_lod_cutoff   = None
        self.lod_cum_in_train = []

    def OnTraversalDone(self, data):
        self.data_manger.SetData(data)
        self.data_manger.NormalizeData()

        # Generate the positive model using the training data and evaluate each variant
        positive_training_data = self.data_manger.GetTrainingData()
        sys.stderr.write('[INFO] Training the goodModel ...\n')
        good_model = self.engine.GenerateModel(positive_training_data, self.vrac.MAX_GAUSSIANS)
        sys.stderr.write(f'[INFO] The converged information of goodModel is: {good_model.converged_}\n')
        sys.stderr.write(f'[INFO] The means of gaussion of goodModel is: {good_model.means_}\n\n')

        self.engine.EvaluateData(self.data_manger.data, good_model, False)
        self.bad_lod_cutoff, self.lod_cum_in_train = self.data_manger.CalculateWorstLodCutoff()

        # Generate the negative model using the worst performing data and 
        # evaluate each variant contrastively
        sys.stderr.write('[INFO] Training the badModel ...\n')
        negative_training_data = self.data_manger.SelectWorstVariants(self.bad_lod_cutoff)
        bad_model = self.engine.GenerateModel(negative_training_data, self.vrac.MAX_GAUSSIANS)
        sys.stderr.write(f'[INFO] The converged information of badModel is: {bad_model.converged_}\n')
        sys.stderr.write(f'[INFO] The means of gaussion of badModel is: {bad_model.means_}\n')
        self.engine.EvaluateData(self.data_manger.data, bad_model, True)

        if (not good_model.converged_) or (not bad_model.converged_): 
            raise ValueError ('[ERROR] NaN LOD value assigned. Clustering '
                              'with these variants and these annotations is '
                              'unsafe. Please consider raising the number of '
                              'variants used to train the negative model or '
                              'lowering the maximum number of Gaussians allowed '
                              'for use in the model.')

        # Find the VQSLOD cutoff values which correspond to the various 
        # tranches of calls requested by the user
        self.engine.CalculateWorstPerformingAnnotation(self.data_manger.data, good_model, bad_model)

    def VisualizationLodVStrainingSet(self, fig_name):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig = plt.figure()
        plt.title('LOD VS Positive training set', fontsize = 14)
        plt.plot(self.lod_cum_in_train[:,0], self.lod_cum_in_train[:,1], 'r-')
        plt.scatter(self.lod_cum_in_train[:,0], self.lod_cum_in_train[:,1],
                    c='r', marker='.', linewidth = 0, alpha = 0.5)
        plt.plot([self.bad_lod_cutoff, self.bad_lod_cutoff], [0,1], 'g--')
        plt.ylim(0, 1.0)
        plt.xlim(-10, 10)
        plt.xlabel('Variant score threshold for the bad model', fontsize = 16)
        plt.ylabel('Rate of Positive->Negative', fontsize = 16)
        fig.savefig(fig_name + '.png')

