import sys
import re
import gzip

import numpy as np
from sklearn.metrics import roc_curve

from . import vcfutils

class VariantRecalibratorArgumentCollection(object):
    def __init__ (self):
        self.NITER                 = 150
        self.NINIT                 = 100
        self.STD_THRESHOLD         = 8.0
        self.TRAIN_SIZE_RATE       = 0.6 # The ratio of Traing date set 
        self.CV_SIZE_RATE          = 0.2 # The ratio of cross validation data set 
        self.TEST_SIZE_RATE        = 0.2 # The ratio of test data set
        self.MAX_GAUSSIANS         = 10
        self.MIN_NUM_BAD_VARIANTS  = 1000
        self.MAX_NUM_TRAINING_DATA = 50000

        # The threshold that the positive trainingset -> negative
        self.POSITIVE_TO_NEGATIVE_RATE = 0.05
        self.MAX_NEG_GAUSSIANS = 10
     
        
class VariantDatum(object):
    def __init__ (self):
        self.annotations           = [] # Will be normalize and use for VQSR
        self.raw_annotations       = [] # Keep the raw value of each variant
        self.prior                 = 2.0
        self.lod                   = None 
        self.at_training_site      = False
        self.at_anti_training_site = False
        self.failing_std_threshold = False
        self.worst_annotation      = None
        self.variant_order         = None
        
        
class VariantDataManager(object):
    def __init__(self, data=None):
        self.vrac = VariantRecalibratorArgumentCollection()
        self.annotation_mean = None
        self.annotation_std  = None
        self.anno_texts      = [['QUAL', 'Float', 'Raw variant quality before VQSR process'],
                                ['DP', 'Integer', 'Total depth of this variant'],
                                ['FS', 'Float', 'Phred-scaled p-value using '
                                                'Fisher\'s exact test to detect strand bias'],
                                ['Indel_SP', 'Integer', 'Indel species around this position.'
                                                        'The less the better.'],
                                ['Indel_TOT', 'Integer', 'Number of Indel around this position.'
                                                         'The less the better.']]

        # The data is a list of VariantDatum
        self.data = self.SetData(data) if data else []

    def SetData(self, data):
        if not isinstance(data[0], VariantDatum):
            raise ValueError(f'[ERROR] The data type should be "VariantDatum" in '
                             f'VariantDataManager(), but found {type(data[0])}\n')

        for i, d in enumerate(data):
            data[i].annotations = np.array(d.annotations)
        
        return data

    def NormalizeData(self):
        data = np.array([d.annotations for d in self.data], dtype=float)
        mean = data.mean(axis=0)
        self.annotation_mean = mean

        std = data.std(axis=0)
        self.annotation_std = std

        # foundZeroVarianceAnnotation
        if any(std < 1e-5):
            raise ValueError('[ERROR] Found annotations with zero variance. '
                             'They must be excluded before proceeding.')

        # Each data now is (x - mean)/std
        for i, d in enumerate(data):
            self.data[i].annotations = (d - mean) / std
            # trim data by standard deviation threshold and mark failing data 
            # for exclusion later
            self.data[i].failing_std_threshold = False
            if any(np.abs(self.data[i].annotations) > self.vrac.STD_THRESHOLD):
                self.data[i].failing_std_threshold = True

    def GetTrainingData(self):
        training_data = [d for d in self.data if ((not d.failing_std_threshold) and d.at_training_site)]
        sys.stderr.write(f'[INFO] Training with {len(training_data)} variants after standard deviation '
                          'thresholding.\n')

        if len(training_data) < self.vrac.MIN_NUM_BAD_VARIANTS:
            sys.stderr.write(('[WARNING] Training with very few variant sites. Please check the model '
                              'reporting PDF to ensure the quality of the model is reliable.\n'))

        if len(training_data) > self.vrac.MAX_NUM_TRAINING_DATA:
            sys.stderr.write(f'[WARING] Very large training set detected. Downsampling to '
                             f'{self.vrac.MAX_NUM_TRAINING_DATA} training variants.\n')

            np.random.shuffle(training_data) # Random shuffling
            return training_data[:self.vrac.MAX_NUM_TRAINING_DATA]

        return training_data 

    def SelectWorstVariants(self, bad_lod):
        training_data = []
        for i,d in enumerate(self.data):
            if(d.lod < bad_lod) and (not d.failing_std_threshold):
                training_data.append(d)
                # I do need: i order to be the same as self.data
                self.data[i].at_anti_training_site = True

        sys.stderr.write(f'[INFO] Training with worst {len(training_data)} scoring '
                         f'variants --> variants with LOD < {bad_lod:.2f}\n')

        if len(training_data) > self.vrac.MAX_NUM_TRAINING_DATA:
            sys.stderr.write(f'[WARING] Very large training set detected. Downsampling to '
                             f'{self.vrac.MAX_NUM_TRAINING_DATA} training variants.\n')

            np.random.shuffle(training_data) # Random shuffling
            return training_data[:self.vrac.MAX_NUM_TRAINING_DATA]

        return training_data

    def CalculateWorstLodCutoff(self):
        lod_threshold, lod_cum = None, []
        if len(self.data) > 0:
            lod_dist = np.array([[d.at_training_site, d.lod] for d in self.data 
                                if(not d.failing_std_threshold)])

            # I just use the 'roc_curve' function to calculate the worst 
            # LOD threshold, not use it to draw ROC curve And 'roc_curve' 
            # function will output the increse order, so that I don't 
            # have to sort it again
            _, tpr, thresholds = roc_curve(lod_dist[:,0], lod_dist[:,1]) 
            lod_cum = [[thresholds[i], 1.0 - r] for i, r in enumerate(tpr)]
            for i, r in enumerate(tpr):
                if r > 1.0 - self.vrac.POSITIVE_TO_NEGATIVE_RATE: 
                    lod_threshold = round(thresholds[i])
                    break

        return lod_threshold, np.array(lod_cum)


def LoadTrainingSiteFromVCF(vcffile):
    """Record the training site
    """
    data_set = set()
    with gzip.open(vcffile, "rt") if vcffile.endswith(".gz") else open(vcffile, "r") as IN_VCF:
        for line in IN_VCF:
            if line.startswith('#'):
                continue

            col = line.strip().split()
            data_set.add(col[0] + ':' + col[1])  # just get the positions

    return data_set


def LoadDataSet(vcf_infile, traning_set):
    if len(traning_set) == 0: 
        raise ValueError('[ERROR] No Training Data found')

    data, h_info = [], vcfutils.Header()
    with gzip.open(vcf_infile, "rt") if vcf_infile.endswith(".gz") else open(vcf_infile, "r") as IN_VCF:
        for line in IN_VCF: # VCF format
            # Record the header information
            if line.startswith('#'):
                h_info.record(line.strip())
                continue

            col = line.strip().split()
            if col[3] in ['N', 'n']:
                continue

            qual = float(col[5])
            dp = re.search(r';?CM_DP=([^;]+)', col[7])
            fs = re.search(r';?FS=([^;]+)', col[7])
            indel_sp = re.search(r';?Indel_SP=([^;]+)', col[7])
            indel_tot = re.search(r';?Indel_TOT=([^;]+)', col[7])

            if any([not dp, not fs, not indel_sp, not indel_tot]):
                continue

            dp = round(float(dp.group(1)), 2)
            fs = round(float(fs.group(1)), 3)
            indel_sp = float(indel_sp.group(1))
            indel_tot = float(indel_tot.group(1))

            datum = VariantDatum()
            datum.raw_annotations = dict(QUAL=qual, DP=dp, FS=fs, Indel_SP=indel_sp, Indel_TOT=indel_tot)
            datum.annotations = [qual, dp, fs, indel_sp, indel_tot]

            datum.variant_order = col[0] + ':' + col[1]
            if datum.variant_order in traning_set:
                datum.at_training_site = True

            data.append(datum)

    return h_info, np.array(data)

