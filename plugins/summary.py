from utils.definitions import *
from utils.utils import *
import pandas as pd

def get_marker_summary(msdata):
    summary = pd.DataFrame(columns = ['marker_name', 
                                        'total_samples',
                                        'missing_count', 
                                        'missing_percentage',
                                        'major_allele',
                                        'major_allele_count',
                                        'major_allale_freq',
                                        'minor_allele',
                                        'minor_allele_count',
                                        'minor_allale_freq',
                                        'het_count',
                                        'het_propotion'])
    for marker in msdata:
        # [0] markers, 
        # [1] total_samples,
        # [2] missing_count, 
        # [3] missing_percentage
        # [4] major_allele
        # [5] major_allele_count
        # [6] major_allale_freq
        # [7] minor_allele
        # [8] minor_allele_count
        # [9] minor_allale_freq
        # [10] het_count
        # [11] het_freq
        counts = initialize_list(12)
        counts[0] = marker
        counts[1] = len(msdata[marker].data.keys())
        gt_calls = list(msdata[marker].data.values())
        gt_call_counts = get_counts(gt_calls)
        for call in MISSING_ALLELES:
            if call in gt_call_counts:
                counts[2] += gt_call_counts[call]
        counts[3] = round(counts[2]/counts[1]*100, 2)
        major_allele = allele_freq(gt_calls, ':')
        counts[4] = major_allele[0]
        counts[5] = major_allele[1]
        counts[6] = round(major_allele[1]/((counts[1]-counts[2])*2), 2)
        counts[7] = major_allele[2]
        counts[8] = major_allele[3]
        counts[9] = round(major_allele[3]/((counts[1]-counts[2])*2), 2)
        counts[10] = major_allele[4]
        counts[11] = round(major_allele[4]/(counts[1]-counts[2]), 2)
        summary.loc[len(summary.index)] = counts

        
    print(summary)
        
    