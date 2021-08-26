#!/usr/bin/env python

from classes.Hapmap import Hapmap
from plugins.bestmarkers import get_best_markers
from plugins.performance import check_performance
from plugins.differences import *
from plugins.plots import make_snp_plots
from classes.Data import Markermetadata
from utils.file_processing import *
from plugins.summary import get_marker_summary, get_sample_summary
from classes.LGC import LGC
from classes.Grid import Grid
from utils.definitions import *
from collections import defaultdict
import sys
import os
import warnings
import numpy as np
from classes.Getoptions import *
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 


def main():
    options = Getoptions(sys.argv[1:])
    print(options.task)
    task = options.task
    arguments = options.arguments
    
    if task == 'convert':
        in_format = arguments['I']
        out_format = arguments['O']
        in_file = arguments['i']
        out_file = arguments['o']
        if in_format == 'lgc':
            gt_data = LGC(in_file)
            msdata = gt_data.msdata
            smdata = gt_data.smdata
            smdata, msdata = fill_gaps_gtdata(smdata, msdata)
            marker_info = get_marker_info()
            markers = defaultdict(Markermetadata)
            markers = get_markers(in_file, markers, marker_info)
            if not markers:
                marker = make_markers(msdata, markers, marker_info)
        if out_format == "grid":
            write_grid_file(out_file, smdata, markers)
        if out_format == "fjk":
            write_flapjack_file(out_file, smdata, markers)
        if out_format == "hmp":
            write_hapmap_file(out_file, smdata, markers)
    
if __name__ == '__main__':
    main()
    

if False:    
        


    in_lgc_file = sys.argv[1]
    gt_data = Grid(in_lgc_file)
    msdata = gt_data.msdata
    smdata = gt_data.smdata

    smdata, msdata = fill_gaps_gtdata(smdata, msdata)

    marker_info = get_marker_info()
    markers = defaultdict(Markermetadata)
    samples = list(smdata.keys())
    markers = get_markers(in_lgc_file, markers, marker_info)

    if not markers:
        marker = make_markers(msdata, markers, marker_info)


    outdir = "../tests/output/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print(f"Created output directory {outdir}")

    # print("Writing Grid file..")
    # write_grid_file("output/out_grid.csv", smdata, markers)

    # print("Writing Flapjack file..")
    # write_flapjack_file("output/out_FJ.data", smdata, markers)

    # print("Writing Hapmap file..")
    # write_hapmap_file("output/out.hmp.txt", smdata, markers)

    # print("Writing SNP plots..")
    # make_snp_plots(gt_data.msdata, markers, outdir, gt_data.name)

    print(gt_data.name, len(smdata.keys()), len(msdata.keys()))

    print("Writing Marker Summary..")
    marker_summary = get_marker_summary(msdata)
    marker_summary.to_csv("output/marker_summary.txt", sep = "\t", index=False)

    print("Writing Sample Summary..")
    sample_summary = get_sample_summary(smdata)
    sample_summary.to_csv("output/sample_summary.txt", sep = "\t", index=False)

    # differences = find_differences(smdata)
    # differences.to_csv("output/differences.txt", sep = "\t", index=False)

    # print("Checking performance..")
    # performance = check_performance(smdata)
    # print(f"Total combinations: {performance[0]}")
    # print(f"Combinations with ZERO polymorphic markers: {performance[1]} ({round((performance[1]/performance[0])*100,2)} %)")
    # print(f"Combinations with >= 1 polymorphic markers: {performance[2]} ({round((performance[2]/performance[0])*100,2)} %)")
    # print(f"Combinations with >= 2 polymorphic markers: {performance[3]} ({round((performance[3]/performance[0])*100,2)} %)")

    # all_summary,ind_summary = get_best_markers(smdata, msdata, marker_summary)
    # all_summary.to_csv("output/BestMarkerSummaryAll.txt", sep = "\t", index=False)
    # ind_summary.to_csv("output/BestMarkerSummaryInd.txt", sep = "\t", index=False)

    # filename = "output/BestMarkerSummaryInd.txt"
    # ind_summary = pd.read_csv(filename, sep="\t", index_col='marker_count')
    # count = 15
    # print(list(ind_summary.loc[count]['marker']))
    # extract_markers = list(ind_summary.loc[count]['marker'])
    # extract_markers = extract_dict(markers, extract_markers)

    flt_sample_summary = sample_summary.loc[sample_summary['missing_percentage'] < SAMPLE_MISSING_CUTOFF]
    flt_sample_summary.to_csv("output/sample_summary_flt.txt", sep = "\t", index=False)
    samples = filter_samples(samples, flt_sample_summary)

    flt_marker_summary = marker_summary.loc[marker_summary['missing_percentage'] < MARKER_MISSING_CUTOFF]
    flt_marker_summary.to_csv("output/marker_summary_flt.txt", sep = "\t", index=False)
    markers = filter_markers(markers, flt_marker_summary)


    smdata = subset_gtdata(smdata, markers, samples, 'samplefast')
    msdata = subset_gtdata(msdata, markers, samples, 'markerfast')

    print("Writing Grid file..")
    write_grid_file("output/out_grid_flt.csv", smdata, markers)

    print("Writing Flapjack file..")
    write_flapjack_file("output/out_FJ_flt.data", smdata, markers)

    print("Writing Hapmap file..")
    write_hapmap_file("output/out_flt.hmp.txt", smdata, markers)


    # print("Checking performance..")
    # performance = check_performance(smdata)
    # print(f"Total combinations: {performance[0]}")
    # print(f"Combinations with ZERO polymorphic markers: {performance[1]} ({round((performance[1]/performance[0])*100,2)} %)")
    # print(f"Combinations with >= 1 polymorphic markers: {performance[2]} ({round((performance[2]/performance[0])*100,2)} %)")
    # print(f"Combinations with >= 2 polymorphic markers: {performance[3]} ({round((performance[3]/performance[0])*100,2)} %)")

    # # smdata, msdata = fill_gaps_gtdata(smdata, msdata)
    all_summary,ind_summary = get_best_markers(smdata, msdata, marker_summary)
    all_summary.to_csv("output/BestMarkerSummaryAll_flt.txt", sep = "\t", index=False)
    ind_summary.to_csv("output/BestMarkerSummaryInd_flt.txt", sep = "\t", index=False)

    # filename = "output/BestMarkerSummaryInd.txt"
    # ind_summary = pd.read_csv(filename, sep="\t", index_col='marker_count')
    # count = 15
    # print(list(ind_summary.loc[count]['marker']))
    # extract_markers = list(ind_summary.loc[count]['marker'])
    # extract_markers = extract_dict(markers, extract_markers)

    # tmp_smdata = subset_gtdata(smdata, extract_markers, samples, 'samplefast')
    # tmp_msdata = subset_gtdata(msdata, extract_markers, samples, 'markerfast')


    # print("Writing Hapmap file..")
    # write_hapmap_file("output/15_out.hmp.txt", tmp_smdata, extract_markers)

    # gt_data = Hapmap("output/15_out.hmp.txt")
    # print("Checking performance..")
    # performance = check_performance(gt_data.smdata)
    # print(f"Total combinations: {performance[0]}")
    # print(f"Combinations with ZERO polymorphic markers: {performance[1]}")
    # print(f"Combinations with >= 1 polymorphic markers: {performance[2]}")
    # print(f"Combinations with >= 2 polymorphic markers: {performance[3]}")





