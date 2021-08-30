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
import warnings
import numpy as np
from classes.Getoptions import *
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


def main():
    # if __name__ == '__main__':
    # main()

    parser = get_opts()
    options = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        exit(1)

    marker_info = get_marker_info()
    markers = defaultdict(Markermetadata)
    out_prefix = options.out_prefix

    if options.lgc_file:
        lgc_file = options.lgc_file
        gt_data = LGC(lgc_file)
        msdata = gt_data.msdata
        smdata = gt_data.smdata
        markers = get_markers(lgc_file, markers, marker_info)
        samples = list(smdata.keys())

    if options.grid_file:
        grid_file = options.grid_file
        gt_data = Grid(grid_file)
        msdata = gt_data.msdata
        smdata = gt_data.smdata
        marker = make_markers(msdata, markers, marker_info)
        samples = list(smdata.keys())

    smdata, msdata = fill_gaps_gtdata(smdata, msdata)

    if options.sample_list_file:
        sample_query = get_list_from_file(options.sample_list_file)
    else:
        sample_query = list(smdata.keys())

    if options.marker_list_file:
        marker_query = get_list_from_file(options.marker_list_file)
    else:
        marker_query = markers.copy()

    if options.sample_list_file or options.marker_list_file:
        print_log(f"Subsetting genotype data...")
        smdata = subset_gtdata(smdata, marker_query,
                               sample_query, 'samplefast')
        msdata = subset_gtdata(msdata, marker_query,
                               sample_query, 'markerfast')

    if options.task_run_summary:
        print_log("Writing Marker Summary..")
        marker_summary = get_marker_summary(msdata)
        marker_summary.to_csv(out_prefix + "_marker_summary.txt",
                              sep="\t", index=False)

        print_log("Writing Sample Summary..")
        sample_summary = get_sample_summary(smdata)
        sample_summary.to_csv(out_prefix + "_sample_summary.txt",
                              sep="\t", index=False)

    if options.task_filter:
        print_log(f"Filtering data")
        flt_sample_summary = sample_summary.loc[(sample_summary['missing_percentage']
                                                < options.cutoff_mx_missing_sample)]
        samples = filter_samples(samples, flt_sample_summary)
        flt_marker_summary = marker_summary.loc[(marker_summary['missing_percentage'] <= options.cutoff_mx_missing_marker)
                                                & (marker_summary['minor_allale_freq'] >= options.cutoff_mn_maf_marker)
                                                & (marker_summary['PIC'] >= options.cutoff_mn_pic_marker)]
        markers = filter_markers(markers, flt_marker_summary)
        smdata = subset_gtdata(smdata, markers, samples, 'samplefast')
        msdata = subset_gtdata(msdata, markers, samples, 'markerfast')
        if options.task_run_summary:
            print_log("Writing Filtered Sample Summary..")
            flt_sample_summary.to_csv(
                out_prefix + "_sample_summary_flt.txt", sep="\t", index=False)
            print_log("Writing Filtered Marker Summary..")
            flt_marker_summary.to_csv(
                out_prefix + "_marker_summary_flt.txt", sep="\t", index=False)

    if options.sample_list_a_file and options.sample_list_b_file:
        sample_list_a = get_list_from_file(options.sample_list_a_file)
        sample_list_b = get_list_from_file(options.sample_list_b_file)
    elif options.sample_list_a_file and not options.sample_list_b_file:
        print_log(f"ERROR: Female parents list file missing.. Quiting...")
        exit(1)
    elif not options.sample_list_a_file and options.sample_list_b_file:
        print_log(f"ERROR: Male parents list file missing.. Quiting...")
        exit(1)
    else:
        sample_list_a = []
        sample_list_b = []

    if options.task_bestmarkers:
        all_summary, ind_summary = get_best_markers(
            smdata, msdata, marker_summary, sample_list_a, sample_list_b)
        all_summary.to_csv(
            out_prefix + "_BestMarkerSummaryAll.txt", sep="\t", index=False)
        ind_summary.to_csv(
            out_prefix + "_BestMarkerSummaryInd.txt", sep="\t", index=False)

    if options.task_find_differences:
        print_log(f"Finding polymorphic markers...")
        differences = find_differences(smdata, sample_list_a, sample_list_b)
        print_log(f"Writing polymorphic markers...")
        differences.to_csv(out_prefix + "_differences.txt",
                           sep="\t", index=False)

    if options.task_checkperformance:
        print_log("Checking performance..")
        performance = check_performance(smdata, sample_list_a, sample_list_b)
        print(f"------------------------------------------------------------")
        print_log(f"Marker Performance on the Data")
        print(f"------------------------------------------------------------")
        print(f"-----------------------------------------------------------------")
        print(
            f"Total combinations                        | {performance[0]}\t\t\t|")
        print(
            f"Combinations with ZERO polymorphic markers| {performance[1]} ({round((performance[1]/performance[0])*100,2)} %)\t\t|")
        print(
            f"Combinations with >= 1 polymorphic markers| {performance[2]} ({round((performance[2]/performance[0])*100,2)} %)\t|")
        print(
            f"Combinations with >= 2 polymorphic markers| {performance[3]} ({round((performance[3]/performance[0])*100,2)} %)\t|")
        print(f"-----------------------------------------------------------------")

    if options.task_write_grid:
        print_log("Writing Grid file..")
        write_grid_file(out_prefix + "_out_grid.csv", smdata, markers)

    if options.task_write_fjk:
        print_log("Writing Flapjack file..")
        write_flapjack_file(out_prefix + "_out_FJ.data", smdata, markers)

    if options.task_write_hapmap:
        print_log("Writing Hapmap file..")
        write_hapmap_file(out_prefix + "_out.hmp.txt", smdata, markers)

    if options.task_make_plots:
        print_log("Writing SNP plots..")
        make_snp_plots(gt_data.msdata, markers, out_prefix, gt_data.name)

    # filename = outdir + "/BestMarkerSummaryInd.txt"
    # ind_summary = pd.read_csv(filename, sep="\t", index_col='marker_count')
    # count = 15
    # print(list(ind_summary.loc[count]['marker']))
    # extract_markers = list(ind_summary.loc[count]['marker'])
    # extract_markers = extract_dict(markers, extract_markers)

    # tmp_smdata = subset_gtdata(smdata, extract_markers, samples, 'samplefast')
    # tmp_msdata = subset_gtdata(msdata, extract_markers, samples, 'markerfast')

    # filename = outdir + "/BestMarkerSummaryInd.txt"
    # ind_summary = pd.read_csv(filename, sep="\t", index_col='marker_count')
    # count = 15
    # print(list(ind_summary.loc[count]['marker']))
    # extract_markers = list(ind_summary.loc[count]['marker'])
    # extract_markers = extract_dict(markers, extract_markers)
    return


if __name__ == '__main__':
    main()
