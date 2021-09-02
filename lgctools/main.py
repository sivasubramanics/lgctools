#!/usr/bin/env python

from plugins.consensus import *
from plugins.bestmarkers import *
from plugins.performance import *
from plugins.differences import *
from plugins.rename import *
from plugins.plots import *
from plugins.summary import *
from classes.Data import Markermetadata
from classes.LGC import LGC
from classes.Grid import Grid
from classes.Hapmap import Hapmap
from utils.file_processing import *
from utils.definitions import *
from collections import defaultdict
import sys
import time
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


def main():

    start_time = time.time()

    parser = get_opts()
    options = parser.parse_args()

    if len(sys.argv) <= 1:
        parser.print_help()
        exit(1)

    marker_info = get_marker_info()
    markers = defaultdict(Markermetadata)
    out_prefix = options.out_prefix
    smdata, msdata = defaultdict(), defaultdict()

    if options.pedigree_file:
        print_log(
            f"Processing complete for pedigree file ({options.pedigree_file})...")
        pedigree_dict = get_pedigree(options.pedigree_file)

    if options.lgc_file:
        print_log(f"Reading input LGC genotype data ({options.lgc_file})...")
        gt_data = LGC(options.lgc_file)
        msdata = gt_data.msdata
        smdata = gt_data.smdata
        markers = get_markers(options.lgc_file, markers, marker_info)
        samples = list(smdata.keys())

    if options.grid_file:
        print_log(f"Reading input Grid genotype data ({options.grid_file})...")
        gt_data = Grid(options.grid_file)
        msdata = gt_data.msdata
        smdata = gt_data.smdata
        marker = make_markers(msdata, markers, marker_info)
        samples = list(smdata.keys())

    if options.hapmap_file:
        print_log(
            f"Reading input Hapmap genotype data ({options.hapmap_file})...")
        gt_data = Hapmap(options.hapmap_file)
        msdata = gt_data.msdata
        smdata = gt_data.smdata
        marker = make_markers(msdata, markers, marker_info)
        samples = list(smdata.keys())

    if not msdata and not smdata:
        print_log(
            f"ERROR: No genotype provided... Quiting...")
        exit(1)

    smdata, msdata = fill_gaps_gtdata(smdata, msdata)
    print_log(
        f"Input data : {len(smdata)} samples x {len(msdata)} markers")

    if options.task_rename:
        if options.samplemap_file:
            print_log(
                f"Processing complete for sample map file ({options.samplemap_file})...")
            samplemap_dict = get_samplemap(options.samplemap_file)
            smdata, msdata = rename_data(smdata, msdata, samplemap_dict)
        else:
            print_log(
                f"ERROR: No sample map file provided (--samplemap-file)... Quiting...")
            exit(1)

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
        print_log(
            f"Subset data : {len(smdata)} samples x {len(msdata)} markers")

    if options.task_run_summary:
        marker_summary = get_marker_summary(msdata)
        print_log(
            f"Writing Marker Summary ({out_prefix + '_marker_summary.txt'})...")
        marker_summary.to_csv(out_prefix + "_marker_summary.txt",
                              sep="\t", index=False)

        sample_summary = get_sample_summary(smdata)
        print_log(
            f"Writing Sample Summary ({out_prefix + '_sample_summary.txt'})..")
        sample_summary.to_csv(out_prefix + "_sample_summary.txt",
                              sep="\t", index=False)

    if options.task_filter:
        if not marker_summary:
            marker_summary = get_marker_summary(msdata)
        if not sample_summary:
            sample_summary = get_sample_summary(smdata)

        print_log(f"Filtering data...")
        flt_sample_summary = sample_summary.loc[(sample_summary['missing_percentage']
                                                < options.cutoff_mx_missing_sample)]
        samples = filter_samples(samples, flt_sample_summary)
        flt_marker_summary = marker_summary.loc[(marker_summary['missing_percentage'] <= options.cutoff_mx_missing_marker)
                                                & (marker_summary['minor_allale_freq'] >= options.cutoff_mn_maf_marker)
                                                & (marker_summary['PIC'] >= options.cutoff_mn_pic_marker)]
        markers = filter_markers(markers, flt_marker_summary)
        smdata = subset_gtdata(smdata, markers, samples, 'samplefast')
        msdata = subset_gtdata(msdata, markers, samples, 'markerfast')
        print_log(
            f"Filtered data : {len(smdata)} samples x {len(msdata)} markers")
        if options.task_run_summary:
            print_log(
                f"Writing Filtered Sample Summary ({out_prefix + '_sample_summary_flt.txt'})...")
            flt_sample_summary.to_csv(
                out_prefix + "_sample_summary_flt.txt", sep="\t", index=False)
            print_log(
                f"Writing Filtered Marker Summary ({out_prefix + '_marker_summary_flt.txt'})...")
            flt_marker_summary.to_csv(
                out_prefix + "_marker_summary_flt.txt", sep="\t", index=False)

    if options.sample_list_a_file and options.sample_list_b_file:
        sample_list_a = get_list_from_file(options.sample_list_a_file)
        sample_list_b = get_list_from_file(options.sample_list_b_file)
    elif options.sample_list_a_file and not options.sample_list_b_file:
        print_log(f"ERROR: Female parents list file missing... Quiting...")
        exit(1)
    elif not options.sample_list_a_file and options.sample_list_b_file:
        print_log(f"ERROR: Male parents list file missing... Quiting...")
        exit(1)
    else:
        sample_list_a = []
        sample_list_b = []

    if options.task_bestmarkers:
        if not marker_summary:
            marker_summary = get_marker_summary(msdata)
        if not sample_summary:
            sample_summary = get_sample_summary(smdata)
        print_log(f"Finding best markers...")
        all_summary, ind_summary = get_best_markers(
            smdata, msdata, marker_summary, sample_list_a, sample_list_b)
        all_summary.to_csv(
            out_prefix + "_BestMarkerSummaryAll.txt", sep="\t", index=False)
        ind_summary.to_csv(
            out_prefix + "_BestMarkerSummaryInd.txt", sep="\t", index=False)
        if options.task_markercloud:
            print_log(f"Making Best Markers cloud Image...")
            make_wordcloud(ind_summary, out_prefix + '_BestMarkerCloud.png')

    if options.task_consensus:
        if not options.designation_file:
            print_log(
                f"ERROR: Designation file is needed to call consensus... Quiting...")
            exit(1)
        msdata, smdata = process_consensus(
            options.designation_file, smdata, msdata, markers, out_prefix)

    if options.task_pedver:
        if not options.pedigree_file:
            print_log(
                f"ERROR: Pedigree file is needed to call consensus and perform pedver... Quiting...")
            exit(1)
        msdata, smdata = process_consensus(
            options.pedigree_file, smdata, msdata, markers, out_prefix)

    if options.task_find_differences:
        print_log(f"Finding polymorphic markers...")
        differences = find_differences(smdata, sample_list_a, sample_list_b)
        print_log(
            f"Writing polymorphic markers ({out_prefix + '_differences.txt'})...")
        differences.to_csv(out_prefix + "_differences.txt",
                           sep="\t", index=False)

    if options.task_checkperformance:
        print_log("Checking performance...")
        performance = check_performance(smdata, sample_list_a, sample_list_b)
        print(f"------------------------------------------------------------")
        print_log(f"Marker Performance on the Data")
        print(f"------------------------------------------------------------")
        print(f"------------------------------------------------------------")
        print(
            f"Total combinations                        : {performance[0]}")
        print(
            f"Combinations with ZERO polymorphic markers: {performance[1]} ({round((performance[1]/performance[0])*100,2)} %)")
        print(
            f"Combinations with >= 1 polymorphic markers: {performance[2]} ({round((performance[2]/performance[0])*100,2)} %)")
        print(
            f"Combinations with >= 2 polymorphic markers: {performance[3]} ({round((performance[3]/performance[0])*100,2)} %)")
        print(f"------------------------------------------------------------")

    if options.task_write_grid:
        print_log(f"Writing Grid file ({out_prefix + '_out_grid.csv'})...")
        write_grid_file(out_prefix + "_out_grid.csv", smdata, markers)

    if options.task_write_fjk:
        print_log(f"Writing Flapjack file ({out_prefix + '_out_FJ.data'})...")
        write_flapjack_file(out_prefix + "_out_FJ.data", smdata, markers)

    if options.task_write_hapmap:
        print_log(f"Writing Hapmap file ({out_prefix + '_out.hmp.txt'})...")
        write_hapmap_file(out_prefix + "_out.hmp.txt", smdata, markers)

    if options.task_make_plots:
        print_log("Writing SNP plots...")
        make_snp_plots(gt_data.msdata, markers, out_prefix, gt_data.name)

    end_time = time.time()
    process_time = round(end_time - start_time, 2)
    print_log(
        f"Total Time taken for the process {secondsToText(process_time)}")
    print()
    return


if __name__ == '__main__':
    main()
