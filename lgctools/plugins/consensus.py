from utils.file_processing import get_replicates, write_grid_file
from utils.utils import *
from classes.Data import *
from collections import defaultdict
import pandas as pd


def get_consensus(calls, cutoff):
    calls = list(map(str, calls))
    counts = get_counts(calls)
    counts = sort_dict(counts, True)
    miss_count = 0
    for miss_call in MISSING_CALLS:
        if miss_call in counts:
            miss_count += counts[miss_call]
    for base in counts:
        if base in MISSING_CALLS:
            continue
        if counts[base] / (len(calls) - miss_count) * 100 > cutoff:
            return base
        else:
            if counts[base] / (len(calls) - miss_count) * 100 == cutoff:
                nextbase = next_key(list(counts.keys()), base)
                if counts[base] == counts[nextbase]:
                    return 'N:N'
                else:
                    return base
            else:
                return 'N:N'


def get_consensus_dict(smdata, msdata, reps_dict, cutoff):
    cons_smdata = defaultdict()
    cons_msdata = defaultdict()
    for cons_name in reps_dict:
        reps = reps_dict[cons_name]
        for marker in msdata:
            calls = msdata[marker].get_data(reps)
            if len(calls) >= 1:
                call = get_consensus(calls, cutoff)
                if not marker in cons_msdata:
                    cons_msdata[marker] = MS(marker)
                if not cons_name in cons_smdata:
                    cons_smdata[cons_name] = SM(cons_name)
                cons_smdata[cons_name].put_data(marker, call)
                cons_msdata[marker].put_data(cons_name, call)
    return cons_smdata, cons_msdata


def get_consensus_summary(cons_rep_smdata, replicates):
    consensus_summary = pd.DataFrame(
        columns=['Rep_Name', 'Match', 'MisMatch', 'Missing', 'PurityScore'])
    for cons_name in replicates:
        if cons_name in cons_rep_smdata:
            cons_calls = cons_rep_smdata[cons_name]
            for rep in replicates[cons_name]:
                rep_calls = cons_rep_smdata[rep]
                counts = compare_two_samples(
                    cons_calls, rep_calls)
                counts.insert(0, rep)
                consensus_summary.loc[len(consensus_summary.index)] = counts
    consensus_summary = consensus_summary.set_index('Rep_Name')
    return consensus_summary


def get_consensus_report(consensus_summary, replicates):
    consensus_report = pd.DataFrame(
        columns=['Name', 'NoReps', 'MeanScore'])
    for cons_name in replicates:
        score = 0
        rep_count = 0
        for rep in replicates[cons_name]:
            if rep in consensus_summary.index:
                rep_count += 1
                score += consensus_summary.at[rep, 'PurityScore']
        if rep_count > 0:
            consensus_report.loc[len(consensus_report.index)] = [
                cons_name, rep_count, round(score/rep_count, 2)]
    return consensus_report


def process_consensus(options, smdata, msdata, markers):
    if options.designation_file:
        designation_file = options.designation_file
    elif options.pedigree_file:
        designation_file = options.pedigree_file
    else:
        print_log(
            f"ERROR: Designation file/Pedigree file is needed to call consensus... Quiting...")
        exit(1)
    out_prefix = options.out_prefix
    # Dictionary with deignation and its replicates
    replicates = get_replicates(designation_file)
    # List of all replicates
    replicates_list = get_list_from_dict(replicates)
    cons_smdata, cons_msdata = get_consensus_dict(
        smdata, msdata, replicates, options.cutoff_consensus)
    rep_msdata = subset_gtdata(
        msdata, markers, replicates_list, 'markerfast')
    rep_smdata = subset_gtdata(
        smdata, markers, replicates_list, 'samplefast')
    cons_rep_msdata = merge_gtdata(cons_msdata, rep_msdata)
    cons_rep_smdata = merge_gtdata(cons_smdata, rep_smdata)
    cons_rep_smdata, cons_rep_msdata = fill_gaps_gtdata(
        cons_rep_smdata, cons_rep_msdata)
    print_log(
        f"Writing consensus calls and replicate calls ({out_prefix + '_consensus_grid.csv'})...")
    write_grid_file(out_prefix + '_consensus_grid.csv',
                    cons_rep_smdata, markers)
    consensus_summary = get_consensus_summary(cons_rep_smdata, replicates)
    print_log(
        f"Writing consensus summary ({out_prefix + '_PuritySummary.txt'})...")
    consensus_summary.to_csv(out_prefix + "_PuritySummary.txt",
                             sep="\t")
    consensus_report = get_consensus_report(consensus_summary, replicates)
    print_log(
        f"Writing consensus report ({out_prefix + '_Purityreport.txt'})...")
    consensus_report.to_csv(out_prefix + "_PurityReport.txt",
                            sep="\t", index=False)

    msdata = merge_gtdata(cons_msdata, msdata)
    smdata = merge_gtdata(cons_smdata, smdata)
    samples = list(smdata.keys())
    samples = remove_elements(samples, replicates_list)
    msdata = subset_gtdata(
        msdata, markers, samples, 'markerfast')
    smdata = subset_gtdata(
        smdata, markers, samples, 'samplefast')
    return msdata, smdata
