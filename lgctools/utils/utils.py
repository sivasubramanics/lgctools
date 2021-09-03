import argparse
import inspect
import time
from collections import defaultdict
from utils.definitions import *
from classes.Data import *


def is_homo(call):
    """
    Check if the genotype call is homozygous call eg: G:G
    """
    calls = call.split(':')
    if call in MISSING_CALLS:
        return False
    elif calls[0] == calls[1]:
        return True
    else:
        return False


def is_hetero(call):
    """
    Check if the genotype call is heterozygous call eg: G:T
    """
    calls = call.split(':')
    if call in MISSING_CALLS:
        return False
    elif calls[0] != calls[1]:
        return True
    else:
        return False


def is_polymorphic(call_a, call_b):
    """
    Check if calls are polymorphic
    """
    if is_homo(call_a) and is_homo(call_b):
        if call_a != call_b:
            return True
    else:
        return False


def unique(in_list):
    """
    Returns list of unique elements
    """
    out_set = set(in_list)
    return list(out_set)


def intersect(list_a, list_b):
    return list(set(list_a) & set(list_b))


def remove_elements(in_list, rm_list):
    out_list = []
    for elem in in_list:
        if not elem in rm_list:
            out_list.append(elem)
    return out_list


def get_counts(in_list):
    """
    Returns dictionary of count for each element from a list
    """
    counts = defaultdict()
    calls = unique(in_list)
    for call in calls:
        counts[call] = in_list.count(call)
    return counts


def initialize_list(count_element):
    """
    Returns list with number of elements initialized to have value ZERO
    """
    out_list = []
    for i in range(count_element):
        out_list.append(0)
    return out_list


def get_pic(freq_a, freq_b):
    pic = 1 - ((freq_a * freq_a) + (freq_b * freq_b)) - \
        (2 * (freq_a * freq_a) * (freq_b * freq_b))
    return round(pic, 2)


def allele_freq(in_dict, seperator):
    """
    Returns list of alleles and its count from a dictionary of allele counts.
    [major allele, major allele count, minor allele, minor allele count, hetero count]
    """
    al_count = defaultdict()
    het_count = 0
    for call in in_dict:
        if not call in MISSING_CALLS:
            alleles = call.split(seperator)
            if not alleles[0] in al_count:
                al_count[alleles[0]] = 0
            if not alleles[1] in al_count:
                al_count[alleles[1]] = 0
            al_count[alleles[0]] += 1
            al_count[alleles[1]] += 1
            if alleles[0] != alleles[1]:
                het_count += 1
    al_count = sort_dict(al_count, True)
    alleles = list(al_count.keys())
    if not alleles:
        return ['na', 0, 'na', 0, 0]
    if len(alleles) == 1:
        return [alleles[0], al_count[alleles[0]], 'na', 0, het_count]
    else:
        return [alleles[0], al_count[alleles[0]], alleles[1], al_count[alleles[1]], het_count]


def to_flapjack(call):
    """
    Converts LGC format genotype call to Flapjack format by replacing ':' to '/'
    """
    if call in MISSING_CALLS:
        return 'N/N'
    else:
        calls = call.split(":")
        return calls[0] + '/' + calls[1]


def to_hmp(call):
    """
    Converts LGC format genotype call to Hapmap format by replacing ':' to ''
    """
    if call in MISSING_CALLS:
        return 'NN'
    else:
        calls = call.split(":")
        return calls[0] + calls[1]


def to_grid(call):
    """
    Converts Hapmap format genotype call to grid format by replacing '' to ':'
    """
    if call in MISSING_CALLS:
        return 'N:N'
    else:
        calls = list(call)
        return calls[0] + ':' + calls[1]


def num(s):
    """
    Returns numeric value of the input
    """
    try:
        return int(s)
    except ValueError:
        return float(s)


def average(lst):
    """
    Returns average (mean) of the list of numbers
    """
    return sum(lst) / len(lst)


def sort_dict(in_dict, reverse):
    """
    Returns sorted dictionary based on value
    """
    if reverse:
        return dict(sorted(in_dict.items(), key=lambda item: item[1], reverse=True))
    else:
        return dict(sorted(in_dict.items(), key=lambda item: item[1], reverse=False))


def subset_gtdata(in_data, markers, samples, data_type):
    if data_type == 'samplefast':
        out_data = defaultdict()
        for sample in in_data:
            if not sample in samples:
                continue
            out_data[sample] = SM(sample)
            for marker in in_data[sample].data:
                if not marker in markers:
                    continue
                out_data[sample].put_data(marker, in_data[sample].data[marker])
    elif data_type == 'markerfast':
        out_data = defaultdict()
        for marker in in_data:
            if not marker in markers:
                continue
            out_data[marker] = MS(marker)
            for sample in in_data[marker].data:
                if not sample in samples:
                    continue
                out_data[marker].put_data(sample, in_data[marker].data[sample])
    else:
        print(
            f"Error: {data_type} is unknown. Should be 'samplefast' or 'markerfast'..")
        exit(1)
    return out_data


def fill_gaps_gtdata(smdata, msdata):
    markers = []
    samples = []
    missing_call = Call('N:N')
    markers = list(msdata.keys())
    samples = list(smdata.keys())

    if not markers:
        print(f"ERROR: Input markerfast dictionary is empty.")
        exit(1)

    if not samples:
        print(f"ERROR: Input samplefast dictionary is empty.")
        exit(1)

    for marker in msdata:
        for sample in samples:
            if not sample in msdata[marker].data:
                msdata[marker].data[sample] = missing_call

    for sample in smdata:
        for marker in markers:
            if not marker in smdata[sample].data:
                smdata[sample].data[marker] = missing_call

    return smdata, msdata


def merge_gtdata(*args):
    out_dict = defaultdict()
    for i in range(0, len(args)):
        out_dict.update(args[i])
    return out_dict


def get_dup_keys(in_dict, q_value):
    return [k for k, v in in_dict.items() if v == q_value]


def extract_dict(in_dict, keys_list):
    out_dict = defaultdict()
    for key in keys_list:
        if key in in_dict:
            out_dict[key] = in_dict[key]
    return out_dict


def filter_markers(markers, flt_marker_summary):
    flt_markers = list(flt_marker_summary['marker_name'])
    new_markers = defaultdict(Markermetadata)
    for marker in markers:
        if marker in flt_markers:
            new_markers[marker] = markers[marker]
    markers = new_markers.copy()
    return markers


def filter_samples(samples, flt_sample_summary):
    flt_samples = list(flt_sample_summary['sample_name'])
    new_samples = []
    for sample in samples:
        if sample in flt_samples:
            new_samples.append(sample)
    samples = new_samples.copy()
    return samples


def get_opts():
    parser = argparse.ArgumentParser(
        description="QC pipeline: Processes the LGC file for tasks involved in purity check")
    parser.add_argument('--summary', dest='task_run_summary', default=False,
                        action='store_true', help="Generate sample and marker summary for the data from genotype file")
    parser.add_argument('--rename', dest='task_rename', default=False,
                        action='store_true', help="Rename the samples provided in the genotype file")
    parser.add_argument('--out-grid', dest='task_write_grid', default=False,
                        action='store_true', help="Write genotype data in grid file format")
    parser.add_argument('--out-flapjack', dest='task_write_fjk', default=False,
                        action='store_true', help="Write genotype data in flapjack file format")
    parser.add_argument('--out-hapmap', dest='task_write_hapmap', default=False,
                        action='store_true', help="Write genotype data in hapmap file format")
    parser.add_argument('--differences', dest='task_find_differences', default=False,
                        action='store_true', help="Find polymorphic markers between genotypes")
    parser.add_argument('--performance', dest='task_checkperformance', default=False,
                        action='store_true', help="Check performance of the provided marker set")
    parser.add_argument('--make-plots', dest='task_make_plots', default=False,
                        action='store_true', help="Create plots based on the LGC data")
    parser.add_argument('--bestmarkers', dest='task_bestmarkers', default=False,
                        action='store_true', help="Find best markers from the given marker set")
    parser.add_argument('--markercloud', dest='task_markercloud', default=False,
                        action='store_true', help="Based on best marker summary make word cloud")
    parser.add_argument('--filter', dest='task_filter', default=False,
                        action='store_true', help="Filter genotype data")
    parser.add_argument('--pedver', dest='task_pedver', default=False,
                        action='store_true', help="Analyze genotype data for F1 Verification")
    parser.add_argument('--consensus', dest='task_consensus', default=False,
                        action='store_true', help="Call consensus and make purity reports")
    parser.add_argument('--fwdbreed', dest='task_pedver', default=False,
                        action='store_true', help="Analyze genotype data for favorable allele propotions")
    parser.add_argument("--lgc-file", dest="lgc_file",
                        metavar="<FILE>", help="LGC raw data File")
    parser.add_argument("--lgc-files", nargs='+', dest="lgc_files",
                        metavar="<FILES>", help="LGC raw data File")
    parser.add_argument("--grid-file", dest="grid_file",
                        metavar="<FILE>", help="LGC Grid Matrix File")
    parser.add_argument("--grid-files", dest="grid_files",
                        metavar="<FILES>", help="Comma seperated LGC Grid Matrix Files")
    parser.add_argument("--hmp-file", dest="hapmap_file",
                        metavar="<FILE>", help="Hapmap genotype File")
    parser.add_argument("--samplemap-file", dest="samplemap_file",
                        metavar="<FILE>", help="Tab seperated sample map file. <SAMPLE_ID> <SAMPLE_NAME>")
    parser.add_argument("--ped-file", dest="pedigree_file",
                        metavar="<FILE>", help="Tab seperated Pedigree File. <F_ONE> <PAR_A> <PAR_B>")
    parser.add_argument("--designation-file", dest="designation_file",
                        metavar="<FILE>", help="Tab seperated Designations parent File. <SAMPLE_NAME> <DESIGNATION>")
    parser.add_argument("--meta-data", dest="metadata_file",
                        metavar="<FILE>", help="Parent information file")
    parser.add_argument("--marker-list", dest="marker_list_file",
                        metavar="<FILE>", help="File with list of snps to analyze")
    parser.add_argument("--sample-list", dest="sample_list_file",
                        metavar="<FILE>", help="File with list of samples to analyze")
    parser.add_argument("--male-parents-list", dest="sample_list_a_file",
                        metavar="<FILE>", help="File with list of male parents to consider")
    parser.add_argument("--female-parents-list", dest="sample_list_b_file",
                        metavar="<FILE>", help="File with list of female parents to consider")
    parser.add_argument("--qtl-file", dest="qtl_file",
                        metavar="<FILE>", help="QTL file in GOBii format")
    parser.add_argument("--out", dest="out_prefix", default="out",
                        metavar="<STR>", help="Output filename prefix")
    parser.add_argument("--f1het-cutOff", dest="cutoff_fone", default=CUTOFF_FONE, type=float,
                        metavar="<INT>", help="Percentage expected heterozygosity for F1 verification")
    parser.add_argument("--consensus-cutOff", dest="cutoff_consensus", default=CUTOFF_CONSENSUS, type=float,
                        metavar="<INT>", help="Percentage propotion to be considered to call consensus")
    parser.add_argument("--max-missing-site", dest="cutoff_mx_missing_marker", default=MARKER_MISSING_MAXIMUM, type=float,
                        metavar="<INT>", help="Percentage propotion to be considered for maximum missing per marker")
    parser.add_argument("--min-pic-site", dest="cutoff_mn_pic_marker", default=MARKER_PIC_MINIMUM, type=float,
                        metavar="<FLOAT>", help="Percentage propotion to be considered for PIC per marker")
    parser.add_argument("--min-maf-site", dest="cutoff_mn_maf_marker", default=MARKER_MAF_MINIMUM, type=float,
                        metavar="<FLOAT>", help="Percentage propotion to be considered for MAF per marker")
    parser.add_argument("--max-missing-sample", dest="cutoff_mx_missing_sample", default=SAMPLE_MISSING_MAXIMUM, type=float,
                        metavar="<INT>", help="Percentage propotion to be considered for maximum missing per sample")
    # Parse commandline arguments
    return parser


def get_list_from_file(in_file):
    func_name = f"{{method: {inspect.stack()[0][3]}}}"
    out_list = []
    with open(in_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.__contains__(TAB) or line.__contains__(CSV):
                print(
                    f"{func_name} ERROR: Input file {in_file} contains more than one column.")
                exit(1)
            out_list.append(line)
    return out_list


def print_log(in_str):
    print(f"[{time.asctime()}] - {in_str}")


def secondsToText(secs):
    days = secs//86400
    hours = (secs - days*86400)//3600
    minutes = (secs - days*86400 - hours*3600)//60
    seconds = secs - days*86400 - hours*3600 - minutes*60
    result = ("{0} day{1}, ".format(days, "s" if days != 1 else "") if days else "") + \
        ("{0} hour{1}, ".format(hours, "s" if hours != 1 else "") if hours else "") + \
        ("{0} minute{1}, ".format(minutes, "s" if minutes != 1 else "") if minutes else "") + \
        ("{0} second{1}".format(
            seconds, "s" if seconds != 1 else "") if seconds else "")
    return result


def next_key(tmpList, current_key):
    # temp = list(test_dict)
    try:
        res = tmpList[tmpList.index(current_key) + 1]
        if res == 'N':
            next_key(tmpList, res)
    except (ValueError, IndexError):
        res = None
    return res


def get_list_from_dict(in_dict):
    out_list = []
    for key in in_dict:
        if type(in_dict[key]) is list:
            out_list += in_dict[key]
        else:
            out_list.append(in_dict[key])
    return unique(out_list)


def compare_two_samples(data_cons, data_rep):
    markers = intersect(list(data_cons.data.keys()),
                        list(data_rep.data.keys()))
    match = 0
    mismatch = 0
    missing = 0
    for marker in markers:
        call_c = data_cons.data[marker].__str__()
        call_r = data_rep.data[marker].__str__()
        if call_c in MISSING_CALLS and call_r in MISSING_CALLS:
            missing += 1
        elif call_c == call_r:
            match += 1
        else:
            if call_r in MISSING_CALLS:
                missing += 1
            else:
                mismatch += 1
    return [match, mismatch, missing, round(match/(match+mismatch)*100, 2)]
