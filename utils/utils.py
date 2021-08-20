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
    pic = 1 - ((freq_a * freq_a) + (freq_b * freq_b)) - (2 * (freq_a * freq_a) * (freq_b * freq_b))
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
        print(f"Error: {data_type} is unknown. Should be 'samplefast' or 'markerfast'..")
        exit(1)
    return out_data
        
def fill_gaps_gtdata(smdata, msdata):
    markers = []
    samples = []
    missing_call = Call('N:N')
    markers = list(msdata.keys())
    samples = list(smdata.keys())
    
    if not markers:
        print(f"ERROR: Input markerfast dictionaty is empty.")
        exit(1)
    
    if not samples:
        print(f"ERROR: Input samplefast dictionaty is empty.")
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

def get_dup_keys(in_dict, q_value):
    return [k for k,v in in_dict.items() if v == q_value]

def extract_dict(in_dict, keys_list):
    out_dict = defaultdict()
    for key in keys_list:
        if key in in_dict:
            out_dict[key] = in_dict[key]
    
    return out_dict