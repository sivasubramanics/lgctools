from collections import defaultdict
from utils.definitions import *

def unique(in_list):
    out_set = set(in_list)
    return list(out_set)

def get_counts(in_list):
    counts = defaultdict()
    calls = unique(in_list)
    for call in calls:
        counts[call] = in_list.count(call)
    return counts
        
def initialize_list(count_element):
    out_list = []
    for i in range(count_element):
        out_list.append(0)
    return out_list

def allele_freq(in_dict, seperator):
    al_count = defaultdict()
    het_count = 0
    for call in in_dict:
        if not call in MISSING_CALLS:
            alleles = call.split(seperator)
            if not alleles[0] in al_count:
                al_count[alleles[0]] = 1
            if not alleles[1] in al_count:
                al_count[alleles[1]] = 1
            al_count[alleles[0]] += 1
            al_count[alleles[1]] += 1
            if alleles[0] != alleles[1]:
                het_count += 1
    al_count = sort_dict(al_count, True)
    alleles = list(al_count.keys())
    if len(alleles) == 1:
        return [alleles[0], al_count[alleles[0]], 'na', 0, het_count]
    else:
        return [alleles[0], al_count[alleles[0]], alleles[1], al_count[alleles[1]], het_count]

def to_flapjack(call):
    if call in MISSING_CALLS:
        return 'N/N'
    else:
        calls = call.split(":")
        return calls[0] + '/' + calls[1]

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def average(lst):
    return sum(lst) / len(lst)


def sort_dict(in_dict, reverse):
    if reverse:
        return dict(sorted(in_dict.items(), key=lambda item: item[1], reverse=True))
    else:
        return dict(sorted(in_dict.items(), key=lambda item: item[1], reverse=False))