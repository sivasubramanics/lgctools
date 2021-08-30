from collections import defaultdict, OrderedDict
from utils.utils import *


def check_performance(smdata, sample_list_a, sample_list_b, *args):
    """
    Returns number of combinations and marker performance as list
    [total combinations, combinations with zero markers, combinations with one marker, combinations with two markers]
    """
    if not sample_list_a and not sample_list_b:
        sample_list = list(smdata.keys())
        no_polymorphic = []
        for a in range(0, len(sample_list)):
            sample_a = sample_list[a]
            for b in range(a+1, len(sample_list)):
                sample_b = sample_list[b]
                if sample_a != sample_b:
                    no_polymorphic.append(check_polymorphic(list(map(
                        str, smdata[sample_a].data.values())), list(map(str, smdata[sample_b].data.values()))))
        counts = dict()
    else:
        # sample_list = list(smdata.keys())
        no_polymorphic = []
        for a in range(0, len(sample_list_a)):
            sample_a = sample_list_a[a]
            for b in range(0, len(sample_list_b)):
                sample_b = sample_list_b[b]
                if sample_a != sample_b:
                    no_polymorphic.append(check_polymorphic(list(map(
                        str, smdata[sample_a].data.values())), list(map(str, smdata[sample_b].data.values()))))
        counts = dict()
    for i in no_polymorphic:
        counts[i] = counts.get(i, 0) + 1
    counts = OrderedDict(sorted(counts.items()))
    total = len(no_polymorphic)
    if not args:
        return [total, counts[0], total-counts[0], total-counts[0]-counts[1]]
    else:
        marker_scores_dict = args[0]
        marker_name = args[1]
        if not marker_name in marker_scores_dict:
            marker_scores_dict[marker_name] = [
                total, counts[0], total-counts[0], total-counts[0]-counts[1]]


def check_polymorphic(list_a, list_b):
    """
    Returns number of polymorphic markers between 2 lists containing genotype calls
    More then 2 it doesn't iterate to save processing time
    """
    len_list_a = len(list_a)
    len_list_b = len(list_b)
    no_polymorphic = 0
    if len_list_a == len_list_b:
        for i in range(0, len_list_a):
            if is_polymorphic(list_a[i], list_b[i]):
                no_polymorphic += 1
            if no_polymorphic == 2:
                return no_polymorphic
    return no_polymorphic
