from utils.utils import *
from classes.Data import *
from collections import defaultdict


def get_consensus(calls):
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
        if counts[base] / (len(calls) - miss_count) * 100 > CONSENSUSCUTOFF:
            return base
        else:
            if counts[base] / (len(calls) - miss_count) * 100 == CONSENSUSCUTOFF:
                nextbase = next_key(list(counts.keys()), base)
                if counts[base] == counts[nextbase]:
                    return 'N:N'
                else:
                    return base
            else:
                return 'N:N'


def get_consensus_dict(smdata, msdata, reps_dict):
    cons_smdata = defaultdict()
    cons_msdata = defaultdict()
    for cons_name in reps_dict:
        reps = reps_dict[cons_name]
        for marker in msdata:
            calls = msdata[marker].get_data(reps)
            if len(calls) >= 1:
                call = get_consensus(calls)
                if not marker in cons_msdata:
                    cons_msdata[marker] = MS(marker)
                if not cons_name in cons_smdata:
                    cons_smdata[cons_name] = SM(cons_name)
                cons_smdata[cons_name].put_data(marker, call)
                cons_msdata[marker].put_data(cons_name, call)
    return cons_smdata, cons_msdata
