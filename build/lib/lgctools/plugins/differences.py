
from ..utils.definitions import CSV, SPACE
import pandas as pd
from ..utils.utils import is_polymorphic


def polymorphic_calls(smdata, sample_a, sample_b):
    call_polymorphic = ['', '', 0, 'na']
    call_polymorphic[0] = sample_a
    call_polymorphic[1] = sample_b
    polymorphic = []
    markers_a = smdata[sample_a].data.keys()
    markers_b = smdata[sample_b].data.keys()
    if markers_a == markers_b:
        for marker in markers_a:
            call_a = smdata[sample_a].data[marker].__str__()
            call_b = smdata[sample_b].data[marker].__str__()
            if is_polymorphic(call_a, call_b):
                polymorphic.append(
                    marker + ' [' + call_a + '|' + call_b + ']')
                call_polymorphic[2] += 1
    if polymorphic:
        call_polymorphic[3] = CSV.join(polymorphic)
    return call_polymorphic


def find_differences(smdata, sample_list_a=[], sample_list_b=[], *args):
    differences = pd.DataFrame(columns=['Genotype_A',
                                        'Genotype_B',
                                        'Count_Polymorphic',
                                        'Markers'])
    total_count = 0
    current_count = 0
    if not sample_list_a and not sample_list_b:
        sample_list = list(smdata.keys())
        total_count = len(sample_list) * (len(sample_list) - 1)/2
        for a in range(0, len(sample_list)):
            sample_a = sample_list[a]
            for b in range(a+1, len(sample_list)):
                current_count += 1
                print(f"\rStatus: {round(current_count/total_count * 100,2)} %",
                      end='   ', flush=True)
                sample_b = sample_list[b]
                if sample_a != sample_b:
                    differences.loc[len(differences.index)] = polymorphic_calls(
                        smdata, sample_a, sample_b)
    else:
        total_count = len(sample_list_a) * len(sample_list_a)
        for a in range(0, len(sample_list_a)):
            sample_a = sample_list_a[a]
            for b in range(0, len(sample_list_b)):
                current_count += 1
                print(f"\rStatus: {round(current_count/total_count * 100,2)} %",
                      end='   ', flush=True)
                sample_b = sample_list_b[b]
                if sample_a != sample_b:
                    if sample_a != sample_b:
                        differences.loc[len(differences.index)] = polymorphic_calls(
                            smdata, sample_a, sample_b)
    print()
    print("\033[A\033[A")
    return differences
