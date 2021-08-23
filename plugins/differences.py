
import pandas as pd
from utils.definitions import TAB
from utils.utils import is_polymorphic


def find_differences(smdata, *args):
    differences = pd.DataFrame(columns = ['Genotype_A', 
                                          'Genotype_B',
                                          'Count_Polymorphic',
                                          'Markers'])
    sample_list = list(smdata.keys())
    for a in range(0, len(sample_list)):
        sample_a = sample_list[a]
        for b in range(a+1, len(sample_list)):
            sample_b = sample_list[b]
            if sample_a != sample_b:
                call_polymorphic = ['', '', 0, []]
                call_polymorphic[0] = sample_a
                call_polymorphic[1] = sample_b
                markers_a = smdata[sample_a].data.keys()
                markers_b = smdata[sample_b].data.keys()
                if markers_a == markers_b:
                    for marker in markers_a:
                        call_a = smdata[sample_a].data[marker].__str__()
                        call_b = smdata[sample_b].data[marker].__str__()
                        if is_polymorphic(call_a, call_b):
                            call_polymorphic[3].append(marker + ' ' + call_a + '|' + call_b)
                            call_polymorphic[2] += 1
                differences.loc[len(differences.index)] = call_polymorphic
    return differences
