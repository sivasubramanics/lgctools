from utils.file_processing import *


def process_pedver(pedigree_file, smdata, msdata, markers, out_prefix):
    pedver_summary = pd.DataFrame(
        columns=['F_One_Name',
                 'Parent_A',
                 'Parent_A_Score',
                 'Parent_B',
                 'Parent_B_Score',
                 'Polymorphic',
                 'Call_Hetero',
                 'Call_ParentA',
                 'Call_ParentB',
                 'Call_Missing',
                 'Percentage_Het',
                 'Comment'])
    pedigree = get_pedigree(pedigree_file)
    consensus_report = pd.read_csv(
        out_prefix + '_PurityReport.txt', sep=TAB, index_col='Name')
    for f_one in pedigree:
        parent_a = pedigree[f_one].get_parent_a()
        parent_b = pedigree[f_one].get_parent_b()
        if f_one in smdata and parent_a in smdata and parent_b in smdata:
            markers_a = list(smdata[parent_a].data.keys())
            markers_b = list(smdata[parent_b].data.keys())
            markers_x = list(smdata[f_one].data.keys())
            count_polymorphic = 0
            count_parent_a = 0
            count_parent_b = 0
            count_het = 0
            count_missing = 0
            verdict = 'failedF1'
            if markers_a == markers_x and markers_b == markers_x:
                for marker in smdata[f_one].data:
                    call_a = smdata[parent_a].data[marker].__str__()
                    call_b = smdata[parent_b].data[marker].__str__()
                    call_x = smdata[f_one].data[marker].__str__()
                    if is_homo(call_a) and is_homo(call_b) and call_a != call_b:
                        count_polymorphic += 1
                        if call_a == call_x:
                            count_parent_a += 1
                        elif call_b == call_x:
                            count_parent_b += 1
                        elif call_x in MISSING_CALLS:
                            count_missing += 1
                        elif call_a.split(":")[0] == call_x.split(":")[0] and call_b.split(":")[0] == call_x.split(":")[1]:
                            count_het += 1
                        elif call_a.split(":")[0] == call_x.split(":")[1] and call_b.split(":")[0] == call_x.split(":")[0]:
                            count_het += 1
                        else:
                            print(call_a, call_b, call_x)
                if count_het == 0:
                    het_percentage = 0.00
                else:
                    het_percentage = round(count_het/count_polymorphic*100, 2)
                if het_percentage == 100:
                    if score_a > 80 and score_b > 80:
                        verdict = 'successfulF1'
                    else:
                        verdict = 'PQF'
                if het_percentage >= CUTOFF_FONE and het_percentage < 100:
                    if score_a > 80 and score_b > 80:
                        verdict = 'possibleF1'
                    else:
                        verdict = 'PQF'
                if het_percentage < CUTOFF_FONE:
                    verdict = 'failedF1'
                score_a = consensus_report.at[parent_a, 'MeanScore']
                score_b = consensus_report.at[parent_b, 'MeanScore']
                pedver_summary.loc[len(pedver_summary.index)] = [f_one, parent_a, score_a, parent_b, score_b, count_polymorphic,
                                                                 count_het, count_parent_a, count_parent_b, count_missing,
                                                                 het_percentage, verdict]
    return pedver_summary
