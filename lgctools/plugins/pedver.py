from utils.file_processing import *
import os


def process_pedver(options, smdata):
    out_prefix = options.out_prefix
    pedigree_file = options.pedigree_file

    pedver_summary = pd.DataFrame(
        columns=['F_One_Name',
                 'Parent_A',
                 'Parent_A_Score',
                 'Parent_B',
                 'Parent_B_Score',
                 'TotalMarkers',
                 'Polymorphic',
                 'Call_Hetero',
                 'Call_ParentA',
                 'Call_ParentB',
                 'Call_Missing',
                 'Percentage_Het',
                 'Comment'])
    pedigree = get_pedigree(pedigree_file)
    if not os.path.isfile(out_prefix + '_PurityReport.txt'):
        print_log(
            f"ERROR: Consensus report file ({out_prefix + '_PurityReport.txt'}) is not accessible. Something is wrong. Quiting...")
        exit(1)
    consensus_report = pd.read_csv(
        out_prefix + '_PurityReport.txt', sep=TAB, index_col='Name')
    for f_one in pedigree:
        parent_a = pedigree[f_one].get_parent_a()
        parent_b = pedigree[f_one].get_parent_b()
        if f_one in smdata:
            if parent_a in smdata and parent_b in smdata:
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
                    count_total = len(markers_x)
                    score_a = consensus_report.at[parent_a, 'MeanScore']
                    score_b = consensus_report.at[parent_b, 'MeanScore']
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
                        het_percentage = round(
                            count_het/count_polymorphic*100, 2)
                    if het_percentage == 100:
                        if score_a > 80 and score_b > 80:
                            verdict = 'successfulF1'
                        else:
                            verdict = 'PQF'
                    if het_percentage >= options.cutoff_fone and het_percentage < 100:
                        if score_a > 80 and score_b > 80:
                            verdict = 'possibleF1'
                        else:
                            verdict = 'PQF'
                    if het_percentage < options.cutoff_fone:
                        verdict = 'failedF1'
                    pedver_summary.loc[len(pedver_summary.index)] = [f_one, parent_a, score_a, parent_b, score_b, count_total, count_polymorphic,
                                                                     count_het, count_parent_a, count_parent_b, count_missing,
                                                                     het_percentage, verdict]
            if not parent_a in smdata or not parent_b in smdata:
                count_parent_a = 'na'
                count_parent_b = 'na'
                score_a = 'na'
                score_b = 'na'
                if parent_a in smdata:
                    count_parent_a = 0
                    score_a = consensus_report.at[parent_a, 'MeanScore']
                if parent_b in smdata:
                    count_parent_b = 0
                    score_b = consensus_report.at[parent_b, 'MeanScore']
                count_polymorphic = 'na'
                count_het = 0
                count_missing = 0
                verdict = 'failedF1'
                count_total = len(list(smdata[f_one].data.keys()))
                for marker in smdata[f_one].data:
                    call_x = smdata[f_one].data[marker].__str__()
                    call_a = ''
                    call_b = ''
                    if parent_a in smdata:
                        call_a = smdata[parent_a].data[marker].__str__()
                    if parent_b in smdata:
                        call_b = smdata[parent_b].data[marker].__str__()
                    if is_hetero(call_x):
                        count_het += 1
                    elif call_x in MISSING_CALLS:
                        count_missing += 1
                    elif call_x == call_a:
                        count_parent_a += 1
                    elif call_x == call_b:
                        count_parent_b += 1
                if count_het == 0:
                    het_percentage = 0.00
                else:
                    het_percentage = round(
                        count_het/(count_total-count_missing)*100, 2)
                if het_percentage == 100:
                    verdict = 'successfulF1'
                if het_percentage >= options.cutoff_fone and het_percentage < 100:
                    verdict = 'possibleF1'
                pedver_summary.loc[len(pedver_summary.index)] = [f_one, parent_a, score_a, parent_b, score_b, count_total, count_polymorphic,
                                                                 count_het, count_parent_a, count_parent_b, count_missing,
                                                                 het_percentage, verdict]
    return pedver_summary
