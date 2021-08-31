from plugins.performance import *
from wordcloud import WordCloud, STOPWORDS
import pandas as pd
import sys
from multiprocessing import Process
import multiprocessing
from utils.definitions import *
from utils.utils import *


def make_wordcloud(best_marker_summary, outfile):
    text = " ".join(cat.split()[0] for cat in best_marker_summary.marker)
    stopwords = set(STOPWORDS)
    wordcloud = WordCloud(width=800, height=800,
                          background_color='white',
                          stopwords=stopwords,
                          prefer_horizontal=0.2,
                          min_font_size=10).generate(text)
    wordcloud.to_file(outfile)
    return


def get_bad_marker(marker_scores):
    tmp_dict_one = defaultdict()
    tmp_dict_two = defaultdict()
    for key in marker_scores:
        tmp_dict_one[key] = marker_scores[key][0]
        tmp_dict_two[key] = marker_scores[key][1]
    tmp_dict_one = sort_dict(tmp_dict_one, True)
    badmarker = list(tmp_dict_one.keys())[0]
    identical_markers = get_dup_keys(tmp_dict_one, tmp_dict_one[badmarker])
    if len(identical_markers) > 1:
        tmp_dict_three = extract_dict(tmp_dict_two, identical_markers)
        tmp_dict_three = sort_dict(tmp_dict_three, False)
        badmarker = list(tmp_dict_three.keys())[0]
    return badmarker


def get_best_markers(smdata, msdata, marker_summary, sample_list_a=[], sample_list_b=[]):
    """
    Returns dataframe with best markers summary
    """
    marker_summary = marker_summary.set_index('marker_name')
    markers = list(msdata.keys())
    samples = list(smdata.keys())
    overall_summary = pd.DataFrame(columns=['marker_count',
                                            'total_combinations',
                                            'zero_count',
                                            'min_one_count',
                                            'min_two_count'])
    best_markers_summary = pd.DataFrame(columns=['marker_count',
                                                 'marker',
                                                 'total_combinations',
                                                 'zero_count',
                                                 'min_one_count',
                                                 'min_two_count'])
    scores = []
    while len(markers) >= 5:
        scores = check_performance(smdata, sample_list_a, sample_list_b)
        scores.insert(0, len(markers))
        overall_summary.loc[len(overall_summary.index)] = scores
        if MULTIPROCESSING:
            manager = multiprocessing.Manager()
            marker_scores_dict = manager.dict()
            df_ms = pd.DataFrame(columns=['marker_name', 'score', 'qual'])
            processess = []
            no_marker = 0
            for marker in markers:
                no_marker += 1
                tmp_markers = markers.copy()
                tmp_markers.remove(marker)
                # print(f"\r{len(markers)} - {no_marker}")
                print(f"\rStatus: {len(markers)} marker - {no_marker} iteration",
                      end='   ', flush=True)
                sys.stdout.flush()
                tmp_smdata = subset_gtdata(
                    smdata, tmp_markers, samples, 'samplefast')
                # tmp_msdata = subset_gtdata(msdata, tmp_markers, samples, 'markerfast')
                process = Process(target=check_performance, args=(
                    tmp_smdata, sample_list_a, sample_list_b, marker_scores_dict, marker))
                processess.append(process)
                process.start()
            for process in processess:
                process.join()
            for marker in marker_scores_dict:
                scores = marker_scores_dict[marker]
                pic = marker_summary.at[marker, 'PIC']
                df_ms.loc[len(df_ms.index)] = [
                    marker, scores[3]/scores[0]*100, pic]
                scores.insert(0, marker)
                scores.insert(0, len(markers))
                best_markers_summary.loc[len(
                    best_markers_summary.index)] = scores
        else:
            for marker in markers:
                tmp_markers = markers.copy()
                tmp_markers.remove(marker)
                print_log(f"\r{len(markers)} - {no_marker}",
                          end='', flush=True)
                sys.stdout.flush()
                tmp_smdata = subset_gtdata(
                    smdata, tmp_markers, samples, 'samplefast')
                # tmp_msdata = subset_gtdata(msdata, tmp_markers, samples, 'markerfast')
                scores = check_performance(
                    tmp_smdata, sample_list_a, sample_list_b)
                pic = marker_summary.at[marker, 'PIC']
                df_ms.loc[len(df_ms.index)] = [
                    marker, scores[3]/scores[0]*100, pic]
                scores.insert(0, marker)
                scores.insert(0, len(markers))
                best_markers_summary.loc[len(
                    best_markers_summary.index)] = scores

        df_ms = df_ms.sort_values(
            by=['score', 'qual'], ascending=[True, False])
        badmarker = df_ms['marker_name'].iloc[-1]
        # print(len(markers), 'dict', badmarker)
        markers.remove(badmarker)
        smdata = subset_gtdata(smdata, markers, samples, 'samplefast')
        msdata = subset_gtdata(msdata, markers, samples, 'markerfast')
    best_markers_summary = best_markers_summary.sort_values(
        by=['marker_count', 'min_two_count'], ascending=[False, True])
    print()
    print("\033[A\033[A")
    return overall_summary, best_markers_summary
