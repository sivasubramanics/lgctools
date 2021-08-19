from plugins.performance import *
import pandas as pd

def get_best_markers(smdata, msdata, marker_summary):
    """
    Returns dataframe with best markers summary
    """
    markers = list(msdata.keys())
    samples = list(smdata.keys())
    overall_summary = pd.DataFrame(columns = ['marker_count', 
                                                'total_combinations',
                                                'zero_count', 
                                                'min_one_count',
                                                'min_two_count'])
    markers_summary = pd.DataFrame(columns = ['marker_count',
                                                'marker', 
                                                'total_combinations',
                                                'zero_count', 
                                                'min_one_count',
                                                'min_two_count'])
    scores = []
    while len(markers) >= 54:
        scores = check_performance(smdata)
        scores.insert(0, len(markers))
        overall_summary.loc[len(overall_summary.index)] = scores
        # print(f"{len(markers)}\tALL\t{scores[3]/scores[0]*100}")
        marker_scores = defaultdict()
        for marker in markers:
            tmp_markers = markers.copy()
            tmp_markers.remove(marker)
            tmp_smdata = subset_gtdata(smdata, tmp_markers, samples, 'samplefast')
            tmp_msdata = subset_gtdata(msdata, tmp_markers, samples, 'markerfast')
            scores = check_performance(tmp_smdata)
            missing_percentage = marker_summary.loc[markers_summary.marker_name == marker,'missing_percentage'].values[0]
            marker_scores[marker] = [scores[3]/scores[0]*100, marker_summary.marker_name[marker]['missing_percentage']]
            scores.insert(0, marker)
            scores.insert(0, len(markers))
            # print(len(scores), scores)
            markers_summary.loc[len(markers_summary.index)] = scores
            
        marker_scores = sort_dict(marker_scores, True)
        badmarker = list(marker_scores.keys())[0]
        markers.remove(badmarker)
        smdata = subset_gtdata(smdata, markers, samples, 'samplefast')
        msdata = subset_gtdata(msdata, markers, samples, 'markerfast')
    return overall_summary, markers_summary


