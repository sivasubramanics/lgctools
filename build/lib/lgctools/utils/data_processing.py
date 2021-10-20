from collections import defaultdict
from ..utils.utils import print_log
from ..classes.Data import MS, SM, Call, Markermetadata


def subset_gtdata(in_data, markers, samples, data_type):
    miss_marker = []
    miss_sample = []
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
        for sample in samples:
            if not sample in out_data:
                miss_sample.append(sample)
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
        for marker in markers:
            if not marker in out_data:
                miss_marker.append(marker)
    else:
        print(
            f"Error: {data_type} is unknown. Should be 'samplefast' or 'markerfast'..")
        exit(1)
    if miss_marker:
        print_log(f"missing markers {miss_marker}")
    if miss_sample:
        print_log(f"missing samples {miss_sample}")
    return out_data


def fill_gaps_gtdata(smdata, msdata):
    markers = []
    samples = []
    missing_call = Call('N:N')
    markers = list(msdata.keys())
    samples = list(smdata.keys())

    if not markers:
        print(f"ERROR: Input markerfast dictionary is empty.")
        exit(1)

    if not samples:
        print(f"ERROR: Input samplefast dictionary is empty.")
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


def merge_dictionary(*args):
    out_dict = defaultdict()
    for i in range(0, len(args)):
        out_dict.update(args[i])
    return out_dict


def merge_gtdata(data_one, data_two, data_type):
    if data_type == "marker-fast":
        msdata = defaultdict()
        for marker in data_one:
            if not marker in msdata:
                msdata[marker] = MS(marker)
            for sample in data_one[marker].data:
                msdata[marker].put_data(sample, data_one[marker].data[sample])
        for marker in data_two:
            if not marker in msdata:
                msdata[marker] = MS(marker)
            for sample in data_two[marker].data:
                msdata[marker].put_data(sample, data_two[marker].data[sample])
        return msdata
    elif data_type == "sample-fast":
        smdata = defaultdict()
        for sample in data_one:
            if not sample in smdata:
                smdata[sample] = SM(sample)
            for marker in data_one[sample].data:
                smdata[sample].put_data(marker, data_one[sample].data[marker])
        for sample in data_two:
            if not sample in smdata:
                smdata[sample] = SM(sample)
            for marker in data_two[sample].data:
                smdata[sample].put_data(marker, data_two[sample].data[marker])
        return smdata
    else:
        print(
            f"Error: {data_type} is unknown. Should be 'sample-fast' or 'marker-fast'..")
        exit(1)


def filter_markers(markers, flt_marker_summary):
    flt_markers = list(flt_marker_summary['marker_name'])
    new_markers = defaultdict(Markermetadata)
    for marker in markers:
        if marker in flt_markers:
            new_markers[marker] = markers[marker]
    markers = new_markers.copy()
    return markers


def filter_samples(samples, flt_sample_summary):
    flt_samples = list(flt_sample_summary['sample_name'])
    new_samples = []
    for sample in samples:
        if sample in flt_samples:
            new_samples.append(sample)
    samples = new_samples.copy()
    return samples
