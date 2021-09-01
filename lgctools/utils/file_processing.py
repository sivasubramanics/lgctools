from utils.definitions import *
from utils.utils import *
from collections import defaultdict
from classes.Data import Markermetadata
from classes.Pedigree import Pedigree
import sys
import os
import pandas as pd


def write_grid_file(outfile, smdata, markers, *args):
    """
    Writes genotype data in grid csv file format
    """
    if args:
        sample_dict = args[0]
    else:
        sample_dict = defaultdict()
    outfile_handle = open(outfile, 'w')
    outfile_handle.write('DNA \\ Assay')
    for marker in markers:
        outfile_handle.write(CSV + marker)
    outfile_handle.write(NEWLINE)
    for sample in smdata:
        if sample in sample_dict:
            outfile_handle.write(sample_dict[sample])
        else:
            outfile_handle.write(sample)
        for marker in markers:
            if marker in smdata[sample].data:
                outfile_handle.write(
                    CSV + smdata[sample].data[marker].__str__())
            else:
                outfile_handle.write(CSV + 'N:N')
        outfile_handle.write(NEWLINE)
    outfile_handle.close()


def write_flapjack_file(outfile, smdata, markers, *args):
    """
    Writes genotype data in flapjack file format
    """
    if args:
        sample_dict = args[0]
    else:
        sample_dict = defaultdict()
    outfile_handle = open(outfile, 'w')
    outfile_handle.write('# fjFile = GENOTYPE' + NEWLINE)
    # outfile_handle.write(TAB)
    for marker in markers:
        outfile_handle.write(TAB + marker)
    outfile_handle.write(NEWLINE)
    for sample in smdata:
        if sample in sample_dict:
            outfile_handle.write(sample_dict[sample])
        else:
            outfile_handle.write(sample)
        for marker in markers:
            if marker in smdata[sample].data:
                outfile_handle.write(
                    TAB + to_flapjack(smdata[sample].data[marker].__str__()))
            else:
                outfile_handle.write(TAB + 'N/N')
        outfile_handle.write(NEWLINE)
    outfile_handle.close()


def write_hapmap_file(outfile, sm_data, markers, *args):
    """
    Writes genotype data in hapmap file format
    """
    hmp_data = defaultdict()
    if args:
        sample_dict = args[0]
    else:
        sample_dict = defaultdict()
    outfile_handle = open(outfile, 'w')
    head_line = list(HMPHEAD)
    for marker in markers:
        hmp_data[marker] = []
        hmp_data[marker].append(marker)
        hmp_data[marker].append(to_flapjack(markers[marker].get_allele_xy()))
        hmp_data[marker].append(markers[marker].get_chr())
        hmp_data[marker].append(markers[marker].get_position())
        hmp_data[marker].append('+')
        hmp_data[marker] += HMPNA
    for sample in sm_data:
        if sample in sample_dict:
            head_line.append(sample_dict[sample].get_sample_name())
        else:
            head_line.append(sample)
        for marker in markers:
            if marker in sm_data[sample].data:
                hmp_data[marker].append(
                    to_hmp(sm_data[sample].data[marker].__str__()))
            else:
                hmp_data[marker].append('N/N')
    # outfile_handle.write(TAB.join(head_line) + NEWLINE)
    hapmap = pd.DataFrame(columns=head_line)
    for marker in hmp_data:
        hapmap.loc[len(hapmap.index)] = hmp_data[marker]
        # outfile_handle.write(TAB.join(hmp_data[marker]) + NEWLINE)
    hapmap.chrom = pd.to_numeric(hapmap.chrom, errors='coerce')
    hapmap.pos = pd.to_numeric(hapmap.pos, errors='coerce')
    hapmap = hapmap.sort_values(by=['chrom', 'pos'])
    hapmap.to_csv(outfile, sep="\t", index=False)
    outfile_handle.close()


def make_markers(msdata, markers, marker_info_dict):
    for marker in msdata:
        markers[marker] = Markermetadata(marker)
        if marker in marker_info_dict:
            markers[marker].put_marker_info(marker_info_dict[marker])
    return markers


def get_markers(in_lgc_file, markers, marker_info_dict):
    """
    Returns list<Markermetadata> from lgc file
    """
    noData = 0
    dataFlag = 0
    noMarkers = 0
    flag_snp = 0
    ms_for_plot = defaultdict()
    markers_for_plot = defaultdict()
    with open(in_lgc_file) as fh:
        noSamples = 0
        noBlankSamples = 0
        list_samples = []
        for line in fh:
            line = line.strip()
            entries = line.split(CSV)
            if line == "":
                flag_snp = 0
                continue
            if entries[0] == 'SNPs':
                flag_snp = 1
                continue
            if flag_snp == 1:
                noMarkers += 1
                if noMarkers == 1:
                    continue
                # lineEntries = line.split(CSV)
                if entries[0] not in markers:
                    markers[entries[0]] = Markermetadata(entries[0])
                    if entries[0] in marker_info_dict:
                        markers[entries[0]].put_marker_info(
                            marker_info_dict[entries[0]])
                if entries[0] not in markers_for_plot:
                    markers_for_plot[entries[0]] = Markermetadata(entries[0])
                markers[entries[0]].put_alleles(line)
                markers_for_plot[entries[0]].put_alleles(line)
    return markers


def get_marker_info():
    """
    Returns dictionary of marker info holding physical positions.
    """
    lineNo = 0
    marker_info_dict = defaultdict()
    tmp_path = os.path.abspath(sys.argv[0])
    cur_dir, file_name = os.path.split(tmp_path)
    with open(cur_dir + '/../tests/data/' + 'KASP_MarkerInfo.txt') as f:
        for line in f:
            line = line.strip()
            lineNo += 1
            if lineNo == 1:
                continue
            lineEntries = line.split(TAB)
            if lineEntries[0] in marker_info_dict:
                continue
            else:
                marker_info_dict[lineEntries[0]] = []
            if len(lineEntries) == 3:
                marker_info_dict[lineEntries[0]].append(lineEntries[1])
                marker_info_dict[lineEntries[0]].append(lineEntries[2])
    return marker_info_dict


def get_pedigree(in_pedigree_file):
    pedigree_dict = defaultdict()
    line_count = 0
    with open(in_pedigree_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            line_count += 1
            if line_count == 1:
                continue
            entries = line.split(TAB)
            if len(entries) < 4:
                continue
            if entries[0] in pedigree_dict:
                continue
            if not entries[0] in pedigree_dict:
                pedigree_dict[entries[0]] = Pedigree(entries[0])
                pedigree_dict[entries[0]].set_designation(entries[1])
                pedigree_dict[entries[0]].set_parent_a(entries[2])
                pedigree_dict[entries[0]].set_parent_b(entries[3])
    fh.close()
    print_log(
        f"Pedigree dictionary contains {len(pedigree_dict)} valid F-One genotypes with parents.")
    return pedigree_dict


def get_samplemap(in_samplemap_file):
    samplemap_dict = defaultdict()
    samples = []
    line_count = 0
    with open(in_samplemap_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            line_count += 1
            if line_count == 1:
                continue
            entries = line.split(TAB)
            if len(entries) < 2:
                continue
            if entries[0] in samplemap_dict:
                continue
            if entries[1] in samples:
                continue
            if not entries[0] in samplemap_dict:
                samplemap_dict[entries[0]] = entries[1]
                samples.append(entries[1])
    fh.close()
    print_log(
        f"Sample map dictionary contains {len(samplemap_dict)} valid genotypes with reference.")
    return samplemap_dict
