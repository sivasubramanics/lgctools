#!/usr/bin/env python

from classes.Hapmap import Hapmap
from plugins.bestmarkers import get_best_markers
from plugins.performance import check_performance
from plugins.plots import make_snp_plots
from classes.Data import Markermetadata
import os
from utils.file_processing import *
from plugins.summary import get_marker_summary, get_sample_summary
import sys
from classes.LGC import LGC
from utils.definitions import *
from collections import defaultdict

in_lgc_file = "data/Genotyping-008.077-08.csv"
gt_data = LGC(in_lgc_file)
msdata = gt_data.msdata
smdata = gt_data.smdata

marker_info = get_marker_info()
markers = defaultdict(Markermetadata)
samples = list(gt_data.smdata.keys())
markers = get_markers(in_lgc_file, markers, marker_info)

outdir = "output/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
print(f"Created output directory {outdir}")

# print("Writing Grid file..")
# write_grid_file("output/out_grid.csv", smdata, markers)

# print("Writing Flapjack file..")
# write_flapjack_file("output/out_FJ.data", smdata, markers)

# print("Writing Hapmap file..")
# write_hapmap_file("output/out.hmp.txt", smdata, markers)

# print("Writing SNP plots..")
# make_snp_plots(gt_data.msdata, markers, outdir, gt_data.name)

# print(gt_data.name, len(smdata.keys()), len(msdata.keys()))

print("Writing Marker Summary..")
marker_summary = get_marker_summary(msdata)
marker_summary.to_csv("output/marker_summary.txt", sep = "\t", index=False)

print("Writing Sample Summary..")
sample_summary = get_sample_summary(smdata)
sample_summary.to_csv("output/sample_summary.txt", sep = "\t", index=False)

# print("Checking performance..")
# performance = check_performance(gt_data.smdata)
# print(f"Total combinations: {performance[0]}")
# print(f"Combinations with ZERO polymorphic markers: {performance[1]}")
# print(f"Combinations with >= 1 polymorphic markers: {performance[2]}")
# print(f"Combinations with >= 2 polymorphic markers: {performance[3]}")

# smdata, msdata = fill_gaps_gtdata(smdata, msdata)
# all_summary,ind_summary = get_best_markers(smdata, msdata, marker_summary)
# all_summary.to_csv("output/BestMarkerSummaryAll.txt", sep = "\t", index=False)
# ind_summary.to_csv("output/BestMarkerSummaryInd.txt", sep = "\t", index=False)

filename = "output/BestMarkerSummaryInd.txt"
ind_summary = pd.read_csv(filename, sep="\t", index_col='marker_count')
count = 15
print(list(ind_summary.loc[count]['marker']))
extract_markers = list(ind_summary.loc[count]['marker'])
extract_markers = extract_dict(markers, extract_markers)

tmp_smdata = subset_gtdata(smdata, extract_markers, samples, 'samplefast')
tmp_msdata = subset_gtdata(msdata, extract_markers, samples, 'markerfast')


print("Writing Hapmap file..")
write_hapmap_file("output/15_out.hmp.txt", tmp_smdata, extract_markers)

gt_data = Hapmap("output/15_out.hmp.txt")
print("Checking performance..")
performance = check_performance(gt_data.smdata)
print(f"Total combinations: {performance[0]}")
print(f"Combinations with ZERO polymorphic markers: {performance[1]}")
print(f"Combinations with >= 1 polymorphic markers: {performance[2]}")
print(f"Combinations with >= 2 polymorphic markers: {performance[3]}")





