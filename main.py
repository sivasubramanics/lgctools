#!/usr/bin/env python

from plugins.plots import make_snp_plots
from classes.Data import Markermetadata
import os
from utils.file_processing import *
from plugins.summary import get_marker_summary, get_sample_summary
import sys
from classes.LGC import LGC
from utils.definitions import *
from collections import defaultdict

in_lgc_file = "data/Genotyping-092.008-01.csv"
gt_data = LGC(in_lgc_file)

marker_info = get_marker_info()
markers = defaultdict(Markermetadata)
samples = list(gt_data.smdata.keys())
markers = get_markers(in_lgc_file, markers, marker_info)

outdir = "output/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
print(f"Created output directory {outdir}")


# write_grid_file("output/out_grid.csv", gt_data.smdata, markers)
# write_flapjack_file("output/out_FJ.data", gt_data.smdata, markers)
# write_hapmap_file("output/out.hmp.txt", gt_data.smdata, markers)
make_snp_plots(gt_data.msdata, markers, outdir, gt_data.name)
print(gt_data.name, len(gt_data.smdata.keys()), len(gt_data.msdata.keys()))
marker_summary = get_marker_summary(gt_data.msdata)
marker_summary.to_csv("output/marker_summary.txt", sep = "\t", index=False)
sample_summary = get_sample_summary(gt_data.smdata)
sample_summary.to_csv("output/sample_summary.txt", sep = "\t", index=False)








