#!/usr/bin/env python

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

write_grid_file("data/out_grid.csv", gt_data.smdata, markers)
write_flapjack_file("data/out_FJ.data", gt_data.smdata, markers)
write_hapmap_file("data/out.hmp.txt", gt_data.smdata, markers)
print(gt_data.name, len(gt_data.smdata.keys()), len(gt_data.msdata.keys()))
get_marker_summary(gt_data.msdata)
get_sample_summary(gt_data.smdata)








