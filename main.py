#!/usr/bin/env python

import os
from plugins.summary import get_marker_summary
import sys
from classes.LGC import LGC
from utils.definitions import *

in_lgc_file = "data/Genotyping-092.008-01.csv"
gt_data = LGC(in_lgc_file)

print(gt_data.name, len(gt_data.smdata.keys()), len(gt_data.msdata.keys()))
get_marker_summary(gt_data.msdata)








