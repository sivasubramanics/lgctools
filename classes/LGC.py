from utils.definitions import *
from collections import defaultdict
from classes.Data import *

class LGC():
    """docstring for ClassName."""
    def __init__(self, in_lgc_file):
        # self.name = name
        self.get_data(in_lgc_file)
        
    def get_data(self, in_lgc_file):
        flag_data = 0
        count_data_points = 0
        smdata = defaultdict()
        msdata = defaultdict()
        with open(in_lgc_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                entries = line.split(CSV)
                if entries[0] == "Title":
                    self.name = entries[1]
                if line == 'Data':
                    flag_data = 1
                    continue
                if flag_data == 1:
                    count_data_points += 1
                    if count_data_points == 1:
                        continue
                    marker = entries[6]
                    sample = entries[7]
                    gt_call = entries[3]
                    if sample in BLANK_SAMPLES:
                        continue
                    if not marker in msdata:
                        msdata[marker] = SM(marker)
                    if not sample in smdata:
                        smdata[sample] = MS(marker)
                    smdata[sample].put_data(marker, gt_call)
                    msdata[marker].put_data(sample, gt_call)
        self.smdata = smdata
        self.msdata = msdata        
                    
                    
